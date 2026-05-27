//! Bit-parallel maximal reverse-complementary match enumeration on
//! cyclic strings.
//!
//! Given a set of cyclic angle sequences ("tiles"), returns every
//! maximal reverse-complementary substring match (length ≥ 1) between
//! a boundary and any tile. Designed for small alphabets (4, 10, 12
//! distinct symbols depending on cyclotomic ring) and the
//! "fresh boundary against a fixed tileset" call pattern: the tileset
//! masks are precomputed once and reused across every boundary stream.
//!
//! # Algorithm
//!
//! For each tile `T_j` of length `L_j`, the matcher targets
//! `P_j = rc(T_j)` (reverse, negate). Forward matching of the boundary
//! against `P_j` is RC matching against `T_j`. Per `(tile, symbol)` we
//! precompute one `L_j`-bit mask `M_j[c]` with bit `p` set iff
//! `P_j[p] == c`.
//!
//! Per tile during a query we maintain:
//!
//! * `alive_prev`: the bitset of `P_j` end-positions where some match
//!   ended at the previous boundary character.
//! * `len[ep]`: length of the *longest* match ending at `(prev_bp, ep)`.
//!
//! Per boundary character `c = B[k]`:
//!
//! ```text
//! alive_curr = M_j[c]
//! shifted    = rotate_left_1(alive_prev)
//! continues  = shifted & alive_curr     // length ≥ 2 ends at ep
//! deaths     = shifted & !alive_curr    // a match died at the previous step
//! starts     = alive_curr & !continues  // fresh length-1 starts
//! ```
//!
//! `deaths` triggers emission of the dying right-maximal match at the
//! previous boundary position with length `len[(ep − 1) mod L_j]`. The
//! `len` array is updated: continuations increment the length of their
//! predecessor; starts set length 1; everything else resets to 0.
//!
//! # Cyclic boundary
//!
//! Cyclic boundaries are handled by streaming `B · B` (length `2N`) and
//! deduplicating by `(pos_a, pos_b)`, keeping the longest length per
//! pair.
//!
//! # Full-tile (cap) matches
//!
//! A match cannot meaningfully exceed `cap = min(N, L_j)`. When a
//! continuation would push length past `cap` we emit at length `cap`
//! and keep the bit alive, so that subsequent cyclic rotations of the
//! same full-tile match are also reported. Dedup canonicalises.

use crate::geom::matches::{EdgeRange, Segment, TileMatch};

/// Fixed-length cyclic bitset over `len` bits, stored as a `Vec<u64>`.
#[derive(Clone)]
struct CyclicBitset {
    bits: Vec<u64>,
    len: usize,
}

impl CyclicBitset {
    fn new(len: usize) -> Self {
        let words = len.div_ceil(64);
        CyclicBitset {
            bits: vec![0u64; words],
            len,
        }
    }

    fn set_bit(&mut self, p: usize) {
        debug_assert!(p < self.len);
        self.bits[p / 64] |= 1u64 << (p % 64);
    }

    /// Mask off any bits at positions ≥ `self.len` in the top word.
    fn mask_top(&mut self) {
        if self.bits.is_empty() {
            return;
        }
        let trailing = self.len % 64;
        if trailing != 0 {
            let mask = (1u64 << trailing) - 1;
            let last = self.bits.len() - 1;
            self.bits[last] &= mask;
        }
    }

    /// Cyclic shift left by 1 over `self.len` bits.
    fn rotate_left_1(&self) -> Self {
        if self.len == 0 {
            return self.clone();
        }
        let mut out = CyclicBitset::new(self.len);
        let mut carry: u64 = 0;
        for i in 0..self.bits.len() {
            out.bits[i] = (self.bits[i] << 1) | carry;
            carry = self.bits[i] >> 63;
        }
        // The top bit at position (len - 1) wraps to position 0.
        let wrap_pos = self.len - 1;
        let wrap_bit = (self.bits[wrap_pos / 64] >> (wrap_pos % 64)) & 1;
        if wrap_bit == 1 {
            out.bits[0] |= 1;
        }
        out.mask_top();
        out
    }

    fn and(&self, other: &Self) -> Self {
        debug_assert_eq!(self.len, other.len);
        let mut out = CyclicBitset::new(self.len);
        for i in 0..self.bits.len() {
            out.bits[i] = self.bits[i] & other.bits[i];
        }
        out
    }

    fn and_not(&self, other: &Self) -> Self {
        debug_assert_eq!(self.len, other.len);
        let mut out = CyclicBitset::new(self.len);
        for i in 0..self.bits.len() {
            out.bits[i] = self.bits[i] & !other.bits[i];
        }
        out.mask_top();
        out
    }

    /// Collect set-bit positions into a Vec (one allocation per call;
    /// the inner streaming loop tolerates this for clarity).
    fn set_bits(&self) -> Vec<usize> {
        let mut out = Vec::new();
        for (i, &word) in self.bits.iter().enumerate() {
            let mut w = word;
            while w != 0 {
                let bit = w.trailing_zeros() as usize;
                out.push(i * 64 + bit);
                w &= w - 1;
            }
        }
        out
    }
}

/// Bit-parallel cyclic RC matcher over a fixed set of tile sequences.
///
/// Indexes a collection of cyclic angle sequences and answers
/// `maximal_rc_matches(i, j)` between any two of them, or
/// `stream_boundary(b, j)` / `stream_boundary_all_tiles(b)` for an
/// arbitrary external boundary against the indexed tileset.
/// Sentinel for "boundary symbol not in this matcher's alphabet".
/// Picked as `u8::MAX` so a single `if idx == NO_SYM` check in the
/// hot loop suffices.
const NO_SYM: u8 = u8::MAX;

pub struct BitParallelMatcher {
    tiles: Vec<Vec<i8>>,
    /// Flat `i8 -> sym_idx` table indexed by `c as u8` (256 entries).
    /// `NO_SYM` marks chars outside the alphabet. The alphabet is
    /// small (4-12 symbols for the cyclotomic rings we care about),
    /// so the 256-byte table is cheaper to allocate and (critically)
    /// faster to lookup in the inner streaming loop than a HashMap.
    sym_to_idx: [u8; 256],
    /// `masks[tile_idx][sym_idx]` = positions in `rc(tile)` holding
    /// that symbol.
    masks: Vec<Vec<CyclicBitset>>,
}

/// Per-tile streaming state, reused across every boundary step.
///
/// One per `(tile, sweep)` pair: holds the bitset of live match
/// end-positions, per-position length arrays, a precomputed empty
/// mask (for boundary symbols not in this tile's alphabet), and the
/// emissions accumulated so far.
struct TileState {
    n_b: usize,
    cap: usize,
    alive_prev: CyclicBitset,
    len_arr: Vec<u8>,
    len_new: Vec<u8>,
    empty: CyclicBitset,
    // (boundary_end_in_stream, tile_end_in_P, len)
    emissions: Vec<(usize, usize, usize)>,
    skip: bool,
}

impl BitParallelMatcher {
    /// Build the matcher by precomputing per-tile symbol masks over
    /// `rc(tile)`.
    pub fn new(tiles: &[Vec<i8>]) -> Self {
        assert!(!tiles.is_empty(), "need at least one tile");

        // Alphabet must include every char that can appear either in a
        // boundary lookup or as a P_j mask value. Boundary chars range
        // over arbitrary inputs, P_j chars are `-chars(tiles)`. Take
        // the union of both so external boundaries with chars in
        // `-chars(tiles)` (e.g. `-1` when only `1` appears in tiles)
        // still resolve to the right mask.
        let mut alpha: Vec<i8> = Vec::new();
        for t in tiles.iter() {
            for &c in t.iter() {
                alpha.push(c);
                alpha.push(-c);
            }
        }
        alpha.sort();
        alpha.dedup();
        assert!(
            alpha.len() < NO_SYM as usize,
            "BitParallelMatcher: alphabet size {} exceeds {}",
            alpha.len(),
            NO_SYM as usize - 1
        );
        let mut sym_to_idx = [NO_SYM; 256];
        for (i, &c) in alpha.iter().enumerate() {
            sym_to_idx[c as u8 as usize] = i as u8;
        }
        let num_syms = alpha.len();

        let masks: Vec<Vec<CyclicBitset>> = tiles
            .iter()
            .map(|t| {
                let l = t.len();
                let mut m: Vec<CyclicBitset> =
                    (0..num_syms).map(|_| CyclicBitset::new(l)).collect();
                // P_j[p] = -t[L - 1 - p]
                for p in 0..l {
                    let c = -t[l - 1 - p];
                    let idx = sym_to_idx[c as u8 as usize];
                    if idx != NO_SYM {
                        m[idx as usize].set_bit(p);
                    }
                }
                m
            })
            .collect();

        BitParallelMatcher {
            tiles: tiles.to_vec(),
            sym_to_idx,
            masks,
        }
    }

    /// Lookup of boundary char `c` in the matcher's alphabet. Returns
    /// the symbol index, or `None` for chars outside the alphabet
    /// (which kill every live match at the streaming step).
    #[inline(always)]
    fn sym_idx_of(&self, c: i8) -> Option<usize> {
        let idx = self.sym_to_idx[c as u8 as usize];
        if idx == NO_SYM {
            None
        } else {
            Some(idx as usize)
        }
    }

    pub fn num_tiles(&self) -> usize {
        self.tiles.len()
    }

    pub fn tile_len(&self, i: usize) -> usize {
        self.tiles[i].len()
    }

    /// All maximal cyclic RC matches between tile `i` and tile `j`,
    /// including length-1 matches. Caller is responsible for any
    /// downstream geometric filtering (e.g. single-edge angle checks).
    pub fn maximal_rc_matches(&self, i: usize, j: usize) -> Vec<TileMatch> {
        let boundary = self.tiles[i].clone();
        let mut matches = self.stream_boundary(&boundary, j);
        for m in &mut matches {
            m.a.tile_id = i;
        }
        matches
    }

    /// Same as [`Self::maximal_rc_matches`] but **post-filtered** to
    /// keep only matches whose A-side starting position is in
    /// `positions`. The cost is unchanged from `maximal_rc_matches` --
    /// the matcher still does a full sweep, the filter only shrinks
    /// the returned `Vec`. Callers with very small `positions` against
    /// large tiles won't see a speedup; this exists purely for
    /// caller convenience.
    pub fn maximal_rc_matches_at_positions(
        &self,
        i: usize,
        j: usize,
        positions: &[usize],
    ) -> Vec<TileMatch> {
        if positions.is_empty() {
            return vec![];
        }
        let n_a = self.tiles[i].len();
        let mut keep = vec![false; n_a];
        for &p in positions {
            if p < n_a {
                keep[p] = true;
            }
        }
        self.maximal_rc_matches(i, j)
            .into_iter()
            .filter(|m| keep[m.a.range.start_offset])
            .collect()
    }

    /// Allocate a fresh per-tile streaming state for `tile_j` against
    /// a boundary of length `n_a`.
    fn init_tile_state(&self, tile_j: usize, n_a: usize) -> TileState {
        let n_b = self.tiles[tile_j].len();
        let cap = n_a.min(n_b);
        TileState {
            n_b,
            cap,
            alive_prev: CyclicBitset::new(n_b),
            len_arr: vec![0u8; n_b],
            len_new: vec![0u8; n_b],
            empty: CyclicBitset::new(n_b),
            emissions: Vec::new(),
            skip: n_b == 0,
        }
    }

    /// Advance one step of the bit-parallel sweep for a single tile.
    ///
    /// `k` is the position in the doubled boundary stream (`0..2*n_a`).
    /// `sym_idx` is the index of the current boundary character in
    /// `self.sym_to_idx`, or `None` if the character is foreign to the
    /// matcher's alphabet (then no `P_j` position can match this char,
    /// so every live bit dies).
    fn step_tile(
        &self,
        ts: &mut TileState,
        tile_j: usize,
        k: usize,
        sym_idx: Option<usize>,
    ) {
        // `alive_curr`: bits of `P_j` that match the current boundary
        // char. When `sym_idx` is `None` we borrow `ts.empty` so the
        // remainder of the algorithm uniformly works against a
        // `&CyclicBitset` (everything dies, as it should).
        let alive_curr: &CyclicBitset = match sym_idx {
            Some(idx) => &self.masks[tile_j][idx],
            None => &ts.empty,
        };

        let shifted = ts.alive_prev.rotate_left_1();
        let continues = shifted.and(alive_curr);
        let deaths = shifted.and_not(alive_curr);
        let starts = alive_curr.and_not(&continues);

        // Deaths: previous-step matches that fail to extend.
        // `ep_new` in `deaths` was reached by rotating `ep_old =
        // (ep_new - 1) mod n_b`; emit the match that ended there at
        // boundary position `k - 1`.
        for ep_new in deaths.set_bits() {
            let ep_old = if ep_new == 0 { ts.n_b - 1 } else { ep_new - 1 };
            let l = ts.len_arr[ep_old] as usize;
            if l > 0 {
                ts.emissions.push((k - 1, ep_old, l));
            }
        }

        // Reset new-length buffer.
        for v in ts.len_new.iter_mut() {
            *v = 0;
        }

        // Continuations: extend length from predecessor. Cap at `cap =
        // min(n_a, n_b)`; when we hit cap we emit at full length but
        // keep the bit alive, so subsequent rotations of the same
        // full-tile match are also reported (the dedup pass keeps the
        // canonical one).
        for ep in continues.set_bits() {
            let ep_pred = if ep == 0 { ts.n_b - 1 } else { ep - 1 };
            let nl = (ts.len_arr[ep_pred] as usize + 1).min(ts.cap);
            ts.len_new[ep] = nl as u8;
            if nl == ts.cap {
                ts.emissions.push((k, ep, ts.cap));
            }
        }

        // Fresh length-1 starts (where no continuation lives). Emit
        // immediately when `cap == 1` (no further extension possible).
        for ep in starts.set_bits() {
            ts.len_new[ep] = 1;
            if ts.cap == 1 {
                ts.emissions.push((k, ep, 1));
            }
        }

        ts.alive_prev = alive_curr.clone();
        std::mem::swap(&mut ts.len_arr, &mut ts.len_new);
    }

    /// Convert raw `(b_end, ep, l)` emissions into `TileMatch`es:
    /// reconstruct `(pos_a, pos_b)` on the original tiles, drop
    /// non-left-maximal matches, sort by `(pos_a, pos_b)` keeping the
    /// longest per pair.
    fn emissions_to_matches(
        &self,
        emissions: Vec<(usize, usize, usize)>,
        boundary: &[i8],
        tile_j: usize,
    ) -> Vec<TileMatch> {
        let n_a = boundary.len();
        let n_b = self.tiles[tile_j].len();
        let cap = n_a.min(n_b);
        let tile_b_seq: &[i8] = &self.tiles[tile_j];

        let mut matches: Vec<TileMatch> = emissions
            .into_iter()
            .filter_map(|(b_end, ep, l)| {
                if l == 0 || l > cap {
                    return None;
                }
                // Decode (b_end, ep, l) in stream/P coordinates back
                // to (pos_a, pos_b) on the original tiles. `b_end` is
                // the boundary position where the match ended; the
                // match covered `l` boundary edges, so its start is
                // `b_end - (l - 1)` (mod n_a), i.e.
                // `(b_end + n_a + 1 - l) % n_a`. Symmetrically on the
                // P-side: tile_start_in_p = ep - (l - 1) (mod n_b).
                // Tile position k in P = (n_b - 1 - k) on the original
                // tile, so pos_b = n_b - 1 - tile_start_in_p (mod n_b).
                let pos_a = (b_end + n_a + 1 - l) % n_a;
                let tile_start_in_p = (ep + n_b + 1 - l) % n_b;
                let pos_b = (n_b + n_b - 1 - tile_start_in_p) % n_b;
                // Left-maximality: drop matches whose left endpoint
                // could extend (= right-suffixes of a longer match
                // starting one boundary char earlier). Matches that
                // wrap a whole boundary or tile are trivially
                // left-maximal.
                if l < n_a && l < n_b {
                    let prev_a = boundary[(pos_a + n_a - 1) % n_a];
                    let prev_b_rc = -tile_b_seq[(pos_b + 1) % n_b];
                    if prev_a == prev_b_rc {
                        return None;
                    }
                }
                Some(TileMatch::new(
                    Segment::new(0, EdgeRange::new(pos_a, l)),
                    Segment::new(tile_j, EdgeRange::new(pos_b, l)),
                ))
            })
            .collect();

        // Sort with `len` descending as the tiebreaker so that the
        // following `dedup_by` (which keeps the FIRST of consecutive
        // equals) keeps the longest match per `(pos_a, pos_b)` pair.
        matches.sort_by(|a, b| {
            a.a.range
                .start_offset
                .cmp(&b.a.range.start_offset)
                .then_with(|| a.b.range.start_offset.cmp(&b.b.range.start_offset))
                .then_with(|| b.len().cmp(&a.len()))
        });
        matches.dedup_by(|a, b| {
            a.a.range.start_offset == b.a.range.start_offset
                && a.b.range.start_offset == b.b.range.start_offset
        });
        matches
    }

    /// Single sweep through `boundary`, updating per-tile state for
    /// **every** tile in this matcher in lockstep. Returns all maximal
    /// cyclic RC matches across every tile. This is the natural BP
    /// usage pattern: one 2N-step pass over the boundary amortises
    /// setup, allocation, and the boundary-char lookup across all
    /// `n_tiles` tiles, instead of paying those costs `n_tiles` times.
    ///
    /// `TileMatch::tile_a` is set to 0; callers can override (the
    /// boundary is opaque to the matcher).
    pub fn stream_boundary_all_tiles(&self, boundary: &[i8]) -> Vec<TileMatch> {
        let n_a = boundary.len();
        let n_tiles = self.tiles.len();
        if n_a == 0 || n_tiles == 0 {
            return vec![];
        }
        let mut state: Vec<TileState> =
            (0..n_tiles).map(|j| self.init_tile_state(j, n_a)).collect();

        for k in 0..(2 * n_a) {
            let sym_idx = self.sym_idx_of(boundary[k % n_a]);
            for (tile_j, ts) in state.iter_mut().enumerate() {
                if ts.skip {
                    continue;
                }
                self.step_tile(ts, tile_j, k, sym_idx);
            }
        }

        let mut out: Vec<TileMatch> = Vec::new();
        for (tile_j, ts) in state.into_iter().enumerate() {
            out.extend(self.emissions_to_matches(ts.emissions, boundary, tile_j));
        }
        out
    }

    /// Stream an arbitrary cyclic boundary against tile `j` and return
    /// all maximal cyclic RC matches (length ≥ 1). The returned
    /// `TileMatch::tile_a` is set to 0 — callers that need a specific
    /// tile_a should set it after the call.
    pub fn stream_boundary(&self, boundary: &[i8], tile_j: usize) -> Vec<TileMatch> {
        let n_a = boundary.len();
        if n_a == 0 || self.tiles[tile_j].is_empty() {
            return vec![];
        }
        let mut ts = self.init_tile_state(tile_j, n_a);

        for k in 0..(2 * n_a) {
            let sym_idx = self.sym_idx_of(boundary[k % n_a]);
            self.step_tile(&mut ts, tile_j, k, sym_idx);
        }

        self.emissions_to_matches(ts.emissions, boundary, tile_j)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn naive_cyclic_rc_matches(a: &[i8], b: &[i8]) -> Vec<(usize, usize, usize)> {
        let n_a = a.len();
        let n_b = b.len();
        if n_a == 0 || n_b == 0 {
            return vec![];
        }
        let max_len = n_a.min(n_b);
        let mut result = Vec::new();

        for ia in 0..n_a {
            for ib_rc in 0..n_b {
                let pos_b = n_b - 1 - ib_rc;
                let mut l = 0;
                while l < max_len {
                    let pa = (ia + l) % n_a;
                    let prcb = (ib_rc + l) % n_b;
                    let rc_val = -b[n_b - 1 - prcb];
                    if a[pa] != rc_val {
                        break;
                    }
                    l += 1;
                }
                if l == 0 {
                    continue;
                }
                let left_maximal = if l >= n_a || l >= n_b {
                    true
                } else {
                    let prev_a = a[(ia + n_a - 1) % n_a];
                    let prev_b_rc = -b[(pos_b + 1) % n_b];
                    prev_a != prev_b_rc
                };
                if !left_maximal {
                    continue;
                }
                result.push((ia, pos_b, l));
            }
        }

        result.sort();
        result.dedup();
        result
    }

    fn check_against_naive(strings: &[Vec<i8>], i: usize, j: usize) {
        let idx = BitParallelMatcher::new(strings);
        let matches = idx.maximal_rc_matches(i, j);
        let expected = naive_cyclic_rc_matches(&strings[i], &strings[j]);

        let mut got: Vec<(usize, usize, usize)> =
            matches.iter().map(|m| (m.a.range.start_offset, m.b.range.start_offset, m.len())).collect();
        got.sort();
        got.dedup();

        assert_eq!(
            got, expected,
            "BitParallel vs naive: tile {i} vs tile {j}\n  got={got:?}\n  exp={expected:?}"
        );
    }

    fn check(strings: &[Vec<i8>], i: usize, j: usize) {
        check_against_naive(strings, i, j);
    }

    /// Independently verify that an emitted `(pos_a, pos_b, len)`
    /// triple actually describes a valid reverse-complementary
    /// substring match in the originals: `a[pos_a..pos_a+len]` walked
    /// forward must equal `b[pos_b..pos_b-len+1]` walked backward,
    /// each character negated.
    fn verify_rc_content(a: &[i8], b: &[i8], pos_a: usize, pos_b: usize, len: usize) {
        let n_a = a.len();
        let n_b = b.len();
        assert!(len > 0, "match length must be positive");
        assert!(len <= n_a && len <= n_b, "match length exceeds tile size");
        for k in 0..len {
            let pa = (pos_a + k) % n_a;
            let pb = (n_b + pos_b - k) % n_b;
            assert_eq!(
                a[pa], -b[pb],
                "RC content mismatch at offset {k}: a[{pa}]={} vs -b[{pb}]={}",
                a[pa], -b[pb],
            );
        }
    }

    /// Structural verifier: for every reported match between tile
    /// `i` and tile `j`, sanity-check the indexing and the actual
    /// substring content. Catches regressions where the emitted
    /// triples drift out of sync with the underlying sequences.
    fn verify_all_matches(strings: &[Vec<i8>], i: usize, j: usize) {
        let bp = BitParallelMatcher::new(strings);
        let matches = bp.maximal_rc_matches(i, j);
        for m in &matches {
            assert_eq!(m.a.tile_id, i);
            assert_eq!(m.b.tile_id, j);
            assert!(m.a.range.start_offset < strings[i].len());
            assert!(m.b.range.start_offset < strings[j].len());
            assert!(!m.is_empty());
            verify_rc_content(&strings[i], &strings[j], m.a.range.start_offset, m.b.range.start_offset, m.len());
        }
    }

    #[test]
    fn identical_tiles_full_match() {
        let s = vec![1i8, 2, 3];
        check(&[s.clone(), s.clone()], 0, 1);
    }

    #[test]
    fn revcomp_tile_full_match() {
        let a = vec![1i8, 2, 3];
        let b: Vec<i8> = a.iter().rev().map(|&x| -x).collect();
        check(&[a, b], 0, 1);
    }

    #[test]
    fn no_match() {
        let a = vec![1i8, 2, 3];
        let b = vec![4i8, 5, 6];
        check(&[a, b], 0, 1);
    }

    #[test]
    fn partial_rc_match() {
        let a = vec![1i8, 2, 3, 4];
        let b = vec![5i8, 6, -3, -4];
        check(&[a, b], 0, 1);
    }

    #[test]
    fn cyclic_wraparound() {
        let a = vec![1i8, 2, 3];
        let b = vec![-1i8, 5, -3];
        check(&[a, b], 0, 1);
    }

    #[test]
    fn self_match() {
        let a = vec![1i8, 2, 1];
        check(&[a], 0, 0);
    }

    #[test]
    fn square_tiles() {
        let sq = vec![2i8, 2, 2, -6];
        check(&[sq.clone(), sq.clone()], 0, 1);
    }

    #[test]
    fn spectre_like_angles() {
        let a: Vec<i8> = vec![3, 2, 0, 2, -3, 2, 3, 2, -3, 2, 3, -2, 3, -2];
        check(&[a.clone(), a.clone()], 0, 1);
    }

    #[test]
    fn empty_tile() {
        // BP must construct successfully when one of the tiles is
        // empty, and return no matches when either side is empty.
        check(&[vec![], vec![1i8, 2]], 0, 1);
        check(&[vec![1i8, 2], vec![]], 0, 1);
        check(&[vec![], vec![1i8, 2]], 1, 0);
    }

    #[test]
    fn e2e_spectre_self_match_max_len_is_three() {
        // The spectre's longest cyclic RC self-match is known to be
        // exactly 3 edges — guard against regressions that would
        // change this invariant.
        let angles: Vec<i8> = vec![3, 2, 0, 2, -3, 2, 3, 2, -3, 2, 3, -2, 3, -2];
        let tiles = vec![angles.clone(), angles];
        check(&tiles, 0, 1);
        verify_all_matches(&tiles, 0, 1);

        let bp = BitParallelMatcher::new(&tiles);
        let matches = bp.maximal_rc_matches(0, 1);
        assert!(!matches.is_empty(), "spectre should have self-matches");
        let max_len = matches.iter().map(|m| m.len()).max().unwrap();
        assert_eq!(max_len, 3, "spectre max RC self-match should be 3");
    }

    #[test]
    fn e2e_spectre_pairwise_rotation() {
        // Spectre vs a rotated copy of itself — different starting
        // angle, same cyclic content. Exercises the cross-tile path
        // with non-identical (but related) inputs.
        let a: Vec<i8> = vec![3, 2, 0, 2, -3, 2, 3, 2, -3, 2, 3, -2, 3, -2];
        let b: Vec<i8> = vec![0, 2, -3, 2, 3, 2, -3, 2, 3, -2, 3, -2, 3, 2];
        let tiles = vec![a, b];
        for i in 0..2 {
            for j in 0..2 {
                check(&tiles, i, j);
                verify_all_matches(&tiles, i, j);
            }
        }
    }

    #[test]
    fn structural_verifier_on_existing_fixtures() {
        // Run the independent RC-content verifier across the
        // larger fixtures, in addition to the naive-oracle and
        // CMI-equivalence checks done by `check`.
        let penrose_wide: Vec<i8> = vec![
            2, 0, -1, 2, -1, 0, 0, 3, 0, 0, 1, -2, 1, 0, 0, 2, 0, 0, -1, 2, -1, 0, 0, 3, 0, 1, -2,
            1, 0,
        ];
        let penrose_narrow: Vec<i8> = vec![
            1, 0, -1, 2, -1, 0, 0, 4, 0, 0, -1, 2, -1, 0, 0, 1, 0, 1, -2, 1, 0, 0, 4, 0, 0, 1, -2,
            1, 0,
        ];
        let tiles = vec![penrose_wide, penrose_narrow];
        for i in 0..2 {
            for j in 0..2 {
                verify_all_matches(&tiles, i, j);
            }
        }

        let s_tet: Vec<i8> = vec![-1, 1, 1, 0, 1, -1, 1, 1, 0, 1];
        verify_all_matches(&[s_tet], 0, 0);
    }

    #[test]
    fn single_angle_tiles() {
        check(&[vec![2i8], vec![-2i8]], 0, 1);
        check(&[vec![2i8], vec![2i8]], 0, 1);
    }

    #[test]
    fn three_tiles_pairwise() {
        let tiles: Vec<Vec<i8>> = vec![vec![1i8, 2, 3], vec![2i8, 3, 4], vec![-3i8, -2, 5]];
        for i in 0..tiles.len() {
            for j in 0..tiles.len() {
                check(&tiles, i, j);
            }
        }
    }

    #[test]
    fn repeated_angles() {
        let a = vec![1i8, 1, 1];
        check(&[a.clone(), a.clone()], 0, 1);
    }

    #[test]
    fn negative_angles() {
        let a = vec![-1i8, -2, -3];
        let b = vec![3i8, 2, 1];
        check(&[a, b], 0, 1);
    }

    #[test]
    fn hexagons_no_rc_match() {
        let hex: Vec<i8> = vec![1, 1, 1, 1, 1, 1];
        check(&[hex.clone(), hex.clone()], 0, 1);
    }

    #[test]
    fn hexagon_vs_reversed_hexagon() {
        let hex: Vec<i8> = vec![1, 1, 1, 1, 1, 1];
        let hex_rev: Vec<i8> = vec![-1, -1, -1, -1, -1, -1];
        check(&[hex, hex_rev], 0, 1);
    }

    #[test]
    fn penrose_rhombs() {
        let wide: Vec<i8> = vec![
            2, 0, -1, 2, -1, 0, 0, 3, 0, 0, 1, -2, 1, 0, 0, 2, 0, 0, -1, 2, -1, 0, 0, 3, 0, 1, -2,
            1, 0,
        ];
        let narrow: Vec<i8> = vec![
            1, 0, -1, 2, -1, 0, 0, 4, 0, 0, -1, 2, -1, 0, 0, 1, 0, 1, -2, 1, 0, 0, 4, 0, 0, 1, -2,
            1, 0,
        ];
        let tiles = vec![wide, narrow];
        for i in 0..2 {
            for j in 0..2 {
                check(&tiles, i, j);
            }
        }
    }

    #[test]
    fn mixed_shapes() {
        let square: Vec<i8> = vec![1, 1, 1, 1];
        let triangle: Vec<i8> = vec![2, 2, 2];
        let rect: Vec<i8> = vec![0, 0, 1, 1, 0, 0, 1, 1];
        let tiles = vec![square, triangle, rect];
        for i in 0..tiles.len() {
            for j in 0..tiles.len() {
                check(&tiles, i, j);
            }
        }
    }

    #[test]
    fn tetrominos() {
        let t_o: Vec<i8> = vec![0, 1, 0, 1, 0, 1, 0, 1];
        let t_i: Vec<i8> = vec![0, 0, 0, 1, 1, 0, 0, 0, 1, 1];
        let t_t: Vec<i8> = vec![-1, 1, 1, -1, 1, 1, 0, 0, 1, 1];
        let tiles = vec![t_o, t_i, t_t];
        for i in 0..tiles.len() {
            for j in 0..tiles.len() {
                check(&tiles, i, j);
            }
        }
    }

    #[test]
    fn s_tetromino_self_match() {
        let t_s: Vec<i8> = vec![-1, 1, 1, 0, 1, -1, 1, 1, 0, 1];
        check(&[t_s], 0, 0);
    }

    #[test]
    fn random_short_pairwise() {
        let mut seed: u64 = 54321;
        for _ in 0..40 {
            let tiles: Vec<Vec<i8>> = (0..4)
                .map(|_| {
                    (0..6)
                        .map(|_| {
                            seed = seed
                                .wrapping_mul(6364136223846793005)
                                .wrapping_add(1442695040888963407);
                            ((seed % 7) as i8) - 3
                        })
                        .collect()
                })
                .collect();
            for i in 0..tiles.len() {
                for j in 0..tiles.len() {
                    check(&tiles, i, j);
                }
            }
        }
    }

    #[test]
    fn random_longer_pairwise() {
        let mut seed: u64 = 99999;
        let mut next_angle = || {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((seed % 7) as i8) - 3
        };
        let tiles: Vec<Vec<i8>> = (0..8)
            .map(|_| (0..10).map(|_| next_angle()).collect())
            .collect();
        for i in 0..tiles.len() {
            for j in 0..tiles.len() {
                check(&tiles, i, j);
            }
        }
    }

    #[test]
    fn long_tile_crosses_word_boundary() {
        // Force tile lengths past 64 bits to exercise multi-word
        // bitset arithmetic and the cyclic wrap-around path.
        let mut seed: u64 = 0xdeadbeefcafebabe;
        let mut next_angle = || {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((seed % 9) as i8) - 4
        };
        let tile_a: Vec<i8> = (0..73).map(|_| next_angle()).collect();
        let tile_b: Vec<i8> = (0..81).map(|_| next_angle()).collect();
        let tiles = vec![tile_a, tile_b];
        for i in 0..2 {
            for j in 0..2 {
                check(&tiles, i, j);
            }
        }
    }

    #[test]
    fn sweep_matches_per_pair_across_fixtures() {
        // For every fixture we use elsewhere, assert that one sweep of
        // `stream_boundary_all_tiles` yields the same TileMatch set
        // (modulo tile_a labelling) as calling `stream_boundary` per
        // tile and unioning.
        type Fixture = (&'static str, Vec<Vec<i8>>, Vec<Vec<i8>>);
        let fixtures: Vec<Fixture> = vec![
            (
                "spectre tile against itself",
                vec![vec![3i8, 2, 0, 2, -3, 2, 3, 2, -3, 2, 3, -2, 3, -2]],
                vec![vec![3i8, 2, 0, 2, -3, 2, 3, 2, -3, 2, 3, -2, 3, -2]],
            ),
            (
                "small mixed boundary vs three tiles",
                vec![vec![1i8, 2, 3, -2, 0, 1, -2]],
                vec![
                    vec![1i8, 2, 3, -2, 0, 1, -2],
                    vec![0i8, 2, -3, 2, 3, -2, 3, -2],
                    vec![1i8, 1, 1, -3, 1, 1, 1, -3],
                ],
            ),
            (
                "external boundary with -1",
                vec![vec![3i8, 2, -3, 1, -2, 0, 2, -1]],
                vec![
                    vec![1i8, 2, 3, -2, 0, 1, -2],
                    vec![0i8, 2, -3, 2, 3, -2, 3, -2],
                ],
            ),
            (
                "penrose pair",
                vec![
                    vec![
                        2i8, 0, -1, 2, -1, 0, 0, 3, 0, 0, 1, -2, 1, 0, 0, 2, 0, 0, -1, 2, -1, 0, 0,
                        3, 0, 1, -2, 1, 0,
                    ],
                    vec![
                        1i8, 0, -1, 2, -1, 0, 0, 4, 0, 0, -1, 2, -1, 0, 0, 1, 0, 1, -2, 1, 0, 0, 4,
                        0, 0, 1, -2, 1, 0,
                    ],
                ],
                vec![
                    vec![
                        2i8, 0, -1, 2, -1, 0, 0, 3, 0, 0, 1, -2, 1, 0, 0, 2, 0, 0, -1, 2, -1, 0, 0,
                        3, 0, 1, -2, 1, 0,
                    ],
                    vec![
                        1i8, 0, -1, 2, -1, 0, 0, 4, 0, 0, -1, 2, -1, 0, 0, 1, 0, 1, -2, 1, 0, 0, 4,
                        0, 0, 1, -2, 1, 0,
                    ],
                ],
            ),
        ];

        for (label, boundaries, tiles) in fixtures {
            let bp = BitParallelMatcher::new(&tiles);
            for boundary in &boundaries {
                let sweep = bp.stream_boundary_all_tiles(boundary);
                let mut per_pair: Vec<TileMatch> = (0..tiles.len())
                    .flat_map(|j| bp.stream_boundary(boundary, j))
                    .collect();

                let mut sweep_tuples: Vec<(usize, usize, usize, usize)> = sweep
                    .iter()
                    .map(|m| (m.b.tile_id, m.a.range.start_offset, m.b.range.start_offset, m.len()))
                    .collect();
                let mut pair_tuples: Vec<(usize, usize, usize, usize)> = per_pair
                    .iter_mut()
                    .map(|m| (m.b.tile_id, m.a.range.start_offset, m.b.range.start_offset, m.len()))
                    .collect();
                sweep_tuples.sort();
                pair_tuples.sort();
                assert_eq!(
                    sweep_tuples, pair_tuples,
                    "{label}: sweep vs per-pair disagrees"
                );
            }
        }
    }

    #[test]
    fn stream_boundary_against_fixed_tileset() {
        // Boundary is not one of the indexed tiles. Verify
        // `stream_boundary` produces the same per-tile results as
        // building a fresh matcher with the boundary as an extra tile
        // and using `maximal_rc_matches(0, j)`.
        let tiles: Vec<Vec<i8>> = vec![
            vec![1, 2, 3, -2, 0, 1, -2],
            vec![0, 2, -3, 2, 3, -2, 3, -2],
            vec![1, 1, 1, -3, 1, 1, 1, -3],
        ];
        let boundary: Vec<i8> = vec![3, 2, -3, 1, -2, 0, 2, -1];

        let bp_fixed = BitParallelMatcher::new(&tiles);
        for j in 0..tiles.len() {
            let from_stream = bp_fixed.stream_boundary(&boundary, j);
            let mut with_boundary = tiles.clone();
            with_boundary.insert(0, boundary.clone());
            let bp_full = BitParallelMatcher::new(&with_boundary);
            let from_pair = bp_full.maximal_rc_matches(0, j + 1);

            let mut a: Vec<(usize, usize, usize)> = from_stream
                .iter()
                .map(|m| (m.a.range.start_offset, m.b.range.start_offset, m.len()))
                .collect();
            let mut b: Vec<(usize, usize, usize)> = from_pair
                .iter()
                .map(|m| (m.a.range.start_offset, m.b.range.start_offset, m.len()))
                .collect();
            a.sort();
            a.dedup();
            b.sort();
            b.dedup();
            assert_eq!(a, b, "stream_boundary vs paired build, j={j}");
        }
    }
}
