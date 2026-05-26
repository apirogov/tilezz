//! Match enumeration between tile boundaries.
//!
//! Two related abstractions live here:
//!
//! - [`MatchFinder`] enumerates legal glue matches between two
//!   tilesets (or a tileset against itself) using cyclic string-match
//!   indices. Use this for one-shot or filtered enumeration when you
//!   need flexibility over which tile pairs to scan and what filter
//!   to apply.
//! - [`MatchTypeIndex`] pre-computes the full match enumeration for a
//!   single tileset and indexes the results by starting `(tile, offset)`.
//!   Use this when you'll query many times against a fixed tileset and
//!   want O(1) candidate lookup.
//!
//! Both share the [`TileMatch`] record describing a single match.

use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

use rustc_hash::FxHashSet;

use crate::cyclotomic::{IsRing};
use crate::geom::glue::junctions_glueable;
use crate::geom::rat::Rat;
use crate::geom::snake::Snake;
use crate::geom::tileset::TileSet;
use crate::geom::matches::{EdgeRange, Segment, TileMatch};
use crate::stringmatch::BitParallelMatcher;

/// Bit-parallel engine state. The B-side masks live in a shared
/// `BitParallelMatcher`; the A-side angle sequences are kept raw and
/// streamed against `b_matcher` per query. The shared masks let the
/// seed-side precompute be built once and shared via `Arc` across
/// many short-lived `MatchFinder`s with varying A-side tilesets —
/// see [`BpSeed`] and [`MatchFinder::crossing_bp_with_seed`].
struct CyclicEngine {
    a_sequences: Vec<Vec<i8>>,
    b_matcher: Arc<BitParallelMatcher>,
    offset_b: usize,
}

impl CyclicEngine {
    fn maximal_rc_matches(&self, i: usize, j: usize) -> Vec<TileMatch> {
        let b_idx = j - self.offset_b;
        let mut matches = self.b_matcher.stream_boundary(&self.a_sequences[i], b_idx);
        for m in &mut matches {
            m.a.tile_id = i;
        }
        matches
    }

    fn maximal_rc_matches_at_positions(
        &self,
        i: usize,
        j: usize,
        positions: &[usize],
    ) -> Vec<TileMatch> {
        if positions.is_empty() {
            return vec![];
        }
        let b_idx = j - self.offset_b;
        let boundary = &self.a_sequences[i];
        let n_a = boundary.len();
        let mut keep = vec![false; n_a];
        for &p in positions {
            if p < n_a {
                keep[p] = true;
            }
        }
        let mut matches = self.b_matcher.stream_boundary(boundary, b_idx);
        matches.retain(|m| keep[m.a.range.start_offset]);
        for m in &mut matches {
            m.a.tile_id = i;
        }
        matches
    }
}

/// Pre-built bit-parallel state for a fixed "seed" (B-side) tileset.
///
/// Construct once, then pass to [`MatchFinder::crossing_with_seed`]
/// for each varying A-side tileset (the typical pattern in
/// `seq_collect` / `seq_explorer` / `patch_enum` where many patch
/// chunks are matched against a fixed tile alphabet). Cheap to clone
/// (`Arc` internals).
pub struct BpSeed<T: IsRing> {
    tileset: Arc<TileSet<T>>,
    matcher: Arc<BitParallelMatcher>,
}

impl<T: IsRing> Clone for BpSeed<T> {
    fn clone(&self) -> Self {
        BpSeed {
            tileset: Arc::clone(&self.tileset),
            matcher: Arc::clone(&self.matcher),
        }
    }
}

impl<T: IsRing> BpSeed<T> {
    /// Pre-compute the bit-parallel state for the given seed tileset.
    pub fn new(tileset: Arc<TileSet<T>>) -> Self {
        let sequences: Vec<Vec<i8>> = tileset.rats().iter().map(|r| r.seq().to_vec()).collect();
        let matcher = Arc::new(BitParallelMatcher::new(&sequences));
        BpSeed { tileset, matcher }
    }

    pub fn tileset(&self) -> &Arc<TileSet<T>> {
        &self.tileset
    }
}

/// Realise a `TileMatch` against the two tilesets, producing the
/// resulting boundary as a `Rat`. Assumes the match has been previously
/// validated (e.g. it came from a `MatchFinder`); will panic otherwise.
pub fn apply_match<T: IsRing>(
    m: &TileMatch,
    rats_a: &[Rat<T>],
    rats_b: &[Rat<T>],
) -> Rat<T> {
    rats_a[m.a.tile_id]
        .try_glue_precomputed(
            (
                m.a.range.start_offset as i64,
                m.len(),
                m.b.range.start_offset as i64,
            ),
            &rats_b[m.b.tile_id],
            true,
        )
        .expect("match was pre-validated")
}

/// The involution of a match: same glue described from the other tile's
/// perspective. Adjusts the asymmetric `start_b` (= first surviving)
/// convention by shifting starts to land on the right edges of the
/// inverted match.
fn match_involution(m: &TileMatch, len_a: usize, len_b: usize) -> TileMatch {
    let len = m.len();
    TileMatch::new(
        Segment::new(
            m.b.tile_id,
            EdgeRange::new((m.b.range.start_offset + len_b - len) % len_b, len),
        ),
        Segment::new(
            m.a.tile_id,
            EdgeRange::new((m.a.range.start_offset + len) % len_a, len),
        ),
    )
}

struct TypeEntry {
    first: TileMatch,
    second: TileMatch,
}

/// Enumerates legal glue matches between two tilesets (or one tileset
/// against itself).
///
/// Use the constructor that matches your scenario:
///
/// * [`MatchFinder::new`] — match a tileset against **itself**. Both
///   `set_a` and `set_b` point at the same `Arc<TileSet>`. The
///   resulting `TileMatch`s describe glues between two tiles drawn
///   from the same tileset.
/// * [`MatchFinder::crossing`] — match **two distinct** tilesets
///   against each other. `set_a` and `set_b` are different. Useful
///   when one tileset represents already-placed boundaries and the
///   other represents candidates to glue.
///
/// Internally backed by a [`BitParallelMatcher`] over the B-side
/// tileset; the A-side angle sequences are streamed against it per
/// query. Build cost is proportional to the B-side only; consider
/// [`Self::crossing_with_seed`] when the B-side is fixed across many
/// `MatchFinder` builds.
pub struct MatchFinder<T: IsRing> {
    set_a: Arc<TileSet<T>>,
    set_b: Arc<TileSet<T>>,
    offset_b: usize,
    engine: CyclicEngine,
}

impl<T: IsRing> MatchFinder<T> {
    /// Build a `MatchFinder` that enumerates self-matches within a
    /// single tileset (`set_a == set_b == tileset`).
    pub fn new(tileset: Arc<TileSet<T>>) -> Self {
        let sequences: Vec<Vec<i8>> = tileset.rats().iter().map(|r| r.seq().to_vec()).collect();
        let b_matcher = Arc::new(BitParallelMatcher::new(&sequences));
        let engine = CyclicEngine {
            a_sequences: sequences,
            b_matcher,
            offset_b: 0,
        };
        MatchFinder {
            set_a: Arc::clone(&tileset),
            set_b: tileset,
            offset_b: 0,
            engine,
        }
    }

    /// Build a `MatchFinder` that enumerates matches between two
    /// distinct tilesets `a` and `b`. The returned `TileMatch` values
    /// have `tile_a` indexed into `a.rats()` and `tile_b` indexed into
    /// `b.rats()`. The B-side BP state is built fresh here — use
    /// [`Self::crossing_with_seed`] to amortise across many builds
    /// with the same B-side.
    pub fn crossing(a: Arc<TileSet<T>>, b: Arc<TileSet<T>>) -> Self {
        let seed = BpSeed::new(b);
        Self::crossing_with_seed(a, seed)
    }

    /// Build a crossing `MatchFinder` reusing a pre-built [`BpSeed`]
    /// for the B-side. Per-call cost is proportional to the A-side
    /// tileset only (cloning the angle sequences); the B-side masks
    /// are shared by `Arc`.
    pub fn crossing_with_seed(a: Arc<TileSet<T>>, seed: BpSeed<T>) -> Self {
        let a_sequences: Vec<Vec<i8>> = a.rats().iter().map(|r| r.seq().to_vec()).collect();
        let offset_b = a.num_tiles();
        let engine = CyclicEngine {
            a_sequences,
            b_matcher: seed.matcher,
            offset_b,
        };
        MatchFinder {
            set_a: a,
            set_b: seed.tileset,
            offset_b,
            engine,
        }
    }

    /// The A-side tileset.
    pub fn set_a(&self) -> &Arc<TileSet<T>> {
        &self.set_a
    }

    /// The B-side tileset (equal to `set_a` when constructed via [`Self::new`]).
    pub fn set_b(&self) -> &Arc<TileSet<T>> {
        &self.set_b
    }

    /// Tile `i` from the A-side tileset.
    pub fn rat_a(&self, i: usize) -> &Rat<T> {
        self.set_a.rat(i)
    }

    /// Tile `j` from the B-side tileset.
    pub fn rat_b(&self, j: usize) -> &Rat<T> {
        self.set_b.rat(j)
    }

    /// Number of tiles in the A-side tileset.
    pub fn num_tiles_a(&self) -> usize {
        self.set_a.num_tiles()
    }

    /// Number of tiles in the B-side tileset.
    pub fn num_tiles_b(&self) -> usize {
        self.set_b.num_tiles()
    }

    /// Realise a `TileMatch` against this finder's tilesets, producing
    /// the resulting boundary `Rat`.
    pub fn apply_match(&self, m: &TileMatch) -> Rat<T> {
        apply_match(m, self.set_a.rats(), self.set_b.rats())
    }

    /// Return the maximal reverse-complementary cyclic substrings
    /// shared between tile `i` (A-side) and tile `j` (B-side).
    /// Includes raw matches that haven't yet been validated as legal
    /// glues (e.g. they might cause self-intersection).
    pub fn shared_boundaries(&self, i: usize, j: usize) -> Vec<TileMatch> {
        self.engine.maximal_rc_matches(i, self.offset_b + j)
    }

    /// Enumerate the legal glue matches between tile `i` (A-side) and
    /// tile `j` (B-side). Matches are validated (no self-intersection,
    /// no degenerate junction angles) and sorted by `(start_a, start_b)`.
    pub fn valid_matches(&self, i: usize, j: usize) -> Vec<TileMatch> {
        let mut m = Self::validated_matches(self.candidates_for_pair(i, j));
        Self::sort_by_interval(&mut m);
        m
    }

    /// Same as [`Self::valid_matches`] but returns the realised
    /// boundary `Rat` alongside each match.
    pub fn valid_matches_with_rats(&self, i: usize, j: usize) -> Vec<(Rat<T>, TileMatch)> {
        let mut m = Self::validated_matches_with_rats(self.candidates_for_pair(i, j));
        m.sort_by(|a, b| {
            a.1.a
                .range
                .start_offset
                .cmp(&b.1.a.range.start_offset)
                .then_with(|| a.1.b.range.start_offset.cmp(&b.1.b.range.start_offset))
        });
        m
    }
}

impl<T: IsRing> MatchFinder<T> {
    /// Like [`Self::valid_matches_with_rats`] but only considers
    /// matches that touch at least one "active" position on the A-side
    /// boundary. `active[k]` selects edge `k` of tile `i`. Useful when
    /// you want to enumerate glues incident with a specific boundary
    /// region (e.g. a partial growth step).
    pub fn valid_matches_filtered(
        &self,
        i: usize,
        j: usize,
        active: &[bool],
    ) -> Vec<(Rat<T>, TileMatch)> {
        let a = self.set_a.rat(i);
        let b = self.set_b.rat(j);
        let cmi_j = self.offset_b + j;
        let n_a = a.len();
        let n_b = b.len();

        if n_a == 0 || n_b == 0 || active.is_empty() {
            return Vec::new();
        }

        let max_match_len = n_a.min(n_b);
        let mut near_active = vec![false; n_a];
        for pos in 0..n_a {
            if active[pos] {
                for d in 0..=max_match_len {
                    near_active[(pos + n_a - d) % n_a] = true;
                    near_active[(pos + d) % n_a] = true;
                }
            }
        }
        let scan_positions: Vec<usize> = (0..n_a).filter(|&p| near_active[p]).collect();

        let mut raw: Vec<(Rat<T>, TileMatch)> = Vec::new();

        let cmi_matches = self
            .engine
            .maximal_rc_matches_at_positions(i, cmi_j, &scan_positions);
        for m in &cmi_matches {
            let (ns, len, ne) = a.get_match((m.a.range.start_offset as i64, m.b.range.start_offset as i64), b);
            if len <= 1 {
                continue;
            }
            if !junctions_glueable(a.seq(), ns as usize, len, b.seq(), ne as usize) {
                continue;
            }
            let ns_u = ns.rem_euclid(n_a as i64) as usize;
            if !match_touches_active(ns_u, len, n_a, active) {
                continue;
            }
            if let Ok(glued) = a.try_glue_precomputed((ns, len, ne), b, true) {
                raw.push((
                    glued,
                    TileMatch::new(
                        Segment::new(i, EdgeRange::new(ns_u, len)),
                        Segment::new(j, EdgeRange::new(ne.rem_euclid(n_b as i64) as usize, len)),
                    ),
                ));
            }
        }

        let seq_a = a.seq();
        let seq_b = b.seq();
        for ia in 0..n_a {
            for ib in 0..n_b {
                if !junctions_glueable(seq_a, ia, 1, seq_b, ib) {
                    continue;
                }
                let (ns, len, ne) = a.get_match((ia as i64, ib as i64), b);
                if len != 1 {
                    continue;
                }
                let ns_u = ns.rem_euclid(n_a as i64) as usize;
                if !match_touches_active(ns_u, len, n_a, active) {
                    continue;
                }
                if let Ok(glued) = a.try_glue_precomputed((ns, len, ne), b, true) {
                    raw.push((
                        glued,
                        TileMatch::new(
                            Segment::new(i, EdgeRange::new(ns_u, len)),
                            Segment::new(j, EdgeRange::new(ne.rem_euclid(n_b as i64) as usize, len)),
                        ),
                    ));
                }
            }
        }

        let mut seen: FxHashSet<Rat<T>> = FxHashSet::default();
        raw.into_iter()
            .filter(|(rat, _)| seen.insert(rat.clone()))
            .collect()
    }

    /// Enumerate all legal glue matches across every `(A-tile, B-tile)`
    /// pair. Sorted into a deterministic global order. Result excludes
    /// matches that are involutions of one another (each glue is listed
    /// once, not twice).
    pub fn all_valid_matches(&self) -> Vec<TileMatch> {
        self.valid_matches_for_pairs(&self.all_pairs())
    }

    /// Like [`Self::all_valid_matches`] but restricted to the given
    /// tile-pair set. `(i, j)` denotes `(A-tile i, B-tile j)`.
    pub fn valid_matches_for_pairs(&self, pairs: &[(usize, usize)]) -> Vec<TileMatch> {
        let mut m = Self::validated_matches(self.collect_candidates(pairs));
        Self::sort_global(&mut m);
        m
    }

    /// Enumerate the distinct *resulting boundaries* (as `Rat`s)
    /// produced by the legal glues across the given tile pairs.
    /// Deduped via the `BTreeSet`.
    pub fn valid_results_for_pairs(&self, pairs: &[(usize, usize)]) -> BTreeSet<Rat<T>> {
        Self::validated_results(self.collect_candidates(pairs))
    }

    fn candidates_for_pair(&self, i: usize, j: usize) -> BTreeMap<Rat<T>, Vec<TileMatch>> {
        let a = self.set_a.rat(i);
        let b = self.set_b.rat(j);
        let cmi_j = self.offset_b + j;
        let n_a = a.len();
        let n_b = b.len();

        let mut groups: BTreeMap<Rat<T>, Vec<TileMatch>> = BTreeMap::new();
        if n_a == 0 || n_b == 0 {
            return groups;
        }

        let cmi_matches = self.engine.maximal_rc_matches(i, cmi_j);
        for m in &cmi_matches {
            let (ns, len, ne) = a.get_match((m.a.range.start_offset as i64, m.b.range.start_offset as i64), b);
            if len <= 1 {
                continue;
            }
            if !junctions_glueable(a.seq(), ns as usize, len, b.seq(), ne as usize) {
                continue;
            }
            if let Ok(glued) = a.try_glue_precomputed((ns, len, ne), b, true) {
                groups.entry(glued).or_default().push(TileMatch::new(
                    Segment::new(i, EdgeRange::new(ns as usize, len)),
                    Segment::new(j, EdgeRange::new(ne as usize, len)),
                ));
            }
        }

        let seq_a = a.seq();
        let seq_b = b.seq();
        for ia in 0..n_a {
            for ib in 0..n_b {
                if !junctions_glueable(seq_a, ia, 1, seq_b, ib) {
                    continue;
                }
                let (ns, len, ne) = a.get_match((ia as i64, ib as i64), b);
                if len != 1 {
                    continue;
                }
                if let Ok(glued) = a.try_glue_precomputed((ns, len, ne), b, true) {
                    groups.entry(glued).or_default().push(TileMatch::new(
                        Segment::new(i, EdgeRange::new(ns as usize, len)),
                        Segment::new(j, EdgeRange::new(ne as usize, len)),
                    ));
                }
            }
        }

        groups
    }

    fn all_pairs(&self) -> Vec<(usize, usize)> {
        let na = self.num_tiles_a();
        let nb = self.num_tiles_b();
        (0..na).flat_map(|i| (0..nb).map(move |j| (i, j))).collect()
    }

    fn collect_candidates(&self, pairs: &[(usize, usize)]) -> BTreeMap<Rat<T>, Vec<TileMatch>> {
        let mut all: BTreeMap<Rat<T>, Vec<TileMatch>> = BTreeMap::new();
        for &(i, j) in pairs {
            for (rat, matches) in self.candidates_for_pair(i, j) {
                all.entry(rat).or_default().extend(matches);
            }
        }
        all
    }

    fn validated_matches(groups: BTreeMap<Rat<T>, Vec<TileMatch>>) -> Vec<TileMatch> {
        groups
            .into_iter()
            .filter(|(rat, _)| Snake::<T>::try_from(rat.seq()).is_ok())
            .flat_map(|(_, matches)| matches)
            .collect()
    }

    fn validated_matches_with_rats(
        groups: BTreeMap<Rat<T>, Vec<TileMatch>>,
    ) -> Vec<(Rat<T>, TileMatch)> {
        groups
            .into_iter()
            .filter(|(rat, _)| Snake::<T>::try_from(rat.seq()).is_ok())
            .flat_map(|(rat, matches)| matches.into_iter().map(move |m| (rat.clone(), m)))
            .collect()
    }

    fn validated_results(groups: BTreeMap<Rat<T>, Vec<TileMatch>>) -> BTreeSet<Rat<T>> {
        groups
            .into_iter()
            .filter(|(rat, _)| Snake::<T>::try_from(rat.seq()).is_ok())
            .map(|(rat, _)| rat)
            .collect()
    }

    fn sort_by_interval(results: &mut [TileMatch]) {
        results.sort_by(|a, b| {
            a.a.range
                .start_offset
                .cmp(&b.a.range.start_offset)
                .then_with(|| a.b.range.start_offset.cmp(&b.b.range.start_offset))
        });
    }

    fn sort_global(results: &mut [TileMatch]) {
        results.sort_by(|a, b| {
            a.a.tile_id
                .cmp(&b.a.tile_id)
                .then_with(|| a.b.tile_id.cmp(&b.b.tile_id))
                .then_with(|| a.a.range.start_offset.cmp(&b.a.range.start_offset))
                .then_with(|| a.b.range.start_offset.cmp(&b.b.range.start_offset))
        });
    }
}

/// Pre-computed index of all legal glue matches for a single tileset,
/// indexed by starting `(tile, offset)`.
///
/// Built by exhaustively enumerating self-matches via [`MatchFinder`]
/// at construction time, then indexed for O(1) lookup. Suited for
/// repeated queries against a fixed tileset (e.g. the candidate-match
/// path used by [`GrowingPatch`](crate::geom::patch::GrowingPatch)).
///
/// Internally deduplicates each match with its involution; the public
/// [`Self::candidates_starting_at`] returns *both* directions
/// (`tile_a` and the inverted `tile_b` lookup are both indexed).
pub struct MatchTypeIndex<T: IsRing> {
    tileset: Arc<TileSet<T>>,
    entries: Vec<TypeEntry>,
    by_start: Vec<Vec<Vec<Segment>>>,
}

impl<T: IsRing> MatchTypeIndex<T> {
    /// Build the index by enumerating all legal self-matches and
    /// indexing them by `(tile_a, start_a)`.
    pub fn new(tileset: Arc<TileSet<T>>) -> Self {
        let mf = MatchFinder::new(Arc::clone(&tileset));
        let all_matches = mf.all_valid_matches();
        let rats = tileset.rats();

        let mut seen: BTreeSet<(TileMatch, TileMatch)> = BTreeSet::new();

        for mt in &all_matches {
            let len_a = rats[mt.a.tile_id].len();
            let len_b = rats[mt.b.tile_id].len();

            let fwd = *mt;
            let bwd = match_involution(mt, len_a, len_b);

            let (first, second) = if fwd <= bwd { (fwd, bwd) } else { (bwd, fwd) };

            seen.insert((first, second));
        }

        let entries: Vec<TypeEntry> = seen
            .into_iter()
            .map(|(first, second)| TypeEntry { first, second })
            .collect();

        let num_tiles = tileset.num_tiles();
        let mut by_start: Vec<Vec<Vec<Segment>>> = (0..num_tiles)
            .map(|_| vec![vec![]; rats[0].len()])
            .collect();
        for ti in 0..num_tiles {
            by_start[ti] = vec![vec![]; rats[ti].len()];
        }
        for mt in &all_matches {
            let cm = Segment::new(mt.b.tile_id, EdgeRange::new(mt.b.range.start_offset, mt.len()));
            by_start[mt.a.tile_id][mt.a.range.start_offset].push(cm);
            let inv =
                match_involution(mt, rats[mt.a.tile_id].len(), rats[mt.b.tile_id].len());
            let cm_inv =
                Segment::new(inv.b.tile_id, EdgeRange::new(inv.b.range.start_offset, inv.len()));
            by_start[inv.a.tile_id][inv.a.range.start_offset].push(cm_inv);
        }

        MatchTypeIndex {
            tileset,
            entries,
            by_start,
        }
    }

    /// The tileset this index was built from.
    pub fn tileset(&self) -> &Arc<TileSet<T>> {
        &self.tileset
    }

    /// Number of distinct match types (each pair of involutions
    /// counted once). IDs run `1..=num_types()`.
    pub fn num_types(&self) -> usize {
        self.entries.len()
    }

    /// Retrieve the canonical-direction `TileMatch` for a 1-based id.
    /// Use [`Self::signed_id`] + [`Self::apply`] if direction matters.
    pub fn get(&self, unsigned_id: usize) -> TileMatch {
        assert!(
            unsigned_id >= 1 && unsigned_id <= self.entries.len(),
            "match type id {unsigned_id} out of range [1, {}]",
            self.entries.len()
        );
        self.entries[unsigned_id - 1].first
    }

    /// Find a `TileMatch`'s id (1-based) in this index, signed by
    /// direction: positive for the canonical direction, negative for
    /// its involution. `None` if not present.
    pub fn signed_id(&self, mt: &TileMatch) -> Option<i32> {
        for (idx, entry) in self.entries.iter().enumerate() {
            if *mt == entry.first {
                return Some((idx + 1) as i32);
            }
            if *mt == entry.second {
                return Some(-((idx + 1) as i32));
            }
        }
        None
    }

    /// Realise the match at the signed id, producing the resulting
    /// boundary as a `Rat`. Positive ids use the canonical direction;
    /// negative ids use the involution.
    pub fn apply(&self, signed_id: i32) -> Rat<T> {
        let unsigned = signed_id.unsigned_abs() as usize;
        assert!(
            unsigned >= 1 && unsigned <= self.entries.len(),
            "signed match type id {signed_id} out of range"
        );
        let entry = &self.entries[unsigned - 1];
        let rats = self.tileset.rats();
        let mt = if signed_id > 0 {
            entry.first
        } else {
            entry.second
        };
        apply_match(&mt, rats, rats)
    }

    /// All candidate matches that start at edge `(tile_id, offset)` on
    /// the A side. O(1) lookup; the returned slice borrows from the
    /// pre-computed index.
    pub fn candidates_starting_at(&self, tile_id: usize, offset: usize) -> &[Segment] {
        &self.by_start[tile_id][offset]
    }

    /// All candidate matches for `tile_id`, grouped by starting
    /// offset. `result[offset]` is the same as `candidates_starting_at(tile_id, offset)`.
    pub fn all_candidates_for_tile(&self, tile_id: usize) -> &[Vec<Segment>] {
        &self.by_start[tile_id]
    }
}

/// True iff the cyclic edge range `[start_a, start_a + mlen)` of a
/// length-`n` boundary touches any active position — either an active
/// edge inside the range, or an active vertex immediately adjacent to
/// either end of the range. Used by [`MatchFinder::valid_matches_filtered`]
/// to restrict enumeration to matches incident with a caller-specified
/// region of the boundary.
fn match_touches_active(start_a: usize, mlen: usize, n: usize, active: &[bool]) -> bool {
    if mlen == 0 || n == 0 {
        return false;
    }
    let before = (start_a + n - 1) % n;
    if active[before] {
        return true;
    }
    let after = (start_a + mlen) % n;
    if active[after] {
        return true;
    }
    for i in 0..mlen {
        if active[(start_a + i) % n] {
            return true;
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic::ZZ12;
    use crate::geom::tiles::{hexagon, spectre, square};
    use std::collections::{BTreeSet, HashSet};

    const MYSTIC_ZZ12: &[i8] = &[
        0, 2, -3, 2, 3, -2, 3, -2, 3, 2, -3, 2, 0, 2, -3, 2, 3, 2, -3, 2,
    ];

    fn brute_force_valid_matches<T: IsRing>(
        a: &Rat<T>,
        b: &Rat<T>,
        i: usize,
        j: usize,
    ) -> Vec<TileMatch> {
        let mut seen: HashSet<(usize, usize, usize)> = HashSet::new();
        let mut results: Vec<TileMatch> = Vec::new();
        for ia in 0..a.len() {
            for ib in 0..b.len() {
                let (ns, len, ne) = a.get_match((ia as i64, ib as i64), b);
                if len == 0 || !seen.insert((ns as usize, len, ne as usize)) {
                    continue;
                }
                if a.try_glue_precomputed((ns, len, ne), b, false).is_ok() {
                    results.push(TileMatch::new(
                        Segment::new(i, EdgeRange::new(ns as usize, len)),
                        Segment::new(j, EdgeRange::new(ne as usize, len)),
                    ));
                }
            }
        }
        results.sort_by(|x, y| {
            x.a.range
                .start_offset
                .cmp(&y.a.range.start_offset)
                .then_with(|| x.b.range.start_offset.cmp(&y.b.range.start_offset))
        });
        results
    }

    fn assert_match_sets_match<T: IsRing>(
        mf_matches: &[TileMatch],
        bf_matches: &[TileMatch],
        rats_a: &[Rat<T>],
        rats_b: &[Rat<T>],
        label: &str,
    ) {
        assert_eq!(
            mf_matches.len(),
            bf_matches.len(),
            "{label}: count mismatch: MatchFinder={}, brute_force={}",
            mf_matches.len(),
            bf_matches.len(),
        );
        for (idx, (mf_m, bf_m)) in mf_matches.iter().zip(bf_matches.iter()).enumerate() {
            assert_eq!(
                (mf_m.a.range.start_offset, mf_m.len(), mf_m.b.range.start_offset),
                (bf_m.a.range.start_offset, bf_m.len(), bf_m.b.range.start_offset),
                "{label}: interval mismatch at index {idx}",
            );
            assert_eq!(
                apply_match(mf_m, rats_a, rats_b),
                apply_match(bf_m, rats_a, rats_b),
                "{label}: result mismatch at ({}, {})",
                mf_m.a.range.start_offset,
                mf_m.b.range.start_offset,
            );
        }
    }

    fn self_match_ts(rats: Vec<Rat<ZZ12>>) -> MatchFinder<ZZ12> {
        let ts = Arc::new(TileSet::new(rats));
        MatchFinder::new(ts)
    }

    #[test]
    fn spectre_self_finds_mystic() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&spectre());
        let mf = self_match_ts(vec![r]);

        let matches = mf.valid_matches(0, 0);
        assert!(
            matches
                .iter()
                .any(|m| mf.apply_match(m).seq() == MYSTIC_ZZ12),
            "mystic not found in {} matches: {:?}",
            matches.len(),
            matches,
        );
    }

    #[test]
    fn hexagon_self_single_edge() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&hexagon());
        let mf = self_match_ts(vec![r]);

        let matches = mf.valid_matches(0, 0);
        assert!(
            !matches.is_empty(),
            "hexagon should have valid self-matches"
        );

        let expected = &[-2, 2, 2, 2, 2, -2, 2, 2, 2, 2];
        assert!(
            matches.iter().any(|m| mf.apply_match(m).seq() == expected),
            "hex+hex match not found in {} matches",
            matches.len(),
        );
    }

    #[test]
    fn self_intersecting_filtered() {
        let r1 = Rat::<ZZ12>::from_slice_unchecked(&[0, 0, 3, 0, 3, 3, -3, -3, 3, 3, 0, 3]);
        let r2 = Rat::<ZZ12>::from_slice_unchecked(&[0, 0, 3, 3, 0, 0, 3, 3]);

        assert!(r1.try_glue((8, 0), &r2).is_err());

        let mf = self_match_ts(vec![r1, r2]);
        let matches = mf.valid_matches(0, 1);
        for m in &matches {
            assert!(
                !(m.a.range.start_offset == 8 && m.b.range.start_offset == 0),
                "self-intersecting glue should be filtered"
            );
        }
    }

    #[test]
    fn vertex_touching_rejected() {
        let r3 = Rat::<ZZ12>::try_from(
            &Snake::try_from(&[
                0, 0, 3, 2, 4, -3, -3, 3, 2, 4, -3, -3, 3, 0, 2, 4, -3, -3, 3, 2, 4, -3, 3, -3,
            ])
            .unwrap(),
        )
        .unwrap();
        let r4 = Rat::<ZZ12>::from_slice_unchecked(&[0, 0, 3, 3, 0, 0, 3, 3]);

        assert!(r3.try_glue((7, 0), &r4).is_err());

        let mf = self_match_ts(vec![r3, r4]);
        let matches = mf.valid_matches(0, 1);
        for m in &matches {
            assert!(
                !(m.a.range.start_offset == 7 && m.b.range.start_offset == 0),
                "vertex-touching glue should be filtered"
            );
        }
    }

    #[test]
    fn deterministic_order() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&hexagon());
        let mf = self_match_ts(vec![r]);
        let matches = mf.valid_matches(0, 0);
        for w in matches.windows(2) {
            assert!(
                w[0].a.range.start_offset <= w[1].a.range.start_offset
                    && (w[0].a.range.start_offset < w[1].a.range.start_offset || w[0].b.range.start_offset <= w[1].b.range.start_offset),
                "matches not in deterministic order"
            );
        }
    }

    #[test]
    fn multi_tile_all_valid_matches() {
        let sq = Rat::<ZZ12>::from_slice_unchecked(&[0, 1, 0, 1, 0, 1, 0, 1]);
        let tri = Rat::<ZZ12>::from_slice_unchecked(&[4, 4, 4]);

        let mf = self_match_ts(vec![sq, tri]);
        let all = mf.all_valid_matches();
        for w in all.windows(2) {
            assert!(
                w[0].a.tile_id < w[1].a.tile_id
                    || (w[0].a.tile_id == w[1].a.tile_id && w[0].b.tile_id < w[1].b.tile_id)
                    || (w[0].a.tile_id == w[1].a.tile_id
                        && w[0].b.tile_id == w[1].b.tile_id
                        && (w[0].a.range.start_offset < w[1].a.range.start_offset
                            || (w[0].a.range.start_offset == w[1].a.range.start_offset && w[0].b.range.start_offset <= w[1].b.range.start_offset))),
                "all_valid_matches not sorted"
            );
        }
    }

    #[test]
    fn spectre_foldback_cases_accepted() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&spectre());
        let cases: Vec<(usize, usize)> = vec![(0, 4), (0, 8), (6, 4), (6, 8)];
        for (ia, ib) in cases {
            assert!(
                r.try_glue((ia as i64, ib as i64), &r).is_ok(),
                "spectre foldback ({ia},{ib}) should be accepted by try_glue",
            );
        }
    }

    #[test]
    fn hexagon_pair_exhaustive() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&hexagon());
        let mf = self_match_ts(vec![r.clone()]);

        let mf_matches = mf.valid_matches(0, 0);
        let bf_matches = brute_force_valid_matches(&r, &r, 0, 0);
        assert_match_sets_match(
            &mf_matches,
            &bf_matches,
            mf.set_a().rats(),
            mf.set_b().rats(),
            "hexagon self",
        );

        for m in &mf_matches {
            assert_eq!(m.len(), 1, "hexagon pair match should be single-edge");
            assert_eq!(
                mf.apply_match(m).len(),
                10,
                "hex+hex result should be 10-gon"
            );
        }

        let canonical = mf.apply_match(&mf_matches[0]);
        for (idx, m) in mf_matches.iter().enumerate() {
            assert_eq!(
                mf.apply_match(m),
                canonical,
                "hex+hex match #{idx} should yield same shape",
            );
        }

        assert_eq!(
            mf_matches.len(),
            36,
            "hexagon pair should have 36 valid matches"
        );
    }

    #[test]
    fn hexamino_pair_exhaustive() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&hexagon());
        let hexamino = r.glue((0, 0), &r);

        let mf = self_match_ts(vec![hexamino.clone()]);
        let mf_matches = mf.valid_matches(0, 0);
        let bf_matches = brute_force_valid_matches(&hexamino, &hexamino, 0, 0);
        assert_match_sets_match(
            &mf_matches,
            &bf_matches,
            mf.set_a().rats(),
            mf.set_b().rats(),
            "hexamino self",
        );

        assert!(
            !mf_matches.is_empty(),
            "hexamino pair should have valid matches"
        );

        let multi_count = mf_matches.iter().filter(|m| m.len() > 1).count();
        let single_count = mf_matches.iter().filter(|m| m.len() == 1).count();
        assert!(
            multi_count > 0,
            "hexamino pair should have multi-edge matches",
        );
        assert!(
            single_count > 0,
            "hexamino pair should have single-edge matches",
        );
    }

    #[test]
    fn spectre_pair_exhaustive() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&spectre());
        let mf = self_match_ts(vec![r.clone()]);

        let mf_matches = mf.valid_matches(0, 0);
        let bf_matches = brute_force_valid_matches(&r, &r, 0, 0);
        assert_match_sets_match(
            &mf_matches,
            &bf_matches,
            mf.set_a().rats(),
            mf.set_b().rats(),
            "spectre self",
        );

        assert!(
            mf_matches
                .iter()
                .any(|m| mf.apply_match(m).seq() == MYSTIC_ZZ12),
            "spectre pair should find mystic, got {} matches",
            mf_matches.len(),
        );

        let multi_count = mf_matches.iter().filter(|m| m.len() > 1).count();
        let single_count = mf_matches.iter().filter(|m| m.len() == 1).count();
        assert!(
            multi_count > 0,
            "spectre pair should have multi-edge matches"
        );
        assert!(
            single_count > 0,
            "spectre pair should have single-edge matches"
        );
    }

    #[test]
    fn mixed_shapes_exhaustive() {
        let sq = Rat::<ZZ12>::from_slice_unchecked(&[0, 1, 0, 1, 0, 1, 0, 1]);
        let tri = Rat::<ZZ12>::from_slice_unchecked(&[4, 4, 4]);
        let mf = self_match_ts(vec![sq, tri]);

        for i in 0..2 {
            for j in 0..2 {
                let mf_matches = mf.valid_matches(i, j);
                let bf_matches = brute_force_valid_matches(mf.rat_a(i), mf.rat_b(j), i, j);
                assert_match_sets_match(
                    &mf_matches,
                    &bf_matches,
                    mf.set_a().rats(),
                    mf.set_b().rats(),
                    &format!("mixed shapes ({i},{j})"),
                );
            }
        }
    }

    #[test]
    fn single_edge_candidate_unit() {
        let hex: &[i8] = &[2, 2, 2, 2, 2, 2];
        assert!(
            junctions_glueable(hex, 0, 1, hex, 0),
            "hex (0,0): 2+2=4 > 0 on both sides"
        );
        assert!(
            junctions_glueable(hex, 0, 1, hex, 3),
            "hex (0,3): 2+2=4 > 0 on both sides"
        );

        let hexamino: &[i8] = &[-2, 2, 2, 2, 2, -2, 2, 2, 2, 2];
        assert!(
            !junctions_glueable(hexamino, 0, 1, hexamino, 0),
            "hexamino (0,0): left=-2+(-2)=-4 < 0, overflow"
        );
        assert!(
            !junctions_glueable(hexamino, 1, 1, hexamino, 1),
            "hexamino (1,1): right=2+(-2)=0, collinear extension"
        );
        assert!(
            junctions_glueable(hexamino, 1, 1, hexamino, 2),
            "hexamino (1,2): left=2+2=4, right=2+2=4, both gaps"
        );

        let spectre: &[i8] = &[3, 2, 0, 2, -3, 2, 3, 2, -3, 2, 3, -2, 3, -2];
        assert!(
            !junctions_glueable(spectre, 3, 1, spectre, 4),
            "spectre (3,4): left=2+(-3)=-1 < 0, overflow"
        );
        assert!(
            junctions_glueable(spectre, 0, 1, spectre, 1),
            "spectre (0,1): left=3+2=5, right=2+3=5, both gaps"
        );
    }

    #[test]
    fn interior_angle_heuristic_sound() {
        let hex = Rat::<ZZ12>::from_unchecked(&hexagon());
        let spec = Rat::<ZZ12>::from_unchecked(&spectre());
        let hexamino = hex.glue((0, 0), &hex);

        let tiles: Vec<(&str, &Rat<ZZ12>)> =
            vec![("hex", &hex), ("spec", &spec), ("hexamino", &hexamino)];

        for (label, rat) in &tiles {
            let n = rat.len();
            let seq = rat.seq();
            for ia in 0..n {
                for ib in 0..n {
                    let (_, len, _) = rat.get_match((ia as i64, ib as i64), rat);
                    if len != 1 {
                        continue;
                    }
                    if rat.try_glue((ia as i64, ib as i64), rat).is_ok() {
                        assert!(
                            junctions_glueable(seq, ia, 1, seq, ib),
                            "{label}: try_glue accepted single-edge ({ia},{ib}) but heuristic rejected",
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn empty_tile_valid_matches() {
        let hex = Rat::<ZZ12>::from_unchecked(&hexagon());
        let mf = self_match_ts(vec![hex]);
        assert!(!mf.valid_matches(0, 0).is_empty());
    }

    #[test]
    fn shared_boundaries_returns_cmi_matches() {
        let h: Snake<ZZ12> = hexagon();
        let r: Rat<ZZ12> = Rat::from_unchecked(&h);
        let mf = self_match_ts(vec![r]);

        let boundaries = mf.shared_boundaries(0, 0);
        assert!(
            boundaries.is_empty(),
            "hexagon self should have no RC matches: {:?}",
            boundaries,
        );

        let s: Snake<ZZ12> = spectre();
        let spec: Rat<ZZ12> = Rat::from_unchecked(&s);
        let mf2 = self_match_ts(vec![spec]);
        let self_bounds = mf2.shared_boundaries(0, 0);
        assert!(
            !self_bounds.is_empty(),
            "spectre self should have RC matches",
        );
        let max_len = self_bounds.iter().map(|m| m.len()).max().unwrap();
        assert_eq!(max_len, 3, "spectre max RC self-match should be 3");
    }

    #[test]
    fn heuristic_sound_on_size4_hexagon_patches() {
        let hex: Rat<ZZ12> = Rat::from_unchecked(&hexagon());

        let ts1 = Arc::new(TileSet::new(vec![hex.clone()]));
        let mf1 = MatchFinder::new(ts1);
        let size2: Vec<Rat<ZZ12>> = mf1
            .all_valid_matches()
            .iter()
            .map(|m| mf1.apply_match(m))
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();

        let ts2 = Arc::new(TileSet::new(size2.clone()));
        let mf2 = MatchFinder::new(ts2);
        let size4: Vec<Rat<ZZ12>> = mf2
            .all_valid_matches()
            .iter()
            .map(|m| mf2.apply_match(m))
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();

        let mut false_rejects = 0;

        for a in &size4 {
            for b in &size4 {
                let seq_a = a.seq();
                let seq_b = b.seq();
                for ia in 0..a.len() {
                    for ib in 0..b.len() {
                        let (_, len, _) = a.get_match((ia as i64, ib as i64), b);
                        if len != 1 {
                            continue;
                        }
                        if a.try_glue((ia as i64, ib as i64), b).is_err() {
                            continue;
                        }
                        if !junctions_glueable(seq_a, ia, 1, seq_b, ib) {
                            false_rejects += 1;
                        }
                    }
                }
            }
        }

        assert_eq!(
            false_rejects, 0,
            "heuristic falsely rejected {false_rejects} valid single-edge glues"
        );
    }

    #[test]
    fn size8_hexagon_valid_matches_matches_brute_force() {
        let hex: Rat<ZZ12> = Rat::from_unchecked(&hexagon());

        let ts1 = Arc::new(TileSet::new(vec![hex.clone()]));
        let mf1 = MatchFinder::new(ts1);
        let size2: Vec<Rat<ZZ12>> = mf1
            .all_valid_matches()
            .iter()
            .map(|m| mf1.apply_match(m))
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();

        let ts2 = Arc::new(TileSet::new(size2.clone()));
        let mf2 = MatchFinder::new(ts2);
        let size4: Vec<Rat<ZZ12>> = mf2
            .all_valid_matches()
            .iter()
            .map(|m| mf2.apply_match(m))
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();

        let ts4 = Arc::new(TileSet::new(size4.clone()));
        let mf4 = MatchFinder::new(ts4);

        let mut index_results: BTreeSet<Rat<ZZ12>> = BTreeSet::new();
        for i in 0..mf4.num_tiles_a() {
            for j in 0..mf4.num_tiles_b() {
                for m in mf4.valid_matches(i, j) {
                    index_results.insert(mf4.apply_match(&m));
                }
            }
        }

        let mut bf_results: BTreeSet<Rat<ZZ12>> = BTreeSet::new();
        for a in &size4 {
            for b in &size4 {
                for ia in 0..a.len() {
                    for ib in 0..b.len() {
                        if let Ok(glued) = a.try_glue((ia as i64, ib as i64), b) {
                            bf_results.insert(glued);
                        }
                    }
                }
            }
        }

        eprintln!(
            "index: {}, brute_force: {}",
            index_results.len(),
            bf_results.len(),
        );

        let missing: Vec<_> = bf_results.difference(&index_results).collect();
        assert_eq!(missing.len(), 0, "index missed {} patches", missing.len(),);
    }

    #[test]
    fn crossing_constructor_distinct_sets() {
        let hex: Rat<ZZ12> = Rat::from_unchecked(&hexagon());
        let spec: Rat<ZZ12> = Rat::from_unchecked(&spectre());

        let ts_a = Arc::new(TileSet::new(vec![hex.clone()]));
        let ts_b = Arc::new(TileSet::new(vec![spec.clone()]));
        let mf = MatchFinder::crossing(ts_a, ts_b);

        assert_eq!(mf.num_tiles_a(), 1);
        assert_eq!(mf.num_tiles_b(), 1);

        let matches = mf.valid_matches(0, 0);
        assert!(
            !matches.is_empty(),
            "hex x spectre crossing should have valid matches"
        );

        for m in &matches {
            assert_eq!(m.a.tile_id, 0);
            assert_eq!(m.b.tile_id, 0);
        }
    }

    #[test]
    fn crossing_self_matches_new_equivalent() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&hexagon());

        let ts = Arc::new(TileSet::new(vec![r.clone()]));
        let mf_self = MatchFinder::new(Arc::clone(&ts));
        let mf_cross = MatchFinder::crossing(Arc::clone(&ts), ts);

        let self_matches = mf_self.valid_matches(0, 0);
        let cross_matches = mf_cross.valid_matches(0, 0);

        assert_eq!(
            self_matches.len(),
            cross_matches.len(),
            "self and crossing should produce same number of matches"
        );

        for (s, c) in self_matches.iter().zip(cross_matches.iter()) {
            assert_eq!(
                (s.a.range.start_offset, s.len(), s.b.range.start_offset),
                (c.a.range.start_offset, c.len(), c.b.range.start_offset),
                "intervals should match"
            );
        }
    }

    #[test]
    fn crossing_multi_patch_against_seed() {
        let hex: Rat<ZZ12> = Rat::from_unchecked(&hexagon());
        let hex2 = hex.glue((0, 0), &hex);

        let ts_patches = Arc::new(TileSet::new(vec![hex.clone(), hex2]));
        let ts_seed = Arc::new(TileSet::new(vec![hex]));
        let mf = MatchFinder::crossing(ts_patches, ts_seed);

        assert_eq!(mf.num_tiles_a(), 2);
        assert_eq!(mf.num_tiles_b(), 1);

        let all = mf.all_valid_matches();
        assert!(
            !all.is_empty(),
            "multi-patch x seed should have valid matches"
        );

        for m in &all {
            assert!(m.a.tile_id < 2, "tile_a should index into patches");
            assert_eq!(m.b.tile_id, 0, "tile_b should index into seed");
        }
    }

    // --- MatchTypeIndex tests ---

    #[test]
    fn hexagon_self_match_types() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&hexagon());
        let ts = Arc::new(TileSet::new(vec![r]));
        let idx = MatchTypeIndex::new(ts);

        assert_eq!(
            idx.num_types(),
            21,
            "hexagon self should have 21 match types"
        );

        for id in 1..=idx.num_types() {
            let mt = idx.get(id);
            assert_eq!((mt.a.tile_id, mt.b.tile_id), (0, 0));
            assert_eq!(mt.len(), 1);
        }
    }

    #[test]
    fn square_self_match_types() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&square());
        let ts = Arc::new(TileSet::new(vec![r]));
        let idx = MatchTypeIndex::new(ts);

        assert_eq!(
            idx.num_types(),
            10,
            "square self should have 10 match types"
        );

        for id in 1..=idx.num_types() {
            let mt = idx.get(id);
            assert_eq!((mt.a.tile_id, mt.b.tile_id), (0, 0));
            assert_eq!(mt.len(), 1);
        }
    }

    #[test]
    fn spectre_self_many_types() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&spectre());
        let ts = Arc::new(TileSet::new(vec![r]));
        let idx = MatchTypeIndex::new(ts);

        assert!(
            idx.num_types() > 10,
            "spectre self should have many match types, got {}",
            idx.num_types()
        );
    }

    #[test]
    fn signed_id_forward_and_backward() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&spectre());
        let ts = Arc::new(TileSet::new(vec![r]));
        let idx = MatchTypeIndex::new(ts);

        let mut found_asymmetric = false;
        for id in 1..=idx.num_types() {
            let mt = idx.get(id);
            if mt.a.range.start_offset == mt.b.range.start_offset {
                continue;
            }
            let len_a = idx.tileset().rat(mt.a.tile_id).len();
            let len_b = idx.tileset().rat(mt.b.tile_id).len();
            let inv = match_involution(&mt, len_a, len_b);
            let sid_fwd = idx.signed_id(&mt).unwrap();
            let sid_inv = idx.signed_id(&inv).unwrap();
            assert_eq!(sid_fwd.abs(), sid_inv.abs(), "same unsigned id");
            if inv == mt {
                assert_eq!(sid_fwd, sid_inv, "self-involutive: both same sign");
                assert!(sid_fwd > 0, "should be positive");
            } else {
                assert_eq!(sid_fwd, -sid_inv, "asymmetric: opposite signs");
                assert!(sid_fwd > 0, "first should be positive");
                found_asymmetric = true;
            }
        }
        assert!(
            found_asymmetric,
            "should find at least one asymmetric match type"
        );
    }

    #[test]
    fn hexagon_self_signed_ids_by_position() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&hexagon());
        let ts = Arc::new(TileSet::new(vec![r]));
        let idx = MatchTypeIndex::new(ts);

        let mut found_ids: HashSet<i32> = HashSet::new();
        for ia in 0..6 {
            for ib in 0..6 {
                let mt = TileMatch::new(
                    Segment::new(0, EdgeRange::new(ia, 1)),
                    Segment::new(0, EdgeRange::new(ib, 1)),
                );
                if let Some(id) = idx.signed_id(&mt) {
                    found_ids.insert(id);
                }
            }
        }
        assert!(
            !found_ids.is_empty(),
            "should find some signed ids for hex self-matches"
        );
    }

    #[test]
    fn apply_produces_valid_rat() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&spectre());
        let ts = Arc::new(TileSet::new(vec![r]));
        let idx = MatchTypeIndex::new(ts);

        for id in 1..=idx.num_types() {
            let result_pos = idx.apply(id as i32);
            assert!(
                Snake::<ZZ12>::try_from(result_pos.seq()).is_ok(),
                "apply(+{id}) should produce valid snake"
            );
            let result_neg = idx.apply(-(id as i32));
            assert!(
                Snake::<ZZ12>::try_from(result_neg.seq()).is_ok(),
                "apply(-{id}) should produce valid snake"
            );
            assert_eq!(
                result_pos, result_neg,
                "forward and backward apply should produce same rat"
            );
        }
    }

    #[test]
    fn bi_hex_bi_square_cross_matches() {
        let hex: Rat<ZZ12> = Rat::from_unchecked(&hexagon());
        let sq: Rat<ZZ12> = Rat::from_unchecked(&square());
        let bi_hex = hex.glue((0, 0), &hex);
        let bi_sq = sq.glue((0, 0), &sq);

        let ts = Arc::new(TileSet::new(vec![bi_hex.clone(), bi_sq.clone()]));
        let idx = MatchTypeIndex::new(Arc::clone(&ts));

        assert_eq!(ts.num_tiles(), 2);

        let n = idx.num_types();
        assert!(n > 0, "should have match types");

        let mut has_self_hex = false;
        let mut has_self_sq = false;
        let mut has_cross = false;
        for id in 1..=n {
            let mt = idx.get(id);
            if mt.a.tile_id == 0 && mt.b.tile_id == 0 {
                has_self_hex = true;
            }
            if mt.a.tile_id == 1 && mt.b.tile_id == 1 {
                has_self_sq = true;
            }
            if mt.a.tile_id != mt.b.tile_id {
                has_cross = true;
            }
        }
        assert!(has_self_hex, "should have bi-hex self-match types");
        assert!(has_self_sq, "should have bi-sq self-match types");
        assert!(has_cross, "should have cross-match types");
    }

    #[test]
    fn signed_id_roundtrip_bi_hex_bi_square() {
        let hex: Rat<ZZ12> = Rat::from_unchecked(&hexagon());
        let sq: Rat<ZZ12> = Rat::from_unchecked(&square());
        let bi_hex = hex.glue((0, 0), &hex);
        let bi_sq = sq.glue((0, 0), &sq);

        let ts = Arc::new(TileSet::new(vec![bi_hex, bi_sq]));
        let idx = MatchTypeIndex::new(Arc::clone(&ts));

        let mf = MatchFinder::new(Arc::clone(&ts));
        let all = mf.all_valid_matches();

        let mut signed_ids: HashSet<i32> = HashSet::new();
        for mt in &all {
            let sid = idx
                .signed_id(mt)
                .unwrap_or_else(|| panic!("match {:?} not found in index", mt));
            signed_ids.insert(sid);

            let result_mf = mf.apply_match(mt);
            let result_idx = idx.apply(sid);
            assert_eq!(
                result_mf, result_idx,
                "MatchFinder and MatchTypeIndex should produce same rat for {:?}",
                mt
            );
        }

        assert!(
            signed_ids.len() >= idx.num_types(),
            "should have at least one signed id per type"
        );

        for &sid in &signed_ids {
            let unsigned = sid.unsigned_abs() as usize;
            assert!(unsigned >= 1 && unsigned <= idx.num_types());
        }
    }

    #[test]
    fn get_match_involution_semantics() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&spectre());
        let n = r.len();

        let mf = self_match_ts(vec![r.clone()]);
        let matches = mf.valid_matches(0, 0);

        let mut checked = 0;
        for mt in &matches {
            if mt.a.range.start_offset == mt.b.range.start_offset {
                continue;
            }

            let fwd = (mt.a.range.start_offset as i64, mt.len(), mt.b.range.start_offset as i64);
            let rat_fwd = r
                .try_glue_precomputed(fwd, &r, true)
                .expect("forward glue should work");

            let ns = mt.a.range.start_offset as i64;
            let ne = mt.b.range.start_offset as i64;
            let len = mt.len() as i64;

            let candidates: Vec<(&str, (i64, i64))> = vec![
                ("(ne, ns)", (ne, ns)),
                ("(ne, ns+len)", (ne, ns + len)),
                ("(ne, ns-len)", (ne, ns - len)),
                ("(ne-len, ns)", (ne - len, ns)),
                ("(ne-len, ns+len)", (ne - len, ns + len)),
                ("(ne+len, ns)", (ne + len, ns)),
                ("(ne+len, ns-len)", (ne + len, ns - len)),
            ];

            let mut found = false;
            for (label, (s, e)) in &candidates {
                let (ns2, len2, ne2) = r.get_match((*s, *e), &r);
                if len2 != mt.len() {
                    continue;
                }
                if let Ok(rat_rev) = r.try_glue_precomputed((ns2, len2, ne2), &r, true) {
                    if rat_rev == rat_fwd {
                        eprintln!(
                            "  MATCH {} works: fwd=({},{},{}) rev=({},{},{}) rat_fwd==rat_rev",
                            label, ns, mt.len(), ne, ns2, len2, ne2
                        );
                        found = true;

                        assert_eq!(
                            ns2,
                            (ne - len).rem_euclid(n as i64),
                            "reverse ext_start should be (ne - len) % n for {}",
                            label
                        );
                        assert_eq!(
                            ne2,
                            (ns + len).rem_euclid(n as i64),
                            "reverse ext_end should be (ns + len) % n for {}",
                            label
                        );
                    }
                }
            }

            assert!(
                found,
                "no reverse candidate produced same rat for fwd=({ns},{len},{ne})"
            );

            checked += 1;
            if checked >= 5 {
                break;
            }
        }
        assert!(
            checked > 0,
            "need at least one match with start_a != start_b"
        );
    }

    #[test]
    fn hexagon_self_apply_both_signs_same_rat() {
        let hex: Rat<ZZ12> = Rat::from_unchecked(&hexagon());
        let ts = Arc::new(TileSet::new(vec![hex]));
        let idx = MatchTypeIndex::new(ts);

        assert!(idx.num_types() > 0, "hex self should have match types");

        for id in 1..=idx.num_types() {
            let r_pos = idx.apply(id as i32);
            let r_neg = idx.apply(-(id as i32));
            assert!(
                Snake::<ZZ12>::try_from(r_pos.seq()).is_ok(),
                "apply(+{id}) should produce valid snake"
            );
            assert!(
                Snake::<ZZ12>::try_from(r_neg.seq()).is_ok(),
                "apply(-{id}) should produce valid snake"
            );
            assert_eq!(
                r_pos, r_neg,
                "apply(+{id}) and apply(-{id}) should produce same rat"
            );
        }
    }

    #[test]
    fn hex_square_cross_apply_both_signs_same_rat() {
        let hex: Rat<ZZ12> = Rat::from_unchecked(&hexagon());
        let sq: Rat<ZZ12> = Rat::from_unchecked(&square());
        let ts = Arc::new(TileSet::new(vec![hex, sq]));
        let idx = MatchTypeIndex::new(ts);

        let cross_ids: Vec<usize> = (1..=idx.num_types())
            .filter(|id| {
                let mt = idx.get(*id);
                mt.a.tile_id != mt.b.tile_id
            })
            .collect();

        assert!(
            !cross_ids.is_empty(),
            "hex+square should have cross-match types"
        );

        for id in cross_ids {
            let r_pos = idx.apply(id as i32);
            let r_neg = idx.apply(-(id as i32));
            assert!(
                Snake::<ZZ12>::try_from(r_pos.seq()).is_ok(),
                "apply(+{id}) should produce valid snake"
            );
            assert!(
                Snake::<ZZ12>::try_from(r_neg.seq()).is_ok(),
                "apply(-{id}) should produce valid snake"
            );
            assert_eq!(
                r_pos, r_neg,
                "apply(+{id}) and apply(-{id}) should produce same rat"
            );
        }
    }

    #[test]
    fn hex_bisq_match_types() {
        let hex: Rat<ZZ12> = Rat::from_unchecked(&hexagon());
        let sq: Rat<ZZ12> = Rat::from_unchecked(&square());
        let bi_sq = sq.glue((0, 0), &sq);
        let ts = Arc::new(TileSet::new(vec![hex, bi_sq]));
        let idx = MatchTypeIndex::new(Arc::clone(&ts));

        let hex_idx = ts
            .index_of(&Rat::<ZZ12>::from_unchecked(&hexagon()))
            .unwrap();
        let bisq_idx = ts.index_of(&sq.glue((0, 0), &sq)).unwrap();
        assert_ne!(hex_idx, bisq_idx, "tiles should be distinct");

        let mut hex_self = 0usize;
        let mut bisq_self = 0usize;
        let mut cross = 0usize;
        for id in 1..=idx.num_types() {
            let mt = idx.get(id);
            if mt.a.tile_id == hex_idx && mt.b.tile_id == hex_idx {
                hex_self += 1;
            } else if mt.a.tile_id == bisq_idx && mt.b.tile_id == bisq_idx {
                bisq_self += 1;
            } else {
                cross += 1;
            }

            let r_pos = idx.apply(id as i32);
            let r_neg = idx.apply(-(id as i32));
            assert!(
                Snake::<ZZ12>::try_from(r_pos.seq()).is_ok(),
                "apply(+{id}) should produce valid snake"
            );
            assert!(
                Snake::<ZZ12>::try_from(r_neg.seq()).is_ok(),
                "apply(-{id}) should produce valid snake"
            );
            assert_eq!(
                r_pos, r_neg,
                "apply(+{id}) and apply(-{id}) should produce same rat"
            );
        }

        assert!(hex_self > 0, "should have hex self-match types");
        assert!(bisq_self > 0, "should have bisq self-match types");
        assert!(cross > 0, "should have cross-match types");
    }

    #[test]
    fn candidates_starting_at_covers_all_valid_matches() {
        let hex_rat = Rat::try_from(&hexagon::<ZZ12>()).unwrap();
        let sq_rat = Rat::try_from(&square::<ZZ12>()).unwrap();
        let ts = Arc::new(TileSet::new(vec![hex_rat, sq_rat]));
        let mf = MatchFinder::new(Arc::clone(&ts));
        let idx = MatchTypeIndex::new(ts.clone());

        let all_matches = mf.all_valid_matches();
        let mut by_start_count = 0usize;
        for tile_a in 0..2 {
            for offset in 0..idx.all_candidates_for_tile(tile_a).len() {
                by_start_count += idx.candidates_starting_at(tile_a, offset).len();
            }
        }

        let all_inv: Vec<TileMatch> = all_matches
            .iter()
            .flat_map(|mt| {
                let inv = match_involution(mt, ts.rat(mt.a.tile_id).len(), ts.rat(mt.b.tile_id).len());
                vec![*mt, inv]
            })
            .collect();
        assert_eq!(
            by_start_count,
            all_inv.len(),
            "secondary index should have same total count as all matches (including involutions)"
        );

        for mt in &all_inv {
            let found = idx
                .candidates_starting_at(mt.a.tile_id, mt.a.range.start_offset)
                .iter()
                .any(|c| c.tile_id == mt.b.tile_id && c.range.start_offset == mt.b.range.start_offset && c.len() == mt.len());
            assert!(
                found,
                "match ({}, {}, {}, {}, {}) not found in secondary index",
                mt.a.tile_id, mt.a.range.start_offset, mt.b.tile_id, mt.b.range.start_offset, mt.len()
            );
        }
    }

    #[test]
    fn match_type_index_matches_brute_force_for_all_tilesets() {
        let single_tile_cases: Vec<(&str, Rat<ZZ12>)> = vec![
            ("spectre self", Rat::try_from(&spectre()).unwrap()),
            ("hexagon self", Rat::try_from(&hexagon()).unwrap()),
            ("square self", Rat::try_from(&square()).unwrap()),
        ];

        for (label, rat) in &single_tile_cases {
            let ts = Arc::new(TileSet::new(vec![rat.clone()]));
            let mf = MatchFinder::new(Arc::clone(&ts));
            let idx = MatchTypeIndex::new(Arc::clone(&ts));

            let mf_matches = mf.valid_matches(0, 0);
            let bf_matches = brute_force_valid_matches(rat, rat, 0, 0);
            assert_match_sets_match(
                &mf_matches,
                &bf_matches,
                mf.set_a().rats(),
                mf.set_b().rats(),
                label,
            );

            let all_mf = mf.all_valid_matches();
            let mut mf_sorted = all_mf.clone();
            mf_sorted.sort_by(|x, y| {
                x.a.range.start_offset
                    .cmp(&y.a.range.start_offset)
                    .then_with(|| x.b.range.start_offset.cmp(&y.b.range.start_offset))
            });

            let mut bf_all = bf_matches.clone();
            bf_all.sort_by(|x, y| {
                x.a.range.start_offset
                    .cmp(&y.a.range.start_offset)
                    .then_with(|| x.b.range.start_offset.cmp(&y.b.range.start_offset))
            });

            assert_eq!(
                mf_sorted, bf_all,
                "{label}: MatchFinder.all_valid_matches differs from brute force"
            );

            let mut idx_all: Vec<TileMatch> = Vec::new();
            for entry in &idx.entries {
                idx_all.push(entry.first);
                idx_all.push(entry.second);
            }
            idx_all.sort_by(|x, y| {
                x.a.range.start_offset
                    .cmp(&y.a.range.start_offset)
                    .then_with(|| x.b.range.start_offset.cmp(&y.b.range.start_offset))
            });
            idx_all.dedup();

            assert_eq!(
                idx_all.len(),
                bf_all.len(),
                "{label}: MatchTypeIndex has {} unique matches but brute force found {}",
                idx_all.len(),
                bf_all.len(),
            );

            for (i, (idx_mt, bf_mt)) in idx_all.iter().zip(bf_all.iter()).enumerate() {
                assert_eq!(
                    (idx_mt.a.range.start_offset, idx_mt.len(), idx_mt.b.range.start_offset),
                    (bf_mt.a.range.start_offset, bf_mt.len(), bf_mt.b.range.start_offset),
                    "{label}: tuple mismatch at index {i}",
                );
            }

            for offset in 0..rat.len() {
                for cm in idx.candidates_starting_at(0, offset) {
                    let found_direct = bf_all.iter().any(|mt| {
                        mt.a.range.start_offset == offset
                            && mt.b.tile_id == 0
                            && mt.b.range.start_offset == cm.range.start_offset
                            && mt.len() == cm.len()
                    });
                    let found_inv = bf_all.iter().any(|mt| {
                        mt.a.range.start_offset == cm.range.start_offset
                            && mt.b.tile_id == 0
                            && mt.b.range.start_offset == offset
                            && mt.len() == cm.len()
                    });
                    assert!(
                        found_direct || found_inv,
                        "{label}: secondary (0,{offset})->(sb={},len={}) not in brute force",
                        cm.range.start_offset,
                        cm.len(),
                    );
                }
            }
        }
    }
}
