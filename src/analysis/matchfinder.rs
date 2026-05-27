//! Match enumeration between two tile boundaries.
//!
//! [`MatchFinder`] is the operational layer between
//! [`crate::stringmatch::BitParallelMatcher`] (raw cyclic RC matching
//! over angle slices) and the geometric validation in
//! [`crate::geom`] (Snake closure, junction angles, glue checks).
//! It enumerates legal glue matches between two tilesets (or a
//! tileset against itself) and filters the BP candidates down to
//! geometrically realisable glues.
//!
//! [`BpSeed`] is the precomputed BP state for a fixed B-side
//! tileset; pass it to [`MatchFinder::crossing_with_seed`] to
//! amortise the matcher build across many A-side batches.
//!
//! The companion module [`crate::analysis::matchtypes`] hosts
//! [`crate::analysis::matchtypes::MatchTypeIndex`] -- the *ontological*
//! catalog of distinct match types pre-computed once for a
//! tileset, indexed for repeated O(1) lookup. `MatchTypeIndex` is
//! built on top of `MatchFinder`; the split here keeps the
//! "find matches" mechanism separate from the "name and catalog
//! the distinct match types" layer above it.

use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

use rustc_hash::FxHashSet;

use crate::analysis::matchtypes::apply_match;
use crate::cyclotomic::IsRing;
use crate::geom::glue::junctions_glueable;
use crate::geom::matches::{EdgeRange, Segment, TileMatch};
use crate::geom::rat::Rat;
use crate::geom::snake::Snake;
use crate::geom::tileset::TileSet;
use crate::stringmatch::BitParallelMatcher;

/// Pre-built bit-parallel state for a fixed "seed" (B-side) tileset.
///
/// Bundles an `Arc<TileSet<T>>` with the [`BitParallelMatcher`] built
/// from its angle sequences. Use this when many short-lived
/// [`MatchFinder`]s will share the same B-side: build a `BpSeed` once,
/// then call [`MatchFinder::crossing_with_seed`] for each A-side
/// without re-running the BP precompute (~`O(num_b_tiles * tile_len)`
/// per call collapses to a single setup).
///
/// Cheap to clone: both fields are `Arc`s. Note: the manual `Clone`
/// impl (rather than `#[derive(Clone)]`) avoids a spurious `T: Clone`
/// bound — the underlying `T: IsRing` data only lives inside the
/// `Arc`s and never needs to be cloned.
///
/// # Canonical usage (patch_enum / seq_explorer pattern)
///
/// ```ignore
/// let seed = BpSeed::new(Arc::clone(&fixed_b_tileset));
/// for batch in /* varying A-side batches */ {
///     let a_ts = Arc::new(TileSet::new(batch));
///     let mf = MatchFinder::crossing_with_seed(a_ts, seed.clone());
///     for m in mf.all_valid_matches() { /* ... */ }
/// }
/// ```
///
/// `seed.clone()` cheaply hands out new `Arc` handles; the matcher
/// masks are shared across all `MatchFinder`s.
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
    /// Pre-compute the bit-parallel masks for the given B-side
    /// tileset. The matcher is wrapped in `Arc` so subsequent
    /// `clone()`s share the masks.
    pub fn new(tileset: Arc<TileSet<T>>) -> Self {
        let sequences: Vec<Vec<i8>> = tileset.rats().iter().map(|r| r.seq().to_vec()).collect();
        let matcher = Arc::new(BitParallelMatcher::new(&sequences));
        BpSeed { tileset, matcher }
    }

    /// The B-side tileset this seed indexes.
    pub fn tileset(&self) -> &Arc<TileSet<T>> {
        &self.tileset
    }
}

/// Enumerates legal glue matches between two tilesets (or one tileset
/// against itself).
///
/// # When to use which constructor
///
/// | Scenario | Constructor | Cost |
/// |---|---|---|
/// | one tileset vs itself | [`Self::new`] | builds BP masks once |
/// | two distinct tilesets, one-shot | [`Self::crossing`] | builds BP masks once for B |
/// | two distinct tilesets, many builds reusing the B-side | [`BpSeed::new`] + [`Self::crossing_with_seed`] | BP masks built once, shared by `Arc` across builds |
///
/// The amortised path matters when you iterate over many A-side
/// tilesets against a fixed B-side, which is the dominant pattern in
/// `patch_enum`, `seq_explorer`, and similar callers (a fixed tile
/// alphabet matched against a stream of grown patches).
///
/// # Tile ids
///
/// The returned [`TileMatch`]es carry **side-local** tile ids:
/// `m.a.tile_id` indexes into `set_a.rats()` (range `[0,
/// num_tiles_a)`) and `m.b.tile_id` indexes into `set_b.rats()`
/// (range `[0, num_tiles_b)`). For `Self::new` (self-match), the two
/// sides are the same tileset so the ranges coincide.
///
/// # Internals
///
/// Backed by a [`BitParallelMatcher`] built once from the B-side
/// tileset's angle sequences; A-side sequences are streamed against
/// it per query, reading directly from `set_a.rat(i).seq()` (no
/// per-MatchFinder clone of the angle data). Build cost is
/// dominated by the B-side BP precompute, which is why the `BpSeed`
/// reuse path exists.
pub struct MatchFinder<T: IsRing> {
    set_a: Arc<TileSet<T>>,
    set_b: Arc<TileSet<T>>,
    b_matcher: Arc<BitParallelMatcher>,
}

impl<T: IsRing> MatchFinder<T> {
    /// Build a `MatchFinder` that enumerates self-matches within a
    /// single tileset (`set_a == set_b == tileset`).
    ///
    /// Both `tile_a` and `tile_b` in the returned `TileMatch`es index
    /// into the same `tileset.rats()`. Use this for the pairwise
    /// "which tiles glue to which" enumeration over one set.
    pub fn new(tileset: Arc<TileSet<T>>) -> Self {
        let sequences: Vec<Vec<i8>> = tileset.rats().iter().map(|r| r.seq().to_vec()).collect();
        let b_matcher = Arc::new(BitParallelMatcher::new(&sequences));
        MatchFinder {
            set_a: Arc::clone(&tileset),
            set_b: tileset,
            b_matcher,
        }
    }

    /// Build a `MatchFinder` that enumerates matches between two
    /// distinct tilesets `a` and `b`. One-shot convenience — builds
    /// the BP masks for `b` inline.
    ///
    /// For callers that iterate many `MatchFinder` builds against
    /// the **same** B-side, use [`BpSeed::new`] +
    /// [`Self::crossing_with_seed`] instead so the precompute is
    /// done once and shared via `Arc`.
    pub fn crossing(a: Arc<TileSet<T>>, b: Arc<TileSet<T>>) -> Self {
        Self::crossing_with_seed(a, BpSeed::new(b))
    }

    /// Build a crossing `MatchFinder` reusing a pre-built [`BpSeed`]
    /// for the B-side.
    ///
    /// Per-call cost is now constant beyond the `Arc` bookkeeping —
    /// `MatchFinder` reads A-side angle sequences directly from the
    /// `set_a` tileset rather than cloning them, so the only state
    /// per build is the two `Arc<TileSet<T>>` handles and the
    /// `Arc<BitParallelMatcher>`. This is the right primitive for
    /// the inner loop of `patch_enum`, `seq_explorer`, and similar
    /// callers — see the [`BpSeed`] doc for the canonical usage
    /// pattern.
    pub fn crossing_with_seed(a: Arc<TileSet<T>>, seed: BpSeed<T>) -> Self {
        MatchFinder {
            set_a: a,
            set_b: seed.tileset,
            b_matcher: seed.matcher,
        }
    }

    /// Enumerate all maximal RC matches between A-tile `i` and B-tile
    /// `j` (both side-local indices). Returned `TileMatch`es have
    /// `m.a.tile_id == i` and `m.b.tile_id == j`. No geometric
    /// validation -- caller is responsible for downstream filters.
    fn maximal_rc_matches(&self, i: usize, j: usize) -> Vec<TileMatch> {
        let mut matches = self.b_matcher.stream_boundary(self.set_a.rat(i).seq(), j);
        for m in &mut matches {
            m.a.tile_id = i;
        }
        matches
    }

    /// Same as [`Self::maximal_rc_matches`] but post-filtered to
    /// matches whose A-side start position is in `positions`.
    fn maximal_rc_matches_at_positions(
        &self,
        i: usize,
        j: usize,
        positions: &[usize],
    ) -> Vec<TileMatch> {
        if positions.is_empty() {
            return vec![];
        }
        let boundary = self.set_a.rat(i).seq();
        let n_a = boundary.len();
        let mut keep = vec![false; n_a];
        for &p in positions {
            if p < n_a {
                keep[p] = true;
            }
        }
        let mut matches = self.b_matcher.stream_boundary(boundary, j);
        matches.retain(|m| keep[m.a.range.start_offset]);
        for m in &mut matches {
            m.a.tile_id = i;
        }
        matches
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
        self.maximal_rc_matches(i, j)
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

        let cmi_matches = self.maximal_rc_matches_at_positions(i, j, &scan_positions);
        for m in &cmi_matches {
            let (ns, len, ne) = a.get_match(
                (m.a.range.start_offset as i64, m.b.range.start_offset as i64),
                b,
            );
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
                            Segment::new(
                                j,
                                EdgeRange::new(ne.rem_euclid(n_b as i64) as usize, len),
                            ),
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
        let n_a = a.len();
        let n_b = b.len();

        let mut groups: BTreeMap<Rat<T>, Vec<TileMatch>> = BTreeMap::new();
        if n_a == 0 || n_b == 0 {
            return groups;
        }

        let cmi_matches = self.maximal_rc_matches(i, j);
        for m in &cmi_matches {
            let (ns, len, ne) = a.get_match(
                (m.a.range.start_offset as i64, m.b.range.start_offset as i64),
                b,
            );
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
