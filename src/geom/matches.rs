//! Layered data types for boundary matches.
//!
//! The match concept decomposes into three layers, each adding caller
//! context to the previous one:
//!
//! 1. **Pure positional** ([`EdgeRange`], [`Match`]). An `EdgeRange` is
//!    just `start_offset + len` on some unspecified boundary. A `Match`
//!    is a pair of equal-length `EdgeRange`s — "these positions line up".
//!    No tile ids, no rings, no anchoring.
//!
//! 2. **Half-bound** ([`Segment`], [`PatchMatch`]). A `Segment` attaches a
//!    tile id to an `EdgeRange`. A `PatchMatch` is asymmetric: an anonymous
//!    A-side `EdgeRange` (the patch's running boundary, no tile id) plus a
//!    B-side `Segment` (the new tile being glued in).
//!
//! 3. **Fully bound** ([`TileMatch`]). A pair of `Segment`s: both sides
//!    have tile ids, drawn from one or two `TileSet`s.
//!
//! Conversions:
//!
//! - `TileMatch::spec()` / `PatchMatch::spec()` → strip tile ids, yielding
//!   the positional `Match`.
//! - `Match::with_tiles(a_tile, b_tile)` → lift to `TileMatch`.
//! - `Match::with_b_tile(b_tile)` → lift to `PatchMatch`.
//!
//! The point of this layering is that function signatures stop carrying
//! information they don't use: `Rat`-level operations work in `Match`
//! (no tile ids), `MatchFinder`-level operations work in `TileMatch`,
//! and the `GrowingPatch` API works in `PatchMatch`. Each layer surfaces
//! exactly the information that layer actually has.
//!
//! Semantics of `start_offset`: this module is agnostic — `EdgeRange`
//! is purely positional. The consumer documents what edges the range
//! refers to (matched edges, surviving edges past the match, etc.).
//! When wiring in the existing callers, we settle on "`start_offset`
//! is the first matched edge" as the universal convention.

use serde::{Deserialize, Serialize};

// ============================================================
// Layer 1: pure positional
// ============================================================

/// Half-open range of edges on some boundary: `[start_offset, start_offset + len)`.
///
/// Indices wrap modulo the boundary length when the consumer is cyclic;
/// this type itself makes no assumption about wrap-around.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub struct EdgeRange {
    pub start_offset: usize,
    pub len: usize,
}

impl EdgeRange {
    pub fn new(start_offset: usize, len: usize) -> Self {
        Self { start_offset, len }
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }
}

/// Two equal-length edge ranges on two unspecified boundaries that line up.
///
/// "Line up" means the consumer treats the `len` edges in `a` and the
/// `len` edges in `b` as compatible (e.g. revcomp-equivalent for gluing).
/// No tile ids, no ring information.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub struct Match {
    pub a: EdgeRange,
    pub b: EdgeRange,
}

impl Match {
    /// Construct from two `EdgeRange`s. In debug builds, asserts that
    /// `a.len == b.len`; otherwise behavior on mismatched lengths is
    /// implementation-defined.
    pub fn new(a: EdgeRange, b: EdgeRange) -> Self {
        debug_assert_eq!(
            a.len, b.len,
            "Match: a and b must have equal length (got a.len={}, b.len={})",
            a.len, b.len,
        );
        Self { a, b }
    }

    /// The match length (same on both sides).
    pub fn len(&self) -> usize {
        debug_assert_eq!(self.a.len, self.b.len);
        self.a.len
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The involution of this match: the same glue described from the
    /// other boundary's perspective.
    ///
    /// **Not a trivial A/B swap.** Under the asymmetric `start_offset`
    /// convention (see [`PatchMatch`]'s docstring), `a.start_offset` is
    /// the first matched A-edge while `b.start_offset` is the first
    /// surviving B-edge past the match. After swapping roles, both
    /// sides need their offsets shifted to land on the correct anchor
    /// under the new role:
    ///
    /// - new A start = (old B start + len_b - len) mod len_b
    /// - new B start = (old A start + len)         mod len_a
    ///
    /// `len_a` / `len_b` are the cyclic boundary lengths of the two
    /// sides. Idempotent under the same input lengths:
    /// `m.involution(la, lb).involution(lb, la) == m`.
    ///
    /// This mirrors [`TileMatch::involution`]; the only difference is
    /// that `Match` is purely positional (no tile ids) so the boundary
    /// lengths are passed in directly rather than looked up via tile
    /// ids.
    pub fn involution(self, len_a: usize, len_b: usize) -> Self {
        let len = self.len();
        Self {
            a: EdgeRange::new((self.b.start_offset + len_b - len) % len_b, len),
            b: EdgeRange::new((self.a.start_offset + len) % len_a, len),
        }
    }

    /// Attach tile ids on both sides, lifting to a [`TileMatch`].
    pub fn with_tiles(self, a_tile: usize, b_tile: usize) -> TileMatch {
        TileMatch {
            a: Segment {
                tile_id: a_tile,
                range: self.a,
            },
            b: Segment {
                tile_id: b_tile,
                range: self.b,
            },
        }
    }

    /// Attach a tile id only on the B side, lifting to a [`PatchMatch`].
    /// The A side stays anonymous (= the patch's own boundary).
    pub fn with_b_tile(self, b_tile: usize) -> PatchMatch {
        PatchMatch {
            a_range: self.a,
            b: Segment {
                tile_id: b_tile,
                range: self.b,
            },
        }
    }
}

// ============================================================
// Layer 2: half-bound
// ============================================================

/// Tile id + edge range — "where + which".
///
/// `tile_id` is meaningful relative to some `TileSet` the caller knows
/// about; this type does not name the tileset.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub struct Segment {
    pub tile_id: usize,
    pub range: EdgeRange,
}

impl Segment {
    pub fn new(tile_id: usize, range: EdgeRange) -> Self {
        Self { tile_id, range }
    }

    pub fn len(&self) -> usize {
        self.range.len
    }

    pub fn is_empty(&self) -> bool {
        self.range.is_empty()
    }
}

/// A maximal contiguous range of edges on a **patch boundary**,
/// anchored to the tile-segment that occupies it.
///
/// Where [`Segment`] is "a range of edges on a tile" (anchored by
/// `tile_id` to a tile in some `TileSet`), `PatchSegment` is the
/// analogous patch-side concept: "a range of edges on the patch's
/// boundary, identified with the tile-segment whose edges occupy
/// those positions".
///
/// # Convention
///
/// Within the patch range, position `k` (counting from
/// `range.start_offset`) corresponds to tile-offset
/// `(tile_seg.range.start_offset + k) mod tile_len`. Forward
/// identification on both sides -- the patch edge at the given
/// position **is** that tile edge (same edge in space, not glued).
/// Contrast [`PatchMatch`] (same shape but different semantic: the
/// edges *would* match anti-parallel under a glue, not currently
/// identical).
///
/// # Cyclic-vs-linear caveat
///
/// The patch boundary is cyclic but `PatchSegment`s are typically
/// produced from the linear array `[0, n)`. If position 0 is not a
/// junction, a single cyclic tile-instance run that straddles
/// position 0 is split into two linear `PatchSegment`s -- one at
/// the start and one at the end. Both have the same
/// `tile_seg.tile_id`, and their tile-offsets are cyclically
/// continuous. Callers that need true cyclic runs must stitch these
/// two halves themselves.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub struct PatchSegment {
    /// Range on the patch boundary.
    pub range: EdgeRange,
    /// The tile-segment whose edges occupy those patch positions.
    /// `tile_seg.range.len` equals `range.len` by construction.
    pub tile_seg: Segment,
}

impl PatchSegment {
    pub fn new(range: EdgeRange, tile_seg: Segment) -> Self {
        debug_assert_eq!(
            range.len, tile_seg.range.len,
            "PatchSegment: range and tile_seg.range must have equal length",
        );
        Self { range, tile_seg }
    }

    pub fn len(&self) -> usize {
        debug_assert_eq!(self.range.len, self.tile_seg.range.len);
        self.range.len
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// Match with a tile id only on the B side. The A side is anonymous
/// because it's typically a patch's running boundary, not a tile drawn
/// from a `TileSet`.
///
/// Asymmetric on purpose: a `PatchMatch` is not just a `TileMatch` with
/// the A-tile forgotten — patches grow over time (boundary changes)
/// while tile ids are fixed lookups into a static set.
///
/// # Asymmetric `start_offset` convention (deliberate)
///
/// The two `start_offset` fields name *different* geometric objects on
/// the two sides:
///
/// - `a_range.start_offset` is the index of the **first matched edge**
///   on A's boundary (tile-forward).
/// - `b.range.start_offset` is the index of the **first surviving edge
///   past the match** on B's boundary (tile-forward), equivalently the
///   index of the boundary *vertex* at which the matched range ends in
///   B's forward direction. This is the **seed-vertex** that
///   [`crate::geom::rat::Rat::get_match`] takes as input.
///
/// This asymmetry is not a convention choice — it reflects the geometry.
/// In an anti-parallel glue, the *same* boundary vertex sits at the
/// matched range's *start* on A and its *end* on B. The shared vertex
/// is what `get_match` and `forward_match_length` actually need as
/// input; storing the B side as the vertex coordinate means
/// `PatchMatch.b.range.start_offset` feeds those functions directly
/// without any shift, and the consumer code stays clear.
///
/// # Why not push symmetric "first matched edge" on both sides?
///
/// We tried. Storing `b.range.start_offset = first matched B edge`
/// looks tidier locally — it makes [`TileMatch::involution`] a trivial
/// swap, makes `glue_raw_angles`'s formulas look symmetric in `a` and
/// `b`, and gives a single answer to "what does `start_offset` mean?".
///
/// But the underlying geometry is *not* symmetric in storage — the
/// anti-parallel glue puts A's first matched edge and B's first matched
/// edge at *opposite ends* of the matched range. Once you store both as
/// "first matched", every consumer that wants to seed a `get_match`
/// (which structurally requires the shared vertex, = `first_matched_b +
/// len`) has to compute `+ len` at the call site. The shifts spread
/// across `~6` call sites in `patch.rs` and `matchtypes.rs`. The
/// `cw_anchor_on_central` / `ccw_anchor_on_central` pair in
/// `NeighborhoodType` (which were aliases under the asymmetric
/// convention) become genuinely different values that need careful
/// distinction at every site. The local tidiness costs more than it
/// buys.
///
/// The asymmetric storage chose its home well: the asymmetry lives in
/// **one place** (the storage convention) where it can be named and
/// documented, and consumers stay shift-free. Pushing the symmetry
/// onto storage just relocates the asymmetry into every consumer.
///
/// If a future refactor wants to revisit this: separate the bug fix
/// (the lossy `tile_offset` encoding in
/// `crate::analysis::vertextypes::bfs_phase`, see the
/// `spectre_old_encoding_*` regression) from the storage convention.
/// They are orthogonal. The tile_offset encoding was fixed without
/// touching `start_offset` semantics.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub struct PatchMatch {
    pub a_range: EdgeRange,
    pub b: Segment,
}

impl PatchMatch {
    pub fn new(a_range: EdgeRange, b: Segment) -> Self {
        debug_assert_eq!(
            a_range.len, b.range.len,
            "PatchMatch: a_range and b.range must have equal length",
        );
        Self { a_range, b }
    }

    /// Strip the B-side tile id, yielding the pure positional match.
    pub fn spec(&self) -> Match {
        Match::new(self.a_range, self.b.range)
    }

    pub fn len(&self) -> usize {
        debug_assert_eq!(self.a_range.len, self.b.range.len);
        self.a_range.len
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

// ============================================================
// Layer 3: fully bound
// ============================================================

/// Match with tile ids on both sides. Used when both endpoints are
/// tiles drawn from (one or two related) `TileSet`s.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub struct TileMatch {
    pub a: Segment,
    pub b: Segment,
}

impl TileMatch {
    pub fn new(a: Segment, b: Segment) -> Self {
        debug_assert_eq!(
            a.range.len, b.range.len,
            "TileMatch: a and b segments must have equal length",
        );
        Self { a, b }
    }

    /// Strip tile ids, yielding the pure positional match.
    pub fn spec(&self) -> Match {
        Match::new(self.a.range, self.b.range)
    }

    pub fn len(&self) -> usize {
        debug_assert_eq!(self.a.range.len, self.b.range.len);
        self.a.range.len
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The involution of this match: the same glue described from the
    /// other tile's perspective. Swaps tile ids and applies the
    /// positional involution to the edge ranges via [`Match::involution`].
    /// See that method for the (asymmetric) offset shift; see
    /// [`PatchMatch`] for why the underlying convention is asymmetric.
    pub fn involution(self, len_a: usize, len_b: usize) -> Self {
        let inv = self.spec().involution(len_a, len_b);
        TileMatch {
            a: Segment {
                tile_id: self.b.tile_id,
                range: inv.a,
            },
            b: Segment {
                tile_id: self.a.tile_id,
                range: inv.b,
            },
        }
    }
}

// ============================================================
// Tests
// ============================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn edge_range_construction() {
        let r = EdgeRange::new(3, 5);
        assert_eq!(r.start_offset, 3);
        assert_eq!(r.len, 5);
        assert!(!r.is_empty());
        assert!(EdgeRange::new(0, 0).is_empty());
    }

    #[test]
    fn match_len_matches_ranges() {
        let m = Match::new(EdgeRange::new(0, 3), EdgeRange::new(4, 3));
        assert_eq!(m.len(), 3);
        assert!(!m.is_empty());
    }

    #[test]
    fn match_empty_when_zero_len() {
        let m = Match::new(EdgeRange::new(0, 0), EdgeRange::new(0, 0));
        assert!(m.is_empty());
    }

    #[test]
    #[should_panic(expected = "equal length")]
    #[cfg(debug_assertions)]
    fn match_unequal_len_panics_in_debug() {
        let _m = Match::new(EdgeRange::new(0, 3), EdgeRange::new(4, 4));
    }

    #[test]
    fn match_involution_shifts_offsets_and_is_idempotent() {
        // len = 2, len_a = 6, len_b = 5. Asymmetric convention:
        //   a = (1, 2): matched A edges 1..3
        //   b = (3, 2): first surviving on B = 3, matched on B = 1..3 cyc
        // Involuted:
        //   new a = (3 + 5 - 2) % 5 = 1 -> (1, 2)
        //   new b = (1 + 2) % 6     = 3 -> (3, 2)
        let m = Match::new(EdgeRange::new(1, 2), EdgeRange::new(3, 2));
        let inv = m.involution(6, 5);
        assert_eq!(inv.a, EdgeRange::new(1, 2));
        assert_eq!(inv.b, EdgeRange::new(3, 2));
        // Idempotent under the swapped-length input.
        assert_eq!(inv.involution(5, 6), m);
    }

    #[test]
    fn match_with_tiles_round_trips_via_spec() {
        let m = Match::new(EdgeRange::new(1, 2), EdgeRange::new(5, 2));
        let tm = m.with_tiles(7, 9);
        assert_eq!(tm.a.tile_id, 7);
        assert_eq!(tm.b.tile_id, 9);
        assert_eq!(tm.spec(), m);
    }

    #[test]
    fn match_with_b_tile_round_trips_via_spec() {
        let m = Match::new(EdgeRange::new(1, 2), EdgeRange::new(5, 2));
        let pm = m.with_b_tile(42);
        assert_eq!(pm.b.tile_id, 42);
        assert_eq!(pm.spec(), m);
    }

    #[test]
    fn segment_basics() {
        let s = Segment::new(3, EdgeRange::new(1, 4));
        assert_eq!(s.tile_id, 3);
        assert_eq!(s.range.start_offset, 1);
        assert_eq!(s.len(), 4);
        assert!(!s.is_empty());
    }

    #[test]
    fn tile_match_spec_strips_tile_ids() {
        let tm = TileMatch::new(
            Segment::new(0, EdgeRange::new(1, 2)),
            Segment::new(1, EdgeRange::new(3, 2)),
        );
        let m = tm.spec();
        assert_eq!(m.a, EdgeRange::new(1, 2));
        assert_eq!(m.b, EdgeRange::new(3, 2));
        assert_eq!(tm.len(), 2);
    }

    #[test]
    fn tile_match_involution_swaps_tile_ids_and_shifts() {
        // Convention reminder: A.start_offset = first matched edge,
        // B.start_offset = first surviving edge just past the match.
        //
        // Take len = 2, len_a = 6, len_b = 5.
        //   a.tile_id=7, a.range = (1, 2) [matched edges 1..3 on A]
        //   b.tile_id=9, b.range = (3, 2) [first surviving = 3 on B, matched = 1..3 cyclically]
        //
        // Involuted:
        //   new a = b's segment, start shifted to "first matched on B"
        //         = (3 + 5 - 2) % 5 = 6 % 5 = 1 -> (9, EdgeRange(1, 2))
        //   new b = a's segment, start shifted to "first surviving on A"
        //         = (1 + 2) % 6 = 3 -> (7, EdgeRange(3, 2))
        let tm = TileMatch::new(
            Segment::new(7, EdgeRange::new(1, 2)),
            Segment::new(9, EdgeRange::new(3, 2)),
        );
        let inv = tm.involution(6, 5);
        assert_eq!(inv.a.tile_id, 9);
        assert_eq!(inv.a.range, EdgeRange::new(1, 2));
        assert_eq!(inv.b.tile_id, 7);
        assert_eq!(inv.b.range, EdgeRange::new(3, 2));

        // Idempotent: involution(la, lb).involution(lb, la) == original.
        assert_eq!(inv.involution(5, 6), tm);
    }

    #[test]
    fn patch_match_spec_strips_b_tile_id() {
        let pm = PatchMatch::new(EdgeRange::new(1, 2), Segment::new(5, EdgeRange::new(9, 2)));
        let m = pm.spec();
        assert_eq!(m.a, EdgeRange::new(1, 2));
        assert_eq!(m.b, EdgeRange::new(9, 2));
        assert_eq!(pm.len(), 2);
    }

    #[test]
    #[should_panic(expected = "equal length")]
    #[cfg(debug_assertions)]
    fn patch_match_unequal_len_panics() {
        let _pm = PatchMatch::new(EdgeRange::new(1, 3), Segment::new(5, EdgeRange::new(9, 2)));
    }

    #[test]
    fn ordering_is_lex_on_fields() {
        let r1 = EdgeRange::new(0, 5);
        let r2 = EdgeRange::new(1, 3);
        assert!(r1 < r2); // earlier start wins, regardless of len

        let r3 = EdgeRange::new(0, 6);
        assert!(r1 < r3); // same start, smaller len comes first
    }

    #[test]
    fn match_via_with_tiles_round_trips_through_tile_match() {
        let m = Match::new(EdgeRange::new(2, 3), EdgeRange::new(7, 3));
        let tm = m.with_tiles(11, 13);
        assert_eq!(tm.spec(), m);
    }

    #[test]
    fn serde_round_trip() {
        let tm = TileMatch::new(
            Segment::new(2, EdgeRange::new(3, 4)),
            Segment::new(5, EdgeRange::new(6, 4)),
        );
        let json = serde_json::to_string(&tm).unwrap();
        let back: TileMatch = serde_json::from_str(&json).unwrap();
        assert_eq!(back, tm);
    }
}
