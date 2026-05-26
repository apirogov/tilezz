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
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize,
)]
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
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize,
)]
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

    /// View the match from the other side: swap `a` and `b`.
    ///
    /// With the universal "`start_offset` = first matched edge"
    /// convention this is a trivial swap. Existing callers that use the
    /// asymmetric "B-side is first surviving edge" convention will need
    /// to shift starts during migration; that shifting is not this
    /// type's responsibility.
    pub fn swap_sides(self) -> Self {
        Self { a: self.b, b: self.a }
    }

    /// Attach tile ids on both sides, lifting to a [`TileMatch`].
    pub fn with_tiles(self, a_tile: usize, b_tile: usize) -> TileMatch {
        TileMatch {
            a: Segment { tile_id: a_tile, range: self.a },
            b: Segment { tile_id: b_tile, range: self.b },
        }
    }

    /// Attach a tile id only on the B side, lifting to a [`PatchMatch`].
    /// The A side stays anonymous (= the patch's own boundary).
    pub fn with_b_tile(self, b_tile: usize) -> PatchMatch {
        PatchMatch {
            a_range: self.a,
            b: Segment { tile_id: b_tile, range: self.b },
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
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize,
)]
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

/// Match with a tile id only on the B side. The A side is anonymous
/// because it's typically a patch's running boundary, not a tile drawn
/// from a `TileSet`.
///
/// Asymmetric on purpose: a `PatchMatch` is not just a `TileMatch` with
/// the A-tile forgotten — patches grow over time (boundary changes)
/// while tile ids are fixed lookups into a static set.
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize,
)]
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
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize,
)]
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

    /// View the match from the other side: swap A and B segments.
    /// See [`Match::swap_sides`] for the convention note.
    pub fn swap_sides(self) -> Self {
        Self { a: self.b, b: self.a }
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
    fn match_swap_sides_is_involutive() {
        let m = Match::new(EdgeRange::new(1, 2), EdgeRange::new(5, 2));
        assert_eq!(m.swap_sides().swap_sides(), m);
        assert_eq!(m.swap_sides().a, m.b);
        assert_eq!(m.swap_sides().b, m.a);
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
    fn tile_match_swap_sides_swaps_tile_ids_too() {
        let tm = TileMatch::new(
            Segment::new(7, EdgeRange::new(1, 2)),
            Segment::new(9, EdgeRange::new(3, 2)),
        );
        let swapped = tm.swap_sides();
        assert_eq!(swapped.a.tile_id, 9);
        assert_eq!(swapped.b.tile_id, 7);
        assert_eq!(swapped.swap_sides(), tm);
    }

    #[test]
    fn patch_match_spec_strips_b_tile_id() {
        let pm = PatchMatch::new(
            EdgeRange::new(1, 2),
            Segment::new(5, EdgeRange::new(9, 2)),
        );
        let m = pm.spec();
        assert_eq!(m.a, EdgeRange::new(1, 2));
        assert_eq!(m.b, EdgeRange::new(9, 2));
        assert_eq!(pm.len(), 2);
    }

    #[test]
    #[should_panic(expected = "equal length")]
    #[cfg(debug_assertions)]
    fn patch_match_unequal_len_panics() {
        let _pm = PatchMatch::new(
            EdgeRange::new(1, 3),
            Segment::new(5, EdgeRange::new(9, 2)),
        );
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
        // and swapping sides commutes with spec()
        assert_eq!(tm.swap_sides().spec(), m.swap_sides());
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
