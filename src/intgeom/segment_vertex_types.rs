//! **Segment vertex types**: the cyclic sequence of boundary
//! [`Segment`]s meeting at one central-tile vertex of a closed
//! corona ([`SurroundedTile`](crate::intgeom::neighborhood::SurroundedTile)
//! with `is_closed == true`).
//!
//! This is a richer counterpart to
//! [`OpenVertexType`](crate::intgeom::patch::OpenVertexType) and
//! [`ClosedVertexType`](crate::intgeom::patch::ClosedVertexType):
//! where those types list incident *edges* (each carrying an
//! `EdgeInfo` = `(tile_id, tile_offset)`), a `SegmentVertexType`
//! lists incident *segments* — runs of consecutive boundary edges on
//! a tile, bounded by junctions. Two segments of the same tile shape
//! at the same `tile_offset` and `length` are considered the same
//! segment intrinsically; instance information (i.e. which copy of
//! the tile in the patch) is not retained.
//!
//! # Cyclic-vec convention
//!
//! At a central-tile vertex with `K + 1` incident tile instances
//! (the central + `K` context tiles), there are `2(K + 1)` segment
//! endpoints meeting (each tile contributes its CW-side and CCW-side
//! segments at this corner). The cyclic vec lists them in CCW order
//! around the vertex; pairs of adjacent entries that belong to the
//! same tile are the tile's two segments meeting at the vertex.
//!
//! # Status (phase A+B plumbing only)
//!
//! This module currently provides the data types and an empty index
//! shell. The BFS-replay-based extraction that populates the catalog
//! lives in phases C and D (not yet implemented). For now,
//! [`SegmentVertexTypeIndex::new_empty`] is the only constructor.

use std::collections::HashMap;
use std::sync::Arc;

use crate::cyclotomic::{IsComplex, IsRingOrField, Units};
use crate::intgeom::neighborhood::NeighborhoodIndex;
use crate::intgeom::patch::canonical_cyclic_rotation;
use crate::intgeom::tileset::TileSet;

/// One contiguous segment on a tile's boundary.
///
/// `tile_id` identifies the tile shape (index into the parent
/// [`TileSet`]); `tile_offset` is the CCW-starting edge offset of the
/// segment within that tile (in `0..tile.len()`); `length` is the
/// segment length measured in **edges** (so a length-3 segment
/// starting at `tile_offset = 2` covers edges 2, 3, 4 of the tile).
///
/// `length` is bounded by the tile's edge count; `tile_offset +
/// length` may wrap modulo the tile length for segments that
/// straddle the tile's edge-0 seam.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Segment {
    pub tile_id: usize,
    pub tile_offset: usize,
    pub length: usize,
}

/// The cyclic CCW sequence of [`Segment`]s meeting at one central-
/// tile vertex of a closed corona. Stored in lex-min cyclic
/// rotation, so two vectors that are cyclic rotations of one another
/// compare equal.
///
/// See the module-level docs for the cyclic-vec convention. The
/// length of the vec is `2 · (number of incident tile instances at
/// the vertex)`.
#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct SegmentVertexType {
    petals: Vec<Segment>,
}

impl SegmentVertexType {
    /// Build from a raw cyclic petal sequence; canonicalises to
    /// lex-min cyclic rotation.
    pub fn from_cyclic(petals: &[Segment]) -> Self {
        SegmentVertexType {
            petals: canonical_cyclic_rotation(petals),
        }
    }

    /// Build from a sequence asserted to already be in canonical
    /// (lex-min cyclic rotation) form. Use [`Self::from_cyclic`] for
    /// raw input.
    #[cfg(test)]
    pub(crate) fn from_canonical_unchecked(petals: Vec<Segment>) -> Self {
        SegmentVertexType { petals }
    }

    /// The petals, in CCW order starting at the lex-min rotation.
    pub fn petals(&self) -> &[Segment] {
        &self.petals
    }

    pub fn len(&self) -> usize {
        self.petals.len()
    }

    pub fn is_empty(&self) -> bool {
        self.petals.is_empty()
    }
}

/// Per-central-vertex catalog of segment vertex types across every
/// closed corona reachable from a [`NeighborhoodIndex`].
///
/// **Status:** the data-type and predecessor-plumbing scaffolding is
/// in place; the BFS-replay-based extraction is not yet implemented.
/// [`Self::new_empty`] is the only working constructor today.
pub struct SegmentVertexTypeIndex<T: IsComplex> {
    tileset: Arc<TileSet<T>>,
    /// All distinct `SegmentVertexType`s discovered, in canonical
    /// (sorted) order. Entry at index `id - 1` has 1-based id `id`.
    entries: Vec<SegmentVertexType>,
    /// Reverse map: `SegmentVertexType -> 1-based id`.
    reverse: HashMap<SegmentVertexType, usize>,
    /// Per closed-SurroundedTile entry id (in the parent
    /// `NeighborhoodIndex`), the `SegmentVertexTypeIndex` ids of its
    /// `central_tile.len()` central-vertex segment configurations,
    /// in CCW order around the central. Empty when the corresponding
    /// `NtEntry` is not a closed phase-2 entry.
    closed_corona_vertex_ids: HashMap<usize, Vec<usize>>,
}

impl<T: IsComplex + IsRingOrField + Units> SegmentVertexTypeIndex<T> {
    /// Construct an empty catalog over the same tileset as `idx`.
    /// Phase-A placeholder: real entries will be populated by the
    /// not-yet-implemented BFS-replay extraction.
    pub fn new_empty(idx: &NeighborhoodIndex<T>) -> Self {
        SegmentVertexTypeIndex {
            tileset: Arc::clone(idx.tileset()),
            entries: Vec::new(),
            reverse: HashMap::new(),
            closed_corona_vertex_ids: HashMap::new(),
        }
    }

    pub fn tileset(&self) -> &Arc<TileSet<T>> {
        &self.tileset
    }

    /// Number of distinct segment vertex types in the catalog.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// All catalog entries in canonical (sorted) order. Entry at
    /// index `id - 1` has id `id`.
    pub fn entries(&self) -> &[SegmentVertexType] {
        &self.entries
    }

    /// 1-based id of `svt` in this catalog, or `None` if absent.
    pub fn get_id(&self, svt: &SegmentVertexType) -> Option<usize> {
        self.reverse.get(svt).copied()
    }

    /// Return the catalog entry for id `id`. Panics if out of range.
    pub fn get(&self, id: usize) -> &SegmentVertexType {
        assert!(id >= 1 && id <= self.entries.len(), "id out of range");
        &self.entries[id - 1]
    }

    /// For closed-SurroundedTile entry id `st_entry_id` (in the
    /// parent `NeighborhoodIndex`), return the per-central-vertex
    /// catalog ids (in CCW order around the central), or `None` if
    /// no extraction has been recorded for that entry.
    pub fn vertex_ids_for_closed_corona(&self, st_entry_id: usize) -> Option<&[usize]> {
        self.closed_corona_vertex_ids
            .get(&st_entry_id)
            .map(|v| v.as_slice())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seg(tile_id: usize, tile_offset: usize, length: usize) -> Segment {
        Segment {
            tile_id,
            tile_offset,
            length,
        }
    }

    #[test]
    fn segment_vertex_type_canonicalises_to_lex_min_rotation() {
        let petals = [seg(2, 0, 1), seg(0, 1, 3), seg(1, 2, 2), seg(0, 3, 1)];
        let canonical = SegmentVertexType::from_cyclic(&petals);
        for shift in 0..petals.len() {
            let rotated: Vec<Segment> =
                (0..petals.len()).map(|i| petals[(shift + i) % petals.len()]).collect();
            assert_eq!(
                SegmentVertexType::from_cyclic(&rotated),
                canonical,
                "rotation by {shift} should canonicalise to the same SVT"
            );
        }
        let petals = canonical.petals();
        for k in 1..petals.len() {
            assert!(petals[0] <= petals[k]);
        }
    }

    #[test]
    fn segment_vertex_type_distinguishes_non_rotation_orderings() {
        let a = SegmentVertexType::from_cyclic(&[seg(0, 0, 1), seg(0, 1, 1), seg(0, 2, 1)]);
        let b = SegmentVertexType::from_cyclic(&[seg(0, 0, 1), seg(0, 2, 1), seg(0, 1, 1)]);
        assert_ne!(a, b);
    }

    #[test]
    fn segment_vertex_type_empty_and_singleton() {
        let empty = SegmentVertexType::from_cyclic(&[]);
        assert!(empty.is_empty());
        let singleton = SegmentVertexType::from_cyclic(&[seg(1, 2, 3)]);
        assert_eq!(singleton.len(), 1);
        assert_eq!(singleton.petals()[0], seg(1, 2, 3));
    }

    #[test]
    fn from_canonical_unchecked_round_trip() {
        let raw = vec![seg(0, 0, 1), seg(1, 0, 1), seg(2, 0, 1)];
        let direct = SegmentVertexType::from_canonical_unchecked(raw.clone());
        let via_cyclic = SegmentVertexType::from_cyclic(&raw);
        assert_eq!(direct, via_cyclic);
    }
}
