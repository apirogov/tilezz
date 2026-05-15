//! Open segment types and the BFS-derived segment type index.
//!
//! A **segment** on a patch open boundary is a maximal run of edges all
//! contributed by a single tile instance, bounded at each end by a
//! [`OpenVertexType`] junction. A **segment type** is the equivalence
//! class of such runs, keyed by the ordered pair of bounding VTs
//! `(vt_cw, vt_ccw)`. Two segments share a type iff their boundary
//! geometry is identical — same tile, same exposed tile offsets, same
//! neighboring tile arrangement at both endpoints.
//!
//! [`OpenSegmentTypeIndex`] is the registry of all segment types that
//! appear on the open boundary of any patch reachable by the BFS in
//! [`crate::intgeom::neighborhood::NeighborhoodIndex`]. It is built
//! once from a `NeighborhoodIndex` by enumerating every consecutive
//! junction pair on every reachable phase-1 ctx boundary, dedup'ing,
//! and assigning stable 0-based ids in canonical sort order.
//!
//! Phase-2 closed coronas are terminal in the BFS and never serve as a
//! ctx, so they do not contribute new segment types.

use std::sync::Arc;

use rustc_hash::FxHashMap;

use crate::cyclotomic::{IsComplex, IsRingOrField, Units};
use crate::intgeom::matchtypes::MatchTypeIndex;
use crate::intgeom::neighborhood::{NeighborhoodIndex, NtEntry};
use crate::intgeom::patch::{GrowingPatch, OpenVertexType};

/// An ordered pair of open vertex types defining one segment type.
/// `cw` is the segment's CW endpoint junction (the one reached by
/// walking CW into the segment from the interior); `ccw` is the CCW
/// endpoint. The exposed tile-offset run between the two endpoints is
/// uniquely determined by the pair.
pub type SegmentTypePair = (OpenVertexType, OpenVertexType);

/// Registry of every open segment type reachable on the open boundary
/// of any phase-1 BFS state in the source `NeighborhoodIndex`.
///
/// Ids are stable: built in canonical sort order over the dedup'd pair
/// set, 0-based.
#[derive(Clone, Debug)]
pub struct OpenSegmentTypeIndex {
    segments: Vec<SegmentTypePair>,
    by_pair: FxHashMap<SegmentTypePair, usize>,
}

impl OpenSegmentTypeIndex {
    /// Build the index from a `NeighborhoodIndex` by reconstructing
    /// every phase-1 entry's ctx and recording every consecutive
    /// junction pair on its full open boundary.
    pub fn build_from_bfs<T: IsComplex + IsRingOrField + Units>(
        idx: &NeighborhoodIndex<T>,
    ) -> Self {
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        let mut seen: FxHashMap<SegmentTypePair, ()> = FxHashMap::default();

        for entry in idx.entries() {
            let NtEntry::Phase1(nt) = entry else {
                continue;
            };
            let Some((ctx, _)) = GrowingPatch::construct_witness_from_vt_sequence(
                &nt.vt_seq,
                Arc::clone(&mi),
            ) else {
                continue;
            };
            let n = ctx.boundary_len();
            // Collect all junctions in CCW order with their positions.
            let juncs: Vec<(usize, OpenVertexType)> = (0..n)
                .filter_map(|i| ctx.junction_vertex_type_at(i).map(|vt| (i, vt)))
                .collect();
            // Cyclic boundary: consecutive pairs include the wrap from
            // the last junction back to the first.
            let k = juncs.len();
            if k < 2 {
                continue;
            }
            for j in 0..k {
                let cw = juncs[j].1.clone();
                let ccw = juncs[(j + 1) % k].1.clone();
                seen.insert((cw, ccw), ());
            }
        }

        let mut segments: Vec<SegmentTypePair> = seen.into_keys().collect();
        segments.sort();
        let by_pair = segments
            .iter()
            .enumerate()
            .map(|(i, p)| (p.clone(), i))
            .collect();
        OpenSegmentTypeIndex { segments, by_pair }
    }

    /// Look up the segment id for a `(cw_vt, ccw_vt)` pair, if known.
    pub fn vertex_types_to_seg(
        &self,
        cw: &OpenVertexType,
        ccw: &OpenVertexType,
    ) -> Option<usize> {
        self.by_pair.get(&(cw.clone(), ccw.clone())).copied()
    }

    /// Look up the VT pair for a segment id, if in range.
    pub fn seg_to_vertex_types(&self, id: usize) -> Option<&SegmentTypePair> {
        self.segments.get(id)
    }

    pub fn len(&self) -> usize {
        self.segments.len()
    }

    pub fn is_empty(&self) -> bool {
        self.segments.is_empty()
    }

    pub fn iter(&self) -> impl Iterator<Item = (usize, &SegmentTypePair)> {
        self.segments.iter().enumerate()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic::ZZ12;
    use crate::intgeom::neighborhood::NeighborhoodIndex;
    use crate::intgeom::rat::Rat;
    use crate::intgeom::tiles;
    use crate::intgeom::tileset::TileSet;
    use std::sync::OnceLock;

    fn hex_idx() -> &'static NeighborhoodIndex<ZZ12> {
        static IDX: OnceLock<NeighborhoodIndex<ZZ12>> = OnceLock::new();
        IDX.get_or_init(|| {
            let rat = Rat::try_from(&tiles::hexagon::<ZZ12>()).unwrap();
            NeighborhoodIndex::new(Arc::new(TileSet::new(vec![rat])))
        })
    }

    fn spectre_idx() -> &'static NeighborhoodIndex<ZZ12> {
        static IDX: OnceLock<NeighborhoodIndex<ZZ12>> = OnceLock::new();
        IDX.get_or_init(|| {
            let rat = Rat::try_from(&tiles::spectre::<ZZ12>()).unwrap();
            NeighborhoodIndex::new(Arc::new(TileSet::new(vec![rat])))
        })
    }

    #[test]
    fn seg_index_roundtrips_hex() {
        let seg_idx = OpenSegmentTypeIndex::build_from_bfs(hex_idx());
        assert!(seg_idx.len() > 0, "expected at least one segment type");
        for (id, pair) in seg_idx.iter() {
            let (cw, ccw) = pair;
            let looked_up = seg_idx
                .vertex_types_to_seg(cw, ccw)
                .expect("pair should be in index");
            assert_eq!(looked_up, id, "id roundtrip failed for {:?}", pair);
        }
    }

    #[test]
    fn seg_index_roundtrips_spectre() {
        let seg_idx = OpenSegmentTypeIndex::build_from_bfs(spectre_idx());
        assert!(seg_idx.len() > 0, "expected at least one segment type");
        for (id, pair) in seg_idx.iter() {
            let (cw, ccw) = pair;
            assert_eq!(
                seg_idx.vertex_types_to_seg(cw, ccw),
                Some(id),
                "id roundtrip failed for id {}",
                id
            );
        }
    }

    #[test]
    fn seg_counts_printed() {
        // Data-only sanity log. Doesn't pin counts (those will drift
        // as BFS structure changes); just records orders of magnitude.
        for (name, idx) in [
            ("hex", hex_idx()),
            ("spectre", spectre_idx()),
        ] {
            let seg_idx = OpenSegmentTypeIndex::build_from_bfs(idx);
            eprintln!("{name}: {} distinct segment types", seg_idx.len());
        }
    }
}
