//! Segment-anchored BFS: enumerate every reachable open-segment type
//! by growing patches one tile at a time.
//!
//! Unlike the central-tile-anchored BFS in
//! [`crate::intgeom::neighborhood`], this BFS doesn't fix a central
//! and grow a corona — it starts from every 2-tile glue ("seed"
//! patch with two junctions), and from each frontier segment tries
//! every legal tile attachment that touches the segment's extended
//! range. Junctions discovered on the resulting boundary become new
//! segment types. By BFS order, the first patch to produce a given
//! seg type is the minimal witness.
//!
//! Goal: produce a seg-type index that is **closed under one further
//! glue** (= the gap diagnosed in `neighborhood`'s module docs).
//! Rule extraction (LHS/RHS pairs) is intentionally out of scope —
//! per-segment classification is conceptually unsound, and the right
//! abstraction is segment-sequence-anchored, deferred to a later
//! phase.

use std::collections::VecDeque;
use std::sync::Arc;

use rustc_hash::FxHashMap;

use crate::cyclotomic::{IsComplex, IsRingOrField, Units};
use crate::intgeom::patch::{EdgeInfo, GrowingPatch, OpenVertexType};
use crate::intgeom::segtypes::OpenSegmentTypeIndex;
use crate::intgeom::tileset::TileSet;

/// Canonical key for a `GrowingPatch` after `normalize()`.
type PatchKey = (Vec<i8>, Vec<EdgeInfo>, Vec<Vec<EdgeInfo>>, Vec<usize>);

/// One witness: a reference into `patches` plus the boundary
/// positions of the two junctions that define this segment on that
/// patch. The first time a seg is seen, the witness is recorded; BFS
/// order guarantees it's minimal.
#[derive(Clone, Debug)]
pub struct SegWitness {
    pub patch_id: usize,
    pub cw_junc_pos: usize,
    pub ccw_junc_pos: usize,
}

#[derive(Debug)]
pub enum BfsError {
    /// Seg count exceeded the safety cap. Reported with the cap.
    SegCapExceeded { cap: usize },
    /// Patch count exceeded the safety cap.
    PatchCapExceeded { cap: usize },
}

pub struct SegmentTypeBFS<T: IsComplex> {
    tileset: Arc<TileSet<T>>,

    // Deduped witness patches. patch_lookup maps PatchKey to index.
    patches: Vec<GrowingPatch<T>>,
    patch_lookup: FxHashMap<PatchKey, usize>,

    // Seg-type index (grows incrementally as new segs are discovered).
    seg_index: OpenSegmentTypeIndex,

    // Witness per seg id: matches the length of seg_index.
    seg_witnesses: Vec<SegWitness>,

    queue: VecDeque<usize>,

    // Safety caps; if either is exceeded, run returns an error.
    seg_cap: usize,
    patch_cap: usize,
}

impl<T: IsComplex + IsRingOrField + Units> SegmentTypeBFS<T> {
    /// Run the segment-anchored BFS to completion. `seg_cap` and
    /// `patch_cap` are safety bounds — if either is exceeded, the
    /// run aborts with a `BfsError`.
    pub fn run(
        tileset: Arc<TileSet<T>>,
        seg_cap: usize,
        patch_cap: usize,
    ) -> Result<Self, BfsError> {
        let mut bfs = Self {
            tileset,
            patches: Vec::new(),
            patch_lookup: FxHashMap::default(),
            seg_index: OpenSegmentTypeIndex::new_empty(),
            seg_witnesses: Vec::new(),
            queue: VecDeque::new(),
            seg_cap,
            patch_cap,
        };
        bfs.seed()?;
        bfs.run_bfs()?;
        Ok(bfs)
    }

    pub fn seg_index(&self) -> &OpenSegmentTypeIndex {
        &self.seg_index
    }

    pub fn patches(&self) -> &[GrowingPatch<T>] {
        &self.patches
    }

    pub fn witness_of(&self, seg_id: usize) -> Option<&SegWitness> {
        self.seg_witnesses.get(seg_id)
    }

    pub fn num_segs(&self) -> usize {
        self.seg_index.len()
    }

    pub fn num_patches(&self) -> usize {
        self.patches.len()
    }

    /// Seed: every 2-tile glue produces 2 segments (= cyclic pair on
    /// the seed patch's 2-junction boundary). Each is registered with
    /// the seed patch as witness.
    fn seed(&mut self) -> Result<(), BfsError> {
        let n_tiles = self.tileset.num_tiles();
        for tile_a in 0..n_tiles {
            let base = GrowingPatch::new(Arc::clone(&self.tileset), tile_a);
            let matches = base.get_all_matches();
            for pm in matches {
                let mut patch = base.clone();
                if !patch.add_tile(&pm) {
                    continue;
                }
                patch.normalize();
                self.harvest_segs_from(patch)?;
            }
        }
        Ok(())
    }

    fn run_bfs(&mut self) -> Result<(), BfsError> {
        while let Some(seg_id) = self.queue.pop_front() {
            self.expand_seg(seg_id)?;
        }
        Ok(())
    }

    /// Expand one frontier seg: enumerate every legal tile attachment
    /// that touches the seg's extended range, apply each, and harvest
    /// the resulting boundary's seg pairs.
    fn expand_seg(&mut self, seg_id: usize) -> Result<(), BfsError> {
        let witness = self.seg_witnesses[seg_id].clone();
        let patch = self.patches[witness.patch_id].clone();
        let n = patch.boundary_len();
        if n == 0 {
            return Ok(());
        }
        let cw_junc = witness.cw_junc_pos;
        let ccw_junc = witness.ccw_junc_pos;
        // Extended range: 1 edge of CW-adjacent seg + own seg edges +
        // 1 edge of CCW-adjacent seg.
        let ext_start = (cw_junc + n - 1) % n;
        let ext_end = ccw_junc;
        let candidates = patch.get_matches_in_edge_range(ext_start, ext_end);
        for pm in candidates {
            let mut trial = patch.clone();
            if !trial.add_tile(&pm) {
                continue;
            }
            trial.normalize();
            self.harvest_segs_from(trial)?;
        }
        Ok(())
    }

    /// Intern the (deduped) patch, enumerate its boundary junction
    /// pairs cyclically, and register any new seg types (each with
    /// the just-interned patch as witness and the corresponding
    /// junction positions).
    fn harvest_segs_from(&mut self, patch: GrowingPatch<T>) -> Result<(), BfsError> {
        let patch_id = self.intern_patch(patch)?;
        let patch_ref = &self.patches[patch_id];
        let n = patch_ref.boundary_len();
        let juncs: Vec<(usize, OpenVertexType)> = (0..n)
            .filter_map(|i| patch_ref.junction_vertex_type_at(i).map(|vt| (i, vt)))
            .collect();
        let k = juncs.len();
        if k < 2 {
            return Ok(());
        }
        for j in 0..k {
            let (cw_pos, cw_vt) = (juncs[j].0, juncs[j].1.clone());
            let (ccw_pos, ccw_vt) = (juncs[(j + 1) % k].0, juncs[(j + 1) % k].1.clone());
            let before = self.seg_index.len();
            let seg_id = self.seg_index.push_pair(cw_vt, ccw_vt);
            if seg_id == before {
                // Newly registered: record witness, enqueue.
                self.seg_witnesses.push(SegWitness {
                    patch_id,
                    cw_junc_pos: cw_pos,
                    ccw_junc_pos: ccw_pos,
                });
                self.queue.push_back(seg_id);
                if self.seg_index.len() > self.seg_cap {
                    return Err(BfsError::SegCapExceeded { cap: self.seg_cap });
                }
            }
        }
        Ok(())
    }

    fn intern_patch(&mut self, patch: GrowingPatch<T>) -> Result<usize, BfsError> {
        let key = patch_key(&patch);
        if let Some(&id) = self.patch_lookup.get(&key) {
            return Ok(id);
        }
        let id = self.patches.len();
        self.patch_lookup.insert(key, id);
        self.patches.push(patch);
        if self.patches.len() > self.patch_cap {
            return Err(BfsError::PatchCapExceeded { cap: self.patch_cap });
        }
        Ok(id)
    }
}

fn patch_key<T: IsComplex + IsRingOrField + Units>(patch: &GrowingPatch<T>) -> PatchKey {
    (
        patch.angles().to_vec(),
        patch.edges().to_vec(),
        patch.inner_chains().to_vec(),
        patch.patch_tile_ids().to_vec(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic::ZZ12;
    use crate::intgeom::rat::Rat;
    use crate::intgeom::tiles;

    fn hex_tileset() -> Arc<TileSet<ZZ12>> {
        Arc::new(TileSet::new(vec![
            Rat::try_from(&tiles::hexagon::<ZZ12>()).unwrap(),
        ]))
    }

    fn spectre_tileset() -> Arc<TileSet<ZZ12>> {
        Arc::new(TileSet::new(vec![
            Rat::try_from(&tiles::spectre::<ZZ12>()).unwrap(),
        ]))
    }

    #[test]
    fn hex_segs_enumerated() {
        let bfs = SegmentTypeBFS::run(hex_tileset(), 100_000, 100_000)
            .expect("hex BFS terminates within caps");
        eprintln!(
            "hex: segs={} patches={}",
            bfs.num_segs(),
            bfs.num_patches()
        );
        assert!(bfs.num_segs() > 0);
    }

    /// Spectre is much larger than hex; gated for cost.
    #[test]
    #[ignore]
    fn spectre_segs_enumerated() {
        let bfs = SegmentTypeBFS::run(spectre_tileset(), 500_000, 500_000)
            .expect("spectre BFS terminates within caps");
        eprintln!(
            "spectre: segs={} patches={}",
            bfs.num_segs(),
            bfs.num_patches()
        );
        assert!(bfs.num_segs() > 0);
    }

    /// Closure: starting from the BFS-produced seg index, take every
    /// witness patch and try every legal further glue — no new segs
    /// should arise.
    fn assert_seg_set_closed<T: IsComplex + IsRingOrField + Units>(bfs: &SegmentTypeBFS<T>) {
        let mut new_segs = 0usize;
        for patch in bfs.patches() {
            for pm in patch.get_all_matches() {
                let mut trial = patch.clone();
                if !trial.add_tile(&pm) {
                    continue;
                }
                trial.normalize();
                let n = trial.boundary_len();
                let juncs: Vec<OpenVertexType> = (0..n)
                    .filter_map(|i| trial.junction_vertex_type_at(i))
                    .collect();
                let k = juncs.len();
                if k < 2 {
                    continue;
                }
                for j in 0..k {
                    let cw = &juncs[j];
                    let ccw = &juncs[(j + 1) % k];
                    if bfs.seg_index().vertex_types_to_seg(cw, ccw).is_none() {
                        new_segs += 1;
                    }
                }
            }
        }
        assert_eq!(new_segs, 0, "BFS seg set should be closed under one more glue");
    }

    #[test]
    fn hex_seg_set_is_closed() {
        let bfs = SegmentTypeBFS::run(hex_tileset(), 100_000, 100_000).unwrap();
        assert_seg_set_closed(&bfs);
    }

    /// Spectre closure: O(patches × glues_per_patch) — large; gated.
    #[test]
    #[ignore]
    fn spectre_seg_set_is_closed() {
        let bfs = SegmentTypeBFS::run(spectre_tileset(), 500_000, 500_000).unwrap();
        assert_seg_set_closed(&bfs);
    }
}
