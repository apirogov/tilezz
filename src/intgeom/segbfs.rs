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
//! Segment types are keyed by pairs of [`CoarseJunction`]s — the
//! boundary-only equivalence (= `cw_edge`, `ccw_edge`, boundary
//! `angle`). The interior-tile arrangement at a junction doesn't
//! affect outgoing tile attachments, so coarser equivalence groups
//! configurations that behave identically for further growth.
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
use crate::intgeom::patch::{CoarseJunction, EdgeInfo, GrowingPatch};
use crate::intgeom::tileset::TileSet;

/// Pair of coarse junctions defining one segment type: the segment
/// stretches CCW from `cw` to `ccw`, contributed by a single tile
/// instance.
pub type SegPair = (CoarseJunction, CoarseJunction);

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

    // Inline seg-type index: SegPair → seg id. Ids are assigned in
    // BFS discovery order. The reverse mapping `seg_pairs[id] -> pair`
    // is kept for lookup-by-id.
    seg_lookup: FxHashMap<SegPair, usize>,
    seg_pairs: Vec<SegPair>,

    // Witness per seg id: matches the length of seg_pairs.
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
            seg_lookup: FxHashMap::default(),
            seg_pairs: Vec::new(),
            seg_witnesses: Vec::new(),
            queue: VecDeque::new(),
            seg_cap,
            patch_cap,
        };
        bfs.seed()?;
        bfs.run_bfs()?;
        Ok(bfs)
    }

    pub fn seg_pairs(&self) -> &[SegPair] {
        &self.seg_pairs
    }

    pub fn lookup_seg(&self, cw: &CoarseJunction, ccw: &CoarseJunction) -> Option<usize> {
        self.seg_lookup.get(&(*cw, *ccw)).copied()
    }

    pub fn patches(&self) -> &[GrowingPatch<T>] {
        &self.patches
    }

    pub fn witness_of(&self, seg_id: usize) -> Option<&SegWitness> {
        self.seg_witnesses.get(seg_id)
    }

    pub fn num_segs(&self) -> usize {
        self.seg_pairs.len()
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
    /// whose absorbed run **overlaps** the seg's extended range (=
    /// seg edges + 1 edge of pad on each side). The witness patch
    /// validates each match's full edge sequence — matches that
    /// extend past the extended range are still legal so long as the
    /// witness's edges accommodate them. Apply each, harvest the
    /// resulting boundary's coarse junction pairs.
    fn expand_seg(&mut self, seg_id: usize) -> Result<(), BfsError> {
        let witness = self.seg_witnesses[seg_id].clone();
        let patch = self.patches[witness.patch_id].clone();
        let n = patch.boundary_len();
        if n == 0 {
            return Ok(());
        }
        let cw_junc = witness.cw_junc_pos;
        let ccw_junc = witness.ccw_junc_pos;
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

    /// Enumerate the patch's boundary coarse-junction pairs, and if
    /// **any** is new, intern the patch and register the new seg
    /// types (each with this patch as witness). If no pair is new,
    /// the patch isn't interned at all — it contributes nothing the
    /// BFS needs.
    fn harvest_segs_from(&mut self, patch: GrowingPatch<T>) -> Result<(), BfsError> {
        let n = patch.boundary_len();
        let juncs: Vec<(usize, CoarseJunction)> = (0..n)
            .filter_map(|i| patch.coarse_junction_at(i).map(|cj| (i, cj)))
            .collect();
        let k = juncs.len();
        if k < 2 {
            return Ok(());
        }
        // Pre-scan: do any of the junction pairs introduce a new seg?
        // If not, skip this patch entirely.
        let any_new = (0..k).any(|j| {
            let pair = (juncs[j].1, juncs[(j + 1) % k].1);
            !self.seg_lookup.contains_key(&pair)
        });
        if !any_new {
            return Ok(());
        }
        let patch_id = self.intern_patch(patch)?;
        for j in 0..k {
            let (cw_pos, cw_cj) = juncs[j];
            let (ccw_pos, ccw_cj) = juncs[(j + 1) % k];
            let pair = (cw_cj, ccw_cj);
            if self.seg_lookup.contains_key(&pair) {
                continue;
            }
            let seg_id = self.seg_pairs.len();
            self.seg_lookup.insert(pair, seg_id);
            self.seg_pairs.push(pair);
            self.seg_witnesses.push(SegWitness {
                patch_id,
                cw_junc_pos: cw_pos,
                ccw_junc_pos: ccw_pos,
            });
            self.queue.push_back(seg_id);
            if self.seg_pairs.len() > self.seg_cap {
                return Err(BfsError::SegCapExceeded { cap: self.seg_cap });
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
    /// witness patch and try every legal further glue — no new
    /// `CoarseJunction` pairs should arise. Skips `normalize()` on
    /// the trial because `CoarseJunction` is normalize-invariant.
    /// Logs progress every `progress_every` patches.
    fn assert_seg_set_closed<T: IsComplex + IsRingOrField + Units>(
        bfs: &SegmentTypeBFS<T>,
        name: &str,
        progress_every: usize,
    ) {
        let total_patches = bfs.patches().len();
        let start = std::time::Instant::now();
        let mut glues_tried = 0usize;
        let mut glues_succeeded = 0usize;
        let mut new_segs = 0usize;
        for (idx, patch) in bfs.patches().iter().enumerate() {
            for pm in patch.get_all_matches() {
                glues_tried += 1;
                let mut trial = patch.clone();
                if !trial.add_tile(&pm) {
                    continue;
                }
                glues_succeeded += 1;
                let n = trial.boundary_len();
                let juncs: Vec<CoarseJunction> = (0..n)
                    .filter_map(|i| trial.coarse_junction_at(i))
                    .collect();
                let k = juncs.len();
                if k < 2 {
                    continue;
                }
                for j in 0..k {
                    let cw = &juncs[j];
                    let ccw = &juncs[(j + 1) % k];
                    if bfs.lookup_seg(cw, ccw).is_none() {
                        new_segs += 1;
                    }
                }
            }
            if progress_every > 0 && (idx + 1) % progress_every == 0 {
                eprintln!(
                    "{name}: closure check patch {}/{} ({:.1}%), glues tried={}, succeeded={}, new segs so far={}, elapsed={:.1}s",
                    idx + 1,
                    total_patches,
                    100.0 * (idx + 1) as f64 / total_patches as f64,
                    glues_tried,
                    glues_succeeded,
                    new_segs,
                    start.elapsed().as_secs_f64(),
                );
            }
        }
        eprintln!(
            "{name}: closure check FINAL: patches={}, glues tried={}, succeeded={}, new segs={}, elapsed={:.1}s",
            total_patches,
            glues_tried,
            glues_succeeded,
            new_segs,
            start.elapsed().as_secs_f64(),
        );
        assert_eq!(new_segs, 0, "{name} BFS seg set should be closed under one more glue");
    }

    #[test]
    fn hex_seg_set_is_closed() {
        let bfs = SegmentTypeBFS::run(hex_tileset(), 100_000, 100_000).unwrap();
        assert_seg_set_closed(&bfs, "hex", 0);
    }

    /// Trace through one missing seg on spectre: find a witness patch
    /// P and a glue M that produces a `CoarseJunction` pair NOT in
    /// the BFS's seg index, then report the geometric details for
    /// analysis.
    #[test]
    #[ignore]
    fn spectre_diagnose_missing_seg() {
        let bfs = SegmentTypeBFS::run(spectre_tileset(), 500_000, 500_000).unwrap();
        eprintln!(
            "spectre BFS: {} segs, {} witness patches",
            bfs.num_segs(),
            bfs.num_patches()
        );
        for (patch_idx, patch) in bfs.patches().iter().enumerate() {
            let patch_n = patch.boundary_len();
            for pm in patch.get_all_matches() {
                let mut trial = patch.clone();
                if !trial.add_tile(&pm) {
                    continue;
                }
                let n = trial.boundary_len();
                let juncs: Vec<(usize, CoarseJunction)> = (0..n)
                    .filter_map(|i| trial.coarse_junction_at(i).map(|cj| (i, cj)))
                    .collect();
                let k = juncs.len();
                if k < 2 {
                    continue;
                }
                for j in 0..k {
                    let (cw_pos, cw_cj) = &juncs[j];
                    let (ccw_pos, ccw_cj) = &juncs[(j + 1) % k];
                    if bfs.lookup_seg(cw_cj, ccw_cj).is_none() {
                        eprintln!("---");
                        eprintln!("MISSING SEG FOUND on patch {}:", patch_idx + 1);
                        eprintln!("  patch boundary_len = {}", patch_n);
                        eprintln!(
                            "  applied match: start_a={}, len={}, start_b={}, tile_id={}",
                            pm.start_a, pm.len, pm.start_b, pm.tile_id
                        );
                        eprintln!("  trial boundary_len = {}", n);
                        eprintln!(
                            "  missing pair position on trial: cw vertex {} → ccw vertex {}",
                            cw_pos, ccw_pos
                        );
                        eprintln!("  missing seg cw junction: {:?}", cw_cj);
                        eprintln!("  missing seg ccw junction: {:?}", ccw_cj);
                        let dist_ccw = (ccw_pos + n - cw_pos) % n;
                        eprintln!(
                            "  CCW distance from cw_pos to ccw_pos on trial: {} edges",
                            dist_ccw
                        );
                        eprintln!(
                            "  match length: {} edges (tile {} perimeter)",
                            pm.len, pm.tile_id
                        );
                        return;
                    }
                }
            }
        }
        panic!("expected to find a missing seg but didn't");
    }

    /// Spectre closure: O(patches × glues_per_patch) — large; gated.
    /// Reports progress every 1000 witness patches.
    #[test]
    #[ignore]
    fn spectre_seg_set_is_closed() {
        let bfs = SegmentTypeBFS::run(spectre_tileset(), 500_000, 500_000).unwrap();
        assert_seg_set_closed(&bfs, "spectre", 1_000);
    }
}
