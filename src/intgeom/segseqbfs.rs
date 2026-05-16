//! Segment-sequence BFS: enumerate all distinct "patch-tile matching
//! contexts" as segment-id sequences.
//!
//! Built on top of a pre-computed [`SegmentTypeBFS`]: each match on
//! any reachable patch is characterized by the smallest seg-sequence
//! on the patch boundary that **fully contains** the match's
//! absorbed run, padded by one segment on each side
//! (= `[pad_cw, seg_with_match_start, ..., seg_with_match_end,
//! pad_ccw]`). The seg ids come from the seg DB; thanks to seg
//! closure, the lookup never fails.
//!
//! The BFS explores all reachable patches (seeds = 2-tile glues,
//! then trials of every legal match). Each patch is deduped
//! canonically. When a patch contributes any new sequence, it is
//! retained as the witness for those sequences; otherwise it stays
//! in the patch pool only to dedup future trials (no new info).
//! Per sequence, we record `(patch_id, applied_match)` so the
//! geometric context is recoverable.
//!
//! Termination: queue empties when every reachable patch has been
//! explored and no new sequences arise. Bounded for finite-seg
//! tilings.

use std::collections::VecDeque;
use std::sync::Arc;

use rustc_hash::FxHashMap;

use crate::cyclotomic::{IsComplex, IsRingOrField, Units};
use crate::intgeom::patch::{EdgeInfo, GrowingPatch, PatchMatch};
use crate::intgeom::segbfs::SegmentTypeBFS;
use crate::intgeom::tileset::TileSet;

/// Canonical key for a `GrowingPatch` after `normalize()`.
type PatchKey = (Vec<i8>, Vec<EdgeInfo>, Vec<Vec<EdgeInfo>>, Vec<usize>);

/// Witness for one discovered sequence: the patch on which it was
/// first observed plus the match that produced it.
#[derive(Clone, Debug)]
pub struct SeqWitness {
    pub patch_id: usize,
    pub pm: PatchMatch,
}

#[derive(Debug)]
pub enum SeqBfsError {
    SeqCapExceeded { cap: usize },
    PatchCapExceeded { cap: usize },
}

pub struct SegSeqBFS<T: IsComplex> {
    tileset: Arc<TileSet<T>>,
    seg_bfs: SegmentTypeBFS<T>,

    // Deduped patches. patch_lookup ensures each unique patch is
    // added at most once. `is_witness[i]` flags whether patch i
    // contributed any new sequences.
    patches: Vec<GrowingPatch<T>>,
    patch_lookup: FxHashMap<PatchKey, usize>,
    is_witness: Vec<bool>,

    // Discovered sequences with witnesses.
    sequences: Vec<Vec<usize>>,
    seq_lookup: FxHashMap<Vec<usize>, usize>,
    seq_witnesses: Vec<SeqWitness>,

    // BFS frontier: patch ids waiting to be explored.
    queue: VecDeque<usize>,

    seq_cap: usize,
    patch_cap: usize,
}

impl<T: IsComplex + IsRingOrField + Units> SegSeqBFS<T> {
    /// Run the segment-sequence BFS to completion.
    ///
    /// First builds a `SegmentTypeBFS` for the same tileset (= the
    /// seg DB), then enumerates sequences as described in the module
    /// doc.
    pub fn run(
        tileset: Arc<TileSet<T>>,
        seq_cap: usize,
        patch_cap: usize,
    ) -> Result<Self, SeqBfsError> {
        let seg_bfs = SegmentTypeBFS::run(Arc::clone(&tileset), seq_cap.max(500_000), patch_cap.max(500_000))
            .map_err(|_| SeqBfsError::SeqCapExceeded { cap: seq_cap })?;
        let mut bfs = Self {
            tileset,
            seg_bfs,
            patches: Vec::new(),
            patch_lookup: FxHashMap::default(),
            is_witness: Vec::new(),
            sequences: Vec::new(),
            seq_lookup: FxHashMap::default(),
            seq_witnesses: Vec::new(),
            queue: VecDeque::new(),
            seq_cap,
            patch_cap,
        };
        bfs.seed()?;
        bfs.run_bfs()?;
        Ok(bfs)
    }

    pub fn num_sequences(&self) -> usize {
        self.sequences.len()
    }

    pub fn num_patches(&self) -> usize {
        self.patches.len()
    }

    pub fn num_witness_patches(&self) -> usize {
        self.is_witness.iter().filter(|&&w| w).count()
    }

    pub fn sequences(&self) -> &[Vec<usize>] {
        &self.sequences
    }

    pub fn seq_witness(&self, seq_id: usize) -> Option<&SeqWitness> {
        self.seq_witnesses.get(seq_id)
    }

    pub fn patches(&self) -> &[GrowingPatch<T>] {
        &self.patches
    }

    pub fn seg_bfs(&self) -> &SegmentTypeBFS<T> {
        &self.seg_bfs
    }

    fn seed(&mut self) -> Result<(), SeqBfsError> {
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
                self.intern_patch(patch)?;
            }
        }
        Ok(())
    }

    /// Intern (dedup) a patch. If new, push onto `queue`. Returns the
    /// patch id (whether new or pre-existing).
    fn intern_patch(&mut self, patch: GrowingPatch<T>) -> Result<usize, SeqBfsError> {
        let key = patch_key(&patch);
        if let Some(&id) = self.patch_lookup.get(&key) {
            return Ok(id);
        }
        let id = self.patches.len();
        self.patch_lookup.insert(key, id);
        self.patches.push(patch);
        self.is_witness.push(false);
        self.queue.push_back(id);
        if self.patches.len() > self.patch_cap {
            return Err(SeqBfsError::PatchCapExceeded {
                cap: self.patch_cap,
            });
        }
        Ok(id)
    }

    fn run_bfs(&mut self) -> Result<(), SeqBfsError> {
        let start = std::time::Instant::now();
        let mut explored = 0usize;
        while let Some(patch_id) = self.queue.pop_front() {
            self.explore_patch(patch_id)?;
            explored += 1;
            if explored % 1000 == 0 {
                eprintln!(
                    "segseqbfs: explored {} patches, total interned={}, witnesses={}, sequences={}, queue={}, elapsed={:.1}s",
                    explored,
                    self.patches.len(),
                    self.num_witness_patches(),
                    self.sequences.len(),
                    self.queue.len(),
                    start.elapsed().as_secs_f64(),
                );
            }
        }
        Ok(())
    }

    fn explore_patch(&mut self, patch_id: usize) -> Result<(), SeqBfsError> {
        let patch = self.patches[patch_id].clone();
        let n = patch.boundary_len();
        if n == 0 {
            return Ok(());
        }
        let juncs: Vec<usize> = (0..n).filter(|&i| patch.is_junction(i)).collect();
        let k = juncs.len();
        if k < 2 {
            return Ok(());
        }
        // Resolve each segment between consecutive junctions to a seg
        // id via the seg DB. (Closure guarantees the lookup never
        // fails for patches reachable from seeds.)
        let mut seg_ids: Vec<usize> = Vec::with_capacity(k);
        for j in 0..k {
            let cw_cj = patch.coarse_junction_at(juncs[j]).expect("junction");
            let ccw_cj = patch
                .coarse_junction_at(juncs[(j + 1) % k])
                .expect("junction");
            let id = self
                .seg_bfs
                .lookup_seg(&cw_cj, &ccw_cj)
                .expect("seg must be in DB (closure verified)");
            seg_ids.push(id);
        }

        let mut any_new_seq = false;
        for pm in patch.get_all_matches() {
            if pm.len == 0 {
                continue;
            }
            let start_edge = pm.start_a;
            let end_edge = (pm.start_a + pm.len - 1) % n;
            let start_seg = seg_containing_edge(&juncs, start_edge, n);
            let end_seg = seg_containing_edge(&juncs, end_edge, n);

            // Build padded match sequence:
            //   [seg_ids[cw_pad], seg_ids[start_seg], ..., seg_ids[end_seg], seg_ids[ccw_pad]]
            let cw_pad = (start_seg + k - 1) % k;
            let ccw_pad = (end_seg + 1) % k;
            let mut seq = Vec::with_capacity(k + 2);
            seq.push(seg_ids[cw_pad]);
            let mut idx = start_seg;
            loop {
                seq.push(seg_ids[idx]);
                if idx == end_seg {
                    break;
                }
                idx = (idx + 1) % k;
                debug_assert!(seq.len() <= k + 2, "infinite loop in seq build");
            }
            seq.push(seg_ids[ccw_pad]);

            if self.seq_lookup.contains_key(&seq) {
                // Sequence already seen — don't enqueue trial; this
                // glue isn't going to teach us anything new from
                // this context.
                continue;
            }

            // New sequence: register, record witness, and apply the
            // match so the resulting trial gets enqueued (only
            // productive matches drive BFS exploration further).
            let seq_id = self.sequences.len();
            self.seq_lookup.insert(seq.clone(), seq_id);
            self.sequences.push(seq);
            self.seq_witnesses.push(SeqWitness {
                patch_id,
                pm: pm.clone(),
            });
            any_new_seq = true;
            if self.sequences.len() > self.seq_cap {
                return Err(SeqBfsError::SeqCapExceeded { cap: self.seq_cap });
            }

            let mut trial = patch.clone();
            if trial.add_tile(&pm) {
                trial.normalize();
                self.intern_patch(trial)?;
            }
        }

        if any_new_seq {
            self.is_witness[patch_id] = true;
        }
        Ok(())
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

/// Given the cyclic-CCW-sorted junction vertex positions `juncs` on
/// a boundary of length `n`, find the index `i` of the segment that
/// owns the edge at position `e`. The segment with index `i` covers
/// edges `[juncs[i], juncs[(i+1) % k])` (mod n).
fn seg_containing_edge(juncs: &[usize], e: usize, n: usize) -> usize {
    let k = juncs.len();
    for i in 0..k {
        let start = juncs[i];
        let end = juncs[(i + 1) % k];
        if start <= end {
            if e >= start && e < end {
                return i;
            }
        } else {
            // Range wraps around.
            if e >= start || e < end {
                return i;
            }
        }
    }
    panic!(
        "edge {} not in any seg (n={}, juncs={:?})",
        e, n, juncs
    );
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

    #[test]
    fn hex_sequences_enumerated() {
        let bfs = SegSeqBFS::run(hex_tileset(), 1_000_000, 1_000_000)
            .expect("hex SegSeqBFS terminates within caps");
        eprintln!(
            "hex: sequences={}, witness_patches={}/{} (total interned)",
            bfs.num_sequences(),
            bfs.num_witness_patches(),
            bfs.num_patches()
        );
        assert!(bfs.num_sequences() > 0);
    }
}
