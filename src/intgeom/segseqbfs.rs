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

use rustc_hash::{FxHashMap, FxHashSet};

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

    // Witness patches only. patch_lookup maps canonical keys to
    // witness slot ids. Patches that don't contribute new sequences
    // are never interned; they're explored from the queue and
    // dropped on the spot.
    patches: Vec<GrowingPatch<T>>,
    patch_lookup: FxHashMap<PatchKey, usize>,

    // Dedup for the BFS: every key we've ever enqueued (= prevents
    // re-processing the same patch). Includes witnesses' keys.
    seen_keys: FxHashSet<PatchKey>,

    // Discovered sequences with witnesses (= patch_id + applied match).
    sequences: Vec<Vec<usize>>,
    seq_lookup: FxHashMap<Vec<usize>, usize>,
    seq_witnesses: Vec<SeqWitness>,

    // BFS frontier: patches awaiting exploration. Owned by value,
    // since we may discard them entirely on exploration.
    queue: VecDeque<GrowingPatch<T>>,

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
            seen_keys: FxHashSet::default(),
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

    /// Number of distinct patches ever enqueued (= seen_keys size).
    pub fn num_patches_seen(&self) -> usize {
        self.seen_keys.len()
    }

    /// Number of patches retained as witnesses.
    pub fn num_witness_patches(&self) -> usize {
        self.patches.len()
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

    pub fn patch(&self, patch_id: usize) -> Option<&GrowingPatch<T>> {
        self.patches.get(patch_id)
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
                self.enqueue_if_unseen(patch)?;
            }
        }
        Ok(())
    }

    /// Compute the patch's key; if not seen before, mark seen and
    /// push onto the BFS queue. Otherwise drop the patch.
    fn enqueue_if_unseen(&mut self, patch: GrowingPatch<T>) -> Result<(), SeqBfsError> {
        let key = patch_key(&patch);
        if self.seen_keys.contains(&key) {
            return Ok(());
        }
        self.seen_keys.insert(key);
        self.queue.push_back(patch);
        if self.seen_keys.len() > self.patch_cap {
            return Err(SeqBfsError::PatchCapExceeded {
                cap: self.patch_cap,
            });
        }
        Ok(())
    }

    /// Register a patch as a witness (only called after we know it
    /// contributes new sequences). Returns the new slot id.
    fn register_witness(&mut self, patch: GrowingPatch<T>) -> usize {
        let key = patch_key(&patch);
        let id = self.patches.len();
        self.patch_lookup.insert(key, id);
        self.patches.push(patch);
        id
    }

    fn run_bfs(&mut self) -> Result<(), SeqBfsError> {
        let start = std::time::Instant::now();
        let mut explored = 0usize;
        while let Some(patch) = self.queue.pop_front() {
            self.explore_patch(patch)?;
            explored += 1;
            if explored % 1000 == 0 {
                eprintln!(
                    "segseqbfs: explored {} patches, seen={}, witnesses={}, sequences={}, queue={}, elapsed={:.1}s",
                    explored,
                    self.seen_keys.len(),
                    self.patches.len(),
                    self.sequences.len(),
                    self.queue.len(),
                    start.elapsed().as_secs_f64(),
                );
            }
        }
        Ok(())
    }

    fn explore_patch(&mut self, patch: GrowingPatch<T>) -> Result<(), SeqBfsError> {
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

        // Enumerate matches. Each match producing a (globally-)new
        // sequence: enqueue its trial eagerly and note (seq, pm) for
        // later registration. Duplicates within this source — same
        // seq from different matches — still enqueue distinct
        // trials (they may be different patches under normalize),
        // but the seq is only noted once.
        let mut local_new_seqs: Vec<(Vec<usize>, PatchMatch)> = Vec::new();
        let mut local_seen: FxHashSet<Vec<usize>> = FxHashSet::default();
        for pm in patch.get_all_matches() {
            if pm.len == 0 {
                continue;
            }
            let start_edge = pm.start_a;
            let end_edge = (pm.start_a + pm.len - 1) % n;
            let start_seg = seg_containing_edge(&juncs, start_edge, n);
            let end_seg = seg_containing_edge(&juncs, end_edge, n);

            // Build padded match sequence.
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
                continue;
            }
            if !local_seen.insert(seq.clone()) {
                // Duplicate seq within this source — skip both the
                // registration and the trial enqueue (= match the
                // original BFS behavior where seq_lookup was updated
                // inline).
                continue;
            }

            // Trial.
            let mut trial = patch.clone();
            if !trial.add_tile(&pm) {
                continue;
            }
            trial.normalize();
            self.enqueue_if_unseen(trial)?;
            local_new_seqs.push((seq, pm.clone()));
        }

        if local_new_seqs.is_empty() {
            return Ok(());
        }

        // Source contributed — register it as a witness and claim
        // the noted sequences.
        let patch_id = self.register_witness(patch);
        for (seq, pm) in local_new_seqs {
            let seq_id = self.sequences.len();
            self.seq_lookup.insert(seq.clone(), seq_id);
            self.sequences.push(seq);
            self.seq_witnesses.push(SeqWitness { patch_id, pm });
            if self.sequences.len() > self.seq_cap {
                return Err(SeqBfsError::SeqCapExceeded { cap: self.seq_cap });
            }
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
            "hex: sequences={}, witness_patches={}, seen={} (total enqueued)",
            bfs.num_sequences(),
            bfs.num_witness_patches(),
            bfs.num_patches_seen()
        );
        assert!(bfs.num_sequences() > 0);
    }
}
