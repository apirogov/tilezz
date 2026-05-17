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

/// The witness-independent metadata for one rewrite: which tile is
/// added, at what edge offset on that tile, how many edges it
/// absorbs, and where in the LHS its CW anchor sits.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub struct MatchMeta {
    pub tile_id: usize,
    pub tile_offset: usize,
    pub match_len: usize,
    pub start_edge_in_lhs: usize,
}

/// One rewriting rule: an LHS seg sequence + MatchMeta + the
/// canonical `(witness_patch_id, witness_pm)` that produced it.
/// Exactly one record per distinct `(LHS, MatchMeta)` pair.
#[derive(Clone, Debug)]
pub struct RuleRecord {
    pub lhs: Vec<usize>,
    pub meta: MatchMeta,
    pub witness_patch_id: usize,
    pub witness_pm: PatchMatch,
}

#[derive(Debug)]
pub enum SeqBfsError {
    SeqCapExceeded { cap: usize },
    PatchCapExceeded { cap: usize },
}

pub struct SegSeqBFS<T: IsComplex> {
    tileset: Arc<TileSet<T>>,
    seg_bfs: SegmentTypeBFS<T>,

    // Witness patches: a patch is interned here iff it is the
    // canonical witness for at least one new rule. `patch_id` indexes
    // both `patches` and `patch_lookup`.
    patches: Vec<GrowingPatch<T>>,
    patch_lookup: FxHashMap<PatchKey, usize>,

    // Dedup for the BFS: every key we've ever enqueued (= prevents
    // re-processing the same patch). Includes witnesses' keys.
    seen_keys: FxHashSet<PatchKey>,

    // Discovered sequences. Each LHS appears exactly once. `seq_witnesses`
    // points to the first patch+match that produced it.
    sequences: Vec<Vec<usize>>,
    seq_lookup: FxHashMap<Vec<usize>, usize>,
    seq_witnesses: Vec<SeqWitness>,

    // Discovered rules. One record per distinct `(LHS, MatchMeta)`
    // pair, with its canonical witness.
    rules: Vec<RuleRecord>,
    rule_lookup: FxHashMap<(Vec<usize>, MatchMeta), usize>,

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
            rules: Vec::new(),
            rule_lookup: FxHashMap::default(),
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

    /// Distinct `(LHS, MatchMeta)` rules observed. Each is backed by
    /// exactly one canonical witness pair `(witness_patch_id, witness_pm)`
    /// stored in `rules[rule_id]`.
    pub fn num_rules(&self) -> usize {
        self.rules.len()
    }

    /// All discovered rules, each with its canonical witness.
    pub fn rules(&self) -> &[RuleRecord] {
        &self.rules
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

        // Enumerate matches. A trial is enqueued the first time its
        // (LHS, MatchMeta) is observed globally; that same first
        // occurrence pins the rule's canonical witness. Subsequent
        // matches producing the same (LHS, MatchMeta) — whether on
        // this witness or any future one — are silently discarded.
        // (LHS, MatchMeta, pm) entries to commit if this patch turns
        // out to be a canonical witness.
        let mut new_rules_local: Vec<(Vec<usize>, MatchMeta, PatchMatch)> =
            Vec::new();
        let mut new_seqs_local: Vec<(Vec<usize>, PatchMatch)> = Vec::new();
        let mut new_seqs_seen: FxHashSet<Vec<usize>> = FxHashSet::default();
        let mut new_rules_seen: FxHashSet<(Vec<usize>, MatchMeta)> =
            FxHashSet::default();
        for pm in patch.get_all_matches() {
            if pm.len == 0 {
                continue;
            }
            let (start_seg, end_seg) = lhs_seg_range(&juncs, n, &pm);

            // Build the LHS sequence: [start_seg, internal..., end_seg].
            let mut seq = Vec::with_capacity(k);
            let mut idx = start_seg;
            loop {
                seq.push(seg_ids[idx]);
                if idx == end_seg {
                    break;
                }
                idx = (idx + 1) % k;
                debug_assert!(seq.len() <= k + 1, "infinite loop in seq build");
            }

            // Compute MatchMeta from the same LHS bracketing.
            let j_outer_cw_src = juncs[start_seg];
            let meta = MatchMeta {
                tile_id: pm.tile_id,
                tile_offset: pm.start_b,
                match_len: pm.len,
                start_edge_in_lhs: (pm.start_a + n - j_outer_cw_src) % n,
            };

            // Rule-based dedup: only do anything if the rule is new.
            let rule_key = (seq.clone(), meta);
            if !new_rules_seen.insert(rule_key.clone()) {
                continue;
            }
            if self.rule_lookup.contains_key(&rule_key) {
                continue;
            }

            // Validate the match — fail means this `(LHS, MatchMeta)`
            // is not realizable on this witness, but might still arise
            // on a later one. Don't insert into rule_lookup yet.
            let mut trial = patch.clone();
            if !trial.add_tile(&pm) {
                continue;
            }
            trial.normalize();
            self.enqueue_if_unseen(trial)?;
            new_rules_local.push((seq.clone(), meta, pm.clone()));
            if new_seqs_seen.insert(seq.clone())
                && !self.seq_lookup.contains_key(&seq)
            {
                new_seqs_local.push((seq, pm.clone()));
            }
        }

        if new_rules_local.is_empty() {
            // No new rules — this patch contributes nothing canonical
            // and is not retained.
            return Ok(());
        }

        // Register this patch as a canonical witness.
        let patch_id = self.register_witness(patch);

        // Allocate new seq slots in discovery order.
        for (seq, primary_pm) in new_seqs_local {
            let seq_id = self.sequences.len();
            self.seq_lookup.insert(seq.clone(), seq_id);
            self.sequences.push(seq);
            self.seq_witnesses.push(SeqWitness {
                patch_id,
                pm: primary_pm,
            });
            if self.sequences.len() > self.seq_cap {
                return Err(SeqBfsError::SeqCapExceeded { cap: self.seq_cap });
            }
        }

        // Register new rules with their canonical witness.
        for (lhs, meta, pm) in new_rules_local {
            let rule_id = self.rules.len();
            self.rule_lookup.insert((lhs.clone(), meta), rule_id);
            self.rules.push(RuleRecord {
                lhs,
                meta,
                witness_patch_id: patch_id,
                witness_pm: pm,
            });
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

/// Compute the LHS segment range `(start_seg, end_seg)` for a match.
///
/// The match touches edges `[start_a, start_a + len)`. Its anchors are
/// the vertices `cw_anchor = start_a` and `ccw_anchor = start_a + len`.
/// The LHS is the CCW arc of segs from `start_seg` to `end_seg`.
///
/// If an anchor is strictly inside a segment, that segment is partially
/// absorbed by the match (= the anchor becomes a new junction). If an
/// anchor lands **on** an existing source junction, the *neighbor*
/// seg on that side is pulled into the LHS instead — because the
/// junction's CJ is modified by the glue and is no longer a "preserved
/// outer junction". The extension keeps the anchor strictly interior
/// to the LHS sweep so the trial walk always finds 4 junctions on the
/// rewrite side (= 3 RHS segs).
pub(crate) fn lhs_seg_range(
    juncs: &[usize],
    n: usize,
    pm: &PatchMatch,
) -> (usize, usize) {
    let k = juncs.len();
    let cw_anchor = pm.start_a;
    let ccw_anchor = (pm.start_a + pm.len) % n;
    let end_edge = (pm.start_a + pm.len - 1) % n;
    let mut start_seg = seg_containing_edge(juncs, pm.start_a, n);
    let mut end_seg = seg_containing_edge(juncs, end_edge, n);
    if juncs[start_seg] == cw_anchor {
        start_seg = (start_seg + k - 1) % k;
    }
    if juncs[(end_seg + 1) % k] == ccw_anchor {
        end_seg = (end_seg + 1) % k;
    }
    // PADDING INVARIANT: the outer junctions of the LHS sweep must be
    // untouched by the match (= not coincide with either anchor). If
    // either outer lands on an anchor, the "preserved outer" promise
    // is broken and the rewrite is meaningless. This can only happen
    // if BOTH anchors land on existing junctions in a patch with too
    // few junctions to extend past them on both sides (e.g. k=2 with
    // a match exactly between the two junctions — geometrically
    // impossible for spectre, per construction).
    let cw_outer = juncs[start_seg];
    let ccw_outer = juncs[(end_seg + 1) % k];
    assert!(
        cw_outer != cw_anchor
            && cw_outer != ccw_anchor
            && ccw_outer != cw_anchor
            && ccw_outer != ccw_anchor,
        "padding invariant violated: outer junctions coincide with a \
         match anchor (juncs={:?}, n={}, pm={:?}, start_seg={}, end_seg={}, \
         cw_outer={}, ccw_outer={}, cw_anchor={}, ccw_anchor={})",
        juncs, n, pm, start_seg, end_seg, cw_outer, ccw_outer, cw_anchor, ccw_anchor,
    );
    (start_seg, end_seg)
}

/// Given the cyclic-CCW-sorted junction vertex positions `juncs` on
/// a boundary of length `n`, find the index `i` of the segment that
/// owns the edge at position `e`. The segment with index `i` covers
/// edges `[juncs[i], juncs[(i+1) % k])` (mod n).
pub(crate) fn seg_containing_edge(juncs: &[usize], e: usize, n: usize) -> usize {
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

    fn spectre_tileset() -> Arc<TileSet<ZZ12>> {
        Arc::new(TileSet::new(vec![
            Rat::try_from(&tiles::spectre::<ZZ12>()).unwrap(),
        ]))
    }

    fn mixed_tileset() -> Arc<TileSet<ZZ12>> {
        Arc::new(TileSet::new(vec![
            Rat::try_from(&tiles::square::<ZZ12>()).unwrap(),
            Rat::try_from(&tiles::hexagon::<ZZ12>()).unwrap(),
        ]))
    }

    #[test]
    #[ignore]
    fn mixed_sequences_enumerated() {
        let bfs = SegSeqBFS::run(mixed_tileset(), 20_000_000, 20_000_000)
            .expect("mixed SegSeqBFS terminates within caps");
        eprintln!(
            "mixed (square+hex): sequences={}, witness_patches={}, seen={} (total enqueued)",
            bfs.num_sequences(),
            bfs.num_witness_patches(),
            bfs.num_patches_seen()
        );
        assert!(bfs.num_sequences() > 0);
    }

    #[test]
    #[ignore]
    fn spectre_sequences_enumerated() {
        let bfs = SegSeqBFS::run(spectre_tileset(), 10_000_000, 10_000_000)
            .expect("spectre SegSeqBFS terminates within caps");
        eprintln!(
            "spectre: sequences={}, witness_patches={}, seen={} (total enqueued)",
            bfs.num_sequences(),
            bfs.num_witness_patches(),
            bfs.num_patches_seen()
        );
        assert!(bfs.num_sequences() > 0);
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
