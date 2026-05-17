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
use serde::{Deserialize, Serialize};

use crate::cyclotomic::{IsComplex, IsRingOrField, Units};
use crate::intgeom::matchtypes::MatchTypeIndex;
use crate::intgeom::patch::{CoarseJunction, EdgeInfo, GrowingPatch, PatchMatch};
use crate::intgeom::segbfs::SegmentTypeBFS;
use crate::intgeom::tileset::TileSet;

/// Canonical key for a `GrowingPatch` after `normalize()`.
type PatchKey = (Vec<i8>, Vec<EdgeInfo>, Vec<Vec<EdgeInfo>>, Vec<usize>);

/// Witness for one discovered sequence: the patch on which it was
/// first observed plus the match that produced it.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct SeqWitness {
    pub patch_id: usize,
    pub pm: PatchMatch,
}

/// The witness-independent metadata for one rewrite: which tile is
/// added, at what edge offset on that tile, how many edges it
/// absorbs, and where in the LHS its CW anchor sits.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, serde::Serialize, serde::Deserialize)]
pub struct MatchMeta {
    pub tile_id: usize,
    pub tile_offset: usize,
    pub match_len: usize,
    pub start_edge_in_lhs: usize,
}

/// One rewriting rule: an LHS seg sequence + MatchMeta + RHS, with
/// the canonical `(witness_patch_id, witness_pm)` that produced it.
/// Exactly one record per distinct `(LHS, MatchMeta)` pair.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct RuleRecord {
    pub lhs: Vec<usize>,
    pub meta: MatchMeta,
    pub rhs: [usize; 3],
    pub witness_patch_id: usize,
    pub witness_pm: PatchMatch,
}

/// All the information needed to reconstruct a `GrowingPatch` via
/// `from_parts`, in a serializable form (no generic `T`).
#[derive(Clone, Debug, Serialize, Deserialize)]
struct SerializedPatch {
    angles: Vec<i8>,
    edges: Vec<EdgeInfo>,
    inner_chains: Vec<Vec<EdgeInfo>>,
    patch_tile_ids: Vec<usize>,
    next_tile_id: usize,
}

impl<T: IsComplex + IsRingOrField + Units> From<&GrowingPatch<T>> for SerializedPatch {
    fn from(p: &GrowingPatch<T>) -> Self {
        SerializedPatch {
            angles: p.angles().to_vec(),
            edges: p.edges().to_vec(),
            inner_chains: p.inner_chains().to_vec(),
            patch_tile_ids: p.patch_tile_ids().to_vec(),
            next_tile_id: p.next_tile_id(),
        }
    }
}

/// On-disk format for `SegSeqBFS::save_to` / `load_from`.
#[derive(Serialize, Deserialize)]
struct Snapshot {
    patches: Vec<SerializedPatch>,
    seg_pairs: Vec<(CoarseJunction, CoarseJunction)>,
    sequences: Vec<Vec<usize>>,
    seq_witnesses: Vec<SeqWitness>,
    rules: Vec<RuleRecord>,
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

    /// Serialize the BFS to a binary file using bincode. Saves
    /// everything an analysis needs: rules, sequences, seq witnesses,
    /// patches (in `from_parts`-compatible form), and the seg DB
    /// (= seg_pairs only — patches of the inner SegmentTypeBFS are
    /// not saved). The tileset is the user's responsibility to
    /// provide on load (same tileset, same tile order).
    pub fn save_to(&self, path: impl AsRef<std::path::Path>) -> std::io::Result<()> {
        let snapshot = Snapshot {
            patches: self
                .patches
                .iter()
                .map(SerializedPatch::from)
                .collect(),
            seg_pairs: self.seg_bfs.seg_pairs().to_vec(),
            sequences: self.sequences.clone(),
            seq_witnesses: self.seq_witnesses.clone(),
            rules: self.rules.clone(),
        };
        let f = std::fs::File::create(path)?;
        let mut w = std::io::BufWriter::new(f);
        bincode::serialize_into(&mut w, &snapshot)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
        Ok(())
    }

    /// Reconstruct a `SegSeqBFS` from a saved snapshot. The user must
    /// provide the same tileset that was used to build the BFS.
    /// Internal lookup maps (`patch_lookup`, `seq_lookup`, `rule_lookup`,
    /// `seen_keys`) are rebuilt from the loaded data. The result has
    /// an empty queue and is not runnable; use it for analysis only.
    pub fn load_from(
        tileset: Arc<TileSet<T>>,
        path: impl AsRef<std::path::Path>,
    ) -> std::io::Result<Self> {
        let f = std::fs::File::open(path)?;
        let r = std::io::BufReader::new(f);
        let snapshot: Snapshot = bincode::deserialize_from(r)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;

        // Reconstruct the inner SegmentTypeBFS from the saved seg_pairs.
        let seg_bfs = SegmentTypeBFS::from_loaded(
            Arc::clone(&tileset),
            snapshot.seg_pairs,
        );

        // Reconstruct patches via from_parts using a shared match_index.
        let match_index = Arc::new(MatchTypeIndex::new(Arc::clone(&tileset)));
        let mut patches: Vec<GrowingPatch<T>> = Vec::with_capacity(snapshot.patches.len());
        for sp in snapshot.patches {
            let patch = GrowingPatch::from_parts(
                Arc::clone(&match_index),
                sp.angles,
                sp.edges,
                sp.inner_chains,
                sp.patch_tile_ids,
                sp.next_tile_id,
            )
            .ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "patch from_parts failed",
                )
            })?;
            patches.push(patch);
        }

        // Rebuild patch_lookup and seen_keys from the loaded patches.
        let mut patch_lookup: FxHashMap<PatchKey, usize> = FxHashMap::default();
        let mut seen_keys: FxHashSet<PatchKey> = FxHashSet::default();
        for (id, p) in patches.iter().enumerate() {
            let key = patch_key(p);
            patch_lookup.insert(key.clone(), id);
            seen_keys.insert(key);
        }

        // Rebuild seq_lookup.
        let mut seq_lookup: FxHashMap<Vec<usize>, usize> = FxHashMap::default();
        for (id, seq) in snapshot.sequences.iter().enumerate() {
            seq_lookup.insert(seq.clone(), id);
        }

        // Rebuild rule_lookup.
        let mut rule_lookup: FxHashMap<(Vec<usize>, MatchMeta), usize> =
            FxHashMap::default();
        for (id, r) in snapshot.rules.iter().enumerate() {
            rule_lookup.insert((r.lhs.clone(), r.meta), id);
        }

        Ok(Self {
            tileset,
            seg_bfs,
            patches,
            patch_lookup,
            seen_keys,
            sequences: snapshot.sequences,
            seq_lookup,
            seq_witnesses: snapshot.seq_witnesses,
            rules: snapshot.rules,
            rule_lookup,
            queue: VecDeque::new(),
            seq_cap: 0,
            patch_cap: 0,
        })
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
        // (LHS, MatchMeta, RHS, pm) entries to commit if this patch
        // turns out to be a canonical witness.
        let mut new_rules_local: Vec<(Vec<usize>, MatchMeta, [usize; 3], PatchMatch)> =
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

            // Compute the RHS by walking the un-normalized trial
            // between the preserved outer junctions. Must happen
            // BEFORE normalize because the position arithmetic
            // depends on the add_tile layout.
            let rhs = match compute_rhs(
                &trial, &juncs, n, &pm, start_seg, end_seg, &self.seg_bfs,
            ) {
                Some(r) => r,
                None => continue,
            };

            trial.normalize();
            self.enqueue_if_unseen(trial)?;
            new_rules_local.push((seq.clone(), meta, rhs, pm.clone()));
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
        for (lhs, meta, rhs, pm) in new_rules_local {
            let rule_id = self.rules.len();
            self.rule_lookup.insert((lhs.clone(), meta), rule_id);
            self.rules.push(RuleRecord {
                lhs,
                meta,
                rhs,
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

/// Compute the RHS of a rewrite by walking the un-normalized trial
/// from the CW preserved-outer junction CCW to the CCW one. Returns
/// the 3 seg ids `[pad_cw', central, pad_ccw']` or `None` if the
/// walk fails / a seg isn't in the seg DB.
pub(crate) fn compute_rhs<T: IsComplex + IsRingOrField + Units>(
    trial: &GrowingPatch<T>,
    src_juncs: &[usize],
    src_n: usize,
    pm: &PatchMatch,
    start_seg: usize,
    end_seg: usize,
    seg_bfs: &SegmentTypeBFS<T>,
) -> Option<[usize; 3]> {
    let k = src_juncs.len();
    let j_outer_cw_src = src_juncs[start_seg];
    let j_outer_ccw_src = src_juncs[(end_seg + 1) % k];
    let aug_n = trial.boundary_len();
    let gap_start = src_n - pm.len;

    // After add_tile, the survivor strip on the trial occupies
    // positions [0, gap_start] with CCW anchor at 0 and CW anchor at
    // gap_start. Source survivor vertex `v` lives on the trial at
    // `(v - (start_a + match_len) + n) % n`.
    let ccw_anchor_src = (pm.start_a + pm.len) % src_n;
    let j_outer_cw_trial = (j_outer_cw_src + src_n - ccw_anchor_src) % src_n;
    let j_outer_ccw_trial = (j_outer_ccw_src + src_n - ccw_anchor_src) % src_n;
    debug_assert!(j_outer_cw_trial <= gap_start);
    debug_assert!(j_outer_ccw_trial <= gap_start);

    let mut walk: Vec<CoarseJunction> = Vec::with_capacity(4);
    let mut pos = j_outer_cw_trial;
    for _ in 0..=aug_n {
        if let Some(cj) = trial.coarse_junction_at(pos) {
            walk.push(cj);
            if pos == j_outer_ccw_trial && walk.len() >= 2 {
                break;
            }
        }
        pos = (pos + 1) % aug_n;
        if walk.len() > aug_n {
            return None;
        }
    }
    if walk.len() != 4 {
        panic!(
            "expected 4 junctions in RHS walk, got {}: pm={:?}, src_n={}, aug_n={}, gap_start={}, src_juncs={:?}, start_seg={}, end_seg={}, walk={:?}",
            walk.len(),
            pm, src_n, aug_n, gap_start,
            src_juncs, start_seg, end_seg,
            walk,
        );
    }
    Some([
        seg_bfs.lookup_seg(&walk[0], &walk[1])?,
        seg_bfs.lookup_seg(&walk[1], &walk[2])?,
        seg_bfs.lookup_seg(&walk[2], &walk[3])?,
    ])
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
