//! Segment-sequence rewrite rules extracted from a finished
//! [`SegSeqBFS`].
//!
//! For each (witness patch, match) pair in the BFS output, we
//! extract a rewrite rule:
//!   - **LHS** = the padded match sequence (= same as the BFS
//!     records).
//!   - **MatchMeta** = the tile + tile_offset + match_len +
//!     `start_edge_in_lhs` describing the rewrite.
//!   - **RHS** = exactly 3 segment ids `[pad_cw', central, pad_ccw']`
//!     describing the post-glue boundary between the unchanged
//!     outer junctions.
//!
//! Witness-independence: for any (LHS, MatchMeta), the RHS is
//! uniquely determined. Different witnesses producing the same
//! (LHS, MatchMeta) yield the same RHS, so dedup is safe and the
//! table holds one entry per distinct rule.

use std::sync::Arc;

use rustc_hash::FxHashMap;

use crate::cyclotomic::{IsComplex, IsRingOrField, Units};
use crate::intgeom::patch::{CoarseJunction, GrowingPatch, PatchMatch};
use crate::intgeom::segseqbfs::{lhs_seg_range, SegSeqBFS};

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

#[derive(Debug)]
pub struct SegRewriteTable {
    /// LHS sequence → all distinct outgoing transitions. Each RHS is
    /// exactly 3 seg ids `[pad_cw', central, pad_ccw']` between the
    /// preserved outer junctions on the trial.
    pub by_lhs: FxHashMap<Vec<usize>, FxHashMap<MatchMeta, [usize; 3]>>,
}

impl SegRewriteTable {
    /// Build the table eagerly from a finished [`SegSeqBFS`].
    /// Returns the table along with the number of `(LHS, MatchMeta)`
    /// collisions observed during the build (= rule re-extracted from
    /// multiple witnesses with the same key). Collisions are
    /// asserted to produce the same RHS (= witness-independence).
    pub fn build_from<T: IsComplex + IsRingOrField + Units>(
        seq_bfs: &SegSeqBFS<T>,
    ) -> Self {
        let (table, _collisions) = Self::build_with_stats(seq_bfs);
        table
    }

    /// Same as `build_from` but also returns the collision count
    /// (= number of times a `(LHS, MatchMeta)` was re-extracted from
    /// a different witness).
    pub fn build_with_stats<T: IsComplex + IsRingOrField + Units>(
        seq_bfs: &SegSeqBFS<T>,
    ) -> (Self, usize) {
        Self::build_filtered(seq_bfs, |_| true)
    }

    /// Build the table from only the witnesses that pass `keep`.
    /// Used for diagnostics (= compare "contributing-only" vs "all").
    pub fn build_filtered<T: IsComplex + IsRingOrField + Units>(
        seq_bfs: &SegSeqBFS<T>,
        keep: impl Fn(usize) -> bool,
    ) -> (Self, usize) {
        let mut by_lhs: FxHashMap<Vec<usize>, FxHashMap<MatchMeta, [usize; 3]>> =
            FxHashMap::default();
        let mut collisions = 0usize;
        let seg_bfs = seq_bfs.seg_bfs();
        for (patch_id, patch) in seq_bfs.patches().iter().enumerate() {
            if !keep(patch_id) {
                continue;
            }
            let n = patch.boundary_len();
            if n == 0 {
                continue;
            }
            let juncs: Vec<usize> =
                (0..n).filter(|&i| patch.is_junction(i)).collect();
            let k = juncs.len();
            if k < 2 {
                continue;
            }
            // Cache seg ids for the patch's k segments (one per
            // consecutive junction pair, CCW).
            let mut seg_ids: Vec<usize> = Vec::with_capacity(k);
            for j in 0..k {
                let cw_cj = patch
                    .coarse_junction_at(juncs[j])
                    .expect("known junction");
                let ccw_cj = patch
                    .coarse_junction_at(juncs[(j + 1) % k])
                    .expect("known junction");
                let id = seg_bfs
                    .lookup_seg(&cw_cj, &ccw_cj)
                    .expect("seg must be in DB");
                seg_ids.push(id);
            }
            for pm in patch.get_all_matches() {
                if pm.len == 0 {
                    continue;
                }
                let (lhs, meta, rhs) = match extract_rule(
                    patch, n, &juncs, &seg_ids, &pm, seg_bfs,
                ) {
                    Some(triple) => triple,
                    None => continue,
                };
                let entries = by_lhs.entry(lhs).or_default();
                if let Some(prev) = entries.insert(meta, rhs) {
                    // Witness-independence: the same (LHS, MatchMeta)
                    // must produce the same RHS regardless of which
                    // witness we extracted from.
                    assert_eq!(
                        prev, rhs,
                        "witness-independence violated for meta {:?}",
                        meta
                    );
                    collisions += 1;
                }
            }
        }
        (SegRewriteTable { by_lhs }, collisions)
    }

    /// Distinct LHS sequences.
    pub fn num_lhs(&self) -> usize {
        self.by_lhs.len()
    }

    /// Total distinct (LHS, MatchMeta, RHS) rules.
    pub fn num_rules(&self) -> usize {
        self.by_lhs.values().map(|m| m.len()).sum()
    }

    /// All rewrites firing on a given LHS, if any.
    pub fn rules_for_lhs(&self, lhs: &[usize]) -> Option<&FxHashMap<MatchMeta, [usize; 3]>> {
        self.by_lhs.get(lhs)
    }
}

/// Try to extract a rule from `(patch, pm)`. Returns `None` if the
/// match is degenerate or the trial fails the layout assumption.
fn extract_rule<T: IsComplex + IsRingOrField + Units>(
    patch: &GrowingPatch<T>,
    n: usize,
    juncs: &[usize],
    seg_ids: &[usize],
    pm: &PatchMatch,
    seg_bfs: &crate::intgeom::segbfs::SegmentTypeBFS<T>,
) -> Option<(Vec<usize>, MatchMeta, [usize; 3])> {
    let k = juncs.len();
    let (start_seg, end_seg) = lhs_seg_range(juncs, n, pm);

    // Build the LHS sequence: [start_seg, internal..., end_seg].
    // Same definition as SegSeqBFS (see `lhs_seg_range`).
    let mut lhs = Vec::with_capacity(k);
    let mut idx = start_seg;
    loop {
        lhs.push(seg_ids[idx]);
        if idx == end_seg {
            break;
        }
        idx = (idx + 1) % k;
        debug_assert!(lhs.len() <= k + 1, "infinite loop building LHS");
    }

    // Preserved outer junctions on the source:
    //   CW outer = start_seg's CW endpoint (= junction CW of start_seg).
    //   CCW outer = end_seg's CCW endpoint (= junction CCW of end_seg).
    // After `lhs_seg_range` extension, both are strictly outside the
    // match's anchor closure.
    let j_outer_cw_src = juncs[start_seg];
    let j_outer_ccw_src = juncs[(end_seg + 1) % k];

    // Match metadata.
    let start_edge_in_lhs = (pm.start_a + n - j_outer_cw_src) % n;
    let meta = MatchMeta {
        tile_id: pm.tile_id,
        tile_offset: pm.start_b,
        match_len: pm.len,
        start_edge_in_lhs,
    };

    // Apply the match.
    let mut trial = patch.clone();
    if !trial.add_tile(pm) {
        return None;
    }
    let aug_n = trial.boundary_len();
    let gap_start = n - pm.len;

    // Map preserved-outer source positions to trial positions.
    // After add_tile, the survivor strip on the trial occupies
    // positions [0, gap_start] with CCW anchor at 0 and CW anchor
    // at gap_start. Source survivor vertex `v` lives on the trial
    // at `(v - (start_a + match_len) + n) % n`.
    let ccw_anchor_src = (pm.start_a + pm.len) % n;
    let j_outer_cw_trial = (j_outer_cw_src + n - ccw_anchor_src) % n;
    let j_outer_ccw_trial = (j_outer_ccw_src + n - ccw_anchor_src) % n;
    debug_assert!(
        j_outer_cw_trial <= gap_start,
        "j_outer_cw_trial={} > gap_start={}",
        j_outer_cw_trial,
        gap_start
    );
    debug_assert!(
        j_outer_ccw_trial <= gap_start,
        "j_outer_ccw_trial={} > gap_start={}",
        j_outer_ccw_trial,
        gap_start
    );

    // Walk junctions CCW from j_outer_cw_trial to j_outer_ccw_trial.
    // We expect exactly 4 junctions: J_outer_cw, new_cw_anchor,
    // new_ccw_anchor, J_outer_ccw → 3 segs in the RHS.
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
        // Unexpected walk length — investigate rather than silently skip.
        let src_juncs_dbg: Vec<(usize, CoarseJunction)> = juncs
            .iter()
            .map(|&i| (i, patch.coarse_junction_at(i).unwrap()))
            .collect();
        panic!(
            "expected 4 junctions in RHS walk, got {}\n\
             walk: {:?}\n\
             pm: {:?}\n\
             source n: {}, aug_n: {}, gap_start: {}\n\
             source juncs ({}): {:?}\n\
             source seg_ids: {:?}\n\
             start_seg={}, end_seg={}, j_outer_cw_src={}, j_outer_ccw_src={}\n\
             ccw_anchor_src={}, j_outer_cw_trial={}, j_outer_ccw_trial={}",
            walk.len(),
            walk,
            pm,
            n, aug_n, gap_start,
            juncs.len(), src_juncs_dbg,
            seg_ids,
            start_seg, end_seg, j_outer_cw_src, j_outer_ccw_src,
            ccw_anchor_src, j_outer_cw_trial, j_outer_ccw_trial,
        );
    }
    let rhs = [
        seg_bfs.lookup_seg(&walk[0], &walk[1])?,
        seg_bfs.lookup_seg(&walk[1], &walk[2])?,
        seg_bfs.lookup_seg(&walk[2], &walk[3])?,
    ];

    Some((lhs, meta, rhs))
}

// Silence unused-import warning when not needed.
#[allow(dead_code)]
fn _suppress(_p: Arc<()>) {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic::ZZ12;
    use crate::intgeom::rat::Rat;
    use crate::intgeom::tiles;
    use crate::intgeom::tileset::TileSet;
    use std::sync::OnceLock;

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

    /// Cached SegSeqBFS for hex — running it once across all tests
    /// avoids the ~25s rebuild.
    fn hex_seq_bfs() -> &'static SegSeqBFS<ZZ12> {
        static BFS: OnceLock<SegSeqBFS<ZZ12>> = OnceLock::new();
        BFS.get_or_init(|| {
            SegSeqBFS::run(hex_tileset(), 1_000_000, 1_000_000).unwrap()
        })
    }

    #[test]
    fn hex_rewrite_rules_extracted() {
        let seq_bfs = hex_seq_bfs();
        let (table, collisions) = SegRewriteTable::build_with_stats(seq_bfs);
        eprintln!(
            "hex rewrites: distinct_lhs={}, total_rules={}, collisions={}, bfs_transitions={}",
            table.num_lhs(),
            table.num_rules(),
            collisions,
            seq_bfs.num_transitions()
        );
        assert!(table.num_rules() > 0);
        assert_eq!(
            seq_bfs.num_transitions(),
            table.num_rules() + collisions,
            "BFS transitions must equal rules + collisions"
        );
        // Sanity: RHS is always 3.
        for metas in table.by_lhs.values() {
            for rhs in metas.values() {
                assert_eq!(rhs.len(), 3);
            }
        }
    }

    /// Witness-independence is verified inside `build_with_stats` via
    /// an `assert_eq` on every (LHS, MatchMeta) collision. The
    /// collision count being >0 confirms the test actually
    /// exercises that path (= multiple witnesses extract the same
    /// rule and agree on the RHS).
    #[test]
    fn hex_witness_independence_exercised() {
        let (_table, collisions) = SegRewriteTable::build_with_stats(hex_seq_bfs());
        assert!(
            collisions > 0,
            "expected at least one (LHS, MatchMeta) collision \
             to exercise witness-independence; got {}",
            collisions
        );
        eprintln!("hex witness-independence collisions: {}", collisions);
    }

    /// Verifies that rules emerging from non-contributing witnesses
    /// are genuinely new `(LHS, MatchMeta)` pairs not present on any
    /// contributing witness — i.e. same LHS reached, but with a
    /// MatchMeta that no contributing witness's match enumerated.
    /// Also confirms witness-independence still holds (= when a
    /// contributing witness also has the rule, the RHS agrees).
    #[test]
    fn hex_extra_rules_are_genuine() {
        let seq_bfs = hex_seq_bfs();
        let (full, _) = SegRewriteTable::build_with_stats(seq_bfs);
        let (contrib_only, _) = SegRewriteTable::build_filtered(seq_bfs, |id| {
            seq_bfs.patch_contributed(id)
        });
        let mut extra: Vec<(&Vec<usize>, &MatchMeta, &[usize; 3])> = Vec::new();
        for (lhs, metas) in &full.by_lhs {
            let contrib_metas = contrib_only.by_lhs.get(lhs);
            for (meta, rhs) in metas {
                let already = contrib_metas
                    .map(|m| m.contains_key(meta))
                    .unwrap_or(false);
                if !already {
                    extra.push((lhs, meta, rhs));
                }
            }
        }
        eprintln!(
            "hex: full_rules={}, contrib_rules={}, extra={}",
            full.num_rules(),
            contrib_only.num_rules(),
            extra.len()
        );
        // Every extra rule's LHS must already exist in BFS sequences
        // (= same LHS, new MatchMeta, not a new LHS).
        let bfs_set: rustc_hash::FxHashSet<&Vec<usize>> =
            seq_bfs.sequences().iter().collect();
        for (lhs, _, _) in &extra {
            assert!(
                bfs_set.contains(*lhs),
                "extra rule has LHS not in BFS sequences: {:?}",
                lhs
            );
        }
        assert_eq!(
            full.num_rules() - contrib_only.num_rules(),
            extra.len(),
            "rule-count delta should equal the set of extra (LHS, MatchMeta) tuples"
        );
    }

    #[test]
    #[ignore]
    fn spectre_extra_rules_are_genuine() {
        let seq_bfs = SegSeqBFS::run(spectre_tileset(), 10_000_000, 10_000_000).unwrap();
        let (full, _) = SegRewriteTable::build_with_stats(&seq_bfs);
        let (contrib_only, _) = SegRewriteTable::build_filtered(&seq_bfs, |id| {
            seq_bfs.patch_contributed(id)
        });
        let mut extra: Vec<(&Vec<usize>, &MatchMeta, &[usize; 3])> = Vec::new();
        for (lhs, metas) in &full.by_lhs {
            let contrib_metas = contrib_only.by_lhs.get(lhs);
            for (meta, rhs) in metas {
                let already = contrib_metas
                    .map(|m| m.contains_key(meta))
                    .unwrap_or(false);
                if !already {
                    extra.push((lhs, meta, rhs));
                }
            }
        }
        eprintln!(
            "spectre: full_rules={}, contrib_rules={}, extra={}",
            full.num_rules(),
            contrib_only.num_rules(),
            extra.len()
        );
        extra.sort_by_key(|(lhs, _, _)| lhs.clone());
        eprintln!("first 5 extra (lhs, meta) tuples:");
        for (lhs, meta, _) in extra.iter().take(5) {
            eprintln!("  lhs={:?} meta={:?}", lhs, meta);
        }
        let bfs_set: rustc_hash::FxHashSet<&Vec<usize>> =
            seq_bfs.sequences().iter().collect();
        for (lhs, _, _) in &extra {
            assert!(
                bfs_set.contains(*lhs),
                "extra rule has LHS not in BFS sequences: {:?}",
                lhs
            );
        }
        assert_eq!(
            full.num_rules() - contrib_only.num_rules(),
            extra.len(),
            "rule-count delta should equal the set of extra (LHS, MatchMeta) tuples"
        );
    }

    #[test]
    #[ignore]
    fn spectre_rewrite_rules_extracted() {
        let seq_bfs = SegSeqBFS::run(spectre_tileset(), 10_000_000, 10_000_000).unwrap();
        let (table, collisions) = SegRewriteTable::build_with_stats(&seq_bfs);
        eprintln!(
            "spectre: sequences(bfs)={}, distinct_lhs(rules)={}, total_rules={}, collisions={}, bfs_transitions={}",
            seq_bfs.num_sequences(),
            table.num_lhs(),
            table.num_rules(),
            collisions,
            seq_bfs.num_transitions()
        );
        assert_eq!(
            seq_bfs.num_transitions(),
            table.num_rules() + collisions,
            "BFS transitions must equal rules + collisions"
        );
        if seq_bfs.num_sequences() != table.num_lhs() {
            let bfs_set: rustc_hash::FxHashSet<&Vec<usize>> =
                seq_bfs.sequences().iter().collect();
            let mut missing_from_bfs: Vec<&Vec<usize>> = table
                .by_lhs
                .keys()
                .filter(|k| !bfs_set.contains(k))
                .collect();
            missing_from_bfs.sort();
            eprintln!(
                "table has {} LHS not in BFS sequences. First 5:",
                missing_from_bfs.len()
            );
            for lhs in missing_from_bfs.iter().take(5) {
                eprintln!("  {:?}", lhs);
            }
            let table_set: rustc_hash::FxHashSet<&Vec<usize>> =
                table.by_lhs.keys().collect();
            let mut missing_from_table: Vec<&Vec<usize>> = seq_bfs
                .sequences()
                .iter()
                .filter(|s| !table_set.contains(s))
                .collect();
            missing_from_table.sort();
            eprintln!(
                "BFS has {} sequences not in table LHS. First 5:",
                missing_from_table.len()
            );
            for seq in missing_from_table.iter().take(5) {
                eprintln!("  {:?}", seq);
            }
            panic!("BFS sequences must match table LHS count");
        }
        for metas in table.by_lhs.values() {
            for rhs in metas.values() {
                assert_eq!(rhs.len(), 3);
            }
        }
    }

    #[test]
    #[ignore]
    fn mixed_rewrite_rules_extracted() {
        let seq_bfs = SegSeqBFS::run(mixed_tileset(), 20_000_000, 20_000_000).unwrap();
        let (table, collisions) = SegRewriteTable::build_with_stats(&seq_bfs);
        eprintln!(
            "mixed: sequences(bfs)={}, distinct_lhs(rules)={}, total_rules={}, collisions={}, bfs_transitions={}",
            seq_bfs.num_sequences(),
            table.num_lhs(),
            table.num_rules(),
            collisions,
            seq_bfs.num_transitions()
        );
        assert_eq!(seq_bfs.num_sequences(), table.num_lhs(),
            "BFS sequences must match table LHS count");
        assert_eq!(
            seq_bfs.num_transitions(),
            table.num_rules() + collisions,
            "BFS transitions must equal rules + collisions"
        );
        for metas in table.by_lhs.values() {
            for rhs in metas.values() {
                assert_eq!(rhs.len(), 3);
            }
        }
    }
}
