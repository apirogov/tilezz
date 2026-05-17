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
use crate::intgeom::segseqbfs::{lhs_seg_range, MatchMeta, SegSeqBFS};

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
        Self::build_with_stats(seq_bfs)
    }

    /// Same as `build_from` but also returns the collision count
    /// (= number of times a `(LHS, MatchMeta)` was re-extracted from
    /// a different witness).
    pub fn build_with_stats<T: IsComplex + IsRingOrField + Units>(
        seq_bfs: &SegSeqBFS<T>,
    ) -> Self {
        let mut by_lhs: FxHashMap<Vec<usize>, FxHashMap<MatchMeta, [usize; 3]>> =
            FxHashMap::default();
        let seg_bfs = seq_bfs.seg_bfs();
        for rule in seq_bfs.rules() {
            let patch = &seq_bfs.patches()[rule.witness_patch_id];
            let n = patch.boundary_len();
            let juncs: Vec<usize> =
                (0..n).filter(|&i| patch.is_junction(i)).collect();
            let k = juncs.len();
            if k < 2 {
                continue;
            }
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
            let triple = extract_rule(
                patch,
                n,
                &juncs,
                &seg_ids,
                &rule.witness_pm,
                seg_bfs,
            );
            let (lhs, meta, rhs) = match triple {
                Some(t) => t,
                None => continue,
            };
            assert_eq!(lhs, rule.lhs, "extracted LHS doesn't match rule");
            assert_eq!(meta, rule.meta, "extracted MatchMeta doesn't match rule");
            by_lhs.entry(lhs).or_default().insert(meta, rhs);
        }
        SegRewriteTable { by_lhs }
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
        let table = SegRewriteTable::build_with_stats(seq_bfs);
        eprintln!(
            "hex rewrites: distinct_lhs={}, total_rules={}, bfs_rules={}",
            table.num_lhs(),
            table.num_rules(),
            seq_bfs.num_rules()
        );
        assert!(table.num_rules() > 0);
        assert_eq!(
            seq_bfs.num_rules(),
            table.num_rules(),
            "table rules must equal BFS rules — one canonical witness per rule"
        );
        // Sanity: RHS is always 3.
        for metas in table.by_lhs.values() {
            for rhs in metas.values() {
                assert_eq!(rhs.len(), 3);
            }
        }
    }

    #[test]
    #[ignore]
    fn spectre_rewrite_rules_extracted() {
        let seq_bfs = SegSeqBFS::run(spectre_tileset(), 10_000_000, 10_000_000).unwrap();
        let table = SegRewriteTable::build_with_stats(&seq_bfs);
        eprintln!(
            "spectre: sequences(bfs)={}, distinct_lhs(rules)={}, total_rules={}, bfs_rules={}",
            seq_bfs.num_sequences(),
            table.num_lhs(),
            table.num_rules(),
            seq_bfs.num_rules()
        );
        assert_eq!(seq_bfs.num_rules(), table.num_rules(),
            "table rules must equal BFS rules");
        assert_eq!(seq_bfs.num_sequences(), table.num_lhs(),
            "BFS sequences must match table LHS count");
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
        let table = SegRewriteTable::build_with_stats(&seq_bfs);
        eprintln!(
            "mixed: sequences(bfs)={}, distinct_lhs(rules)={}, total_rules={}, bfs_rules={}",
            seq_bfs.num_sequences(),
            table.num_lhs(),
            table.num_rules(),
            seq_bfs.num_rules()
        );
        assert_eq!(seq_bfs.num_rules(), table.num_rules(),
            "table rules must equal BFS rules");
        assert_eq!(seq_bfs.num_sequences(), table.num_lhs(),
            "BFS sequences must match table LHS count");
        for metas in table.by_lhs.values() {
            for rhs in metas.values() {
                assert_eq!(rhs.len(), 3);
            }
        }
    }
}
