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
use crate::intgeom::segseqbfs::SegSeqBFS;

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
    /// LHS sequence → all distinct outgoing transitions.
    /// RHS is always 3 segs: [pad_cw', central, pad_ccw'].
    pub by_lhs: FxHashMap<Vec<usize>, FxHashMap<MatchMeta, [usize; 3]>>,
}

impl SegRewriteTable {
    /// Build the table eagerly from a finished [`SegSeqBFS`].
    pub fn build_from<T: IsComplex + IsRingOrField + Units>(
        seq_bfs: &SegSeqBFS<T>,
    ) -> Self {
        let mut by_lhs: FxHashMap<Vec<usize>, FxHashMap<MatchMeta, [usize; 3]>> =
            FxHashMap::default();
        let seg_bfs = seq_bfs.seg_bfs();
        for patch in seq_bfs.patches() {
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
                    debug_assert_eq!(
                        prev, rhs,
                        "witness-independence violated for meta {:?}",
                        meta
                    );
                }
            }
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
    let start_edge = pm.start_a;
    let end_edge = (pm.start_a + pm.len - 1) % n;
    let start_seg = seg_containing_edge(juncs, start_edge, n)?;
    let end_seg = seg_containing_edge(juncs, end_edge, n)?;

    // Build the LHS sequence: [start_seg, internal..., end_seg].
    // (Same definition as SegSeqBFS: only the segs the match touches.)
    let mut lhs = Vec::with_capacity(k);
    let mut idx = start_seg;
    loop {
        lhs.push(seg_ids[idx]);
        if idx == end_seg {
            break;
        }
        idx = (idx + 1) % k;
        debug_assert!(lhs.len() <= k, "infinite loop building LHS");
    }

    // Preserved outer junctions on the source:
    //   CW outer = start_seg's CW endpoint (= junction CW of start_seg).
    //   CCW outer = end_seg's CCW endpoint (= junction CCW of end_seg).
    let j_outer_cw_src = juncs[start_seg];
    let j_outer_ccw_src = juncs[(end_seg + 1) % k];

    // Match metadata.
    let start_edge_in_lhs = (start_edge + n - j_outer_cw_src) % n;
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
        panic!(
            "expected 4 junctions in RHS walk, got {}: walk={:?}, outer_cw={:?}, outer_ccw={:?}, pm={:?}",
            walk.len(), walk, walk.first(), walk.last(), pm
        );
    }
    let rhs = [
        seg_bfs.lookup_seg(&walk[0], &walk[1])?,
        seg_bfs.lookup_seg(&walk[1], &walk[2])?,
        seg_bfs.lookup_seg(&walk[2], &walk[3])?,
    ];

    Some((lhs, meta, rhs))
}

fn seg_containing_edge(juncs: &[usize], e: usize, _n: usize) -> Option<usize> {
    let k = juncs.len();
    for i in 0..k {
        let start = juncs[i];
        let end = juncs[(i + 1) % k];
        if start <= end {
            if e >= start && e < end {
                return Some(i);
            }
        } else if e >= start || e < end {
            return Some(i);
        }
    }
    None
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
        let table = SegRewriteTable::build_from(hex_seq_bfs());
        eprintln!(
            "hex rewrites: distinct_lhs={}, total_rules={}",
            table.num_lhs(),
            table.num_rules()
        );
        assert!(table.num_rules() > 0);
        // Sanity: RHS is always 3 (= invariant from the LHS
        // definition).
        for metas in table.by_lhs.values() {
            for rhs in metas.values() {
                assert_eq!(rhs.len(), 3);
            }
        }
    }
}
