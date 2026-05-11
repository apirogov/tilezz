use std::collections::VecDeque;
use std::sync::Arc;

use rustc_hash::FxHashMap;

use crate::cyclotomic::{IsComplex, IsRingOrField, Units};
use crate::intgeom::matchtypes::MatchTypeIndex;
use crate::intgeom::patch::{EdgeInfo, GrowingPatch, PatchMatch, VertexType};
use crate::intgeom::tileset::TileSet;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct NeighborhoodType {
    pub central_tile_id: usize,
    pub cw_anchor_on_central: usize,
    pub cw_anchor_on_context: usize,
    pub vt_seq: Vec<VertexType>,
}

pub const NT_CLOSED_ID: usize = 0;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum NtKind {
    Dead,
    Undead,
    Blessed,
    Free,
}

#[derive(Clone, Debug)]
pub struct NtTransition {
    pub src_id: usize,
    pub dst_id: usize,
    pub tile_id: usize,
    pub tile_offset: usize,
}

pub struct NeighborhoodIndex<T: IsComplex> {
    tileset: Arc<TileSet<T>>,
    entries: Vec<NeighborhoodType>,
    transitions: Vec<NtTransition>,
}

impl<T: IsComplex + IsRingOrField + Units> NeighborhoodIndex<T> {
    pub fn new(tileset: Arc<TileSet<T>>) -> Self {
        let match_index = Arc::new(MatchTypeIndex::new(Arc::clone(&tileset)));
        let mut entries: Vec<NeighborhoodType> = Vec::new();
        let mut transitions: Vec<NtTransition> = Vec::new();
        // IDs are 1-based so NT_CLOSED_ID = 0 is unambiguous: entries[i] has
        // id i + 1. `seen` and transitions store ids, not indices.
        let mut seen: FxHashMap<NeighborhoodType, usize> = FxHashMap::default();
        let mut queue: VecDeque<NeighborhoodType> = VecDeque::new();

        for id in 1..=match_index.num_types() {
            let mt = match_index.get(id);
            let mut patch = GrowingPatch::new(Arc::clone(match_index.tileset()), mt.tile_a);
            let pm = PatchMatch {
                start_a: mt.start_a,
                len: mt.len,
                start_b: mt.start_b,
                tile_id: mt.tile_b,
            };
            if patch.add_tile(&pm).is_none() {
                continue;
            }
            patch.normalize();
            let patch_n = patch.boundary_len();

            // Collect all ctx-ctx junctions with their positions, sorted by
            // boundary position. We then filter per third_pm to only those
            // incident with the match (= within the covered segment).
            let all_juncs: Vec<(usize, VertexType)> = (0..patch_n)
                .filter_map(|i| patch.junction_vertex_type_at(i).map(|vt| (i, vt)))
                .collect();

            for third_pm in patch.get_all_matches() {
                let central_tile_id = third_pm.tile_id;
                let cw_anchor_on_central = third_pm.start_b;

                // vt_seq is the junctions incident with the central match:
                // those at positions in [anchor, anchor + match_len], i.e.
                // CCW distance from the anchor is at most match_len.
                let filtered: Vec<&(usize, VertexType)> = all_juncs
                    .iter()
                    .filter(|(pos, _)| (pos + patch_n - third_pm.start_a) % patch_n <= third_pm.len)
                    .collect();
                if filtered.is_empty() {
                    continue;
                }
                let first_junc = filtered[0].0;
                let vt_seq: Vec<VertexType> = filtered.iter().map(|(_, vt)| vt.clone()).collect();
                // cw_anchor_on_context = CW distance from first_junc back to
                // the anchor (equivalently, CCW distance from anchor to
                // first_junc). The first VT sits at or CCW of the anchor, so
                // this distance is small and invariant under BFS growth
                // (which only adds edges CCW of the last VT, never between
                // the anchor and first_junc).
                let cw_anchor_on_context = (first_junc + patch_n - third_pm.start_a) % patch_n;
                let nt = NeighborhoodType {
                    central_tile_id,
                    cw_anchor_on_central,
                    cw_anchor_on_context,
                    vt_seq,
                };
                if seen.contains_key(&nt) {
                    continue;
                }
                let id = entries.len() + 1;
                seen.insert(nt.clone(), id);
                entries.push(nt.clone());
                queue.push_back(nt);
            }
        }

        let mut explored = 0usize;
        while let Some(state) = queue.pop_front() {
            explored += 1;
            let src_id = *seen.get(&state).expect("dequeued state must be in seen");
            for outcome in explore_one(&state, &match_index) {
                let dst_id = match outcome.kind {
                    OutcomeKind::Closed => NT_CLOSED_ID,
                    OutcomeKind::Open { nt: new_nt, .. } => {
                        if let Some(&id) = seen.get(&new_nt) {
                            id
                        } else {
                            let id = entries.len() + 1;
                            seen.insert(new_nt.clone(), id);
                            entries.push(new_nt.clone());
                            queue.push_back(new_nt);
                            id
                        }
                    }
                };
                transitions.push(NtTransition {
                    src_id,
                    dst_id,
                    tile_id: outcome.petal_pm.tile_id,
                    tile_offset: outcome.petal_pm.start_b,
                });
            }
        }

        eprintln!(
            "  total: {} types, {} transitions, {} explored",
            entries.len(),
            transitions.len(),
            explored,
        );

        NeighborhoodIndex {
            tileset,
            entries,
            transitions,
        }
    }

    pub fn entries(&self) -> &[NeighborhoodType] {
        &self.entries
    }

    pub fn num_types(&self) -> usize {
        self.entries.len()
    }

    pub fn transitions(&self) -> &[NtTransition] {
        &self.transitions
    }

    pub fn tileset(&self) -> &Arc<TileSet<T>> {
        &self.tileset
    }

    pub fn classify_all(&self) -> Vec<NtKind> {
        // Transition ids are 1-based; convert with `- 1` when indexing into
        // these per-entry vectors. NT_CLOSED_ID (= 0) is the closed sentinel
        // and never appears as a src_id.
        let n = self.entries.len();
        let mut succ_sets: Vec<Vec<usize>> = vec![vec![]; n];
        let mut has_outgoing = vec![false; n];

        for t in &self.transitions {
            debug_assert!(t.src_id != NT_CLOSED_ID, "closed cannot be a src");
            let src_idx = t.src_id - 1;
            has_outgoing[src_idx] = true;
            if t.dst_id != NT_CLOSED_ID {
                succ_sets[src_idx].push(t.dst_id - 1);
            }
        }

        let mut cursed = vec![false; n];
        for i in 0..n {
            cursed[i] = !has_outgoing[i];
        }

        let mut changed = true;
        while changed {
            changed = false;
            for i in 0..n {
                if cursed[i] {
                    continue;
                }
                let succs = &succ_sets[i];
                if succs.is_empty() {
                    continue;
                }
                if succs.iter().all(|&s| cursed[s]) {
                    cursed[i] = true;
                    changed = true;
                }
            }
        }

        let mut blessed = vec![false; n];
        let mut has_closed_or_blessed_succ = vec![false; n];

        for t in &self.transitions {
            let src_idx = t.src_id - 1;
            if t.dst_id == NT_CLOSED_ID {
                has_closed_or_blessed_succ[src_idx] = true;
            }
        }

        changed = true;
        while changed {
            changed = false;
            for i in 0..n {
                if blessed[i] || cursed[i] {
                    continue;
                }
                let mut any = false;
                let src_id_i = i + 1;
                let all_good = self
                    .transitions
                    .iter()
                    .filter(|t| t.src_id == src_id_i)
                    .all(|t| {
                        any = true;
                        t.dst_id == NT_CLOSED_ID || blessed[t.dst_id - 1]
                    });
                if any && all_good {
                    blessed[i] = true;
                    has_closed_or_blessed_succ[i] = true;
                    changed = true;
                }
            }
        }

        self.entries
            .iter()
            .enumerate()
            .map(|(i, _)| {
                if cursed[i] && !has_outgoing[i] {
                    NtKind::Dead
                } else if cursed[i] {
                    NtKind::Undead
                } else if blessed[i] {
                    NtKind::Blessed
                } else {
                    NtKind::Free
                }
            })
            .collect()
    }

    pub fn validate(&self) -> Vec<usize> {
        todo!("validate")
    }

    pub fn write_collection(&self, _out: &mut impl std::io::Write) -> std::io::Result<()> {
        todo!("write_collection")
    }

    pub fn parse_file(_tileset: Arc<TileSet<T>>, _input: &str) -> Result<Self, String> {
        todo!("parse_file")
    }
}

#[allow(dead_code)]
struct AugmentedContext<T: IsComplex> {
    augmented: GrowingPatch<T>,
    gap_start: usize,
    gap_len: usize,
    central_ptid: usize,
}

#[allow(dead_code)]
struct FrontierInfo {
    dist_to_frontier: usize,
}

#[allow(dead_code)]
struct ExploreOutcome<T: IsComplex> {
    petal_pm: PatchMatch,
    kind: OutcomeKind<T>,
}

#[allow(dead_code, clippy::large_enum_variant)]
enum OutcomeKind<T: IsComplex> {
    Closed,
    Open {
        nt: NeighborhoodType,
        trial: GrowingPatch<T>,
        central_ptid: usize,
    },
}

#[allow(dead_code)]
fn attach_central<T: IsComplex + IsRingOrField + Units>(
    nt: &NeighborhoodType,
    match_index: &Arc<MatchTypeIndex<T>>,
) -> Option<AugmentedContext<T>> {
    let (mut context, junc_positions) =
        GrowingPatch::construct_witness_from_vt_sequence(&nt.vt_seq, Arc::clone(match_index))?;
    let first_junc = junc_positions[0];
    let ctx_n = context.boundary_len();
    let anchor_pos = (first_junc + ctx_n - nt.cw_anchor_on_context) % ctx_n;

    let pm = PatchMatch {
        start_a: anchor_pos,
        len: {
            let tileset = match_index.tileset();
            let central_seq = tileset.rat(nt.central_tile_id).seq();
            crate::intgeom::patch::forward_match_length(
                context.angles(),
                anchor_pos,
                central_seq,
                nt.cw_anchor_on_central,
            )
        },
        start_b: nt.cw_anchor_on_central,
        tile_id: nt.central_tile_id,
    };

    let tileset = match_index.tileset();
    let central_len = tileset.rat(nt.central_tile_id).seq().len();
    let ctx_n = context.boundary_len();

    context.add_tile(&pm)?;

    let aug_n = context.boundary_len();
    let gap_len = central_len - pm.len;
    let seg_len_old = ctx_n - pm.len;
    let gap_start = seg_len_old % aug_n;
    let central_ptid = context.patch_tile_ids()[gap_start];

    Some(AugmentedContext {
        augmented: context,
        gap_start,
        gap_len,
        central_ptid,
    })
}

#[allow(dead_code)]
fn find_gap_frontier<T: IsComplex + IsRingOrField + Units>(
    aug: &AugmentedContext<T>,
) -> Option<FrontierInfo> {
    let patch = &aug.augmented;
    let n = patch.boundary_len();

    let gap_end = (aug.gap_start + aug.gap_len) % n;

    let mut frontier_pos = gap_end;
    for _ in 0..n {
        if patch.is_junction(frontier_pos) {
            break;
        }
        frontier_pos = (frontier_pos + 1) % n;
    }
    if !patch.is_junction(frontier_pos) {
        return None;
    }

    let mut dist_to_frontier = aug.gap_start;
    for i in (0..aug.gap_start).rev() {
        if patch.is_junction(i) {
            dist_to_frontier = aug.gap_start - i - 1;
            break;
        }
    }

    Some(FrontierInfo { dist_to_frontier })
}

#[allow(dead_code)]
fn explore_step<T: IsComplex + IsRingOrField + Units>(
    nt: &NeighborhoodType,
    match_index: &Arc<MatchTypeIndex<T>>,
) -> Option<(AugmentedContext<T>, FrontierInfo, Vec<ExploreOutcome<T>>)> {
    let aug = attach_central(nt, match_index)?;
    let frontier = find_gap_frontier(&aug)?;
    let patch = &aug.augmented;
    let n = patch.boundary_len();
    let central_ptid = aug.central_ptid;

    let gap_end = (aug.gap_start + aug.gap_len) % n;
    let mut frontier_pos = gap_end;
    for _ in 0..n {
        if patch.is_junction(frontier_pos) {
            break;
        }
        frontier_pos = (frontier_pos + 1) % n;
    }

    let candidates = GrowingPatch::compute_candidates_covering_position(
        patch.match_index(),
        patch.angles(),
        patch.edges(),
        frontier_pos,
    );

    let mut outcomes = Vec::new();
    for petal_pm in &candidates {
        let mut trial = aug.augmented.clone_for_mutation();
        if trial.add_tile(petal_pm).is_none() {
            continue;
        }

        let kind = if find_remaining_gap(&trial, central_ptid).is_none() {
            OutcomeKind::Closed
        } else {
            OutcomeKind::Open {
                nt: nt.clone(),
                trial,
                central_ptid,
            }
        };
        outcomes.push(ExploreOutcome {
            petal_pm: petal_pm.clone(),
            kind,
        });
    }

    Some((aug, frontier, outcomes))
}

#[allow(dead_code)]
fn find_remaining_gap<T: IsComplex + IsRingOrField + Units>(
    trial: &GrowingPatch<T>,
    central_ptid: usize,
) -> Option<(usize, usize)> {
    let n = trial.boundary_len();
    let ptids = trial.patch_tile_ids();
    let first = (0..n).find(|&i| ptids[i] == central_ptid)?;
    let mut len = 1;
    while len < n && ptids[(first + len) % n] == central_ptid {
        len += 1;
    }
    Some((first, len))
}

#[allow(dead_code)]
struct AttachedContext<T: IsComplex> {
    aug: GrowingPatch<T>,
    central_ptid: usize,
    frontier_pos_on_aug: usize,
    frontier_is_junction_in_ctx: bool,
    last_covered_ctx_edge: EdgeInfo,
}

#[allow(dead_code)]
fn build_attached_context<T: IsComplex + IsRingOrField + Units>(
    nt: &NeighborhoodType,
    match_index: &Arc<MatchTypeIndex<T>>,
) -> Option<AttachedContext<T>> {
    let (ctx, junc_positions) =
        GrowingPatch::construct_witness_from_vt_sequence(&nt.vt_seq, Arc::clone(match_index))?;
    let ctx_n = ctx.boundary_len();
    if ctx_n == 0 {
        return None;
    }
    let first_junc = junc_positions[0];
    let anchor_pos = (first_junc + ctx_n - nt.cw_anchor_on_context) % ctx_n;

    let tileset = match_index.tileset();
    let central_seq = tileset.rat(nt.central_tile_id).seq();
    let central_n = central_seq.len();
    let match_len = crate::intgeom::patch::forward_match_length(
        ctx.angles(),
        anchor_pos,
        central_seq,
        nt.cw_anchor_on_central,
    );
    if match_len == 0 || match_len >= central_n {
        return None;
    }

    let frontier_on_ctx = (anchor_pos + match_len) % ctx_n;
    let frontier_is_junction_in_ctx = ctx.is_junction(frontier_on_ctx);
    let last_covered_ctx_edge = ctx.edges()[(anchor_pos + match_len - 1) % ctx_n];

    let central_pm = PatchMatch {
        start_a: anchor_pos,
        len: match_len,
        start_b: nt.cw_anchor_on_central,
        tile_id: nt.central_tile_id,
    };
    let mut aug = ctx.clone_for_mutation();
    aug.add_tile(&central_pm)?;

    let aug_n = aug.boundary_len();
    let gap_start = ctx_n - match_len;
    let central_ptid = aug.patch_tile_ids()[gap_start];

    // gap_end == 0 by construction; advance CCW until a junction is found.
    let mut frontier_pos_on_aug = 0usize;
    let mut found = false;
    for _ in 0..aug_n {
        if aug.is_junction(frontier_pos_on_aug) {
            found = true;
            break;
        }
        frontier_pos_on_aug = (frontier_pos_on_aug + 1) % aug_n;
    }
    if !found {
        return None;
    }

    Some(AttachedContext {
        aug,
        central_ptid,
        frontier_pos_on_aug,
        frontier_is_junction_in_ctx,
        last_covered_ctx_edge,
    })
}

#[allow(dead_code)]
fn explore_one<T: IsComplex + IsRingOrField + Units>(
    nt: &NeighborhoodType,
    match_index: &Arc<MatchTypeIndex<T>>,
) -> Vec<ExploreOutcome<T>> {
    let ac = match build_attached_context(nt, match_index) {
        Some(a) => a,
        None => return Vec::new(),
    };

    let candidates = GrowingPatch::compute_candidates_covering_position(
        ac.aug.match_index(),
        ac.aug.angles(),
        ac.aug.edges(),
        ac.frontier_pos_on_aug,
    );

    let mut results = Vec::new();
    for petal_pm in candidates {
        let mut trial = ac.aug.clone_for_mutation();
        if trial.add_tile(&petal_pm).is_none() {
            continue;
        }

        // Closed: petal absorbed all central edges.
        if !trial.patch_tile_ids().contains(&ac.central_ptid) {
            results.push(ExploreOutcome {
                petal_pm,
                kind: OutcomeKind::Closed,
            });
            continue;
        }

        // The first surviving petal edge has tile_offset == petal_pm.start_b
        // (see add_tile_growing in patch.rs). That edge starts at the OLD
        // frontier vertex on ctx, and is the ccw side of the new VT.
        let petal_edge = EdgeInfo {
            tile_id: petal_pm.tile_id,
            tile_offset: petal_pm.start_b,
        };

        let mut new_vt_seq = nt.vt_seq.clone();
        if ac.frontier_is_junction_in_ctx {
            // Frontier vertex was already a ctx-ctx junction (last VT in vt_seq).
            // The petal becomes the new ccw-most tile at this junction; the old
            // ccw tile is pushed into `inner`.
            let last = new_vt_seq.last_mut().expect("vt_seq is non-empty");
            let old_ccw = last.ccw;
            last.inner.push(old_ccw);
            last.ccw = petal_edge;
        } else {
            // Frontier vertex was not a ctx-ctx junction; the petal creates a
            // new junction there. cw side is the last covered ctx edge.
            new_vt_seq.push(VertexType {
                cw: ac.last_covered_ctx_edge,
                inner: vec![],
                ccw: petal_edge,
            });
        }

        let new_nt = NeighborhoodType {
            central_tile_id: nt.central_tile_id,
            cw_anchor_on_central: nt.cw_anchor_on_central,
            cw_anchor_on_context: nt.cw_anchor_on_context,
            vt_seq: new_vt_seq,
        };

        results.push(ExploreOutcome {
            petal_pm,
            kind: OutcomeKind::Open {
                nt: new_nt,
                trial,
                central_ptid: ac.central_ptid,
            },
        });
    }
    results
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic::ZZ12;
    use crate::intgeom::matchtypes::MatchTypeIndex;
    use crate::intgeom::patch::GrowingPatch;
    use crate::intgeom::rat::Rat;
    use crate::intgeom::tiles;

    fn square_tileset() -> Arc<TileSet<ZZ12>> {
        let sq = tiles::square::<ZZ12>();
        let rat = Rat::try_from(&sq).unwrap();
        Arc::new(TileSet::new(vec![rat]))
    }

    fn hex_tileset() -> Arc<TileSet<ZZ12>> {
        let hex = tiles::hexagon::<ZZ12>();
        let rat = Rat::try_from(&hex).unwrap();
        Arc::new(TileSet::new(vec![rat]))
    }

    #[test]
    fn square_seed_count() {
        let idx = NeighborhoodIndex::new(square_tileset());
        assert!(idx.num_types() > 0, "expected non-empty seed collection");
        for nt in idx.entries() {
            assert_eq!(nt.central_tile_id, 0, "single-tile tileset");
            assert!(!nt.vt_seq.is_empty(), "seeds must have non-empty vt_seq");
        }
    }

    #[test]
    fn hex_seed_count() {
        let idx = NeighborhoodIndex::new(hex_tileset());
        assert!(idx.num_types() > 0, "expected non-empty seed collection");
        for nt in idx.entries() {
            assert_eq!(nt.central_tile_id, 0, "single-tile tileset");
            assert!(!nt.vt_seq.is_empty(), "seeds must have non-empty vt_seq");
        }
    }

    #[test]
    fn seeds_have_no_duplicate_keys() {
        let idx = NeighborhoodIndex::new(hex_tileset());
        let mut keys = std::collections::HashSet::new();
        for (i, nt) in idx.entries().iter().enumerate() {
            assert!(
                keys.insert((
                    nt.central_tile_id,
                    nt.cw_anchor_on_central,
                    nt.cw_anchor_on_context,
                    nt.vt_seq.clone(),
                )),
                "[{}] duplicate seed key: central={}, anchor={}, ctx_anchor={}, vt_seq={:?}",
                i,
                nt.central_tile_id,
                nt.cw_anchor_on_central,
                nt.cw_anchor_on_context,
                nt.vt_seq,
            );
        }
    }

    /// Group NTs by their reconstructed context's normalized boundary
    /// (angle-sequence rotated to lex-min). Many NTs share the same context
    /// shape (different anchors / different central-tile orientations against
    /// the same context). For square / hex, the number of distinct shapes
    /// should be small even though the NT count is large.
    #[test]
    fn context_shape_diversity() {
        for (name, idx) in [
            ("SQUARE", NeighborhoodIndex::new(square_tileset())),
            ("HEX", NeighborhoodIndex::new(hex_tileset())),
        ] {
            let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
            let mut by_shape: std::collections::HashMap<Vec<i8>, usize> =
                std::collections::HashMap::new();
            let mut by_tile_count: std::collections::BTreeMap<usize, usize> =
                std::collections::BTreeMap::new();
            for nt in idx.entries() {
                let (ctx, _) =
                    GrowingPatch::construct_witness_from_vt_sequence(&nt.vt_seq, Arc::clone(&mi))
                        .expect("reconstruct");
                let angles = ctx.angles();
                let rot = crate::intgeom::rat::lex_min_rot(angles);
                let mut norm = angles.to_vec();
                norm.rotate_left(rot);
                *by_shape.entry(norm).or_insert(0) += 1;
                let n_tiles = ctx
                    .patch_tile_ids()
                    .iter()
                    .collect::<std::collections::HashSet<_>>()
                    .len();
                *by_tile_count.entry(n_tiles).or_insert(0) += 1;
            }
            eprintln!(
                "{}: {} NTs -> {} distinct context shapes",
                name,
                idx.num_types(),
                by_shape.len()
            );
            for (n, count) in &by_tile_count {
                eprintln!("  {} ctx tiles: {} NTs", n, count);
            }
            // Print a sample of the smallest shapes (by boundary length).
            let mut shapes: Vec<(&Vec<i8>, &usize)> = by_shape.iter().collect();
            shapes.sort_by_key(|(angles, _)| angles.len());
            eprintln!("  smallest shapes (by boundary length):");
            for (angles, count) in shapes.iter().take(8) {
                eprintln!(
                    "    len={} count={} angles={:?}",
                    angles.len(),
                    count,
                    angles
                );
            }
        }
    }

    #[test]
    fn hex_level4_outcomes() {
        let idx = NeighborhoodIndex::new(hex_tileset());
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        // Find a few len-4 NTs and check what explore_one returns.
        let len4: Vec<&NeighborhoodType> = idx
            .entries()
            .iter()
            .filter(|n| n.vt_seq.len() == 4)
            .take(5)
            .collect();
        assert!(!len4.is_empty(), "expected some len-4 NTs");
        for (i, nt) in len4.iter().enumerate() {
            let outcomes = explore_one(nt, &mi);
            let n_open = outcomes
                .iter()
                .filter(|o| matches!(o.kind, OutcomeKind::Open { .. }))
                .count();
            let n_closed = outcomes
                .iter()
                .filter(|o| matches!(o.kind, OutcomeKind::Closed))
                .count();
            eprintln!(
                "len-4 nt[{}]: ac={} ax={} -> {} open, {} closed outcomes",
                i, nt.cw_anchor_on_central, nt.cw_anchor_on_context, n_open, n_closed
            );
        }
    }

    #[test]
    fn bfs_vt_seq_length_distribution() {
        let sq = NeighborhoodIndex::new(square_tileset());
        let mut hist: std::collections::BTreeMap<usize, usize> = std::collections::BTreeMap::new();
        let mut max_inner = 0;
        for nt in sq.entries() {
            *hist.entry(nt.vt_seq.len()).or_insert(0) += 1;
            for vt in &nt.vt_seq {
                max_inner = max_inner.max(vt.inner.len());
            }
        }
        eprintln!(
            "SQUARE: entries={} transitions={} max_inner_per_vt={}",
            sq.num_types(),
            sq.transitions().len(),
            max_inner
        );
        for (len, count) in &hist {
            eprintln!("  vt_seq_len={}: {} entries", len, count);
        }

        let hex = NeighborhoodIndex::new(hex_tileset());
        let mut hist: std::collections::BTreeMap<usize, usize> = std::collections::BTreeMap::new();
        let mut max_inner = 0;
        for nt in hex.entries() {
            *hist.entry(nt.vt_seq.len()).or_insert(0) += 1;
            for vt in &nt.vt_seq {
                max_inner = max_inner.max(vt.inner.len());
            }
        }
        eprintln!(
            "HEX: entries={} transitions={} max_inner_per_vt={}",
            hex.num_types(),
            hex.transitions().len(),
            max_inner
        );
        for (len, count) in &hist {
            eprintln!("  vt_seq_len={}: {} entries", len, count);
        }
    }

    #[test]
    fn bfs_produces_transitions() {
        let sq_idx = NeighborhoodIndex::new(square_tileset());
        assert!(sq_idx.num_types() > 0);
        assert!(
            !sq_idx.transitions().is_empty(),
            "BFS should produce transitions"
        );
        // src_id is a 1-based entry id; dst_id is either 1-based or
        // NT_CLOSED_ID (= 0).
        for t in sq_idx.transitions() {
            assert!(t.src_id >= 1 && t.src_id <= sq_idx.num_types());
            assert!(t.dst_id == NT_CLOSED_ID || t.dst_id <= sq_idx.num_types());
        }

        let hex_idx = NeighborhoodIndex::new(hex_tileset());
        assert!(hex_idx.num_types() > 0);
        assert!(!hex_idx.transitions().is_empty());
        for t in hex_idx.transitions() {
            assert!(t.src_id >= 1 && t.src_id <= hex_idx.num_types());
            assert!(t.dst_id == NT_CLOSED_ID || t.dst_id <= hex_idx.num_types());
        }
    }

    fn validate_seeds<T: IsComplex + IsRingOrField + Units>(
        idx: &NeighborhoodIndex<T>,
    ) -> Vec<String> {
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        let mut errors = Vec::new();
        for nt in idx.entries() {
            let Some((mut ctx, junc_positions)) =
                GrowingPatch::construct_witness_from_vt_sequence(&nt.vt_seq, Arc::clone(&mi))
            else {
                errors.push(format!(
                    "nt c={} ac={} ax={}: reconstruct failed",
                    nt.central_tile_id, nt.cw_anchor_on_central, nt.cw_anchor_on_context
                ));
                continue;
            };
            let ctx_n = ctx.boundary_len();
            let first_junc = junc_positions[0];
            let anchor_pos = (first_junc + ctx_n - nt.cw_anchor_on_context) % ctx_n;
            let found = ctx.get_all_matches().into_iter().find(|pm| {
                pm.tile_id == nt.central_tile_id
                    && pm.start_a == anchor_pos
                    && pm.start_b == nt.cw_anchor_on_central
            });
            let Some(pm) = found else {
                errors.push(format!(
                    "nt c={} ac={} ax={}: no match at anchor positions",
                    nt.central_tile_id, nt.cw_anchor_on_central, nt.cw_anchor_on_context
                ));
                continue;
            };
            if ctx.add_tile(&pm).is_none() {
                errors.push(format!(
                    "nt c={} ac={} ax={}: add_tile failed",
                    nt.central_tile_id, nt.cw_anchor_on_central, nt.cw_anchor_on_context
                ));
            }
        }
        errors
    }

    #[test]
    fn square_seeds_validate() {
        let idx = NeighborhoodIndex::new(square_tileset());
        let errors = validate_seeds(&idx);
        assert!(
            errors.is_empty(),
            "validation errors:\n{}",
            errors.join("\n")
        );
    }

    #[test]
    fn hex_seeds_validate() {
        let idx = NeighborhoodIndex::new(hex_tileset());
        let errors = validate_seeds(&idx);
        assert!(
            errors.is_empty(),
            "validation errors:\n{}",
            errors.join("\n")
        );
    }

    #[test]
    fn attach_central_square_seeds() {
        let idx = NeighborhoodIndex::new(square_tileset());
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        for (i, nt) in idx.entries().iter().enumerate() {
            let aug = attach_central(nt, &mi)
                .unwrap_or_else(|| panic!("seed {}: attach_central failed for nt {:?}", i, nt));
            assert!(aug.gap_len > 0, "seed {}: gap_len should be > 0", i);
            assert!(
                aug.gap_start < aug.augmented.boundary_len(),
                "seed {}: gap_start out of range",
                i
            );
            let tileset = mi.tileset();
            let central_len = tileset.rat(nt.central_tile_id).seq().len();
            assert!(
                aug.gap_len < central_len,
                "seed {}: gap_len {} should be < central_len {}",
                i,
                aug.gap_len,
                central_len
            );
        }
    }

    #[test]
    fn attach_central_hex_seeds() {
        let idx = NeighborhoodIndex::new(hex_tileset());
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        for (i, nt) in idx.entries().iter().enumerate() {
            let aug = attach_central(nt, &mi)
                .unwrap_or_else(|| panic!("seed {}: attach_central failed for nt {:?}", i, nt));
            assert!(aug.gap_len > 0, "seed {}: gap_len should be > 0", i);
        }
    }

    #[test]
    fn find_gap_frontier_square_seeds() {
        let idx = NeighborhoodIndex::new(square_tileset());
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        for (i, nt) in idx.entries().iter().enumerate() {
            let (aug, frontier, petal_outcomes) =
                explore_step(nt, &mi).unwrap_or_else(|| panic!("seed {}: explore_step failed", i));
            assert!(
                frontier.dist_to_frontier < aug.augmented.boundary_len(),
                "seed {}: dist_to_frontier out of range",
                i
            );
            assert!(
                !petal_outcomes.is_empty(),
                "seed {}: should have at least one successful petal",
                i
            );
        }
    }

    #[test]
    fn find_gap_frontier_hex_seeds() {
        let idx = NeighborhoodIndex::new(hex_tileset());
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        for (i, nt) in idx.entries().iter().enumerate() {
            let (aug, frontier, petal_outcomes) =
                explore_step(nt, &mi).unwrap_or_else(|| panic!("seed {}: explore_step failed", i));
            assert!(
                frontier.dist_to_frontier < aug.augmented.boundary_len(),
                "seed {}: dist_to_frontier out of range",
                i
            );
            assert!(
                !petal_outcomes.is_empty(),
                "seed {}: should have at least one successful petal",
                i
            );
        }
    }

    #[test]
    fn gap_edges_are_central_tile() {
        let idx = NeighborhoodIndex::new(square_tileset());
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        for (i, nt) in idx.entries().iter().enumerate() {
            let aug = attach_central(nt, &mi)
                .unwrap_or_else(|| panic!("seed {}: attach_central failed", i));
            let edges = aug.augmented.edges();
            let n = aug.augmented.boundary_len();
            let ptids = aug.augmented.patch_tile_ids();
            let central_ptid = ptids[aug.gap_start];
            for k in 0..aug.gap_len {
                let pos = (aug.gap_start + k) % n;
                assert_eq!(
                    edges[pos].tile_id, nt.central_tile_id,
                    "seed {}: gap edge at pos {} should be central tile",
                    i, pos
                );
                assert_eq!(
                    ptids[pos], central_ptid,
                    "seed {}: gap edge at pos {} should have same ptid",
                    i, pos
                );
            }
        }
    }

    #[test]
    fn frontier_adjacent_to_gap() {
        let idx = NeighborhoodIndex::new(hex_tileset());
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        for (i, nt) in idx.entries().iter().enumerate() {
            let (aug, _frontier, _petal_outcomes) =
                explore_step(nt, &mi).unwrap_or_else(|| panic!("seed {}: explore_step failed", i));
            let n = aug.augmented.boundary_len();
            let gap_end = (aug.gap_start + aug.gap_len) % n;
            assert!(
                aug.augmented.is_junction(gap_end),
                "seed {}: gap_end {} should be a junction",
                i,
                gap_end
            );
        }
    }

    #[test]
    fn explore_step_square_has_petals() {
        let idx = NeighborhoodIndex::new(square_tileset());
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        let nt = &idx.entries()[0];
        let (_aug, _frontier, petal_outcomes) =
            explore_step(nt, &mi).expect("explore_step on first seed");
        assert!(
            !petal_outcomes.is_empty(),
            "should have at least one petal outcome"
        );
    }

    #[test]
    fn find_remaining_gap_after_petal() {
        let idx = NeighborhoodIndex::new(square_tileset());
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        for (i, nt) in idx.entries().iter().enumerate() {
            let (aug, _frontier, petal_outcomes) =
                explore_step(nt, &mi).unwrap_or_else(|| panic!("seed {}: explore_step failed", i));
            let central_ptid = aug.central_ptid;
            for outcome in &petal_outcomes {
                match &outcome.kind {
                    OutcomeKind::Closed => {}
                    OutcomeKind::Open {
                        trial,
                        central_ptid: cptid,
                        ..
                    } => {
                        let gap = find_remaining_gap(trial, *cptid);
                        let n = trial.boundary_len();
                        let ptids = trial.patch_tile_ids();
                        let central_count = ptids.iter().filter(|&&id| id == central_ptid).count();
                        if central_count == 0 {
                            assert!(
                                gap.is_none(),
                                "seed {}: gap should be None when fully consumed",
                                i
                            );
                        } else {
                            let (gs, gl) = gap.unwrap_or_else(|| {
                                panic!(
                                    "seed {}: gap should exist ({} central edges)",
                                    i, central_count
                                )
                            });
                            assert_eq!(gl, central_count, "seed {}: gap length mismatch", i);
                            for k in 0..gl {
                                assert_eq!(
                                    ptids[(gs + k) % n],
                                    central_ptid,
                                    "seed {}: gap edge at offset {} should be central_ptid",
                                    i,
                                    k
                                );
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn explore_one_hex_seeds_reconstruct() {
        let idx = NeighborhoodIndex::new(hex_tileset());
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        let mut checked = 0usize;
        let mut errors = Vec::new();
        for (i, nt) in idx.entries().iter().enumerate() {
            for outcome in explore_one(nt, &mi) {
                let new_nt = match &outcome.kind {
                    OutcomeKind::Closed => continue,
                    OutcomeKind::Open { nt, .. } => nt.clone(),
                };
                let (mut ctx, junc_pos) = GrowingPatch::construct_witness_from_vt_sequence(
                    &new_nt.vt_seq,
                    Arc::clone(&mi),
                )
                .unwrap_or_else(|| {
                    panic!(
                        "seed {} result: reconstruct failed for new_nt {:?}",
                        i, new_nt
                    )
                });
                let first_junc = junc_pos[0];
                let ctx_n = ctx.boundary_len();
                let anchor_pos = (first_junc + ctx_n - new_nt.cw_anchor_on_context) % ctx_n;
                let found = ctx.get_all_matches().into_iter().find(|pm| {
                    pm.tile_id == new_nt.central_tile_id
                        && pm.start_a == anchor_pos
                        && pm.start_b == new_nt.cw_anchor_on_central
                });
                let Some(pm) = found else {
                    errors.push(format!(
                        "seed {}: no match c={} ac={} ax={} (ctx_n={})",
                        i,
                        new_nt.central_tile_id,
                        new_nt.cw_anchor_on_central,
                        new_nt.cw_anchor_on_context,
                        ctx.boundary_len()
                    ));
                    continue;
                };
                if ctx.add_tile(&pm).is_none() {
                    errors.push(format!(
                        "seed {}: add_tile failed for new_nt {:?}",
                        i, new_nt
                    ));
                    continue;
                }
                checked += 1;
            }
        }
        eprintln!("validated {} new NTs from hex seeds", checked);
        assert!(checked > 0);
        if !errors.is_empty() {
            for e in &errors[..errors.len().min(10)] {
                eprintln!("  {}", e);
            }
            panic!("{} reattach errors", errors.len());
        }
    }

    #[test]
    fn dist_to_frontier_computation() {
        let idx = NeighborhoodIndex::new(square_tileset());
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        for (i, nt) in idx.entries().iter().enumerate() {
            let (aug, frontier, _petal_outcomes) =
                explore_step(nt, &mi).unwrap_or_else(|| panic!("seed {}: explore_step failed", i));
            assert!(
                frontier.dist_to_frontier < aug.augmented.boundary_len(),
                "seed {}: dist_to_frontier out of range",
                i
            );
        }
    }

    #[test]
    fn explore_one_square_seeds_reconstruct() {
        let idx = NeighborhoodIndex::new(square_tileset());
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        let mut checked = 0usize;
        let mut errors = Vec::new();
        for (i, nt) in idx.entries().iter().enumerate() {
            for outcome in explore_one(nt, &mi) {
                let new_nt = match &outcome.kind {
                    OutcomeKind::Closed => continue,
                    OutcomeKind::Open { nt, .. } => nt.clone(),
                };
                let (mut ctx, junc_pos) = GrowingPatch::construct_witness_from_vt_sequence(
                    &new_nt.vt_seq,
                    Arc::clone(&mi),
                )
                .unwrap_or_else(|| {
                    panic!(
                        "seed {} result: reconstruct failed for new_nt {:?}",
                        i, new_nt
                    )
                });
                let first_junc = junc_pos[0];
                let ctx_n = ctx.boundary_len();
                let anchor_pos = (first_junc + ctx_n - new_nt.cw_anchor_on_context) % ctx_n;
                let found = ctx.get_all_matches().into_iter().find(|pm| {
                    pm.tile_id == new_nt.central_tile_id
                        && pm.start_a == anchor_pos
                        && pm.start_b == new_nt.cw_anchor_on_central
                });
                let Some(pm) = found else {
                    errors.push(format!(
                        "seed {}: no match c={} ac={} ax={} (ctx_n={})",
                        i,
                        new_nt.central_tile_id,
                        new_nt.cw_anchor_on_central,
                        new_nt.cw_anchor_on_context,
                        ctx.boundary_len()
                    ));
                    continue;
                };
                if ctx.add_tile(&pm).is_none() {
                    errors.push(format!(
                        "seed {}: add_tile failed for new_nt {:?}",
                        i, new_nt
                    ));
                    continue;
                }
                checked += 1;
            }
        }
        eprintln!("validated {} new NTs from square seeds", checked);
        assert!(checked > 0);
        if !errors.is_empty() {
            for e in &errors[..errors.len().min(10)] {
                eprintln!("  {}", e);
            }
            panic!("{} reattach errors", errors.len());
        }
    }
}
