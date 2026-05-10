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
        let transitions: Vec<NtTransition> = Vec::new();
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
            let vt_seq = extract_context_vt_seq(&patch);
            if vt_seq.is_empty() {
                continue;
            }

            for third_pm in patch.get_all_matches() {
                let central_tile_id = third_pm.tile_id;
                let cw_anchor_on_central = third_pm.start_b;
                let cw_anchor_on_context = third_pm.start_a;
                let nt = NeighborhoodType {
                    central_tile_id,
                    cw_anchor_on_central,
                    cw_anchor_on_context,
                    vt_seq: vt_seq.clone(),
                };
                if seen.contains_key(&nt) {
                    continue;
                }
                let idx = entries.len();
                seen.insert(nt.clone(), idx);
                entries.push(nt.clone());
                queue.push_back(nt);
            }
        }

        let mut explored = 0usize;
        while let Some(state) = queue.pop_front() {
            explored += 1;
            let _state_id = seen.get(&state).expect("dequeued state must be in seen");
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
        let n = self.entries.len();
        let mut succ_sets: Vec<Vec<usize>> = vec![vec![]; n];
        let mut has_outgoing = vec![false; n];

        for t in &self.transitions {
            if t.src_id == NT_CLOSED_ID {
                continue;
            }
            let src_idx = t.src_id;
            has_outgoing[src_idx] = true;
            if t.dst_id != NT_CLOSED_ID {
                succ_sets[src_idx].push(t.dst_id);
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
            if t.src_id == NT_CLOSED_ID {
                continue;
            }
            let src_idx = t.src_id;
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
                let all_good = self.transitions.iter().filter(|t| t.src_id == i).all(|t| {
                    any = true;
                    t.dst_id == NT_CLOSED_ID || blessed[t.dst_id]
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

fn extract_context_vt_seq<T: IsComplex + IsRingOrField + Units>(
    patch: &GrowingPatch<T>,
) -> Vec<VertexType> {
    let n = patch.edges().len();
    if n == 0 {
        return Vec::new();
    }
    let mut vts = Vec::new();
    for pos in 0..n {
        if patch.is_junction(pos) {
            if let Some(vt) = patch.junction_vertex_type_at(pos) {
                vts.push(vt);
            }
        }
    }
    vts
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
enum PetalOutcome<T: IsComplex> {
    Closed {
        tile_id: usize,
        tile_offset: usize,
    },
    Open {
        tile_id: usize,
        tile_offset: usize,
        trial: GrowingPatch<T>,
        central_ptid: usize,
    },
}

#[allow(dead_code)]
struct ExploreResult<T: IsComplex> {
    outcome: ExploreOutcome<T>,
}

#[allow(dead_code)]
enum ExploreOutcome<T: IsComplex> {
    Closed,
    Open(NeighborhoodType, #[expect(dead_code)] GrowingPatch<T>),
}

#[allow(dead_code)]
fn attach_central<T: IsComplex + IsRingOrField + Units>(
    nt: &NeighborhoodType,
    match_index: &Arc<MatchTypeIndex<T>>,
) -> Option<AugmentedContext<T>> {
    let (mut context, _) =
        GrowingPatch::construct_witness_from_vt_sequence(&nt.vt_seq, Arc::clone(match_index))?;

    let pm = PatchMatch {
        start_a: nt.cw_anchor_on_context,
        len: {
            let tileset = match_index.tileset();
            let central_seq = tileset.rat(nt.central_tile_id).seq();
            crate::intgeom::patch::forward_match_length(
                context.angles(),
                nt.cw_anchor_on_context,
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
) -> Option<(AugmentedContext<T>, FrontierInfo, Vec<PetalOutcome<T>>)> {
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

    let mut petal_outcomes = Vec::new();
    for petal_pm in &candidates {
        let mut trial = aug.augmented.clone_for_mutation();
        if trial.add_tile(petal_pm).is_none() {
            continue;
        }

        let has_gap = find_remaining_gap(&trial, central_ptid).is_some();
        if !has_gap {
            petal_outcomes.push(PetalOutcome::Closed {
                tile_id: petal_pm.tile_id,
                tile_offset: petal_pm.start_b,
            });
        } else {
            petal_outcomes.push(PetalOutcome::Open {
                tile_id: petal_pm.tile_id,
                tile_offset: petal_pm.start_b,
                trial,
                central_ptid,
            });
        }
    }

    Some((aug, frontier, petal_outcomes))
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
fn find_petal_edge_adjacent_to_gap<T: IsComplex + IsRingOrField + Units>(
    trial: &GrowingPatch<T>,
    new_gap_end: usize,
    central_ptid: usize,
    petal_ptid: usize,
) -> Option<EdgeInfo> {
    let n = trial.boundary_len();
    let ptids = trial.patch_tile_ids();
    let edges = trial.edges();
    let pos_before_gap = (new_gap_end + n - 1) % n;
    if ptids[pos_before_gap] == petal_ptid {
        Some(edges[pos_before_gap])
    } else {
        None
    }
}

#[allow(dead_code)]
fn explore_one<T: IsComplex + IsRingOrField + Units>(
    nt: &NeighborhoodType,
    match_index: &Arc<MatchTypeIndex<T>>,
) -> Vec<ExploreOutcome<T>> {
    let (aug, frontier, petal_outcomes) = match explore_step(nt, match_index) {
        Some(r) => r,
        None => return Vec::new(),
    };

    let dist_to_frontier = frontier.dist_to_frontier;

    let mut results = Vec::new();
    for outcome in petal_outcomes {
        match outcome {
            PetalOutcome::Closed { .. } => {
                results.push(ExploreOutcome::Closed);
            }
            PetalOutcome::Open {
                trial,
                central_ptid,
                ..
            } => {
                let gap = match find_remaining_gap(&trial, central_ptid) {
                    Some(g) => g,
                    None => {
                        results.push(ExploreOutcome::Closed);
                        continue;
                    }
                };

                let trial_n = trial.boundary_len();
                let new_gap_end = (gap.0 + gap.1) % trial_n;
                let petal_ptid = *trial.patch_tile_ids().iter().max().unwrap();

                let petal_edge_info = match find_petal_edge_adjacent_to_gap(
                    &trial,
                    new_gap_end,
                    central_ptid,
                    petal_ptid,
                ) {
                    Some(e) => e,
                    None => continue,
                };

                let mut new_vt_seq = nt.vt_seq.clone();
                if dist_to_frontier > 0 {
                    let new_junc_aug_pos = aug.gap_start - dist_to_frontier;
                    let cw_edge = aug.augmented.edges()[(new_junc_aug_pos
                        + aug.augmented.boundary_len()
                        - 1)
                        % aug.augmented.boundary_len()];
                    new_vt_seq.push(VertexType {
                        cw: cw_edge,
                        inner: vec![],
                        ccw: petal_edge_info,
                    });
                } else {
                    let last = new_vt_seq.last_mut().expect("vt_seq is non-empty");
                    last.inner.push(last.ccw);
                    last.ccw = petal_edge_info;
                }

                let new_nt = NeighborhoodType {
                    central_tile_id: nt.central_tile_id,
                    cw_anchor_on_central: nt.cw_anchor_on_central,
                    cw_anchor_on_context: nt.cw_anchor_on_context,
                    vt_seq: new_vt_seq,
                };

                results.push(ExploreOutcome::Open(new_nt, trial));
            }
        }
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

    #[test]
    fn bfs_explores_all_seeds() {
        let sq_idx = NeighborhoodIndex::new(square_tileset());
        assert!(sq_idx.num_types() > 0);
        assert!(sq_idx.transitions().is_empty());

        let hex_idx = NeighborhoodIndex::new(hex_tileset());
        assert!(hex_idx.num_types() > 0);
        assert!(hex_idx.transitions().is_empty());
    }

    fn validate_seeds<T: IsComplex + IsRingOrField + Units>(
        idx: &NeighborhoodIndex<T>,
    ) -> Vec<String> {
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        let mut errors = Vec::new();
        for nt in idx.entries() {
            let Some((mut ctx, _)) =
                GrowingPatch::construct_witness_from_vt_sequence(&nt.vt_seq, Arc::clone(&mi))
            else {
                errors.push(format!(
                    "nt c={} ac={} ax={}: reconstruct failed",
                    nt.central_tile_id, nt.cw_anchor_on_central, nt.cw_anchor_on_context
                ));
                continue;
            };
            let found = ctx.get_all_matches().into_iter().find(|pm| {
                pm.tile_id == nt.central_tile_id
                    && pm.start_a == nt.cw_anchor_on_context
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
                match outcome {
                    PetalOutcome::Closed { .. } => {}
                    PetalOutcome::Open {
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
    fn petal_edge_found_at_gap_boundary() {
        let idx = NeighborhoodIndex::new(square_tileset());
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        for (i, nt) in idx.entries().iter().enumerate() {
            let (aug, _frontier, petal_outcomes) =
                explore_step(nt, &mi).unwrap_or_else(|| panic!("seed {}: explore_step failed", i));
            for outcome in &petal_outcomes {
                let (trial, petal_tile_id) = match outcome {
                    PetalOutcome::Open { trial, tile_id, .. } => (trial, *tile_id),
                    _ => continue,
                };
                let central_ptid = aug.central_ptid;
                let gap = match find_remaining_gap(trial, central_ptid) {
                    Some(g) => g,
                    None => continue,
                };
                let n = trial.boundary_len();
                let gap_end = (gap.0 + gap.1) % n;
                let petal_ptid = *trial.patch_tile_ids().iter().max().unwrap();
                let petal_edge =
                    find_petal_edge_adjacent_to_gap(trial, gap_end, central_ptid, petal_ptid);
                let pei = petal_edge
                    .unwrap_or_else(|| panic!("seed {}: no petal edge adjacent to gap found", i));
                assert_eq!(
                    pei.tile_id, petal_tile_id,
                    "seed {}: petal edge tile_id mismatch",
                    i
                );
            }
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
    fn diag_correct_vt_seq_from_trial() {
        let idx = NeighborhoodIndex::new(square_tileset());
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        let nt = &idx.entries()[0];

        let (aug, _frontier, petal_outcomes) = explore_step(nt, &mi).expect("explore_step");
        let central_ptid = aug.central_ptid;
        let petal_outcome = &petal_outcomes[0];
        let (trial, petal_tile_id) = match petal_outcome {
            PetalOutcome::Open { trial, tile_id, .. } => (trial, *tile_id),
            _ => panic!("expected open"),
        };
        let trial_n = trial.boundary_len();

        // Show what explore_one produces for this candidate
        let results = explore_one(nt, &mi);
        let new_nt = match &results[0] {
            ExploreOutcome::Open(nt, _) => nt.clone(),
            _ => panic!("expected open"),
        };
        eprintln!("explore_one vt_seq:");
        for (vi, vt) in new_nt.vt_seq.iter().enumerate() {
            eprintln!(
                "  vt[{}]: cw=({},{}) inner={:?} ccw=({},{})",
                vi,
                vt.cw.tile_id,
                vt.cw.tile_offset,
                vt.inner
                    .iter()
                    .map(|e| (e.tile_id, e.tile_offset))
                    .collect::<Vec<_>>(),
                vt.ccw.tile_id,
                vt.ccw.tile_offset
            );
        }

        // Now extract the CORRECT vt_seq from the trial boundary
        // by reading junction VTs at context-only positions
        eprintln!("\ncorrect vt_seq (from trial context-only junctions):");
        for i in 0..trial_n {
            if !trial.is_junction(i) {
                continue;
            }
            let ptid_cw = trial.patch_tile_ids()[(i + trial_n - 1) % trial_n];
            let ptid_ccw = trial.patch_tile_ids()[i];
            let is_context = ptid_cw != central_ptid && ptid_ccw != central_ptid;
            if is_context {
                let vt = trial.junction_vertex_type_at(i).unwrap();
                eprintln!(
                    "  trial[{}]: cw=({},{}) inner={:?} ccw=({},{})",
                    i,
                    vt.cw.tile_id,
                    vt.cw.tile_offset,
                    vt.inner
                        .iter()
                        .map(|e| (e.tile_id, e.tile_offset))
                        .collect::<Vec<_>>(),
                    vt.ccw.tile_id,
                    vt.ccw.tile_offset
                );
            }
        }

        // Reconstruct from the correct vt_seq and check
        let mut correct_vts = Vec::new();
        for i in 0..trial_n {
            if !trial.is_junction(i) {
                continue;
            }
            let ptid_cw = trial.patch_tile_ids()[(i + trial_n - 1) % trial_n];
            let ptid_ccw = trial.patch_tile_ids()[i];
            let is_context = ptid_cw != central_ptid && ptid_ccw != central_ptid;
            if is_context {
                correct_vts.push(trial.junction_vertex_type_at(i).unwrap());
            }
        }

        let (ctx, _) =
            GrowingPatch::construct_witness_from_vt_sequence(&correct_vts, Arc::clone(&mi))
                .expect("correct vt_seq should reconstruct");
        eprintln!("\ncorrect reconstruction: ctx_n={}", ctx.boundary_len());

        // Check if central can attach
        let matches: Vec<_> = ctx
            .get_all_matches()
            .into_iter()
            .filter(|pm| pm.tile_id == nt.central_tile_id && pm.start_b == nt.cw_anchor_on_central)
            .collect();
        eprintln!(
            "matches at start_b={}: {}",
            nt.cw_anchor_on_central,
            matches.len()
        );
        for m in &matches {
            eprintln!("  start_a={} len={}", m.start_a, m.len);
        }

        // Check if cw_anchor_on_context works
        let at_anchor: Vec<_> = matches
            .iter()
            .filter(|pm| pm.start_a == nt.cw_anchor_on_context)
            .collect();
        eprintln!(
            "matches at start_a={}: {}",
            nt.cw_anchor_on_context,
            at_anchor.len()
        );
    }

    #[test]
    fn explore_one_square_seeds_reconstruct() {
        let idx = NeighborhoodIndex::new(square_tileset());
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        let mut checked = 0usize;
        let mut errors = Vec::new();
        for (i, nt) in idx.entries().iter().enumerate() {
            for outcome in explore_one(nt, &mi) {
                let new_nt = match &outcome {
                    ExploreOutcome::Closed => continue,
                    ExploreOutcome::Open(nt, _) => nt.clone(),
                };
                let (mut ctx, _) = GrowingPatch::construct_witness_from_vt_sequence(
                    &new_nt.vt_seq,
                    Arc::clone(&mi),
                )
                .unwrap_or_else(|| {
                    panic!(
                        "seed {} result: reconstruct failed for new_nt {:?}",
                        i, new_nt
                    )
                });
                let found = ctx.get_all_matches().into_iter().find(|pm| {
                    pm.tile_id == new_nt.central_tile_id
                        && pm.start_a == new_nt.cw_anchor_on_context
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
        assert!(checked > 0, "should have checked at least one result");
        if !errors.is_empty() {
            for e in &errors[..errors.len().min(10)] {
                eprintln!("  {}", e);
            }
            panic!("{} reattach errors", errors.len());
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
                let new_nt = match &outcome {
                    ExploreOutcome::Closed => continue,
                    ExploreOutcome::Open(nt, _) => nt.clone(),
                };
                let (mut ctx, _) = GrowingPatch::construct_witness_from_vt_sequence(
                    &new_nt.vt_seq,
                    Arc::clone(&mi),
                )
                .unwrap_or_else(|| {
                    panic!(
                        "seed {} result: reconstruct failed for new_nt {:?}",
                        i, new_nt
                    )
                });
                let found = ctx.get_all_matches().into_iter().find(|pm| {
                    pm.tile_id == new_nt.central_tile_id
                        && pm.start_a == new_nt.cw_anchor_on_context
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
}
