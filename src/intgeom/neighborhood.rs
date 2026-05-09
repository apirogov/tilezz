use std::collections::{BTreeSet, VecDeque};
use std::sync::Arc;

use rustc_hash::FxHashMap;

use crate::cyclotomic::{IsComplex, IsRingOrField, Units};
use crate::intgeom::matchtypes::MatchTypeIndex;
use crate::intgeom::patch::{
    compute_glue_angles, forward_match_length, glue_match_to_raw_boundary,
    glue_tile_to_raw_boundary, vertex_type_raw_from, EdgeInfo, GrowingPatch, PatchMatch,
    RawBoundary, VertexType,
};
use crate::intgeom::rat::Rat;
use crate::intgeom::snake::Snake;
use crate::intgeom::tileset::TileSet;

const MAX_BOUNDARY_SIZE: usize = 80;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NeighborhoodType {
    pub id: usize,
    pub central_tile_id: usize,
    pub central_anchor_edge: usize,
    pub gap_len: u8,
    pub context_vertex_types: Vec<VertexType>,
    pub angles: Vec<i8>,
    pub edges: Vec<EdgeInfo>,
    pub inner_chains: Vec<Vec<EdgeInfo>>,
    pub gap_start: usize,
}

impl NeighborhoodType {
    pub fn gap_len(&self) -> usize {
        self.gap_len as usize
    }

    pub fn validate<T: IsComplex + IsRingOrField + Units>(
        &self,
        tileset: &TileSet<T>,
    ) -> Result<(), String> {
        let n = self.angles.len();
        if n == 0 {
            return Err("empty boundary".into());
        }
        if self.edges.len() != n || self.inner_chains.len() != n {
            return Err(format!(
                "length mismatch: angles={} edges={} inner={}",
                n,
                self.edges.len(),
                self.inner_chains.len()
            ));
        }
        if self.gap_len() >= n {
            return Err(format!("gap_len {} >= boundary_len {}", self.gap_len(), n));
        }
        if self.gap_len() == 0 {
            return Err("gap_len is zero".into());
        }

        let tile_rat = tileset.rat(self.central_tile_id);
        let tile_seq = tile_rat.seq();
        let m = tile_seq.len();

        for i in 0..self.gap_len() {
            let pos = (self.gap_start + i) % n;
            let expected_offset = (self.central_anchor_edge + i) % m;
            let edge = &self.edges[pos];
            if edge.tile_id != self.central_tile_id {
                return Err(format!(
                    "gap edge at pos={} has tile_id={} expected={}",
                    pos, edge.tile_id, self.central_tile_id
                ));
            }
            if edge.tile_offset != expected_offset {
                return Err(format!(
                    "gap edge at pos={} has tile_offset={} expected={}",
                    pos, edge.tile_offset, expected_offset
                ));
            }
        }

        Ok(())
    }
}

fn in_consumed_range(pos: usize, start_a: usize, match_len: usize, n: usize) -> bool {
    if match_len == 0 || n == 0 {
        return false;
    }
    if match_len >= n {
        return true;
    }
    let end = (start_a + match_len) % n;
    if end == 0 {
        pos >= start_a
    } else if start_a < end {
        pos >= start_a && pos < end
    } else {
        pos >= start_a || pos < end
    }
}

fn remap_surviving_position(
    old_pos: usize,
    start_a: usize,
    match_len: usize,
    old_n: usize,
) -> Option<usize> {
    if in_consumed_range(old_pos, start_a, match_len, old_n) {
        return None;
    }
    let ccw_pos = (start_a + match_len) % old_n;
    let old_survivor_len = old_n - match_len;
    let new_pos = (old_pos + old_n - ccw_pos) % old_n;
    if new_pos < old_survivor_len {
        Some(new_pos)
    } else {
        None
    }
}

fn extract_covered_vertex_types(
    edges: &[EdgeInfo],
    inner_chains: &[Vec<EdgeInfo>],
    gap_start: usize,
    gap_len: usize,
) -> Vec<VertexType> {
    let n = edges.len();
    let covered_len = n - gap_len;
    if covered_len < 2 {
        return Vec::new();
    }
    let num_inner = covered_len - 1;
    let mut result = Vec::with_capacity(num_inner);
    for k in 1..=num_inner {
        let pos = (gap_start + n - k) % n;
        result.push(vertex_type_raw_from(edges, inner_chains, pos));
    }
    result
}

fn tile_boundary<T: IsComplex + IsRingOrField + Units>(
    tileset: &TileSet<T>,
    tile_id: usize,
) -> RawBoundary {
    let seq = tileset.rat(tile_id).seq();
    RawBoundary {
        angles: seq.to_vec(),
        edges: (0..seq.len())
            .map(|tile_offset| EdgeInfo {
                tile_id,
                tile_offset,
            })
            .collect(),
        inner_chains: vec![vec![]; seq.len()],
        patch_tile_ids: vec![0; seq.len()],
    }
}

fn raw_is_junction<T: IsComplex + IsRingOrField + Units>(
    boundary: &RawBoundary,
    tileset: &TileSet<T>,
    pos: usize,
) -> bool {
    let edge = boundary.edges[pos];
    tileset.rat(edge.tile_id).seq()[edge.tile_offset] != boundary.angles[pos]
}

fn validate_raw_boundary<T: IsComplex + IsRingOrField + Units>(boundary: &RawBoundary) -> bool {
    if boundary.angles.iter().any(|a| a.abs() == T::hturn()) {
        return false;
    }
    Snake::<T>::try_from(boundary.angles.as_slice())
        .map(|snake| snake.is_closed())
        .unwrap_or(false)
}

#[derive(Clone)]
struct NtBfsState {
    central_tile_id: usize,
    cw_anchor_on_central: usize,
    vt_seq: Vec<VertexType>,
    seed_tile_id: usize,
    cw_junc_on_seed: usize,
    seed_match_len: usize,
    cw_junc_on_ctx: usize,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct NtStateKey {
    central_tile_id: usize,
    cw_anchor_on_central: usize,
    vt_seq: Vec<VertexType>,
}

impl NtBfsState {
    fn key(&self) -> NtStateKey {
        NtStateKey {
            central_tile_id: self.central_tile_id,
            cw_anchor_on_central: self.cw_anchor_on_central,
            vt_seq: self.vt_seq.clone(),
        }
    }
}

struct AttachResult {
    augmented: RawBoundary,
    covered_len: usize,
    gap_len: usize,
    gap_start: usize,
    frontier: usize,
}

fn frontier_distance(ctx_n: usize, cw_junc: usize, ctx_frontier: usize) -> usize {
    (cw_junc + ctx_n - ctx_frontier) % ctx_n
}

fn reconstruct_context<T: IsComplex + IsRingOrField + Units>(
    state: &NtBfsState,
    match_index: &Arc<MatchTypeIndex<T>>,
) -> Option<(RawBoundary, usize)> {
    if state.vt_seq.is_empty() {
        return Some((
            tile_boundary(match_index.tileset().as_ref(), state.seed_tile_id),
            state.cw_junc_on_ctx,
        ));
    }

    let (patch, _) =
        GrowingPatch::construct_witness_from_vt_sequence(&state.vt_seq, Arc::clone(match_index))?;

    Some((
        RawBoundary {
            angles: patch.angles().to_vec(),
            edges: patch.edges().to_vec(),
            inner_chains: patch.inner_chains().to_vec(),
            patch_tile_ids: patch.patch_tile_ids().to_vec(),
        },
        state.cw_junc_on_ctx,
    ))
}

fn find_gap_frontier<T: IsComplex + IsRingOrField + Units>(
    augmented: &RawBoundary,
    tileset: &TileSet<T>,
    gap_start: usize,
    gap_len: usize,
) -> usize {
    let n = augmented.angles.len();
    for step in 1..=gap_len {
        let pos = (gap_start + step) % n;
        if raw_is_junction(augmented, tileset, pos) {
            return pos;
        }
    }
    (gap_start + gap_len) % n
}

fn attach_central<T: IsComplex + IsRingOrField + Units>(
    context: &RawBoundary,
    cw_junc_on_ctx: usize,
    central_tile_id: usize,
    cw_anchor_on_central: usize,
    tileset: &Arc<TileSet<T>>,
    first_step: bool,
) -> Option<Option<AttachResult>> {
    let ctx_n = context.angles.len();
    if ctx_n == 0 || cw_junc_on_ctx >= ctx_n {
        return None;
    }

    let central_rat = tileset.rat(central_tile_id);
    let central_n = central_rat.len();
    let ctx_rat = Rat::from_slice_unchecked(&context.angles);
    let (ctx_match_start, covered_len, central_match_end) = ctx_rat.get_match(
        (cw_junc_on_ctx as i64, cw_anchor_on_central as i64),
        central_rat,
    );

    if covered_len == central_n {
        return Some(None);
    }

    if ctx_match_start.rem_euclid(ctx_n as i64) as usize != cw_junc_on_ctx
        || central_match_end.rem_euclid(central_n as i64) as usize != cw_anchor_on_central
        || covered_len == 0
        || covered_len > central_n
    {
        return None;
    }

    let gap_len = central_n - covered_len;

    let pm = PatchMatch {
        start_a: cw_junc_on_ctx,
        len: covered_len,
        start_b: cw_anchor_on_central,
        tile_id: central_tile_id,
    };
    let glue = glue_match_to_raw_boundary::<T>(context, &pm, tileset.as_ref(), first_step, 0)?;
    if !validate_raw_boundary::<T>(&glue.boundary) {
        return None;
    }

    let gap_start = glue.old_survivor_len;
    let frontier = find_gap_frontier(&glue.boundary, tileset.as_ref(), gap_start, gap_len);
    Some(Some(AttachResult {
        augmented: glue.boundary,
        covered_len,
        gap_len,
        gap_start,
        frontier,
    }))
}

fn make_neighborhood_entry(
    id: usize,
    state: &NtBfsState,
    attachment: AttachResult,
) -> Option<NeighborhoodType> {
    let gap_len = u8::try_from(attachment.gap_len).ok()?;
    Some(NeighborhoodType {
        id,
        central_tile_id: state.central_tile_id,
        central_anchor_edge: state.cw_anchor_on_central,
        gap_len,
        context_vertex_types: state.vt_seq.clone(),
        angles: attachment.augmented.angles,
        edges: attachment.augmented.edges,
        inner_chains: attachment.augmented.inner_chains,
        gap_start: attachment.gap_start,
    })
}

fn valid_two_tile_glue<T: IsComplex + IsRingOrField + Units>(
    seed_angles: &[i8],
    pm: &PatchMatch,
    tileset: &Arc<TileSet<T>>,
) -> bool {
    compute_glue_angles::<T>(seed_angles, pm, tileset)
        .and_then(|angles| Snake::<T>::try_from(angles.as_slice()).ok())
        .map(|snake| snake.is_closed())
        .unwrap_or(false)
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
    pub match_start: usize,
    pub match_len: usize,
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
        let mut seen: FxHashMap<NtStateKey, usize> = FxHashMap::default();
        let mut queue: VecDeque<NtBfsState> = VecDeque::new();
        let mut seen_seed_glues = BTreeSet::new();

        for seed_tile_id in 0..tileset.num_tiles() {
            let seed_rat = tileset.rat(seed_tile_id);
            assert!(seed_rat.len() <= 32, "tile edge count exceeds 32");
            let seed_angles = seed_rat.seq().to_vec();
            let seed = GrowingPatch::new(Arc::clone(&tileset), seed_tile_id);

            for pm in seed.get_all_matches() {
                if !valid_two_tile_glue::<T>(&seed_angles, &pm, &tileset) {
                    continue;
                }

                let other_tile_id = pm.tile_id;
                let oriented = (seed_tile_id, pm.start_a, other_tile_id, pm.start_b, pm.len);
                let reversed = (other_tile_id, pm.start_b, seed_tile_id, pm.start_a, pm.len);
                if !seen_seed_glues.insert(oriented.min(reversed)) {
                    continue;
                }

                let seed_specs = [
                    (seed_tile_id, other_tile_id, pm.start_a, pm.start_b),
                    (other_tile_id, seed_tile_id, pm.start_b, pm.start_a),
                ];

                for (central_tile_id, context_tile_id, central_anchor, context_anchor) in seed_specs
                {
                    if pm.len >= tileset.rat(central_tile_id).len() {
                        continue;
                    }

                    let state = NtBfsState {
                        central_tile_id,
                        cw_anchor_on_central: central_anchor,
                        vt_seq: Vec::new(),
                        seed_tile_id: context_tile_id,
                        cw_junc_on_seed: context_anchor,
                        seed_match_len: pm.len,
                        cw_junc_on_ctx: context_anchor,
                    };

                    let context = tile_boundary(tileset.as_ref(), context_tile_id);
                    let attachment = match attach_central(
                        &context,
                        context_anchor,
                        central_tile_id,
                        central_anchor,
                        &tileset,
                        true,
                    ) {
                        Some(Some(a)) => a,
                        _ => continue,
                    };

                    let key = state.key();
                    if seen.contains_key(&key) {
                        continue;
                    }

                    let id = entries.len() + 1;
                    let Some(entry) = make_neighborhood_entry(id, &state, attachment) else {
                        continue;
                    };
                    seen.insert(key, id);
                    entries.push(entry);
                    queue.push_back(state);
                }
            }
        }

        let mut explored = 0usize;
        while let Some(state) = queue.pop_front() {
            explored += 1;
            if explored.is_multiple_of(5000) {
                eprintln!(
                    "  explored {}, queue {}, types {}, transitions {}",
                    explored,
                    queue.len(),
                    entries.len(),
                    transitions.len(),
                );
            }

            let state_id = *seen
                .get(&state.key())
                .expect("dequeued state must be in seen");

            let Some((context, cw_junc_on_ctx)) = reconstruct_context(&state, &match_index) else {
                continue;
            };

            let Some(Some(source_att)) = attach_central(
                &context,
                cw_junc_on_ctx,
                state.central_tile_id,
                state.cw_anchor_on_central,
                &tileset,
                state.vt_seq.is_empty(),
            ) else {
                continue;
            };

            if source_att.augmented.angles.len() > MAX_BOUNDARY_SIZE {
                continue;
            }

            let ctx_n = context.angles.len();
            let ctx_frontier = (cw_junc_on_ctx + source_att.covered_len) % ctx_n;
            let src_frontier_dist = frontier_distance(ctx_n, cw_junc_on_ctx, ctx_frontier);

            for tile_id in 0..tileset.num_tiles() {
                let tile_rat = tileset.rat(tile_id);
                let tile_seq = tile_rat.seq();
                let tile_len = tile_rat.len();
                for tile_offset in 0..tile_len {
                    let fml = forward_match_length(
                        &source_att.augmented.angles,
                        source_att.frontier,
                        tile_seq,
                        tile_offset,
                    );
                    if fml == 0 {
                        continue;
                    }

                    let target = EdgeInfo {
                        tile_id,
                        tile_offset,
                    };

                    let Some(aug_glue) = glue_tile_to_raw_boundary::<T>(
                        &source_att.augmented,
                        source_att.frontier,
                        target,
                        tileset.as_ref(),
                        false,
                        0,
                    ) else {
                        continue;
                    };
                    if !validate_raw_boundary::<T>(&aug_glue.boundary) {
                        continue;
                    }

                    if !aug_glue
                        .boundary
                        .edges
                        .iter()
                        .any(|e| e.tile_id == state.central_tile_id)
                    {
                        transitions.push(NtTransition {
                            src_id: state_id,
                            dst_id: NT_CLOSED_ID,
                            tile_id,
                            tile_offset,
                            match_start: source_att.frontier,
                            match_len: aug_glue.match_len,
                        });
                        continue;
                    }

                    let Some(replay_glue) = glue_tile_to_raw_boundary::<T>(
                        &context,
                        ctx_frontier,
                        target,
                        tileset.as_ref(),
                        state.vt_seq.is_empty(),
                        0,
                    ) else {
                        continue;
                    };
                    if replay_glue.boundary.angles.len() > MAX_BOUNDARY_SIZE
                        || !validate_raw_boundary::<T>(&replay_glue.boundary)
                    {
                        continue;
                    }

                    let Some(new_cw_junc) = remap_surviving_position(
                        cw_junc_on_ctx,
                        ctx_frontier,
                        replay_glue.match_len,
                        ctx_n,
                    ) else {
                        continue;
                    };

                    let new_context = replay_glue.boundary;

                    let new_att = match attach_central(
                        &new_context,
                        new_cw_junc,
                        state.central_tile_id,
                        state.cw_anchor_on_central,
                        &tileset,
                        false,
                    ) {
                        Some(Some(a)) => a,
                        Some(None) => {
                            transitions.push(NtTransition {
                                src_id: state_id,
                                dst_id: NT_CLOSED_ID,
                                tile_id,
                                tile_offset,
                                match_start: source_att.frontier,
                                match_len: aug_glue.match_len,
                            });
                            continue;
                        }
                        None => continue,
                    };

                    if new_att.covered_len < source_att.covered_len {
                        continue;
                    }

                    let new_frontier_vt = vertex_type_raw_from(
                        &new_att.augmented.edges,
                        &new_att.augmented.inner_chains,
                        new_att.frontier,
                    );

                    let new_ctx_n = new_context.angles.len();
                    let new_ctx_frontier = (new_cw_junc + new_att.covered_len) % new_ctx_n;
                    let new_frontier_dist =
                        frontier_distance(new_ctx_n, new_cw_junc, new_ctx_frontier);

                    let mut new_vt_seq = state.vt_seq.clone();
                    if new_vt_seq.is_empty() {
                        new_vt_seq.push(new_frontier_vt);
                    } else if new_frontier_dist == src_frontier_dist {
                        *new_vt_seq.last_mut().expect("non-empty") = new_frontier_vt;
                    } else {
                        new_vt_seq.push(new_frontier_vt);
                    }
                    if new_vt_seq.len() > MAX_BOUNDARY_SIZE {
                        continue;
                    }

                    let new_state = NtBfsState {
                        central_tile_id: state.central_tile_id,
                        cw_anchor_on_central: state.cw_anchor_on_central,
                        vt_seq: new_vt_seq,
                        seed_tile_id: state.seed_tile_id,
                        cw_junc_on_seed: state.cw_junc_on_seed,
                        seed_match_len: state.seed_match_len,
                        cw_junc_on_ctx: new_cw_junc,
                    };

                    let key = new_state.key();
                    let dst_id = if let Some(&id) = seen.get(&key) {
                        id
                    } else {
                        let id = entries.len() + 1;
                        let Some(entry) = make_neighborhood_entry(id, &new_state, new_att) else {
                            continue;
                        };
                        seen.insert(key, id);
                        entries.push(entry);
                        queue.push_back(new_state);
                        id
                    };

                    transitions.push(NtTransition {
                        src_id: state_id,
                        dst_id,
                        tile_id,
                        tile_offset,
                        match_start: source_att.frontier,
                        match_len: aug_glue.match_len,
                    });
                }
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

    pub fn classify_all(&self) -> Vec<NtKind> {
        let n = self.entries.len();
        let mut succ_sets: Vec<Vec<usize>> = vec![vec![]; n];
        let mut has_outgoing = vec![false; n];

        for t in &self.transitions {
            if t.src_id == NT_CLOSED_ID {
                continue;
            }
            let src_idx = t.src_id - 1;
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
                if succs.iter().all(|&s| cursed[s - 1]) {
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
                let id = i + 1;
                let mut any = false;
                let all_good = self.transitions.iter().filter(|t| t.src_id == id).all(|t| {
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

    pub fn tileset(&self) -> &Arc<TileSet<T>> {
        &self.tileset
    }

    pub fn validate(&self) -> Vec<usize> {
        let match_index = Arc::new(MatchTypeIndex::new(Arc::clone(&self.tileset)));
        let mut invalid = Vec::new();
        for nhood in &self.entries {
            if !self.validate_one(nhood, &match_index) {
                invalid.push(nhood.id);
            }
        }
        invalid
    }

    fn validate_one(&self, nhood: &NeighborhoodType, match_index: &Arc<MatchTypeIndex<T>>) -> bool {
        if nhood.gap_len == 0 {
            return false;
        }
        if nhood.angles.len() != nhood.edges.len() || nhood.angles.is_empty() {
            return false;
        }
        if nhood.inner_chains.len() != nhood.angles.len() {
            return false;
        }

        let patch = match GrowingPatch::from_parts(
            Arc::clone(match_index),
            nhood.angles.clone(),
            nhood.edges.clone(),
            nhood.inner_chains.clone(),
            vec![0; nhood.angles.len()],
            1,
        ) {
            Some(p) => p,
            None => return false,
        };

        if patch.boundary_len() != nhood.angles.len() {
            return false;
        }
        if patch.edges().len() != nhood.edges.len() {
            return false;
        }
        for (i, (a, b)) in patch.angles().iter().zip(nhood.angles.iter()).enumerate() {
            if a != b {
                eprintln!("  angle mismatch at {}: {} vs {}", i, a, b);
                return false;
            }
        }
        for (i, (a, b)) in patch.edges().iter().zip(nhood.edges.iter()).enumerate() {
            if a != b {
                eprintln!(
                    "  edge mismatch at {}: {}.{} vs {}.{}",
                    i, a.tile_id, a.tile_offset, b.tile_id, b.tile_offset
                );
                return false;
            }
        }

        if let Err(err) = nhood.validate(self.tileset.as_ref()) {
            eprintln!("  neighborhood validation failed: {err}");
            return false;
        }

        let gs = nhood.gap_start;

        if nhood.edges[gs].tile_offset != nhood.central_anchor_edge {
            eprintln!(
                "  central_anchor_edge mismatch: found {} vs claimed {}",
                nhood.edges[gs].tile_offset, nhood.central_anchor_edge
            );
            return false;
        }

        true
    }

    pub fn write_collection(&self, out: &mut impl std::io::Write) -> std::io::Result<()> {
        for nhood in &self.entries {
            let n = nhood.angles.len();
            let angles_str = nhood
                .angles
                .iter()
                .map(|a| a.to_string())
                .collect::<Vec<_>>()
                .join(" ");
            let edges_str = nhood
                .edges
                .iter()
                .map(|e| format!("{}.{}", e.tile_id, e.tile_offset))
                .collect::<Vec<_>>()
                .join(" ");
            let inner_chains_str = nhood
                .inner_chains
                .iter()
                .map(|chain| {
                    if chain.is_empty() {
                        "-".to_string()
                    } else {
                        chain
                            .iter()
                            .map(|e| format!("{}.{}", e.tile_id, e.tile_offset))
                            .collect::<Vec<_>>()
                            .join(",")
                    }
                })
                .collect::<Vec<_>>()
                .join(" ");
            let vts_str = if nhood.context_vertex_types.is_empty() {
                "-".to_string()
            } else {
                nhood
                    .context_vertex_types
                    .iter()
                    .map(format_vtype)
                    .collect::<Vec<_>>()
                    .join(" ")
            };

            writeln!(
                out,
                "NTYPE {} {} {} {} {} {} {} {} {} {}",
                nhood.id,
                nhood.central_tile_id,
                nhood.gap_start,
                nhood.gap_len,
                nhood.central_anchor_edge,
                n,
                angles_str,
                edges_str,
                inner_chains_str,
                vts_str,
            )?;
        }
        Ok(())
    }

    pub fn parse_file(tileset: Arc<TileSet<T>>, input: &str) -> Result<Self, String> {
        let mut entries = Vec::new();
        for line in input.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.is_empty() || parts[0] != "NTYPE" {
                continue;
            }
            if parts.len() < 7 {
                return Err(format!("NTYPE line too short: {}", line));
            }

            let id: usize = parts[1]
                .parse()
                .map_err(|_| format!("bad id: {}", parts[1]))?;
            let central_tile_id: usize = parts[2]
                .parse()
                .map_err(|_| format!("bad central_tile_id: {}", parts[2]))?;
            let gap_start: usize = parts[3]
                .parse()
                .map_err(|_| format!("bad gap_start: {}", parts[3]))?;
            let gap_len: u8 = parts[4]
                .parse()
                .map_err(|_| format!("bad gap_len: {}", parts[4]))?;
            let central_anchor_edge: usize = parts[5]
                .parse()
                .map_err(|_| format!("bad central_anchor_edge: {}", parts[5]))?;
            let n: usize = parts[6]
                .parse()
                .map_err(|_| format!("bad n: {}", parts[6]))?;

            if parts.len() < 7 + n {
                return Err(format!(
                    "NTYPE {} expected {} angles, got {} fields",
                    id,
                    n,
                    parts.len() - 7
                ));
            }
            let angles: Vec<i8> = parts[7..7 + n]
                .iter()
                .map(|s| s.parse::<i8>())
                .collect::<Result<Vec<_>, _>>()
                .map_err(|e| format!("bad angle: {}", e))?;

            let edges_start = 7 + n;
            if parts.len() < edges_start + n {
                return Err(format!(
                    "NTYPE {} expected {} edges, got {} fields",
                    id,
                    n,
                    parts.len() - edges_start
                ));
            }
            let edges: Vec<EdgeInfo> = parts[edges_start..edges_start + n]
                .iter()
                .map(|s| parse_edge_info(s))
                .collect::<Result<Vec<_>, _>>()?;

            let inner_start = edges_start + n;
            let inner_chains: Vec<Vec<EdgeInfo>> = if parts.len() >= inner_start + n {
                parts[inner_start..inner_start + n]
                    .iter()
                    .map(|s| {
                        if *s == "-" {
                            Ok(vec![])
                        } else {
                            s.split(',')
                                .map(parse_edge_info)
                                .collect::<Result<Vec<_>, _>>()
                        }
                    })
                    .collect::<Result<Vec<_>, _>>()?
            } else {
                vec![vec![]; n]
            };

            let vts_start = inner_start + n;
            let context_vertex_types = if parts.len() > vts_start {
                if parts[vts_start..] == ["-"] {
                    Vec::new()
                } else {
                    parts[vts_start..]
                        .iter()
                        .map(|s| parse_vtype(s))
                        .collect::<Result<Vec<_>, _>>()?
                }
            } else {
                extract_covered_vertex_types(&edges, &inner_chains, gap_start, gap_len as usize)
            };

            entries.push(NeighborhoodType {
                id,
                central_tile_id,
                central_anchor_edge,
                gap_len,
                context_vertex_types,
                angles,
                edges,
                inner_chains,
                gap_start,
            });
        }
        Ok(NeighborhoodIndex {
            tileset,
            entries,
            transitions: Vec::new(),
        })
    }
}

fn format_vtype(vt: &VertexType) -> String {
    let inner_str = if vt.inner.is_empty() {
        "-".to_string()
    } else {
        vt.inner
            .iter()
            .map(|e| format!("{}.{}", e.tile_id, e.tile_offset))
            .collect::<Vec<_>>()
            .join(",")
    };
    format!(
        "{}.{};{};{}.{}",
        vt.cw.tile_id, vt.cw.tile_offset, inner_str, vt.ccw.tile_id, vt.ccw.tile_offset
    )
}

fn parse_vtype(s: &str) -> Result<VertexType, String> {
    let parts: Vec<&str> = s.split(';').collect();
    if parts.len() != 3 {
        return Err(format!("bad vertex type '{}', expected cw;inner;ccw", s));
    }
    let cw = parse_edge_info(parts[0])?;
    let inner = if parts[1] == "-" {
        vec![]
    } else {
        parts[1]
            .split(',')
            .map(parse_edge_info)
            .collect::<Result<Vec<_>, _>>()?
    };
    let ccw = parse_edge_info(parts[2])?;
    Ok(VertexType { cw, inner, ccw })
}

fn parse_edge_info(s: &str) -> Result<EdgeInfo, String> {
    let parts: Vec<&str> = s.split('.').collect();
    if parts.len() != 2 {
        return Err(format!("bad edge info '{}', expected tile_id.offset", s));
    }
    let tile_id: usize = parts[0]
        .parse()
        .map_err(|_| format!("bad tile_id: {}", parts[0]))?;
    let tile_offset: usize = parts[1]
        .parse()
        .map_err(|_| format!("bad tile_offset: {}", parts[1]))?;
    Ok(EdgeInfo {
        tile_id,
        tile_offset,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic::{ZZ12, ZZ4};
    use crate::intgeom::rat::Rat;
    use crate::intgeom::tiles;

    fn square_index_zz4() -> NeighborhoodIndex<ZZ4> {
        let sq: crate::intgeom::snake::Snake<ZZ4> = tiles::square();
        let rat = Rat::try_from(&sq).unwrap();
        NeighborhoodIndex::new(Arc::new(TileSet::new(vec![rat])))
    }

    #[test]
    fn square_has_valid_open_types() {
        let idx = square_index_zz4();
        assert!(idx.num_types() > 0, "expected non-empty NT collection");
        let invalid = idx.validate();
        assert!(invalid.is_empty(), "invalid NT ids: {invalid:?}");
    }

    #[test]
    fn square_roundtrip_collection() {
        let idx = square_index_zz4();
        let mut buf = Vec::new();
        idx.write_collection(&mut buf).unwrap();
        let serialized = String::from_utf8(buf).unwrap();

        let sq: crate::intgeom::snake::Snake<ZZ4> = tiles::square();
        let rat = Rat::try_from(&sq).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));

        let parsed = NeighborhoodIndex::parse_file(ts, &serialized).unwrap();
        assert_eq!(parsed.num_types(), idx.num_types());
    }

    #[test]
    fn classify_vector_matches_entry_count() {
        let idx = square_index_zz4();
        let kinds = idx.classify_all();
        assert_eq!(kinds.len(), idx.num_types());
    }

    #[test]
    fn spectre_collection_validates() {
        let sp: crate::intgeom::snake::Snake<ZZ12> = tiles::spectre();
        let rat = Rat::try_from(&sp).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = NeighborhoodIndex::new(ts);
        assert!(idx.num_types() > 0, "expected spectre NTs");
        let invalid = idx.validate();
        assert!(invalid.is_empty(), "invalid spectre NT ids: {invalid:?}");
    }
}
