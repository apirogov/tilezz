use rustc_hash::FxHashSet;
use std::sync::Arc;

use crate::cyclotomic::geometry::intersect;
use crate::cyclotomic::{IsComplex, IsRingOrField, Units};
use crate::intgeom::angles;
use crate::intgeom::grid::UnitSquareGrid;
use crate::intgeom::matchtypes::{
    is_single_edge_candidate, junction_gap_nonnegative, CandidateMatch, MatchTypeIndex,
};
use crate::intgeom::rat::Rat;
use crate::intgeom::snake::Snake;
use crate::intgeom::tileset::TileSet;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct EdgeInfo {
    pub tile_id: usize,
    pub tile_offset: usize,
}

/// Auxiliary per-vertex info on a patch boundary.
///
/// Each boundary position has an angle and two adjacent edges (cw, ccw).
/// This is instrumental infrastructure — the interesting structure lives
/// in [`VertexType`] at actual junction vertices.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct PatchVertexInfo {
    pub angle: i8,
    pub cw: EdgeInfo,
    pub ccw: EdgeInfo,
}

pub struct TileSegment {
    pub patch_start: usize,
    pub patch_end: usize,
    pub tile_id: usize,
    pub offset_start: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PatchMatch {
    pub start_a: usize,
    pub len: usize,
    pub start_b: usize,
    pub tile_id: usize,
}

/// Junction vertex type: the arrangement of tiles meeting at a single junction vertex.
///
/// At a junction (where the boundary angle differs from the tile's internal angle),
/// multiple tile instances converge. `cw` is the edge arriving from the CW direction,
/// `ccw` is the edge departing in the CCW direction, and `inner` lists any enclosed
/// tile edges between them. Only meaningful at junction positions.
#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct VertexType {
    pub cw: EdgeInfo,
    pub inner: Vec<EdgeInfo>,
    pub ccw: EdgeInfo,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum TransitionSide {
    Cw,
    Ccw,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Transition {
    pub src_id: usize,
    pub dst_id: Option<usize>,
    pub side: TransitionSide,
    pub tile_id: usize,
    pub tile_offset: usize,
}

pub struct AddTileDiff;

#[derive(Clone)]
pub struct GrowingPatch<T: IsComplex> {
    match_index: Arc<MatchTypeIndex<T>>,
    state: PatchState<T>,
}

#[derive(Clone)]
enum PatchState<T: IsComplex> {
    Seed {
        tile_id: usize,
        cached_matches: Vec<PatchMatch>,
    },
    Growing {
        angles: Vec<i8>,
        edges: Vec<EdgeInfo>,
        candidates_by_start: Option<Vec<Vec<PatchMatch>>>,
        inner_chains: Vec<Vec<EdgeInfo>>,
        positions: Vec<T>,
        grid: UnitSquareGrid,
        seg_data: Vec<(T, T)>,
        boundary_edge_ids: Vec<usize>,
        next_edge_id: usize,
    },
    _Phantom(std::marker::PhantomData<T>),
}

/// Raw vertex type extraction at a boundary position (no junction check).
///
/// Constructs the (cw, inner, ccw) triple from boundary edge/inner-chain data.
/// Used internally; prefer [`GrowingPatch::junction_vertex_type_at`] which
/// returns `None` at non-junction positions.
pub(crate) fn vertex_type_raw_from(
    edges: &[EdgeInfo],
    inner_chains: &[Vec<EdgeInfo>],
    pos: usize,
) -> VertexType {
    let n = edges.len();
    VertexType {
        cw: edges[(pos + n - 1) % n],
        inner: inner_chains[pos].clone(),
        ccw: edges[pos],
    }
}

/// Extract the junction vertex type at a boundary position.
///
/// Returns `None` if `pos` is not a junction vertex (i.e. the boundary angle
/// equals the tile's internal angle at that edge).
pub fn junction_vertex_type_from<T: IsComplex + IsRingOrField + Units>(
    edges: &[EdgeInfo],
    inner_chains: &[Vec<EdgeInfo>],
    angles: &[i8],
    tileset: &TileSet<T>,
    pos: usize,
) -> Option<VertexType> {
    let ei = edges[pos];
    let tile = tileset.rat(ei.tile_id);
    if tile.seq()[ei.tile_offset] == angles[pos] {
        return None;
    }
    Some(vertex_type_raw_from(edges, inner_chains, pos))
}

fn forward_match_length(
    self_angles: &[i8],
    self_start: usize,
    other_angles: &[i8],
    other_junction: usize,
) -> usize {
    let n = self_angles.len();
    let m = other_angles.len();
    let max_len = n.min(m);
    let mut len = 1;
    for i in 1..max_len {
        let xi = self_angles[(self_start + i) % n];
        let yi = -other_angles[(other_junction + m - i) % m];
        if xi != yi {
            break;
        }
        len += 1;
    }
    len
}

pub fn update_inner_chains(
    old_inner: &[Vec<EdgeInfo>],
    old_edges: &[EdgeInfo],
    pm: &PatchMatch,
    new_n: usize,
) -> Vec<Vec<EdgeInfo>> {
    let n = old_edges.len();
    let seg_len_old = n - pm.len;
    let ccw_pos = (pm.start_a + pm.len) % n;
    let cw_end_matched = (pm.start_a + pm.len - 1) % n;

    let mut new_inner = vec![Vec::new(); new_n];

    for i in 1..seg_len_old {
        new_inner[i] = old_inner[(ccw_pos + i) % n].clone();
    }

    new_inner[0] = old_inner[ccw_pos]
        .iter()
        .chain(std::iter::once(&old_edges[cw_end_matched]))
        .cloned()
        .collect();

    new_inner[seg_len_old] = old_inner[pm.start_a]
        .iter()
        .chain(std::iter::once(&old_edges[pm.start_a]))
        .cloned()
        .collect();

    new_inner
}

impl<T: IsComplex + IsRingOrField + Units> GrowingPatch<T> {
    pub fn new(tileset: Arc<TileSet<T>>, seed_tile_id: usize) -> Self {
        let match_index = Arc::new(MatchTypeIndex::new(Arc::clone(&tileset)));
        let seed_matches =
            Self::compute_seed_matches(&match_index, match_index.tileset(), seed_tile_id);
        GrowingPatch {
            match_index,
            state: PatchState::Seed {
                tile_id: seed_tile_id,
                cached_matches: seed_matches,
            },
        }
    }

    pub fn from_parts(
        match_index: Arc<MatchTypeIndex<T>>,
        angles: Vec<i8>,
        edges: Vec<EdgeInfo>,
        inner_chains: Vec<Vec<EdgeInfo>>,
    ) -> Option<Self> {
        if angles.is_empty() || angles.len() != edges.len() || inner_chains.len() != angles.len() {
            return None;
        }
        let n = angles.len();
        let positions = trace_boundary_positions::<T>(&angles);
        let (grid, seg_data, boundary_edge_ids, next_edge_id) =
            build_boundary_grid::<T>(&positions, n);
        Some(GrowingPatch {
            match_index,
            state: PatchState::Growing {
                angles,
                edges,
                candidates_by_start: None,
                inner_chains,
                positions,
                grid,
                seg_data,
                boundary_edge_ids,
                next_edge_id,
            },
        })
    }

    pub fn construct_minimal_witness(
        vtype: &VertexType,
        match_index: Arc<MatchTypeIndex<T>>,
    ) -> Option<(Self, usize)> {
        let tileset = match_index.tileset();
        let seed_id = vtype.cw.tile_id;
        let seed_seq = tileset.rat(seed_id).seq();
        let n_seed = seed_seq.len();

        let mut boundary_angles = seed_seq.to_vec();
        let mut boundary_edges: Vec<EdgeInfo> = (0..n_seed)
            .map(|i| EdgeInfo {
                tile_id: seed_id,
                tile_offset: i,
            })
            .collect();
        let mut boundary_inner: Vec<Vec<EdgeInfo>> = vec![vec![]; n_seed];

        let mut junc_pos = (vtype.cw.tile_offset + 1) % n_seed;

        let mut targets = vtype.inner.clone();
        targets.push(vtype.ccw);

        let mut first_step = true;

        for target in &targets {
            let tile_seq = tileset.rat(target.tile_id).seq();
            let m = tile_seq.len();
            let tile_junc = target.tile_offset;
            let n = boundary_angles.len();

            let mlen = forward_match_length(&boundary_angles, junc_pos, tile_seq, tile_junc);
            if mlen == 0 || mlen > n || mlen > m {
                return None;
            }

            let ccw_pos = (junc_pos + mlen) % n;
            let seg_len_old = n - mlen;
            let seg_len_new = m - mlen;
            let new_n = seg_len_old + seg_len_new;

            let mut new_edges = Vec::with_capacity(new_n);
            for i in 0..seg_len_old {
                new_edges.push(boundary_edges[(ccw_pos + i) % n]);
            }
            new_edges.push(EdgeInfo {
                tile_id: target.tile_id,
                tile_offset: tile_junc,
            });
            for k in 1..seg_len_new {
                new_edges.push(EdgeInfo {
                    tile_id: target.tile_id,
                    tile_offset: (tile_junc + mlen + k) % m,
                });
            }

            let new_inner = if first_step {
                first_step = false;
                vec![vec![]; new_n]
            } else {
                let pm = PatchMatch {
                    start_a: junc_pos,
                    len: mlen,
                    start_b: tile_junc,
                    tile_id: target.tile_id,
                };
                update_inner_chains(&boundary_inner, &boundary_edges, &pm, new_n)
            };

            let gr = angles::glue_raw_angles::<T>(
                &boundary_angles,
                tile_seq,
                junc_pos,
                mlen,
                tile_junc,
            )?;
            boundary_angles = gr.angles;
            boundary_edges = new_edges;
            boundary_inner = new_inner;

            junc_pos = seg_len_old;
        }

        let patch = GrowingPatch::from_parts(
            Arc::clone(&match_index),
            boundary_angles,
            boundary_edges,
            boundary_inner,
        )?;

        let final_edges = patch.edges();
        let final_inner = patch.inner_chains();
        for pos in 0..patch.boundary_len() {
            let vt = vertex_type_raw_from(final_edges, final_inner, pos);
            if vt == *vtype {
                return Some((patch, pos));
            }
        }

        None
    }

    pub fn compute_candidates_covering_position(
        match_index: &Arc<MatchTypeIndex<T>>,
        angles: &[i8],
        edges: &[EdgeInfo],
        target: usize,
    ) -> Vec<PatchMatch> {
        let n = angles.len();
        let rat = Rat::from_slice_unchecked(angles);
        let tileset = match_index.tileset();
        let max_tile_len = tileset.rats().iter().map(|r| r.len()).max().unwrap_or(0);
        let mut global_seen: FxHashSet<(usize, usize, usize, usize)> = FxHashSet::default();
        let mut result: Vec<PatchMatch> = Vec::new();

        let segments = compute_segments(angles, edges, tileset);
        let junctions_set: FxHashSet<usize> = compute_junctions(angles, edges, tileset)
            .into_iter()
            .collect();

        for offset in 0..max_tile_len.min(n) {
            let pos = (target + n - offset) % n;

            if let Some(segment) = segments
                .iter()
                .find(|s| pos >= s.patch_start && pos < s.patch_end)
            {
                let tile_id = segment.tile_id;
                let tile_offset = (segment.offset_start + (pos - segment.patch_start))
                    % tileset.rat(tile_id).len();
                let mut tmp_by_start: Vec<Vec<PatchMatch>> = vec![Vec::new(); n];
                compute_candidates_at_position(
                    pos,
                    angles,
                    &rat,
                    n,
                    tile_id,
                    tile_offset,
                    match_index,
                    tileset,
                    &mut global_seen,
                    &mut tmp_by_start,
                );
                for pm in tmp_by_start.into_iter().flatten() {
                    if cyclic_range_contains(pm.start_a, pm.len, target, n) {
                        result.push(pm);
                    }
                }
            }

            if junctions_set.contains(&pos) {
                for tile_id_b in 0..tileset.num_tiles() {
                    let tile_b = tileset.rat(tile_id_b);
                    let b_seq = tile_b.seq();
                    let m = b_seq.len();
                    for ib in 0..m {
                        if !is_single_edge_candidate(angles, pos, b_seq, ib) {
                            continue;
                        }
                        let (ns, len, ne) = rat.get_match((pos as i64, ib as i64), tile_b);
                        if len != 1 {
                            continue;
                        }
                        let ns_u = ns.rem_euclid(n as i64) as usize;
                        let ne_u = ne.rem_euclid(m as i64) as usize;
                        let key = (ns_u, len, ne_u, tile_id_b);
                        if !global_seen.insert(key) {
                            continue;
                        }
                        if cyclic_range_contains(ns_u, len, target, n)
                            && rat
                                .try_glue_precomputed((ns, len, ne), tile_b, true)
                                .is_ok()
                        {
                            result.push(PatchMatch {
                                start_a: ns_u,
                                len,
                                start_b: ne_u,
                                tile_id: tile_id_b,
                            });
                        }
                    }
                }
            }
        }

        result
    }

    pub fn compute_all_candidates(
        match_index: &Arc<MatchTypeIndex<T>>,
        angles: &[i8],
        edges: &[EdgeInfo],
    ) -> Vec<Vec<PatchMatch>> {
        let n = angles.len();
        let rat = Rat::from_slice_unchecked(angles);
        let tileset = match_index.tileset();
        let mut result = vec![Vec::new(); n];

        let segments = compute_segments(angles, edges, tileset);
        let junctions = compute_junctions(angles, edges, tileset);

        let mut global_seen: FxHashSet<(usize, usize, usize, usize)> = FxHashSet::default();

        for segment in &segments {
            let tile_id = segment.tile_id;
            let seg_len = segment.patch_end - segment.patch_start;
            for local_k in 0..seg_len {
                let tile_offset = (segment.offset_start + local_k) % tileset.rat(tile_id).len();
                let patch_pos = segment.patch_start + local_k;
                compute_candidates_at_position(
                    patch_pos,
                    angles,
                    &rat,
                    n,
                    tile_id,
                    tile_offset,
                    match_index,
                    tileset,
                    &mut global_seen,
                    &mut result,
                );
            }
        }

        for &junc_idx in &junctions {
            for tile_id_b in 0..tileset.num_tiles() {
                let tile_b = tileset.rat(tile_id_b);
                let b_seq = tile_b.seq();
                let m = b_seq.len();
                for ib in 0..m {
                    if !is_single_edge_candidate(angles, junc_idx, b_seq, ib) {
                        continue;
                    }
                    let (ns, len, ne) = rat.get_match((junc_idx as i64, ib as i64), tile_b);
                    if len != 1 {
                        continue;
                    }
                    let ns_u = ns.rem_euclid(n as i64) as usize;
                    let ne_u = ne.rem_euclid(m as i64) as usize;
                    let key = (ns_u, len, ne_u, tile_id_b);
                    if !global_seen.insert(key) {
                        continue;
                    }
                    if rat
                        .try_glue_precomputed((ns, len, ne), tile_b, true)
                        .is_ok()
                    {
                        result[ns_u].push(PatchMatch {
                            start_a: ns_u,
                            len,
                            start_b: ne_u,
                            tile_id: tile_id_b,
                        });
                    }
                }
            }
        }

        result
    }

    pub fn is_growing(&self) -> bool {
        matches!(self.state, PatchState::Growing { .. })
    }

    pub(crate) fn clone_for_mutation(&self) -> Self {
        GrowingPatch {
            match_index: Arc::clone(&self.match_index),
            state: match &self.state {
                PatchState::Seed { tile_id, .. } => PatchState::Seed {
                    tile_id: *tile_id,
                    cached_matches: Vec::new(),
                },
                PatchState::Growing {
                    angles,
                    edges,
                    candidates_by_start: _,
                    inner_chains,
                    positions,
                    grid,
                    seg_data,
                    boundary_edge_ids,
                    next_edge_id,
                } => PatchState::Growing {
                    angles: angles.clone(),
                    edges: edges.clone(),
                    candidates_by_start: None,
                    inner_chains: inner_chains.clone(),
                    positions: positions.clone(),
                    grid: grid.clone(),
                    seg_data: seg_data.clone(),
                    boundary_edge_ids: boundary_edge_ids.clone(),
                    next_edge_id: *next_edge_id,
                },
                PatchState::_Phantom(p) => PatchState::_Phantom(*p),
            },
        }
    }

    pub fn get_all_matches(&self) -> Vec<PatchMatch> {
        match &self.state {
            PatchState::Seed { cached_matches, .. } => cached_matches.clone(),
            PatchState::Growing {
                angles,
                edges,
                candidates_by_start,
                ..
            } => match candidates_by_start {
                Some(cbs) => cbs.iter().flatten().cloned().collect(),
                None => Self::compute_all_candidates(&self.match_index, angles, edges)
                    .into_iter()
                    .flatten()
                    .collect(),
            },
            PatchState::_Phantom(_) => Vec::new(),
        }
    }

    pub fn get_matches_for_tile(&self, tile_id: usize) -> impl Iterator<Item = PatchMatch> + '_ {
        self.get_all_matches()
            .into_iter()
            .filter(move |m| m.tile_id == tile_id)
    }

    pub fn get_matches_touching_vertex(&self, vertex_index: usize) -> Vec<PatchMatch> {
        match &self.state {
            PatchState::Seed { cached_matches, .. } => {
                let n = self.match_index.tileset().rat(0).len();
                cached_matches
                    .iter()
                    .filter(|pm| cyclic_range_contains(pm.start_a, pm.len, vertex_index, n))
                    .cloned()
                    .collect()
            }
            PatchState::Growing {
                candidates_by_start,
                angles,
                edges,
                ..
            } => {
                let n = angles.len();
                match candidates_by_start {
                    Some(cbs) => {
                        let k = self
                            .match_index
                            .tileset()
                            .rats()
                            .iter()
                            .map(|r| r.len())
                            .max()
                            .unwrap_or(0)
                            .min(n);
                        let mut result = Vec::new();
                        for offset in 0..=k {
                            let start = (vertex_index + n - offset) % n;
                            for pm in &cbs[start] {
                                if cyclic_range_contains(pm.start_a, pm.len, vertex_index, n) {
                                    result.push(pm.clone());
                                }
                            }
                        }
                        result
                    }
                    None => {
                        let all = Self::compute_all_candidates(&self.match_index, angles, edges);
                        let k = self
                            .match_index
                            .tileset()
                            .rats()
                            .iter()
                            .map(|r| r.len())
                            .max()
                            .unwrap_or(0)
                            .min(n);
                        let mut result = Vec::new();
                        for offset in 0..=k {
                            let start = (vertex_index + n - offset) % n;
                            for pm in &all[start] {
                                if cyclic_range_contains(pm.start_a, pm.len, vertex_index, n) {
                                    result.push(pm.clone());
                                }
                            }
                        }
                        result
                    }
                }
            }
            PatchState::_Phantom(_) => Vec::new(),
        }
    }

    pub fn ensure_candidates_materialized(&mut self) {
        if let PatchState::Growing {
            angles,
            edges,
            candidates_by_start,
            ..
        } = &mut self.state
        {
            if candidates_by_start.is_none() {
                *candidates_by_start = Some(Self::compute_all_candidates(
                    &self.match_index,
                    angles,
                    edges,
                ));
            }
        }
    }

    pub fn num_tiles(&self) -> usize {
        self.match_index.tileset().num_tiles()
    }

    pub fn boundary_len(&self) -> usize {
        match &self.state {
            PatchState::Growing { angles, .. } => angles.len(),
            _ => 0,
        }
    }

    pub fn angles(&self) -> &[i8] {
        match &self.state {
            PatchState::Growing { angles, .. } => angles,
            _ => &[],
        }
    }

    pub fn to_rat(&self) -> Rat<T> {
        match &self.state {
            PatchState::Seed { tile_id, .. } => self.match_index.tileset().rat(*tile_id).clone(),
            PatchState::Growing { angles, .. } => Rat::from_slice_unchecked(angles),
            _ => panic!("patch has no boundary"),
        }
    }

    pub fn tileset(&self) -> &Arc<TileSet<T>> {
        self.match_index.tileset()
    }

    pub fn match_index(&self) -> &Arc<MatchTypeIndex<T>> {
        &self.match_index
    }

    pub fn edges(&self) -> &[EdgeInfo] {
        match &self.state {
            PatchState::Growing { edges, .. } => edges,
            _ => &[],
        }
    }

    pub fn inner_chains(&self) -> &[Vec<EdgeInfo>] {
        match &self.state {
            PatchState::Growing { inner_chains, .. } => inner_chains,
            _ => &[],
        }
    }

    pub fn candidates_by_start(&self) -> &[Vec<PatchMatch>] {
        match &self.state {
            PatchState::Growing {
                candidates_by_start,
                ..
            } => candidates_by_start.as_deref().unwrap_or(&[]),
            _ => &[],
        }
    }

    pub fn vertex_type_at(&self, i: usize) -> Option<PatchVertexInfo> {
        let n = self.boundary_len();
        if n == 0 || i >= n {
            return None;
        }
        let edges = match &self.state {
            PatchState::Growing { edges, .. } => edges,
            _ => return None,
        };
        let angles = self.angles();
        Some(PatchVertexInfo {
            angle: angles[i],
            cw: edges[(i + n - 1) % n],
            ccw: edges[i],
        })
    }

    /// Return the junction vertex type at position `i`, or `None` if not a junction.
    pub fn junction_vertex_type_at(&self, i: usize) -> Option<VertexType> {
        let n = self.boundary_len();
        if n == 0 || i >= n {
            return None;
        }
        match &self.state {
            PatchState::Growing {
                edges,
                inner_chains,
                ..
            } => {
                if !self.is_junction(i) {
                    return None;
                }
                let n = edges.len();
                Some(VertexType {
                    cw: edges[(i + n - 1) % n],
                    inner: inner_chains[i].clone(),
                    ccw: edges[i],
                })
            }
            _ => None,
        }
    }

    pub fn is_junction(&self, i: usize) -> bool {
        let edges = match &self.state {
            PatchState::Growing { edges, .. } => edges,
            _ => return false,
        };
        let n = edges.len();
        if n == 0 {
            return false;
        }
        let ei = edges[i];
        let tileset = self.match_index.tileset();
        let tile = tileset.rat(ei.tile_id);
        tile.seq()[ei.tile_offset] != self.angles()[i]
    }

    pub fn neighbor_junction_offsets(&self, pos: usize) -> Option<(usize, usize)> {
        let (edges, n) = match &self.state {
            PatchState::Growing { edges, .. } => (edges, edges.len()),
            _ => return None,
        };
        if n == 0 || pos >= n {
            return None;
        }

        let mut j_cw = (pos + n - 1) % n;
        while j_cw != pos && !self.is_junction(j_cw) {
            j_cw = (j_cw + n - 1) % n;
        }

        let mut j_ccw = (pos + 1) % n;
        while j_ccw != pos && !self.is_junction(j_ccw) {
            j_ccw = (j_ccw + 1) % n;
        }

        let cw_offset = edges[j_cw].tile_offset;
        let ccw_prev = (j_ccw + n - 1) % n;
        let ccw_edge = edges[ccw_prev];
        let tile_len = self.match_index.tileset().rat(ccw_edge.tile_id).len();
        let ccw_offset = (ccw_edge.tile_offset + 1) % tile_len;

        Some((cw_offset, ccw_offset))
    }

    pub fn junction_vertices(&self) -> Vec<(usize, PatchVertexInfo)> {
        let n = self.boundary_len();
        let mut result = Vec::new();
        for i in 0..n {
            if self.is_junction(i) {
                result.push((i, self.vertex_type_at(i).unwrap()));
            }
        }
        result
    }

    pub fn tile_segments(&self) -> Vec<TileSegment> {
        let edges = match &self.state {
            PatchState::Growing { edges, .. } => edges,
            _ => return vec![],
        };
        let n = edges.len();
        if n == 0 {
            return vec![];
        }
        compute_segments(self.angles(), edges, self.match_index.tileset())
    }

    pub fn add_tile(&mut self, pm: &PatchMatch) -> Option<AddTileDiff> {
        match &self.state {
            PatchState::Seed { tile_id, .. } => {
                let seed_id = *tile_id;
                self.init_from_first_add(seed_id, pm)
            }
            PatchState::Growing { .. } => self.add_tile_growing(pm),
            PatchState::_Phantom(_) => None,
        }
    }

    fn init_from_first_add(&mut self, seed_id: usize, pm: &PatchMatch) -> Option<AddTileDiff> {
        let tileset = self.match_index.tileset();
        let seed_rat = tileset.rat(seed_id);
        let n = seed_rat.seq().len();
        let m = tileset.rat(pm.tile_id).seq().len();
        let mlen = pm.len;

        if mlen == 0 || mlen > n || mlen > m {
            return None;
        }

        let (_ns, seed_len, _ne) = seed_rat.get_match(
            (pm.start_a as i64, pm.start_b as i64),
            tileset.rat(pm.tile_id),
        );
        if seed_len == 0 {
            return None;
        }

        let seed_angles = seed_rat.seq().to_vec();
        let new_angles = compute_glue_angles::<T>(&seed_angles, pm, tileset)?;

        let seg_len_old = n - mlen;
        let seg_len_new = m - mlen;
        let new_len = seg_len_old + seg_len_new;
        if new_len == 0 {
            return None;
        }

        if Snake::<T>::try_from(new_angles.as_slice()).is_err() {
            return None;
        }

        let positions = trace_boundary_positions::<T>(&new_angles);
        let (grid, seg_data, boundary_edge_ids, next_edge_id) =
            build_boundary_grid::<T>(&positions, new_len);

        let ccw_pos = (pm.start_a + mlen) % n;

        let mut edges = Vec::with_capacity(new_len);
        for i in 0..seg_len_old {
            edges.push(EdgeInfo {
                tile_id: seed_id,
                tile_offset: (ccw_pos + i) % n,
            });
        }
        edges.push(EdgeInfo {
            tile_id: pm.tile_id,
            tile_offset: pm.start_b,
        });
        for k in 1..seg_len_new {
            edges.push(EdgeInfo {
                tile_id: pm.tile_id,
                tile_offset: (pm.start_b + mlen + k) % m,
            });
        }

        debug_assert_eq!(edges.len(), new_len);

        let inner_chains = vec![vec![]; new_len];

        self.state = PatchState::Growing {
            angles: new_angles,
            edges,
            candidates_by_start: None,
            inner_chains,
            positions,
            grid,
            seg_data,
            boundary_edge_ids,
            next_edge_id,
        };

        Some(AddTileDiff)
    }

    fn add_tile_growing(&mut self, pm: &PatchMatch) -> Option<AddTileDiff> {
        let (
            angles,
            edges,
            old_inner,
            positions,
            mut grid,
            seg_data,
            boundary_edge_ids,
            next_edge_id,
        ) = match &mut self.state {
            PatchState::Growing {
                angles,
                edges,
                candidates_by_start: _,
                inner_chains,
                positions,
                grid,
                seg_data,
                boundary_edge_ids,
                next_edge_id,
            } => (
                std::mem::take(angles),
                std::mem::take(edges),
                std::mem::take(inner_chains),
                std::mem::take(positions),
                std::mem::take(grid),
                std::mem::take(seg_data),
                std::mem::take(boundary_edge_ids),
                *next_edge_id,
            ),
            _ => return None,
        };

        let n = angles.len();
        let mlen = pm.len;
        let tileset = self.match_index.tileset();
        let m = tileset.rat(pm.tile_id).seq().len();

        if mlen == 0 || mlen > n || mlen > m {
            self.restore_growing(
                angles,
                edges,
                old_inner,
                positions,
                grid,
                seg_data,
                boundary_edge_ids,
                next_edge_id,
            );
            return None;
        }

        let new_angles = match compute_glue_angles::<T>(&angles, pm, tileset) {
            Some(a) => a,
            None => {
                self.restore_growing(
                    angles,
                    edges,
                    old_inner,
                    positions,
                    grid,
                    seg_data,
                    boundary_edge_ids,
                    next_edge_id,
                );
                return None;
            }
        };

        let seg_len_old = n - mlen;
        let seg_len_new = m - mlen;
        let new_len = seg_len_old + seg_len_new;

        if new_len == 0 {
            return None;
        }

        let ccw_pos = (pm.start_a + mlen) % n;

        let removed_ids: Vec<usize> = (0..mlen)
            .map(|i| boundary_edge_ids[(pm.start_a + i) % n])
            .collect();

        for &id in &removed_ids {
            unregister_segment::<T>(&mut grid, &seg_data, id);
        }

        let abs_dir_at_cw: i8 = {
            let prev = (pm.start_a + n - 1) % n;
            dir_of_edge::<T>(positions[prev], positions[pm.start_a])
        };

        let cw_junction = positions[pm.start_a];
        let ccw_junction = positions[ccw_pos];
        let mut new_tile_positions = vec![cw_junction];
        let mut dir = abs_dir_at_cw;
        for &a in &new_angles[seg_len_old..] {
            dir = (dir as i64 + a as i64).rem_euclid(T::turn() as i64) as i8;
            let last = *new_tile_positions.last().unwrap();
            new_tile_positions.push(last + T::unit(dir));
        }

        debug_assert_eq!(
            *new_tile_positions.last().unwrap(),
            ccw_junction,
            "new tile trace should end at CCW junction"
        );

        let mut all_clear = true;
        for i in 0..seg_len_new {
            let p1 = new_tile_positions[i];
            let p2 = new_tile_positions[i + 1];
            let allowed = if i == seg_len_new - 1 {
                Some(ccw_junction)
            } else {
                None
            };
            if !check_segment_clear::<T>(&grid, &seg_data, p1, p2, allowed) {
                all_clear = false;
                break;
            }
        }

        if !all_clear {
            for &id in &removed_ids {
                let (p1, p2) = seg_data[id];
                register_segment::<T>(&mut grid, p1, p2, id);
            }
            self.restore_growing(
                angles,
                edges,
                old_inner,
                positions,
                grid,
                seg_data,
                boundary_edge_ids,
                next_edge_id,
            );
            return None;
        }

        let mut new_positions = Vec::with_capacity(new_len + 1);
        for i in 0..=seg_len_old {
            new_positions.push(positions[(ccw_pos + i) % n]);
        }
        for p in new_tile_positions.iter().take(seg_len_new + 1).skip(1) {
            new_positions.push(*p);
        }

        let mut new_boundary_edge_ids = Vec::with_capacity(new_len);
        let mut new_seg_data = seg_data;
        let mut next_id = next_edge_id;

        for i in 0..seg_len_old {
            new_boundary_edge_ids.push(boundary_edge_ids[(ccw_pos + i) % n]);
        }

        for i in 0..seg_len_new {
            let id = next_id;
            next_id += 1;
            new_boundary_edge_ids.push(id);
            let p1 = new_tile_positions[i];
            let p2 = new_tile_positions[i + 1];
            new_seg_data.push((p1, p2));
            register_segment::<T>(&mut grid, p1, p2, id);
        }

        let mut new_edges = Vec::with_capacity(new_len);
        for i in 0..seg_len_old {
            new_edges.push(edges[(ccw_pos + i) % n]);
        }
        new_edges.push(EdgeInfo {
            tile_id: pm.tile_id,
            tile_offset: pm.start_b,
        });
        for k in 1..seg_len_new {
            new_edges.push(EdgeInfo {
                tile_id: pm.tile_id,
                tile_offset: (pm.start_b + mlen + k) % m,
            });
        }

        debug_assert_eq!(new_edges.len(), new_len);

        let new_inner = update_inner_chains(&old_inner, &edges, pm, new_len);

        self.state = PatchState::Growing {
            angles: new_angles,
            edges: new_edges,
            candidates_by_start: None,
            inner_chains: new_inner,
            positions: new_positions,
            grid,
            seg_data: new_seg_data,
            boundary_edge_ids: new_boundary_edge_ids,
            next_edge_id: next_id,
        };

        Some(AddTileDiff)
    }

    #[allow(clippy::too_many_arguments)]
    fn restore_growing(
        &mut self,
        angles: Vec<i8>,
        edges: Vec<EdgeInfo>,
        inner_chains: Vec<Vec<EdgeInfo>>,
        positions: Vec<T>,
        grid: UnitSquareGrid,
        seg_data: Vec<(T, T)>,
        boundary_edge_ids: Vec<usize>,
        next_edge_id: usize,
    ) {
        self.state = PatchState::Growing {
            angles,
            edges,
            candidates_by_start: None,
            inner_chains,
            positions,
            grid,
            seg_data,
            boundary_edge_ids,
            next_edge_id,
        };
    }

    fn compute_seed_matches(
        match_index: &MatchTypeIndex<T>,
        tileset: &TileSet<T>,
        seed_tile_id: usize,
    ) -> Vec<PatchMatch> {
        let seed = tileset.rat(seed_tile_id);
        let seed_seq = seed.seq();
        let n = seed_seq.len();
        let mut seen: FxHashSet<(usize, usize, usize, usize)> = FxHashSet::default();
        let mut matches = Vec::new();

        for offset in 0..n {
            for cand in match_index.candidates_starting_at(seed_tile_id, offset) {
                let tile_b = tileset.rat(cand.tile_b);
                let (ns, len, ne) = seed.get_match((offset as i64, cand.start_b as i64), tile_b);
                if len == 0 {
                    continue;
                }
                let ns_u = ns.rem_euclid(n as i64) as usize;
                let ne_u = ne.rem_euclid(tile_b.len() as i64) as usize;
                if !junction_gap_nonnegative(seed_seq, ns_u, len, tile_b.seq(), ne_u) {
                    continue;
                }
                if let Ok(_glued) = seed.try_glue_precomputed((ns, len, ne), tile_b, true) {
                    let key = (ns_u, len, ne_u, cand.tile_b);
                    if seen.insert(key) {
                        matches.push(PatchMatch {
                            start_a: ns_u,
                            len,
                            start_b: ne_u,
                            tile_id: cand.tile_b,
                        });
                    }
                }
            }
        }

        matches
    }
}

fn compute_segments<T: IsComplex + IsRingOrField + Units>(
    angles: &[i8],
    edges: &[EdgeInfo],
    tileset: &TileSet<T>,
) -> Vec<TileSegment> {
    let n = edges.len();
    if n == 0 {
        return vec![];
    }
    let mut segments: Vec<TileSegment> = Vec::new();
    let mut seg_start = 0;
    for i in 1..=n {
        let at_end = i == n;
        let tile_change = !at_end
            && (edges[i].tile_id != edges[seg_start].tile_id
                || edges[i].tile_offset != edges[seg_start].tile_offset + (i - seg_start));
        let junction_break = !at_end && {
            let ei = edges[i];
            tileset.rat(ei.tile_id).seq()[ei.tile_offset] != angles[i]
        };
        if at_end || tile_change || junction_break {
            segments.push(TileSegment {
                patch_start: seg_start,
                patch_end: i,
                tile_id: edges[seg_start].tile_id,
                offset_start: edges[seg_start].tile_offset,
            });
            seg_start = i;
        }
    }
    segments
}

fn compute_junctions<T: IsComplex + IsRingOrField + Units>(
    angles: &[i8],
    edges: &[EdgeInfo],
    tileset: &TileSet<T>,
) -> Vec<usize> {
    let n = edges.len();
    let mut juncs = Vec::new();
    for i in 0..n {
        let ei = edges[i];
        if tileset.rat(ei.tile_id).seq()[ei.tile_offset] != angles[i] {
            juncs.push(i);
        }
    }
    juncs
}

#[allow(clippy::too_many_arguments)]
fn compute_candidates_at_position<T: IsComplex + IsRingOrField + Units>(
    pos: usize,
    angles: &[i8],
    rat: &Rat<T>,
    n: usize,
    tile_id: usize,
    tile_offset: usize,
    match_index: &MatchTypeIndex<T>,
    tileset: &TileSet<T>,
    global_seen: &mut FxHashSet<(usize, usize, usize, usize)>,
    result: &mut [Vec<PatchMatch>],
) {
    for cand in match_index.candidates_starting_at(tile_id, tile_offset) {
        append_match_candidate(pos, angles, rat, n, cand, tileset, global_seen, result);
    }
}

#[allow(clippy::too_many_arguments)]
fn append_match_candidate<T: IsComplex + IsRingOrField + Units>(
    pos: usize,
    angles: &[i8],
    rat: &Rat<T>,
    n: usize,
    cand: &CandidateMatch,
    tileset: &TileSet<T>,
    seen: &mut FxHashSet<(usize, usize, usize, usize)>,
    result: &mut [Vec<PatchMatch>],
) {
    let tile_b = tileset.rat(cand.tile_b);
    let (ns, len, ne) = rat.get_match((pos as i64, cand.start_b as i64), tile_b);
    if len == 0 {
        return;
    }
    let ns_u = ns.rem_euclid(n as i64) as usize;
    let ne_u = ne.rem_euclid(tile_b.len() as i64) as usize;
    if !junction_gap_nonnegative(angles, ns_u, len, tile_b.seq(), ne_u) {
        return;
    }
    let key = (ns_u, len, ne_u, cand.tile_b);
    if !seen.insert(key) {
        return;
    }
    if rat
        .try_glue_precomputed((ns, len, ne), tile_b, true)
        .is_ok()
    {
        result[ns_u].push(PatchMatch {
            start_a: ns_u,
            len,
            start_b: ne_u,
            tile_id: cand.tile_b,
        });
    }
}

pub fn candidates_from_flat(n: usize, flat: Vec<PatchMatch>) -> Vec<Vec<PatchMatch>> {
    let mut result = vec![Vec::new(); n];
    for pm in flat {
        if pm.start_a < n {
            result[pm.start_a].push(pm);
        }
    }
    result
}

pub fn cyclic_range_contains(start: usize, len: usize, index: usize, n: usize) -> bool {
    if len == 0 || n == 0 {
        return false;
    }
    let end = start + len;
    if end <= n {
        index >= start && index <= end
    } else {
        index >= start || index <= end % n
    }
}

pub fn compute_glue_angles<T: IsComplex + IsRingOrField + Units>(
    angles: &[i8],
    pm: &PatchMatch,
    tileset: &Arc<TileSet<T>>,
) -> Option<Vec<i8>> {
    let other_seq = tileset.rat(pm.tile_id).seq();
    let gr = angles::glue_raw_angles::<T>(angles, other_seq, pm.start_a, pm.len, pm.start_b)?;
    if let (Some(a_yx), Some(a_xy)) = (gr.a_yx, gr.a_xy) {
        if a_yx.abs() == T::hturn() || a_xy.abs() == T::hturn() {
            return None;
        }
    }
    Some(gr.angles)
}

fn trace_boundary_positions<T: IsComplex + IsRingOrField + Units>(angles: &[i8]) -> Vec<T> {
    let n = angles.len();
    let mut positions = Vec::with_capacity(n + 1);
    positions.push(T::zero());
    let mut dir: i8 = 0;
    for &a in angles {
        dir = (dir as i64 + a as i64).rem_euclid(T::turn() as i64) as i8;
        let last = *positions.last().unwrap();
        positions.push(last + T::unit(dir));
    }
    positions
}

fn build_boundary_grid<T: IsComplex + IsRingOrField + Units>(
    positions: &[T],
    n: usize,
) -> (UnitSquareGrid, Vec<(T, T)>, Vec<usize>, usize) {
    let mut grid = UnitSquareGrid::new();
    let mut seg_data = Vec::with_capacity(n);
    let boundary_edge_ids: Vec<usize> = (0..n).collect();
    for i in 0..n {
        let p1 = positions[i];
        let p2 = positions[i + 1];
        seg_data.push((p1, p2));
        let cells = UnitSquareGrid::seg_neighborhood_of(p1, p2);
        for &cell in &cells {
            grid.add(cell, i);
        }
    }
    (grid, seg_data, boundary_edge_ids, n)
}

fn register_segment<T: IsComplex + IsRingOrField + Units>(
    grid: &mut UnitSquareGrid,
    p1: T,
    p2: T,
    id: usize,
) {
    let cells = UnitSquareGrid::seg_neighborhood_of(p1, p2);
    for &cell in &cells {
        grid.add(cell, id);
    }
}

fn unregister_segment<T: IsComplex + IsRingOrField + Units>(
    grid: &mut UnitSquareGrid,
    seg_data: &[(T, T)],
    id: usize,
) {
    let (p1, p2) = seg_data[id];
    let cells = UnitSquareGrid::seg_neighborhood_of(p1, p2);
    for &cell in &cells {
        grid.remove(cell, id);
    }
}

fn check_segment_clear<T: IsComplex + IsRingOrField + Units>(
    grid: &UnitSquareGrid,
    seg_data: &[(T, T)],
    p1: T,
    p2: T,
    allowed_endpoint: Option<T>,
) -> bool {
    let cells = UnitSquareGrid::seg_neighborhood_of(p1, p2);
    let near_ids = grid.get_cells(&cells);
    for &id in &near_ids {
        let (x, y) = seg_data[id];
        let is_allowed = allowed_endpoint == Some(p2);
        if !is_allowed && (p2 == x || p2 == y) {
            return false;
        }
        if intersect(&(p1, p2), &(x, y)) {
            return false;
        }
    }
    true
}

fn dir_of_edge<T: IsComplex + IsRingOrField + Units>(from: T, to: T) -> i8 {
    let d = to - from;
    for dir in 0..T::turn() {
        if T::unit(dir) == d {
            return dir;
        }
    }
    panic!("edge is not a unit vector");
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic::{ZZ12, ZZ4};
    use crate::intgeom::snake::Snake;
    use crate::intgeom::tiles;
    use std::collections::BTreeMap;

    fn square_patch() -> GrowingPatch<ZZ4> {
        let sq: Snake<ZZ4> = tiles::square();
        let rat = Rat::try_from(&sq).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        GrowingPatch::new(ts, 0)
    }

    fn hex_patch() -> GrowingPatch<ZZ12> {
        let hex: Snake<ZZ12> = tiles::hexagon();
        let rat = Rat::try_from(&hex).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        GrowingPatch::new(ts, 0)
    }

    #[test]
    fn seed_patch_has_no_boundary() {
        let gp = square_patch();
        assert!(!gp.is_growing());
        assert_eq!(gp.boundary_len(), 0);
        assert_eq!(gp.edges().len(), 0);
    }

    #[test]
    fn seed_patch_has_matches() {
        let gp = square_patch();
        let matches = gp.get_all_matches();
        assert!(!matches.is_empty(), "seed should have matches");
        for pm in &matches {
            assert_eq!(pm.tile_id, 0, "only tile 0 in single-tile set");
            assert!(pm.len > 0, "match length must be positive");
        }
    }

    #[test]
    fn first_add_produces_growing() {
        let mut gp = hex_patch();
        let pm = gp.get_all_matches()[0].clone();
        let diff = gp.add_tile(&pm).expect("first add");

        assert!(gp.is_growing());
        assert_eq!(gp.boundary_len(), 12 - 2 * pm.len);
        assert_eq!(gp.edges().len(), gp.boundary_len());
        assert_eq!(gp.angles().len(), gp.boundary_len());
        let _ = diff;
    }

    #[test]
    fn segment_cyclic_invariant() {
        let mut gp = hex_patch();
        let matches = gp.get_all_matches();
        let pm = matches[0].clone();
        gp.add_tile(&pm).expect("first add");
        assert_eq!(gp.angles().len(), gp.edges().len(), "angles == edges");

        let matches2 = gp.get_all_matches();
        if let Some(pm2) = matches2.first() {
            if gp.add_tile(pm2).is_some() {
                assert_eq!(
                    gp.angles().len(),
                    gp.edges().len(),
                    "angles == edges after second add"
                );
            }
        }
    }

    #[test]
    fn junction_vertex_ids_nonempty_after_each_add() {
        let mut gp = hex_patch();
        let matches = gp.get_all_matches();

        for (step, pm) in matches.iter().enumerate() {
            if gp.add_tile(pm).is_none() {
                break;
            }
            let edges = gp.edges();
            assert!(
                !edges.is_empty(),
                "step {}: edges should not be empty",
                step
            );
            if step < 3 {
                let junctions = gp.junction_vertices();
                assert!(
                    !junctions.is_empty(),
                    "step {}: should have junction vertices",
                    step
                );
            }
        }
    }

    #[test]
    fn hexagon_all_36_matches_produce_valid_bi_hexes() {
        let gp = hex_patch();
        let matches = gp.get_all_matches();
        assert_eq!(matches.len(), 36, "hex self-matches = 36");

        for pm in &matches {
            let mut gp2 = hex_patch();
            assert!(
                gp2.add_tile(pm).is_some(),
                "first add should succeed for pm {:?}",
                pm
            );
            assert!(gp2.is_growing());
            assert_eq!(gp2.boundary_len(), 12 - 2 * pm.len);
            assert_eq!(gp2.edges().len(), gp2.boundary_len());

            let rat = gp2.to_rat();
            assert!(
                Snake::<ZZ12>::try_from(rat.seq()).is_ok(),
                "valid snake for pm {:?}",
                pm
            );
        }
    }

    #[test]
    fn square_all_16_matches_produce_valid_bi_squares() {
        let gp = square_patch();
        let matches = gp.get_all_matches();
        assert_eq!(matches.len(), 16, "square self-matches = 16");

        for pm in &matches {
            let mut gp2 = square_patch();
            assert!(
                gp2.add_tile(pm).is_some(),
                "first add should succeed for pm {:?}",
                pm
            );
            assert!(gp2.is_growing());
            assert_eq!(gp2.boundary_len(), 8 - 2 * pm.len);

            let rat = gp2.to_rat();
            assert!(
                Snake::<ZZ4>::try_from(rat.seq()).is_ok(),
                "valid snake for pm {:?}",
                pm
            );
        }
    }

    #[test]
    fn to_rat_matches_direct_glue_for_all_matches() {
        let gp = hex_patch();
        let matches = gp.get_all_matches();
        let ts = gp.tileset().clone();

        for pm in &matches {
            let mut gp2 = GrowingPatch::<ZZ12>::new(Arc::clone(&ts), 0);
            gp2.add_tile(pm).expect("first add");
            let rat = gp2.to_rat();

            let seed_rat = ts.rat(0);
            let new_rat = ts.rat(pm.tile_id);
            let glued = seed_rat.try_glue((pm.start_a as i64, pm.start_b as i64), new_rat);
            match glued {
                Ok(g) => assert_eq!(rat.seq(), g.seq(), "mismatch for pm {:?}", pm),
                Err(e) => panic!("glue failed for pm {:?}: {}", pm, e),
            }
        }
    }

    #[test]
    fn mixed_hex_square_add_tile() {
        let hex_snake: Snake<ZZ12> = tiles::hexagon();
        let hex_rat = Rat::try_from(&hex_snake).unwrap();
        let sq_snake: Snake<ZZ12> = tiles::square();
        let sq_rat = Rat::try_from(&sq_snake).unwrap();
        let ts = Arc::new(TileSet::new(vec![hex_rat, sq_rat]));
        let mut gp = GrowingPatch::<ZZ12>::new(Arc::clone(&ts), 0);

        let matches = gp.get_all_matches();
        assert!(!matches.is_empty());

        if let Some(pm) = matches.first() {
            let result = gp.add_tile(pm);
            assert!(result.is_some());
            assert!(gp.is_growing());
        }
    }

    #[test]
    fn edges_self_consistent() {
        let gp_sq: GrowingPatch<ZZ4> = square_patch();
        for pm in gp_sq.get_all_matches() {
            let mut gp2 = gp_sq.clone();
            if gp2.add_tile(&pm).is_none() || !gp2.is_growing() {
                continue;
            }
            verify_edges_consistency(&gp2, gp2.tileset(), &format!("bi-sq pm {:?}", pm));
        }
        let gp_hex: GrowingPatch<ZZ12> = hex_patch();
        for pm in gp_hex.get_all_matches() {
            let mut gp2 = gp_hex.clone();
            if gp2.add_tile(&pm).is_none() || !gp2.is_growing() {
                continue;
            }
            verify_edges_consistency(&gp2, gp2.tileset(), &format!("bi-hex pm {:?}", pm));
            for pm2 in gp2.get_all_matches() {
                let mut gp3 = gp2.clone();
                if gp3.add_tile(&pm2).is_some() && gp3.is_growing() {
                    verify_edges_consistency(&gp3, gp3.tileset(), "3-hex");
                }
            }
        }
    }

    fn verify_edges_consistency<T: IsComplex + IsRingOrField + Units>(
        gp: &GrowingPatch<T>,
        ts: &Arc<TileSet<T>>,
        label: &str,
    ) {
        let n = gp.boundary_len();
        assert!(n > 0, "[{}] patch should be growing", label);
        let edges = gp.edges();
        assert_eq!(edges.len(), n, "[{}] edges length", label);

        for (i, edge) in edges.iter().enumerate().take(n) {
            assert!(
                edge.tile_id < ts.num_tiles(),
                "[{}] pos {}: invalid tile_id {}",
                label,
                i,
                edge.tile_id
            );
            let tile_len = ts.rat(edge.tile_id).len();
            assert!(
                edge.tile_offset < tile_len,
                "[{}] pos {}: invalid offset {} for tile {} (len {})",
                label,
                i,
                edge.tile_offset,
                edge.tile_id,
                tile_len
            );
        }

        for i in 0..n {
            let j = (i + 1) % n;
            if edges[i].tile_id == edges[j].tile_id && !gp.is_junction(i) && !gp.is_junction(j) {
                let tile_len = ts.rat(edges[i].tile_id).len();
                let expected_next = (edges[i].tile_offset + 1) % tile_len;
                assert_eq!(
                    edges[j].tile_offset, expected_next,
                    "[{}] pos {}->{}: same-tile continuation expected offset {} got {}",
                    label, i, j, expected_next, edges[j].tile_offset
                );
            }
        }

        let angles = gp.angles();
        assert_eq!(angles.len(), n, "[{}] angles length", label);
    }

    #[test]
    fn edges_bi_hex_consistency() {
        let gp = hex_patch();
        for pm in gp.get_all_matches() {
            let mut gp2 = hex_patch();
            if gp2.add_tile(&pm).is_some() && gp2.is_growing() {
                verify_edges_consistency(&gp2, gp2.tileset(), &format!("bi-hex pm {:?}", pm));
            }
        }
    }

    #[test]
    fn edges_bi_square_consistency() {
        let gp = square_patch();
        for pm in gp.get_all_matches() {
            let mut gp2 = square_patch();
            if gp2.add_tile(&pm).is_some() && gp2.is_growing() {
                verify_edges_consistency(&gp2, gp2.tileset(), &format!("bi-sq pm {:?}", pm));
            }
        }
    }

    #[test]
    fn edges_multi_hex_consistency() {
        let gp = hex_patch();
        let matches = gp.get_all_matches();
        for pm1 in &matches {
            let mut gp2 = hex_patch();
            if gp2.add_tile(pm1).is_none() || !gp2.is_growing() {
                continue;
            }
            verify_edges_consistency(&gp2, gp2.tileset(), "2-hex");
            for pm2 in gp2.get_all_matches() {
                let mut gp3 = gp2.clone();
                if gp3.add_tile(&pm2).is_some() && gp3.is_growing() {
                    verify_edges_consistency(&gp3, gp3.tileset(), "3-hex");
                }
            }
        }
    }

    #[test]
    fn edges_mixed_consistency() {
        let hex_snake: Snake<ZZ12> = tiles::hexagon();
        let sq_snake: Snake<ZZ12> = tiles::square();
        let hex_rat = Rat::try_from(&hex_snake).unwrap();
        let sq_rat = Rat::try_from(&sq_snake).unwrap();
        let ts = Arc::new(TileSet::new(vec![hex_rat, sq_rat]));

        for seed_id in 0..ts.num_tiles() {
            let gp = GrowingPatch::<ZZ12>::new(Arc::clone(&ts), seed_id);
            for pm in gp.get_all_matches() {
                let mut gp2 = GrowingPatch::new(Arc::clone(&ts), seed_id);
                if gp2.add_tile(&pm).is_some() && gp2.is_growing() {
                    verify_edges_consistency(
                        &gp2,
                        &ts,
                        &format!("mixed seed={} pm {:?}", seed_id, pm),
                    );
                }
            }
        }
    }

    #[test]
    fn brute_force_squares_up_to_4_tiles() {
        let sq: Snake<ZZ4> = tiles::square();
        let rat = Rat::try_from(&sq).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let patches = brute_force_zz4(&ts, 4);

        let mut by_tiles: BTreeMap<usize, (usize, usize)> = BTreeMap::new();
        for ways in patches.values() {
            let n = ways[0].len() + 1;
            let e = by_tiles.entry(n).or_insert((0, 0));
            e.0 += 1;
            e.1 += ways.len();
        }

        assert_eq!(
            by_tiles.get(&1).map(|(s, _)| *s).unwrap_or(0),
            1,
            "1 mono-square"
        );
        assert_eq!(by_tiles.get(&2), Some(&(1, 16)), "1 bi-square, 16 ways");
        assert!(
            by_tiles.get(&3).map(|(s, _)| *s).unwrap_or(0) >= 2,
            "at least 2 tri-squares"
        );
    }

    #[test]
    fn brute_force_hexagons_up_to_3_tiles() {
        let hex: Snake<ZZ12> = tiles::hexagon();
        let rat = Rat::try_from(&hex).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let patches = brute_force_zz12(&ts, 3);

        let mut by_tiles: BTreeMap<usize, usize> = BTreeMap::new();
        for ways in patches.values() {
            let n = ways[0].len() + 1;
            by_tiles.entry(n).and_modify(|c| *c += 1).or_insert(1);
        }

        assert_eq!(by_tiles.get(&1).copied().unwrap_or(0), 1, "1 mono-hex");
        assert_eq!(by_tiles.get(&2).copied().unwrap_or(0), 1, "1 bi-hex");
        assert!(
            by_tiles.get(&3).copied().unwrap_or(0) >= 1,
            "at least 1 tri-hex"
        );
    }

    fn brute_force_zz4(
        ts: &Arc<TileSet<ZZ4>>,
        max_tiles: usize,
    ) -> BTreeMap<Rat<ZZ4>, Vec<Vec<PatchMatch>>> {
        let mut results: BTreeMap<Rat<ZZ4>, Vec<Vec<PatchMatch>>> = BTreeMap::new();

        fn recurse(
            gp: &mut GrowingPatch<ZZ4>,
            history: &mut Vec<PatchMatch>,
            max_tiles: usize,
            results: &mut BTreeMap<Rat<ZZ4>, Vec<Vec<PatchMatch>>>,
        ) {
            let num_tiles = history.len() + 1;
            let rat = gp.to_rat();
            results.entry(rat).or_default().push(history.clone());

            if num_tiles >= max_tiles || !gp.is_growing() {
                return;
            }

            let matches = gp.get_all_matches();
            for pm in &matches {
                let mut gp2 = gp.clone();
                if gp2.add_tile(pm).is_some() {
                    history.push(pm.clone());
                    recurse(&mut gp2, history, max_tiles, results);
                    history.pop();
                }
            }
        }

        results
            .entry(ts.rat(0).clone())
            .or_default()
            .push(Vec::new());

        let seed_matches = GrowingPatch::new(Arc::clone(ts), 0).get_all_matches();
        for pm in &seed_matches {
            let mut gp = GrowingPatch::new(Arc::clone(ts), 0);
            gp.add_tile(pm).expect("first add");
            let mut history = vec![pm.clone()];
            recurse(&mut gp, &mut history, max_tiles, &mut results);
        }

        results
    }

    fn brute_force_zz12(
        ts: &Arc<TileSet<ZZ12>>,
        max_tiles: usize,
    ) -> BTreeMap<Rat<ZZ12>, Vec<Vec<PatchMatch>>> {
        let mut results: BTreeMap<Rat<ZZ12>, Vec<Vec<PatchMatch>>> = BTreeMap::new();

        fn recurse(
            gp: &mut GrowingPatch<ZZ12>,
            history: &mut Vec<PatchMatch>,
            max_tiles: usize,
            results: &mut BTreeMap<Rat<ZZ12>, Vec<Vec<PatchMatch>>>,
        ) {
            let num_tiles = history.len() + 1;
            if gp.is_growing() {
                let rat = gp.to_rat();
                results.entry(rat).or_default().push(history.clone());
            }

            if num_tiles >= max_tiles || !gp.is_growing() {
                return;
            }

            let matches = gp.get_all_matches();
            for pm in &matches {
                let mut gp2 = gp.clone();
                if gp2.add_tile(pm).is_some() {
                    history.push(pm.clone());
                    recurse(&mut gp2, history, max_tiles, results);
                    history.pop();
                }
            }
        }

        results
            .entry(ts.rat(0).clone())
            .or_default()
            .push(Vec::new());

        let seed_matches = GrowingPatch::new(Arc::clone(ts), 0).get_all_matches();
        for pm in &seed_matches {
            let mut gp = GrowingPatch::new(Arc::clone(ts), 0);
            gp.add_tile(pm).expect("first add");
            let mut history = vec![pm.clone()];
            recurse(&mut gp, &mut history, max_tiles, &mut results);
        }

        results
    }

    #[test]
    fn inner_chains_empty_after_first_glue() {
        let gp = hex_patch();
        let pm = gp.get_all_matches()[0].clone();
        let mut gp2 = gp.clone();
        gp2.add_tile(&pm).expect("first add");
        let inner = match &gp2.state {
            PatchState::Growing { inner_chains, .. } => inner_chains.clone(),
            _ => panic!("expected Growing"),
        };
        for (i, chain) in inner.iter().enumerate() {
            assert!(
                chain.is_empty(),
                "inner chain at position {i} should be empty after first glue, got {chain:?}"
            );
        }
    }

    #[test]
    fn inner_chains_grow_on_second_glue() {
        let gp = hex_patch();
        let first_match = gp.get_all_matches()[0].clone();
        let mut gp2 = gp.clone();
        gp2.add_tile(&first_match).expect("first add");

        let _n = gp2.boundary_len();
        let candidates = gp2.get_all_matches();
        let second = candidates
            .iter()
            .find(|pm| pm.len == 1)
            .expect("need len-1 match");
        let pos_of_match = second.start_a;
        let mut gp3 = gp2.clone();
        gp3.add_tile(second).expect("second add");

        let inner = match &gp3.state {
            PatchState::Growing { inner_chains, .. } => inner_chains.clone(),
            _ => panic!("expected Growing"),
        };

        let new_n = gp3.boundary_len();
        let mut found_nonempty = false;
        for (i, chain) in inner.iter().enumerate() {
            if !chain.is_empty() {
                found_nonempty = true;
                assert_eq!(
                    chain.len(),
                    1,
                    "inner chain at position {i} should have exactly 1 entry after second glue"
                );
                assert_eq!(
                    chain[0].tile_id, 0,
                    "inner entry should be the original hex tile"
                );
            }
        }
        assert!(found_nonempty, "at least one inner chain should be non-empty after second glue on a {new_n}-edge boundary, pos_of_match={pos_of_match}");
    }

    #[test]
    fn junction_vertex_type_roundtrip_after_first_glue() {
        let gp = hex_patch();
        let pm = gp.get_all_matches()[0].clone();
        let mut gp2 = gp.clone();
        gp2.add_tile(&pm).expect("first add");
        let n = gp2.boundary_len();
        let mut junction_count = 0;
        for i in 0..n {
            if let Some(vt) = gp2.junction_vertex_type_at(i) {
                assert!(vt.inner.is_empty(), "inner should be empty at pos {i}");
                junction_count += 1;
            }
        }
        assert!(junction_count > 0, "should have at least one junction");
    }

    #[test]
    fn construct_minimal_witness_hex_roundtrip() {
        let gp = hex_patch();
        let matches = gp.get_all_matches();
        let ts = gp.match_index.clone();

        for pm in &matches {
            let mut glued = gp.clone();
            glued.add_tile(pm).expect("glue should succeed");
            let n = glued.boundary_len();
            for pos in 0..n {
                let vt = match glued.junction_vertex_type_at(pos) {
                    Some(vt) => vt,
                    None => continue,
                };
                let result = GrowingPatch::construct_minimal_witness(&vt, Arc::clone(&ts));
                assert!(
                    result.is_some(),
                    "construct_minimal_witness should succeed for vt={vt:?} from pm={pm:?} pos={pos}"
                );
                let (witness, wpos) = result.unwrap();
                let w_edges = witness.edges();
                let w_inner = witness.inner_chains();
                let reconstructed = vertex_type_raw_from(w_edges, w_inner, wpos);
                assert_eq!(
                    reconstructed, vt,
                    "roundtrip failed for vt={vt:?} from pm={pm:?} pos={pos}"
                );
            }
        }
    }

    #[test]
    fn construct_minimal_witness_square_roundtrip() {
        let gp = square_patch();
        let matches = gp.get_all_matches();
        let ts = gp.match_index.clone();

        for pm in &matches {
            let mut glued = gp.clone();
            glued.add_tile(pm).expect("glue should succeed");
            let n = glued.boundary_len();
            for pos in 0..n {
                let vt = match glued.junction_vertex_type_at(pos) {
                    Some(vt) => vt,
                    None => continue,
                };
                let result = GrowingPatch::construct_minimal_witness(&vt, Arc::clone(&ts));
                assert!(
                    result.is_some(),
                    "construct_minimal_witness should succeed for vt={vt:?} from pm={pm:?} pos={pos}"
                );
                let (witness, wpos) = result.unwrap();
                let w_edges = witness.edges();
                let w_inner = witness.inner_chains();
                let reconstructed = vertex_type_raw_from(w_edges, w_inner, wpos);
                assert_eq!(
                    reconstructed, vt,
                    "roundtrip failed for vt={vt:?} from pm={pm:?} pos={pos}"
                );
            }
        }
    }

    #[test]
    fn construct_minimal_witness_hex_with_inner() {
        let gp = hex_patch();
        let ts = gp.match_index.clone();
        let first = gp.get_all_matches()[0].clone();
        let mut gp2 = gp.clone();
        gp2.add_tile(&first).expect("first add");

        let second_candidates = gp2.get_all_matches();
        let len1_match = second_candidates
            .iter()
            .find(|pm| pm.len == 1)
            .expect("need len-1 match")
            .clone();
        let mut gp3 = gp2.clone();
        gp3.add_tile(&len1_match).expect("second add");

        for pos in 0..gp3.boundary_len() {
            let vt = match gp3.junction_vertex_type_at(pos) {
                Some(vt) => vt,
                None => continue,
            };
            let result = GrowingPatch::construct_minimal_witness(&vt, Arc::clone(&ts));
            assert!(
                result.is_some(),
                "construct_minimal_witness should succeed for vt={vt:?} with inner"
            );
            let (witness, wpos) = result.unwrap();
            let w_edges = witness.edges();
            let w_inner = witness.inner_chains();
            let reconstructed = vertex_type_raw_from(w_edges, w_inner, wpos);
            assert_eq!(reconstructed, vt, "roundtrip failed for vt={vt:?}");
        }
    }

    #[test]
    fn get_matches_touching_vertex_lazy_matches_eager() {
        let ts: Arc<TileSet<ZZ12>> =
            Arc::new(TileSet::new(
                vec![Rat::try_from(&tiles::spectre()).unwrap()],
            ));
        let mut gp = GrowingPatch::new(Arc::clone(&ts), 0);
        let pm = gp.get_all_matches().into_iter().next().unwrap();
        gp.add_tile(&pm).unwrap();
        gp.ensure_candidates_materialized();

        let n = gp.boundary_len();
        for target in 0..n {
            let lazy = {
                let gp2 = gp.clone_for_mutation();
                gp2.get_matches_touching_vertex(target)
            };
            let eager = gp.get_matches_touching_vertex(target);
            let mut lazy_sorted = lazy;
            lazy_sorted.sort_by_key(|pm| (pm.start_a, pm.len, pm.start_b, pm.tile_id));
            let mut eager_sorted = eager;
            eager_sorted.sort_by_key(|pm| (pm.start_a, pm.len, pm.start_b, pm.tile_id));
            assert_eq!(
                lazy_sorted,
                eager_sorted,
                "Mismatch at target={target}: lazy={}, eager={}",
                lazy_sorted.len(),
                eager_sorted.len()
            );
        }
    }

    #[test]
    fn compute_candidates_covering_position_is_subset_of_all_candidates() {
        let ts: Arc<TileSet<ZZ12>> =
            Arc::new(TileSet::new(
                vec![Rat::try_from(&tiles::spectre()).unwrap()],
            ));
        let mi: Arc<MatchTypeIndex<ZZ12>> = Arc::new(MatchTypeIndex::new(Arc::clone(&ts)));
        let mut gp = GrowingPatch::new(Arc::clone(&ts), 0);
        let pm = gp.get_all_matches().into_iter().next().unwrap();
        gp.add_tile(&pm).unwrap();

        let all_cands = GrowingPatch::compute_all_candidates(&mi, gp.angles(), gp.edges());
        let n = gp.angles().len();
        let max_tile_len = ts.rats().iter().map(|r| r.len()).max().unwrap_or(0);
        let k = max_tile_len.min(n);

        for target in 0..n {
            let covering = GrowingPatch::compute_candidates_covering_position(
                &mi,
                gp.angles(),
                gp.edges(),
                target,
            );

            let mut touching_from_all: Vec<_> = Vec::new();
            for offset in 0..=k {
                let start = (target + n - offset) % n;
                for pm in &all_cands[start] {
                    if cyclic_range_contains(pm.start_a, pm.len, target, n) {
                        touching_from_all.push(pm.clone());
                    }
                }
            }

            for pm in &covering {
                assert!(
                    touching_from_all.contains(pm),
                    "covering returned extra match at target={target}: {:?}",
                    pm
                );
            }
        }
    }

    #[test]
    fn construct_minimal_witness_hex_boundary_matches_brute_force() {
        let gp = hex_patch();
        let ts = gp.match_index.clone();

        for pm in &gp.get_all_matches() {
            let mut brute = gp.clone();
            brute.add_tile(pm).expect("brute glue");
            let brute_angles = brute.angles().to_vec();
            let brute_edges = brute.edges().to_vec();
            let brute_inner = brute.inner_chains().to_vec();

            for pos in 0..brute.boundary_len() {
                let vt = match brute.junction_vertex_type_at(pos) {
                    Some(vt) => vt,
                    None => continue,
                };

                let (witness, _wpos) =
                    GrowingPatch::construct_minimal_witness(&vt, Arc::clone(&ts))
                        .expect("reconstruction should succeed");

                let w_edges = witness.edges();
                let w_inner = witness.inner_chains();

                let mut found = false;
                for wpos in 0..witness.boundary_len() {
                    let wvt = vertex_type_raw_from(w_edges, w_inner, wpos);
                    if wvt == vt {
                        let brute_vt = vertex_type_raw_from(&brute_edges, &brute_inner, pos);
                        assert_eq!(
                            wvt, brute_vt,
                            "reconstructed VT should match brute-force VT"
                        );
                        found = true;
                        break;
                    }
                }
                assert!(found, "should find matching position in witness");

                let mut wa: Vec<i8> = witness.angles().to_vec();
                let mut ba = brute_angles.clone();
                wa.sort();
                ba.sort();
                assert_eq!(
                    wa, ba,
                    "reconstructed angles should be a rotation of brute-force angles"
                );
            }
        }
    }

    #[test]
    fn construct_minimal_witness_square_boundary_matches_brute_force() {
        let gp = square_patch();
        let ts = gp.match_index.clone();

        for pm in &gp.get_all_matches() {
            let mut brute = gp.clone();
            brute.add_tile(pm).expect("brute glue");
            let brute_edges = brute.edges().to_vec();
            let brute_inner = brute.inner_chains().to_vec();

            for pos in 0..brute.boundary_len() {
                let vt = match brute.junction_vertex_type_at(pos) {
                    Some(vt) => vt,
                    None => continue,
                };

                let (witness, _wpos) =
                    GrowingPatch::construct_minimal_witness(&vt, Arc::clone(&ts))
                        .expect("reconstruction should succeed");

                let w_edges = witness.edges();
                let w_inner = witness.inner_chains();

                let mut found = false;
                for wpos in 0..witness.boundary_len() {
                    let wvt = vertex_type_raw_from(w_edges, w_inner, wpos);
                    if wvt == vt {
                        let brute_vt = vertex_type_raw_from(&brute_edges, &brute_inner, pos);
                        assert_eq!(wvt, brute_vt);
                        found = true;
                        break;
                    }
                }
                assert!(
                    found,
                    "should find matching position in witness for vt={vt:?}"
                );
            }
        }
    }

    #[test]
    fn construct_minimal_witness_spectre_roundtrip() {
        let ts: Arc<TileSet<ZZ12>> =
            Arc::new(TileSet::new(
                vec![Rat::try_from(&tiles::spectre()).unwrap()],
            ));
        let mut gp = GrowingPatch::new(Arc::clone(&ts), 0);

        let pm = gp.get_all_matches().into_iter().next().unwrap();
        gp.add_tile(&pm).expect("first spectre glue");

        let mi = gp.match_index.clone();
        for pos in 0..gp.boundary_len() {
            let vt = match gp.junction_vertex_type_at(pos) {
                Some(vt) => vt,
                None => continue,
            };
            let result = GrowingPatch::construct_minimal_witness(&vt, Arc::clone(&mi));
            assert!(
                result.is_some(),
                "spectre witness should succeed for vt={vt:?} at pos={pos}"
            );
            let (witness, wpos) = result.unwrap();
            let w_edges = witness.edges();
            let w_inner = witness.inner_chains();
            let reconstructed = vertex_type_raw_from(w_edges, w_inner, wpos);
            assert_eq!(reconstructed, vt, "spectre roundtrip failed for vt={vt:?}");
        }
    }

    #[test]
    fn forward_match_length_hex_basic() {
        let hex: Snake<ZZ12> = tiles::hexagon();
        let rat = Rat::try_from(&hex).unwrap();
        let seq = rat.seq();

        assert_eq!(forward_match_length(seq, 0, seq, 0), 1);
        assert_eq!(forward_match_length(seq, 3, seq, 3), 1);
        assert_eq!(forward_match_length(seq, 0, seq, 1), 1);

        let boundary: Vec<i8> = vec![-2, 2, 2, 2, 2, -2, 2, 2, 2, 2];
        assert_eq!(forward_match_length(&boundary, 5, seq, 0), 1);
        assert_eq!(forward_match_length(&boundary, 0, seq, 0), 1);
    }

    #[test]
    fn forward_match_length_square_basic() {
        let sq: Snake<ZZ4> = tiles::square();
        let rat = Rat::try_from(&sq).unwrap();
        let seq = rat.seq();

        assert_eq!(forward_match_length(seq, 0, seq, 0), 1);
        assert_eq!(forward_match_length(seq, 2, seq, 2), 1);
    }

    #[test]
    fn glue_raw_angles_hex_self_glue() {
        let hex: Snake<ZZ12> = tiles::hexagon();
        let rat = Rat::try_from(&hex).unwrap();
        let seq = rat.seq().to_vec();

        let result = angles::glue_raw_angles::<ZZ12>(&seq, &seq, 0, 1, 0);
        assert!(result.is_some());
        let gr = result.unwrap();
        assert_eq!(gr.angles.len(), 10);
        assert_eq!(gr.a_yx, Some(-2));
        assert_eq!(gr.a_xy, Some(-2));
    }

    #[test]
    fn glue_raw_angles_matches_rat_glue() {
        let hex: Snake<ZZ12> = tiles::hexagon();
        let rat = Rat::try_from(&hex).unwrap();
        let seq = rat.seq();

        let rat_result = rat.try_glue((0, 0), &rat).expect("rat glue");
        let raw_result = angles::glue_raw_angles::<ZZ12>(seq, seq, 0, 1, 0).expect("raw glue");

        let rat_seq = rat_result.seq();
        let raw_seq = &raw_result.angles;

        let mut rat_sorted: Vec<i8> = rat_seq.to_vec();
        let mut raw_sorted = raw_seq.clone();
        rat_sorted.sort();
        raw_sorted.sort();
        assert_eq!(rat_sorted, raw_sorted);
    }
}
