use std::collections::BTreeSet;
use std::sync::Arc;

use crate::cyclotomic::{IsComplex, IsRingOrField, Units};
use crate::intgeom::angles::normalize_angle;
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

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct PatchVertexType {
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
        candidates_by_start: Vec<Vec<PatchMatch>>,
    },
    _Phantom(std::marker::PhantomData<T>),
}

impl<T: IsComplex + IsRingOrField + Units> GrowingPatch<T> {
    pub fn new(tileset: Arc<TileSet<T>>, seed_tile_id: usize) -> Self {
        let match_index = Arc::new(MatchTypeIndex::new(Arc::clone(&tileset)));
        let seed_matches =
            Self::compute_seed_matches(&match_index, &match_index.tileset(), seed_tile_id);
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
        candidates_by_start: Vec<Vec<PatchMatch>>,
    ) -> Option<Self> {
        if angles.is_empty()
            || angles.len() != edges.len()
            || candidates_by_start.len() != angles.len()
        {
            return None;
        }
        Some(GrowingPatch {
            match_index,
            state: PatchState::Growing {
                angles,
                edges,
                candidates_by_start,
            },
        })
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

        let mut global_seen: BTreeSet<(usize, usize, usize, usize)> = BTreeSet::new();

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
                        result[junc_idx].push(PatchMatch {
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
                    candidates_by_start,
                    ..
                } => PatchState::Growing {
                    angles: angles.clone(),
                    edges: edges.clone(),
                    candidates_by_start: candidates_by_start.clone(),
                },
                PatchState::_Phantom(p) => PatchState::_Phantom(p.clone()),
            },
        }
    }

    pub fn get_all_matches(&self) -> Vec<PatchMatch> {
        match &self.state {
            PatchState::Seed { cached_matches, .. } => cached_matches.clone(),
            PatchState::Growing {
                candidates_by_start,
                ..
            } => candidates_by_start.iter().flatten().cloned().collect(),
            PatchState::_Phantom(_) => Vec::new(),
        }
    }

    pub fn get_matches_for_tile(&self, tile_id: usize) -> impl Iterator<Item = PatchMatch> + '_ {
        self.get_all_matches()
            .into_iter()
            .filter(move |m| m.tile_id == tile_id)
    }

    pub fn get_matches_touching_vertex(
        &self,
        vertex_index: usize,
    ) -> impl Iterator<Item = &PatchMatch> {
        match &self.state {
            PatchState::Seed { cached_matches, .. } => {
                let n = self.match_index.tileset().rat(0).len();
                cached_matches
                    .iter()
                    .filter(move |pm| cyclic_range_contains(pm.start_a, pm.len, vertex_index, n))
                    .collect::<Vec<_>>()
                    .into_iter()
            }
            PatchState::Growing {
                candidates_by_start,
                angles,
                ..
            } => {
                let n = angles.len();
                let k = self.max_tile_len().min(n);
                let mut result = Vec::new();
                for offset in 0..=k {
                    let start = (vertex_index + n - offset) % n;
                    for pm in &candidates_by_start[start] {
                        if cyclic_range_contains(pm.start_a, pm.len, vertex_index, n) {
                            result.push(pm);
                        }
                    }
                }
                result.into_iter()
            }
            PatchState::_Phantom(_) => Vec::new().into_iter(),
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

    pub fn candidates_by_start(&self) -> &[Vec<PatchMatch>] {
        match &self.state {
            PatchState::Growing {
                candidates_by_start,
                ..
            } => candidates_by_start,
            _ => &[],
        }
    }

    pub fn vertex_type_at(&self, i: usize) -> Option<PatchVertexType> {
        let n = self.boundary_len();
        if n == 0 || i >= n {
            return None;
        }
        let edges = match &self.state {
            PatchState::Growing { edges, .. } => edges,
            _ => return None,
        };
        let angles = self.angles();
        Some(PatchVertexType {
            angle: angles[i],
            cw: edges[(i + n - 1) % n],
            ccw: edges[i],
        })
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

    pub fn junction_vertices(&self) -> Vec<(usize, PatchVertexType)> {
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

        let candidates_by_start =
            Self::compute_all_candidates(&self.match_index, &new_angles, &edges);

        self.state = PatchState::Growing {
            angles: new_angles,
            edges,
            candidates_by_start,
        };

        Some(AddTileDiff)
    }

    fn add_tile_growing(&mut self, pm: &PatchMatch) -> Option<AddTileDiff> {
        let (angles, edges, old_candidates) = match &mut self.state {
            PatchState::Growing {
                angles,
                edges,
                candidates_by_start,
            } => (
                std::mem::take(angles),
                std::mem::take(edges),
                std::mem::take(candidates_by_start),
            ),
            _ => return None,
        };

        let n = angles.len();
        let mlen = pm.len;
        let tileset = self.match_index.tileset();
        let m = tileset.rat(pm.tile_id).seq().len();

        if mlen == 0 || mlen > n || mlen > m {
            self.restore_growing(angles, edges, old_candidates);
            return None;
        }

        let new_angles = match compute_glue_angles::<T>(&angles, pm, tileset) {
            Some(a) => a,
            None => {
                self.restore_growing(angles, edges, old_candidates);
                return None;
            }
        };

        let seg_len_old = n - mlen;
        let seg_len_new = m - mlen;
        let new_len = seg_len_old + seg_len_new;

        if new_len == 0 {
            return None;
        }

        if Snake::<T>::try_from(new_angles.as_slice()).is_err() {
            self.restore_growing(angles, edges, old_candidates);
            return None;
        }

        let ccw_pos = (pm.start_a + mlen) % n;

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

        let new_candidates =
            Self::compute_all_candidates(&self.match_index, &new_angles, &new_edges);

        self.state = PatchState::Growing {
            angles: new_angles,
            edges: new_edges,
            candidates_by_start: new_candidates,
        };

        Some(AddTileDiff)
    }

    fn restore_growing(
        &mut self,
        angles: Vec<i8>,
        edges: Vec<EdgeInfo>,
        candidates_by_start: Vec<Vec<PatchMatch>>,
    ) {
        self.state = PatchState::Growing {
            angles,
            edges,
            candidates_by_start,
        };
    }

    fn max_tile_len(&self) -> usize {
        self.match_index
            .tileset()
            .rats()
            .iter()
            .map(|r| r.len())
            .max()
            .unwrap_or(0)
    }

    fn incremental_candidates(
        old_angles: &[i8],
        _old_edges: &[EdgeInfo],
        old_candidates: &[Vec<PatchMatch>],
        new_angles: &[i8],
        new_edges: &[EdgeInfo],
        start_a: usize,
        mlen: usize,
        match_index: &Arc<MatchTypeIndex<T>>,
    ) -> Vec<Vec<PatchMatch>> {
        let n = old_angles.len();
        let n_new = new_angles.len();
        let k = match_index
            .tileset()
            .rats()
            .iter()
            .map(|r| r.len())
            .max()
            .unwrap_or(0);

        let mut new_candidates = vec![Vec::new(); n_new];

        let seg_len_old = n - mlen;

        for new_pos in 0..n_new {
            let is_old_part = new_pos < seg_len_old;
            let is_junction_pos = new_pos == 0 || new_pos == seg_len_old;

            if is_old_part && !is_junction_pos {
                let old_pos = (start_a + mlen + new_pos) % n;
                let offset_from_junction_start = new_pos;
                let offset_from_junction_end = seg_len_old - new_pos;

                if offset_from_junction_start > k && offset_from_junction_end > k {
                    new_candidates[new_pos] = old_candidates[old_pos].clone();
                    continue;
                }
            }

            new_candidates[new_pos] =
                Self::compute_candidates_at(new_pos, new_angles, new_edges, match_index);
        }

        new_candidates
    }

    fn compute_candidates_at(
        pos: usize,
        angles: &[i8],
        edges: &[EdgeInfo],
        match_index: &MatchTypeIndex<T>,
    ) -> Vec<PatchMatch> {
        let n = angles.len();
        let rat = Rat::from_slice_unchecked(angles);
        let tileset = match_index.tileset();
        let tile_id = edges[pos].tile_id;
        let tile_offset = edges[pos].tile_offset;

        let mut seen: BTreeSet<(usize, usize, usize, usize)> = BTreeSet::new();
        let mut candidates = Vec::new();

        for cand in match_index.candidates_starting_at(tile_id, tile_offset) {
            append_match_candidate(
                pos,
                angles,
                &rat,
                n,
                cand,
                match_index,
                tileset,
                &mut seen,
                &mut candidates,
            );
        }

        let ei = edges[pos];
        if tileset.rat(ei.tile_id).seq()[ei.tile_offset] != angles[pos] {
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
                    if !seen.insert(key) {
                        continue;
                    }
                    if rat
                        .try_glue_precomputed((ns, len, ne), tile_b, true)
                        .is_ok()
                    {
                        candidates.push(PatchMatch {
                            start_a: ns_u,
                            len,
                            start_b: ne_u,
                            tile_id: tile_id_b,
                        });
                    }
                }
            }
        }

        candidates
    }

    fn compute_seed_matches(
        match_index: &MatchTypeIndex<T>,
        tileset: &TileSet<T>,
        seed_tile_id: usize,
    ) -> Vec<PatchMatch> {
        let seed = tileset.rat(seed_tile_id);
        let seed_seq = seed.seq();
        let n = seed_seq.len();
        let mut seen: BTreeSet<(usize, usize, usize, usize)> = BTreeSet::new();
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
                if let Ok(_glued) =
                    seed.try_glue_precomputed((ns as i64, len, ne as i64), tile_b, true)
                {
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

fn compute_candidates_at_position<T: IsComplex + IsRingOrField + Units>(
    pos: usize,
    angles: &[i8],
    rat: &Rat<T>,
    n: usize,
    tile_id: usize,
    tile_offset: usize,
    match_index: &MatchTypeIndex<T>,
    tileset: &TileSet<T>,
    global_seen: &mut BTreeSet<(usize, usize, usize, usize)>,
    result: &mut [Vec<PatchMatch>],
) {
    for cand in match_index.candidates_starting_at(tile_id, tile_offset) {
        append_match_candidate(
            pos,
            angles,
            rat,
            n,
            cand,
            match_index,
            tileset,
            global_seen,
            &mut result[pos],
        );
    }
}

fn append_match_candidate<T: IsComplex + IsRingOrField + Units>(
    pos: usize,
    angles: &[i8],
    rat: &Rat<T>,
    n: usize,
    cand: &CandidateMatch,
    _match_index: &MatchTypeIndex<T>,
    tileset: &TileSet<T>,
    seen: &mut BTreeSet<(usize, usize, usize, usize)>,
    candidates: &mut Vec<PatchMatch>,
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
        candidates.push(PatchMatch {
            start_a: ns_u,
            len,
            start_b: ne_u,
            tile_id: cand.tile_b,
        });
    }
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

pub(crate) fn compute_glue_angles<T: IsComplex + IsRingOrField + Units>(
    angles: &[i8],
    pm: &PatchMatch,
    tileset: &Arc<TileSet<T>>,
) -> Option<Vec<i8>> {
    let n = angles.len();
    let mlen = pm.len;
    let seed_seq = tileset.rat(pm.tile_id).seq();
    let m = seed_seq.len();

    let x_raw_len = n - mlen + 1;
    let y_raw_len = m - mlen + 1;

    if x_raw_len <= 1 && y_raw_len <= 1 {
        return None;
    }

    let x: Vec<i8> = (0..x_raw_len)
        .map(|i| angles[(pm.start_a + mlen + i) % n])
        .collect();
    let y: Vec<i8> = (0..y_raw_len)
        .map(|i| seed_seq[(pm.start_b + i) % m])
        .collect();

    let seg_len_old = x_raw_len - 1;
    let seg_len_new = y_raw_len - 1;

    let mut new_angles = Vec::with_capacity(seg_len_old + seg_len_new);
    if seg_len_old > 0 {
        new_angles.extend_from_slice(&x[..seg_len_old]);
    }
    if seg_len_new > 0 {
        new_angles.extend_from_slice(&y[..seg_len_new]);
    }

    if new_angles.is_empty() {
        return None;
    }

    let a_yx = if x_raw_len > 1 && y_raw_len > 1 {
        normalize_angle::<T>(x[0] + y[y_raw_len - 1] - T::hturn())
    } else {
        new_angles[0]
    };
    let a_xy = if x_raw_len > 1 && y_raw_len > 1 {
        normalize_angle::<T>(y[0] + x[x_raw_len - 1] - T::hturn())
    } else {
        new_angles[0]
    };

    if a_yx.abs() == T::hturn() || a_xy.abs() == T::hturn() {
        return None;
    }

    if seg_len_old > 0 {
        new_angles[0] = a_yx;
    }
    if seg_len_old > 0 && seg_len_new > 0 {
        new_angles[seg_len_old] = a_xy;
    } else if seg_len_new > 0 {
        new_angles[0] = a_xy;
    }

    Some(new_angles)
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
            verify_edges_consistency(&gp2, &gp2.tileset(), &format!("bi-sq pm {:?}", pm));
        }
        let gp_hex: GrowingPatch<ZZ12> = hex_patch();
        for pm in gp_hex.get_all_matches() {
            let mut gp2 = gp_hex.clone();
            if gp2.add_tile(&pm).is_none() || !gp2.is_growing() {
                continue;
            }
            verify_edges_consistency(&gp2, &gp2.tileset(), &format!("bi-hex pm {:?}", pm));
            for pm2 in gp2.get_all_matches() {
                let mut gp3 = gp2.clone();
                if gp3.add_tile(&pm2).is_some() && gp3.is_growing() {
                    verify_edges_consistency(&gp3, &gp3.tileset(), "3-hex");
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

        for i in 0..n {
            assert!(
                edges[i].tile_id < ts.num_tiles(),
                "[{}] pos {}: invalid tile_id {}",
                label,
                i,
                edges[i].tile_id
            );
            let tile_len = ts.rat(edges[i].tile_id).len();
            assert!(
                edges[i].tile_offset < tile_len,
                "[{}] pos {}: invalid offset {} for tile {} (len {})",
                label,
                i,
                edges[i].tile_offset,
                edges[i].tile_id,
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
                verify_edges_consistency(&gp2, &gp2.tileset(), &format!("bi-hex pm {:?}", pm));
            }
        }
    }

    #[test]
    fn edges_bi_square_consistency() {
        let gp = square_patch();
        for pm in gp.get_all_matches() {
            let mut gp2 = square_patch();
            if gp2.add_tile(&pm).is_some() && gp2.is_growing() {
                verify_edges_consistency(&gp2, &gp2.tileset(), &format!("bi-sq pm {:?}", pm));
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
            verify_edges_consistency(&gp2, &gp2.tileset(), "2-hex");
            for pm2 in gp2.get_all_matches() {
                let mut gp3 = gp2.clone();
                if gp3.add_tile(&pm2).is_some() && gp3.is_growing() {
                    verify_edges_consistency(&gp3, &gp3.tileset(), "3-hex");
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
}
