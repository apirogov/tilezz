use std::sync::Arc;

use crate::cyclotomic::{IsComplex, IsRingOrField, Units};
use crate::intgeom::angles::normalize_angle;
use crate::intgeom::matchtypes::{MatchFinder, MatchType, MatchTypeIndex};
use crate::intgeom::rat::Rat;
use crate::intgeom::tileset::TileSet;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BoundarySegment {
    pub tile_id: usize,
    pub tile_start: usize,
    pub tile_len: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VertexInfo {
    pub matches: Vec<i32>,
    pub is_open: bool,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct VertexType(pub Vec<i32>);

#[derive(Debug, Clone)]
pub struct PatchMatch {
    pub start_a: usize,
    pub len: usize,
    pub start_b: usize,
    pub tile_id: usize,
}

pub struct AddTileDiff {
    pub closed_vertices: Vec<VertexInfo>,
}

#[derive(Clone)]
pub struct GrowingPatch<T: IsComplex> {
    tileset: Arc<TileSet<T>>,
    match_type_index: Arc<MatchTypeIndex<T>>,
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
        segments: Vec<BoundarySegment>,
        vertices: Vec<VertexInfo>,
        cached_matches: Vec<PatchMatch>,
    },
    Closed,
    _Phantom(std::marker::PhantomData<T>),
}

impl<T: IsComplex + IsRingOrField + Units> GrowingPatch<T> {
    pub fn new(
        tileset: Arc<TileSet<T>>,
        match_type_index: Arc<MatchTypeIndex<T>>,
        seed_tile_id: usize,
    ) -> Self {
        let seed_matches = Self::compute_seed_matches(&tileset, seed_tile_id);
        GrowingPatch {
            tileset,
            match_type_index,
            state: PatchState::Seed {
                tile_id: seed_tile_id,
                cached_matches: seed_matches,
            },
        }
    }

    pub fn seed_tile_id(&self) -> usize {
        match &self.state {
            PatchState::Seed { tile_id, .. } => *tile_id,
            PatchState::Growing { segments, .. } => {
                segments.first().map(|s| s.tile_id).unwrap_or(0)
            }
            PatchState::Closed | PatchState::_Phantom(_) => 0,
        }
    }

    pub fn is_growing(&self) -> bool {
        matches!(self.state, PatchState::Growing { .. })
    }

    pub fn is_closed(&self) -> bool {
        matches!(self.state, PatchState::Closed)
    }

    pub fn get_all_matches(&self) -> &[PatchMatch] {
        match &self.state {
            PatchState::Seed { cached_matches, .. } => cached_matches,
            PatchState::Growing { cached_matches, .. } => cached_matches,
            PatchState::Closed | PatchState::_Phantom(_) => &[],
        }
    }

    pub fn get_matches_for_tile(&self, tile_id: usize) -> impl Iterator<Item = &PatchMatch> {
        self.get_all_matches()
            .iter()
            .filter(move |m| m.tile_id == tile_id)
    }

    pub fn num_tiles(&self) -> usize {
        self.tileset.num_tiles()
    }

    pub fn boundary_len(&self) -> usize {
        match &self.state {
            PatchState::Growing { angles, .. } => angles.len(),
            _ => 0,
        }
    }

    pub fn segments(&self) -> &[BoundarySegment] {
        match &self.state {
            PatchState::Growing { segments, .. } => segments,
            _ => &[],
        }
    }

    pub fn vertices(&self) -> &[VertexInfo] {
        match &self.state {
            PatchState::Growing { vertices, .. } => vertices,
            _ => &[],
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

    pub fn to_vertices(&self) -> Vec<VertexInfo> {
        match &self.state {
            PatchState::Growing { vertices, .. } => vertices.clone(),
            _ => Vec::new(),
        }
    }

    pub fn to_segments(&self) -> Vec<BoundarySegment> {
        match &self.state {
            PatchState::Growing { segments, .. } => segments.clone(),
            _ => Vec::new(),
        }
    }

    pub fn tileset(&self) -> &Arc<TileSet<T>> {
        &self.tileset
    }

    pub fn add_tile(&mut self, pm: &PatchMatch) -> Option<AddTileDiff> {
        match &self.state {
            PatchState::Seed { tile_id, .. } => {
                let seed_id = *tile_id;
                self.init_from_first_add(seed_id, pm)?;
                Some(AddTileDiff {
                    closed_vertices: Vec::new(),
                })
            }
            PatchState::Growing { .. } => self.add_tile_growing(pm),
            PatchState::Closed | PatchState::_Phantom(_) => None,
        }
    }

    fn init_from_first_add(&mut self, seed_id: usize, pm: &PatchMatch) -> Option<()> {
        let seed_rat = self.tileset.rat(seed_id);
        let n = seed_rat.seq().len();
        let m = self.tileset.rat(pm.tile_id).seq().len();
        let mlen = pm.len;

        if mlen == 0 || mlen > n || mlen > m {
            return None;
        }

        let (ns, seed_len, ne) = seed_rat.get_match(
            (pm.start_a as i64, pm.start_b as i64),
            self.tileset.rat(pm.tile_id),
        );
        if seed_len == 0 {
            return None;
        }
        let signed_id = self.match_type_index.signed_id(&MatchType {
            tile_a: seed_id,
            start_a: ns as usize,
            tile_b: pm.tile_id,
            start_b: ne as usize,
            len: seed_len,
        });

        let seed_angles = seed_rat.seq().to_vec();
        let new_angles = compute_glue_angles::<T>(&seed_angles, pm, &self.tileset)?;

        let seg_len_old = n - mlen;
        let seg_len_new = m - mlen;

        let mut segments = Vec::new();
        let mut vertices = Vec::new();

        if seg_len_old > 0 {
            segments.push(BoundarySegment {
                tile_id: seed_id,
                tile_start: (pm.start_a + mlen) % n,
                tile_len: seg_len_old,
            });
        }
        if seg_len_new > 0 {
            segments.push(BoundarySegment {
                tile_id: pm.tile_id,
                tile_start: (pm.start_b + mlen) % m,
                tile_len: seg_len_new,
            });
        }

        if segments.is_empty() {
            return None;
        }

        for i in 0..segments.len() {
            let is_cw = i == 0 || segments.len() == 1;
            let is_ccw = i == segments.len() - 1 || segments.len() == 1;
            let mut matches = Vec::new();
            if let Some(sid) = signed_id {
                if is_cw {
                    matches.push(-sid);
                }
                if is_ccw && !(is_cw && sid == -sid) {
                    matches.push(sid);
                }
            }
            vertices.push(VertexInfo {
                matches,
                is_open: true,
            });
        }

        debug_assert_eq!(segments.len(), vertices.len());

        self.state = PatchState::Growing {
            angles: new_angles,
            segments,
            vertices,
            cached_matches: Vec::new(),
        };
        self.recompute_matches();
        Some(())
    }

    fn add_tile_growing(&mut self, pm: &PatchMatch) -> Option<AddTileDiff> {
        let (angles, segments, vertices) = match &mut self.state {
            PatchState::Growing {
                angles,
                segments,
                vertices,
                ..
            } => (
                std::mem::take(angles),
                std::mem::take(segments),
                std::mem::take(vertices),
            ),
            _ => return None,
        };

        let n = angles.len();
        let mlen = pm.len;
        let m = self.tileset.rat(pm.tile_id).seq().len();

        if mlen == 0 || mlen > n || mlen > m {
            self.restore_growing(angles, segments, vertices);
            return None;
        }

        let signed_id = resolve_signed_id(
            &segments,
            &angles,
            pm,
            &self.tileset,
            &self.match_type_index,
        );

        let new_angles = match compute_glue_angles::<T>(&angles, pm, &self.tileset) {
            Some(a) => a,
            None => {
                self.restore_growing(angles, segments, vertices);
                return None;
            }
        };

        let result = apply_match_locally(
            &segments,
            &vertices,
            pm,
            &self.tileset,
            &self.match_type_index,
            signed_id,
        );

        match result {
            Some((new_segments, new_vertices, closed_vertices)) => {
                let has_open = new_vertices.iter().any(|v| v.is_open);
                self.state = if has_open {
                    PatchState::Growing {
                        angles: new_angles,
                        segments: new_segments,
                        vertices: new_vertices,
                        cached_matches: Vec::new(),
                    }
                } else {
                    PatchState::Closed
                };
                self.recompute_matches();
                Some(AddTileDiff { closed_vertices })
            }
            None => {
                self.restore_growing(angles, segments, vertices);
                None
            }
        }
    }

    fn restore_growing(
        &mut self,
        angles: Vec<i8>,
        segments: Vec<BoundarySegment>,
        vertices: Vec<VertexInfo>,
    ) {
        self.state = PatchState::Growing {
            angles,
            segments,
            vertices,
            cached_matches: Vec::new(),
        };
        self.recompute_matches();
    }

    fn compute_seed_matches(tileset: &Arc<TileSet<T>>, seed_tile_id: usize) -> Vec<PatchMatch> {
        let seed = tileset.rat(seed_tile_id);
        let patch_ts = Arc::new(TileSet::new(vec![seed.clone()]));
        let mf = MatchFinder::crossing(Arc::clone(&patch_ts), Arc::clone(tileset));

        let mut matches = Vec::new();
        for j in 0..mf.num_tiles_b() {
            for mt in mf.valid_matches(0, j) {
                matches.push(PatchMatch {
                    start_a: mt.start_a,
                    len: mt.len,
                    start_b: mt.start_b,
                    tile_id: j,
                });
            }
        }
        matches
    }

    fn recompute_matches(&mut self) {
        let (angles, cached_matches) = match &mut self.state {
            PatchState::Growing {
                angles,
                cached_matches,
                ..
            } => (angles.clone(), cached_matches),
            _ => return,
        };

        cached_matches.clear();
        if angles.is_empty() {
            return;
        }

        let patch_ts = Arc::new(TileSet::new(vec![Rat::from_slice_unchecked(&angles)]));
        let mf = MatchFinder::crossing(Arc::clone(&patch_ts), Arc::clone(&self.tileset));

        for j in 0..mf.num_tiles_b() {
            for mt in mf.valid_matches(0, j) {
                cached_matches.push(PatchMatch {
                    start_a: mt.start_a,
                    len: mt.len,
                    start_b: mt.start_b,
                    tile_id: j,
                });
            }
        }
    }
}

fn compute_glue_angles<T: IsComplex + IsRingOrField + Units>(
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

struct SubMatchSpan {
    seg_idx: usize,
    first_pos: usize,
    last_pos: usize,
    signed_id: Option<i32>,
}

#[allow(clippy::too_many_arguments)]
fn apply_match_locally<T: IsComplex + IsRingOrField + Units>(
    old_segments: &[BoundarySegment],
    old_vertices: &[VertexInfo],
    pm: &PatchMatch,
    tileset: &Arc<TileSet<T>>,
    match_type_index: &Arc<MatchTypeIndex<T>>,
    signed_id: Option<i32>,
) -> Option<(Vec<BoundarySegment>, Vec<VertexInfo>, Vec<VertexInfo>)> {
    let n: usize = old_segments.iter().map(|s| s.tile_len).sum();
    let mlen = pm.len;
    let m = tileset.rat(pm.tile_id).seq().len();

    let offsets = compute_seg_offsets(old_segments);
    let spans = compute_sub_match_spans(old_segments, &offsets, pm, tileset, match_type_index, n);

    let cw_pos = pm.start_a;
    let ccw_pos = (pm.start_a + mlen) % n;

    let seg_len_old = n - mlen;
    let seg_len_new = m - mlen;

    // Map each boundary position in match range to its span index
    let mut pos_to_span: Vec<Option<usize>> = vec![None; n];
    for (si, span) in spans.iter().enumerate() {
        let mut p = span.first_pos;
        loop {
            pos_to_span[p] = Some(si);
            if p == span.last_pos {
                break;
            }
            p = (p + 1) % n;
        }
    }

    // Classify old vertices
    let mut interior_vtx: Vec<usize> = Vec::new();
    for (vi, &vp) in offsets.iter().enumerate() {
        if vp == cw_pos || vp == ccw_pos {
            continue;
        }
        if in_cyclic_range_strict(vp, pm.start_a, mlen, n) {
            interior_vtx.push(vi);
        }
    }

    // Process closed vertices
    let mut closed_vertices: Vec<VertexInfo> = Vec::new();
    for &vi in &interior_vtx {
        let vp = offsets[vi];
        let mut matches = Vec::new();
        if vi < old_vertices.len() {
            matches.extend(old_vertices[vi].matches.iter().copied());
        }
        // CW sub-match: span covering the edge just before the vertex
        let cw_edge = (vp + n - 1) % n;
        if let Some(si) = pos_to_span[cw_edge] {
            if let Some(sid) = spans[si].signed_id {
                matches.push(sid);
            }
        }
        // CCW sub-match: span covering the edge at the vertex
        if let Some(si) = pos_to_span[vp] {
            if let Some(sid) = spans[si].signed_id {
                matches.push(sid);
            }
        }
        canonicalize_closed(&mut matches);
        closed_vertices.push(VertexInfo {
            matches,
            is_open: false,
        });
    }

    // Track which old vertices are NOT consumed (still on boundary, outside match range)
    let consumed: Vec<bool> = (0..old_vertices.len())
        .map(|vi| {
            let vp = offsets[vi];
            vp == cw_pos || vp == ccw_pos || in_cyclic_range_strict(vp, pm.start_a, mlen, n)
        })
        .collect();

    // Build a mapping from old boundary position to old vertex index (for surviving vertices)
    // A vertex survives if it's at a position in the unmatched old portion AND it's not consumed
    // Also: CW junction at cw_pos, CCW junction at ccw_pos are special

    let cw_vtx_idx = offsets.iter().position(|&p| p == cw_pos);
    let ccw_vtx_idx = offsets.iter().position(|&p| p == ccw_pos);

    // Walk the new boundary and build segments + vertices
    // New boundary: unmatched old (ccw_pos to cw_pos, exclusive of match) then new tile unmatched
    let mut new_segments: Vec<BoundarySegment> = Vec::new();
    let mut new_vertices: Vec<VertexInfo> = Vec::new();

    if seg_len_old > 0 {
        let mut p = ccw_pos;
        let mut last_old_seg_idx: Option<usize> = None;
        for _ in 0..seg_len_old {
            let seg_idx = find_segment_at_offsets(&offsets, p);
            let seg = &old_segments[seg_idx];
            let offset_in_seg = (p + n - offsets[seg_idx]) % n;
            let tile_pos = (seg.tile_start + offset_in_seg) % tileset.rat(seg.tile_id).len();

            let same_old_seg = last_old_seg_idx == Some(seg_idx);
            let has_surviving_vertex = offsets
                .iter()
                .enumerate()
                .any(|(vi, &off)| off == p && !consumed[vi] && vi < old_vertices.len());

            if same_old_seg && !has_surviving_vertex {
                new_segments.last_mut().unwrap().tile_len += 1;
            } else {
                if !new_segments.is_empty() {
                    emit_boundary_vertex(
                        p,
                        &offsets,
                        old_vertices,
                        cw_pos,
                        ccw_pos,
                        signed_id,
                        cw_vtx_idx,
                        ccw_vtx_idx,
                        &consumed,
                        &mut new_vertices,
                        n,
                    );
                }
                new_segments.push(BoundarySegment {
                    tile_id: seg.tile_id,
                    tile_start: tile_pos,
                    tile_len: 1,
                });
                last_old_seg_idx = Some(seg_idx);
            }
            p = (p + 1) % n;
        }
    }

    // Transition from old unmatched to new tile unmatched (or wrap)
    // Never fuse here: the CW junction vertex must always be emitted
    if seg_len_new > 0 {
        let new_tile_start = (pm.start_b + mlen) % m;
        if !new_segments.is_empty() {
            emit_boundary_vertex(
                cw_pos,
                &offsets,
                old_vertices,
                cw_pos,
                ccw_pos,
                signed_id,
                cw_vtx_idx,
                ccw_vtx_idx,
                &consumed,
                &mut new_vertices,
                n,
            );
        }
        new_segments.push(BoundarySegment {
            tile_id: pm.tile_id,
            tile_start: new_tile_start,
            tile_len: seg_len_new,
        });
    }

    // Close the cyclic boundary: vertex between last segment and first
    if !new_segments.is_empty() {
        // This vertex is at the CCW junction (start of new boundary = ccw_pos on old)
        let ccw_junction_pos = ccw_pos;
        emit_boundary_vertex(
            ccw_junction_pos,
            &offsets,
            old_vertices,
            cw_pos,
            ccw_pos,
            signed_id,
            cw_vtx_idx,
            ccw_vtx_idx,
            &consumed,
            &mut new_vertices,
            n,
        );
    }

    if new_segments.is_empty() {
        return None;
    }

    debug_assert_eq!(
        new_segments.len(),
        new_vertices.len(),
        "cyclic invariant: segs={} verts={}",
        new_segments.len(),
        new_vertices.len()
    );

    Some((new_segments, new_vertices, closed_vertices))
}

#[allow(clippy::too_many_arguments)]
fn emit_boundary_vertex(
    old_boundary_pos: usize,
    offsets: &[usize],
    old_vertices: &[VertexInfo],
    cw_pos: usize,
    ccw_pos: usize,
    signed_id: Option<i32>,
    cw_vtx_idx: Option<usize>,
    ccw_vtx_idx: Option<usize>,
    consumed: &[bool],
    new_vertices: &mut Vec<VertexInfo>,
    _n: usize,
) {
    let mut matches = Vec::new();

    if old_boundary_pos == cw_pos {
        if let Some(vi) = cw_vtx_idx {
            if vi < old_vertices.len() {
                matches.extend(old_vertices[vi].matches.iter().copied());
            }
        }
        if let Some(sid) = signed_id {
            let inv = -sid;
            if !matches.contains(&inv) {
                matches.push(inv);
            }
        }
        new_vertices.push(VertexInfo {
            matches,
            is_open: true,
        });
        return;
    }

    if old_boundary_pos == ccw_pos {
        if let Some(vi) = ccw_vtx_idx {
            if vi < old_vertices.len() {
                matches.extend(old_vertices[vi].matches.iter().copied());
            }
        }
        if let Some(sid) = signed_id {
            if !matches.contains(&sid) {
                matches.push(sid);
            }
        }
        new_vertices.push(VertexInfo {
            matches,
            is_open: true,
        });
        return;
    }

    // Look for a surviving old vertex at this position
    for (vi, &off) in offsets.iter().enumerate() {
        if off == old_boundary_pos && !consumed[vi] && vi < old_vertices.len() {
            new_vertices.push(old_vertices[vi].clone());
            return;
        }
    }

    // New vertex (segment split): compute from the match at this position
    // This is a vertex inside a segment that got split by a junction endpoint
    // It gets the signed_id from the match type
    if let Some(sid) = signed_id {
        matches.push(sid);
    }
    new_vertices.push(VertexInfo {
        matches,
        is_open: true,
    });
}

#[allow(clippy::too_many_arguments)]
fn compute_sub_match_spans<T: IsComplex + IsRingOrField + Units>(
    old_segments: &[BoundarySegment],
    offsets: &[usize],
    pm: &PatchMatch,
    tileset: &Arc<TileSet<T>>,
    match_type_index: &Arc<MatchTypeIndex<T>>,
    n: usize,
) -> Vec<SubMatchSpan> {
    let mlen = pm.len;
    let mut spans: Vec<SubMatchSpan> = Vec::new();
    let mut p = pm.start_a;

    for _ in 0..mlen {
        let seg_idx = find_segment_at_offsets(offsets, p);
        let continued = spans
            .last()
            .is_some_and(|s| s.seg_idx == seg_idx && (s.last_pos + 1) % n == p);

        if continued {
            spans.last_mut().unwrap().last_pos = p;
        } else {
            spans.push(SubMatchSpan {
                seg_idx,
                first_pos: p,
                last_pos: p,
                signed_id: None,
            });
        }
        p = (p + 1) % n;
    }

    for span in &mut spans {
        let seg = &old_segments[span.seg_idx];
        let offset_in_seg = (span.first_pos + n - offsets[span.seg_idx]) % n;
        let orig_pos = (seg.tile_start + offset_in_seg) % tileset.rat(seg.tile_id).len();
        let new_offset = (span.first_pos + n - pm.start_a) % n;
        let new_pos = (pm.start_b + new_offset) % tileset.rat(pm.tile_id).len();

        let seed_a = tileset.rat(seg.tile_id);
        let seed_b = tileset.rat(pm.tile_id);
        let (ns, slen, ne) = seed_a.get_match((orig_pos as i64, new_pos as i64), seed_b);

        if slen > 0 {
            span.signed_id = match_type_index.signed_id(&MatchType {
                tile_a: seg.tile_id,
                start_a: ns as usize,
                tile_b: pm.tile_id,
                start_b: ne as usize,
                len: slen,
            });
        }
    }

    spans
}

fn in_cyclic_range_strict(pos: usize, start: usize, len: usize, n: usize) -> bool {
    if len <= 1 {
        return false;
    }
    (1..len).any(|i| (start + i) % n == pos)
}

fn find_segment_at_offsets(offsets: &[usize], pos: usize) -> usize {
    for i in (0..offsets.len()).rev() {
        if offsets[i] <= pos {
            return i;
        }
    }
    0
}

fn compute_seg_offsets(segments: &[BoundarySegment]) -> Vec<usize> {
    let mut offsets = Vec::with_capacity(segments.len());
    let mut off = 0usize;
    for seg in segments {
        offsets.push(off);
        off += seg.tile_len;
    }
    offsets
}

fn resolve_signed_id<T: IsComplex + IsRingOrField + Units>(
    segments: &[BoundarySegment],
    angles: &[i8],
    pm: &PatchMatch,
    tileset: &Arc<TileSet<T>>,
    match_type_index: &Arc<MatchTypeIndex<T>>,
) -> Option<i32> {
    if segments.is_empty() {
        return None;
    }
    let n = if angles.is_empty() {
        segments.iter().map(|s| s.tile_len).sum()
    } else {
        angles.len()
    };
    let offsets = compute_seg_offsets(segments);
    let seg_idx = find_segment_at_offsets(&offsets, pm.start_a);
    let seg = &segments[seg_idx];
    let offset_in_seg = (pm.start_a + n - offsets[seg_idx]) % n;
    let orig_pos = (seg.tile_start + offset_in_seg) % tileset.rat(seg.tile_id).len();

    let seed_a = tileset.rat(seg.tile_id);
    let seed_b = tileset.rat(pm.tile_id);
    let (ns, seed_len, ne) = seed_a.get_match((orig_pos as i64, pm.start_b as i64), seed_b);

    if seed_len == 0 {
        return None;
    }

    match_type_index.signed_id(&MatchType {
        tile_a: seg.tile_id,
        start_a: ns as usize,
        tile_b: pm.tile_id,
        start_b: ne as usize,
        len: seed_len,
    })
}

fn canonicalize_closed(matches: &mut [i32]) {
    let rot = lex_min_rotation(matches);
    if rot > 0 {
        matches.rotate_left(rot);
    }
}

fn lex_min_rotation(seq: &[i32]) -> usize {
    if seq.len() <= 1 {
        return 0;
    }
    let n = seq.len();
    let mut best = 0;
    for i in 1..n {
        for j in 0..n {
            let a = seq[(best + j) % n];
            let b = seq[(i + j) % n];
            if b < a {
                best = i;
                break;
            }
            if b > a {
                break;
            }
        }
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic::{ZZ12, ZZ4};
    use crate::intgeom::snake::Snake;
    use crate::intgeom::tiles;

    fn hex_patch() -> GrowingPatch<ZZ12> {
        let hex: Rat<ZZ12> = Rat::from_unchecked(&tiles::hexagon());
        let ts = Arc::new(TileSet::new(vec![hex]));
        let mti = Arc::new(MatchTypeIndex::new(Arc::clone(&ts)));
        GrowingPatch::new(ts, mti, 0)
    }

    fn mixed_patch() -> GrowingPatch<ZZ12> {
        let hex: Rat<ZZ12> = Rat::from_unchecked(&tiles::hexagon());
        let sq: Rat<ZZ12> = Rat::from_unchecked(&tiles::square());
        let ts = Arc::new(TileSet::new(vec![hex, sq]));
        let mti = Arc::new(MatchTypeIndex::new(Arc::clone(&ts)));
        GrowingPatch::new(ts, mti, 0)
    }

    #[test]
    fn seed_patch_has_no_boundary() {
        let gp = hex_patch();
        assert!(!gp.is_growing());
        assert!(!gp.is_closed());
        assert_eq!(gp.boundary_len(), 0);
        assert_eq!(gp.segments().len(), 0);
        assert_eq!(gp.vertices().len(), 0);
    }

    #[test]
    fn seed_patch_has_matches() {
        let gp = hex_patch();
        let matches = gp.get_all_matches();
        assert!(!matches.is_empty(), "seed hex should have matches");
        for m in matches {
            assert_eq!(m.tile_id, 0, "self-match only");
        }
    }

    #[test]
    fn first_add_initializes_patch() {
        let mut gp = hex_patch();
        let matches = gp.get_all_matches().to_vec();
        assert!(!matches.is_empty());

        let pm = &matches[0];
        let _diff = gp.add_tile(pm).expect("first add_tile should succeed");

        assert!(gp.is_growing());
        assert_eq!(gp.boundary_len(), 10);
        assert_eq!(gp.segments().len(), gp.vertices().len(), "cyclic invariant");

        let seg_sum: usize = gp.segments().iter().map(|s| s.tile_len).sum();
        assert_eq!(seg_sum, gp.boundary_len());

        for vtx in gp.vertices() {
            assert!(vtx.is_open);
            assert!(!vtx.matches.is_empty(), "vertex should have match IDs");
        }
    }

    #[test]
    fn first_add_to_rat_matches_glue() {
        let hex: Rat<ZZ12> = Rat::from_unchecked(&tiles::hexagon());
        let ts = Arc::new(TileSet::new(vec![hex.clone()]));
        let mti = Arc::new(MatchTypeIndex::new(Arc::clone(&ts)));
        let mut gp = GrowingPatch::new(ts, mti, 0);

        let matches = gp.get_all_matches().to_vec();
        let pm = &matches[0];
        gp.add_tile(pm).expect("add_tile");

        let rat = gp.to_rat();
        let expected = hex.glue((pm.start_a as i64, pm.start_b as i64), &hex);
        assert_eq!(rat, expected);
    }

    fn verify_bi_tile_patch<T: IsComplex + IsRingOrField + Units>(
        gp: &GrowingPatch<T>,
        seed_len: usize,
        label: &str,
    ) {
        let expected_boundary = 2 * seed_len - 2;
        assert_eq!(gp.boundary_len(), expected_boundary, "{label}: boundary");

        let seg_sum: usize = gp.segments().iter().map(|s| s.tile_len).sum();
        assert_eq!(seg_sum, gp.boundary_len(), "{label}: seg sum");

        assert_eq!(
            gp.segments().len(),
            gp.vertices().len(),
            "{label}: segments == vertices"
        );

        let last = gp.vertices().len() - 1;
        assert!(
            !gp.vertices()[0].matches.is_empty(),
            "{label}: first vertex should have match ID"
        );
        assert!(
            !gp.vertices()[last].matches.is_empty(),
            "{label}: last vertex should have match ID"
        );

        for (i, vtx) in gp.vertices().iter().enumerate() {
            assert!(vtx.is_open, "{label}: vertex {i} should be open");
            for &sid in &vtx.matches {
                assert_ne!(sid, 0, "{label}: vertex {i} signed id nonzero");
            }
        }

        let rat = gp.to_rat();
        assert!(
            Snake::<T>::try_from(rat.seq()).is_ok(),
            "{label}: valid snake"
        );
    }

    #[test]
    fn square_all_16_matches_produce_valid_bi_squares() {
        let sq: Rat<ZZ4> = Rat::from_unchecked(&tiles::square());
        let ts = Arc::new(TileSet::new(vec![sq]));
        let mti = Arc::new(MatchTypeIndex::new(Arc::clone(&ts)));

        let mut seen_signed_ids: std::collections::HashSet<i32> = std::collections::HashSet::new();
        for (mi, pm) in GrowingPatch::new(Arc::clone(&ts), Arc::clone(&mti), 0)
            .get_all_matches()
            .iter()
            .enumerate()
        {
            let label = format!("sq-match-{mi}");
            let mut gp = GrowingPatch::new(Arc::clone(&ts), Arc::clone(&mti), 0);
            gp.add_tile(pm).expect("add_tile");
            verify_bi_tile_patch(&gp, 4, &label);

            let last = gp.vertices().len() - 1;
            seen_signed_ids.insert(gp.vertices()[0].matches[0]);
            seen_signed_ids.insert(gp.vertices()[last].matches[0]);
        }

        assert_eq!(
            seen_signed_ids.len(),
            20,
            "16 raw matches should produce 20 distinct signed IDs"
        );
    }

    #[test]
    fn hexagon_all_36_matches_produce_valid_bi_hexes() {
        let hex: Rat<ZZ12> = Rat::from_unchecked(&tiles::hexagon());
        let ts = Arc::new(TileSet::new(vec![hex]));
        let mti = Arc::new(MatchTypeIndex::new(Arc::clone(&ts)));

        let mut seen_signed_ids: std::collections::HashSet<i32> = std::collections::HashSet::new();
        for (mi, pm) in GrowingPatch::new(Arc::clone(&ts), Arc::clone(&mti), 0)
            .get_all_matches()
            .iter()
            .enumerate()
        {
            let label = format!("hex-match-{mi}");
            let mut gp = GrowingPatch::new(Arc::clone(&ts), Arc::clone(&mti), 0);
            gp.add_tile(pm).expect("add_tile");
            verify_bi_tile_patch(&gp, 6, &label);

            let last = gp.vertices().len() - 1;
            seen_signed_ids.insert(gp.vertices()[0].matches[0]);
            seen_signed_ids.insert(gp.vertices()[last].matches[0]);
        }

        assert_eq!(
            seen_signed_ids.len(),
            42,
            "36 raw matches should produce 42 distinct signed IDs"
        );
    }

    #[test]
    fn to_rat_matches_direct_glue_for_all_matches() {
        let hex: Rat<ZZ12> = Rat::from_unchecked(&tiles::hexagon());
        let ts = Arc::new(TileSet::new(vec![hex.clone()]));
        let mti = Arc::new(MatchTypeIndex::new(Arc::clone(&ts)));

        for pm in GrowingPatch::new(Arc::clone(&ts), Arc::clone(&mti), 0)
            .get_all_matches()
            .iter()
        {
            let mut gp = GrowingPatch::new(Arc::clone(&ts), Arc::clone(&mti), 0);
            gp.add_tile(pm).expect("add_tile");
            let rat = gp.to_rat();
            let expected = hex.glue((pm.start_a as i64, pm.start_b as i64), &hex);
            assert_eq!(rat, expected, "mismatch for {:?}", pm);
        }
    }

    #[test]
    fn segment_cyclic_invariant() {
        let mut gp = hex_patch();
        let matches = gp.get_all_matches().to_vec();
        gp.add_tile(&matches[0]).expect("first add");

        assert_eq!(
            gp.segments().len(),
            gp.vertices().len(),
            "segments == vertices after init"
        );

        for step in 0..5 {
            let matches = gp.get_all_matches().to_vec();
            if matches.is_empty() {
                break;
            }
            gp.add_tile(&matches[0]).expect("add_tile");
            assert_eq!(
                gp.segments().len(),
                gp.vertices().len(),
                "step {}: segments == vertices",
                step + 1
            );
            assert!(
                gp.segments().iter().map(|s| s.tile_len).sum::<usize>() == gp.boundary_len(),
                "step {}: seg sum == boundary",
                step + 1
            );
        }
    }

    #[test]
    fn mixed_hex_square_add_tile() {
        let hex: Rat<ZZ12> = Rat::from_unchecked(&tiles::hexagon());
        let sq: Rat<ZZ12> = Rat::from_unchecked(&tiles::square());
        let ts = Arc::new(TileSet::new(vec![hex.clone(), sq.clone()]));
        let mti = Arc::new(MatchTypeIndex::new(Arc::clone(&ts)));

        let gp = GrowingPatch::new(Arc::clone(&ts), Arc::clone(&mti), 0);
        let hex_sq_matches: Vec<PatchMatch> = gp.get_matches_for_tile(1).cloned().collect();
        assert!(
            !hex_sq_matches.is_empty(),
            "hex seed should have matches with square tile"
        );

        for (i, pm) in hex_sq_matches.iter().enumerate() {
            let label = format!("hex+sq-{i}");
            let mut gp = GrowingPatch::new(Arc::clone(&ts), Arc::clone(&mti), 0);
            gp.add_tile(pm).expect("add_tile");

            let seg_sum: usize = gp.segments().iter().map(|s| s.tile_len).sum();
            assert_eq!(seg_sum, gp.boundary_len(), "{label}: seg sum");

            assert_eq!(
                gp.segments().len(),
                gp.vertices().len(),
                "{label}: cyclic invariant"
            );

            let rat = gp.to_rat();
            let expected = hex.glue((pm.start_a as i64, pm.start_b as i64), &sq);
            assert_eq!(rat, expected, "{label}: to_rat matches glue");
            assert!(
                Snake::<ZZ12>::try_from(rat.seq()).is_ok(),
                "{label}: valid snake"
            );
        }
    }

    #[test]
    fn junction_vertex_ids_nonempty_after_each_add() {
        let mut gp = hex_patch();
        let matches = gp.get_all_matches().to_vec();
        gp.add_tile(&matches[0]).expect("first add");

        for step in 0..4 {
            let matches = gp.get_all_matches().to_vec();
            if matches.is_empty() {
                break;
            }
            gp.add_tile(&matches[0]).expect("add_tile");

            assert_eq!(
                gp.segments().len(),
                gp.vertices().len(),
                "step {}: cyclic",
                step + 1
            );

            if !gp.is_growing() {
                break;
            }

            let last = gp.vertices().len() - 1;
            assert!(
                !gp.vertices()[0].matches.is_empty(),
                "step {}: first vertex has IDs",
                step + 1
            );
            assert!(
                !gp.vertices()[last].matches.is_empty(),
                "step {}: last vertex has IDs",
                step + 1
            );
        }
    }

    use std::collections::{BTreeMap, BTreeSet};

    struct PatchCensus<T: IsComplex + IsRingOrField + Units> {
        ts: Arc<TileSet<T>>,
        mti: Arc<MatchTypeIndex<T>>,
        patches: BTreeMap<Vec<i8>, Vec<Vec<PatchMatch>>>,
    }

    impl PatchCensus<ZZ4> {
        fn for_squares(max_tiles: usize) -> Self {
            let sq: Rat<ZZ4> = Rat::from_unchecked(&tiles::square());
            let ts = Arc::new(TileSet::new(vec![sq]));
            let mti = Arc::new(MatchTypeIndex::new(Arc::clone(&ts)));
            let patches = brute_force_zz4(&ts, &mti, max_tiles);
            Self { ts, mti, patches }
        }
    }

    impl PatchCensus<ZZ12> {
        fn for_hexagons(max_tiles: usize) -> Self {
            let hex: Rat<ZZ12> = Rat::from_unchecked(&tiles::hexagon());
            let ts = Arc::new(TileSet::new(vec![hex]));
            let mti = Arc::new(MatchTypeIndex::new(Arc::clone(&ts)));
            let patches = brute_force_zz12(&ts, &mti, max_tiles);
            Self { ts, mti, patches }
        }
    }

    impl<T: IsComplex + IsRingOrField + Units> PatchCensus<T> {
        fn replay(&self, history: &[PatchMatch]) -> Option<GrowingPatch<T>> {
            self.replay_with_diffs(history).map(|(gp, _)| gp)
        }

        fn replay_with_diffs(
            &self,
            history: &[PatchMatch],
        ) -> Option<(GrowingPatch<T>, Vec<AddTileDiff>)> {
            if history.is_empty() {
                return None;
            }
            let mut gp = GrowingPatch::new(Arc::clone(&self.ts), Arc::clone(&self.mti), 0);
            let mut diffs = Vec::new();
            for pm in history {
                let diff = gp.add_tile(pm)?;
                diffs.push(diff);
            }
            Some((gp, diffs))
        }

        fn replay_with_snapshots(
            &self,
            history: &[PatchMatch],
        ) -> Option<(GrowingPatch<T>, Vec<StepSnapshot<T>>)> {
            if history.is_empty() {
                return None;
            }
            let mut gp = GrowingPatch::new(Arc::clone(&self.ts), Arc::clone(&self.mti), 0);
            let mut snapshots = Vec::new();
            for pm in history {
                let old_segs = gp.to_segments();
                let old_verts = gp.to_vertices();
                let old_offsets = compute_seg_offsets(&old_segs);
                let old_n = gp.boundary_len();
                let diff = gp.add_tile(pm)?;
                let new_segs = gp.to_segments();
                let new_verts = gp.to_vertices();
                snapshots.push(StepSnapshot {
                    pm: pm.clone(),
                    old_segs,
                    old_verts,
                    old_offsets,
                    old_n,
                    diff,
                    new_segs,
                    new_verts,
                    _phantom: std::marker::PhantomData,
                });
            }
            Some((gp, snapshots))
        }

        fn verify_all_steps(&self, label: &str) {
            let mut checked_steps = 0usize;
            let mut checked_closed = 0usize;
            for ways in self.patches.values() {
                for history in ways {
                    if history.is_empty() {
                        continue;
                    }
                    if let Some((gp, snapshots)) = self.replay_with_snapshots(history) {
                        verify_patch_invariants(&gp, label);
                        for snap in &snapshots {
                            verify_step(snap, &self.ts, &self.mti, label);
                            checked_steps += 1;
                            checked_closed += snap.diff.closed_vertices.len();
                        }
                    }
                }
            }
            eprintln!(
                "  [{}] verified {} steps, {} closed vertices",
                label, checked_steps, checked_closed
            );
        }

        fn tile_count(&self, canonical: &[i8]) -> usize {
            self.patches
                .get(canonical)
                .map(|ways| ways[0].len() + 1)
                .unwrap_or(0)
        }

        fn collect_all_vertex_types(&self) -> (BTreeSet<Vec<i32>>, BTreeSet<Vec<i32>>) {
            let mut open: BTreeSet<Vec<i32>> = BTreeSet::new();
            let mut closed: BTreeSet<Vec<i32>> = BTreeSet::new();
            for (canonical, ways) in &self.patches {
                let _ = canonical;
                for history in ways {
                    if let Some((gp, diffs)) = self.replay_with_diffs(history) {
                        for v in gp.to_vertices() {
                            open.insert(v.matches.clone());
                        }
                        for diff in &diffs {
                            for v in &diff.closed_vertices {
                                closed.insert(v.matches.clone());
                            }
                        }
                    }
                }
            }
            (open, closed)
        }

        fn collect_vertex_sequences(&self) -> BTreeMap<Vec<i8>, BTreeSet<Vec<Vec<i32>>>> {
            let mut seqs: BTreeMap<Vec<i8>, BTreeSet<Vec<Vec<i32>>>> = BTreeMap::new();
            for (canonical, ways) in &self.patches {
                let entry = seqs.entry(canonical.clone()).or_default();
                for history in ways {
                    if let Some(gp) = self.replay(history) {
                        let vseq: Vec<Vec<i32>> =
                            gp.to_vertices().iter().map(|v| v.matches.clone()).collect();
                        entry.insert(vseq);
                    }
                }
            }
            seqs
        }

        fn report(&self, label: &str) {
            let mut by_tiles: BTreeMap<usize, (usize, usize)> = BTreeMap::new();
            for ways in self.patches.values() {
                let n = ways[0].len() + 1;
                let e = by_tiles.entry(n).or_insert((0, 0));
                e.0 += 1;
                e.1 += ways.len();
            }

            eprintln!("\n=== {} census ===", label);
            eprintln!("Distinct patches by tile count:");
            for (n, (shapes, ways)) in &by_tiles {
                eprintln!("  {} tiles: {} distinct, {} ways", n, shapes, ways);
            }
            eprintln!("  total distinct patches: {}", self.patches.len());

            let (open_vtypes, closed_vtypes) = self.collect_all_vertex_types();
            eprintln!("\nOpen vertex types ({} total):", open_vtypes.len());
            for vt in &open_vtypes {
                eprintln!("  arity={}: {:?}", vt.len(), vt);
            }
            eprintln!("\nClosed vertex types ({} total):", closed_vtypes.len());
            for vt in &closed_vtypes {
                eprintln!("  arity={}: {:?}", vt.len(), vt);
            }

            let vseqs = self.collect_vertex_sequences();
            eprintln!("\nVertex sequences per rat:");
            for (canonical, seq_set) in &vseqs {
                let n = self.tile_count(canonical);
                eprintln!(
                    "  n={} boundary={:?} ({} distinct vertex sequences):",
                    n,
                    canonical,
                    seq_set.len()
                );
                for vseq in seq_set {
                    eprintln!("    {:?}", vseq);
                }
            }
        }
    }

    struct StepSnapshot<T: IsComplex> {
        pm: PatchMatch,
        old_segs: Vec<BoundarySegment>,
        old_verts: Vec<VertexInfo>,
        old_offsets: Vec<usize>,
        old_n: usize,
        diff: AddTileDiff,
        new_segs: Vec<BoundarySegment>,
        new_verts: Vec<VertexInfo>,
        _phantom: std::marker::PhantomData<T>,
    }

    fn verify_patch_invariants<T: IsComplex + IsRingOrField + Units>(
        gp: &GrowingPatch<T>,
        label: &str,
    ) {
        if !gp.is_growing() {
            return;
        }
        assert_eq!(
            gp.to_segments().len(),
            gp.to_vertices().len(),
            "[{}] cyclic invariant: segs={} verts={}",
            label,
            gp.to_segments().len(),
            gp.to_vertices().len()
        );
        assert_eq!(
            gp.to_segments().iter().map(|s| s.tile_len).sum::<usize>(),
            gp.boundary_len(),
            "[{}] seg sum == boundary",
            label
        );
    }

    fn verify_step<T: IsComplex + IsRingOrField + Units>(
        snap: &StepSnapshot<T>,
        ts: &Arc<TileSet<T>>,
        mti: &Arc<MatchTypeIndex<T>>,
        label: &str,
    ) {
        let n = snap.old_n;
        let mlen = snap.pm.len;

        if snap.old_segs.is_empty() {
            verify_first_add(snap, ts, mti, label);
            return;
        }

        let cw_pos = snap.pm.start_a;
        let ccw_pos = (snap.pm.start_a + mlen) % n;

        // Test A: closed vertex count matches geometric prediction
        let expected_closed: Vec<usize> = snap
            .old_offsets
            .iter()
            .enumerate()
            .filter(|(_, &vp)| {
                vp != cw_pos
                    && vp != ccw_pos
                    && in_cyclic_range_strict(vp, snap.pm.start_a, mlen, n)
            })
            .map(|(vi, _)| vi)
            .collect();
        assert_eq!(
            snap.diff.closed_vertices.len(),
            expected_closed.len(),
            "[{}] closed vertex count: got {} expected {}",
            label,
            snap.diff.closed_vertices.len(),
            expected_closed.len()
        );

        // Test B: no closed vertices when match doesn't cross any vertex
        if expected_closed.is_empty() {
            assert!(
                snap.diff.closed_vertices.is_empty(),
                "[{}] no closed vertices expected for match not crossing any vertex",
                label
            );
        }

        // For each closed vertex, verify content
        for (ci, &vi) in expected_closed.iter().enumerate() {
            let closed_vtx = &snap.diff.closed_vertices[ci];
            let vp = snap.old_offsets[vi];

            // Test C: sub-match signed IDs verified independently
            let mut expected_matches = Vec::new();
            if vi < snap.old_verts.len() {
                expected_matches.extend(snap.old_verts[vi].matches.iter().copied());
            }
            let cw_edge = (vp + n - 1) % n;
            if let Some(sid) = independent_sub_match_signed_id(
                &snap.old_segs,
                &snap.old_offsets,
                cw_edge,
                &snap.pm,
                ts,
                mti,
                n,
            ) {
                expected_matches.push(sid);
            }
            if let Some(sid) = independent_sub_match_signed_id(
                &snap.old_segs,
                &snap.old_offsets,
                vp,
                &snap.pm,
                ts,
                mti,
                n,
            ) {
                expected_matches.push(sid);
            }

            // Test D: canonical rotation correctness
            let expected_canonical = lex_min_rotation_clone(&expected_matches);

            assert_eq!(
                closed_vtx.matches, expected_canonical,
                "[{}] closed vertex {}: got {:?} expected {:?} (raw {:?})",
                label, ci, closed_vtx.matches, expected_canonical, expected_matches
            );
            assert!(
                !closed_vtx.is_open,
                "[{}] closed vertex should have is_open=false",
                label
            );
        }

        // Test E: junction vertex updates correct
        let signed_id = resolve_signed_id(&snap.old_segs, &[], &snap.pm, ts, mti);

        if let Some(sid) = signed_id {
            let inv = -sid;
            let cw_vtx_idx = snap.old_offsets.iter().position(|&p| p == cw_pos);
            if let Some(cwvi) = cw_vtx_idx {
                let old_matches_cw = &snap.old_verts[cwvi].matches;
                let already_has_inv = old_matches_cw.contains(&inv);
                let expected_len = snap.old_verts.get(cwvi).map_or(1, |ov| {
                    ov.matches.len() + if already_has_inv { 0 } else { 1 }
                });
                let has_cw = snap
                    .new_verts
                    .iter()
                    .any(|v| v.matches.contains(&inv) && v.matches.len() == expected_len);
                assert!(
                    has_cw,
                    "[{}] CW junction should have -signed_id={}",
                    label, inv
                );
            }

            let ccw_vtx_idx = snap.old_offsets.iter().position(|&p| p == ccw_pos);
            if let Some(ccwvi) = ccw_vtx_idx {
                let old_matches_ccw = &snap.old_verts[ccwvi].matches;
                let already_has_sid = old_matches_ccw.contains(&sid);
                let expected_len = snap.old_verts.get(ccwvi).map_or(1, |ov| {
                    ov.matches.len() + if already_has_sid { 0 } else { 1 }
                });
                let has_ccw = snap
                    .new_verts
                    .iter()
                    .any(|v| v.matches.contains(&sid) && v.matches.len() == expected_len);
                assert!(
                    has_ccw,
                    "[{}] CCW junction should have +signed_id={}",
                    label, sid
                );
            }
        }

        // Test F: surviving vertices unchanged (excluding junction vertices which get updated)
        for (vi, nv) in snap.old_verts.iter().enumerate() {
            let vp = snap.old_offsets[vi];
            let is_junction = vp == cw_pos || vp == ccw_pos;
            let consumed = is_junction || in_cyclic_range_strict(vp, snap.pm.start_a, mlen, n);
            if !consumed {
                let survived = snap
                    .new_verts
                    .iter()
                    .any(|new_v| new_v.matches == nv.matches);
                assert!(
                    survived,
                    "[{}] surviving vertex {} with matches {:?} should appear in new boundary",
                    label, vi, nv.matches
                );
            }
        }
    }

    fn verify_first_add<T: IsComplex + IsRingOrField + Units>(
        snap: &StepSnapshot<T>,
        ts: &Arc<TileSet<T>>,
        mti: &Arc<MatchTypeIndex<T>>,
        label: &str,
    ) {
        // Test G: bi-tile vertices match match type index
        let seed_id = 0;
        let seed_rat = ts.rat(seed_id);
        let (ns, slen, ne) = seed_rat.get_match(
            (snap.pm.start_a as i64, snap.pm.start_b as i64),
            ts.rat(snap.pm.tile_id),
        );
        assert!(slen > 0, "[{}] first add should have valid match", label);
        let signed_id = mti.signed_id(&MatchType {
            tile_a: seed_id,
            start_a: ns as usize,
            tile_b: snap.pm.tile_id,
            start_b: ne as usize,
            len: slen,
        });

        assert!(
            snap.diff.closed_vertices.is_empty(),
            "[{}] first add should have no closed vertices",
            label
        );

        assert_eq!(
            snap.new_segs.len(),
            snap.new_verts.len(),
            "[{}] bi-tile cyclic invariant",
            label
        );

        if snap.new_segs.len() == 2 {
            let expected_cw: Vec<i32> = signed_id.map(|sid| -sid).into_iter().collect();
            let expected_ccw: Vec<i32> = signed_id.into_iter().collect();
            assert!(
                snap.new_verts[0].matches == expected_cw
                    || snap.new_verts[1].matches == expected_cw,
                "[{}] CW junction should be {:?}, got {:?} and {:?}",
                label,
                expected_cw,
                snap.new_verts[0].matches,
                snap.new_verts[1].matches
            );
            assert!(
                snap.new_verts[0].matches == expected_ccw
                    || snap.new_verts[1].matches == expected_ccw,
                "[{}] CCW junction should be {:?}, got {:?} and {:?}",
                label,
                expected_ccw,
                snap.new_verts[0].matches,
                snap.new_verts[1].matches
            );
        }
    }

    fn independent_sub_match_signed_id<T: IsComplex + IsRingOrField + Units>(
        old_segs: &[BoundarySegment],
        old_offsets: &[usize],
        boundary_pos: usize,
        pm: &PatchMatch,
        ts: &Arc<TileSet<T>>,
        mti: &Arc<MatchTypeIndex<T>>,
        n: usize,
    ) -> Option<i32> {
        let seg_idx = find_segment_at_offsets(old_offsets, boundary_pos);
        let seg = &old_segs[seg_idx];
        let offset_in_seg = (boundary_pos + n - old_offsets[seg_idx]) % n;
        let orig_pos = (seg.tile_start + offset_in_seg) % ts.rat(seg.tile_id).len();
        let new_offset = (boundary_pos + n - pm.start_a) % n;
        let new_pos = (pm.start_b + new_offset) % ts.rat(pm.tile_id).len();

        let seed_a = ts.rat(seg.tile_id);
        let seed_b = ts.rat(pm.tile_id);
        let (ns, slen, ne) = seed_a.get_match((orig_pos as i64, new_pos as i64), seed_b);

        if slen == 0 {
            return None;
        }

        mti.signed_id(&MatchType {
            tile_a: seg.tile_id,
            start_a: ns as usize,
            tile_b: pm.tile_id,
            start_b: ne as usize,
            len: slen,
        })
    }

    fn lex_min_rotation_clone(matches: &[i32]) -> Vec<i32> {
        let rot = lex_min_rotation(matches);
        let mut result = matches.to_vec();
        if rot > 0 {
            result.rotate_left(rot);
        }
        result
    }

    fn brute_force_zz4(
        ts: &Arc<TileSet<ZZ4>>,
        mti: &Arc<MatchTypeIndex<ZZ4>>,
        max_tiles: usize,
    ) -> BTreeMap<Vec<i8>, Vec<Vec<PatchMatch>>> {
        let mut results: BTreeMap<Vec<i8>, Vec<Vec<PatchMatch>>> = BTreeMap::new();

        fn recurse(
            gp: &mut GrowingPatch<ZZ4>,
            history: &mut Vec<PatchMatch>,
            max_tiles: usize,
            results: &mut BTreeMap<Vec<i8>, Vec<Vec<PatchMatch>>>,
        ) {
            let num_tiles = history.len() + 1;
            let canonical = gp.to_rat().seq().to_vec();
            results.entry(canonical).or_default().push(history.clone());

            if num_tiles >= max_tiles || !gp.is_growing() {
                return;
            }

            let matches = gp.get_all_matches().to_vec();
            for pm in &matches {
                let mut gp2 = gp.clone();
                if gp2.add_tile(pm).is_some() {
                    history.push(pm.clone());
                    recurse(&mut gp2, history, max_tiles, results);
                    history.pop();
                }
            }
        }

        let seed_seq = ts.rat(0).seq().to_vec();
        results.entry(seed_seq).or_default().push(Vec::new());

        let seed_matches = GrowingPatch::new(Arc::clone(ts), Arc::clone(mti), 0)
            .get_all_matches()
            .to_vec();
        for pm in &seed_matches {
            let mut gp = GrowingPatch::new(Arc::clone(ts), Arc::clone(mti), 0);
            gp.add_tile(pm).expect("first add");
            let mut history = vec![pm.clone()];
            recurse(&mut gp, &mut history, max_tiles, &mut results);
        }

        results
    }

    fn brute_force_zz12(
        ts: &Arc<TileSet<ZZ12>>,
        mti: &Arc<MatchTypeIndex<ZZ12>>,
        max_tiles: usize,
    ) -> BTreeMap<Vec<i8>, Vec<Vec<PatchMatch>>> {
        let mut results: BTreeMap<Vec<i8>, Vec<Vec<PatchMatch>>> = BTreeMap::new();

        fn recurse(
            gp: &mut GrowingPatch<ZZ12>,
            history: &mut Vec<PatchMatch>,
            max_tiles: usize,
            results: &mut BTreeMap<Vec<i8>, Vec<Vec<PatchMatch>>>,
        ) {
            let num_tiles = history.len() + 1;
            if gp.is_growing() {
                let canonical = gp.to_rat().seq().to_vec();
                results.entry(canonical).or_default().push(history.clone());
            }

            if num_tiles >= max_tiles || !gp.is_growing() {
                return;
            }

            let matches = gp.get_all_matches().to_vec();
            for pm in &matches {
                let mut gp2 = gp.clone();
                if gp2.add_tile(pm).is_some() {
                    history.push(pm.clone());
                    recurse(&mut gp2, history, max_tiles, results);
                    history.pop();
                }
            }
        }

        let seed_seq = ts.rat(0).seq().to_vec();
        results.entry(seed_seq).or_default().push(Vec::new());

        let seed_matches = GrowingPatch::new(Arc::clone(ts), Arc::clone(mti), 0)
            .get_all_matches()
            .to_vec();
        for pm in &seed_matches {
            let mut gp = GrowingPatch::new(Arc::clone(ts), Arc::clone(mti), 0);
            gp.add_tile(pm).expect("first add");
            let mut history = vec![pm.clone()];
            recurse(&mut gp, &mut history, max_tiles, &mut results);
        }

        results
    }

    #[test]
    fn square_brute_force_up_to_4_tiles() {
        let census = PatchCensus::for_squares(4);
        census.report("square");
        census.verify_all_steps("square");

        let mut by_tiles: BTreeMap<usize, (usize, usize)> = BTreeMap::new();
        for ways in census.patches.values() {
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

        for ways in census.patches.values() {
            for history in ways {
                if history.is_empty() {
                    continue;
                }
                let mut gp = GrowingPatch::new(Arc::clone(&census.ts), Arc::clone(&census.mti), 0);
                for pm in history {
                    assert!(gp.add_tile(pm).is_some(), "replay should succeed");
                }
                assert_eq!(
                    gp.to_segments().len(),
                    gp.to_vertices().len(),
                    "cyclic invariant"
                );
                assert_eq!(
                    gp.to_segments().iter().map(|s| s.tile_len).sum::<usize>(),
                    gp.boundary_len(),
                    "seg sum == boundary"
                );
                if gp.is_growing() {
                    assert!(
                        Snake::<ZZ4>::try_from(gp.to_rat().seq()).is_ok(),
                        "valid snake"
                    );
                }
            }
        }
    }

    #[test]
    fn hex_brute_force_up_to_3_tiles() {
        let census = PatchCensus::for_hexagons(3);
        census.report("hexagon");
        census.verify_all_steps("hexagon");

        let mut by_tiles: BTreeMap<usize, usize> = BTreeMap::new();
        for ways in census.patches.values() {
            let n = ways[0].len() + 1;
            by_tiles.entry(n).and_modify(|c| *c += 1).or_insert(1);
        }

        assert_eq!(by_tiles.get(&1).copied().unwrap_or(0), 1, "1 mono-hex");
        assert_eq!(by_tiles.get(&2).copied().unwrap_or(0), 1, "1 bi-hex");
        assert!(
            by_tiles.get(&3).copied().unwrap_or(0) >= 1,
            "at least 1 tri-hex"
        );

        for ways in census.patches.values() {
            for history in ways {
                if history.is_empty() {
                    continue;
                }
                let mut gp = GrowingPatch::new(Arc::clone(&census.ts), Arc::clone(&census.mti), 0);
                for pm in history {
                    gp.add_tile(pm).expect("add_tile");
                }
                assert!(gp.is_growing(), "patch should be growing");
                assert_eq!(
                    gp.to_segments().len(),
                    gp.to_vertices().len(),
                    "cyclic invariant"
                );
                assert!(
                    Snake::<ZZ12>::try_from(gp.to_rat().seq()).is_ok(),
                    "valid snake"
                );
            }
        }
    }
}
