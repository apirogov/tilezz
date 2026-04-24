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
        let seed_seq = seed_rat.seq();
        let n = seed_seq.len();
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

        let seed_angles = seed_seq.to_vec();
        let fake_old_angles = seed_angles.clone();
        let fake_seg = BoundarySegment {
            tile_id: seed_id,
            tile_start: 0,
            tile_len: n,
        };

        let new_angles = compute_glue_angles::<T>(&fake_old_angles, pm, &self.tileset)?;

        let seg_len_old = n - mlen;
        let seg_len_new = m - mlen;

        let (segments, vertices, closed) = rebuild_by_walk(
            &[fake_seg],
            &[],
            pm,
            &fake_old_angles,
            &self.tileset,
            signed_id,
            seg_len_old,
            seg_len_new,
        );

        if segments.is_empty() {
            return None;
        }

        self.state = PatchState::Growing {
            angles: new_angles,
            segments,
            vertices,
            cached_matches: Vec::new(),
        };
        self.recompute_matches();

        let _ = closed;
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

        let signed_id = resolve_signed_id_from_segments(
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

        let seg_len_old = n - mlen;
        let seg_len_new = m - mlen;

        let (new_segments, new_vertices, closed_vertices) = rebuild_by_walk(
            &segments,
            &vertices,
            pm,
            &angles,
            &self.tileset,
            signed_id,
            seg_len_old,
            seg_len_new,
        );

        if new_segments.is_empty() {
            self.restore_growing(angles, segments, vertices);
            return None;
        }

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

#[allow(clippy::too_many_arguments)]
fn rebuild_by_walk<T: IsComplex + IsRingOrField + Units>(
    old_segments: &[BoundarySegment],
    old_vertices: &[VertexInfo],
    pm: &PatchMatch,
    old_angles: &[i8],
    tileset: &Arc<TileSet<T>>,
    signed_id: Option<i32>,
    seg_len_old: usize,
    seg_len_new: usize,
) -> (Vec<BoundarySegment>, Vec<VertexInfo>, Vec<VertexInfo>) {
    let n = old_angles.len();
    let mlen = pm.len;
    let m = tileset.rat(pm.tile_id).seq().len();
    let match_end = (pm.start_a + mlen) % n;
    let new_boundary_len = seg_len_old + seg_len_new;

    let mut old_offsets = Vec::with_capacity(old_segments.len());
    let mut off = 0usize;
    for seg in old_segments {
        old_offsets.push(off);
        off += seg.tile_len;
    }

    let mut segments: Vec<BoundarySegment> = Vec::new();
    let mut inherited: Vec<Option<usize>> = Vec::new();

    let mut cur_tile_id = usize::MAX;
    let mut cur_tile_start = 0usize;
    let mut cur_len = 0usize;
    let mut cur_old_vertex: Option<usize> = None;

    for p in 0..new_boundary_len {
        let (tile_id, tile_pos_in_tile) = if p < seg_len_old {
            let old_pos = (match_end + p) % n;
            let mut found = None;
            for (i, seg) in old_segments.iter().enumerate() {
                if old_offsets[i] <= old_pos && old_pos < old_offsets[i] + seg.tile_len {
                    let offset_in_seg = old_pos - old_offsets[i];
                    let tp = (seg.tile_start + offset_in_seg) % tileset.rat(seg.tile_id).len();
                    found = Some((seg.tile_id, tp));
                    break;
                }
            }
            match found {
                Some(v) => v,
                None => return (Vec::new(), Vec::new(), Vec::new()),
            }
        } else {
            let tp = (pm.start_b + (p - seg_len_old)) % m;
            (pm.tile_id, tp)
        };

        let tile_len_total = tileset.rat(tile_id).len();
        let contiguous = cur_len > 0
            && tile_id == cur_tile_id
            && tile_pos_in_tile == (cur_tile_start + cur_len) % tile_len_total;

        if !contiguous && cur_len > 0 {
            segments.push(BoundarySegment {
                tile_id: cur_tile_id,
                tile_start: cur_tile_start,
                tile_len: cur_len,
            });
            inherited.push(cur_old_vertex.take());
            cur_len = 0;
        }

        if cur_len == 0 {
            cur_tile_id = tile_id;
            cur_tile_start = tile_pos_in_tile;
            cur_len = 0;
            if p < seg_len_old {
                let old_pos = (match_end + p) % n;
                cur_old_vertex = old_offsets
                    .iter()
                    .enumerate()
                    .find(|&(_, &off)| off == old_pos)
                    .map(|(i, _)| i);
            } else {
                cur_old_vertex = None;
            }
        }
        cur_len += 1;
    }

    if cur_len > 0 {
        segments.push(BoundarySegment {
            tile_id: cur_tile_id,
            tile_start: cur_tile_start,
            tile_len: cur_len,
        });
        inherited.push(cur_old_vertex);
    }

    let num_old_out = segments
        .iter()
        .enumerate()
        .take_while(|(i, _)| {
            let seg_end: usize = segments[..=*i].iter().map(|s| s.tile_len).sum();
            seg_end <= seg_len_old
        })
        .count();

    if segments.is_empty() {
        return (Vec::new(), Vec::new(), Vec::new());
    }

    let inv_sid = signed_id.map(|id| -id);
    let num_segs = segments.len();

    let cw_junction_vtx = 0;
    let ccw_junction_vtx = if num_old_out < num_segs {
        num_old_out
    } else {
        0
    };

    let mut open_vertices = Vec::with_capacity(num_segs);
    let closed_vertices = Vec::new();

    for i in 0..num_segs {
        let mut matches = Vec::new();

        if let Some(oi) = inherited.get(i).and_then(|x| *x) {
            if oi < old_vertices.len() {
                matches.extend(old_vertices[oi].matches.iter().copied());
            }
        }

        if i == cw_junction_vtx {
            if let Some(sid) = inv_sid {
                matches.push(sid);
            }
        }
        if ccw_junction_vtx != cw_junction_vtx && i == ccw_junction_vtx {
            if let Some(sid) = signed_id {
                matches.push(sid);
            }
        }
        if ccw_junction_vtx == cw_junction_vtx {
            if let Some(sid) = inv_sid {
                if !matches.contains(&sid) {
                    matches.push(sid);
                }
            }
            if let Some(sid) = signed_id {
                if !matches.contains(&sid) {
                    matches.push(sid);
                }
            }
        }

        let is_open = true;

        open_vertices.push(VertexInfo { matches, is_open });
    }

    (segments, open_vertices, closed_vertices)
}

fn resolve_signed_id_from_segments<T: IsComplex + IsRingOrField + Units>(
    segments: &[BoundarySegment],
    angles: &[i8],
    pm: &PatchMatch,
    tileset: &Arc<TileSet<T>>,
    match_type_index: &Arc<MatchTypeIndex<T>>,
) -> Option<i32> {
    if segments.is_empty() {
        return None;
    }

    let seg_idx = find_segment_at(segments, pm.start_a);
    let seg = &segments[seg_idx];

    let offsets = compute_seg_offsets(segments);
    let offset_in_seg = (pm.start_a + angles.len() - offsets[seg_idx]) % angles.len();
    let orig_pos = (seg.tile_start + offset_in_seg) % tileset.rat(seg.tile_id).len();

    let seed_a = tileset.rat(seg.tile_id);
    let seed_b = tileset.rat(pm.tile_id);
    let (ns, seed_len, ne) = seed_a.get_match((orig_pos as i64, pm.start_b as i64), seed_b);

    if seed_len == 0 {
        return None;
    }

    let orig_mt = MatchType {
        tile_a: seg.tile_id,
        start_a: ns as usize,
        tile_b: pm.tile_id,
        start_b: ne as usize,
        len: seed_len,
    };

    match_type_index.signed_id(&orig_mt)
}

fn find_segment_at(segments: &[BoundarySegment], pos: usize) -> usize {
    let offsets = compute_seg_offsets(segments);
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
}
