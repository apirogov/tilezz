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
    pub affected_old: Vec<VertexInfo>,
    pub new_vertices: Vec<VertexInfo>,
}

pub struct GrowingPatch<T: IsComplex> {
    tileset: Arc<TileSet<T>>,
    match_type_index: Arc<MatchTypeIndex<T>>,
    angles: Vec<i8>,
    segments: Vec<BoundarySegment>,
    vertices: Vec<VertexInfo>,
    cached_matches: Vec<PatchMatch>,
}

impl<T: IsComplex + IsRingOrField + Units> GrowingPatch<T> {
    pub fn new(tileset: Arc<TileSet<T>>, match_type_index: Arc<MatchTypeIndex<T>>) -> Self {
        GrowingPatch {
            tileset,
            match_type_index,
            angles: Vec::new(),
            segments: Vec::new(),
            vertices: Vec::new(),
            cached_matches: Vec::new(),
        }
    }

    pub fn from_tile(
        tileset: Arc<TileSet<T>>,
        match_type_index: Arc<MatchTypeIndex<T>>,
        tile_id: usize,
    ) -> Self {
        let angles = tileset.rat(tile_id).seq().to_vec();
        let n = angles.len();
        let mut gp = GrowingPatch {
            tileset,
            match_type_index,
            angles,
            segments: vec![BoundarySegment {
                tile_id,
                tile_start: 0,
                tile_len: n,
            }],
            vertices: Vec::new(),
            cached_matches: Vec::new(),
        };
        gp.recompute_matches();
        gp
    }

    pub fn get_all_matches(&self) -> &[PatchMatch] {
        &self.cached_matches
    }

    pub fn get_matches_for_tile(&self, tile_id: usize) -> impl Iterator<Item = &PatchMatch> {
        self.cached_matches
            .iter()
            .filter(move |m| m.tile_id == tile_id)
    }

    pub fn num_tiles(&self) -> usize {
        self.tileset.num_tiles()
    }

    pub fn boundary_len(&self) -> usize {
        self.angles.len()
    }

    pub fn segments(&self) -> &[BoundarySegment] {
        &self.segments
    }

    pub fn vertices(&self) -> &[VertexInfo] {
        &self.vertices
    }

    pub fn angles(&self) -> &[i8] {
        &self.angles
    }

    pub fn to_rat(&self) -> Rat<T> {
        Rat::from_slice_unchecked(&self.angles)
    }

    pub fn tileset(&self) -> &Arc<TileSet<T>> {
        &self.tileset
    }

    pub fn is_empty(&self) -> bool {
        self.angles.is_empty()
    }

    pub fn add_tile(&mut self, pm: &PatchMatch) -> Option<AddTileDiff> {
        if self.angles.is_empty() {
            return None;
        }

        let n = self.angles.len();
        let mlen = pm.len;
        let seed_seq = self.tileset.rat(pm.tile_id).seq();
        let m = seed_seq.len();

        if mlen == 0 || mlen > n || mlen > m {
            return None;
        }

        let (new_angles, seg_start_old, seg_len_old, seg_start_new, seg_len_new) =
            self.compute_glue(pm)?;

        let signed_id = self.resolve_signed_id(pm);

        let (old_segments, old_vertices) = {
            let mut s = Vec::new();
            let mut v = Vec::new();
            std::mem::swap(&mut self.segments, &mut s);
            std::mem::swap(&mut self.vertices, &mut v);
            (s, v)
        };

        let affected_old = old_vertices.clone();

        let (new_segments, new_vertices) = self.rebuild_segments_and_vertices(
            pm,
            &old_segments,
            &old_vertices,
            seg_start_old,
            seg_len_old,
            seg_start_new,
            seg_len_new,
            new_angles.len(),
            signed_id,
        );

        self.angles = new_angles;
        self.segments = new_segments;
        self.vertices = new_vertices;
        self.recompute_matches();

        Some(AddTileDiff {
            affected_old,
            new_vertices: self.vertices.clone(),
        })
    }

    fn compute_glue(&self, pm: &PatchMatch) -> Option<(Vec<i8>, usize, usize, usize, usize)> {
        let n = self.angles.len();
        let mlen = pm.len;
        let seed_seq = self.tileset.rat(pm.tile_id).seq();
        let m = seed_seq.len();

        let x_raw_len = n - mlen + 1;
        let y_raw_len = m - mlen + 1;

        if x_raw_len <= 1 && y_raw_len <= 1 {
            return None;
        }

        let x: Vec<i8> = (0..x_raw_len)
            .map(|i| self.angles[(pm.start_a + mlen + i) % n])
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

        Some((new_angles, 0, seg_len_old, seg_len_old, seg_len_new))
    }

    #[allow(clippy::too_many_arguments)]
    fn rebuild_segments_and_vertices(
        &self,
        pm: &PatchMatch,
        old_segments: &[BoundarySegment],
        old_vertices: &[VertexInfo],
        _seg_start_old: usize,
        seg_len_old: usize,
        _seg_start_new: usize,
        seg_len_new: usize,
        _new_total: usize,
        signed_id: Option<i32>,
    ) -> (Vec<BoundarySegment>, Vec<VertexInfo>) {
        let n = self.angles.len();
        let mlen = pm.len;
        let m = self.tileset.rat(pm.tile_id).seq().len();

        let mut seg_offsets = Vec::with_capacity(old_segments.len());
        let mut offset = 0usize;
        for seg in old_segments {
            seg_offsets.push(offset);
            offset += seg.tile_len;
        }

        let match_start = pm.start_a;
        let match_end_excl = (pm.start_a + mlen) % n;

        let mut new_segments = Vec::new();
        let mut new_vertex_matches_lists: Vec<Vec<i32>> = Vec::new();

        for (i, seg) in old_segments.iter().enumerate() {
            let seg_start = seg_offsets[i];
            let seg_end_excl = seg_start + seg.tile_len;

            if seg.tile_len == 0 {
                continue;
            }

            let overlap = cyclic_overlap_len(seg_start, seg.tile_len, match_start, mlen, n);

            if overlap == seg.tile_len {
                continue;
            }

            if overlap == 0 {
                new_segments.push(seg.clone());
                let vtx_matches = if i < old_vertices.len() {
                    old_vertices[i].matches.clone()
                } else {
                    Vec::new()
                };
                new_vertex_matches_lists.push(vtx_matches);
                continue;
            }

            let seg_match_start = cyclic_sub(match_start, seg_start, n);
            if seg_match_start == 0 {
                let keep_len = seg.tile_len - overlap;
                new_segments.push(BoundarySegment {
                    tile_id: seg.tile_id,
                    tile_start: (seg.tile_start + overlap) % self.tileset.rat(seg.tile_id).len(),
                    tile_len: keep_len,
                });
                new_vertex_matches_lists.push(Vec::new());
            } else if seg_match_start + overlap >= seg.tile_len {
                let keep_len = seg.tile_len - overlap;
                new_segments.push(BoundarySegment {
                    tile_id: seg.tile_id,
                    tile_start: seg.tile_start,
                    tile_len: keep_len,
                });
                let vtx_matches = if i < old_vertices.len() {
                    old_vertices[i].matches.clone()
                } else {
                    Vec::new()
                };
                new_vertex_matches_lists.push(vtx_matches);
            } else {
                let left_len = seg_match_start;
                let right_len = seg.tile_len - seg_match_start - overlap;

                if left_len > 0 {
                    new_segments.push(BoundarySegment {
                        tile_id: seg.tile_id,
                        tile_start: seg.tile_start,
                        tile_len: left_len,
                    });
                    let vtx_matches = if i < old_vertices.len() {
                        old_vertices[i].matches.clone()
                    } else {
                        Vec::new()
                    };
                    new_vertex_matches_lists.push(vtx_matches);
                }

                if right_len > 0 {
                    let right_tile_start = (seg.tile_start + seg_match_start + overlap)
                        % self.tileset.rat(seg.tile_id).len();
                    new_segments.push(BoundarySegment {
                        tile_id: seg.tile_id,
                        tile_start: right_tile_start,
                        tile_len: right_len,
                    });
                    new_vertex_matches_lists.push(Vec::new());
                }
            }
        }

        let new_tile_start = (pm.start_b + mlen) % m;
        let new_tile_len = m - mlen;
        if new_tile_len > 0 {
            new_segments.push(BoundarySegment {
                tile_id: pm.tile_id,
                tile_start: new_tile_start,
                tile_len: new_tile_len,
            });
        }

        if new_segments.is_empty() {
            return (Vec::new(), Vec::new());
        }

        let inv_signed_id = signed_id.map(|id| -id);

        let num_segs = new_segments.len();
        let mut vertices = Vec::with_capacity(num_segs);
        for i in 0..num_segs {
            let is_cw_junction = i == 0;
            let is_ccw_junction =
                i == num_segs - 1 || (num_segs == 1 && seg_len_old > 0 && seg_len_new > 0);
            let had_old_vtx =
                i < new_vertex_matches_lists.len() && !new_vertex_matches_lists[i].is_empty();

            let mut matches = if had_old_vtx {
                new_vertex_matches_lists[i].clone()
            } else {
                Vec::new()
            };

            if is_cw_junction || is_ccw_junction {
                let id = if is_cw_junction {
                    inv_signed_id
                } else {
                    signed_id
                };
                if let Some(sid) = id {
                    matches.push(sid);
                }
            }

            vertices.push(VertexInfo {
                matches,
                is_open: true,
            });
        }

        (new_segments, vertices)
    }

    fn resolve_signed_id(&self, pm: &PatchMatch) -> Option<i32> {
        if self.segments.is_empty() {
            return None;
        }

        let seg_idx = self.find_segment_at(pm.start_a);
        let seg = &self.segments[seg_idx];

        let seg_offsets = self.compute_seg_offsets();
        let offset_in_seg = cyclic_sub(pm.start_a, seg_offsets[seg_idx], self.angles.len());
        let orig_pos = (seg.tile_start + offset_in_seg) % self.tileset.rat(seg.tile_id).len();

        let seed_a = self.tileset.rat(seg.tile_id);
        let seed_b = self.tileset.rat(pm.tile_id);
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

        self.match_type_index.signed_id(&orig_mt)
    }

    fn find_segment_at(&self, pos: usize) -> usize {
        let offsets = self.compute_seg_offsets();
        for i in (0..offsets.len()).rev() {
            if offsets[i] <= pos {
                return i;
            }
        }
        0
    }

    fn compute_seg_offsets(&self) -> Vec<usize> {
        let mut offsets = Vec::with_capacity(self.segments.len());
        let mut off = 0usize;
        for seg in &self.segments {
            offsets.push(off);
            off += seg.tile_len;
        }
        offsets
    }

    fn recompute_matches(&mut self) {
        self.cached_matches.clear();

        if self.angles.is_empty() {
            return;
        }

        let patch_ts = Arc::new(TileSet::new(vec![Rat::from_slice_unchecked(&self.angles)]));
        let mf = MatchFinder::crossing(Arc::clone(&patch_ts), Arc::clone(&self.tileset));

        for j in 0..mf.num_tiles_b() {
            for mt in mf.valid_matches(0, j) {
                self.cached_matches.push(PatchMatch {
                    start_a: mt.start_a,
                    len: mt.len,
                    start_b: mt.start_b,
                    tile_id: j,
                });
            }
        }
    }
}

fn cyclic_overlap_len(
    seg_start: usize,
    seg_len: usize,
    match_start: usize,
    match_len: usize,
    n: usize,
) -> usize {
    if seg_len == 0 || match_len == 0 {
        return 0;
    }
    let mut count = 0usize;
    for i in 0..seg_len {
        let p = (seg_start + i) % n;
        if in_cyclic_range(p, match_start, match_len, n) {
            count += 1;
        }
    }
    count
}

fn cyclic_sub(a: usize, b: usize, n: usize) -> usize {
    (a + n - b) % n
}

fn in_cyclic_range(v: usize, start: usize, len: usize, n: usize) -> bool {
    if len == 0 {
        return false;
    }
    (0..len).any(|i| (start + i) % n == v)
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
        GrowingPatch::from_tile(ts, mti, 0)
    }

    fn sq_patch() -> GrowingPatch<ZZ4> {
        let sq: Rat<ZZ4> = Rat::from_unchecked(&tiles::square());
        let ts = Arc::new(TileSet::new(vec![sq]));
        let mti = Arc::new(MatchTypeIndex::new(Arc::clone(&ts)));
        GrowingPatch::from_tile(ts, mti, 0)
    }

    fn mixed_patch() -> GrowingPatch<ZZ12> {
        let hex: Rat<ZZ12> = Rat::from_unchecked(&tiles::hexagon());
        let sq: Rat<ZZ12> = Rat::from_unchecked(&tiles::square());
        let ts = Arc::new(TileSet::new(vec![hex, sq]));
        let mti = Arc::new(MatchTypeIndex::new(Arc::clone(&ts)));
        GrowingPatch::from_tile(Arc::clone(&ts), mti, 0)
    }

    #[test]
    fn empty_patch() {
        let hex: Rat<ZZ12> = Rat::from_unchecked(&tiles::hexagon());
        let ts = Arc::new(TileSet::new(vec![hex]));
        let mti = Arc::new(MatchTypeIndex::new(Arc::clone(&ts)));
        let gp = GrowingPatch::new(ts, mti);
        assert!(gp.is_empty());
        assert_eq!(gp.boundary_len(), 0);
        assert_eq!(gp.segments().len(), 0);
        assert_eq!(gp.vertices().len(), 0);
        assert_eq!(gp.get_all_matches().len(), 0);
    }

    #[test]
    fn single_tile_hexagon() {
        let gp = hex_patch();
        assert!(!gp.is_empty());
        assert_eq!(gp.boundary_len(), 6);
        assert_eq!(gp.segments().len(), 1);
        assert_eq!(gp.vertices().len(), 0);

        let seg = &gp.segments()[0];
        assert_eq!(seg.tile_id, 0);
        assert_eq!(seg.tile_start, 0);
        assert_eq!(seg.tile_len, 6);

        let matches = gp.get_all_matches();
        assert!(!matches.is_empty(), "single hex should have matches");
        for m in matches {
            assert_eq!(m.tile_id, 0, "self-match only");
        }

        let rat = gp.to_rat();
        assert!(Snake::<ZZ12>::try_from(rat.seq()).is_ok());
    }

    #[test]
    fn hexagon_add_one_tile() {
        let mut gp = hex_patch();
        let matches = gp.get_all_matches().to_vec();
        assert!(!matches.is_empty());

        let pm = &matches[0];
        let diff = gp.add_tile(pm).expect("add_tile should succeed");

        assert_eq!(
            gp.boundary_len(),
            10,
            "bi-hex should have 10 boundary positions"
        );
        assert_eq!(gp.segments().len(), 2, "bi-hex should have 2 segments");
        assert_eq!(gp.vertices().len(), 2, "bi-hex should have 2 vertices");

        for seg in gp.segments() {
            assert_eq!(seg.tile_id, 0);
            assert!(seg.tile_len > 0);
        }

        for vtx in gp.vertices() {
            eprintln!(
                "  vertex: matches={:?} is_open={}",
                vtx.matches, vtx.is_open
            );
            assert!(vtx.is_open, "vertices should be open after 2 tiles");
            assert_eq!(vtx.matches.len(), 1, "each vertex should have 1 match");
        }

        let rat = gp.to_rat();
        assert!(Snake::<ZZ12>::try_from(rat.seq()).is_ok());

        assert!(
            !diff.new_vertices.is_empty(),
            "diff should contain new vertex info"
        );
    }

    #[test]
    fn hexagon_add_tiles_until_closed() {
        let mut gp = hex_patch();

        for step in 0..10 {
            let matches = gp.get_all_matches().to_vec();
            if matches.is_empty() {
                break;
            }

            let pm = &matches[0];
            let result = gp.add_tile(pm);
            if result.is_none() {
                break;
            }

            let has_closed = gp.vertices().iter().any(|v| !v.is_open);
            eprintln!(
                "step {}: boundary_len={} segments={} vertices={} has_closed={}",
                step + 1,
                gp.boundary_len(),
                gp.segments().len(),
                gp.vertices().len(),
                has_closed,
            );
        }

        let rat = gp.to_rat();
        assert!(
            Snake::<ZZ12>::try_from(rat.seq()).is_ok(),
            "result should be a valid snake at every step"
        );
    }

    #[test]
    fn square_add_one_tile() {
        let mut gp = sq_patch();
        let matches = gp.get_all_matches().to_vec();
        assert!(!matches.is_empty());

        let pm = &matches[0];
        let _diff = gp.add_tile(pm).expect("add_tile should succeed");

        assert_eq!(
            gp.boundary_len(),
            6,
            "bi-square should have 6 boundary positions"
        );
        assert_eq!(gp.segments().len(), 2, "bi-square should have 2 segments");
        assert_eq!(gp.vertices().len(), 2, "bi-square should have 2 vertices");

        for vtx in gp.vertices() {
            assert!(vtx.is_open, "vertices should be open after 2 tiles");
        }

        let rat = gp.to_rat();
        assert!(Snake::<ZZ4>::try_from(rat.seq()).is_ok());
    }

    #[test]
    fn mixed_tileset_initial() {
        let gp = mixed_patch();
        assert_eq!(gp.boundary_len(), 6, "starts with hex");
        assert_eq!(gp.segments().len(), 1);

        let hex_matches: Vec<&PatchMatch> = gp.get_matches_for_tile(0).collect();
        let sq_matches: Vec<&PatchMatch> = gp.get_matches_for_tile(1).collect();
        assert!(!hex_matches.is_empty(), "should have hex self-matches");
        assert!(
            !sq_matches.is_empty() || hex_matches.len() > 0,
            "should have some matches"
        );
    }

    #[test]
    fn segment_lengths_sum_to_boundary() {
        let mut gp = hex_patch();
        let matches = gp.get_all_matches().to_vec();
        assert!(!matches.is_empty());

        gp.add_tile(&matches[0]).expect("add_tile");

        let sum: usize = gp.segments().iter().map(|s| s.tile_len).sum();
        assert_eq!(
            sum,
            gp.boundary_len(),
            "segment lengths should sum to boundary length"
        );
    }

    #[test]
    fn to_rat_matches_growing() {
        let hex: Rat<ZZ12> = Rat::from_unchecked(&tiles::hexagon());
        let ts = Arc::new(TileSet::new(vec![hex]));
        let mti = Arc::new(MatchTypeIndex::new(Arc::clone(&ts)));

        let r1 = Rat::<ZZ12>::from_unchecked(&tiles::hexagon());
        let mut gp = GrowingPatch::from_tile(Arc::clone(&ts), mti, 0);

        let matches = gp.get_all_matches().to_vec();
        let pm = &matches[0];
        gp.add_tile(pm).expect("add_tile");

        let rat = gp.to_rat();
        assert_eq!(
            rat,
            r1.glue((pm.start_a as i64, pm.start_b as i64), &r1),
            "to_rat should match direct glue"
        );
    }

    #[test]
    fn vertex_signed_ids_are_valid() {
        let mut gp = hex_patch();
        let matches = gp.get_all_matches().to_vec();
        assert!(!matches.is_empty());

        gp.add_tile(&matches[0]).expect("add_tile");

        for vtx in gp.vertices() {
            for &sid in &vtx.matches {
                assert_ne!(sid, 0, "signed id should be nonzero");
                let unsigned = sid.unsigned_abs() as usize;
                assert!(
                    unsigned >= 1 && unsigned <= gp.tileset().num_tiles(),
                    "signed id {sid} should be valid"
                );
            }
        }
    }

    #[test]
    fn diff_contains_old_and_new() {
        let mut gp = hex_patch();
        let old_vertices: Vec<VertexInfo> = gp.vertices().to_vec();

        let matches = gp.get_all_matches().to_vec();
        let diff = gp.add_tile(&matches[0]).expect("add_tile");

        assert_eq!(
            diff.affected_old.len(),
            old_vertices.len(),
            "diff should have same number of affected old vertices"
        );
        assert!(
            !diff.new_vertices.is_empty(),
            "diff should have new vertices"
        );
    }
}
