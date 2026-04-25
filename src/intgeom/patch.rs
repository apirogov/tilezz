use std::collections::BTreeSet;
use std::sync::Arc;

use crate::cyclotomic::{IsComplex, IsRingOrField, Units};
use crate::intgeom::angles::normalize_angle;
use crate::intgeom::matchtypes::MatchFinder;
use crate::intgeom::rat::Rat;
use crate::intgeom::tileset::TileSet;

pub type VertexType = Vec<(usize, usize)>;

#[derive(Debug, Clone)]
pub struct PatchMatch {
    pub start_a: usize,
    pub len: usize,
    pub start_b: usize,
    pub tile_id: usize,
}

pub struct AddTileDiff {
    pub old_local: Vec<VertexType>,
    pub new_local: Vec<VertexType>,
    pub closed_vertex_types: Vec<VertexType>,
}

#[derive(Clone)]
pub struct GrowingPatch<T: IsComplex> {
    tileset: Arc<TileSet<T>>,
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
        vertex_types: Vec<VertexType>,
        cached_matches: Vec<PatchMatch>,
    },
    Closed,
    _Phantom(std::marker::PhantomData<T>),
}

impl<T: IsComplex + IsRingOrField + Units> GrowingPatch<T> {
    pub fn new(tileset: Arc<TileSet<T>>, seed_tile_id: usize) -> Self {
        let seed_matches = Self::compute_seed_matches(&tileset, seed_tile_id);
        GrowingPatch {
            tileset,
            state: PatchState::Seed {
                tile_id: seed_tile_id,
                cached_matches: seed_matches,
            },
        }
    }

    pub fn seed_tile_id(&self) -> usize {
        match &self.state {
            PatchState::Seed { tile_id, .. } => *tile_id,
            PatchState::Growing { vertex_types, .. } => vertex_types
                .first()
                .and_then(|v| v.first())
                .map(|(id, _)| *id)
                .unwrap_or(0),
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

    pub fn get_matches_touching_vertex(
        &self,
        vertex_index: usize,
    ) -> impl Iterator<Item = &PatchMatch> {
        let n = self.boundary_len();
        self.get_all_matches()
            .iter()
            .filter(move |pm| cyclic_range_contains(pm.start_a, pm.len, vertex_index, n))
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

    pub fn vertex_types(&self) -> &[VertexType] {
        match &self.state {
            PatchState::Growing { vertex_types, .. } => vertex_types,
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
                self.init_from_first_add(seed_id, pm)
            }
            PatchState::Growing { .. } => self.add_tile_growing(pm),
            PatchState::Closed | PatchState::_Phantom(_) => None,
        }
    }

    fn init_from_first_add(&mut self, seed_id: usize, pm: &PatchMatch) -> Option<AddTileDiff> {
        let seed_rat = self.tileset.rat(seed_id);
        let n = seed_rat.seq().len();
        let m = self.tileset.rat(pm.tile_id).seq().len();
        let mlen = pm.len;

        if mlen == 0 || mlen > n || mlen > m {
            return None;
        }

        let (_ns, seed_len, _ne) = seed_rat.get_match(
            (pm.start_a as i64, pm.start_b as i64),
            self.tileset.rat(pm.tile_id),
        );
        if seed_len == 0 {
            return None;
        }

        let seed_angles = seed_rat.seq().to_vec();
        let new_angles = compute_glue_angles::<T>(&seed_angles, pm, &self.tileset)?;

        let seg_len_old = n - mlen;
        let seg_len_new = m - mlen;
        let new_len = seg_len_old + seg_len_new;
        if new_len == 0 {
            return None;
        }

        let cw_pos = pm.start_a;
        let ccw_pos = (pm.start_a + mlen) % n;

        let seed_vtx: Vec<VertexType> = (0..n).map(|i| vec![(seed_id, i)]).collect();

        let mut closed = Vec::new();
        for k in 1..mlen {
            let p = (cw_pos + k) % n;
            let mut vt = seed_vtx[p].clone();
            let new_offset = (pm.start_b + k) % m;
            vt.push((pm.tile_id, new_offset));
            canonicalize_vtx(&mut vt);
            closed.push(vt);
        }

        let mut old_local = Vec::new();
        for k in 0..=mlen {
            let p = (cw_pos + k) % n;
            old_local.push(seed_vtx[p].clone());
        }

        let mut new_vtx = Vec::with_capacity(new_len);

        let mut ccw_junction = seed_vtx[ccw_pos].clone();
        let ccw_new_offset = (pm.start_b + mlen) % m;
        ccw_junction.insert(0, (pm.tile_id, ccw_new_offset));
        new_vtx.push(ccw_junction);

        for i in 1..seg_len_old {
            let old_p = (ccw_pos + i) % n;
            new_vtx.push(seed_vtx[old_p].clone());
        }

        let mut cw_junction = seed_vtx[cw_pos].clone();
        cw_junction.push((pm.tile_id, pm.start_b));
        new_vtx.push(cw_junction);

        for k in 1..seg_len_new {
            let offset = (pm.start_b + mlen + k) % m;
            new_vtx.push(vec![(pm.tile_id, offset)]);
        }

        debug_assert_eq!(new_vtx.len(), new_len);

        let mut new_local = Vec::new();
        new_local.push(new_vtx[0].clone());
        for k in 1..seg_len_new {
            new_local.push(new_vtx[seg_len_old + k].clone());
        }
        new_local.push(new_vtx[seg_len_old].clone());

        self.state = PatchState::Growing {
            angles: new_angles,
            vertex_types: new_vtx,
            cached_matches: Vec::new(),
        };
        self.recompute_matches();

        Some(AddTileDiff {
            old_local,
            new_local,
            closed_vertex_types: closed,
        })
    }

    fn add_tile_growing(&mut self, pm: &PatchMatch) -> Option<AddTileDiff> {
        let (angles, vertex_types) = match &mut self.state {
            PatchState::Growing {
                angles,
                vertex_types,
                ..
            } => (std::mem::take(angles), std::mem::take(vertex_types)),
            _ => return None,
        };

        let n = angles.len();
        let mlen = pm.len;
        let m = self.tileset.rat(pm.tile_id).seq().len();

        if mlen == 0 || mlen > n || mlen > m {
            self.restore_growing(angles, vertex_types);
            return None;
        }

        let new_angles = match compute_glue_angles::<T>(&angles, pm, &self.tileset) {
            Some(a) => a,
            None => {
                self.restore_growing(angles, vertex_types);
                return None;
            }
        };

        let seg_len_old = n - mlen;
        let seg_len_new = m - mlen;
        let new_len = seg_len_old + seg_len_new;

        if new_len == 0 {
            self.state = PatchState::Closed;
            return Some(AddTileDiff {
                old_local: Vec::new(),
                new_local: Vec::new(),
                closed_vertex_types: Vec::new(),
            });
        }

        let cw_pos = pm.start_a;
        let ccw_pos = (pm.start_a + mlen) % n;

        let mut closed = Vec::new();
        for k in 1..mlen {
            let p = (cw_pos + k) % n;
            let mut vt = vertex_types[p].clone();
            let new_offset = (pm.start_b + k) % m;
            vt.push((pm.tile_id, new_offset));
            canonicalize_vtx(&mut vt);
            closed.push(vt);
        }

        let mut old_local = Vec::new();
        for k in 0..=mlen {
            let p = (cw_pos + k) % n;
            old_local.push(vertex_types[p].clone());
        }

        let mut new_vtx = Vec::with_capacity(new_len);

        let mut ccw_junction = vertex_types[ccw_pos].clone();
        let ccw_new_offset = (pm.start_b + mlen) % m;
        ccw_junction.insert(0, (pm.tile_id, ccw_new_offset));
        new_vtx.push(ccw_junction);

        for i in 1..seg_len_old {
            let old_p = (ccw_pos + i) % n;
            new_vtx.push(vertex_types[old_p].clone());
        }

        let mut cw_junction = vertex_types[cw_pos].clone();
        cw_junction.push((pm.tile_id, pm.start_b));
        new_vtx.push(cw_junction);

        for k in 1..seg_len_new {
            let offset = (pm.start_b + mlen + k) % m;
            new_vtx.push(vec![(pm.tile_id, offset)]);
        }

        debug_assert_eq!(new_vtx.len(), new_len);

        let mut new_local = Vec::new();
        new_local.push(new_vtx[0].clone());
        for k in 1..seg_len_new {
            new_local.push(new_vtx[seg_len_old + k].clone());
        }
        new_local.push(new_vtx[seg_len_old].clone());

        self.state = PatchState::Growing {
            angles: new_angles,
            vertex_types: new_vtx,
            cached_matches: Vec::new(),
        };
        self.recompute_matches();

        Some(AddTileDiff {
            old_local,
            new_local,
            closed_vertex_types: closed,
        })
    }

    fn restore_growing(&mut self, angles: Vec<i8>, vertex_types: Vec<VertexType>) {
        self.state = PatchState::Growing {
            angles,
            vertex_types,
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

fn cyclic_range_contains(start: usize, len: usize, index: usize, n: usize) -> bool {
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

fn canonicalize_vtx(vt: &mut Vec<(usize, usize)>) {
    let n = vt.len();
    if n <= 1 {
        return;
    }
    let best = lex_min_rotation_clone(vt);
    *vt = best;
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

fn lex_min_rotation_clone(vt: &[(usize, usize)]) -> VertexType {
    let n = vt.len();
    if n <= 1 {
        return vt.to_vec();
    }
    let mut best = vt.to_vec();
    for r in 1..n {
        let candidate: Vec<(usize, usize)> = (0..n).map(|i| vt[(r + i) % n]).collect();
        if candidate < best {
            best = candidate;
        }
    }
    best
}

fn synthetic_closed_vtypes(
    vertex_len: usize,
    tile_id: usize,
    tile_edges: usize,
) -> BTreeSet<VertexType> {
    let seed: Vec<VertexType> = (0..tile_edges).map(|i| vec![(tile_id, i)]).collect();
    let mut result: BTreeSet<VertexType> = BTreeSet::new();
    let mut combo: Vec<(usize, usize)> = Vec::with_capacity(vertex_len);
    fn recurse(
        depth: usize,
        max_depth: usize,
        seed: &[VertexType],
        combo: &mut Vec<(usize, usize)>,
        result: &mut BTreeSet<VertexType>,
    ) {
        if depth == max_depth {
            result.insert(lex_min_rotation_clone(combo));
            return;
        }
        for s in seed {
            combo.push(s[0]);
            recurse(depth + 1, max_depth, seed, combo, result);
            combo.pop();
        }
    }
    recurse(0, vertex_len, &seed, &mut combo, &mut result);
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic::{ZZ12, ZZ4};
    use crate::intgeom::snake::Snake;
    use crate::intgeom::tiles;
    use std::collections::BTreeMap;
    use std::collections::BTreeSet;
    use std::collections::VecDeque;

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
        assert!(!gp.is_closed());
        assert_eq!(gp.boundary_len(), 0);
        assert_eq!(gp.vertex_types().len(), 0);
    }

    #[test]
    fn seed_patch_has_matches() {
        let gp = square_patch();
        let matches = gp.get_all_matches();
        assert!(!matches.is_empty(), "seed should have matches");
        for pm in matches {
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
        assert_eq!(gp.vertex_types().len(), gp.boundary_len());

        assert!(diff.closed_vertex_types.is_empty() || pm.len > 2);
        assert!(!diff.old_local.is_empty());
        assert!(!diff.new_local.is_empty());
    }

    #[test]
    fn segment_cyclic_invariant() {
        let mut gp = hex_patch();
        let matches = gp.get_all_matches().to_vec();
        let pm = matches[0].clone();
        gp.add_tile(&pm).expect("first add");
        assert_eq!(
            gp.angles().len(),
            gp.vertex_types().len(),
            "angles == vertex_types"
        );

        let matches2 = gp.get_all_matches().to_vec();
        if let Some(pm2) = matches2.first() {
            if gp.add_tile(pm2).is_some() {
                assert_eq!(
                    gp.angles().len(),
                    gp.vertex_types().len(),
                    "angles == vertex_types after second add"
                );
            }
        }
    }

    #[test]
    fn junction_vertex_ids_nonempty_after_each_add() {
        let mut gp = hex_patch();
        let matches = gp.get_all_matches().to_vec();

        for (step, pm) in matches.iter().enumerate() {
            let diff = gp.add_tile(pm);
            if diff.is_none() {
                break;
            }
            let vt = gp.vertex_types();
            assert!(
                !vt.is_empty(),
                "step {}: vertex types should not be empty",
                step
            );
            if step < 3 {
                assert!(
                    !vt[0].is_empty(),
                    "step {}: first vertex type should be non-empty",
                    step
                );
                let last = vt.len() - 1;
                assert!(
                    !vt[last].is_empty(),
                    "step {}: last vertex type should be non-empty",
                    step
                );
            }
        }
    }

    #[test]
    fn hexagon_all_36_matches_produce_valid_bi_hexes() {
        let gp = hex_patch();
        let matches = gp.get_all_matches().to_vec();
        assert_eq!(matches.len(), 36, "hex self-matches = 36");

        for pm in &matches {
            let mut gp2 = hex_patch();
            let diff = gp2.add_tile(pm).expect("first add");
            assert!(gp2.is_growing());
            assert_eq!(gp2.boundary_len(), 12 - 2 * pm.len);
            assert_eq!(gp2.vertex_types().len(), gp2.boundary_len());

            let rat = gp2.to_rat();
            assert!(
                Snake::<ZZ12>::try_from(rat.seq()).is_ok(),
                "valid snake for pm {:?}",
                pm
            );

            assert_eq!(diff.old_local.len(), pm.len + 1);
            assert!(diff.new_local.len() >= 2);
        }
    }

    #[test]
    fn square_all_16_matches_produce_valid_bi_squares() {
        let gp = square_patch();
        let matches = gp.get_all_matches().to_vec();
        assert_eq!(matches.len(), 16, "square self-matches = 16");

        for pm in &matches {
            let mut gp2 = square_patch();
            let diff = gp2.add_tile(pm).expect("first add");
            assert!(gp2.is_growing());
            assert_eq!(gp2.boundary_len(), 8 - 2 * pm.len);

            let rat = gp2.to_rat();
            assert!(
                Snake::<ZZ4>::try_from(rat.seq()).is_ok(),
                "valid snake for pm {:?}",
                pm
            );

            assert_eq!(diff.old_local.len(), pm.len + 1);
            assert!(diff.new_local.len() >= 2);
        }
    }

    #[test]
    fn to_rat_matches_direct_glue_for_all_matches() {
        let gp = hex_patch();
        let matches = gp.get_all_matches().to_vec();
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

        let matches = gp.get_all_matches().to_vec();
        assert!(!matches.is_empty());

        if let Some(pm) = matches.first() {
            let diff = gp.add_tile(pm);
            assert!(diff.is_some());
            assert!(gp.is_growing());
        }
    }

    struct PatchCensus<T: IsComplex> {
        ts: Arc<TileSet<T>>,
        patches: BTreeMap<Rat<T>, Vec<Vec<PatchMatch>>>,
    }

    impl PatchCensus<ZZ4> {
        fn for_squares(max_tiles: usize) -> Self {
            let sq: Snake<ZZ4> = tiles::square();
            let rat = Rat::try_from(&sq).unwrap();
            let ts = Arc::new(TileSet::new(vec![rat]));
            let patches = brute_force_zz4(&ts, max_tiles);
            PatchCensus { ts, patches }
        }

        fn synthetic_closed_vtypes_len4(&self) -> BTreeSet<VertexType> {
            synthetic_closed_vtypes(4, 0, 4)
        }
    }

    impl PatchCensus<ZZ12> {
        fn for_hexagons(max_tiles: usize) -> Self {
            let hex: Snake<ZZ12> = tiles::hexagon();
            let rat = Rat::try_from(&hex).unwrap();
            let ts = Arc::new(TileSet::new(vec![rat]));
            let patches = brute_force_zz12(&ts, max_tiles);
            PatchCensus { ts, patches }
        }

        fn synthetic_closed_vtypes_len3(&self) -> BTreeSet<VertexType> {
            synthetic_closed_vtypes(3, 0, 6)
        }
    }

    impl<T: IsComplex + IsRingOrField + Units> PatchCensus<T> {
        fn report(&self, label: &str) {
            let mut by_tiles: BTreeMap<usize, usize> = BTreeMap::new();
            for ways in self.patches.values() {
                let n = ways[0].len() + 1;
                by_tiles.entry(n).and_modify(|c| *c += 1).or_insert(1);
            }
            eprintln!("[{}] Patch census:", label);
            for (n, count) in &by_tiles {
                eprintln!("  {}-tile patches: {}", n, count);
            }

            eprintln!("\nDetailed patches by tile count:");
            for (rat, ways) in &self.patches {
                let n = ways[0].len() + 1;
                if n <= 3 {
                    eprintln!("  n={} boundary={:?} ({} ways)", n, rat.seq(), ways.len());
                }
            }

            let mut open_vtypes: BTreeSet<VertexType> = BTreeSet::new();
            let mut closed_vtypes: BTreeSet<VertexType> = BTreeSet::new();

            for ways in self.patches.values() {
                for history in ways {
                    let mut gp = GrowingPatch::new(Arc::clone(&self.ts), 0);
                    let mut all_closed: Vec<VertexType> = Vec::new();
                    for pm in history {
                        if let Some(diff) = gp.add_tile(pm) {
                            all_closed.extend(diff.closed_vertex_types);
                        }
                    }
                    if gp.is_growing() {
                        for vt in gp.vertex_types() {
                            open_vtypes.insert(vt.clone());
                        }
                    }
                    closed_vtypes.extend(all_closed);
                }
            }

            let mut open_by_n: BTreeMap<usize, usize> = BTreeMap::new();
            for vt in &open_vtypes {
                *open_by_n.entry(vt.len()).or_insert(0) += 1;
            }
            eprintln!("\n  Open vertex types by touched-tile count:");
            for (n, count) in &open_by_n {
                eprintln!("    touching {} tiles: {}", n, count);
            }

            let mut closed_by_n: BTreeMap<usize, usize> = BTreeMap::new();
            for vt in &closed_vtypes {
                *closed_by_n.entry(vt.len()).or_insert(0) += 1;
            }
            eprintln!("\n  Closed vertex types by touched-tile count:");
            for (n, count) in &closed_by_n {
                eprintln!("    touching {} tiles: {}", n, count);
            }
        }

        fn verify_all_steps(&self, label: &str) {
            let mut checked_steps = 0usize;
            let mut checked_closed = 0usize;

            for ways in self.patches.values() {
                for history in ways {
                    if history.is_empty() {
                        continue;
                    }
                    let mut gp = GrowingPatch::new(Arc::clone(&self.ts), 0);

                    for (step, pm) in history.iter().enumerate() {
                        let old_n = gp.boundary_len();
                        let old_vtypes = if old_n > 0 {
                            gp.vertex_types().to_vec()
                        } else {
                            Vec::new()
                        };

                        let diff = gp.add_tile(pm);
                        assert!(diff.is_some(), "[{}] step {} add_tile failed", label, step);
                        let diff = diff.unwrap();

                        if old_n > 0 {
                            assert_eq!(
                                diff.old_local.len(),
                                pm.len + 1,
                                "[{}] step {}: old_local len should be mlen+1={}",
                                label,
                                step,
                                pm.len + 1
                            );

                            let cw_pos = pm.start_a;
                            for k in 0..=pm.len {
                                let p = (cw_pos + k) % old_n;
                                assert_eq!(
                                    diff.old_local[k], old_vtypes[p],
                                    "[{}] step {}: old_local[{}] mismatch at pos {}",
                                    label, step, k, p
                                );
                            }

                            let m = self.ts.rat(pm.tile_id).seq().len();

                            if pm.len > 1 {
                                assert_eq!(
                                    diff.closed_vertex_types.len(),
                                    pm.len - 1,
                                    "[{}] step {}: closed count",
                                    label,
                                    step
                                );
                                for (ck, cvt) in diff.closed_vertex_types.iter().enumerate() {
                                    let k = ck + 1;
                                    let p = (cw_pos + k) % old_n;
                                    let mut expected = old_vtypes[p].clone();
                                    expected.push((pm.tile_id, (pm.start_b + k) % m));
                                    let expected_canon = lex_min_rotation_clone(&expected);
                                    assert_eq!(
                                        *cvt, expected_canon,
                                        "[{}] step {}: closed vtx {} at pos {}",
                                        label, step, ck, p
                                    );
                                }
                            } else {
                                assert!(
                                    diff.closed_vertex_types.is_empty(),
                                    "[{}] step {}: no closed for mlen=1",
                                    label,
                                    step
                                );
                            }

                            assert_eq!(
                                gp.angles().len(),
                                gp.vertex_types().len(),
                                "[{}] step {}: angles==vtx_types",
                                label,
                                step
                            );
                        }

                        checked_steps += 1;
                        checked_closed += diff.closed_vertex_types.len();
                    }

                    assert_eq!(
                        gp.angles().len(),
                        gp.vertex_types().len(),
                        "[{}] final cyclic invariant",
                        label
                    );
                }
            }
            eprintln!(
                "  [{}] verified {} steps, {} closed vertices",
                label, checked_steps, checked_closed
            );
        }
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

        results
            .entry(ts.rat(0).clone())
            .or_default()
            .push(Vec::new());

        let seed_matches = GrowingPatch::new(Arc::clone(ts), 0)
            .get_all_matches()
            .to_vec();
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

        results
            .entry(ts.rat(0).clone())
            .or_default()
            .push(Vec::new());

        let seed_matches = GrowingPatch::new(Arc::clone(ts), 0)
            .get_all_matches()
            .to_vec();
        for pm in &seed_matches {
            let mut gp = GrowingPatch::new(Arc::clone(ts), 0);
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

        let synthetic = census.synthetic_closed_vtypes_len4();
        eprintln!(
            "\n  Synthetic closed vtypes of length 4: {}",
            synthetic.len()
        );
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

        let synthetic = census.synthetic_closed_vtypes_len3();
        eprintln!(
            "\n  Synthetic closed vtypes of length 3: {}",
            synthetic.len()
        );

        for ways in census.patches.values() {
            for history in ways {
                if history.is_empty() {
                    continue;
                }
                let mut gp = GrowingPatch::new(Arc::clone(&census.ts), 0);
                for pm in history {
                    gp.add_tile(pm).expect("add_tile");
                }
                assert!(gp.is_growing(), "patch should be growing");
                assert_eq!(
                    gp.angles().len(),
                    gp.vertex_types().len(),
                    "angles == vertex_types"
                );
                assert!(
                    Snake::<ZZ12>::try_from(gp.to_rat().seq()).is_ok(),
                    "valid snake"
                );
            }
        }
    }

    struct VertexTypeCollector<T: IsComplex> {
        tileset: Arc<TileSet<T>>,
        open_types: BTreeSet<VertexType>,
        closed_types: BTreeSet<VertexType>,
        dead_types: BTreeSet<VertexType>,
    }

    impl<T: IsComplex + IsRingOrField + Units> VertexTypeCollector<T> {
        fn collect(tileset: Arc<TileSet<T>>) -> Self {
            let mut open_types: BTreeSet<VertexType> = BTreeSet::new();
            let mut closed_types: BTreeSet<VertexType> = BTreeSet::new();
            let mut dead_types: BTreeSet<VertexType> = BTreeSet::new();
            let mut visited: BTreeSet<VertexType> = BTreeSet::new();
            let mut queue: VecDeque<(GrowingPatch<T>, usize)> = VecDeque::new();

            for seed_id in 0..tileset.num_tiles() {
                let seed = GrowingPatch::new(Arc::clone(&tileset), seed_id);
                for pm in seed.get_all_matches() {
                    let mut gp = seed.clone();
                    if let Some(diff) = gp.add_tile(pm) {
                        closed_types.extend(diff.closed_vertex_types);
                        if gp.is_growing() {
                            for pos in 0..gp.boundary_len() {
                                let vt = gp.vertex_types()[pos].clone();
                                if visited.insert(vt.clone()) {
                                    queue.push_back((gp.clone(), pos));
                                }
                            }
                        }
                    }
                }
            }

            while let Some((gp, pos)) = queue.pop_front() {
                let vt = gp.vertex_types()[pos].clone();
                if dead_types.contains(&vt) {
                    continue;
                }

                let touching: Vec<&PatchMatch> = gp.get_matches_touching_vertex(pos).collect();
                if touching.is_empty() {
                    dead_types.insert(vt);
                    continue;
                }

                open_types.insert(vt.clone());

                for pm in touching {
                    let mut gp2 = gp.clone();
                    if let Some(diff) = gp2.add_tile(pm) {
                        closed_types.extend(diff.closed_vertex_types);
                        if gp2.is_growing() {
                            for new_pos in 0..gp2.boundary_len() {
                                let new_vt = gp2.vertex_types()[new_pos].clone();
                                if new_vt.len() >= vt.len() && visited.insert(new_vt.clone()) {
                                    queue.push_back((gp2.clone(), new_pos));
                                }
                            }
                        }
                    }
                }
            }

            VertexTypeCollector {
                tileset,
                open_types,
                closed_types,
                dead_types,
            }
        }

        fn report(&self, label: &str) {
            let mut open_by_n: BTreeMap<usize, usize> = BTreeMap::new();
            for vt in &self.open_types {
                *open_by_n.entry(vt.len()).or_insert(0) += 1;
            }
            let mut closed_by_n: BTreeMap<usize, usize> = BTreeMap::new();
            for vt in &self.closed_types {
                *closed_by_n.entry(vt.len()).or_insert(0) += 1;
            }
            let mut dead_by_n: BTreeMap<usize, usize> = BTreeMap::new();
            for vt in &self.dead_types {
                *dead_by_n.entry(vt.len()).or_insert(0) += 1;
            }

            eprintln!("[{}] Vertex type collection:", label);
            eprintln!("  Open (non-terminal) by touched-tile count:");
            for (n, c) in &open_by_n {
                eprintln!("    {} tiles: {}", n, c);
            }
            eprintln!("  Open total: {}", self.open_types.len());
            eprintln!("  Closed (terminal) by touched-tile count:");
            for (n, c) in &closed_by_n {
                eprintln!("    {} tiles: {}", n, c);
            }
            eprintln!("  Closed total: {}", self.closed_types.len());
            eprintln!("  Dead (open, terminal) by touched-tile count:");
            for (n, c) in &dead_by_n {
                eprintln!("    {} tiles: {}", n, c);
            }
            eprintln!("  Dead total: {}", self.dead_types.len());
        }
    }

    #[test]
    fn square_systematic_vtypes() {
        let sq: Snake<ZZ4> = tiles::square();
        let rat = Rat::try_from(&sq).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let collector = VertexTypeCollector::collect(ts);
        collector.report("square");

        let synthetic_closed = synthetic_closed_vtypes(4, 0, 4);
        assert_eq!(
            collector
                .closed_types
                .iter()
                .filter(|vt| vt.len() == 4)
                .count(),
            synthetic_closed.len(),
            "closed vtypes of length 4 should match synthetic"
        );
    }

    #[test]
    fn hex_systematic_vtypes() {
        let hex: Snake<ZZ12> = tiles::hexagon();
        let rat = Rat::try_from(&hex).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let collector = VertexTypeCollector::collect(ts);
        collector.report("hexagon");

        let synthetic_closed = synthetic_closed_vtypes(3, 0, 6);
        assert_eq!(
            collector
                .closed_types
                .iter()
                .filter(|vt| vt.len() == 3)
                .count(),
            synthetic_closed.len(),
            "closed vtypes of length 3 should match synthetic"
        );
    }
}
