use std::collections::{BTreeSet, HashMap, VecDeque};
use std::sync::Arc;

use crate::cyclotomic::{IsComplex, IsRingOrField, SymNum, Units};
use crate::intgeom::matchtypes::MatchFinder;
use crate::intgeom::patch::{compute_glue_angles, cyclic_range_contains, PatchMatch, VertexType};
use crate::intgeom::rat::{prepare_seq, Rat};
use crate::intgeom::tileset::TileSet;
use crate::intgeom::vertextypes::{compute_gap_angle, VertexTypeKind};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct BoundaryVertexType {
    cw: (usize, usize),
    ccw: (usize, usize),
    angle: i8,
}

impl BoundaryVertexType {
    pub fn from_vertex_type<T: IsComplex + IsRingOrField + Units + SymNum>(
        vtype: &VertexType,
        tileset: &Arc<TileSet<T>>,
    ) -> Self {
        assert!(!vtype.is_empty());
        let cw = vtype[0];
        let ccw = *vtype.last().unwrap();
        let angle = compute_gap_angle(tileset, vtype);
        BoundaryVertexType { cw, ccw, angle }
    }

    pub fn cw(&self) -> (usize, usize) {
        self.cw
    }

    pub fn ccw(&self) -> (usize, usize) {
        self.ccw
    }

    pub fn angle(&self) -> i8 {
        self.angle
    }
}

pub struct BoundaryVertexTypeInfo {
    btype: BoundaryVertexType,
    kind: VertexTypeKind,
    successors: Vec<usize>,
    predecessors: Vec<usize>,
    is_cursed: bool,
}

impl BoundaryVertexTypeInfo {
    pub fn btype(&self) -> &BoundaryVertexType {
        &self.btype
    }

    pub fn kind(&self) -> VertexTypeKind {
        self.kind
    }

    pub fn successors(&self) -> &[usize] {
        &self.successors
    }

    pub fn predecessors(&self) -> &[usize] {
        &self.predecessors
    }

    pub fn is_cursed(&self) -> bool {
        self.is_cursed
    }
}

struct FrontierPatch {
    angles: Vec<i8>,
    vertex_types: Vec<VertexType>,
}

impl FrontierPatch {
    fn canonicalize(self) -> Self {
        let n = self.angles.len();
        if n == 0 {
            return self;
        }
        let (_, offset) = prepare_seq(&self.angles);
        if offset == 0 {
            return self;
        }
        let mut angles = self.angles;
        let mut vertex_types = self.vertex_types;
        angles.rotate_left(offset);
        vertex_types.rotate_left(offset);
        FrontierPatch {
            angles,
            vertex_types,
        }
    }
}

fn apply_match_to_frontier<T: IsComplex + IsRingOrField + Units>(
    fp: &FrontierPatch,
    pm: &PatchMatch,
    tileset: &Arc<TileSet<T>>,
) -> Option<FrontierPatch> {
    let n = fp.angles.len();
    let m = tileset.rat(pm.tile_id).seq().len();
    let mlen = pm.len;

    if mlen == 0 || mlen > n || mlen > m {
        return None;
    }

    let new_angles = compute_glue_angles::<T>(&fp.angles, pm, tileset)?;

    let seg_len_old = n - mlen;
    let seg_len_new = m - mlen;
    if seg_len_old + seg_len_new == 0 {
        return None;
    }

    let cw_pos = pm.start_a;
    let ccw_pos = (pm.start_a + mlen) % n;

    let mut new_vtx = Vec::with_capacity(seg_len_old + seg_len_new);

    let mut ccw_junction = fp.vertex_types[ccw_pos].clone();
    let ccw_new_offset = (pm.start_b + mlen) % m;
    ccw_junction.insert(0, (pm.tile_id, ccw_new_offset));
    new_vtx.push(ccw_junction);

    for i in 1..seg_len_old {
        let old_p = (ccw_pos + i) % n;
        new_vtx.push(fp.vertex_types[old_p].clone());
    }

    let mut cw_junction = fp.vertex_types[cw_pos].clone();
    cw_junction.push((pm.tile_id, pm.start_b));
    new_vtx.push(cw_junction);

    for k in 1..seg_len_new {
        let offset = (pm.start_b + mlen + k) % m;
        new_vtx.push(vec![(pm.tile_id, offset)]);
    }

    debug_assert_eq!(new_vtx.len(), seg_len_old + seg_len_new);

    Some(
        FrontierPatch {
            angles: new_angles,
            vertex_types: new_vtx,
        }
        .canonicalize(),
    )
}

fn make_seed_frontier<T: IsComplex + IsRingOrField + Units>(
    tileset: &Arc<TileSet<T>>,
    seed_id: usize,
    pm: &PatchMatch,
) -> Option<FrontierPatch> {
    let seed_rat = tileset.rat(seed_id);
    let n = seed_rat.seq().len();
    let m = tileset.rat(pm.tile_id).seq().len();
    let mlen = pm.len;

    if mlen == 0 || mlen > n || mlen > m {
        return None;
    }

    let seed_angles = seed_rat.seq().to_vec();
    let new_angles = compute_glue_angles::<T>(&seed_angles, pm, tileset)?;

    let seg_len_old = n - mlen;
    let seg_len_new = m - mlen;
    if seg_len_old + seg_len_new == 0 {
        return None;
    }

    let cw_pos = pm.start_a;
    let ccw_pos = (pm.start_a + mlen) % n;

    let seed_vtx: Vec<VertexType> = (0..n).map(|i| vec![(seed_id, i)]).collect();

    let mut new_vtx = Vec::with_capacity(seg_len_old + seg_len_new);

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

    debug_assert_eq!(new_vtx.len(), seg_len_old + seg_len_new);

    Some(
        FrontierPatch {
            angles: new_angles,
            vertex_types: new_vtx,
        }
        .canonicalize(),
    )
}

fn extract_bvtypes(fp: &FrontierPatch) -> Vec<BoundaryVertexType> {
    fp.vertex_types
        .iter()
        .zip(fp.angles.iter())
        .map(|(vt, &angle)| BoundaryVertexType {
            cw: vt[0],
            ccw: *vt.last().unwrap(),
            angle,
        })
        .collect()
}

pub struct BoundaryVertexTypeIndex<T: IsComplex> {
    tileset: Arc<TileSet<T>>,
    entries: Vec<BoundaryVertexTypeInfo>,
    reverse: HashMap<BoundaryVertexType, usize>,
}

impl<T: IsComplex + IsRingOrField + Units + SymNum> BoundaryVertexTypeIndex<T> {
    pub fn new(tileset: Arc<TileSet<T>>) -> Self {
        Self::new_with_max_layers(tileset, usize::MAX)
    }

    pub fn new_with_max_layers(tileset: Arc<TileSet<T>>, max_layers: usize) -> Self {
        let mut all_types: BTreeSet<BoundaryVertexType> = BTreeSet::new();
        let mut kind_map: HashMap<BoundaryVertexType, VertexTypeKind> = HashMap::new();
        let mut transitions: Vec<(BoundaryVertexType, BoundaryVertexType)> = Vec::new();

        let mut visited: BTreeSet<BoundaryVertexType> = BTreeSet::new();
        let mut queue: VecDeque<(FrontierPatch, usize)> = VecDeque::new();

        let seed_ts = Arc::new(TileSet::new(tileset.rats().to_vec()));
        eprintln!("  Building seed MatchFinder...");
        let seed_mf = MatchFinder::crossing(Arc::clone(&seed_ts), Arc::clone(&tileset));
        eprintln!("  Seed MatchFinder ready, exploring seeds...");

        let mut seed_count = 0usize;
        for seed_id in 0..tileset.num_tiles() {
            for j in 0..seed_mf.num_tiles_b() {
                for mt in seed_mf.valid_matches(seed_id, j) {
                    seed_count += 1;
                    let pm = PatchMatch {
                        start_a: mt.start_a,
                        len: mt.len,
                        start_b: mt.start_b,
                        tile_id: j,
                    };
                    if let Some(fp) = make_seed_frontier(&tileset, seed_id, &pm) {
                        let bvtypes = extract_bvtypes(&fp);
                        for (pos, bvt) in bvtypes.iter().enumerate() {
                            all_types.insert(*bvt);
                            if visited.insert(*bvt) {
                                queue.push_back((
                                    FrontierPatch {
                                        angles: fp.angles.clone(),
                                        vertex_types: fp.vertex_types.clone(),
                                    },
                                    pos,
                                ));
                            }
                        }
                    }
                }
            }
        }
        eprintln!(
            "  Seeds done: {} matches, {} visited, {} queued",
            seed_count,
            visited.len(),
            queue.len()
        );

        let mut layer = 0;
        while !queue.is_empty() {
            if layer >= max_layers {
                eprintln!(
                    "  BFS stopped early at layer {} (max_layers={})",
                    layer, max_layers
                );
                break;
            }
            let batch: Vec<_> = queue.drain(..).collect();
            layer += 1;

            let mut unique_angles: BTreeSet<Vec<i8>> = BTreeSet::new();
            for (fp, _) in &batch {
                unique_angles.insert(fp.angles.clone());
            }

            let rat_vec: Vec<Rat<T>> = unique_angles
                .iter()
                .map(|a| Rat::from_slice_unchecked(a))
                .collect();

            let batch_ts = Arc::new(TileSet::new(rat_vec));
            let mf = MatchFinder::crossing(Arc::clone(&batch_ts), Arc::clone(&tileset));

            let angles_to_tile_id: HashMap<&[i8], usize> = (0..batch_ts.num_tiles())
                .map(|id| (batch_ts.rat(id).seq(), id))
                .collect();

            let mut batch_matches: Vec<Vec<PatchMatch>> = vec![Vec::new(); batch_ts.num_tiles()];
            for rat_idx in 0..batch_ts.num_tiles() {
                for j in 0..mf.num_tiles_b() {
                    for mt in mf.valid_matches(rat_idx, j) {
                        batch_matches[rat_idx].push(PatchMatch {
                            start_a: mt.start_a,
                            len: mt.len,
                            start_b: mt.start_b,
                            tile_id: j,
                        });
                    }
                }
            }

            for (fp, pos) in &batch {
                let bvt = BoundaryVertexType {
                    cw: fp.vertex_types[*pos][0],
                    ccw: *fp.vertex_types[*pos].last().unwrap(),
                    angle: fp.angles[*pos],
                };

                if kind_map.get(&bvt) == Some(&VertexTypeKind::Dead) {
                    continue;
                }

                let rat_idx = angles_to_tile_id[fp.angles.as_slice()];
                let n = fp.angles.len();

                let touching: Vec<&PatchMatch> = batch_matches[rat_idx]
                    .iter()
                    .filter(|pm| cyclic_range_contains(pm.start_a, pm.len, *pos, n))
                    .collect();

                if touching.is_empty() {
                    kind_map.insert(bvt, VertexTypeKind::Dead);
                    continue;
                }

                kind_map.entry(bvt).or_insert(VertexTypeKind::Open);

                for pm in touching {
                    if let Some(new_fp) = apply_match_to_frontier(fp, pm, &tileset) {
                        let junction_pos = if pm.start_a == *pos {
                            n - pm.len
                        } else {
                            0
                        };
                        let new_bvtypes = extract_bvtypes(&new_fp);
                        let new_bvt = new_bvtypes[junction_pos];
                        transitions.push((bvt, new_bvt));
                        for (new_pos, nbvt) in new_bvtypes.into_iter().enumerate() {
                            all_types.insert(nbvt);
                            if visited.insert(nbvt) {
                                queue.push_back((
                                    FrontierPatch {
                                        angles: new_fp.angles.clone(),
                                        vertex_types: new_fp.vertex_types.clone(),
                                    },
                                    new_pos,
                                ));
                            }
                        }
                    }
                }
            }

            if layer % 1 == 0 {
                eprintln!(
                    "  layer {}: batch={} unique={} visited={} queue={}",
                    layer,
                    batch.len(),
                    unique_angles.len(),
                    visited.len(),
                    queue.len(),
                );
            }
        }

        eprintln!(
            "  BFS done: {} layers, {} boundary vertex types",
            layer,
            visited.len()
        );

        let entries: Vec<BoundaryVertexType> = all_types.into_iter().collect();
        let reverse: HashMap<BoundaryVertexType, usize> = entries
            .iter()
            .enumerate()
            .map(|(i, bvt)| (*bvt, i + 1))
            .collect();

        let mut succ_sets: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); entries.len()];
        let mut pred_sets: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); entries.len()];

        for (src, dst) in &transitions {
            let src_id = reverse[src];
            let dst_id = reverse[dst];
            succ_sets[src_id - 1].insert(dst_id);
            pred_sets[dst_id - 1].insert(src_id);
        }

        let is_cursed = compute_bvtype_cursed(&entries, &kind_map, &succ_sets);

        let info_entries: Vec<BoundaryVertexTypeInfo> = entries
            .into_iter()
            .enumerate()
            .map(|(i, bvt)| {
                let id = i + 1;
                BoundaryVertexTypeInfo {
                    kind: kind_map.get(&bvt).copied().unwrap_or(VertexTypeKind::Dead),
                    successors: succ_sets[i].iter().copied().collect(),
                    predecessors: pred_sets[i].iter().copied().collect(),
                    is_cursed: is_cursed[&id],
                    btype: bvt,
                }
            })
            .collect();

        BoundaryVertexTypeIndex {
            tileset,
            entries: info_entries,
            reverse,
        }
    }

    pub fn num_types(&self) -> usize {
        self.entries.len()
    }

    pub fn get_id(&self, btype: &BoundaryVertexType) -> Option<usize> {
        self.reverse.get(btype).copied()
    }

    pub fn get_info(&self, id: usize) -> &BoundaryVertexTypeInfo {
        assert!(
            id >= 1 && id <= self.entries.len(),
            "boundary vertex type id out of range"
        );
        &self.entries[id - 1]
    }

    pub fn tileset(&self) -> &Arc<TileSet<T>> {
        &self.tileset
    }

    pub fn to_rat(&self, ids: &[usize]) -> Option<Rat<T>> {
        if ids.is_empty() {
            return None;
        }
        let angles: Vec<i8> = ids.iter().map(|&id| self.entries[id - 1].btype.angle).collect();
        Some(Rat::from_slice_unchecked(&angles))
    }
}

fn compute_bvtype_cursed(
    entries: &[BoundaryVertexType],
    kind_map: &HashMap<BoundaryVertexType, VertexTypeKind>,
    succ_sets: &[BTreeSet<usize>],
) -> HashMap<usize, bool> {
    let n = entries.len();
    let mut cursed: HashMap<usize, bool> = HashMap::with_capacity(n);

    for (i, bvt) in entries.iter().enumerate() {
        let id = i + 1;
        cursed.insert(id, kind_map.get(bvt) == Some(&VertexTypeKind::Dead));
    }

    let mut changed = true;
    while changed {
        changed = false;
        for (i, _bvt) in entries.iter().enumerate() {
            let id = i + 1;
            if cursed[&id] {
                continue;
            }
            let succs = &succ_sets[i];
            if succs.is_empty() {
                continue;
            }
            if succs.iter().all(|s| cursed[s]) {
                cursed.insert(id, true);
                changed = true;
            }
        }
    }

    cursed
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic::{ZZ10, ZZ12};
    use crate::intgeom::snake::Snake;
    use crate::intgeom::tiles;
    use crate::intgeom::vertextypes::VertexTypeIndex;
    use std::collections::BTreeMap;
    use std::sync::Arc;

    fn hex_ts() -> Arc<TileSet<ZZ12>> {
        let hex: Snake<ZZ12> = tiles::hexagon();
        let rat = Rat::try_from(&hex).unwrap();
        Arc::new(TileSet::new(vec![rat]))
    }

    fn sq_ts() -> Arc<TileSet<ZZ12>> {
        let sq: Snake<ZZ12> = tiles::square();
        let rat = Rat::try_from(&sq).unwrap();
        Arc::new(TileSet::new(vec![rat]))
    }

    fn mixed_ts() -> Arc<TileSet<ZZ12>> {
        let sq: Snake<ZZ12> = tiles::square();
        let hex: Snake<ZZ12> = tiles::hexagon();
        let sq_rat = Rat::try_from(&sq).unwrap();
        let hex_rat = Rat::try_from(&hex).unwrap();
        Arc::new(TileSet::new(vec![sq_rat, hex_rat]))
    }

    fn penrose_ts() -> Arc<TileSet<ZZ10>> {
        let narrow: Snake<ZZ10> = tiles::penrose_p3_narrow();
        let wide: Snake<ZZ10> = tiles::penrose_p3_wide();
        let narrow_rat = Rat::try_from(&narrow).unwrap();
        let wide_rat = Rat::try_from(&wide).unwrap();
        Arc::new(TileSet::new(vec![narrow_rat, wide_rat]))
    }

    fn report<T: IsComplex + IsRingOrField + Units + SymNum>(
        idx: &BoundaryVertexTypeIndex<T>,
        label: &str,
    ) {
        let mut open = 0usize;
        let mut dead = 0usize;
        let mut cursed = 0usize;
        let mut initial = 0usize;
        for id in 1..=idx.num_types() {
            let info = idx.get_info(id);
            match info.kind() {
                VertexTypeKind::Open => open += 1,
                VertexTypeKind::Dead => dead += 1,
                VertexTypeKind::Closed => {}
            }
            if info.predecessors().is_empty() {
                initial += 1;
            }
            if info.is_cursed() {
                cursed += 1;
            }
        }
        eprintln!(
            "[{}] BoundaryVertexTypeIndex: {} types, open={} dead={} initial={} cursed={}",
            label, idx.num_types(), open, dead, initial, cursed
        );
    }

    #[test]
    fn hex_boundary_vtype_index() {
        let idx = BoundaryVertexTypeIndex::new(hex_ts());
        report(&idx, "hex");

        let full_idx = VertexTypeIndex::new(hex_ts());
        eprintln!(
            "  full VertexTypeIndex: {} types (for comparison)",
            full_idx.num_types()
        );
        assert_eq!(idx.num_types(), 42);
        assert!(idx.num_types() <= full_idx.num_types());
    }

    #[test]
    fn square_boundary_vtype_index() {
        let idx = BoundaryVertexTypeIndex::new(sq_ts());
        report(&idx, "square");

        let full_idx = VertexTypeIndex::new(sq_ts());
        eprintln!(
            "  full VertexTypeIndex: {} types (for comparison)",
            full_idx.num_types()
        );
        assert_eq!(idx.num_types(), 36);
        assert!(idx.num_types() <= full_idx.num_types());
    }

    #[test]
    fn mixed_boundary_vtype_index() {
        let idx = BoundaryVertexTypeIndex::new(mixed_ts());
        report(&idx, "square+hex");

        let full_idx = VertexTypeIndex::new(mixed_ts());
        eprintln!(
            "  full VertexTypeIndex: {} types (for comparison)",
            full_idx.num_types()
        );

        let mut by_kind: BTreeMap<VertexTypeKind, usize> = BTreeMap::new();
        for id in 1..=idx.num_types() {
            let info = idx.get_info(id);
            *by_kind.entry(info.kind()).or_insert(0) += 1;
        }
        eprintln!("  By kind: {:?}", by_kind);

        assert_eq!(idx.num_types(), 274);
        assert!(idx.num_types() <= full_idx.num_types());
    }

    #[test]
    fn penrose_p3_boundary_vtype_index() {
        let idx = BoundaryVertexTypeIndex::new(penrose_ts());
        report(&idx, "penrose P3");

        let mut by_kind: BTreeMap<VertexTypeKind, usize> = BTreeMap::new();
        for id in 1..=idx.num_types() {
            let info = idx.get_info(id);
            *by_kind.entry(info.kind()).or_insert(0) += 1;
        }
        eprintln!("  By kind: {:?}", by_kind);
    }
}
