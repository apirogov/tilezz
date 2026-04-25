use std::collections::{BTreeSet, HashMap, VecDeque};
use std::sync::Arc;

use crate::cyclotomic::{IsComplex, IsRingOrField, SymNum, Units};
use crate::intgeom::patch::{GrowingPatch, VertexType};
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

    pub fn from_patch_position<T: IsComplex + IsRingOrField + Units + SymNum>(
        gp: &GrowingPatch<T>,
        pos: usize,
    ) -> Self {
        let vtype = &gp.vertex_types()[pos];
        let angle = gp.angles()[pos];
        assert!(!vtype.is_empty());
        BoundaryVertexType {
            cw: vtype[0],
            ccw: *vtype.last().unwrap(),
            angle,
        }
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

pub struct BoundaryVertexTypeIndex<T: IsComplex> {
    tileset: Arc<TileSet<T>>,
    entries: Vec<BoundaryVertexTypeInfo>,
    reverse: HashMap<BoundaryVertexType, usize>,
}

impl<T: IsComplex + IsRingOrField + Units + SymNum> BoundaryVertexTypeIndex<T> {
    pub fn new(tileset: Arc<TileSet<T>>) -> Self {
        let mut all_types: BTreeSet<BoundaryVertexType> = BTreeSet::new();
        let mut kind_map: HashMap<BoundaryVertexType, VertexTypeKind> = HashMap::new();
        let mut transitions: Vec<(BoundaryVertexType, BoundaryVertexType)> = Vec::new();

        let mut visited: BTreeSet<BoundaryVertexType> = BTreeSet::new();
        let mut queue: VecDeque<(GrowingPatch<T>, usize)> = VecDeque::new();

        for seed_id in 0..tileset.num_tiles() {
            let seed = GrowingPatch::new(Arc::clone(&tileset), seed_id);
            for pm in seed.get_all_matches() {
                let mut gp = seed.clone();
                if let Some(_diff) = gp.add_tile(pm) {
                    if gp.is_growing() {
                        for pos in 0..gp.boundary_len() {
                            let bvt = BoundaryVertexType::from_patch_position(&gp, pos);
                            all_types.insert(bvt);
                            if visited.insert(bvt) {
                                queue.push_back((gp.clone(), pos));
                            }
                        }
                    }
                }
            }
        }

        while let Some((gp, pos)) = queue.pop_front() {
            let bvt = BoundaryVertexType::from_patch_position(&gp, pos);
            if kind_map.get(&bvt) == Some(&VertexTypeKind::Dead) {
                continue;
            }

            if visited.len().is_multiple_of(1000) {
                eprintln!("  BFS: visited={} queue={}", visited.len(), queue.len(),);
            }

            if visited.len().is_multiple_of(1000) {
                eprintln!("  BFS: visited={} queue={}", visited.len(), queue.len(),);
            }

            let touching: Vec<_> = gp.get_matches_touching_vertex(pos).collect();
            if touching.is_empty() {
                kind_map.insert(bvt, VertexTypeKind::Dead);
                continue;
            }

            kind_map.entry(bvt).or_insert(VertexTypeKind::Open);

            for pm in touching {
                let mut gp2 = gp.clone();
                if let Some(_diff) = gp2.add_tile(pm) {
                    if gp2.is_growing() {
                        let n = gp.boundary_len();
                        let junction_pos = if pm.start_a == pos { n - pm.len } else { 0 };
                        let new_bvt = BoundaryVertexType::from_patch_position(&gp2, junction_pos);
                        transitions.push((bvt, new_bvt));
                        for new_pos in 0..gp2.boundary_len() {
                            let nbvt = BoundaryVertexType::from_patch_position(&gp2, new_pos);
                            if visited.insert(nbvt) {
                                all_types.insert(nbvt);
                                queue.push_back((gp2.clone(), new_pos));
                            }
                        }
                    }
                }
            }
        }

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

    pub fn to_rat(&self, ids: &[usize]) -> Option<crate::intgeom::rat::Rat<T>> {
        if ids.is_empty() {
            return None;
        }
        let angles: Vec<i8> = ids
            .iter()
            .map(|&id| self.entries[id - 1].btype.angle)
            .collect();
        Some(crate::intgeom::rat::Rat::from_slice_unchecked(&angles))
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
    use crate::intgeom::rat::Rat;
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

    fn report(idx: &BoundaryVertexTypeIndex<ZZ12>, label: &str) {
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
            label,
            idx.num_types(),
            open,
            dead,
            initial,
            cursed
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

        assert!(idx.num_types() <= full_idx.num_types());
    }

    #[test]
    fn boundary_vtype_angle_reconstructs_rat() {
        let ts = hex_ts();
        let mut gp = GrowingPatch::new(Arc::clone(&ts), 0);
        let pm = gp.get_all_matches()[0].clone();
        gp.add_tile(&pm).expect("first add");
        let expected = gp.to_rat();

        let idx = BoundaryVertexTypeIndex::new(Arc::clone(&ts));
        let bvtypes: Vec<BoundaryVertexType> = (0..gp.boundary_len())
            .map(|pos| BoundaryVertexType::from_patch_position(&gp, pos))
            .collect();

        let ids: Vec<usize> = bvtypes
            .iter()
            .map(|bvt| idx.get_id(bvt).expect("should be found"))
            .collect();

        let rat = idx.to_rat(&ids).expect("reconstruct");
        assert_eq!(rat.seq(), expected.seq(), "reconstructed rat should match");
    }

    #[test]
    #[ignore]
    fn penrose_p3_boundary_vtype_index() {
        let narrow: Snake<ZZ10> = tiles::penrose_p3_narrow();
        let wide: Snake<ZZ10> = tiles::penrose_p3_wide();
        let narrow_rat = Rat::try_from(&narrow).unwrap();
        let wide_rat = Rat::try_from(&wide).unwrap();
        let ts = Arc::new(TileSet::new(vec![narrow_rat, wide_rat]));

        for seed_id in 0..ts.num_tiles() {
            let seed = GrowingPatch::new(Arc::clone(&ts), seed_id);
            let matches = seed.get_all_matches();
            eprintln!(
                "seed {} (len={}): {} matches",
                seed_id,
                ts.rat(seed_id).len(),
                matches.len()
            );
        }

        let idx = BoundaryVertexTypeIndex::new(ts);

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
            "[penrose P3] BoundaryVertexTypeIndex: {} types, open={} dead={} initial={} cursed={}",
            idx.num_types(),
            open,
            dead,
            initial,
            cursed
        );

        let mut by_kind: BTreeMap<VertexTypeKind, usize> = BTreeMap::new();
        for id in 1..=idx.num_types() {
            let info = idx.get_info(id);
            *by_kind.entry(info.kind()).or_insert(0) += 1;
        }
        eprintln!("  By kind: {:?}", by_kind);
    }
}
