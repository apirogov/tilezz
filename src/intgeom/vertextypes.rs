use std::collections::{BTreeSet, HashMap, VecDeque};
use std::sync::Arc;

use crate::cyclotomic::{IsComplex, IsRingOrField, SymNum, Units};
use crate::intgeom::patch::{GrowingPatch, PatchMatch, VertexType};
use crate::intgeom::tileset::TileSet;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum VertexTypeKind {
    Open,
    Closed,
    Dead,
}

pub struct VertexTypeInfo {
    vtype: VertexType,
    kind: VertexTypeKind,
    successors: Vec<usize>,
    predecessors: Vec<usize>,
    is_cursed: bool,
}

impl VertexTypeInfo {
    pub fn vtype(&self) -> &VertexType {
        &self.vtype
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

    pub fn is_initial(&self) -> bool {
        self.predecessors.is_empty()
    }

    pub fn is_terminal(&self) -> bool {
        self.successors.is_empty()
    }

    pub fn is_closed(&self) -> bool {
        self.kind == VertexTypeKind::Closed
    }

    pub fn is_cursed(&self) -> bool {
        self.is_cursed
    }
}

pub struct VertexTypeIndex<T: IsComplex> {
    tileset: Arc<TileSet<T>>,
    entries: Vec<VertexTypeInfo>,
    reverse: HashMap<VertexType, usize>,
}

impl<T: IsComplex + IsRingOrField + Units + SymNum> VertexTypeIndex<T> {
    pub fn new(tileset: Arc<TileSet<T>>) -> Self {
        let mut all_types: BTreeSet<VertexType> = BTreeSet::new();
        let mut kind_map: HashMap<VertexType, VertexTypeKind> = HashMap::new();
        let mut transitions: Vec<(VertexType, VertexType)> = Vec::new();

        let mut visited: BTreeSet<VertexType> = BTreeSet::new();
        let mut queue: VecDeque<(GrowingPatch<T>, usize)> = VecDeque::new();

        for seed_id in 0..tileset.num_tiles() {
            let seed = GrowingPatch::new(Arc::clone(&tileset), seed_id);
            for pm in seed.get_all_matches() {
                let mut gp = seed.clone();
                if let Some(diff) = gp.add_tile(pm) {
                    for cvt in &diff.closed_vertex_types {
                        all_types.insert(cvt.clone());
                        kind_map.insert(cvt.clone(), VertexTypeKind::Closed);
                    }
                    if gp.is_growing() {
                        for pos in 0..gp.boundary_len() {
                            let vt = gp.vertex_types()[pos].clone();
                            if visited.insert(vt.clone()) {
                                all_types.insert(vt.clone());
                                queue.push_back((gp.clone(), pos));
                            }
                        }
                    }
                }
            }
        }

        while let Some((gp, pos)) = queue.pop_front() {
            let vt = gp.vertex_types()[pos].clone();
            if kind_map.get(&vt) == Some(&VertexTypeKind::Dead) {
                continue;
            }

            let touching: Vec<&PatchMatch> = gp.get_matches_touching_vertex(pos).collect();
            if touching.is_empty() {
                kind_map.insert(vt, VertexTypeKind::Dead);
                continue;
            }

            kind_map.entry(vt.clone()).or_insert(VertexTypeKind::Open);

            for pm in touching {
                let mut gp2 = gp.clone();
                if let Some(diff) = gp2.add_tile(pm) {
                    for cvt in &diff.closed_vertex_types {
                        all_types.insert(cvt.clone());
                        kind_map.insert(cvt.clone(), VertexTypeKind::Closed);
                    }
                    if gp2.is_growing() {
                        let n = gp.boundary_len();
                        let junction_pos = if pm.start_a == pos { n - pm.len } else { 0 };
                        let new_vt = gp2.vertex_types()[junction_pos].clone();
                        transitions.push((vt.clone(), new_vt.clone()));
                        for new_pos in 0..gp2.boundary_len() {
                            let nv = gp2.vertex_types()[new_pos].clone();
                            if nv.len() >= vt.len() && visited.insert(nv) {
                                all_types.insert(gp2.vertex_types()[new_pos].clone());
                                queue.push_back((gp2.clone(), new_pos));
                            }
                        }
                    }
                }
            }
        }

        validate_closed_types(&all_types, &kind_map, &tileset);

        let entries: Vec<VertexType> = all_types.into_iter().collect();
        let reverse: HashMap<VertexType, usize> = entries
            .iter()
            .enumerate()
            .map(|(i, vt)| (vt.clone(), i + 1))
            .collect();

        let mut succ_sets: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); entries.len()];
        let mut pred_sets: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); entries.len()];

        for (src, dst) in &transitions {
            let src_id = reverse[src];
            let dst_id = reverse[dst];
            succ_sets[src_id - 1].insert(dst_id);
            pred_sets[dst_id - 1].insert(src_id);
        }

        let is_cursed = compute_cursed(&entries, &kind_map, &succ_sets, &reverse);

        let info_entries: Vec<VertexTypeInfo> = entries
            .into_iter()
            .enumerate()
            .map(|(i, vt)| {
                let id = i + 1;
                VertexTypeInfo {
                    kind: kind_map[&vt],
                    successors: succ_sets[i].iter().copied().collect(),
                    predecessors: pred_sets[i].iter().copied().collect(),
                    is_cursed: is_cursed[&id],
                    vtype: vt,
                }
            })
            .collect();

        VertexTypeIndex {
            tileset,
            entries: info_entries,
            reverse,
        }
    }

    pub fn num_types(&self) -> usize {
        self.entries.len()
    }

    pub fn get_id(&self, vtype: &VertexType) -> Option<usize> {
        self.reverse.get(vtype).copied()
    }

    pub fn get_type(&self, id: usize) -> &VertexType {
        assert!(
            id >= 1 && id <= self.entries.len(),
            "vertex type id out of range"
        );
        &self.entries[id - 1].vtype
    }

    pub fn get_info(&self, id: usize) -> &VertexTypeInfo {
        assert!(
            id >= 1 && id <= self.entries.len(),
            "vertex type id out of range"
        );
        &self.entries[id - 1]
    }

    pub fn tileset(&self) -> &Arc<TileSet<T>> {
        &self.tileset
    }
}

fn validate_closed_types<T: IsComplex + IsRingOrField + Units + SymNum>(
    all_types: &BTreeSet<VertexType>,
    kind_map: &HashMap<VertexType, VertexTypeKind>,
    tileset: &Arc<TileSet<T>>,
) {
    let hturn = T::hturn() as i64;
    for vt in all_types {
        if kind_map[vt] != VertexTypeKind::Closed {
            continue;
        }
        let angle_sum: i64 = vt
            .iter()
            .map(|&(tile_id, offset)| tileset.rat(tile_id).seq()[offset] as i64)
            .sum();
        let expected = (vt.len() as i64 - 2) * hturn;
        assert_eq!(
            angle_sum, expected,
            "closed vertex type {:?}: angle sum {} != expected {}",
            vt, angle_sum, expected
        );
    }
}

fn compute_cursed(
    entries: &[VertexType],
    kind_map: &HashMap<VertexType, VertexTypeKind>,
    succ_sets: &[BTreeSet<usize>],
    _reverse: &HashMap<VertexType, usize>,
) -> HashMap<usize, bool> {
    let n = entries.len();
    let mut cursed: HashMap<usize, bool> = HashMap::with_capacity(n);

    for (i, vt) in entries.iter().enumerate() {
        let id = i + 1;
        cursed.insert(id, kind_map[vt] == VertexTypeKind::Dead);
    }

    let mut changed = true;
    while changed {
        changed = false;
        for (i, vt) in entries.iter().enumerate() {
            let id = i + 1;
            if cursed[&id] {
                continue;
            }
            if kind_map[vt] == VertexTypeKind::Closed {
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
