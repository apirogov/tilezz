use std::collections::{BTreeSet, HashMap, VecDeque};
use std::sync::Arc;

use crate::cyclotomic::{IsComplex, IsRingOrField, Units};
use crate::intgeom::patch::{GrowingPatch, PatchMatch, PatchVertexType};
use crate::intgeom::rat::Rat;
use crate::intgeom::tileset::TileSet;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum VertexTypeKind {
    Open,
    Dead,
}

pub struct VertexTypeInfo<T: IsComplex> {
    vtype: PatchVertexType,
    kind: VertexTypeKind,
    successors: Vec<usize>,
    predecessors: Vec<usize>,
    is_cursed: bool,
    realizing_rat: Rat<T>,
    witness: GrowingPatch<T>,
    witness_pos: usize,
}

impl<T: IsComplex> VertexTypeInfo<T> {
    pub fn vtype(&self) -> &PatchVertexType {
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

    pub fn is_cursed(&self) -> bool {
        self.is_cursed
    }

    pub fn realizing_rat(&self) -> &Rat<T> {
        &self.realizing_rat
    }

    pub fn gap_angle(&self) -> i8 {
        self.vtype.angle
    }

    pub fn witness(&self) -> &GrowingPatch<T> {
        &self.witness
    }

    pub fn witness_pos(&self) -> usize {
        self.witness_pos
    }
}

pub struct TransitionInfo {
    pub src_id: usize,
    pub dst_id: usize,
    pub patch_match: PatchMatch,
}

pub struct VertexTypeIndex<T: IsComplex> {
    tileset: Arc<TileSet<T>>,
    entries: Vec<VertexTypeInfo<T>>,
    transitions: Vec<TransitionInfo>,
    reverse: HashMap<PatchVertexType, usize>,
}

impl<T: IsComplex + IsRingOrField + Units> VertexTypeIndex<T> {
    pub fn new(tileset: Arc<TileSet<T>>) -> Self {
        let mut all_types: BTreeSet<PatchVertexType> = BTreeSet::new();
        let mut kind_map: HashMap<PatchVertexType, VertexTypeKind> = HashMap::new();
        let mut transitions: Vec<(PatchVertexType, PatchVertexType, PatchMatch)> = Vec::new();
        let mut witnesses: HashMap<PatchVertexType, (GrowingPatch<T>, usize)> = HashMap::new();

        let mut visited: BTreeSet<PatchVertexType> = BTreeSet::new();
        let mut queue: VecDeque<(GrowingPatch<T>, usize)> = VecDeque::new();

        for seed_id in 0..tileset.num_tiles() {
            let seed = GrowingPatch::new(Arc::clone(&tileset), seed_id);
            for pm in seed.get_all_matches() {
                let mut gp = seed.clone();
                if gp.add_tile(pm).is_some() && gp.is_growing() {
                    for pos in 0..gp.boundary_len() {
                        if let Some(pvt) = gp.vertex_type_at(pos) {
                            all_types.insert(pvt);
                            witnesses.entry(pvt).or_insert_with(|| (gp.clone(), pos));
                            if visited.insert(pvt) {
                                queue.push_back((gp.clone(), pos));
                            }
                        }
                    }
                }
            }
        }

        while let Some((gp, pos)) = queue.pop_front() {
            let pvt = gp.vertex_type_at(pos).unwrap();
            if kind_map.get(&pvt) == Some(&VertexTypeKind::Dead) {
                continue;
            }

            let touching: Vec<&PatchMatch> = gp.get_matches_touching_vertex(pos).collect();
            if touching.is_empty() {
                kind_map.insert(pvt, VertexTypeKind::Dead);
                continue;
            }

            kind_map.entry(pvt).or_insert(VertexTypeKind::Open);

            for pm in touching {
                let mut gp2 = gp.clone();
                if gp2.add_tile(pm).is_some() && gp2.is_growing() {
                    let old_n = gp.boundary_len();
                    let junction_pos = if pm.start_a == pos { old_n - pm.len } else { 0 };
                    if let Some(new_pvt) = gp2.vertex_type_at(junction_pos) {
                        transitions.push((pvt, new_pvt, pm.clone()));
                    }
                    for new_pos in 0..gp2.boundary_len() {
                        if let Some(nv) = gp2.vertex_type_at(new_pos) {
                            all_types.insert(nv);
                            witnesses
                                .entry(nv)
                                .or_insert_with(|| (gp2.clone(), new_pos));
                            if visited.insert(nv) {
                                queue.push_back((gp2.clone(), new_pos));
                            }
                        }
                    }
                }
            }
        }

        let entries: Vec<PatchVertexType> = all_types.into_iter().collect();
        let reverse: HashMap<PatchVertexType, usize> = entries
            .iter()
            .enumerate()
            .map(|(i, vt)| (*vt, i + 1))
            .collect();

        let mut succ_sets: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); entries.len()];
        let mut pred_sets: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); entries.len()];

        let mut transition_infos: Vec<TransitionInfo> = Vec::new();
        for (src, dst, pm) in &transitions {
            if let (Some(&src_id), Some(&dst_id)) = (reverse.get(src), reverse.get(dst)) {
                succ_sets[src_id - 1].insert(dst_id);
                pred_sets[dst_id - 1].insert(src_id);
                transition_infos.push(TransitionInfo {
                    src_id,
                    dst_id,
                    patch_match: pm.clone(),
                });
            }
        }

        let is_cursed = compute_cursed(&entries, &kind_map, &succ_sets);

        let info_entries: Vec<VertexTypeInfo<T>> = entries
            .into_iter()
            .enumerate()
            .map(|(i, vt)| {
                let id = i + 1;
                let (witness, witness_pos) = witnesses.remove(&vt).expect("witness missing");
                let rat = witness.to_rat();
                VertexTypeInfo {
                    kind: kind_map.get(&vt).copied().unwrap_or(VertexTypeKind::Open),
                    successors: succ_sets[i].iter().copied().collect(),
                    predecessors: pred_sets[i].iter().copied().collect(),
                    is_cursed: is_cursed[&id],
                    realizing_rat: rat,
                    vtype: vt,
                    witness,
                    witness_pos,
                }
            })
            .collect();

        VertexTypeIndex {
            tileset,
            entries: info_entries,
            transitions: transition_infos,
            reverse,
        }
    }

    pub fn num_types(&self) -> usize {
        self.entries.len()
    }

    pub fn transitions(&self) -> &[TransitionInfo] {
        &self.transitions
    }

    pub fn entries(&self) -> &[VertexTypeInfo<T>] {
        &self.entries
    }

    pub fn get_id(&self, vtype: &PatchVertexType) -> Option<usize> {
        self.reverse.get(vtype).copied()
    }

    pub fn get_type(&self, id: usize) -> &PatchVertexType {
        assert!(
            id >= 1 && id <= self.entries.len(),
            "vertex type id out of range"
        );
        &self.entries[id - 1].vtype
    }

    pub fn get_info(&self, id: usize) -> &VertexTypeInfo<T> {
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

fn compute_cursed(
    entries: &[PatchVertexType],
    kind_map: &HashMap<PatchVertexType, VertexTypeKind>,
    succ_sets: &[BTreeSet<usize>],
) -> HashMap<usize, bool> {
    let n = entries.len();
    let mut cursed: HashMap<usize, bool> = HashMap::with_capacity(n);

    for (i, vt) in entries.iter().enumerate() {
        let id = i + 1;
        cursed.insert(id, kind_map.get(vt) == Some(&VertexTypeKind::Dead));
    }

    let mut changed = true;
    while changed {
        changed = false;
        for (i, _vt) in entries.iter().enumerate() {
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
    use crate::cyclotomic::{ZZ12, ZZ4};
    use crate::intgeom::tiles;
    use crate::intgeom::tileset::TileSet;
    use std::sync::Arc;

    #[test]
    fn hexagon_vertex_type_counts() {
        let hex: crate::intgeom::snake::Snake<ZZ12> = tiles::hexagon();
        let rat = Rat::try_from(&hex).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = VertexTypeIndex::new(ts);
        assert!(idx.num_types() > 0, "should discover some vertex types");
        let open = idx
            .entries
            .iter()
            .filter(|e| e.kind == VertexTypeKind::Open)
            .count();
        let dead = idx
            .entries
            .iter()
            .filter(|e| e.kind == VertexTypeKind::Dead)
            .count();
        assert!(open > 0, "should have open types");
        eprintln!(
            "hex: {} types ({} open, {} dead)",
            idx.num_types(),
            open,
            dead
        );
    }

    #[test]
    fn square_vertex_type_counts() {
        let sq: crate::intgeom::snake::Snake<ZZ4> = tiles::square();
        let rat = Rat::try_from(&sq).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = VertexTypeIndex::new(ts);
        assert!(idx.num_types() > 0);
        let open = idx
            .entries
            .iter()
            .filter(|e| e.kind == VertexTypeKind::Open)
            .count();
        eprintln!("square: {} types ({} open)", idx.num_types(), open);
    }

    #[test]
    fn hexagon_no_duplicate_witnesses() {
        let hex: crate::intgeom::snake::Snake<ZZ12> = tiles::hexagon();
        let rat = Rat::try_from(&hex).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = VertexTypeIndex::new(ts);
        for info in &idx.entries {
            let pvt = info.witness().vertex_type_at(info.witness_pos());
            assert_eq!(
                pvt,
                Some(info.vtype),
                "witness should realize its claimed type"
            );
        }
    }
}
