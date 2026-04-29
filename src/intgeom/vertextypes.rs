use std::collections::{BTreeSet, HashMap, VecDeque};
use std::sync::Arc;

use crate::cyclotomic::{IsComplex, IsRingOrField, Units};
use crate::intgeom::patch::{
    cyclic_range_contains, GrowingPatch, PatchMatch, TransitionSide, VertexType,
};
use crate::intgeom::rat::Rat;
use crate::intgeom::tileset::TileSet;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
enum HasTransitions {
    Yes,
    No,
}

pub struct VertexTypeInfo<T: IsComplex> {
    vtype: VertexType,
    has_transitions: HasTransitions,
    successors: Vec<usize>,
    predecessors: Vec<usize>,
    is_cursed: bool,
    realizing_rat: Rat<T>,
    witness: GrowingPatch<T>,
    witness_pos: usize,
    gap_angle: i8,
}

impl<T: IsComplex> VertexTypeInfo<T> {
    pub fn vtype(&self) -> &VertexType {
        &self.vtype
    }

    pub fn is_alive(&self) -> bool {
        !self.is_cursed
    }

    pub fn is_dead(&self) -> bool {
        matches!(self.has_transitions, HasTransitions::No)
    }

    pub fn is_cursed(&self) -> bool {
        self.is_cursed
    }

    pub fn realizing_rat(&self) -> &Rat<T> {
        &self.realizing_rat
    }

    pub fn gap_angle(&self) -> i8 {
        self.gap_angle
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
    pub side: TransitionSide,
    pub patch_match: PatchMatch,
}

pub struct VertexTypeIndex<T: IsComplex> {
    tileset: Arc<TileSet<T>>,
    entries: Vec<VertexTypeInfo<T>>,
    transitions: Vec<TransitionInfo>,
    reverse: HashMap<VertexType, usize>,
}

impl<T: IsComplex + IsRingOrField + Units> VertexTypeIndex<T> {
    pub fn new(tileset: Arc<TileSet<T>>) -> Self {
        let mut all_types: BTreeSet<VertexType> = BTreeSet::new();
        let mut transition_map: HashMap<VertexType, HasTransitions> = HashMap::new();
        let mut raw_transitions: Vec<(VertexType, VertexType, TransitionSide, PatchMatch)> =
            Vec::new();

        let mut visited: BTreeSet<VertexType> = BTreeSet::new();
        let mut queue: VecDeque<VertexType> = VecDeque::new();
        let mut witness_store: HashMap<VertexType, (GrowingPatch<T>, usize, i8)> = HashMap::new();

        for seed_id in 0..tileset.num_tiles() {
            let seed = GrowingPatch::new(Arc::clone(&tileset), seed_id);
            for pm in seed.get_all_matches() {
                let mut gp = seed.clone_for_mutation();
                if gp.add_tile(&pm).is_none() || !gp.is_growing() {
                    continue;
                }
                for pos in 0..gp.boundary_len() {
                    if let Some(vt) = gp.full_vertex_type_at(pos) {
                        all_types.insert(vt.clone());
                        witness_store
                            .entry(vt.clone())
                            .or_insert_with(|| (gp.clone(), pos, gp.angles()[pos]));
                        if visited.insert(vt.clone()) {
                            queue.push_back(vt);
                        }
                    }
                }
            }
        }

        while let Some(vt) = queue.pop_front() {
            if transition_map.get(&vt) == Some(&HasTransitions::No) {
                continue;
            }

            let touching = {
                let (gp, pos, _gap) = &witness_store[&vt];
                let n = gp.boundary_len();
                let t: Vec<PatchMatch> = gp.get_matches_touching_vertex(*pos).cloned().collect();
                if t.is_empty() {
                    transition_map.insert(vt, HasTransitions::No);
                    continue;
                }
                transition_map
                    .entry(vt.clone())
                    .or_insert(HasTransitions::Yes);
                (t, n)
            };
            let (touching, n) = touching;
            let pos = witness_store[&vt].1;

            for pm in touching {
                let gp = &witness_store[&vt].0;
                let mut gp2 = gp.clone_for_mutation();
                if gp2.add_tile(&pm).is_none() || !gp2.is_growing() {
                    continue;
                }

                let covers_ccw = cyclic_range_contains(pm.start_a, pm.len, pos, n);
                let covers_cw = cyclic_range_contains(pm.start_a, pm.len, (pos + n - 1) % n, n);

                let junction_pos = if pm.start_a == pos { n - pm.len } else { 0 };

                if let Some(new_vt) = gp2.full_vertex_type_at(junction_pos) {
                    let side = if covers_cw && !covers_ccw {
                        TransitionSide::Cw
                    } else {
                        TransitionSide::Ccw
                    };
                    raw_transitions.push((vt.clone(), new_vt, side, pm));
                }

                for new_pos in 0..gp2.boundary_len() {
                    if let Some(nv) = gp2.full_vertex_type_at(new_pos) {
                        all_types.insert(nv.clone());
                        if visited.insert(nv.clone()) {
                            witness_store
                                .insert(nv.clone(), (gp2.clone(), new_pos, gp2.angles()[new_pos]));
                            queue.push_back(nv);
                        }
                    }
                }
            }
        }

        let entries: Vec<VertexType> = all_types.into_iter().collect();
        let reverse: HashMap<VertexType, usize> = entries
            .iter()
            .enumerate()
            .map(|(i, vt)| (vt.clone(), i + 1))
            .collect();

        let mut succ_sets: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); entries.len()];
        let mut pred_sets: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); entries.len()];

        let mut transition_infos: Vec<TransitionInfo> = Vec::new();
        for (src, dst, side, pm) in &raw_transitions {
            if let (Some(&src_id), Some(&dst_id)) = (reverse.get(src), reverse.get(dst)) {
                succ_sets[src_id - 1].insert(dst_id);
                pred_sets[dst_id - 1].insert(src_id);
                transition_infos.push(TransitionInfo {
                    src_id,
                    dst_id,
                    side: *side,
                    patch_match: pm.clone(),
                });
            }
        }

        let is_cursed = compute_cursed(&entries, &transition_map, &succ_sets);

        let info_entries: Vec<VertexTypeInfo<T>> = entries
            .into_iter()
            .enumerate()
            .map(|(i, vt)| {
                let id = i + 1;
                let (witness, witness_pos, gap_angle) = witness_store.remove(&vt).unwrap();
                let rat = witness.to_rat();
                VertexTypeInfo {
                    has_transitions: transition_map
                        .get(&vt)
                        .copied()
                        .unwrap_or(HasTransitions::Yes),
                    successors: succ_sets[i].iter().copied().collect(),
                    predecessors: pred_sets[i].iter().copied().collect(),
                    is_cursed: is_cursed[&id],
                    realizing_rat: rat,
                    vtype: vt,
                    witness,
                    witness_pos,
                    gap_angle,
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
    entries: &[VertexType],
    transition_map: &HashMap<VertexType, HasTransitions>,
    succ_sets: &[BTreeSet<usize>],
) -> HashMap<usize, bool> {
    let n = entries.len();
    let mut cursed: HashMap<usize, bool> = HashMap::with_capacity(n);

    for (i, vt) in entries.iter().enumerate() {
        let id = i + 1;
        cursed.insert(id, transition_map.get(vt) == Some(&HasTransitions::No));
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
        let alive = idx.entries.iter().filter(|e| e.is_alive()).count();
        let dead = idx.entries.iter().filter(|e| e.is_dead()).count();
        let cursed = idx.entries.iter().filter(|e| e.is_cursed).count();
        eprintln!(
            "hex: {} types ({} alive, {} dead, {} cursed)",
            idx.num_types(),
            alive,
            dead,
            cursed
        );
    }

    #[test]
    fn square_vertex_type_counts() {
        let sq: crate::intgeom::snake::Snake<ZZ4> = tiles::square();
        let rat = Rat::try_from(&sq).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = VertexTypeIndex::new(ts);
        assert!(idx.num_types() > 0);
        let alive = idx.entries.iter().filter(|e| e.is_alive()).count();
        eprintln!("square: {} types ({} alive)", idx.num_types(), alive);
    }

    #[test]
    fn hexagon_witness_consistency() {
        let hex: crate::intgeom::snake::Snake<ZZ12> = tiles::hexagon();
        let rat = Rat::try_from(&hex).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = VertexTypeIndex::new(ts);
        for info in &idx.entries {
            let vt = info.witness().full_vertex_type_at(info.witness_pos());
            assert_eq!(
                vt,
                Some(info.vtype.clone()),
                "witness should realize its claimed type"
            );
        }
    }
}
