use std::collections::VecDeque;
use std::sync::Arc;

use rustc_hash::FxHashMap;

use crate::cyclotomic::{IsComplex, IsRingOrField, Units};
use crate::intgeom::matchtypes::MatchTypeIndex;
use crate::intgeom::patch::{GrowingPatch, PatchMatch, VertexType};
use crate::intgeom::tileset::TileSet;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct NeighborhoodType {
    pub central_tile_id: usize,
    pub cw_anchor_on_central: usize,
    pub cw_anchor_on_context: usize,
    pub vt_seq: Vec<VertexType>,
}

pub const NT_CLOSED_ID: usize = 0;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum NtKind {
    Dead,
    Undead,
    Blessed,
    Free,
}

#[derive(Clone, Debug)]
pub struct NtTransition {
    pub src_id: usize,
    pub dst_id: usize,
    pub tile_id: usize,
    pub tile_offset: usize,
}

pub struct NeighborhoodIndex<T: IsComplex> {
    tileset: Arc<TileSet<T>>,
    entries: Vec<NeighborhoodType>,
    transitions: Vec<NtTransition>,
}

impl<T: IsComplex + IsRingOrField + Units> NeighborhoodIndex<T> {
    pub fn new(tileset: Arc<TileSet<T>>) -> Self {
        let match_index = Arc::new(MatchTypeIndex::new(Arc::clone(&tileset)));
        let mut entries: Vec<NeighborhoodType> = Vec::new();
        let transitions: Vec<NtTransition> = Vec::new();
        let mut seen: FxHashMap<NeighborhoodType, usize> = FxHashMap::default();
        let mut queue: VecDeque<NeighborhoodType> = VecDeque::new();

        for id in 1..=match_index.num_types() {
            let mt = match_index.get(id);
            let contexts = [
                (mt.tile_a, mt.start_a, mt.tile_b, mt.start_b),
                (mt.tile_b, mt.start_b, mt.tile_a, mt.start_a),
            ];
            for (seed_tile, seed_start, other_tile, other_start) in contexts {
                let mut patch = GrowingPatch::new(Arc::clone(match_index.tileset()), seed_tile);
                let pm = PatchMatch {
                    start_a: seed_start,
                    len: mt.len,
                    start_b: other_start,
                    tile_id: other_tile,
                };
                if patch.add_tile(&pm).is_none() {
                    continue;
                }

                for third_pm in patch.get_all_matches() {
                    let mut three = patch.clone();
                    if three.add_tile(&third_pm).is_none() {
                        continue;
                    }
                    let central_tile_id = third_pm.tile_id;
                    let cw_anchor_on_central = third_pm.start_b;
                    let cw_anchor_on_context = third_pm.start_a;
                    let vt_seq = extract_context_vt_seq(&three, central_tile_id);
                    if vt_seq.is_empty() {
                        continue;
                    }
                    let nt = NeighborhoodType {
                        central_tile_id,
                        cw_anchor_on_central,
                        cw_anchor_on_context,
                        vt_seq,
                    };
                    if seen.contains_key(&nt) {
                        continue;
                    }
                    let idx = entries.len();
                    seen.insert(nt.clone(), idx);
                    entries.push(nt.clone());
                    queue.push_back(nt);
                }
            }
        }

        let mut explored = 0usize;
        while let Some(state) = queue.pop_front() {
            explored += 1;
            let _state_id = seen.get(&state).expect("dequeued state must be in seen");
        }

        eprintln!(
            "  total: {} types, {} transitions, {} explored",
            entries.len(),
            transitions.len(),
            explored,
        );

        NeighborhoodIndex {
            tileset,
            entries,
            transitions,
        }
    }

    pub fn entries(&self) -> &[NeighborhoodType] {
        &self.entries
    }

    pub fn num_types(&self) -> usize {
        self.entries.len()
    }

    pub fn transitions(&self) -> &[NtTransition] {
        &self.transitions
    }

    pub fn tileset(&self) -> &Arc<TileSet<T>> {
        &self.tileset
    }

    pub fn classify_all(&self) -> Vec<NtKind> {
        let n = self.entries.len();
        let mut succ_sets: Vec<Vec<usize>> = vec![vec![]; n];
        let mut has_outgoing = vec![false; n];

        for t in &self.transitions {
            if t.src_id == NT_CLOSED_ID {
                continue;
            }
            let src_idx = t.src_id;
            has_outgoing[src_idx] = true;
            if t.dst_id != NT_CLOSED_ID {
                succ_sets[src_idx].push(t.dst_id);
            }
        }

        let mut cursed = vec![false; n];
        for i in 0..n {
            cursed[i] = !has_outgoing[i];
        }

        let mut changed = true;
        while changed {
            changed = false;
            for i in 0..n {
                if cursed[i] {
                    continue;
                }
                let succs = &succ_sets[i];
                if succs.is_empty() {
                    continue;
                }
                if succs.iter().all(|&s| cursed[s]) {
                    cursed[i] = true;
                    changed = true;
                }
            }
        }

        let mut blessed = vec![false; n];
        let mut has_closed_or_blessed_succ = vec![false; n];

        for t in &self.transitions {
            if t.src_id == NT_CLOSED_ID {
                continue;
            }
            let src_idx = t.src_id;
            if t.dst_id == NT_CLOSED_ID {
                has_closed_or_blessed_succ[src_idx] = true;
            }
        }

        changed = true;
        while changed {
            changed = false;
            for i in 0..n {
                if blessed[i] || cursed[i] {
                    continue;
                }
                let mut any = false;
                let all_good = self.transitions.iter().filter(|t| t.src_id == i).all(|t| {
                    any = true;
                    t.dst_id == NT_CLOSED_ID || blessed[t.dst_id]
                });
                if any && all_good {
                    blessed[i] = true;
                    has_closed_or_blessed_succ[i] = true;
                    changed = true;
                }
            }
        }

        self.entries
            .iter()
            .enumerate()
            .map(|(i, _)| {
                if cursed[i] && !has_outgoing[i] {
                    NtKind::Dead
                } else if cursed[i] {
                    NtKind::Undead
                } else if blessed[i] {
                    NtKind::Blessed
                } else {
                    NtKind::Free
                }
            })
            .collect()
    }

    pub fn validate(&self) -> Vec<usize> {
        todo!("validate")
    }

    pub fn write_collection(&self, _out: &mut impl std::io::Write) -> std::io::Result<()> {
        todo!("write_collection")
    }

    pub fn parse_file(_tileset: Arc<TileSet<T>>, _input: &str) -> Result<Self, String> {
        todo!("parse_file")
    }
}

fn extract_context_vt_seq<T: IsComplex + IsRingOrField + Units>(
    patch: &GrowingPatch<T>,
    central_tile_id: usize,
) -> Vec<VertexType> {
    let edges = patch.edges();
    let n = edges.len();
    if n == 0 {
        return Vec::new();
    }
    let mut vts = Vec::new();
    for (pos, edge) in edges.iter().enumerate() {
        if patch.is_junction(pos) && edge.tile_id == central_tile_id {
            if let Some(vt) = patch.junction_vertex_type_at(pos) {
                vts.push(vt);
            }
        }
    }
    vts
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic::ZZ12;
    use crate::intgeom::rat::Rat;
    use crate::intgeom::tiles;

    fn square_tileset() -> Arc<TileSet<ZZ12>> {
        let sq = tiles::square::<ZZ12>();
        let rat = Rat::try_from(&sq).unwrap();
        Arc::new(TileSet::new(vec![rat]))
    }

    fn hex_tileset() -> Arc<TileSet<ZZ12>> {
        let hex = tiles::hexagon::<ZZ12>();
        let rat = Rat::try_from(&hex).unwrap();
        Arc::new(TileSet::new(vec![rat]))
    }

    #[test]
    fn square_seed_count() {
        let idx = NeighborhoodIndex::new(square_tileset());
        assert!(idx.num_types() > 0, "expected non-empty seed collection");
        for nt in idx.entries() {
            assert_eq!(nt.central_tile_id, 0, "single-tile tileset");
            assert!(!nt.vt_seq.is_empty(), "seeds must have non-empty vt_seq");
        }
    }

    #[test]
    fn hex_seed_count() {
        let idx = NeighborhoodIndex::new(hex_tileset());
        assert!(idx.num_types() > 0, "expected non-empty seed collection");
        for nt in idx.entries() {
            assert_eq!(nt.central_tile_id, 0, "single-tile tileset");
            assert!(!nt.vt_seq.is_empty(), "seeds must have non-empty vt_seq");
        }
    }

    #[test]
    fn seeds_have_no_duplicate_keys() {
        let idx = NeighborhoodIndex::new(hex_tileset());
        let mut keys = std::collections::HashSet::new();
        for (i, nt) in idx.entries().iter().enumerate() {
            assert!(
                keys.insert((
                    nt.central_tile_id,
                    nt.cw_anchor_on_central,
                    nt.cw_anchor_on_context,
                    nt.vt_seq.clone(),
                )),
                "[{}] duplicate seed key: central={}, anchor={}, ctx_anchor={}, vt_seq={:?}",
                i,
                nt.central_tile_id,
                nt.cw_anchor_on_central,
                nt.cw_anchor_on_context,
                nt.vt_seq,
            );
        }
    }

    #[test]
    fn bfs_explores_all_seeds() {
        let sq_idx = NeighborhoodIndex::new(square_tileset());
        assert!(sq_idx.num_types() > 0);
        assert!(sq_idx.transitions().is_empty());

        let hex_idx = NeighborhoodIndex::new(hex_tileset());
        assert!(hex_idx.num_types() > 0);
        assert!(hex_idx.transitions().is_empty());
    }
}
