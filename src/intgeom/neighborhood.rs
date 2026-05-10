use std::collections::VecDeque;
use std::sync::Arc;

use rustc_hash::FxHashMap;

use crate::cyclotomic::{IsComplex, IsRingOrField, Units};
use crate::intgeom::matchtypes::MatchTypeIndex;
use crate::intgeom::patch::VertexType;
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
            let specs = [
                (mt.tile_a, mt.start_a, mt.start_b),
                (mt.tile_b, mt.start_b, mt.start_a),
            ];
            for (central_tile_id, cw_anchor_on_central, cw_anchor_on_context) in specs {
                let nt = NeighborhoodType {
                    central_tile_id,
                    cw_anchor_on_central,
                    cw_anchor_on_context,
                    vt_seq: Vec::new(),
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

        eprintln!(
            "  seeds: {} types, {} transitions",
            entries.len(),
            transitions.len(),
        );

        let _ = match_index;

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
            assert!(nt.vt_seq.is_empty(), "seeds have empty vt_seq");
        }
    }

    #[test]
    fn hex_seed_count() {
        let idx = NeighborhoodIndex::new(hex_tileset());
        assert!(idx.num_types() > 0, "expected non-empty seed collection");
        for nt in idx.entries() {
            assert_eq!(nt.central_tile_id, 0, "single-tile tileset");
            assert!(nt.vt_seq.is_empty(), "seeds have empty vt_seq");
        }
    }

    #[test]
    fn seeds_have_no_duplicate_keys() {
        let idx = NeighborhoodIndex::new(hex_tileset());
        let mut keys = std::collections::HashSet::new();
        for nt in idx.entries() {
            assert!(
                keys.insert((
                    nt.central_tile_id,
                    nt.cw_anchor_on_central,
                    nt.cw_anchor_on_context,
                )),
                "duplicate seed key: central={}, anchor={}, ctx_anchor={}",
                nt.central_tile_id,
                nt.cw_anchor_on_central,
                nt.cw_anchor_on_context,
            );
        }
    }
}
