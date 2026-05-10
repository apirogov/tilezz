use std::sync::Arc;

use crate::cyclotomic::{IsComplex, IsRingOrField, Units};
use crate::intgeom::patch::VertexType;
use crate::intgeom::tileset::TileSet;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct NeighborhoodType {
    pub central_tile_id: usize,
    pub cw_anchor_on_central: usize,
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
    pub fn new(_tileset: Arc<TileSet<T>>) -> Self {
        todo!("NT BFS enumeration")
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
