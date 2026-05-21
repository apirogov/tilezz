//! Cyclic reverse-complementary substring matching.
//!
//! The only engine in this module is [`BitParallelMatcher`], a
//! streaming Shift-And-style matcher tuned for the small-alphabet
//! cyclic-RC use case (see its module docs). It precomputes per-tile,
//! per-symbol bitsets once and streams arbitrary boundaries against
//! them.
//!
//! `period::repetition_factor` is unrelated string-matching utility
//! left available because it stands on its own (KMP prefix function).

mod bitparallel;
mod period;

pub use bitparallel::BitParallelMatcher;
pub use period::repetition_factor;

/// A maximal reverse-complementary match between two cyclic angle sequences.
///
/// Represents a shared boundary segment between two polygonal tiles:
/// tile A's subsequence starting at `pos_a` matches tile B's subsequence
/// traced in reverse (CW direction) starting at `pos_b`, for `len` edges.
///
/// The match is maximal: it cannot be extended in either direction cyclically.
pub struct CyclicMatch {
    pub tile_a: usize,
    pub pos_a: usize,
    pub tile_b: usize,
    pub pos_b: usize,
    pub len: usize,
}

impl std::fmt::Debug for CyclicMatch {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "CyclicMatch(tile_a={}, pos_a={}, tile_b={}, pos_b={}, len={})",
            self.tile_a, self.pos_a, self.tile_b, self.pos_b, self.len
        )
    }
}
