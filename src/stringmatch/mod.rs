//! Cyclic reverse-complementary substring matching.
//!
//! The only engine in this module is [`BitParallelMatcher`], a
//! streaming Shift-And-style matcher tuned for the small-alphabet
//! cyclic-RC use case (see its module docs). It precomputes per-tile,
//! per-symbol bitsets once and streams arbitrary boundaries against
//! them.
//!
//! `period::repetition_factor` and `cyclic::cyclic_contains` are
//! standalone string-matching utilities (both built on the same KMP
//! prefix function) — unrelated to `BitParallelMatcher` but exposed
//! here because they belong to the same family of primitives.

mod bitparallel;
mod cyclic;
mod period;

pub use bitparallel::BitParallelMatcher;
pub use cyclic::cyclic_contains;
pub use period::repetition_factor;

// A maximal reverse-complementary match between two cyclic angle sequences
// is the crate-wide `TileMatch` type; re-exported here for ergonomic
// `stringmatch::TileMatch` access at the BitParallelMatcher API boundary.
pub use crate::geom::matches::TileMatch;
