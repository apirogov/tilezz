//! Cyclic reverse-complementary substring matching, plus the DAFSA
//! storage stack used to ship large enumerated rat sets.
//!
//! Two independent pieces share this module:
//!
//! 1. **String-matching primitives.** [`BitParallelMatcher`] is the
//!    bulk-matching engine: a streaming Shift-And-style matcher tuned
//!    for the small-alphabet cyclic-RC use case (see its module
//!    docs). The smaller primitives serve different niches:
//!    * [`match_length`] / [`forward_match_length`] (in `extend.rs`)
//!      extend a single match outward from a known anchor pair.
//!      `match_length` is the generic two-slice matcher used by
//!      `Rat::get_match`; `forward_match_length` is a zero-alloc
//!      anti-parallel cyclic specialization for the patch /
//!      neighborhood paths.
//!    * `period::repetition_factor` and `cyclic::cyclic_contains` are
//!      standalone KMP-based utilities -- unrelated to
//!      `BitParallelMatcher` but exposed here because they belong to
//!      the same family of primitives.
//!
//! 2. **DAFSA storage.** See [`dafsa`] for the three-layer stack:
//!    [`Dafsa`] (the automaton), [`RatDafsa`] (length-prefixed rat
//!    wrapper, the single-file format), and [`LazyRatDafsa`] /
//!    [`LazyRatDafsaAsync`] (block-served readers for incremental
//!    loading from disk / HTTP).

mod bitparallel;
mod cyclic;
pub mod dafsa;
mod extend;
mod period;

pub use bitparallel::BitParallelMatcher;
pub use cyclic::cyclic_contains;
pub use dafsa::{Dafsa, LazyRatDafsa, LazyRatDafsaAsync, RatDafsa};
pub use extend::{forward_match_length, match_length};
pub use period::repetition_factor;

// A maximal reverse-complementary match between two cyclic angle sequences
// is the crate-wide `TileMatch` type; re-exported here for ergonomic
// `stringmatch::TileMatch` access at the BitParallelMatcher API boundary.
pub use crate::geom::matches::TileMatch;
