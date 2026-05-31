//! Cyclic reverse-complementary substring matching.
//!
//! The bulk-matching engine is [`BitParallelMatcher`], a streaming
//! Shift-And-style matcher tuned for the small-alphabet cyclic-RC use
//! case (see its module docs). It precomputes per-tile, per-symbol
//! bitsets once and streams arbitrary boundaries against them.
//!
//! The smaller primitives serve different niches:
//!   * [`match_length`] / [`forward_match_length`] (in `extend.rs`)
//!     extend a single match outward from a known anchor pair.
//!     `match_length` is the generic two-slice matcher used by
//!     `Rat::get_match`; `forward_match_length` is a zero-alloc
//!     anti-parallel cyclic specialization for the patch /
//!     neighborhood paths.
//!   * `period::repetition_factor` and `cyclic::cyclic_contains` are
//!     standalone KMP-based utilities — unrelated to
//!     `BitParallelMatcher` but exposed here because they belong to
//!     the same family of primitives.
//!   * [`Dafsa`] is a minimum-state DAFSA with `from_seqs(sorted)` /
//!     `iter()` round-trip plus `contains` / `get(i)` / `walk_prefix`
//!     queries. Independent of the rest of the module; used as a
//!     compact storage format for large enumerated sequence sets
//!     (e.g. dumping `rat_enum` output for a browser-side explorer).

mod bitparallel;
mod cyclic;
mod dafsa;
mod extend;
mod period;

pub use bitparallel::BitParallelMatcher;
pub use cyclic::cyclic_contains;
pub use dafsa::{
    Dafsa, DafsaBuilder, DafsaCursor, DafsaIter, JSON_SCHEMA_DOC as DAFSA_JSON_SCHEMA_DOC,
};
pub use extend::{forward_match_length, match_length};
pub use period::repetition_factor;

// A maximal reverse-complementary match between two cyclic angle sequences
// is the crate-wide `TileMatch` type; re-exported here for ergonomic
// `stringmatch::TileMatch` access at the BitParallelMatcher API boundary.
pub use crate::geom::matches::TileMatch;
