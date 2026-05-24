//! Standalone patch-enumeration utilities.
//!
//! Independent explorers and growth algorithms that build on the
//! core [`crate::intgeom`] geometry primitives (`Rat`, `Snake`,
//! `TileSet`, `GrowingPatch`, `MatchFinder`) but are not part of
//! the central feature-type / catalog story. Grouped here so they
//! can evolve independently and so consumers of `intgeom` don't
//! pay the cost of compiling them when unwanted.
//!
//! - [`seq_explorer`]: BFS over patches that enumerates every
//!   cyclic boundary subsequence reachable from a tileset.
//! - [`vtseq_explorer`]: similar in spirit to `seq_explorer`, built
//!   directly on `GrowingPatch` instead of `MatchFinder`.
//! - [`redelmeier`]: Redelmeier-style polyomino-like enumeration
//!   over a single polygonal tile.

pub mod redelmeier;
pub mod seq_explorer;
pub mod vtseq_explorer;
