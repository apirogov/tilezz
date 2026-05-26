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
//! - [`patch_grow`]: layer-BFS patch enumeration over an arbitrary
//!   tileset (any combination of tile shapes), returning all
//!   distinct patches up to a given size.

pub mod patch_grow;
pub mod profile;
pub mod seq_explorer;
