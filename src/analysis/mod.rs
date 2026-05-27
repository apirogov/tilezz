//! Analysis layer: enumeration and classification of structure built on
//! top of [`crate::geom`] primitives.
//!
//! Where `geom` provides the data shapes for tiles, patches, and
//! boundary glues (`Rat`, `Snake`, `TileSet`, `GrowingPatch`,
//! `glue_raw_angles`), the analysis layer asks **questions** about those
//! shapes -- enumerating legal glues, classifying vertex configurations,
//! finding fixed points of boundary subsequences, growing patches up to
//! a maximum size.
//!
//! Modules:
//!
//! - [`matchfinder`]: legal-glue enumeration (`MatchFinder`, `BpSeed`)
//!   -- the operational layer that bridges raw cyclic RC matching
//!   (`stringmatch::BitParallelMatcher`) and geometric validation
//!   (Snake / glue / junction checks).
//! - [`matchtypes`]: catalog layer above `MatchFinder`
//!   (`MatchTypeIndex`) -- pre-computes the full match enumeration
//!   for a fixed tileset and indexes it for O(1) lookup.
//! - [`vertextypes`]: BFS over open-vertex configurations of a tileset
//!   (`OpenVertexTypeIndex`, `Collection`).
//! - [`neighborhood`]: corona / phase-2 classification of local tile
//!   neighborhoods (`NeighborhoodIndex`, `NtKind`, `Collection`).
//! - [`patch_enum`]: layer-BFS enumeration of all distinct patches of a
//!   tileset up to a given size.
//! - [`seq_explorer`]: fixed-point enumeration of cyclic boundary
//!   subsequences reachable from a tileset.

pub mod matchfinder;
pub mod matchtypes;
pub mod neighborhood;
pub mod patch_enum;
pub mod seq_explorer;
pub mod vertextypes;
