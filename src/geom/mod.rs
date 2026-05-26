//! Data structures for intrinsic, integer-based geometry.
//!
//! Tiles, patches, snakes, rats, glue mechanics, and the angle / edge
//! primitives they're built on. Higher-level enumeration / classification
//! over these shapes lives in [`crate::analysis`].

pub mod angles;

pub mod grid;

pub mod snake;

pub mod rat;

pub mod tiles;

pub mod tileset;

pub mod matches;

pub mod glue;

pub mod patch;
