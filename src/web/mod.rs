//! Browser-facing entry points compiled into a WASM bundle by
//! `wasm-pack build --no-default-features --features rat_explorer`.
//!
//! Currently hosts the [`rat_explorer`] static page; future static
//! sites would live as sibling modules.

pub mod rat_explorer;
