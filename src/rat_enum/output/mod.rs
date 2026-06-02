//! Post-enumeration output helpers.
//!
//! Once the DFS has produced its `Vec<Vec<i8>>` of canonical
//! sequences, the binary's various `--mode` arms turn them into
//! - polylines (for `--mode render` / GIF output),
//! - DAFSA blobs (handled by `tilezz::stringmatch::RatDafsa`),
//! - human-readable per-length statistics (`--stats`).
//!
//! This module collects the polyline + statistics helpers. The
//! DAFSA writers themselves live in `tilezz::stringmatch`; the
//! binary calls them directly in its mode handlers.

pub mod render;
pub mod stats_summary;

pub use render::{polygons, run_rat_enum_polylines};
pub use stats_summary::{is_achiral, print_stats};
