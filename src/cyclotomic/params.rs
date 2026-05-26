//! Per-ring `ZZParams` constants (just `full_turn_steps`) plus the `f64`
//! constants used by the various `complex64_fn` projections.
//!
//! All other ring-specific data (reduction rules, decomp tables, real-sign
//! helpers, Cartesian projections) live in `rings.rs` next to each ring's
//! `define_integral_zz!` invocation.

use std::f64::consts::SQRT_2;

use super::symnum::ZZParams;

// ---- f64 nested-radical constants used by the per-ring complex64_fn ----

/// `2 + sqrt(2)` -- used by ZZ16 (and ZZ32 as the inner nested radical).
pub const ZZ16_Y: f64 = 2.0 + SQRT_2;

/// `2 + sqrt(2 + sqrt(2))` -- used by ZZ32. The literal `1.847...` is
/// `sqrt(2 + sqrt(2))` to 16 digits.
pub const ZZ32_Z: f64 = 2.0 + 1.847_759_065_022_573_5;

/// `sqrt(5)` to 16 digits.
pub const SQRT_5: f64 = 2.236_067_977_499_79;

/// `2 * (5 - sqrt(5))` -- used by ZZ10, ZZ20, ZZ60.
pub const ZZ10_Y: f64 = 2.0 * (5.0 - SQRT_5);

// ---- Per-ring `ZZParams` ----

pub const ZZ4_PARAMS: ZZParams = ZZParams { full_turn_steps: 4 };
pub const ZZ8_PARAMS: ZZParams = ZZParams { full_turn_steps: 8 };
pub const ZZ10_PARAMS: ZZParams = ZZParams { full_turn_steps: 10 };
pub const ZZ12_PARAMS: ZZParams = ZZParams { full_turn_steps: 12 };
pub const ZZ16_PARAMS: ZZParams = ZZParams { full_turn_steps: 16 };
pub const ZZ20_PARAMS: ZZParams = ZZParams { full_turn_steps: 20 };
pub const ZZ24_PARAMS: ZZParams = ZZParams { full_turn_steps: 24 };
pub const ZZ32_PARAMS: ZZParams = ZZParams { full_turn_steps: 32 };
pub const ZZ60_PARAMS: ZZParams = ZZParams { full_turn_steps: 60 };
