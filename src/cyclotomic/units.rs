//! Cached lookup of roots of unity (unit vectors) per cyclotomic complex type.
//!
//! This is intentionally *not* part of `SymNum` to avoid conflicting impls
//! and to allow per-concrete-type caching. Every concrete ring picks one of
//! two impls in its `define_integral_zz!` block:
//!
//! * `impl_integral_units_via_basis!` -- one `OnceLock<Vec<[i64; PHI]>>`
//!   built from `derive_units_lookup` at first call.
//! * A hand-rolled `static UNIT_TABLE: [[i64; PHI]; N] = [...]` overriding
//!   `unit(angle)` for the rat_enum hot path (ZZ12).

use super::numtraits::Ccw;
use super::traits::IsRingOrField;
use num_traits::One;

/// Provides cached lookup for unit vectors (roots of unity) for a concrete
/// complex type.
pub trait Units: Copy + One + Ccw + IsRingOrField {
    /// Return unit length vector pointing in direction of given angle.
    fn unit(angle: i8) -> Self;

    /// Build the full unit table of size `turn()` by walking `ccw()^k`.
    /// Used by helpers that need an iterator over all roots of unity.
    #[inline]
    fn build_unit_table() -> Vec<Self> {
        let n = Self::turn();
        (0..n).map(|k| Self::ccw().zz_pow(k as u8)).collect()
    }
}
