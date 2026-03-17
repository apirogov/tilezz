//! Fast lookup of roots of unity (unit vectors) per cyclotomic complex type.
//!
//! This is intentionally *not* part of `SymNum` to avoid conflicting impls and
//! to allow per-concrete-type caching via `OnceLock`.

use super::numtraits::Ccw;
use super::traits::{IsComplex, IsRingOrField};
use num_traits::One;

/// Provides cached lookup for unit vectors (roots of unity) for a concrete complex type.
pub trait Units: Copy + One + Ccw + IsComplex + IsRingOrField {
    /// Return unit length vector pointing in direction of given angle.
    fn unit(angle: i8) -> Self;

    /// Build the full unit table of size `turn()`.
    #[inline]
    fn build_unit_table() -> Vec<Self> {
        let n = Self::turn();
        (0..n).map(|k| Self::ccw().zz_pow(k as u8)).collect()
    }
}

macro_rules! impl_units_for {
    ($t:ty) => {
        impl Units for $t {
            #[inline]
            fn unit(angle: i8) -> Self {
                static TABLE: std::sync::OnceLock<Vec<$t>> = std::sync::OnceLock::new();
                let tab = TABLE.get_or_init(<$t as Units>::build_unit_table);
                let idx = angle.rem_euclid(<$t as super::symnum::SymNum>::turn()) as usize;
                tab[idx]
            }
        }
    };
}

pub(crate) use impl_units_for;
