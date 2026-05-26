//! Cyclotomic rings.
//!
//! Every `ZZn` ring lives on the **integer-basis** representation in
//! `rings.rs`: storage `[i64; phi(n)]` over the power basis
//! `{1, zeta, ..., zeta^(phi(n)-1)}`, with closed-form exact sign
//! extraction (no f64 fallback). The shared engine lives in
//! `integral_basis.rs` (free functions + `define_integral_zz!` macro).

pub(crate) mod numtraits;

pub mod gaussint;

pub(crate) mod symnum;

mod sign;

#[cfg(test)]
mod unit_tests;

pub(crate) mod params;

pub(crate) mod traits;

mod rings;

mod units;

pub mod integral_basis;

// --------

pub mod constants;

pub mod linalg;

pub mod geometry;

pub use num_traits::{One, Pow, Zero};

pub use numtraits::{
    Ccw, Conj, IntField, IntRing, IntersectUnitSegments, OneImag, ReImSign, WithinRadius,
};

pub use symnum::{SymNum, ZZComplex};
pub use units::Units;

pub use traits::{ComplexTraits, RingTraits};
pub use traits::{HasZZ10, HasZZ12, HasZZ4, HasZZ6, HasZZ8};
pub use traits::{IsRing, IsRingOrField, IsZZ4};
pub use traits::ZZType;

pub use rings::{ZZ10, ZZ12, ZZ16, ZZ20, ZZ24, ZZ32, ZZ4, ZZ60, ZZ8};
