//! Cyclotomic rings and fields.

pub(crate) mod numtraits;

pub mod gaussint;

#[macro_use]
mod symnum;

mod sign;

#[cfg(test)]
mod unit_tests;

mod params;

mod traits;

#[macro_use]
mod types;

mod zz12;

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

pub use types::{ZZ10, ZZ16, ZZ20, ZZ24, ZZ32, ZZ4, ZZ60, ZZ8};
pub use zz12::ZZ12;
