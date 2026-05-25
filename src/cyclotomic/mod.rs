//! Cyclotomic rings and fields.

mod numtraits;

pub mod gaussint;

#[macro_use]
mod symnum;

mod sign;

#[cfg(test)]
mod sign_extra_tests;
#[cfg(test)]
mod unit_tests;

mod params;

mod traits;

#[macro_use]
mod types;

mod units;

// --------

pub mod constants;

pub mod linalg;

pub mod geometry;

pub use num_traits::{One, Pow, Zero};

pub use numtraits::{Ccw, Conj, IntField, IntRing, OneImag};

pub use symnum::{SymNum, ZZComplex};
pub use units::Units;

pub use traits::{ComplexTraits, RealTraits, RingTraits};
pub use traits::{HasZZ10, HasZZ12, HasZZ4, HasZZ6, HasZZ8};
pub use traits::{IsComplex, IsReal, IsRing, IsRingOrField, IsZZ4};
pub use traits::{ZType, ZZType};

pub use types::{Z10, Z12, Z16, Z20, Z24, Z30, Z32, Z4, Z6, Z60, Z8};
pub use types::{ZZ10, ZZ12, ZZ16, ZZ20, ZZ24, ZZ30, ZZ32, ZZ4, ZZ6, ZZ60, ZZ8};
