//! Cyclotomic rings and fields.

mod numtraits;

mod gaussint;

#[macro_use]
mod symnum;

mod sign;

mod div;

mod params;

mod traits;

#[macro_use]
mod types;

// --------

pub mod constants;

pub mod geometry;

pub use num_traits::{One, Zero};

pub use numtraits::{Ccw, Conj, IntField, IntRing, OneImag};

pub use symnum::{SymNum, ZZComplex};

pub use traits::{ComplexTraits, FieldTraits, RealTraits, RingTraits};
pub use traits::{HasZZ10, HasZZ12, HasZZ4, HasZZ6, HasZZ8};
pub use traits::{IsComplex, IsField, IsReal, IsRealOrComplex, IsRing, IsRingOrField};
pub use traits::{QQType, QType, ZType, ZZType};

pub use types::{Z10, Z12, Z16, Z20, Z24, Z30, Z32, Z4, Z6, Z60, Z8};
pub use types::{ZZ10, ZZ12, ZZ16, ZZ20, ZZ24, ZZ30, ZZ32, ZZ4, ZZ6, ZZ60, ZZ8};

// NOTE: division does not work for the other fields yet (more complicated case)
pub use types::{Q10, Q12, Q24, Q4, Q6, Q8};
pub use types::{QQ10, QQ12, QQ24, QQ4, QQ6, QQ8};
