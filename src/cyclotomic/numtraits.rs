//! Foundational numeric traits that are *not* cyclotomic-specific.
//!
//! Cyclotomic-ring traits (`SymNum`, `Ccw`, `Conj`, `OneImag`, `ReImSign`,
//! `IntersectUnitSegments`, `WithinRadius`, `Units`, `IsRingOrField`, ...)
//! live in `traits.rs`. This file holds only the building blocks that are
//! sensible for general scalar types (i32, i64, Ratio<T>).

use num_integer::Integer;
use num_rational::Ratio;
use num_traits::{One, Signed, Zero};
use std::marker::Copy;
use std::ops::{Add, Div, Mul, Neg, Sub};

/// Bag of the trait obligations a reasonable integer ring should satisfy.
pub trait IntRing:
    Copy
    + Clone
    + Zero
    + One
    + Neg<Output = Self>
    + Add<Self, Output = Self>
    + Sub<Self, Output = Self>
    + Mul<Self, Output = Self>
    + PartialEq
    + Eq
{
}

/// A field is a ring that supports division.
pub trait IntField: IntRing + Div<Self, Output = Self> {}

impl IntRing for i32 {}
impl IntField for i32 {}
impl IntRing for i64 {}
impl IntField for i64 {}
impl<T: Integer + IntRing> IntRing for Ratio<T> {}
impl<T: Integer + IntField> IntField for Ratio<T> {}

/// Like `num_traits::Signed`, but without the `Num`-like constraints. Lets
/// us define `signum`-style sign-of-an-expression helpers generically over
/// `T: IntRing + ZSigned` without dragging in division-like requirements.
pub trait ZSigned: IntRing {
    fn signum(&self) -> Self;

    fn abs(&self) -> Self {
        *self * self.signum()
    }
    fn is_positive(&self) -> bool {
        self.signum() == Self::one()
    }
    fn is_negative(&self) -> bool {
        self.signum() == -Self::one()
    }
    fn abs_sub(&self, other: &Self) -> Self {
        <Self as ZSigned>::abs(&(*self - *other))
    }
}

impl<T: IntRing + Signed> ZSigned for T {
    fn signum(&self) -> Self {
        self.signum()
    }
}

/// Utility trait to access the integer type on top of which a numeric type
/// is built. Used to bridge concrete-int rings (i32, i64) and their
/// rational extensions.
pub trait InnerIntType {
    type IntType;
}

impl InnerIntType for i32 {
    type IntType = i32;
}
impl InnerIntType for i64 {
    type IntType = i64;
}
impl<T: Integer + IntRing + InnerIntType> InnerIntType for Ratio<T> {
    type IntType = T::IntType;
}
