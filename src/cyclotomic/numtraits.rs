//! General traits that are needed or convenient for technical reasons to simplify implementation.

use num_integer::Integer;
use num_rational::Ratio;
use num_traits::{One, Signed, Zero};
use std::marker::Copy;
use std::ops::{Add, Div, Mul, Neg, Sub};

/// Trait with all traits a reasonable integer ring should provide.
pub trait IntRing:
    Copy
    + Clone
    + Eq
    + PartialEq
    + Neg<Output = Self>
    + Add<Self, Output = Self>
    + Sub<Self, Output = Self>
    + Mul<Self, Output = Self>
    + Zero
    + One
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

// NOTE: this only makes sense for the complex fields.
pub trait Ccw {
    /// Return unit of (counterclockwise) rotation for a complex integer ring,
    /// i.e. unit vector rotated by the smallest possible rotation step.
    ///
    /// Technically, this is the actual single generator of the cyclotomic ring,
    /// as we can recover the normal unit by taking a power of the rotation unit.
    fn ccw() -> Self;

    /// Return true if the value is equal to the rotation unit.
    fn is_ccw(&self) -> bool;
}

// NOTE: Cyclotomic rings are closed under complex conjugation, so this is ok to do.
pub trait Conj {
    /// For representations of complex numbers, return conjugated complex number.
    /// For one-dimensional numbers, return same number.
    fn conj(&self) -> Self;
}

impl<T: Integer + Clone + std::ops::Neg<Output = T>> Conj for Ratio<T> {
    fn conj(&self) -> Self {
        self.clone()
    }
}

/// Trait to be implemented by rings that support the imaginary unit i directly.
// NOTE: Only to be implemented for types that are a superset of the Gaussian integers.
pub trait OneImag {
    /// Return imaginary unit.
    fn one_i() -> Self;

    /// Return true if the value is equal to the imaginary unit.
    fn is_one_i(&self) -> bool;
}

/// Utility trait to access the integer type on top of which
/// more complex numeric types are built.
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

/// ZSigned is like num_traits::Signed, but without the Num-like constraints.
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
