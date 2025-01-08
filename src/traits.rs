use num_integer::Integer;
use num_rational::Ratio;
use num_traits::{One, Zero};
use std::marker::Copy;
use std::ops::{Add, Mul, Neg, Sub};

/// Trait with all traits a reasonable integer ring should provide
pub trait IntRing:
    Zero
    + One
    + Eq
    + PartialEq
    + Neg<Output = Self>
    + Add<Self, Output = Self>
    + Sub<Self, Output = Self>
    + Mul<Self, Output = Self>
    + Copy
    + Clone
{
}
impl IntRing for i32 {}
impl IntRing for i64 {}
impl<T: Integer + IntRing> IntRing for Ratio<T> {}

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

/// A complex integer ring is a cyclotomic ring, i.e.
/// the ring of integers extended by a root of unity.
pub trait ComplexIntRing: IntRing + Ccw {}

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
