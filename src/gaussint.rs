use crate::traits::{Ccw, InnerIntType, IntField, IntRing, ZSigned};
use num_integer::Integer;
use num_rational::{Ratio, Rational32, Rational64};
use num_traits::{FromPrimitive, One, Zero};
use std::fmt;
use std::fmt::{Debug, Display};
use std::marker::Copy;
use std::ops::{Add, Div, Mul, Neg, Sub};

/// Gaussian Integer (complex number with real and imaginary part both integers).
///
/// This type mainly exists for bootstrapping, as a scalar coefficient type
/// that is used for constructing other rings and fields.
///
/// NOTE: even though complex numbers cannot really be ordered,
/// for convenience, we derive Ord instance, which is "lexicographic".
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct GaussInt<T> {
    pub real: T,
    pub imag: T,
}

impl<T: IntRing + InnerIntType> InnerIntType for GaussInt<T> {
    type IntType = T::IntType;
}

impl<T: Copy + Neg<Output = T>> GaussInt<T> {
    pub const fn new(re: T, im: T) -> Self {
        Self { real: re, imag: im }
    }

    /// Return complex conjugate (i.e. with negated real value).
    pub fn conj(&self) -> Self {
        Self {
            real: self.real,
            imag: -self.imag,
        }
    }
}

impl<T: IntRing> Zero for GaussInt<T> {
    fn zero() -> Self {
        Self::new(T::zero(), T::zero())
    }

    fn is_zero(&self) -> bool {
        self.real == T::zero() && self.imag == T::zero()
    }
}

impl<T: IntRing> One for GaussInt<T>
where
    GaussInt<T>: Mul<Output = Self>,
{
    fn one() -> Self {
        Self::new(T::one(), T::zero())
    }

    fn is_one(&self) -> bool {
        self.real == T::one() && self.imag == T::zero()
    }
}

impl<T: IntRing + FromPrimitive> FromPrimitive for GaussInt<T> {
    fn from_i64(n: i64) -> Option<Self> {
        T::from_i64(n).and_then(|num| Some(Self::new(num, T::zero())))
    }
    fn from_u64(n: u64) -> Option<Self> {
        T::from_u64(n).and_then(|num| Some(Self::new(num, T::zero())))
    }
}

impl<T: IntRing> Ccw for GaussInt<T> {
    fn ccw() -> Self {
        Self::new(T::zero(), T::one())
    }

    fn is_ccw(&self) -> bool {
        self.real == T::zero() && self.imag == T::one()
    }
}

impl<T: Copy + Neg<Output = T>> Neg for GaussInt<T> {
    type Output = Self;

    fn neg(self) -> Self {
        Self::new(-self.real, -self.imag)
    }
}

impl<T: Add<Output = T>> Add<GaussInt<T>> for GaussInt<T> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            real: self.real + other.real,
            imag: self.imag + other.imag,
        }
    }
}

impl<T: Sub<Output = T>> Sub<GaussInt<T>> for GaussInt<T> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            real: self.real - other.real,
            imag: self.imag - other.imag,
        }
    }
}

impl<T: Copy + Add<Output = T> + Mul<Output = T> + Sub<Output = T>> Mul<GaussInt<T>>
    for GaussInt<T>
{
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        Self {
            real: self.real * other.real - self.imag * other.imag,
            imag: self.real * other.imag + self.imag * other.real,
        }
    }
}

impl<T: IntField> GaussInt<T> {
    /// Return the multiplicative inverse
    pub fn inv(&self) -> Self {
        // 1/(a + b*i) = (a - i*b)/(a^2 + b^2)
        let denom = self.real * self.real + self.imag * self.imag;
        Self {
            real: self.real / denom,
            imag: -self.imag / denom,
        }
    }
}

impl<T: IntField> Div<GaussInt<T>> for GaussInt<T> {
    type Output = Self;

    fn div(self, other: Self) -> Self {
        self * other.inv()
    }
}

impl<T: Display + IntRing + PartialOrd> Display for GaussInt<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut terms: Vec<String> = Vec::new();
        if !self.real.is_zero() {
            terms.push(format!("{}", self.real));
        }
        let neg1: T = -T::one();
        if !self.imag.is_zero() {
            terms.push(if self.imag == T::one() {
                "i".to_string()
            } else if self.imag == neg1 {
                "-i".to_string()
            } else {
                format!("{}i", self.imag)
            });
        }
        return write!(
            f,
            "{}",
            if terms.is_empty() {
                "0".to_string()
            } else {
                terms.join(if terms.len() == 2 && self.imag < T::zero() {
                    ""
                } else {
                    "+"
                })
            }
        );
    }
}

// implement from(...)
impl<T: IntRing> From<T> for GaussInt<T> {
    fn from(value: T) -> Self {
        GaussInt::new(value, T::zero())
    }
}
impl<T: IntRing + Integer> From<T> for GaussInt<Ratio<T>> {
    fn from(value: T) -> Self {
        GaussInt::new(Ratio::new(value, T::one()), Ratio::zero())
    }
}
impl<T: IntRing + Integer> From<(T,)> for GaussInt<Ratio<T>> {
    fn from((value,): (T,)) -> Self {
        GaussInt::from(value)
    }
}
impl<T: IntRing + Integer> From<(T, T)> for GaussInt<Ratio<T>> {
    fn from((re, im): (T, T)) -> Self {
        GaussInt::new(Ratio::new(re, T::one()), Ratio::new(im, T::one()))
    }
}

impl<T: IntRing + ZSigned> ZSigned for GaussInt<T> {
    fn signum(&self) -> Self {
        GaussInt::new(self.real.signum(), T::zero())
    }
}

// for scalars on the left, use singleton tuples as wrapper
impl<T: IntRing + Integer> Add<GaussInt<Ratio<T>>> for (T,) {
    type Output = GaussInt<Ratio<T>>;
    fn add(self, other: GaussInt<Ratio<T>>) -> GaussInt<Ratio<T>> {
        return GaussInt::<Ratio<T>>::from(self) + other;
    }
}
impl<T: IntRing + Integer> Mul<GaussInt<Ratio<T>>> for (T,) {
    type Output = GaussInt<Ratio<T>>;
    fn mul(self, other: GaussInt<Ratio<T>>) -> GaussInt<Ratio<T>> {
        return GaussInt::<Ratio<T>>::from(self) * other;
    }
}
// for scalars on the right, no need for the tuple wrapper
impl<T: IntRing + Integer> Add<T> for GaussInt<Ratio<T>> {
    type Output = GaussInt<Ratio<T>>;
    fn add(self, other: T) -> GaussInt<Ratio<T>> {
        return self + GaussInt::<Ratio<T>>::from(other);
    }
}
impl<T: IntRing + Integer> Mul<T> for GaussInt<Ratio<T>> {
    type Output = GaussInt<Ratio<T>>;
    fn mul(self, other: T) -> GaussInt<Ratio<T>> {
        return self * GaussInt::<Ratio<T>>::from(other);
    }
}
impl<T: IntField + Integer> Div<T> for GaussInt<Ratio<T>> {
    type Output = GaussInt<Ratio<T>>;
    fn div(self, other: T) -> GaussInt<Ratio<T>> {
        return self / GaussInt::<Ratio<T>>::from(other);
    }
}

impl<T: IntRing> IntRing for GaussInt<T> {}
impl<T: IntField> IntField for GaussInt<T> {}

/// The types we will typically want to use.
pub type GaussInt32 = GaussInt<Rational32>;
pub type GaussInt64 = GaussInt<Rational64>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from() {
        assert_eq!(GaussInt32::from(2).real, Rational32::new(2, 1));
        assert_eq!(GaussInt64::from(3).real, Rational64::new(3, 1));
        assert_eq!(GaussInt64::from((4,)).real, Rational64::new(4, 1));
        assert_eq!(GaussInt64::from((5,)).imag, Rational64::new(0, 1));

        let x = GaussInt32::from((3, 5));
        assert_eq!(x.real, Rational32::from(3));
        assert_eq!(x.imag, Rational32::from(5));

        let x_conj = x.conj();
        assert_eq!(x_conj.real, Rational32::from(3));
        assert_eq!(x_conj.imag, Rational32::from(-5));
    }

    #[test]
    fn test_units() {
        type GInt = GaussInt32;
        assert_eq!(GInt::zero() + GInt::zero(), GInt::zero());
        assert_eq!(GInt::one() + GInt::zero(), GInt::one());
        assert_eq!(GInt::zero() + GInt::one(), GInt::one());
        assert_eq!(-GInt::one() + GInt::one(), GInt::zero());
        assert_eq!(GInt::one() - GInt::one(), GInt::zero());

        assert_eq!(GInt::zero() * GInt::zero(), GInt::zero());
        assert_eq!(GInt::one() * GInt::zero(), GInt::zero());
        assert_eq!(GInt::zero() * GInt::one(), GInt::zero());
        assert_eq!(GInt::one() * GInt::one(), GInt::one());
        assert_eq!(-GInt::one() * GInt::one(), -GInt::one());
        assert_eq!(GInt::one() * (-GInt::one()), -GInt::one());
        assert_eq!((-GInt::one()) * (-GInt::one()), GInt::one());

        // NOTE: for Gaussian integers, ccw() is just the complex unit i
        assert_eq!(GInt::ccw() * GInt::ccw(), -GInt::one());
        assert_eq!(-(GInt::ccw()) * -(GInt::ccw()), -GInt::one());
        assert_eq!((-GInt::ccw()) * GInt::ccw(), GInt::one());
        assert_eq!(GInt::one() * GInt::ccw(), GInt::ccw());
        assert_eq!(-(GInt::one()) * GInt::ccw(), -GInt::ccw());
    }

    #[test]
    fn test_ring_ops() {
        type GInt = GaussInt<i32>; // NOTE: this is != GaussInt32
        assert_eq!(-GInt::new(1, 2), GInt::new(-1, -2));
        assert_eq!(GInt::new(1, 2) + GInt::new(3, 4), GInt::new(4, 6));
        assert_eq!(GInt::new(1, 2) - GInt::new(3, 4), GInt::new(-2, -2));
        assert_eq!(GInt::new(1, 2) * GInt::new(3, 4), GInt::new(-5, 10));
    }

    #[test]
    fn test_div() {
        type GInt = GaussInt<Ratio<i32>>;
        let x: GInt = GaussInt::new(2.into(), 3.into());
        let y: GInt = GaussInt::new(5.into(), (-7).into());
        let z: GInt = GaussInt::new((-11, 74).into(), (29, 74).into());
        assert_eq!(x / y, z);
    }

    #[test]
    fn test_display() {
        type GInt = GaussInt64;
        assert_eq!(format!("{}", GInt::zero()), "0");
        assert_eq!(format!("{}", GInt::one()), "1");
        assert_eq!(format!("{}", -GInt::one()), "-1");
        assert_eq!(format!("{}", GInt::ccw()), "i");
        assert_eq!(format!("{}", -GInt::ccw()), "-i");

        assert_eq!(format!("{}", GInt::from((2, 0))), "2");
        assert_eq!(format!("{}", GInt::from((0, 3))), "3i");
        assert_eq!(format!("{}", GInt::from((-4, 0))), "-4");
        assert_eq!(format!("{}", GInt::from((5, -1))), "5-i");
        assert_eq!(format!("{}", GInt::from((6, -2))), "6-2i");
        assert_eq!(format!("{}", GInt::from((-7, -3))), "-7-3i");
        assert_eq!(format!("{}", GInt::from((-8, 9))), "-8+9i");
        assert_eq!(format!("{}", GInt::from((0, -10))), "-10i");
    }
}
