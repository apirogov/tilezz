use super::traits::{Ccw, InnerIntType, IntRing};
use num_integer::Integer;
use num_rational::{Ratio, Rational32, Rational64};
use num_traits::{One, Zero};
use std::fmt;
use std::fmt::{Debug, Display};
use std::marker::Copy;
use std::ops::{Add, Mul, Neg, Sub};

/// Gaussian Integer (complex number with real and imaginary part both integers).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Hash)]
pub struct GaussInt<T> {
    pub real: T,
    pub imag: T,
}

impl<T: IntRing + InnerIntType> InnerIntType for GaussInt<T> {
    type IntType = T::IntType;
}

impl<T: IntRing> GaussInt<T> {
    pub fn new(re: T, im: T) -> Self {
        Self { real: re, imag: im }
    }

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

impl<T: IntRing> Ccw for GaussInt<T> {
    fn ccw() -> Self {
        Self::new(T::zero(), T::one())
    }

    fn is_ccw(&self) -> bool {
        self.real == T::zero() && self.imag == T::one()
    }
}

impl<T: IntRing> Neg for GaussInt<T> {
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

impl<T: Display + Zero + One + Eq + PartialOrd + Neg<Output = T>> Display for GaussInt<T> {
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

impl<T: Integer + IntRing> IntRing for GaussInt<T> {}

/// Lift an integer into a GaussInt of a rational of the original type.
pub fn to_gint<T: IntRing + Integer>(re: T) -> GaussInt<Ratio<T>> {
    GaussInt::new(Ratio::new(re, T::one()), Ratio::zero())
}

pub fn to_gint2<T: IntRing + Integer>(re: T, im: T) -> GaussInt<Ratio<T>> {
    GaussInt::new(Ratio::new(re, T::one()), Ratio::new(im, T::one()))
}

// convenience impls (so we can write (n,) instead of to_gint(n) to lift n)
impl<T: IntRing + Integer> Add<GaussInt<Ratio<T>>> for (T,) {
    type Output = GaussInt<Ratio<T>>;
    fn add(self, other: GaussInt<Ratio<T>>) -> GaussInt<Ratio<T>> {
        return to_gint(T::from(self.0)) + other;
    }
}
impl<T: IntRing + Integer> Mul<GaussInt<Ratio<T>>> for (T,) {
    type Output = GaussInt<Ratio<T>>;
    fn mul(self, other: GaussInt<Ratio<T>>) -> GaussInt<Ratio<T>> {
        return to_gint(T::from(self.0)) * other;
    }
}

// The types we will typically want to use
pub type GaussInt32 = GaussInt<Rational32>;
pub type GaussInt64 = GaussInt<Rational64>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_to_gint_to_gint2_conj() {
        assert_eq!(to_gint(2).real, Rational32::new(2, 1));
        assert_eq!(to_gint(3).real, Rational64::new(3, 1));
        assert_eq!(to_gint(4).imag, Rational64::new(0, 1));

        let x = to_gint2(3, 5);
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
    fn test_display() {
        type GInt = GaussInt64;
        let x = GInt::zero();
        assert_eq!(format!("{x}"), "0");

        let x = GInt::one();
        assert_eq!(format!("{x}"), "1");
        let x = -GInt::one();
        assert_eq!(format!("{x}"), "-1");
        let x = GInt::ccw();
        assert_eq!(format!("{x}"), "i");
        let x = -GInt::ccw();
        assert_eq!(format!("{x}"), "-i");

        let x = to_gint2(2, 0);
        assert_eq!(format!("{x}"), "2");
        let x = to_gint2(0, 3);
        assert_eq!(format!("{x}"), "3i");
        let x = to_gint2(-4, 0);
        assert_eq!(format!("{x}"), "-4");
        let x = to_gint2(5, -1);
        assert_eq!(format!("{x}"), "5-i");
        let x = to_gint2(6, -2);
        assert_eq!(format!("{x}"), "6-2i");
        let x = to_gint2(-7, -3);
        assert_eq!(format!("{x}"), "-7-3i");
        let x = to_gint2(-8, 9);
        assert_eq!(format!("{x}"), "-8+9i");
        let x = to_gint2(0, -10);
        assert_eq!(format!("{x}"), "-10i");
    }
}
