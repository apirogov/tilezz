//! Full cyclotomic fields (allowing division)
use std::fmt;
use std::fmt::Display;
use std::ops::{Add, Div, Mul, Neg, Sub};

use num_complex::Complex64;
use num_traits::{One, Zero};

use crate::traits::{Ccw, Conj, InnerIntType, IntRing};
use crate::zz::{HasZZ10, HasZZ12, HasZZ4, HasZZ6, HasZZ8};
use crate::zz::{Z10, Z12, Z4, Z6, Z8, ZZ10, ZZ12, ZZ4, ZZ6, ZZ8};
use crate::zzbase::{ZCommon, ZNum, ZZBase, ZZComplex, ZZNum, ZZParams};
use crate::zzsigned::ZSigned;

/// Compute inverse by repeated rationalizing of the denominator.
///
/// IMPORTANT: This only works correctly if there are no nested roots in
/// the representation of the cyclotomic field
// FIXME
#[inline]
pub fn zz_inv<Z: ZCommon>(val: &Z) -> Z {
    // for x/y where y = a + b, b being a single (scaled) square root,
    // we compute y' = a - b and produce x*y'/(a+b)(a-b) = x*y'/a^2-b^2
    // where a^2 is a simpler term and b is rational.
    // we repeat this for all square roots in the denominator.

    // println!("------------------------\n{val}");

    let num_terms = Z::zz_params().sym_roots_num;

    let mut numer = Z::one();
    let mut denom = val.clone();

    // first ensure that the denominator is real
    let denom_conj = denom.conj();
    numer = numer * denom_conj;
    denom = denom * denom_conj;

    let mut root_ix = 0;
    let mut non_root_ix = 0;
    let mut root_found = true;
    while root_found {
        // println!("\n{numer}\n----\n{denom}\n");

        root_found = false;
        for i in 0..num_terms {
            if Z::zz_params().sym_roots_sqs[i].is_one() {
                non_root_ix = i; // we need the non-irrational part later
                continue; // non-irrational term
            }
            let c = denom.zz_coeffs()[i];
            if !c.is_zero() {
                root_found = true;
                root_ix = i;
                break;
            }
        }

        if !root_found {
            // println!("NO MORE ROOTS");
            break; // rational denominator
        }
        // println!(
        //     "ROOT FOUND: (#{root_ix}), coeff: {}",
        //     denom.zz_coeffs()[root_ix]
        // );

        // compute y' = a - b
        let denom_conj = {
            let curr_root_coeff = denom.zz_coeffs()[root_ix];
            let mut yy = denom.clone();
            yy.zz_coeffs_mut()[root_ix] = -curr_root_coeff;
            yy
        };
        // println!("CONJUGATED DENOM: {}", denom_conj);

        // update numerator (= x * y')
        numer = numer * denom_conj;
        // update denominator (= (a + b)(a - b) = a^2 - b^2)
        denom = denom * denom_conj;
    }

    // now we have a rational denominator (i.e. no square root terms)

    // just flip it and multiply with the numerator to get the final result
    let mut inv_denom = Z::zero();
    inv_denom.zz_coeffs_mut()[non_root_ix] = Z::Scalar::one() / denom.zz_coeffs()[non_root_ix];
    let inv_denom = inv_denom;

    return numer * inv_denom;
}

/// cyclotomic field or its real part
pub trait QCommon {}
/// real part of cyclotomic field
pub trait QNum {}
/// cyclotomic field
pub trait QQNum {}

macro_rules! qq_impl_struct {
    ($q:ident, $z:ident) => {
        #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
        pub struct $q($z);

        impl QCommon for $q {}

        impl From<$z> for $q {
            /// Lift ring value into the corresponding cyclomatic field
            fn from(value: $z) -> Self {
                Self(value)
            }
        }

        impl Add<$q> for $q {
            type Output = Self;
            fn add(self, other: Self) -> Self {
                Self(self.0.add(other.0))
            }
        }
        impl Sub<$q> for $q {
            type Output = Self;
            fn sub(self, other: Self) -> Self {
                Self(self.0.sub(other.0))
            }
        }
        impl Neg for $q {
            type Output = Self;
            fn neg(self) -> Self {
                Self(self.0.neg())
            }
        }
        impl Mul<$q> for $q {
            type Output = Self;
            fn mul(self, other: Self) -> Self {
                Self(self.0.mul(other.0))
            }
        }

        // NOW THIS IS WHY THIS TYPE EXISTS
        // --------
        impl Div<$q> for $q {
            type Output = Self;
            fn div(self, other: Self) -> Self {
                Self(self.0.mul(zz_inv(&other.0)))
            }
        }
        // --------

        impl Zero for $q {
            fn zero() -> Self {
                Self($z::zero())
            }
            fn is_zero(&self) -> bool {
                self.0.is_zero()
            }
        }
        impl One for $q {
            fn one() -> Self {
                Self($z::one())
            }
            fn is_one(&self) -> bool {
                self.0.is_one()
            }
        }
        impl Conj for $q {
            fn conj(&self) -> Self {
                Self(self.0.conj())
            }
            fn co_conj(&self) -> Self {
                Self(self.0.co_conj())
            }
        }

        impl InnerIntType for $q {
            type IntType = <$z as InnerIntType>::IntType;
        }
        impl IntRing for $q {}

        impl Display for $q {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                self.0.fmt(f)
            }
        }
    };
}

macro_rules! qq_impl_signed {
    ($($q:ty)*) => ($(
        impl ZSigned for $q {
            fn signum(&self) -> Self {
                Self(self.0.signum())
            }
        }
    )*)
}

macro_rules! qq_impl_qnum {
    ($($q:ty)*) => ($(
        impl ZNum for $q {}
        impl QNum for $q {}
    )*)
}

macro_rules! qq_impl_ccw {
    ($qq:ident, $zz:ident) => {
        impl Ccw for $qq {
            fn ccw() -> Self {
                Self($zz::ccw())
            }
            fn is_ccw(&self) -> bool {
                *self == Self::ccw()
            }
        }
    };
}

macro_rules! qq_impl_qqnum {
    ($($q:ty)*) => ($(
        impl ZZNum for $q {}
        impl QQNum for $q {}
    )*)
}

macro_rules! qq_impl_complex {
    ($q:ident, $q_real:ident) => {
        impl ZZComplex for $q {
            #[inline]
            fn is_real(&self) -> bool {
                self.0.is_real()
            }
            #[inline]
            fn is_imag(&self) -> bool {
                self.0.is_imag()
            }
            #[inline]
            fn re(&self) -> <Self as ZZBase>::Real {
                $q_real(self.0.re())
            }
            #[inline]
            fn im(&self) -> <Self as ZZBase>::Real {
                $q_real(self.0.im())
            }
        }
    };
}

macro_rules! qq_impl_common {
    ($q:ident, $q_real: ident, $z:ident) => {
        impl ZZBase for $q {
            type Scalar = <$z as ZZBase>::Scalar;
            type Real = $q_real;

            #[inline]
            fn new(coeffs: &[Self::Scalar]) -> Self {
                Self($z::new(coeffs))
            }
            #[inline]
            fn complex64(&self) -> Complex64 {
                self.0.complex64()
            }
            #[inline]
            fn zz_coeffs(&self) -> &[Self::Scalar] {
                self.0.zz_coeffs()
            }
            #[inline]
            fn zz_coeffs_mut(&mut self) -> &mut [Self::Scalar] {
                self.0.zz_coeffs_mut()
            }
            #[inline]
            fn zz_params() -> &'static ZZParams<'static> {
                $z::zz_params()
            }
            #[inline]
            fn zz_mul_arrays(x: &[Self::Scalar], y: &[Self::Scalar]) -> Vec<Self::Scalar> {
                $z::zz_mul_arrays(x, y)
            }
            #[inline]
            fn zz_mul_scalar(x: &[Self::Scalar], scalar: i64) -> Vec<Self::Scalar> {
                $z::zz_mul_scalar(x, scalar)
            }
        }

        impl ZCommon for $q {}
    };
}

qq_impl_struct!(Q4, Z4);
qq_impl_struct!(Q6, Z6);
qq_impl_struct!(Q8, Z8);
qq_impl_struct!(Q10, Z10);
qq_impl_struct!(Q12, Z12);
qq_impl_signed!(Q4 Q6 Q8 Q10 Q12);
qq_impl_common!(Q4, Q4, Z4);
qq_impl_common!(Q6, Q6, Z6);
qq_impl_common!(Q8, Q8, Z8);
qq_impl_common!(Q10, Q10, Z10);
qq_impl_common!(Q12, Q12, Z12);
qq_impl_qnum!(Q4 Q6 Q8 Q10 Q12);

qq_impl_struct!(QQ4, ZZ4);
qq_impl_struct!(QQ6, ZZ6);
qq_impl_struct!(QQ8, ZZ8);
qq_impl_struct!(QQ10, ZZ10);
qq_impl_struct!(QQ12, ZZ12);
qq_impl_common!(QQ4, Q4, ZZ4);
qq_impl_common!(QQ6, Q6, ZZ6);
qq_impl_common!(QQ8, Q8, ZZ8);
qq_impl_common!(QQ10, Q10, ZZ10);
qq_impl_common!(QQ12, Q12, ZZ12);
qq_impl_complex!(QQ4, Q4);
qq_impl_complex!(QQ6, Q6);
qq_impl_complex!(QQ8, Q8);
qq_impl_complex!(QQ10, Q10);
qq_impl_complex!(QQ12, Q12);
qq_impl_ccw!(QQ4, ZZ4);
qq_impl_ccw!(QQ6, ZZ6);
qq_impl_ccw!(QQ8, ZZ8);
qq_impl_ccw!(QQ10, ZZ10);
qq_impl_ccw!(QQ12, ZZ12);
qq_impl_qqnum!(QQ4 QQ6 QQ8 QQ10 QQ12);

pub trait FieldRingPair {
    type Field;
    type Ring;
}

macro_rules! impl_ring_field_typealias {
    ($q:ident, $z:ident) => {
        impl FieldRingPair for $q {
            type Field = $q;
            type Ring = $z;
        }
        impl FieldRingPair for $z {
            type Field = $q;
            type Ring = $z;
        }
    };
}
impl_ring_field_typealias!(Q4, Z4);
impl_ring_field_typealias!(Q6, Z6);
impl_ring_field_typealias!(Q8, Z8);
impl_ring_field_typealias!(Q10, Z10);
impl_ring_field_typealias!(Q12, Z12);

impl_ring_field_typealias!(QQ4, ZZ4);
impl_ring_field_typealias!(QQ6, ZZ6);
impl_ring_field_typealias!(QQ8, ZZ8);
impl_ring_field_typealias!(QQ10, ZZ10);
impl_ring_field_typealias!(QQ12, ZZ12);

impl<T: QCommon + FieldRingPair> HasZZ4 for T where T: HasZZ4 {}
impl<T: QCommon + FieldRingPair> HasZZ6 for T where T: HasZZ6 {}
impl<T: QCommon + FieldRingPair> HasZZ8 for T where T: HasZZ8 {}
impl<T: QCommon + FieldRingPair> HasZZ10 for T where T: HasZZ10 {}
impl<T: QCommon + FieldRingPair> HasZZ12 for T where T: HasZZ12 {}

#[cfg(test)]
mod tests {
    macro_rules! zz_div_tests {
    ($($name:ident: $type:ty,)*) => {$(
mod $name {
    use super::*;
    use crate::zzbase::ZZBase;

    type QQi = <$type as FieldRingPair>::Field;
    type ZZi = <$type as FieldRingPair>::Ring;

    #[test]
    fn test_inv() {
        // try to invert combinations of units
        for k in 0..ZZi::turn() {
            for l in 0..ZZi::turn() {
                let val = ZZi::unit(k) + ZZi::unit(l);
                if val.is_zero() {
                    continue;
                }
                let inv = zz_inv(&val);
                assert_eq!(val * inv, ZZi::one());
            }
        }
    }

    #[test]
    fn test_div() {
        use crate::zz::constants::zz_units_sum;

        // test with worst case - all roots are there
        let z = zz_units_sum::<QQi>();

        // symmetry of dividing by rotation unit
        // (and testing some operations via QQi types)
        assert_eq!(z / ZZi::ccw().into(), z * QQi::unit(-1));
    }
}
    )*}
}

    use super::*;

    zz_div_tests! {
        zz4: ZZ4,
        zz6: ZZ6,
        zz8: ZZ8,
        zz10: ZZ10,
        zz12: ZZ12,
        // zz16: ZZ16,
        // zz20: ZZ20,
        // zz24: ZZ24,
        // zz30: ZZ30,
        // zz32: ZZ32,
        // zz60: ZZ60,
    }
}
