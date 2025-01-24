//! Full cyclotomic fields (allowing division)
use paste::paste;
use std::fmt;
use std::fmt::Display;
use std::ops::{Add, Div, Mul, Neg, Sub};

use num_complex::Complex64;
use num_integer::Integer;
use num_traits::{FromPrimitive, One, ToPrimitive, Zero};

use crate::traits::{Ccw, Conj, InnerIntType, IntRing, OneImag};
use crate::zz::{HasZZ10, HasZZ12, HasZZ4, HasZZ6, HasZZ8};
use crate::zz::{Z10, Z12, Z24, Z4, Z6, Z8, ZZ10, ZZ12, ZZ24, ZZ4, ZZ6, ZZ8};
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

/// Utility trait to connect related field and ring types.
pub trait FieldRingPair {
    type Field;
    type Ring;
}

impl<T: QCommon + FieldRingPair> HasZZ4 for T where T::Ring: HasZZ4 {}
impl<T: QCommon + FieldRingPair> HasZZ6 for T where T::Ring: HasZZ6 {}
impl<T: QCommon + FieldRingPair> HasZZ8 for T where T::Ring: HasZZ8 {}
impl<T: QCommon + FieldRingPair> HasZZ10 for T where T::Ring: HasZZ10 {}
impl<T: QCommon + FieldRingPair> HasZZ12 for T where T::Ring: HasZZ12 {}

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

        impl From<<$q as InnerIntType>::IntType> for $q {
            /// lift an integer
            fn from(value: <<Self as ZZBase>::Scalar as InnerIntType>::IntType) -> Self {
                Self::one().scale(value)
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

macro_rules! qq_impl_qq_from_q {
    ($q:ident, $q_real: ident) => {
        impl From<$q_real> for $q {
            /// Lift real-valued ring value into the corresponding cyclomatic ring
            fn from(value: $q_real) -> Self {
                let cs: Vec<<$q as ZZBase>::Scalar> =
                    value.zz_coeffs().iter().map(|v| (*v).into()).collect();
                Self::new(cs.as_slice())
            }
        }
    };
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

macro_rules! qq_derive {
    ($($n:expr)*) => ($(
        paste! {
            qq_impl_struct!([<Q $n>], [<Z $n>]);
            qq_impl_common!([<Q $n>], [<Q $n>] , [<Z $n>]);
            qq_impl_signed!([<Q $n>]);
            qq_impl_qnum!([<Q $n>]);

            qq_impl_struct!([<QQ $n>], [<ZZ $n>]);
            qq_impl_common!([<QQ $n>], [<Q $n>] , [<ZZ $n>]);
            qq_impl_complex!([<QQ $n>], [<Q $n>]);
            qq_impl_ccw!([<QQ $n>], [<ZZ $n>]);
            qq_impl_qqnum!([<QQ $n>]);
            qq_impl_qq_from_q!([<QQ $n>], [<Q $n>]);

            impl_ring_field_typealias!([<Q $n>], [<Z $n>]);
            impl_ring_field_typealias!([<QQ $n>], [<ZZ $n>]);
        }
    )*)
}

// generate all the boilerplate code
qq_derive!(4 6 8 10 12 24);

// NOTE: list should be equal to members of HasZZ4
impl_from_int_pair!(QQ4 QQ8 QQ12 /* QQ16 QQ20 */ QQ24 /* QQ32 QQ60 */);

#[cfg(test)]
mod tests {
    use super::*;

    use crate::zz::constants::zz_units_sum;
    use crate::zzbase::GInt;

    use crate::zz_test;

    macro_rules! qq_div_test {
        ($name:ident, $type:ty) => {
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
        };
    }
    macro_rules! qq_div_tests {
    ($($n:expr)*) => ($(
        paste! {
            zz_test!([<qq_ring_ops $n>], [<QQ $n>]);
            qq_div_test!([<qq_div $n>], [<QQ $n>]);
         }
    )*)
}
    qq_div_tests!(4 6 8 10 12 /* 16 20 */ 24 /* 30 32 60 */);
}
