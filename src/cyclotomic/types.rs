//! In this module the different ring types are generated using macros.

use std::cmp::Ordering;
use std::fmt;
use std::fmt::Display;
use std::ops::{Add, Div, Mul, Neg, Sub};

use num_complex::Complex64;
use num_traits::{One, Pow, ToPrimitive, Zero};

use paste::paste;

use super::div::zz_inv;
use super::gaussint::GaussInt;
use super::numtraits::{Ccw, Conj, InnerIntType, IntField, IntRing, OneImag, ZSigned};
use super::params::*;
use super::sign::{
    zz_partial_signum_1_sym, zz_partial_signum_2_sym, zz_partial_signum_4_sym,
    zz_partial_signum_fallback,
};
use super::symnum::{GIntT, RatioT, SymNum, ZZComplex, ZZParams};
use super::traits::{
    ComplexTraits, FieldTraits, HasZZ10Impl, HasZZ12Impl, HasZZ4Impl, HasZZ6Impl, HasZZ8Impl,
    IsComplex, IsField, IsReal, IsRealOrComplex, IsRing, IsRingOrField, QQType, QType, RealTraits,
    RingTraits, ZType, ZZType,
};

// generate boilerplate implementations

macro_rules! impl_symnum_struct_for {
    ($($n:expr)*) => ($(
        paste! {
            impl_symnum!([<ZZ $n>], [<Z $n>], [<QQ $n>], [<Q $n>], [<ZZ $n _PARAMS>], [<zz $n _mul>]);
        }
    )*)
}
macro_rules! impl_conversions_for {
    ($($n:expr)*) => ($(
        paste! {
            impl_conversions!([<ZZ $n>], [<Z $n>], [<QQ $n>], [<Q $n>]);
        }
    )*)
}

macro_rules! zz_from_int_pair {
    ($t:ty) => {
        impl From<(<$t as InnerIntType>::IntType, <$t as InnerIntType>::IntType)> for $t {
            fn from(
                (re, im): (
                    <Self as InnerIntType>::IntType,
                    <Self as InnerIntType>::IntType,
                ),
            ) -> Self {
                Self::one().scale(re) + Self::one_i().scale(im)
            }
        }
    };
}
macro_rules! impl_from_int_pair {
    ($($n:expr)*) => ($(
        paste! {
            zz_from_int_pair!([<ZZ $n>]);
            zz_from_int_pair!([<QQ $n>]);
        }
    )*)
}
macro_rules! impl_from_real_pair {
    ($($t:ty)*) => ($(
        /// Convert from result of projection (which is only closed in the field)
        impl From<(<Self as IsField>::Real, <Self as IsField>::Real)> for $t {
            fn from((re, im): (<$t as QQType>::Real, <$t as QQType>::Real)) -> Self {
                Self::from(re) + Self::one_i() * Self::from(im)
            }
        }

        /// Convert from result of projection (which is only closed in the field)
        impl From<(<<Self as IsComplex>::Ring as IsRing>::Real, <<Self as IsComplex>::Ring as IsRing>::Real)> for $t {
            fn from((re, im): (<<Self as IsComplex>::Ring as IsRing>::Real, <<Self as IsComplex>::Ring as IsRing>::Real)) -> Self {
                Self::from(re) + Self::one_i() * Self::from(im)
            }
        }
    )*)
}

impl_symnum_struct_for!(4 6 8 10 12 16 20 24 30 32 60);
impl_functional_traits!(ZZ4, Z4, QQ4, Q4, zz_partial_signum_1_sym);
impl_functional_traits!(ZZ6, Z6, QQ6, Q6, zz_partial_signum_2_sym);
impl_functional_traits!(ZZ8, Z8, QQ8, Q8, zz_partial_signum_2_sym);
impl_functional_traits!(ZZ10, Z10, QQ10, Q10, zz_partial_signum_fallback);
impl_functional_traits!(ZZ12, Z12, QQ12, Q12, zz_partial_signum_2_sym);
impl_functional_traits!(ZZ16, Z16, QQ16, Q16, zz_partial_signum_fallback);
impl_functional_traits!(ZZ20, Z20, QQ20, Q20, zz_partial_signum_fallback);
impl_functional_traits!(ZZ24, Z24, QQ24, Q24, zz_partial_signum_4_sym);
impl_functional_traits!(ZZ30, Z30, QQ30, Q30, zz_partial_signum_fallback);
impl_functional_traits!(ZZ32, Z32, QQ32, Q32, zz_partial_signum_fallback);
impl_functional_traits!(ZZ60, Z60, QQ60, Q60, zz_partial_signum_fallback);
impl_conversions_for!(4 6 8 10 12 16 20 24 30 32 60);

zz_triv_impl!(HasZZ4Impl, ZZ4 ZZ8 ZZ12 ZZ16 ZZ20 ZZ24 ZZ32 ZZ60);
impl_from_int_pair!(4 8 12 16 20 24 32 60);
impl_from_real_pair!(QQ4 QQ8 QQ12 QQ16 QQ20 QQ24 QQ32 QQ60);

zz_triv_impl!(HasZZ6Impl, ZZ6 ZZ12 ZZ24 ZZ30 ZZ60);
zz_triv_impl!(HasZZ8Impl, ZZ8 ZZ16 ZZ24 ZZ32);
zz_triv_impl!(HasZZ10Impl, ZZ10 ZZ20 ZZ30 ZZ60);
zz_triv_impl!(HasZZ12Impl, ZZ12 ZZ24 ZZ60);

#[cfg(test)]
mod tests {
    use super::*;

    use super::super::constants::zz_units_sum;
    use super::super::symnum::{GIntT, SymNum};

    macro_rules! zz_tests {
        ($mod: ident, $type: ident, $($n:expr)*) => ($(
            paste! { zz_test!([<$mod $n>], [<$type $n>]); }
        )*)
    }

    #[test]
    fn test_display() {
        let x = ZZ24::zero();
        assert_eq!(format!("{x}"), "0");

        let x = ZZ24::one();
        assert_eq!(format!("{x}"), "1");

        let x = ZZ24::one() + ZZ24::one();
        assert_eq!(format!("{x}"), "2");

        let x = -ZZ24::one();
        assert_eq!(format!("{x}"), "-1");

        let x = ZZ24::one() + (ZZ24::ccw()).pow(2i8);
        assert_eq!(format!("{x}"), "1+1/2i + (1/2)*sqrt(3)");

        let x: ZZ10 = zz_units_sum();
        assert_eq!(
            format!("{x}"),
            "-5 + (-15/4i)*sqrt(2(5-sqrt(5))) + (-5/4i)*sqrt(10(5-sqrt(5)))"
        );
    }

    macro_rules! qq_div_test {
        ($name:ident, $type:ty) => {
            mod $name {
                use super::*;

                type QQi = $type;
                type ZZi = <$type as IsRealOrComplex>::Ring;

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
    zz_tests!(zz, ZZ, 4 6 8 10 12 16 20 24 30 32 60);
    qq_div_tests!(4 6 8 10 12 /* 16 20 */ 24 /* 30 32 60 */);
}
