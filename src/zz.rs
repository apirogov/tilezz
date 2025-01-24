//! In this module the different ring types are generated using macros.

use std::fmt;
use std::fmt::Display;
use std::ops::{Add, Mul, Neg, Sub};

use num_complex::Complex64;
use num_integer::Integer;
use num_traits::{FromPrimitive, One, ToPrimitive, Zero};

use paste::paste;

use crate::gaussint::GaussInt;
use crate::traits::{Ccw, Conj, InnerIntType, IntRing, OneImag};
use crate::zzbase::{Frac, GInt, ZCommon, ZNum, ZZBase, ZZComplex, ZZNum, ZZParams};

use crate::zzparams::*;

use crate::zzsigned::{
    zz_partial_signum_1_sym, zz_partial_signum_2_sym, zz_partial_signum_4_sym,
    zz_partial_signum_fallback, ZSigned,
};
use crate::{zz_base_impl, zz_ops_impl, zz_triv_impl};

// generate boilerplate implementations
zz_base_impl!(ZZ4, Z4, ZZ4_PARAMS, zz4_mul, zz_partial_signum_1_sym);
zz_base_impl!(ZZ6, Z6, ZZ6_PARAMS, zz6_mul, zz_partial_signum_2_sym);
zz_base_impl!(ZZ8, Z8, ZZ8_PARAMS, zz8_mul, zz_partial_signum_2_sym);
zz_base_impl!(ZZ10, Z10, ZZ10_PARAMS, zz10_mul, zz_partial_signum_fallback);
zz_base_impl!(ZZ12, Z12, ZZ12_PARAMS, zz12_mul, zz_partial_signum_2_sym);
zz_base_impl!(ZZ16, Z16, ZZ16_PARAMS, zz16_mul, zz_partial_signum_fallback);
zz_base_impl!(ZZ20, Z20, ZZ20_PARAMS, zz20_mul, zz_partial_signum_fallback);
zz_base_impl!(ZZ24, Z24, ZZ24_PARAMS, zz24_mul, zz_partial_signum_4_sym);
zz_base_impl!(ZZ30, Z30, ZZ30_PARAMS, zz30_mul, zz_partial_signum_fallback);
zz_base_impl!(ZZ32, Z32, ZZ32_PARAMS, zz32_mul, zz_partial_signum_fallback);
zz_base_impl!(ZZ60, Z60, ZZ60_PARAMS, zz60_mul, zz_partial_signum_fallback);

macro_rules! z_zz_ops_impl {
    ($($n:expr)*) => ($(
        paste! {
                zz_ops_impl!([<Z $n>]);
                zz_ops_impl!([<ZZ $n>]);
        }
    )*)
}
z_zz_ops_impl!(4 6 8 10 12 16 20 24 30 32 60);

macro_rules! zz_from_int_pair {
    ($t:ty) => {
        impl<I: Integer + FromPrimitive + ToPrimitive> From<(I, I)> for $t {
            /// Convert a Gaussian integer into a cyclotomic integer (only for fields including them)
            fn from((re, im): (I, I)) -> Self {
                let r = I::from_i64(re.to_i64().unwrap()).unwrap();
                let i = I::from_i64(im.to_i64().unwrap()).unwrap();
                Self::one().scale(r) + Self::one_i().scale(i)
            }
        }
    };
}

/// rings containing ZZ4
pub trait HasZZ4 {}
/// rings containing ZZ6
pub trait HasZZ6 {}
/// rings containing ZZ8
pub trait HasZZ8 {}
/// rings containing ZZ10
pub trait HasZZ10 {}
/// rings containing ZZ12
pub trait HasZZ12 {}

zz_triv_impl!(HasZZ4, ZZ4 ZZ8 ZZ12 ZZ16 ZZ20 ZZ24 ZZ32 ZZ60);
zz_triv_impl!(HasZZ6, ZZ6 ZZ12 ZZ24 ZZ30 ZZ60);
zz_triv_impl!(HasZZ8, ZZ8 ZZ16 ZZ24 ZZ32);
zz_triv_impl!(HasZZ10, ZZ10 ZZ20 ZZ30 ZZ60);
zz_triv_impl!(HasZZ12, ZZ12 ZZ24 ZZ60);

macro_rules! impl_from_int_pair {
    ($($t:ty)*) => ($(
            zz_from_int_pair!($t);
    )*)
}

// NOTE: list should be equal to HasZZ4 members
impl_from_int_pair!(ZZ4 ZZ8 ZZ12 ZZ16 ZZ20 ZZ24 ZZ32 ZZ60);

impl<T: ZZNum + HasZZ4> OneImag for T {
    fn one_i() -> Self {
        Self::unit(Self::qturn())
    }

    fn is_one_i(&self) -> bool {
        *self == Self::one_i()
    }
}

#[cfg(test)]
mod tests {
    use constants::zz_units_sum;

    use super::*;

    macro_rules! zz_tests {
        ($mod: ident, $type: ident, $($n:expr)*) => ($(
            paste! { zz_test!([<$mod $n>], [<$type $n>]); }
        )*)
    }
    zz_tests!(zz, ZZ, 4 6 8 10 12 16 20 24 30 32 60);

    #[test]
    fn test_signum() {
        use super::constants::*;
        let sq2: ZZ24 = sqrt2();
        let sq3: ZZ24 = sqrt3();
        let sq6: ZZ24 = sqrt6();

        let z = Z24::zero();
        let p = Z24::one();
        let m = -p;

        // use same test as above
        let sign_zz24 = |a, b, c, d| {
            (ZZ24::from(a) + ZZ24::from(b) * sq2 + ZZ24::from(c) * sq3 + ZZ24::from(d) * sq6)
                .re()
                .signum()
        };

        let (a, b, c, d) = (485, 343, 280, 198);
        assert_eq!(sign_zz24(0, 0, 0, 0), z);
        assert_eq!(sign_zz24(-a, -b, c, d), m);
        assert_eq!(sign_zz24(-a, b, -c, d), p);
        assert_eq!(sign_zz24(-a, b, c, -d), p);
        assert_eq!(sign_zz24(a, -b, -c, d), m);
        assert_eq!(sign_zz24(a, -b, c, -d), m);
        assert_eq!(sign_zz24(a, b, -c, -d), p);
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

        let x = ZZ24::one() + (ZZ24::ccw()).pow(2);
        assert_eq!(format!("{x}"), "1+1/2i + (1/2)*sqrt(3)");

        let x: ZZ10 = zz_units_sum();
        assert_eq!(
            format!("{x}"),
            "-5 + (-15/4i)*sqrt(2(5-sqrt(5))) + (-5/4i)*sqrt(10(5-sqrt(5)))"
        );
    }
}

pub mod constants {
    use super::*;

    // Returns the sum of all units of a complex integer ring.
    pub fn zz_units_sum<T: ZZNum>() -> T {
        let mut p = T::zero();
        for i in 0..T::turn() {
            p = p + T::unit(i).scale(i as i64);
        }
        p
    }

    // NOTE: as we can get the real-valued square roots represented,
    // it means that we can represent any linear combination
    // in a ring that supports quarter turn rotation (i.e. ZZDiv4).

    pub fn sqrt2<T: ZZNum + HasZZ8>() -> T {
        let sc = T::zz_params().full_turn_steps / 8;
        T::unit(sc) + T::unit(-sc)
    }

    pub fn sqrt3<T: ZZNum + HasZZ12>() -> T {
        let sc = T::zz_params().full_turn_steps / 12;
        T::unit(sc) + T::unit(-sc)
    }

    pub fn sqrt5<T: ZZNum + HasZZ10>() -> T {
        let sc = T::zz_params().full_turn_steps / 10;
        (T::unit(sc) + T::unit(-sc)) * T::one().scale(2) - T::one()
    }

    pub fn sqrt6<T: ZZNum + HasZZ8 + HasZZ12>() -> T {
        let sc = T::zz_params().full_turn_steps / 24;
        (T::unit(sc) + T::unit(-sc)) * T::one().scale(2) - sqrt2::<T>()
    }

    // misc. irregular
    pub fn isqrt3() -> ZZ6 {
        ZZ6::unit(1) + ZZ6::unit(2)
    }
    pub fn zz10_isqrt_penta() -> ZZ10 {
        ZZ10::unit(1) * ZZ10::from(4) - ZZ10::one() - sqrt5()
    }
    pub fn zz20_half_sqrt_penta() -> ZZ20 {
        ZZ20::unit(3) + ZZ20::unit(-3)
    }

    #[cfg(test)]
    mod tests {
        #[test]
        fn test_constants() {
            use super::*;
            use std::f64::consts::SQRT_2;

            let sq2 = SQRT_2;
            let sq3 = 3.0_f64.sqrt();
            let sq_penta = ZZ10_Y.sqrt();
            let hsq_penta = 0.5 * ZZ10_Y.sqrt();
            let sq5 = 5.0_f64.sqrt();
            let sq6 = 6.0_f64.sqrt();

            assert_eq!(sqrt2::<ZZ8>().complex64().re, sq2);
            assert_eq!(sqrt2::<ZZ16>().complex64().re, sq2);
            assert_eq!(sqrt2::<ZZ24>().complex64().re, sq2);
            assert_eq!(sqrt2::<ZZ32>().complex64().re, sq2);

            assert_eq!(sqrt3::<ZZ12>().complex64().re, sq3);
            assert_eq!(sqrt3::<ZZ24>().complex64().re, sq3);
            assert_eq!(sqrt3::<ZZ60>().complex64().re, sq3);

            assert_eq!(sqrt5::<ZZ10>().complex64().re, sq5);
            assert_eq!(sqrt5::<ZZ20>().complex64().re, sq5);
            assert_eq!(sqrt5::<ZZ30>().complex64().re, sq5);
            assert_eq!(sqrt5::<ZZ60>().complex64().re, sq5);

            assert_eq!(sqrt6::<ZZ24>().complex64().re, sq6);
            // assert_eq!(sqrt6::<ZZ120>().complex().re, sq6);
            // assert_eq!(sqrt6::<ZZ240>().complex().re, sq6);

            assert_eq!(isqrt3().complex64().im, sq3);
            assert_eq!(zz10_isqrt_penta().complex64().im, sq_penta);
            assert_eq!(zz20_half_sqrt_penta().complex64().re, hsq_penta);
        }
    }
}
