use super::gaussint::to_gint;
use super::traits::{Ccw, InnerIntType, IntRing};
use super::zzbase::{Frac, GInt, ZZBase, ZZNum, ZZParams};
use crate::traits::ComplexIntRing;
use crate::{zz_base_impl, zz_ops_impl};

use num_traits::{One, Zero};
use std::fmt;
use std::fmt::Display;
use std::marker::PhantomData;
use std::ops::{Add, Mul, Neg, Sub};

// definitions needed to derive different ZZn types

// numeric constants
const SQRT_5: f64 = 2.23606797749978969;
const PENTA: f64 = 2.0 * (5.0 - SQRT_5);
// const PENTA_POS: f64 = 2.0 * (5.0 + SQRT_5);

/// Gauss integers
pub const ZZ4_PARAMS: ZZParams<Frac> = ZZParams {
    phantom: PhantomData,
    full_turn_steps: 4,
    sym_roots_num: 1,
    sym_roots_sqs: &[1.0],
    scaling_fac: 1,
    ccw_unit_coeffs: &[[0, 1]],
};
fn zz4_mul(x: &[GInt], y: &[GInt]) -> Vec<GInt> {
    match [*array_ref!(x, 0, 1), *array_ref!(y, 0, 1)] {
        [[a], [b]] => vec![a * b],
    }
}
/// Eisenstein integers
pub const ZZ6_PARAMS: ZZParams<Frac> = ZZParams {
    phantom: PhantomData,
    full_turn_steps: 6,
    sym_roots_num: 2,
    sym_roots_sqs: &[1.0, 3.0],
    scaling_fac: 2,
    ccw_unit_coeffs: &[[1, 0], [0, 1]],
};
/// Dimension multiplication matrix (for Z[i]-valued vectors):
///    c d
/// a [1 s]
/// b [s 3]
/// where s = sqrt(3)
fn zz6_mul(x: &[GInt], y: &[GInt]) -> Vec<GInt> {
    match [*array_ref!(x, 0, 2), *array_ref!(y, 0, 2)] {
        [[a, b], [c, d]] => vec![a * c + ((3,) * b * d), a * d + b * c],
    }
}
/// Compass integers
pub const ZZ8_PARAMS: ZZParams<Frac> = ZZParams {
    phantom: PhantomData,
    full_turn_steps: 8,
    sym_roots_num: 2,
    sym_roots_sqs: &[1.0, 2.0],
    scaling_fac: 2,
    ccw_unit_coeffs: &[[0, 0], [1, 1]],
};
/// Dimension multiplication matrix (for Z[i]-valued vectors):
///    c d
/// a [1 s]
/// b [s 2]
/// where s = sqrt(2)
fn zz8_mul(x: &[GInt], y: &[GInt]) -> Vec<GInt> {
    match [*array_ref!(x, 0, 2), *array_ref!(y, 0, 2)] {
        [[a, b], [c, d]] => vec![a * c + ((2,) * b * d), a * d + b * c],
    }
}
/// Halfrose integers (impractical, has no quarter turn)
pub const ZZ10_PARAMS: ZZParams<Frac> = ZZParams {
    phantom: PhantomData,
    full_turn_steps: 10,
    sym_roots_num: 4,
    sym_roots_sqs: &[1.0, 5.0, PENTA, 5.0 * PENTA],
    scaling_fac: 8,
    ccw_unit_coeffs: &[[2, 0], [2, 0], [0, 2], [0, 0]],
};
/// Dimension multiplication matrix (for Z[i]-valued vectors):
///    e     f      g        h
/// a [1   ,  x   ,   y    ,   xy   ]
/// b [ x  , 5    ,  xy    ,  5 y   ]
/// c [  y ,  xy  , 10-2x  , 10(x-1)]
/// d [ xy , 5 y  , 10(x-1), 10(5-x)]
/// where x = sqrt(5), y = sqrt(2*(5-sqrt(5)))
fn zz10_mul(x: &[GInt], y: &[GInt]) -> Vec<GInt> {
    match [*array_ref!(x, 0, 4), *array_ref!(y, 0, 4)] {
        [[a, b, c, d], [e, f, g, h]] => {
            let w = a * e + (5,) * b * f + (10,) * (c * g + (5,) * d * h - c * h - d * g);
            let x = a * f + b * e - (2,) * c * g + (10,) * (c * h + d * g - d * h);
            let y = a * g + (5,) * (b * h + d * f) + c * e;
            let z = a * h + b * g + c * f + d * e;
            vec![w, x, y, z]
        }
    }
}
/// Clock integers
pub const ZZ12_PARAMS: ZZParams<Frac> = ZZParams {
    phantom: PhantomData,
    full_turn_steps: 12,
    sym_roots_num: 2,
    sym_roots_sqs: &[1.0, 3.0],
    scaling_fac: 2,
    ccw_unit_coeffs: &[[0, 1], [1, 0]],
};
fn zz12_mul(x: &[GInt], y: &[GInt]) -> Vec<GInt> {
    return zz6_mul(x, y);
}
/// Penrose integers
pub const ZZ20_PARAMS: ZZParams<Frac> = ZZParams {
    phantom: PhantomData,
    full_turn_steps: 20,
    sym_roots_num: 4,
    sym_roots_sqs: &[1.0, 5.0, PENTA, 5.0 * PENTA],
    scaling_fac: 8,
    ccw_unit_coeffs: &[[0, -2], [0, 2], [1, 0], [1, 0]],
};
/// Let x = sqrt(5), y = sqrt(2*(5-sqrt(5))), z = sqrt(2*(5+sqrt(5)))
/// We have that e^(i*pi/10) = 1/4(-i + ix + z)
/// But as 1/2(xz - z) = y ^ 1/2(xy + y) = z,
/// we can reuse ZZ10 logic (which has y = PENTA), as 2z = xy + y
fn zz20_mul(x: &[GInt], y: &[GInt]) -> Vec<GInt> {
    zz10_mul(x, y)
}
/// Digiclock integers
pub const ZZ24_PARAMS: ZZParams<Frac> = ZZParams {
    phantom: PhantomData,
    full_turn_steps: 24,
    sym_roots_num: 4,
    sym_roots_sqs: &[1.0, 2.0, 3.0, 6.0],
    scaling_fac: 4,
    ccw_unit_coeffs: &[[0, 0], [1, -1], [0, 0], [1, 1]],
};
/// Dimension multiplication matrix (for Z[i]-valued vectors):
///     e    f     g    h
/// a [1   ,  x  ,   y ,  xy ]
/// b [ x  , 2   ,  xy , 2 y ]
/// c [  y ,  xy , 3   , 3x  ]
/// d [ xy , 2 y , 3x  , 6   ]
fn zz24_mul(x: &[GInt], y: &[GInt]) -> Vec<GInt> {
    match [*array_ref!(x, 0, 4), *array_ref!(y, 0, 4)] {
        [[a, b, c, d], [e, f, g, h]] => {
            let w = a * e + (2,) * b * f + (3,) * c * g + (6,) * d * h;
            let x = a * f + b * e + (3,) * (c * h + d * g);
            let y = a * g + c * e + (2,) * (b * h + d * f);
            let z = a * h + b * g + c * f + d * e;
            vec![w, x, y, z]
        }
    }
}
// --------

// generate boilerplate implementations
// TODO: is there a practical representation of ZZ16?
zz_base_impl!(ZZ4, ZZ4_PARAMS, zz4_mul);
zz_base_impl!(ZZ6, ZZ6_PARAMS, zz6_mul);
zz_base_impl!(ZZ8, ZZ8_PARAMS, zz8_mul);
zz_base_impl!(ZZ10, ZZ10_PARAMS, zz10_mul);
zz_base_impl!(ZZ12, ZZ12_PARAMS, zz12_mul);
zz_base_impl!(ZZ20, ZZ20_PARAMS, zz20_mul);
zz_base_impl!(ZZ24, ZZ24_PARAMS, zz24_mul);
zz_ops_impl!(ZZ4 ZZ6 ZZ8 ZZ10 ZZ12 ZZ20 ZZ24);

// implementations for different complex integer rings
// (using default underlying integer type on the innermost level)
impl ZZNum for ZZ4 {}
impl ZZNum for ZZ6 {}
impl ZZNum for ZZ8 {}
impl ZZNum for ZZ10 {}
impl ZZNum for ZZ12 {}
impl ZZNum for ZZ20 {}
impl ZZNum for ZZ24 {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gaussint::to_gint2;
    use num_complex::Complex64;
    use std::collections::HashSet;

    // TODO: make macro to generate the tests for all instances
    // https://eli.thegreenplace.net/2021/testing-multiple-implementations-of-a-trait-in-rust/
    type ZZi = ZZ20;

    #[test]
    fn test_basic() {
        // # of full turn steps is an even natural (so a half turn is possible)
        assert!(ZZi::turn() > 0);
        assert!(ZZi::turn() % 2 == 0);

        // check vector sizes
        let roots_num = ZZi::zz_params().sym_roots_num;
        assert_eq!(ZZi::zz_params().sym_roots_sqs.len(), roots_num);
        assert_eq!(ZZi::zz_params().ccw_unit_coeffs.len(), roots_num);

        // check zero-vector
        let z = ZZi::zero();
        let cs = z.zz_coeffs();
        assert_eq!(cs.len(), roots_num);
        for i in 0..roots_num {
            assert_eq!(cs[i], GInt::zero());
        }

        // check one-vector
        let o = ZZi::one();
        let cs = o.zz_coeffs();
        assert_eq!(cs.len(), roots_num);
        assert_eq!(cs[0], GInt::one());
        for i in 1..roots_num {
            assert_eq!(cs[i], GInt::zero());
        }
    }

    #[test]
    fn test_add_sub() {
        // test addition / subtraction
        assert_eq!(ZZi::zero() + ZZi::zero(), ZZi::zero());
        assert_eq!(ZZi::one() + ZZi::zero(), ZZi::one());
        assert_eq!(ZZi::zero() + ZZi::one(), ZZi::one());
        assert_eq!(-ZZi::one() + ZZi::one(), ZZi::zero());
        assert_eq!(ZZi::one() - ZZi::one(), ZZi::zero());
    }

    #[test]
    fn test_mul() {
        // test scalar multiplication
        assert_eq!(ZZi::zero().scale(2), ZZi::zero());
        assert_eq!(ZZi::ccw().scale(3), ZZi::ccw() + ZZi::ccw() + ZZi::ccw());
        assert_eq!(ZZi::one().scale(-1), -ZZi::one());
        assert_eq!(ZZi::one().scale(-42), ZZi::from_int(-42));

        // test multiplication
        assert_eq!(ZZi::zero() * ZZi::zero(), ZZi::zero());
        assert_eq!(ZZi::one() * ZZi::zero(), ZZi::zero());
        assert_eq!(ZZi::zero() * ZZi::one(), ZZi::zero());
        assert_eq!(ZZi::one() * ZZi::one(), ZZi::one());
        assert_eq!(-ZZi::one() * ZZi::one(), -ZZi::one());
        assert_eq!(ZZi::one() * (-ZZi::one()), -ZZi::one());
        assert_eq!((-ZZi::one()) * (-ZZi::one()), ZZi::one());
    }

    #[test]
    fn test_rotations() {
        // test ccw()
        assert_eq!(ZZi::ccw() * ZZi::ccw().conj(), ZZi::one());
        assert_eq!(-(-(ZZi::one()) * ZZi::ccw()), ZZi::ccw());

        // test going around the unit circle step by step
        let mut x = ZZi::one();
        for _ in 0..ZZi::turn() {
            x = x * ZZi::ccw();
        }
        assert_eq!(x, ZZi::one());

        // test unit()
        assert_eq!(ZZi::unit(0), ZZi::one());
        assert_eq!(ZZi::unit(-1), ZZi::unit(ZZi::turn() - 1));
        assert_eq!(ZZi::unit(1), ZZi::unit(ZZi::turn() + 1));
        assert_eq!(ZZi::unit(-ZZi::hturn()), ZZi::unit(ZZi::hturn()));
        assert_eq!(ZZi::unit(ZZi::hturn()), -ZZi::one());
        if ZZi::turn() % 4 == 0 {
            assert_eq!(ZZi::one_i().zz_coeffs()[0], to_gint2(0, 1));
        }

        // test powi()
        assert_eq!(ZZi::ccw().powi(ZZi::hturn()), -ZZi::one());
        assert_eq!(ZZi::ccw().powi(ZZi::turn()), ZZi::one());
        assert_eq!(ZZi::ccw().powi(ZZi::hturn()).powi(2), ZZi::one());
    }

    #[test]
    #[should_panic]
    fn test_neg_powi() {
        ZZi::one().powi(-1);
    }

    #[test]
    fn test_scaling_fac() {
        // test scaling fac is correct by checking denom. of coeffs of all units
        // (that the denom. always can be expressed as multple of scaling factor)
        // and that the chosen constant factor is indeed minimal
        let sc_fac = ZZi::zz_params().scaling_fac;
        let mut max_fac: i64 = 0;
        for i in 0..ZZi::turn() {
            let x = ZZi::unit(i);
            println!("{x}");
            for c in x.coeffs {
                assert_eq!(sc_fac % c.real.denom(), 0);
                assert_eq!(sc_fac % c.imag.denom(), 0);
                max_fac = max_fac.max(*c.real.denom());
                max_fac = max_fac.max(*c.imag.denom());
            }
        }
        assert_eq!(sc_fac, max_fac);
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

        let x = ZZ24::one() + (ZZ24::ccw()).powi(2);
        assert_eq!(format!("{x}"), "1+1/2i + (1/2)*sqrt(3)");
    }

    #[test]
    fn test_complex() {
        let x = ZZi::zero();
        assert_eq!(x.complex(), Complex64::zero());
        let x = ZZi::one();
        assert_eq!(x.complex(), Complex64::one());
        let x = -ZZi::one();
        assert_eq!(x.complex(), -Complex64::one());
        let x = ZZi::one() + ZZi::one();
        assert_eq!(x.complex(), Complex64::new(2.0, 0.0));

        let x = ZZi::ccw();
        let c = x.complex();
        println!("{c} = {x}");
    }

    #[test]
    fn test_hashable() {
        let mut s: HashSet<ZZi> = HashSet::new();
        s.insert(ZZi::zero());
        s.insert(ZZi::one());
        s.insert(ZZi::ccw());
        assert!(s.contains(&ZZi::ccw()));
        assert!(s.contains(&(ZZi::ccw() + ZZi::zero())));
        assert!(!s.contains(&(ZZi::ccw() + ZZi::one())));
    }

    // Test point layout:
    // -------
    // E F
    //
    // A B C D
    // -------
    // static A: ZZi = ZZi::zero();
    // static B: GaussInt<i64> = ZZi::one();
    // static C: GaussInt<i64> = ZZi::from_int(2);
    // static D: GaussInt<i64> = GaussInt::new(3, 0);
    // static E: GaussInt<i64> = GaussInt::new(0, 1);
    // static F: GaussInt<i64> = GaussInt::new(1, 1);

    fn get_non_triv_point() -> ZZ12 {
        // -0.035898384862245614-0.401923788646684i
        let tmp1 = (ZZ12::one().scale(2) - ZZ12::unit(2) - ZZ12::unit(-1).scale(2)) * ZZ12::unit(1);
        let tmp2 = (ZZ12::one().scale(3) - ZZ12::unit(2) - ZZ12::unit(-1).scale(3)) * ZZ12::unit(2);
        (tmp1 - tmp2) * ZZ12::unit(-2)
    }

    #[test]
    fn test_xy() {
        // imaginary unit (note that we only have that
        // in ZZi where i is divisible by 4, so in general
        // the splitting operation may not be reversible
        // with ring operations on the two parts.
        let i = ZZ12::unit(3);

        // test correctness of splitting and reconstruction
        let p = get_non_triv_point();
        let (x, y) = p.xy();
        assert_eq!(p, x + y * i);
    }

    #[test]
    fn test_is_real_imag_complex() {
        assert!(ZZ12::zero().is_real());
        assert!(ZZ12::zero().is_imag());

        assert!(ZZ12::one().is_real());
        assert!((-ZZ12::one()).is_real());
        assert!(!ZZ12::one().is_imag());
        assert!(!ZZ12::one().is_complex());

        assert!(!ZZ12::unit(1).is_real());
        assert!(!ZZ12::unit(2).is_imag());
        assert!(ZZ12::unit(1).is_complex());
        assert!(ZZ12::unit(2).is_complex());

        assert!(!ZZ12::unit(3).is_real());
        assert!(ZZ12::unit(3).is_imag());
        assert!((-ZZ12::unit(3)).is_imag());
        assert!(!ZZ12::unit(3).is_complex());
    }

    #[test]
    fn test_dot() {
        let p1 = ZZ12::one();
        let p2 = ZZ12::from_int(2);
        let p3 = ZZ12::from_int(3);
        let pi = ZZ12::unit(3); // i
        let p60 = ZZ12::unit(2);
        let pm60 = ZZ12::unit(-2);

        assert_eq!(p1.dot(&pi), ZZ12::zero());
        assert_eq!(ZZ12::zero().dot(&pi), ZZ12::zero());
        assert_eq!(p2.dot(&p3), ZZ12::from_int(6));

        // {0, 1} dot {1/2, sqrt(3)/2} = sqrt(3) / 2
        // => dot^2 = 3/4
        let d1 = p60.dot(&pi).powi(2).complex();
        assert_eq!(d1.re, 0.75);
        assert_eq!(d1.im, 0.0);

        // same but with negative sign (check indirectly)
        let d2 = pm60.dot(&pi).complex();
        assert!(d2.re < 0.0);
        assert_eq!(d2.im, 0.0);
        assert_eq!(pm60.dot(&pi).powi(2).complex().re, 0.75);

        let p = get_non_triv_point();
        let q = p.norm_sq();

        // all rotations of the same point around origin
        // have the same squared distance, i.e. quadrance
        for i in 1..ZZ12::turn() {
            let pi = p * ZZ12::unit(i);
            let qi = pi.norm_sq();
            let ci = pi.complex();
            let ni = ci.norm();
            println!("{ci} {ni} {qi}");
            assert_eq!(qi, q);
        }
    }

    #[test]
    fn test_is_between() {
        let a: ZZi = ZZi::zero();
        let b: ZZi = ZZi::one();
        let c: ZZi = ZZi::from_int(2);
        let e: ZZi = ZZi::unit(ZZi::hturn() / 2);
        let f: ZZi = b + e;
        let g: ZZi = ZZi::unit(1) + ZZi::unit(-1) - ZZi::one();

        // is actually in between
        assert!(b.is_between(&a, &c));
        assert!(g.is_between(&a, &b));
        // is an endpoint
        assert!(!b.is_between(&a, &b));
        assert!(!a.is_between(&a, &b));
        // colinear, but not between
        assert!(!c.is_between(&a, &b));
        // not colinear
        assert!(!f.is_between(&a, &b));
    }

    #[test]
    fn test_colinear() {
        let a: ZZi = ZZi::zero();
        let b: ZZi = ZZi::one();
        let c: ZZi = ZZi::from_int(2);
        let d: ZZi = ZZi::from_int(3);
        let e: ZZi = ZZi::unit(ZZi::hturn() / 2);
        let f: ZZi = b + e;

        let l_ab = a.line_through(&b);
        let l_ac = a.line_through(&c);
        let l_af = a.line_through(&f);

        // colinear, overlap
        assert!(b.is_colinear(&l_ac));
        assert!(d.is_colinear(&l_ac));
        // colinear, no overlap
        assert!(c.is_colinear(&l_ab));
        assert!(d.is_colinear(&l_ab));
        // parallel (not colinear)
        assert!(!e.is_colinear(&l_ab));
        assert!(!f.is_colinear(&l_ab));
        // perpendicular (touches in one point, not in the other)
        assert!(!(a.is_colinear(&l_ab) && e.is_colinear(&l_ab)));
        // general case
        assert!(!b.is_colinear(&l_af));
        assert!(!d.is_colinear(&l_af));
    }
}
