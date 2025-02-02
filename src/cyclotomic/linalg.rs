//! Core linear algebra utils
use super::traits::{IsComplex, IsRingOrField};
use num_traits::Zero;

use super::numtraits::ZSigned;

/// Dot product between two values, i.e. a dot b = |a||b|cos(ab).
pub fn dot<ZZ: IsRingOrField + IsComplex>(p1: &ZZ, p2: &ZZ) -> ZZ::Real {
    let (p1x, p1y) = p1.re_im();
    let (p2x, p2y) = p2.re_im();
    (p1x * p2x) + (p1y * p2y)
}

/// 2D "cross" product (wedge is the general name),
/// i.e. a wedge b = |a|b|sin(ab).
pub fn wedge<ZZ: IsRingOrField + IsComplex>(p1: &ZZ, p2: &ZZ) -> ZZ::Real {
    let (p1x, p1y) = p1.re_im();
    let (p2x, p2y) = p2.re_im();
    (p1x * p2y) - (p1y * p2x)
}

/// Squared norm square of a value, i.e. dot product of the value with itself.
pub fn norm_sq<ZZ: IsRingOrField + IsComplex>(p: &ZZ) -> ZZ::Real {
    dot::<ZZ>(&p, &p)
}

/// Projection of v onto u, i.e. dot(u,v)/|u|^2 * u
pub fn project<ZZ: IsRingOrField + IsComplex + From<ZZ::Real>>(v: &ZZ, u: &ZZ) -> ZZ::Field
where
    ZZ::Field: From<ZZ>,
{
    let v_dot_u = ZZ::Field::from(dot(v, u).into());
    let norm_u_sq = ZZ::Field::from(norm_sq(u).into());
    let u_val = ZZ::Field::from(*u);
    v_dot_u / norm_u_sq * u_val
}

/// Return true if angle (wrt. positive real line) is in closed interval `[a,b]`,
/// assuming that a and b are in counterclockwise order and their ccw angle
/// is less than a half turn.
pub fn angle_between<ZZ: IsRingOrField + IsComplex>(p: &ZZ, (a, b): (&ZZ, &ZZ)) -> bool {
    !wedge(a, p).is_negative() && !wedge(p, b).is_negative()
}

/// Return whether this point is strictly between the other two.
///
/// NOTE: we already assume all three involved points are colinear.
pub fn is_between<ZZ: IsRingOrField + IsComplex>(p: &ZZ, (a, b): (&ZZ, &ZZ)) -> bool {
    let v = *a - *p;
    let w = *p - *b;
    wedge(&v, &w).is_zero() && dot(&v, &w).is_positive()
}

/// Return whether the segments pa and pb (in that order) have a ccw angle,
/// i.e. b is ccw of a with respect to rotation around p.
pub fn is_ccw<ZZ: IsRingOrField + IsComplex>(p: &ZZ, (a, b): (&ZZ, &ZZ)) -> bool {
    wedge(&(*a - *p), &(*b - *p)).is_positive()
    // NOTE: wedge(*a - *p, *p - *b).is_zero() is subexp. of is_between... interesting symmetry
}

#[cfg(test)]
mod tests {
    use num_traits::{One, Pow, Zero};

    use super::super::constants::zz_units_sum;
    use super::super::symnum::SymNum;
    use super::super::types::{Z12, ZZ12};
    use super::*;

    type ZZi = ZZ12;

    #[test]
    fn test_dot() {
        // get a non-trivial point
        let p: ZZi = zz_units_sum();

        // all rotations of the same point around origin
        // have the same squared distance, i.e. quadrance
        // and it is real-valued.
        let q: ZZi = norm_sq(&p).into();
        for i in 1..ZZi::turn() {
            let pi = p * ZZi::unit(i);
            let qi: ZZi = norm_sq(&pi).into();
            assert_eq!(qi, q);
        }
    }

    #[test]
    fn test_dot_extra() {
        let p1 = ZZ12::one();
        let p2 = ZZ12::from(2);
        let p3 = ZZ12::from(3);
        let pi = ZZ12::unit(3); // i
        let p60 = ZZ12::unit(2);
        let pm60 = ZZ12::unit(-2);

        assert_eq!(dot(&p1, &pi), Z12::zero());
        assert_eq!(dot(&0.into(), &pi), Z12::zero());
        assert_eq!(dot(&p2, &p3), Z12::from(6));

        // {0, 1} dot {1/2, sqrt(3)/2} = sqrt(3) / 2
        // => dot^2 = 3/4
        let d1 = ZZ12::from(dot(&p60, &pi).pow(2_u8)).complex64();
        assert_eq!(d1.re, 0.75);
        assert_eq!(d1.im, 0.0);

        // same but with negative sign (check indirectly)
        let d2 = ZZ12::from(dot(&pm60, &pi)).complex64();
        assert!(d2.re < 0.0);
        assert_eq!(d2.im, 0.0);
        assert_eq!(ZZ12::from(dot(&pm60, &pi).pow(2_u8)).complex64().re, 0.75);
    }

    #[test]
    fn test_is_between() {
        let a: ZZi = ZZi::zero();
        let b: ZZi = ZZi::one();
        let c: ZZi = ZZi::from(2);
        let e: ZZi = ZZi::unit(ZZi::hturn() / 2);
        let f: ZZi = b + e;
        let g: ZZi = ZZi::unit(1) + ZZi::unit(-1) - ZZi::one();

        // is actually in between
        assert!(is_between(&b, (&a, &c)));
        assert!(is_between(&g, (&a, &b)));
        // is an endpoint
        assert!(!is_between(&b, (&a, &b)));
        assert!(!is_between(&a, (&a, &b)));
        // colinear, but not between
        assert!(!is_between(&c, (&a, &b)));
        // not colinear
        assert!(!is_between(&f, (&a, &b)));
    }

    #[test]
    fn test_is_ccw() {
        type ZZ = ZZ12;
        assert!(is_ccw(&ZZ::zero(), (&ZZ::unit(0), &ZZ::unit(1))));
        assert!(is_ccw(&ZZ::zero(), (&ZZ::unit(2), &ZZ::unit(3))));
        assert!(!is_ccw(&ZZ::zero(), (&ZZ::unit(0), &ZZ::unit(0))));
        assert!(!is_ccw(&ZZ::zero(), (&ZZ::unit(1), &ZZ::unit(1))));
        assert!(!is_ccw(&ZZ::zero(), (&ZZ::unit(1), &ZZ::unit(0))));
        assert!(!is_ccw(&ZZ::zero(), (&ZZ::unit(6), &ZZ::unit(5))));
    }

    #[test]
    fn test_angle_between() {
        type ZZ = ZZ12;
        // test inside
        assert!(angle_between(&ZZ::unit(1), (&ZZ::unit(0), &ZZ::unit(2))));
        assert!(angle_between(&ZZ::unit(-2), (&ZZ::unit(-3), &ZZ::unit(1))));
        // test outside
        assert!(!angle_between(&ZZ::unit(1), (&ZZ::unit(-1), &ZZ::unit(0))));
        // bound inclusive
        assert!(angle_between(&ZZ::unit(1), (&ZZ::unit(1), &ZZ::unit(2))));
        assert!(angle_between(&ZZ::unit(2), (&ZZ::unit(1), &ZZ::unit(2))));
    }

    #[test]
    fn test_project() {
        use super::super::numtraits::OneImag;
        use super::super::traits::ZZType;
        type ZZ = ZZ12;
        type QQ = <ZZ as ZZType>::Field;
        let r: ZZ = zz_units_sum();
        let re = project(&r.into(), &QQ::one());
        let im = project(&r.into(), &QQ::one_i());
        println!("{r}\n{re}\n{im}\n{}", re + im);
        assert_eq!(re + im, r.into());
    }
}
