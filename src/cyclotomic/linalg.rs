//! Core linear algebra utils
use num_traits::{One, Zero};

use super::numtraits::ZSigned;
use super::symnum::SymNum;
use super::traits::{IsComplex, IsReal, IsRingOrField};

/// Dot product between two values, i.e. a dot b = |a||b|cos(ab).
pub fn dot<ZZ: IsComplex>(p1: &ZZ, p2: &ZZ) -> ZZ::Real {
    let (p1x, p1y) = p1.re_im();
    let (p2x, p2y) = p2.re_im();
    (p1x * p2x) + (p1y * p2y)
}

/// 2D "cross" product (wedge is the general name),
/// i.e. a wedge b = |a|b|sin(ab).
pub fn wedge<ZZ: IsComplex>(p1: &ZZ, p2: &ZZ) -> ZZ::Real {
    let (p1x, p1y) = p1.re_im();
    let (p2x, p2y) = p2.re_im();
    (p1x * p2y) - (p1y * p2x)
}

/// Squared norm square of a value, i.e. dot product of the value with itself.
pub fn norm_sq<ZZ: IsComplex>(p: &ZZ) -> ZZ::Real {
    dot::<ZZ>(&p, &p)
}

/// Projection of v onto u, i.e. dot(u,v)/|u|^2 * u
pub fn project<ZZ: IsComplex + From<ZZ::Real>>(v: &ZZ, u: &ZZ) -> ZZ::Field
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

pub type Angle<ZZ> = <<ZZ as IsRingOrField>::Real as IsReal>::Field;

/// Pseudo-sinus (in Wilderberger, this quantity is called "spread").
/// The value is between [-1, 1] for |p|>=1 and is 0 iff p is on the real axis.
/// For any rational n, pseudo_sin(n*p) = pseudo_sin(p) / n
pub fn pseudo_sin<ZZ: IsComplex>(p: &ZZ) -> Angle<ZZ>
where
    Angle<ZZ>: From<<ZZ as IsRingOrField>::Real>,
{
    type RealField<ZZ> = <<ZZ as IsRingOrField>::Real as IsReal>::Field;
    let numer: RealField<ZZ> = wedge(&ZZ::one(), p).into();
    let denom: RealField<ZZ> = norm_sq(p).into();
    numer / denom
}

/// Pseudo-cosinus (in Wilderberger, this quantity is called "cross").
/// The value is between [-1, 1] for |p|>=1 and is 0 iff p is on the imaginary axis.
/// For any rational n, pseudo_cos(n*p) = pseudo_cos(p) / n
pub fn pseudo_cos<ZZ: IsComplex>(p: &ZZ) -> Angle<ZZ>
where
    Angle<ZZ>: From<<ZZ as IsRingOrField>::Real>,
{
    type RealField<ZZ> = <<ZZ as IsRingOrField>::Real as IsReal>::Field;
    let numer: RealField<ZZ> = dot(&ZZ::one(), p).into();
    let denom: RealField<ZZ> = norm_sq(p).into();
    numer / denom
}

/// Returns information about the quadrant where a point is located in,
/// based on its pseudo-angles s, c.
///       3
///     4   2
///   5   0   1
///     6   8
///       7
/// where 0 indicates location of the origin.
fn quadrant<ZZ: IsComplex>(s: &Angle<ZZ>, c: &Angle<ZZ>) -> i8
where
    Angle<ZZ>: From<<ZZ as IsRingOrField>::Real>,
{
    // if p.is_zero() {
    //     return 0;
    // }
    // s = 0 <=> parallel to x axis ^ c = 0 <=> perp to x axis
    // s > 0 <=> positive y axis ^ c > 0 <=> positive x axis
    let (s_sign, c_sign) = (s.signum(), c.signum());
    if s_sign.is_zero() {
        if c_sign.is_one() {
            1
        } else {
            5
        }
    } else if c_sign.is_zero() {
        if s_sign.is_one() {
            3
        } else {
            7
        }
    } else if s_sign.is_one() {
        if c_sign.is_one() {
            2
        } else {
            4
        }
    } else {
        // s_sign = -1
        if c_sign.is_one() {
            8
        } else {
            6
        }
    }
}

/// Pseudo-atan2, returns a rational pseudo-angle of p between [0, 4),
/// where
/// * 0 means that p is on the positive real axis,
/// * 1 means that p is on the positive imaginary axis,
/// * 2 means that p is on the negative real axis,
/// * 3 means that p is on the negative imaginary axis.
///
/// **NOTE:** The pseudo-angle is **not norm-invariant**, i.e. it can only be
/// used to compare angles of points with the *same* (squared) norm!
pub fn pseudo_angle<ZZ: IsComplex>(p: &ZZ) -> Angle<ZZ>
where
    Angle<ZZ>: From<<ZZ as IsRingOrField>::Real>,
{
    let (s, c) = (pseudo_sin::<ZZ>(&p), pseudo_cos::<ZZ>(&p));
    let one = Angle::<ZZ>::one();
    match quadrant::<ZZ>(&s, &c) {
        1 => Angle::<ZZ>::zero(),
        2 => s,
        3 => one,
        4 => one + (one - s),
        5 => one.scale(2),
        6 => one.scale(2) - s,
        7 => one.scale(3),
        8 => one.scale(3) + (s + one),
        _ => panic!("Something unexpected happened!"),
    }
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
    fn test_norm_sq() {
        let one_half: <ZZ12 as IsComplex>::Field =
            <ZZ12 as IsComplex>::Field::from(1) / <ZZ12 as IsComplex>::Field::from(2);
        assert_eq!(ZZ12::from(norm_sq(&ZZ12::unit(1))), 1.into());
        assert_eq!(ZZ12::from(norm_sq(&ZZ12::unit(1).scale(2))), 4.into());
        assert_eq!(
            <ZZ12 as IsComplex>::Field::from(norm_sq(&(one_half * ZZ12::unit(1).into()))),
            one_half * one_half
        );
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

    #[test]
    fn test_sin_cos() {
        use crate::prelude::constants::sqrt3;

        type Field<ZZ> = <ZZ as IsComplex>::Field;

        let one_half = Angle::<ZZ12>::from(1) / Angle::<ZZ12>::from(2);
        let one_half_sqrt3 = Field::<ZZ12>::from(sqrt3::<ZZ12>()) / Field::<ZZ12>::from(2);

        // pseudo-sin behaves like sin for unit vectors
        assert_eq!(pseudo_sin(&ZZ12::unit(0)), 0.into());
        assert_eq!(pseudo_sin(&ZZ12::unit(1)), one_half);
        assert_eq!(
            Field::<ZZ12>::from(pseudo_sin(&ZZ12::unit(2))),
            one_half_sqrt3
        );
        assert_eq!(pseudo_sin(&ZZ12::unit(3)), 1.into());
        assert_eq!(pseudo_sin(&ZZ12::unit(6)), 0.into());
        assert_eq!(pseudo_sin(&ZZ12::unit(9)), (-1).into());

        // pseudo-cos behaves like cos for unit vectors
        assert_eq!(pseudo_cos(&ZZ12::unit(0)), 1.into());
        assert_eq!(
            Field::<ZZ12>::from(pseudo_cos(&ZZ12::unit(1))),
            one_half_sqrt3
        );
        assert_eq!(pseudo_cos(&ZZ12::unit(2)), one_half);
        assert_eq!(pseudo_cos(&ZZ12::unit(3)), 0.into());
        assert_eq!(pseudo_cos(&ZZ12::unit(6)), (-1).into());
        assert_eq!(pseudo_cos(&ZZ12::unit(9)), 0.into());

        // check invariant between pseudo sin/cos of different scalings of same point
        for i in 1..ZZ12::turn() {
            for j in 1..ZZ12::turn() {
                let p = ZZ12::unit(i).scale(5) + ZZ12::unit(j).scale(7);
                let p3 = p.scale(3);
                let (sp, s3p) = (pseudo_sin::<ZZ12>(&p), pseudo_sin::<ZZ12>(&p3));
                let (cp, c3p) = (pseudo_cos::<ZZ12>(&p), pseudo_cos::<ZZ12>(&p3));
                assert_eq!(s3p, sp / Angle::<ZZ12>::from(3));
                assert_eq!(c3p, cp / Angle::<ZZ12>::from(3));
            }
        }
    }

    #[test]
    fn test_pseudo_angle() {
        // check behavior on units
        let angles: Vec<Angle<ZZ12>> = (0..ZZ12::turn())
            .into_iter()
            .map(|i| pseudo_angle(&ZZ12::unit(i)))
            .collect();
        assert_eq!(angles[0], 0.into());
        assert_eq!(angles[3], 1.into());
        assert_eq!(angles[6], 2.into());
        assert_eq!(angles[9], 3.into());

        for i in 1..ZZ12::turn() {
            let ii = i as usize;
            assert!(angles[ii - 1] < angles[ii]);
            assert!(angles[ii] >= 0.into());
            assert!(angles[ii] < 4.into());
        }

        // check behavior on a non-unit that is rotated
        let angles: Vec<Angle<ZZ12>> = (0..ZZ12::turn())
            .into_iter()
            .map(|i| pseudo_angle(&((ZZ12::one().scale(2) + ZZ12::unit(1)) * ZZ12::unit(i))))
            .collect();
        for i in 1..ZZ12::turn() {
            let ii = i as usize;
            assert!(angles[ii - 1] < angles[ii]);
            assert!(angles[ii] >= 0.into());
            assert!(angles[ii] < 4.into());
        }

        // check behavior on field element
        // TODO: FIXME - make pseudoangle also work correctly for points with abs < 1
        // (compute sin/cos of the multiplicative inverse of the point instead)
        /*
        type Field<ZZ> = <ZZ as IsComplex>::Field;
        let points: Vec<_> = (0..ZZ12::turn())
            .into_iter()
            .map(|i| {
                Field::<ZZ12>::from((ZZ12::one().scale(2) + ZZ12::unit(1)) * ZZ12::unit(i))
                    / Field::<ZZ12>::from(100)
            })
            .collect();
        let angles: Vec<_> = points
            .iter()
            .map(|p: &Field<ZZ12>| pseudo_angle(p))
            .collect();
        for i in 1..ZZ12::turn() {
            let ii = i as usize;
            println!("{} -> {}", points[ii], angles[ii]);
            assert!(angles[ii - 1] < angles[ii]);
            assert!(angles[ii] >= 0.into());
            assert!(angles[ii] < 4.into());
        }
        */
    }
}
