use crate::traits::RealSigned;
use crate::zz::ZZNum;
use num_traits::{One, Zero};

// Core linear algebra utils
// -----------------

/// Dot product between two values, i.e. a dot b = |a||b|cos(ab).
pub fn dot<ZZ: ZZNum>(p1: &ZZ, p2: &ZZ) -> ZZ::Real {
    let (p1x, p1y) = p1.re_im();
    let (p2x, p2y) = p2.re_im();
    (p1x * p2x) + (p1y * p2y)
}

/// 2D "cross" product (wedge is the general name),
/// i.e. a wedge b = |a|b|sin(ab).
pub fn wedge<ZZ: ZZNum>(p1: &ZZ, p2: &ZZ) -> ZZ::Real {
    let (p1x, p1y) = p1.re_im();
    let (p2x, p2y) = p2.re_im();
    (p1x * p2y) - (p1y * p2x)
}

/// norm of a value, i.e. dot product of the value with itself.
pub fn norm<ZZ: ZZNum>(p: &ZZ) -> ZZ::Real {
    dot::<ZZ>(&p, &p)
}

/// Returns 1 if |self| < |other|, 0 on equality and -1 otherwise.
pub fn abs_signum<ZZ: ZZNum>(p1: &ZZ, p2: &ZZ) -> ZZ::Real {
    (norm(p2) - norm(p1)).re_signum()
}

/// Sign of the relative ccw-positive angle between this and the other number.
///
/// Result is always real-valued (Zn).
pub fn angle_signum<ZZ: ZZNum>(p1: &ZZ, p2: &ZZ) -> ZZ::Real {
    wedge(p1, p2).re_signum()
}

/// Return true if angle is in closed interval [a,b],
/// assuming that a and b are in counterclockwise order and their ccw angle
/// is less than a half turn.
pub fn angle_between<ZZ: ZZNum>(p: &ZZ, (a, b): (&ZZ, &ZZ)) -> bool {
    angle_signum(a, p) != -ZZ::Real::one() && angle_signum(p, b) != -ZZ::Real::one()
}

/// Return whether this point is strictly between the other two.
///
/// NOTE: we already assume all three involved points are colinear.
pub fn is_between<ZZ: ZZNum>(p: &ZZ, (a, b): (&ZZ, &ZZ)) -> bool {
    let v = *a - *p;
    let w = *p - *b;
    wedge(&v, &w).is_zero() && dot(&v, &w).re_signum().is_one()
}

/// Return whether the segments ab and ac (in that order) have a ccw angle,
/// i.e. c is ccw of b with respect to rotation around a.
pub fn is_ccw<ZZ: ZZNum>(p: &ZZ, (a, b): (&ZZ, &ZZ)) -> bool {
    angle_signum(&(*a - *p), &(*b - *p)).is_one()
    // NOTE: wedge(*a - *p, *p - *b).is_zero() is subexp. of is_between... interesting symmetry
}

// 2D geometry utils
// -----------------

/// Get (a, b, c) for line ax + by + c = 0 based on two points.
/// y(x) = -(ax + c)/b = -a/b x - c/b = mx + b' with m = -a/b and b' = -c/b
fn line_through<ZZ: ZZNum>(p1: &ZZ, p2: &ZZ) -> (ZZ::Real, ZZ::Real, ZZ::Real) {
    // (wedge(&ZZ::one(), &(*p1 - *p2)), dot(&ZZ::one(), &(*p2 - *p1)), wedge(&p1, &p2))
    ((*p1 - *p2).im(), (*p2 - *p1).re(), wedge::<ZZ>(&p1, &p2))
}

/// Return whether the point is on the given line.
fn is_colinear<ZZ: ZZNum>(p: &ZZ, (a, b, c): &(ZZ::Real, ZZ::Real, ZZ::Real)) -> bool {
    // ((*a) * dot(&ZZ::one(), p) + (*b) * wedge(&ZZ::one(), p) + *c).is_zero()
    (*a * p.re() + *b * p.im() + *c).is_zero()
}

/// Return whether line segments AB and CD intersect.
/// Note that touching in only endpoints does not count as intersection.
/// Based on: https://stackoverflow.com/a/9997374/432908
pub fn intersect<ZZ: ZZNum>(&(a, b): &(ZZ, ZZ), &(c, d): &(ZZ, ZZ)) -> bool {
    if a == c || a == d || b == c || b == d {
        // we ignore touching endpoints
        // (not counting it as proper intersection of _disjoint_ segments)
        // NOTE: if true, there also cannot be any other intersection
        return false;
    }
    let l_ab = line_through(&a, &b);
    if is_colinear(&c, &l_ab) && is_colinear(&d, &l_ab) {
        is_between(&c, (&a, &b)) || is_between(&d, (&a, &b))
    } else {
        // TODO: understand logic why this works
        is_ccw(&a, (&c, &d)) != is_ccw(&b, (&c, &d)) && is_ccw(&a, (&b, &c)) != is_ccw(&a, (&b, &d))
    }
}

// TODO: test and use
pub fn point_in_rect<ZZ: ZZNum>(p: &ZZ, (pos_min, pos_max): (&ZZ, &ZZ)) -> bool {
    let (px, py) = p.re_im();
    let (x_min, y_min) = pos_min.re_im();
    let (x_max, y_max) = pos_max.re_im();
    let x_in_open_bound = (px - x_min).re_signum().is_one() && (x_max - px).re_signum().is_one();
    let y_in_open_bound = (py - y_min).re_signum().is_one() && (y_max - py).re_signum().is_one();
    x_in_open_bound && y_in_open_bound
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::zz::ZZBase;
    use crate::zz::{Z12, ZZ12};
    use num_traits::{One, Zero};

    type ZZi = ZZ12;

    #[test]
    fn test_dot() {
        // get a non-trivial point

        let mut p = ZZi::zero();
        for i in 0..ZZi::turn() {
            p = p + ZZi::unit(i).scale(i as i64);
        }
        //p = zz_units_sum();

        // all rotations of the same point around origin
        // have the same squared distance, i.e. quadrance
        // and it is real-valued.
        let q: ZZi = norm(&p).into();
        for i in 1..ZZi::turn() {
            let pi = p * ZZi::unit(i);
            let qi: ZZi = norm(&pi).into();
            assert_eq!(qi, q);
        }
    }

    #[test]
    fn test_colinear() {
        if !ZZi::has_qturn() {
            return;
        }

        // Test point layout:
        // -------
        // E F
        //
        // A B C D
        // -------
        let a: ZZi = ZZi::zero();
        let b: ZZi = ZZi::one();
        let c: ZZi = ZZi::from(2);
        let d: ZZi = ZZi::from(3);
        let e: ZZi = ZZi::one_i();
        let f: ZZi = b + e;

        let l_ab = line_through(&a, &b);
        let l_ac = line_through(&a, &c);
        let l_af = line_through(&a, &f);

        // colinear, overlap
        assert!(is_colinear(&b, &l_ac));
        assert!(is_colinear(&d, &l_ac));
        // colinear, no overlap
        assert!(is_colinear(&c, &l_ab));
        assert!(is_colinear(&d, &l_ab));
        // parallel (not colinear)
        assert!(!is_colinear(&e, &l_ab));
        assert!(!is_colinear(&f, &l_ab));
        // perpendicular (touches in one point, not in the other)
        assert!(!(is_colinear(&a, &l_ab) && is_colinear(&e, &l_ab)));
        // general case
        assert!(!is_colinear(&b, &l_af));
        assert!(!is_colinear(&d, &l_af));
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
        let d1 = ZZ12::from(dot(&p60, &pi).pow(2)).complex64();
        assert_eq!(d1.re, 0.75);
        assert_eq!(d1.im, 0.0);

        // same but with negative sign (check indirectly)
        let d2 = ZZ12::from(dot(&pm60, &pi)).complex64();
        assert!(d2.re < 0.0);
        assert_eq!(d2.im, 0.0);
        assert_eq!(ZZ12::from(dot(&pm60, &pi).pow(2)).complex64().re, 0.75);
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
    fn test_intersect() {
        let a: ZZ12 = ZZ12::zero();
        let b: ZZ12 = ZZ12::one();
        let c: ZZ12 = ZZ12::from(2);
        let d: ZZ12 = ZZ12::from(3);
        let e: ZZ12 = ZZ12::one_i();
        let f: ZZ12 = b + e;

        // colinear cases:
        // ----
        assert!(!intersect(&(a, b), &(c, d))); // no touch
        assert!(!intersect(&(d, c), &(a, b))); // same, permutated
        assert!(!intersect(&(a, b), &(b, c))); // touch in an endpoint
        assert!(intersect(&(a, c), &(b, d))); // overlap
        assert!(intersect(&(a, d), &(b, c))); // overlap (subsuming)

        // non-colinear cases:
        // ----
        // parallel
        assert!(!intersect(&(a, b), &(e, f)));
        assert!(!intersect(&(a, e), &(b, f)));

        // no touch
        assert!(!intersect(&(a, e), &(b, c))); // perp
        assert!(!intersect(&(a, f), &(b, c))); // non-perp

        // touch in start/end-points
        assert!(!intersect(&(a, e), &(a, b))); // perp
        assert!(!intersect(&(a, f), &(a, b))); // non-perp
        assert!(!intersect(&(a, e), &(e, e - b))); // non-perp

        // endpoint of one segment intersects a non-endpoint of other seg
        assert!(intersect(&(b, f), &(a, c))); // perp
        assert!(intersect(&(b, e), &(a, c))); // non-perp

        // proper intersection of open segments
        assert!(intersect(&(a, f), &(b, e))); // perp
        assert!(intersect(&(a, f), &(c, e))); // non-perp
        assert!(intersect(&(a, f), &(d, e))); // non-perp
    }

    #[test]
    fn test_point_in_rect() {
        use crate::traits::Ccw;

        let min: ZZ12 = 0.into();
        let max: ZZ12 = (2, 2).into();

        // inside
        assert!(point_in_rect(&(1, 1).into(), (&min, &max)));
        assert!(point_in_rect(&ZZ12::ccw(), (&min, &max)));

        // boundary does not count
        assert!(!point_in_rect(&min, (&min, &max)));
        assert!(!point_in_rect(&max, (&min, &max)));

        // outside
        assert!(!point_in_rect(&(3, 1).into(), (&min, &max)));
    }
}
