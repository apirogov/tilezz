//! Utility functions to use cyclotomic rings for 2D geometry
use num_traits::Zero;
use std::ops::Range;

use super::linalg::{is_between, is_ccw, wedge};
use super::numtraits::{OneImag, ZSigned};
use super::traits::{HasZZ4, IsComplex, IsReal, IsRingOrField};

pub type Line<T> = (T, T, T);

/// Get (a, b) for equivalence class of lines line ax + by = c for all c.
///
/// This function is mostly interesting if the constant c does not need to be computed.
fn line_class_through<ZZ: IsComplex>(p1: &ZZ, p2: &ZZ) -> (ZZ::Real, ZZ::Real) {
    ((*p1 - *p2).im(), (*p2 - *p1).re())
}

/// Get (a, b, c) for line ax + by + c = 0 based on two points.
/// y(x) = -(ax + c)/b = -a/b x - c/b = mx + b' with m = -a/b and b' = -c/b
pub fn line_through<ZZ: IsComplex>(p1: &ZZ, p2: &ZZ) -> Line<ZZ::Real> {
    // (wedge(&ZZ::one(), &(*p1 - *p2)), dot(&ZZ::one(), &(*p2 - *p1)), wedge(&p1, &p2))
    let (a, b) = line_class_through::<ZZ>(p1, p2);
    (a, b, wedge::<ZZ>(p1, p2))
}

/// Return whether the point is on the given line.
pub fn is_colinear<ZZ: IsComplex>(p: &ZZ, (a, b, c): &Line<ZZ::Real>) -> bool {
    // ((*a) * dot(&ZZ::one(), p) + (*b) * wedge(&ZZ::one(), p) + *c).is_zero()
    (*a * p.re() + *b * p.im() + *c).is_zero()
}

/// Return whether two lines are parallel (which includes colinearity).
pub fn lines_parallel<ZZ: IsComplex>(l1: &Line<ZZ::Real>, l2: &Line<ZZ::Real>) -> bool {
    let ((a1, b1, _), (a2, b2, _)) = (l1, l2);
    (*a1 * *b2 - *a2 * *b1).is_zero()
}

/// Return whether two lines are perpendicular.
pub fn lines_perp<ZZ: IsComplex>(l1: &Line<ZZ::Real>, l2: &Line<ZZ::Real>) -> bool {
    let ((a1, b1, _), (a2, b2, _)) = (l1, l2);
    (*a1 * *a2 + *b1 * *b2).is_zero()
}

/// Return whether line segments AB and CD intersect.
/// Note that touching in only endpoints does not count as intersection.
/// Based on: <https://stackoverflow.com/a/9997374/432908>
pub fn intersect<ZZ: IsComplex>(&(a, b): &(ZZ, ZZ), &(c, d): &(ZZ, ZZ)) -> bool {
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

/// Return whether a point is inside a rectangle or on its boundary.
/// If strict is true, will not consider a point on a boundary as inside.
pub fn point_in_rect<ZZ: IsComplex>(p: &ZZ, (pos_min, pos_max): &(ZZ, ZZ), strict: bool) -> bool {
    let (px, py) = p.re_im();
    let (x_min, y_min) = pos_min.re_im();
    let (x_max, y_max) = pos_max.re_im();
    let x_in_open_bound = if strict {
        (px - x_min).is_positive() && (x_max - px).is_positive()
    } else {
        !(px - x_min).is_negative() && !(x_max - px).is_negative()
    };
    let y_in_open_bound = if strict {
        (py - y_min).is_positive() && (y_max - py).is_positive()
    } else {
        !(py - y_min).is_negative() && !(y_max - py).is_negative()
    };
    x_in_open_bound && y_in_open_bound
}

/// Put a number into a closed interval by repeated addition or subtraction.
pub fn mod_bound<Z: IsReal>(p: &Z, rng: Range<&Z>) -> Z {
    let (l, r) = (rng.start, rng.end);
    let w = *r - *l;
    let mut ret = *p;

    let mut was_less = false;
    while ret < *l {
        was_less = true;
        ret = ret + w;
    }
    if was_less {
        return ret;
    }
    while ret > *r {
        ret = ret - w;
    }
    ret
}

/// Return number modulo rect by repeated addition or subtraction of the x/y components.
pub fn point_mod_rect<ZZ>(p: &ZZ, (pos_min, pos_max): &(ZZ, ZZ)) -> ZZ
where
    ZZ: IsComplex + HasZZ4 + OneImag + From<<ZZ as IsRingOrField>::Real>,
{
    let (p_x, p_y) = p.re_im();
    let (x_min, y_min) = pos_min.re_im();
    let (x_max, y_max) = pos_max.re_im();
    let ret_x: <ZZ as IsRingOrField>::Real = mod_bound(&p_x, &x_min..&x_max);
    let ret_y: <ZZ as IsRingOrField>::Real = mod_bound(&p_y, &y_min..&y_max);
    ZZ::from(ret_x) + ZZ::one_i() * ZZ::from(ret_y)
}

#[cfg(test)]
mod tests {
    use num_traits::{One, Zero};

    use super::super::numtraits::{Ccw, OneImag};
    use super::super::types::ZZ12;
    use super::*;

    type ZZi = ZZ12;

    /// Return collection of test points.
    ///
    /// Test point layout:
    /// -------
    /// E F
    ///
    /// A B C D
    /// -------
    fn get_test_points() -> (ZZi, ZZi, ZZi, ZZi, ZZi, ZZi) {
        let a: ZZi = ZZi::zero();
        let b: ZZi = ZZi::one();
        let c: ZZi = ZZi::from(2);
        let d: ZZi = ZZi::from(3);
        let e: ZZi = ZZi::one_i();
        let f: ZZi = b + e;
        (a, b, c, d, e, f)
    }

    #[test]
    fn test_colinear_parallel_perp() {
        let (a, b, c, d, e, f) = get_test_points();

        let l_ab = line_through(&a, &b);
        let l_ac = line_through(&a, &c);
        let l_ae = line_through(&a, &e);
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

        assert!(lines_parallel::<ZZi>(&l_ab, &l_ac));
        assert!(!lines_parallel::<ZZi>(&l_ab, &l_af));
        assert!(lines_perp::<ZZi>(&l_ab, &l_ae));
        assert!(!lines_perp::<ZZi>(&l_ab, &l_af));
    }

    #[test]
    fn test_intersect() {
        let (a, b, c, d, e, f) = get_test_points();

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
        let rect @ (min, max): (ZZ12, ZZ12) = (0.into(), (2, 2).into());

        // inside
        assert!(point_in_rect(&(1, 1).into(), &rect, true));
        assert!(point_in_rect(&ZZ12::ccw(), &rect, true));

        // boundary
        assert!(!point_in_rect(&min, &rect, true));
        assert!(!point_in_rect(&max, &rect, true));
        assert!(point_in_rect(&min, &rect, false));
        assert!(point_in_rect(&max, &rect, false));

        // outside
        assert!(!point_in_rect(&(3, 1).into(), &rect, true));
    }
}
