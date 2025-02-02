//! Utility functions to use cycotomic rings and fields for 2D geometry
use num_traits::Zero;
use std::ops::Range;

use super::linalg::{is_between, is_ccw, wedge};
use super::numtraits::ZSigned;
use super::traits::{HasZZ4, IsComplex, IsReal, IsRingOrField};

// -----------------

/// Get (a, b, c) for line ax + by + c = 0 based on two points.
/// y(x) = -(ax + c)/b = -a/b x - c/b = mx + b' with m = -a/b and b' = -c/b
fn line_through<ZZ: IsRingOrField + IsComplex>(p1: &ZZ, p2: &ZZ) -> (ZZ::Real, ZZ::Real, ZZ::Real) {
    // (wedge(&ZZ::one(), &(*p1 - *p2)), dot(&ZZ::one(), &(*p2 - *p1)), wedge(&p1, &p2))
    ((*p1 - *p2).im(), (*p2 - *p1).re(), wedge::<ZZ>(&p1, &p2))
}

/// Return whether the point is on the given line.
fn is_colinear<ZZ: IsRingOrField + IsComplex>(
    p: &ZZ,
    (a, b, c): &(ZZ::Real, ZZ::Real, ZZ::Real),
) -> bool {
    // ((*a) * dot(&ZZ::one(), p) + (*b) * wedge(&ZZ::one(), p) + *c).is_zero()
    (*a * p.re() + *b * p.im() + *c).is_zero()
}

/// Return whether line segments AB and CD intersect.
/// Note that touching in only endpoints does not count as intersection.
/// Based on: <https://stackoverflow.com/a/9997374/432908>
pub fn intersect<ZZ: IsRingOrField + IsComplex>(&(a, b): &(ZZ, ZZ), &(c, d): &(ZZ, ZZ)) -> bool {
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
pub fn point_in_rect<ZZ: IsRingOrField + IsComplex>(
    p: &ZZ,
    (pos_min, pos_max): &(ZZ, ZZ),
    strict: bool,
) -> bool {
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
pub fn point_mod_rect<ZZ: IsComplex + HasZZ4>(
    p: &ZZ,
    (pos_min, pos_max): &(ZZ, ZZ),
) -> <ZZ as IsComplex>::Field
where
    // FIXME: figure out how to make this extra constraint not necessary
    <ZZ as IsComplex>::Field: From<(<ZZ as IsRingOrField>::Real, <ZZ as IsRingOrField>::Real)>,
{
    let (p_x, p_y) = p.re_im();
    let (x_min, y_min) = pos_min.re_im();
    let (x_max, y_max) = pos_max.re_im();
    let ret_x: <ZZ as IsRingOrField>::Real = mod_bound(&p_x, &x_min..&x_max).into();
    let ret_y: <ZZ as IsRingOrField>::Real = mod_bound(&p_y, &y_min..&y_max).into();
    (ret_x, ret_y).into()
}

#[cfg(test)]
mod tests {
    use num_traits::{One, Zero};

    use super::super::numtraits::{Ccw, OneImag};
    use super::super::symnum::SymNum;
    use super::super::types::ZZ12;
    use super::*;

    type ZZi = ZZ12;

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
