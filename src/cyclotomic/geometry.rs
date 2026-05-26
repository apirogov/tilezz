//! Utility functions to use cyclotomic rings for 2D geometry
use super::linalg::{dot_sign, is_between, is_ccw, wedge_sign};
use super::traits::{HasZZ4, IntersectUnitSegments, IsRing, OneImag};

/// Return whether the point `p` lies on the line through `a` and `b`.
///
/// Uses the identity `wedge(b - a, p - a) == 0` (i.e. the imaginary part of
/// `conj(b - a) * (p - a)` is zero), so no real-subring intermediate is materialized.
pub fn is_colinear<ZZ: IsRing>(p: &ZZ, a: &ZZ, b: &ZZ) -> bool {
    wedge_sign(&(*b - *a), &(*p - *a)) == 0
}

/// Return whether the line through `(a, b)` is parallel to the line through `(c, d)`
/// (which includes colinearity).
///
/// Uses the identity `wedge(b - a, d - c) == 0`.
pub fn lines_parallel<ZZ: IsRing>(
    (a, b): (&ZZ, &ZZ),
    (c, d): (&ZZ, &ZZ),
) -> bool {
    wedge_sign(&(*b - *a), &(*d - *c)) == 0
}

/// Return whether the line through `(a, b)` is perpendicular to the line through `(c, d)`.
///
/// Uses the identity `dot(b - a, d - c) == 0`.
pub fn lines_perp<ZZ: IsRing>(
    (a, b): (&ZZ, &ZZ),
    (c, d): (&ZZ, &ZZ),
) -> bool {
    dot_sign(&(*b - *a), &(*d - *c)) == 0
}

/// Return whether line segments AB and CD intersect.
/// Note that touching in only endpoints does not count as intersection.
/// Based on: <https://stackoverflow.com/a/9997374/432908>
pub fn intersect<ZZ: IsRing + PartialEq>(
    &(a, b): &(ZZ, ZZ),
    &(c, d): &(ZZ, ZZ),
) -> bool {
    if a == c || a == d || b == c || b == d {
        // we ignore touching endpoints
        // (not counting it as proper intersection of _disjoint_ segments)
        // NOTE: if true, there also cannot be any other intersection
        return false;
    }
    if is_colinear(&c, &a, &b) && is_colinear(&d, &a, &b) {
        is_between(&c, (&a, &b)) || is_between(&d, (&a, &b))
    } else {
        // TODO: understand logic why this works
        is_ccw(&a, (&c, &d)) != is_ccw(&b, (&c, &d)) && is_ccw(&a, (&b, &c)) != is_ccw(&a, (&b, &d))
    }
}

/// Return whether two unit-length segments intersect.
///
/// **Precondition**: both `s1.1 - s1.0` and `s2.1 - s2.0` are unit vectors of
/// the ring's CCW unit group. Behaviour is unspecified if either segment is
/// not unit-length; callers that cannot guarantee unit-length should use the
/// generic [`intersect`] instead.
///
/// Touching only in endpoints does **not** count as intersection.
///
/// Per-ring overrides may exploit the unit-length precondition for speed
/// (see `IntersectUnitSegments` and the ZZ12 specialization).
#[inline]
pub fn intersect_unit_segments<ZZ: IntersectUnitSegments>(
    s1: &(ZZ, ZZ),
    s2: &(ZZ, ZZ),
) -> bool {
    ZZ::intersect_unit_segments(s1, s2)
}

/// Return the four signs that characterize where `p` sits relative to the
/// axis-aligned rectangle `[pos_min, pos_max]`:
///
/// ```text
///   [ sign(Re(p) - Re(pos_min)),    // -1 left of rect, 0 on left edge, +1 inside re bounds
///     sign(Im(p) - Im(pos_min)),    // -1 below rect,  0 on bottom edge, +1 inside im bounds
///     sign(Re(pos_max) - Re(p)),    // -1 right of rect, 0 on right edge, +1 inside re bounds
///     sign(Im(pos_max) - Im(p)) ]   // -1 above rect, 0 on top edge, +1 inside im bounds
/// ```
///
/// All sign tests stay in pure ring arithmetic via `re_sign`/`im_sign`; no
/// real-subring intermediate is materialized. The four signs are the
/// primitive containment query: every callable predicate ("is `p` strictly
/// inside", "is `p` in the closed rect", "which way is `p` outside")
/// reduces to a Boolean of these signs.
pub fn rect_signs<ZZ: IsRing>(p: &ZZ, pos_min: &ZZ, pos_max: &ZZ) -> [i8; 4] {
    let dlo = *p - *pos_min;
    let dhi = *pos_max - *p;
    [dlo.re_sign(), dlo.im_sign(), dhi.re_sign(), dhi.im_sign()]
}

/// Return whether a point is inside a rectangle or on its boundary.
/// If `strict` is true, points on a boundary do not count as inside.
pub fn point_in_rect<ZZ: IsRing>(p: &ZZ, (pos_min, pos_max): &(ZZ, ZZ), strict: bool) -> bool {
    let signs = rect_signs(p, pos_min, pos_max);
    if strict {
        signs.iter().all(|&s| s > 0)
    } else {
        signs.iter().all(|&s| s >= 0)
    }
}

/// Return number modulo rect by repeated addition or subtraction of the x/y components.
///
/// The rectangle must be a lattice rectangle: `pos_max - pos_min` must be
/// of the form `width + height*i` for integer `width, height >= 0`. This
/// is true when `pos_min` and `pos_max` are built from `From<(i64, i64)>`.
/// `LatticePair::to_lattice_pair` extracts the integer width/height in
/// exact ring arithmetic (debug-asserts the difference is a lattice
/// point).
///
/// Mods `p` into the half-open box `[pos_min, pos_max)` by repeated
/// lattice shifts of `width` along the real axis and `height` along the
/// imaginary axis, matching the prior `mod_bound`-on-closed-range
/// semantics.
pub fn point_mod_rect<ZZ>(p: &ZZ, (pos_min, pos_max): &(ZZ, ZZ)) -> ZZ
where
    ZZ: HasZZ4,
{
    // Extract integer width and height in exact ring arithmetic.
    let (width, height) = (*pos_max - *pos_min).to_lattice_pair();

    // Build lattice shift steps as ZZ values: w_step = width * 1, h_step = height * i.
    let w_step: ZZ = ZZ::from(width);
    let h_step: ZZ = ZZ::one_i() * ZZ::from(height);

    let mut ret = *p;

    // Real axis: bring ret.re into [pos_min.re, pos_max.re), matching the
    // closed-range mod_bound semantics (< 0 vs > 0 sign tests).
    let mut was_less = false;
    while (ret - *pos_min).re_sign() < 0 {
        was_less = true;
        ret = ret + w_step;
    }
    if !was_less {
        while (*pos_max - ret).re_sign() < 0 {
            ret = ret - w_step;
        }
    }

    // Imaginary axis: same shape, using h_step and im_sign.
    let mut was_less = false;
    while (ret - *pos_min).im_sign() < 0 {
        was_less = true;
        ret = ret + h_step;
    }
    if !was_less {
        while (*pos_max - ret).im_sign() < 0 {
            ret = ret - h_step;
        }
    }

    ret
}

#[cfg(test)]
mod tests {
    use num_traits::{One, Zero};

    use super::super::rings::ZZ12;
    use super::super::traits::{Ccw, OneImag, Units};
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

        // colinear, overlap
        assert!(is_colinear(&b, &a, &c));
        assert!(is_colinear(&d, &a, &c));
        // colinear, no overlap
        assert!(is_colinear(&c, &a, &b));
        assert!(is_colinear(&d, &a, &b));
        // parallel (not colinear)
        assert!(!is_colinear(&e, &a, &b));
        assert!(!is_colinear(&f, &a, &b));
        // perpendicular (touches in one point, not in the other)
        assert!(!(is_colinear(&a, &a, &b) && is_colinear(&e, &a, &b)));
        // general case
        assert!(!is_colinear(&b, &a, &f));
        assert!(!is_colinear(&d, &a, &f));

        assert!(lines_parallel::<ZZi>((&a, &b), (&a, &c)));
        assert!(!lines_parallel::<ZZi>((&a, &b), (&a, &f)));
        assert!(lines_perp::<ZZi>((&a, &b), (&a, &e)));
        assert!(!lines_perp::<ZZi>((&a, &b), (&a, &f)));
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
    fn test_intersect_unit_segments_matches_intersect() {
        // Sanity check: for every pair of unit-length segments we can build
        // from points on the ZZ12 unit grid (`+/-` real axis, +/- imag axis,
        // and a few diagonals), the specialized fast path agrees with the
        // generic `intersect`.
        let zero: ZZi = ZZi::zero();
        let one: ZZi = ZZi::one();
        let mi_one: ZZi = -one;
        let one_i: ZZi = ZZi::one_i();
        let mi_i: ZZi = -one_i;
        let one_plus_i: ZZi = one + one_i;

        // Each point with its 4 axis-unit neighbors gives unit segments.
        let centers = [zero, one, mi_one, one_i, mi_i, one_plus_i];

        let mut unit_segs: Vec<(ZZi, ZZi)> = Vec::new();
        for &c in &centers {
            for k in 0..12i8 {
                let u = <ZZi as Units>::unit(k);
                unit_segs.push((c, c + u));
            }
        }

        for s1 in &unit_segs {
            for s2 in &unit_segs {
                let general = intersect(s1, s2);
                let specialized = intersect_unit_segments(s1, s2);
                assert_eq!(
                    general, specialized,
                    "mismatch for s1={s1:?}, s2={s2:?}: general={general} specialized={specialized}"
                );
            }
        }
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
