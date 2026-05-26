//! Core linear algebra utils
use super::traits::IsRingOrField;

/// Sign of the wedge product of `p1` and `p2`.
///
/// By the identity `wedge(p1, p2) = Im(conj(p1) * p2)`, this is the same as
/// the sign of the imaginary component of `conj(p1) * p2`. Returns -1, 0, or 1.
///
/// Stays inside `ZZ` for the multiplication and only extracts a sign at the end.
pub fn wedge_sign<ZZ: IsRingOrField>(p1: &ZZ, p2: &ZZ) -> i8 {
    (p1.conj() * *p2).im_sign()
}

/// Sign of the dot product of `p1` and `p2`.
///
/// By the identity `dot(p1, p2) = Re(conj(p1) * p2)`, this is the same as
/// the sign of the real component of `conj(p1) * p2`. Returns -1, 0, or 1.
///
/// Stays inside `ZZ` for the multiplication and only extracts a sign at the end.
pub fn dot_sign<ZZ: IsRingOrField>(p1: &ZZ, p2: &ZZ) -> i8 {
    (p1.conj() * *p2).re_sign()
}

/// Return true if angle (wrt. positive real line) is in closed interval `[a,b]`,
/// assuming that a and b are in counterclockwise order and their ccw angle
/// is less than a half turn.
pub fn angle_between<ZZ: IsRingOrField>(p: &ZZ, (a, b): (&ZZ, &ZZ)) -> bool {
    wedge_sign(a, p) >= 0 && wedge_sign(p, b) >= 0
}

/// Return whether this point is strictly between the other two.
///
/// NOTE: we already assume all three involved points are colinear.
pub fn is_between<ZZ: IsRingOrField>(p: &ZZ, (a, b): (&ZZ, &ZZ)) -> bool {
    // wedge(a - p, p - b) == 0 AND dot(a - p, p - b) > 0
    let v = *a - *p;
    let w = *p - *b;
    wedge_sign(&v, &w) == 0 && dot_sign(&v, &w) > 0
}

/// Return whether the segments pa and pb (in that order) have a ccw angle,
/// i.e. b is ccw of a with respect to rotation around p.
pub fn is_ccw<ZZ: IsRingOrField>(p: &ZZ, (a, b): (&ZZ, &ZZ)) -> bool {
    // wedge(a - p, b - p) > 0
    wedge_sign(&(*a - *p), &(*b - *p)) > 0
    // NOTE: wedge(*a - *p, *p - *b).is_zero() is subexp. of is_between... interesting symmetry
}

#[cfg(test)]
mod tests {
    use num_traits::{One, Zero};

    use super::super::rings::ZZ12;
    use super::super::traits::{SymNum, Units};
    use super::*;

    type ZZi = ZZ12;

    #[test]
    fn test_is_between() {
        let a: ZZi = ZZi::zero();
        let b: ZZi = ZZi::one();
        let c: ZZi = ZZi::from(2);
        let e: ZZi = <ZZi as Units>::unit(ZZi::hturn() / 2);
        let f: ZZi = b + e;
        let g: ZZi = <ZZi as Units>::unit(1) + <ZZi as Units>::unit(-1) - ZZi::one();

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
        assert!(is_ccw(
            &ZZ::zero(),
            (&<ZZ as Units>::unit(0), &<ZZ as Units>::unit(1))
        ));
        assert!(is_ccw(
            &ZZ::zero(),
            (&<ZZ as Units>::unit(2), &<ZZ as Units>::unit(3))
        ));
        assert!(!is_ccw(
            &ZZ::zero(),
            (&<ZZ as Units>::unit(0), &<ZZ as Units>::unit(0))
        ));
        assert!(!is_ccw(
            &ZZ::zero(),
            (&<ZZ as Units>::unit(1), &<ZZ as Units>::unit(1))
        ));
        assert!(!is_ccw(
            &ZZ::zero(),
            (&<ZZ as Units>::unit(1), &<ZZ as Units>::unit(0))
        ));
        assert!(!is_ccw(
            &ZZ::zero(),
            (&<ZZ as Units>::unit(6), &<ZZ as Units>::unit(5))
        ));
    }

    #[test]
    fn test_angle_between() {
        type ZZ = ZZ12;
        // test inside
        assert!(angle_between(
            &<ZZ as Units>::unit(1),
            (&<ZZ as Units>::unit(0), &<ZZ as Units>::unit(2))
        ));
        assert!(angle_between(
            &<ZZ as Units>::unit(-2),
            (&<ZZ as Units>::unit(-3), &<ZZ as Units>::unit(1))
        ));
        // test outside
        assert!(!angle_between(
            &<ZZ as Units>::unit(1),
            (&<ZZ as Units>::unit(-1), &<ZZ as Units>::unit(0))
        ));
        // bound inclusive
        assert!(angle_between(
            &<ZZ as Units>::unit(1),
            (&<ZZ as Units>::unit(1), &<ZZ as Units>::unit(2))
        ));
        assert!(angle_between(
            &<ZZ as Units>::unit(2),
            (&<ZZ as Units>::unit(1), &<ZZ as Units>::unit(2))
        ));
    }
}
