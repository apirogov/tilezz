use super::zzbase::ZZNum;
use num_traits::One;

// Misc utils
// ----------

/// Get reverse complement of an angle sequence,
/// i.e. with reversed order and sign.
pub fn revcomp(angles: &[i8]) -> Vec<i8> {
    angles.iter().rev().map(|a| -a).collect()
}

/// Normalize an angle to (-half-turn, half-turn). This is
/// used to have a unique symbolic snake representation.
pub fn normalize_angle<T: ZZNum>(angle: i8) -> i8 {
    let a = angle % T::turn();
    if a.abs() <= T::hturn() {
        a
    } else {
        -(a.signum() * T::turn() - a)
    }
}

/// Rescale the angle sequence from one complex integer ring to another.
/// Assumes that the target ring contains the original ring of the sequence.
pub fn upscale_angles<T: ZZNum>(src_ring: i8, angles: &[i8]) -> Vec<i8> {
    // NOTE: using assertion here because using
    // incompatible rings here is an implementation error.
    assert_eq!(T::zz_params().full_turn_steps % src_ring, 0);

    let scale = T::zz_params().full_turn_steps / src_ring;
    angles.iter().map(|x| x * scale).collect()
}

// 2D geometry utils
// -----------------

fn is_ccw<T: ZZNum>(a: T, b: T, c: T) -> bool {
    (b - a).wedge(&(c - a)).re_signum().is_one()
}

/// Return whether line segments AB and CD intersect.
/// Note that touching in only endpoints does not count as intersection.
/// Based on: https://stackoverflow.com/a/9997374/432908
pub fn intersect<T: ZZNum>(&(a, b): &(T, T), &(c, d): &(T, T)) -> bool {
    let l_ab = a.line_through(&b);
    if c.is_colinear(&l_ab) && d.is_colinear(&l_ab) {
        c.is_between(&a, &b) || d.is_between(&a, &b)
    } else {
        if a == c || a == d || b == c || b == d {
            false // ignore touching endpoints in non-colinear case
        } else {
            is_ccw(a, c, d) != is_ccw(b, c, d) && is_ccw(a, b, c) != is_ccw(a, b, d)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::zz::{ZZ12, ZZ24, ZZ6};
    use crate::zzbase::ZZBase;
    use num_traits::{One, Zero};

    #[test]
    fn test_revcomp() {
        assert_eq!(revcomp(vec![].as_slice()), vec![]);
        assert_eq!(revcomp(vec![6].as_slice()), vec![-6]);
        let v = vec![1, 3, -2, 4];
        assert_eq!(revcomp(v.as_slice()), vec![-4, 2, -3, -1]);
        assert_eq!(revcomp(revcomp(v.as_slice()).as_slice()).as_slice(), v);
    }

    #[test]
    fn test_normalize_angles() {
        assert_eq!(normalize_angle::<ZZ12>(0), 0);
        assert_eq!(normalize_angle::<ZZ12>(1), 1);
        assert_eq!(normalize_angle::<ZZ12>(-1), -1);
        assert_eq!(normalize_angle::<ZZ12>(5), 5);
        assert_eq!(normalize_angle::<ZZ12>(-5), -5);
        assert_eq!(normalize_angle::<ZZ12>(6), 6);
        assert_eq!(normalize_angle::<ZZ12>(-6), -6);
        assert_eq!(normalize_angle::<ZZ12>(7), -5);
        assert_eq!(normalize_angle::<ZZ12>(-7), 5);
        assert_eq!(normalize_angle::<ZZ12>(11), -1);
        assert_eq!(normalize_angle::<ZZ12>(-11), 1);
        assert_eq!(normalize_angle::<ZZ12>(12), 0);
        assert_eq!(normalize_angle::<ZZ12>(-12), 0);
        assert_eq!(normalize_angle::<ZZ12>(13), 1);
        assert_eq!(normalize_angle::<ZZ12>(-13), -1);
    }

    #[test]
    fn test_upscale_angles() {
        let exp: &[i8] = &[-4, 8, 12];
        assert_eq!(upscale_angles::<ZZ24>(ZZ6::turn(), &[-1, 2, 3]), exp);
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
}
