//! Angle-related utilities

use crate::cyclotomic::SymNum;

/// Get complement of an angle sequence, i.e. with flipped sign.
pub fn comp(angles: &[i8]) -> Vec<i8> {
    angles.iter().map(|a| -a).collect()
}

/// Get reverse complement of an angle sequence,
/// i.e. with reversed order and sign.
pub fn revcomp(angles: &[i8]) -> Vec<i8> {
    angles.iter().rev().map(|a| -a).collect()
}

/// Convert a fraction of a turn into the corresponding unit index.
pub fn unit_idx<T: SymNum>(angle: i8) -> i8 {
    let a = normalize_angle::<T>(angle);
    let half = T::hturn();
    if a == 0 {
        return 1;
    }
    if a == half {
        return -1;
    }
    if a > 0 {
        return a + 1;
    }
    -(half + a) - 1
}

/// Convert a normalized relative polygon angle sequence to an absolute angle sequence.
///
/// The directions are in ccw order starting from the pos real line,
/// the other half-circle has the the dual opposite directions.
pub fn to_abs_seq<T: SymNum>(angles: &[i8]) -> Vec<i8> {
    let mut currdir: i8 = 0;
    let mut result: Vec<i8> = Vec::new();
    for a in angles {
        currdir += a;
        result.push(unit_idx::<T>(currdir));
    }
    result
}

/// Normalize an angle to the closed interval `[-H, H]`, where `H` is the
/// half-turn of the ring. This is used to have a unique symbolic snake
/// representation.
pub fn normalize_angle<T: SymNum>(angle: i8) -> i8 {
    let a = angle % T::turn();
    if a.abs() <= T::hturn() {
        a
    } else {
        -(a.signum() * T::turn() - a)
    }
}

/// Result of gluing two boundary angle sequences.
pub struct GlueResult {
    /// The new boundary angle sequence.
    pub angles: Vec<i8>,
    /// Index in `angles` of the CW junction (where old boundary meets new tile).
    /// `None` if one side is empty (no junction to fix).
    pub cw_junction_idx: Option<usize>,
    /// Index in `angles` of the CCW junction (where new tile meets old boundary).
    /// `None` if one side is empty (no junction to fix).
    pub ccw_junction_idx: Option<usize>,
    /// The CW junction angle (already written into `angles`).
    pub a_xy: Option<i8>,
    /// The CCW junction angle (already written into `angles`).
    pub a_yx: Option<i8>,
}

/// Compute the angle sequence resulting from gluing two boundary segments.
///
/// Given two cyclic angle sequences (`self_angles` and `other_angles`) and a
/// match described by `(start_a, mlen, start_b)`, computes the new boundary
/// formed by the non-matched portions of both sequences.
///
/// The match covers `mlen` consecutive positions starting at `start_a` on
/// `self_angles` (cyclically) and ends at `start_b` on `other_angles` (cyclically).
/// The surviving portions are concatenated with corrected junction angles.
///
/// Returns `None` if the result would be empty. Does **not** check for ±hturn
/// degeneracy — the caller is responsible for rejecting degenerate junctions.
pub fn glue_raw_angles<T: SymNum>(
    self_angles: &[i8],
    other_angles: &[i8],
    start_a: usize,
    mlen: usize,
    start_b: usize,
) -> Option<GlueResult> {
    let n = self_angles.len();
    let m = other_angles.len();

    let x_raw_len = n - mlen + 1;
    let y_raw_len = m - mlen + 1;

    if x_raw_len <= 1 && y_raw_len <= 1 {
        return None;
    }

    let seg_len_old = x_raw_len - 1;
    let seg_len_new = y_raw_len - 1;

    let mut result = Vec::with_capacity(seg_len_old + seg_len_new);
    for i in 0..seg_len_old {
        result.push(self_angles[(start_a + mlen + i) % n]);
    }
    for i in 0..seg_len_new {
        result.push(other_angles[(start_b + i) % m]);
    }

    if result.is_empty() {
        return None;
    }

    let (a_yx, a_xy) = if x_raw_len > 1 && y_raw_len > 1 {
        let x_first = self_angles[(start_a + mlen) % n];
        let x_last = self_angles[(start_a + mlen + seg_len_old) % n];
        let y_first = other_angles[start_b % m];
        let y_last = other_angles[(start_b + seg_len_new) % m];
        let ayx = normalize_angle::<T>(x_first + y_last - T::hturn());
        let axy = normalize_angle::<T>(y_first + x_last - T::hturn());
        if seg_len_old > 0 {
            result[0] = ayx;
        }
        if seg_len_old > 0 && seg_len_new > 0 {
            result[seg_len_old] = axy;
        } else if seg_len_new > 0 {
            result[0] = axy;
        }
        (Some(ayx), Some(axy))
    } else {
        (None, None)
    };

    let cw_junction_idx = if seg_len_old > 0 && seg_len_new > 0 {
        Some(seg_len_old)
    } else {
        None
    };
    let ccw_junction_idx = if seg_len_old > 0 { Some(0) } else { None };

    Some(GlueResult {
        angles: result,
        cw_junction_idx,
        ccw_junction_idx,
        a_yx,
        a_xy,
    })
}

/// Rescale the angle sequence from one complex integer ring to another.
/// Assumes that the target ring contains the original ring of the sequence.
pub fn upscale_angles<T: SymNum>(src_ring: i8, angles: &[i8]) -> Vec<i8> {
    // NOTE: using assertion here because using
    // incompatible rings here is an implementation error.
    assert_eq!(T::zz_params().full_turn_steps % src_ring, 0);

    let scale = T::zz_params().full_turn_steps / src_ring;
    angles.iter().map(|x| x * scale).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic::{SymNum, ZZ12, ZZ24, ZZ6};

    #[test]
    fn test_revcomp() {
        assert_eq!(revcomp(vec![].as_slice()), vec![]);
        assert_eq!(revcomp(vec![6].as_slice()), vec![-6]);
        let v = vec![1, 3, -2, 4];
        assert_eq!(revcomp(v.as_slice()), vec![-4, 2, -3, -1]);
        assert_eq!(revcomp(revcomp(v.as_slice()).as_slice()).as_slice(), v);
    }

    #[test]
    fn test_abs_seq() {
        assert_eq!(
            to_abs_seq::<ZZ12>(&[0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]),
            &[1, 1, 2, 3, 4, 5, 6, -1, -2, -3, -4, -5, -6, 1]
        );
        assert_eq!(
            to_abs_seq::<ZZ12>(&[0, 4, 2, 1, -3, -5]),
            &[1, 5, -1, -2, 5, -6]
        );
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
}
