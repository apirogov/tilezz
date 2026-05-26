//! Boundary-glue operation: given two angle sequences and a known-valid
//! match interval between them, produce the angle sequence of the
//! combined boundary.
//!
//! This is the *operation* side of gluing. The data shapes describing
//! "where to glue" (positional spec, tile ids, etc.) live in
//! [`super::matches`]; pure-sequence angle utilities (revcomp, abs/rel
//! conversion, normalization) live in [`super::angles`].

use crate::cyclotomic::IsRing;

use super::angles::normalize_angle;

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

/// Compute the angle sequence resulting from gluing two boundary segments
/// along a **known-valid** match interval.
///
/// Given two cyclic angle sequences (`self_angles` and `other_angles`) and a
/// match described by `(start_a, mlen, start_b)`, computes the new boundary
/// formed by the non-matched portions of both sequences.
///
/// # Match-interval convention
///
/// The match covers `mlen` consecutive edges starting at index `start_a` on
/// `self_angles` (so edges `start_a, start_a+1, …, start_a+mlen-1` are
/// consumed, cyclically). On `other_angles` the match **ends** at index
/// `start_b` — that is, `start_b` is the first *surviving* index on
/// `other_angles` (the convention matches what [`super::rat::Rat::get_match`]
/// returns as `ext_end`). Surviving portions are concatenated and the two
/// junction angles `a_yx` / `a_xy` are written in.
///
/// # Caller contract (preconditions)
///
/// Callers must guarantee:
///
/// 1. **Valid match interval.** The `mlen - 1` interior angle pairs along
///    the match must satisfy the revcomp relation
///    `self_angles[(start_a + i) % n] == -other_angles[(start_b + m - i) % m]`
///    for every `i ∈ 1..mlen`. The canonical way to obtain such a tuple is
///    [`super::rat::Rat::get_match`] or
///    [`super::patch::forward_match_length`] /
///    [`super::patch::GrowingPatch::get_all_matches`]. Hand-constructed
///    `PatchMatch`es that don't go through one of these enumerators are
///    *not* validated by this function (or by `add_tile` downstream — the
///    geometric self-intersection check downstream of this function may
///    happen to accept a bogus interval if the resulting polyline is
///    non-self-intersecting, producing a geometrically meaningless glue).
///    A `debug_assert!` enforces the contract in test builds.
/// 2. **No ±hturn junctions.** Caller must reject results whose `a_yx` or
///    `a_xy` equals ±hturn (degenerate / pinched junction). This function
///    does not perform that check.
///
/// # Non-local validity (post-glue self-intersection)
///
/// `glue_raw_angles` only builds the new angle sequence. Whether the
/// resulting boundary is a valid simple polygon (no self-intersection,
/// hole-free) is the caller's responsibility, enforced by one of:
///
/// - [`super::snake::Snake::try_from`] on the returned `angles` (used by
///   `Rat::try_glue` and `GrowingPatch::init_from_first_add`);
/// - incremental grid-segment checks (`check_segment_clear`) used by
///   `GrowingPatch::add_tile_growing`;
/// - upstream canonical enumeration (`forward_match_length` chains in
///   `flower_petal_glue`) guarantees the result is a valid simple polygon
///   when fed back into `GrowingPatch::from_parts`.
///
/// # Returns
///
/// `Some(GlueResult)` on success; `None` if the surviving boundary would
/// be empty.
pub fn glue_raw_angles<T: IsRing>(
    self_angles: &[i8],
    other_angles: &[i8],
    start_a: usize,
    mlen: usize,
    start_b: usize,
) -> Option<GlueResult> {
    debug_assert!(
        {
            let n = self_angles.len();
            let m = other_angles.len();
            mlen == 0
                || mlen > n
                || mlen > m
                || (1..mlen)
                    .all(|i| self_angles[(start_a + i) % n] == -other_angles[(start_b + m - i) % m])
        },
        "glue_raw_angles: claimed match interval (start_a={start_a}, mlen={mlen}, \
         start_b={start_b}) violates the revcomp relation at some interior offset. \
         The caller must pass a real match — use Rat::get_match, \
         forward_match_length, or GrowingPatch::get_all_matches to obtain one."
    );
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
    } else if y_raw_len == 1 && x_raw_len > 1 {
        // Keystone glue: the petal's full perimeter is matched, so it
        // contributes no surviving boundary edges. The two normal-case
        // junction vertices (`pm.start_a` and `ccw_pos`) collapse to a
        // single new boundary vertex at `result[0]`. The merged angle
        // there absorbs both adjustments — equivalently, `ayx + axy`
        // less the first petal angle that would have lived at position
        // `seg_len_old` in the normal case:
        //     merged = ayx + axy - y
        //            = (x_first + y - hturn) + (y + x_last - hturn) - y
        //            = x_first + x_last + y - turn
        let x_first = self_angles[(start_a + mlen) % n];
        let x_last = self_angles[(start_a + mlen + seg_len_old) % n];
        let y = other_angles[start_b % m];
        let merged = normalize_angle::<T>(x_first + x_last + y - T::turn());
        if merged.abs() == T::hturn() {
            return None;
        }
        result[0] = merged;
        (Some(merged), Some(merged))
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
