//! Maximal-match extension from a known seed position.
//!
//! Two primitives that share the same task -- "given an anchor pair,
//! how far does the match extend?" -- but differ in what setup the
//! caller is expected to have already done:
//!
//! - [`match_length`] is a **generic two-slice matcher**: it walks
//!   forward looking for `x[i] == y[i]`, and (when the wrap-around
//!   anchor `x[0] == y[0]` also holds) attempts a cyclic
//!   left-extension. Knows nothing about revcomp / anti-parallel
//!   semantics -- it's the caller's job to set up `x` and `y` so
//!   that positive equality is the right test (typically `y =
//!   revcomp(...)`). Returns `(len, left_offset)` so the caller can
//!   recover both the matched range and how far the wrap extends.
//!
//! - [`forward_match_length`] is a **zero-alloc cyclic specialization**:
//!   takes the underlying angle arrays directly with seed positions on
//!   each side, walks forward only, with the anti-parallel index map
//!   `self[(self_start + i) % n] == -other[(other_junction + m - i) %
//!   m]` baked in. No revcomp allocation; no left-extension. Used by
//!   the patch / neighborhood paths where the seed is pinned to a
//!   junction and only forward extension matters.
//!
//! Equivalently: `forward_match_length` is what `match_length` reduces
//! to when you can guarantee the seed doesn't extend backwards and
//! you don't want to materialize the revcomp slice.

/// Length of the longest match between two equal-rooted slices,
/// plus the size of the optional cyclic left-extension.
///
/// Starting from index 0, walk forward until `x[i] != y[i]`; record
/// `len`. Then if `x[0] == y[0]`, attempt to extend cyclically to
/// the left as well (treating `x` and `y` as cyclic), returning a
/// non-zero `offset` for how far the match wraps back. The caller is
/// expected to have arranged `x` and `y` so that positive equality
/// is the right matching condition (e.g. `y = revcomp(...)` for
/// anti-parallel matching of angle sequences).
///
/// Returns `(total_len, left_offset)` where `total_len` includes both
/// the forward run and the left extension.
pub fn match_length(x: &[i8], y: &[i8]) -> (usize, usize) {
    let min_len = x.len().min(y.len());
    if min_len < 2 {
        // at least one sequence has no nodes or only single node ->
        // no segment
        return (0, 0);
    }

    // we start checking at the 2nd point and look for the right side end
    let mut len = min_len;
    for i in 1..min_len {
        if x[i] != y[i] {
            // first non-complementary angle -> right end of match
            len = i;
            break;
        }
    }
    if x[0] != y[0] {
        // cannot extend to the left -> we are done
        return (len, 0);
    }

    // go cyclically into the other direction to try extending the interval
    let remaining = min_len - len;
    if remaining == 0 {
        return (len, 0);
    }
    let mut offset = remaining;
    for i in 1..remaining {
        if x[x.len() - i] != y[y.len() - i] {
            // first non-complementary angle -> left end of match
            offset = i;
            len += i;
            return (len, offset);
        }
    }
    len += offset;
    (len, offset)
}

/// Length of the longest anti-parallel match anchored at the seed
/// vertex pair `(self_start, other_junction)`, walking forward only.
///
/// Walks outward from a shared anchor: `self` advances CCW from
/// `self_start`, `other` advances CW from `other_junction` (the
/// anti-parallel index map). Edges are compatible when
/// `self_angles[k] == -other_angles[mirror(k)]` (head-to-tail
/// anti-parallel meet).
///
/// Length is at least 1 (the anchor edge always matches by
/// construction) and stops at the first disagreement, capped at
/// `min(self_len, other_len)`. See [`crate::geom::rat::Rat::get_match`]
/// for the seed-vertex convention; see [`match_length`] for the
/// related generic variant that handles cyclic left-extension at the
/// cost of requiring a pre-built revcomp slice.
pub fn forward_match_length(
    self_angles: &[i8],
    self_start: usize,
    other_angles: &[i8],
    other_junction: usize,
) -> usize {
    let n = self_angles.len();
    let m = other_angles.len();
    let max_len = n.min(m);
    let mut len = 1;
    for i in 1..max_len {
        let xi = self_angles[(self_start + i) % n];
        let yi = -other_angles[(other_junction + m - i) % m];
        if xi != yi {
            break;
        }
        len += 1;
    }
    len
}
