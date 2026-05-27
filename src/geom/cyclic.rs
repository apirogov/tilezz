//! Cyclic-index arithmetic on a boundary of length `n`.
//!
//! Tiny standalone helpers -- no types, no rings, no references to
//! patches or rats. Used by [`crate::geom::patch`] (for "is this
//! match incident with a given vertex / edge range?") and elsewhere.
//! Kept separate so the patch module can stay focused on the live
//! state machine.

/// True iff vertex `index` is touched by a match of `len` edges
/// starting at edge position `start` on a cyclic boundary of length
/// `n`.
///
/// A match of `len` edges starting at edge `start` covers edges
/// `[start, start+1, ..., start+len-1]` and therefore touches the
/// `len + 1` boundary vertices `[start, start+1, ..., start+len]`
/// (the CW endpoint of the first edge through the CCW endpoint of
/// the last edge). Inclusive on both ends -- the CW vertex (`start`)
/// and the CCW vertex (`start + len`) are both counted as touched.
///
/// The forward-cyclic-distance formulation is correct regardless of
/// whether `start + len` wraps around the boundary; a naive
/// `end <= n` / `index <= end` split has an off-by-one bug exactly
/// when `start + len == n` (the wrap vertex at `0` is missed).
pub fn cyclic_range_contains(start: usize, len: usize, index: usize, n: usize) -> bool {
    if len == 0 || n == 0 {
        return false;
    }
    let cyclic_diff = (index + n - start % n) % n;
    cyclic_diff <= len
}

/// True iff two **edge-inclusive** cyclic arcs on a boundary of
/// length `n` share at least one edge.
///
/// Arc `A` = edges `[a, a+1, ..., a+l_a - 1]` (mod `n`); arc `B` =
/// edges `[b, b+1, ..., b+l_b - 1]` (mod `n`). Either arc may wrap.
/// Empty arcs (`l_a == 0` or `l_b == 0`) never overlap.
pub fn cyclic_arcs_overlap(a: usize, l_a: usize, b: usize, l_b: usize, n: usize) -> bool {
    if l_a == 0 || l_b == 0 || n == 0 {
        return false;
    }
    // Two cyclic arcs overlap iff either's start is "in" the other.
    // "Start of B in A": forward cyclic distance from `a` to `b` is
    // strictly less than `l_a`. Symmetric for "start of A in B".
    let b_in_a = (b + n - a % n) % n < l_a;
    let a_in_b = (a + n - b % n) % n < l_b;
    b_in_a || a_in_b
}
