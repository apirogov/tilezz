//! Boundary-glue operation: given two angle sequences and a known-valid
//! match interval between them, produce the angle sequence of the
//! combined boundary.
//!
//! This is the *operation* side of gluing. The data shapes describing
//! "where to glue" (positional spec, tile ids, etc.) live in
//! [`crate::matches`]; pure-sequence angle utilities (revcomp, abs/rel
//! conversion, normalization) live in [`super::angles`].
//!
//! ## Match-interval convention
//!
//! Every glue operation here uses the asymmetric `(ns, mlen, ne)`
//! convention: the match consumes A's edges `[ns, ns + mlen)` cyclically,
//! and on B the match **ends** at index `ne` — i.e. `ne` is the first
//! *surviving* edge on B. Equivalently, the consumed B-edges run
//! `[ne - mlen, ne)` cyclically (walked in reverse since B is glued
//! anti-parallel to A). This matches what `Rat::get_match` returns and
//! what [`crate::matches::TileMatch`] / [`crate::matches::PatchMatch`]
//! store.
//!
//! ## Junction angles
//!
//! After a glue, the new boundary has (up to) two junction vertices
//! where the surviving A-edges meet the surviving B-edges:
//!
//! - **CW junction** at position `ns` (where surviving B meets
//!   surviving A on the right side of the match).
//! - **CCW junction** at position `(ns + mlen) mod n` (where
//!   surviving A meets surviving B on the left side).
//!
//! Each junction's "raw angle sum" is the sum of the two boundary
//! angles meeting there; the actual junction angle after gluing is
//! `sum - hturn` (normalized to `[-hturn, hturn]`). [`junction_sums`]
//! returns both sums; [`junctions_glueable`]
//! are cheap pre-filters.

use crate::cyclotomic::IsRing;

use super::angles::normalize_angle;

// ============================================================
// Junction-validity primitives (cheap pre-filters).
// ============================================================

/// Raw angle sums at the two junction vertices of a candidate match.
///
/// Returns `(cw, ccw)` where:
///
///   `cw  = a[ns]              + b[ne]                  `
///   `ccw = a[(ns + mlen) % n] + b[(ne + m - mlen) % m] `
///
/// The actual signed junction angle after gluing is `sum - hturn`,
/// normalized to `[-hturn, hturn]`. Each sum has three regimes with
/// distinct geometric meanings:
///
/// - `sum >  0`: the junction is a **valid, non-degenerate corner** of
///   a glue that is **maximal** at this endpoint (cannot be extended).
/// - `sum == 0`: the revcomp identity `a[idx] == -b[idx']` holds at
///   this endpoint — the match could be **extended past this endpoint**
///   into a longer match. The current `mlen` is *not maximal* on this side.
/// - `sum <  0`: the would-be junction angle is `< -hturn`. After
///   normalization the boundary would **fold back / self-intersect** at
///   this junction; no valid simple-polygon glue is possible.
pub fn junction_sums(a: &[i8], ns: usize, mlen: usize, b: &[i8], ne: usize) -> (i32, i32) {
    let n = a.len();
    let m = b.len();
    let cw = a[ns] as i32 + b[ne] as i32;
    let ccw = a[(ns + mlen) % n] as i32 + b[(ne + m - mlen) % m] as i32;
    (cw, ccw)
}

/// Pre-filter: the candidate match's two new junctions are
/// **glueable** — valid corners of a **maximal** glue at this endpoint.
///
/// Specifically, both `cw` and `ccw` from [`junction_sums`] are strictly
/// positive. The strict `> 0` simultaneously rejects two distinct bad
/// cases at each junction:
///
/// 1. `sum < 0`: the glue would self-intersect (junction angle below
///    fold-back). The full glue cannot exist.
/// 2. `sum == 0`: this candidate is **not maximal** at this endpoint —
///    a longer match exists with the same edges. Whoever is enumerating
///    matches should use that longer one instead of this candidate.
///
/// Necessary but not sufficient: the full glue still has to pass
/// [`glue_raw_angles`]'s `±hturn` check at each junction (which catches
/// the upper-end degeneracy `sum == 2*hturn`) and the downstream
/// self-intersection check.
///
/// Use this filter:
///
/// - **Before** `get_match`, when brute-force-enumerating single-edge
///   matches by sweeping `(ia, ib)` pairs (e.g.
///   `matchtypes::MatchFinder::candidates_for_pair`): rejects pairs that
///   are non-maximal as length-1 matches.
/// - **After** `get_match`, for multi-edge matches: maximality is
///   already guaranteed by `get_match`'s construction, so the `sum == 0`
///   case won't fire in practice; the check effectively reduces to
///   "no self-intersection at the junctions".
pub fn junctions_glueable(a: &[i8], ns: usize, mlen: usize, b: &[i8], ne: usize) -> bool {
    let (cw, ccw) = junction_sums(a, ns, mlen, b, ne);
    cw > 0 && ccw > 0
}

// ============================================================
// Glue operation.
// ============================================================

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

/// One junction angle: `(angle_a + angle_b - hturn)` normalised to the
/// canonical `[-hturn, hturn]` range. The two arguments are the two
/// boundary angles meeting at the new junction (one from each side).
#[inline]
fn junction_angle<T: IsRing>(angle_a: i8, angle_b: i8) -> i8 {
    normalize_angle::<T>(angle_a + angle_b - T::hturn())
}

/// Compute the angle sequence resulting from gluing two boundary
/// segments along a **known-valid** match interval.
///
/// This is a pure algebraic transform on angle sequences: it does not
/// check or know about planar geometry, simplicity, or self-intersection
/// of the resulting boundary. The caller is responsible for ensuring
/// the input is geometrically sensible and for validating the output
/// (see "Non-local validity" below).
///
/// # Output shape per case
///
/// Three branches, distinguished by which side(s) retain surviving edges:
///
/// ## Normal glue (`surv_a > 0 && surv_b > 0`)
///
/// Both sides contribute. The new boundary is `[A-survivors,
/// B-survivors]`. Two distinct junction vertices are created -- `ccw`
/// (`a_yx`) at position `0` and `cw` (`a_xy`) at position `surv_a`.
/// This is the configuration of every "normal" tile-onto-patch glue.
///
/// ## Keystone glue (`surv_b == 0`, `mlen == m`)
///
/// All of B's edges are matched; B contributes zero surviving edges. **
/// This is a hole-filling operation, not a normal glue.** The
/// pre-glue patch A is *weakly simple*: its boundary has a tile-shaped
/// notch enclosing an empty cavity, with the notch's CW and CCW
/// endpoints coinciding at a single pinch-point vertex. B is the tile
/// that fits exactly into the cavity; placing it heals the pinch into
/// a strictly simple polygon.
///
/// Geometrically the CW and CCW junctions of the normal case collapse
/// into one **seam vertex** at position `0` whose merged angle absorbs
/// both junction contributions plus the B-side angle that would have
/// lived between them:
///
/// ```text
/// merged = a_yx + a_xy - y_first
///        = (x_first + y_first - hturn) + (y_first + x_last - hturn) - y_first
///        = x_first + x_last + y_first - turn
/// ```
///
/// Returns `None` if `merged.abs() == hturn` (degenerate pinched
/// junction).
///
/// The only production caller is [`crate::geom::patch::flower_petal_glue`]
/// during Penrose-style flower constructions where a central tile +
/// petals leave a tile-shaped cavity that the keystone tile closes.
/// `Snake` would reject the weakly-simple intermediate state; this is
/// why the keystone path operates at the `RawBoundary` level rather
/// than going through Snake validation at every step.
///
/// ## Degenerate / unreachable cases
///
/// - `surv_a == 0 && surv_b > 0` (**A-keystone**): mirror image of the
///   keystone case with A and B swapped -- A is the tile filling a
///   tile-shaped cavity in B's boundary. By convention every caller
///   has A = patch (large) and B = single tile, so this never fires;
///   `glue_raw_angles` panics via `unreachable!`. (Symmetric merged
///   formula could be added if ever needed.)
/// - `surv_a == 0 && surv_b == 0` (**both consumed**, `mlen == n == m`):
///   would require A and B to be the same tile shape glued along
///   their entire perimeters, producing a topological sphere with no
///   planar polygon realisation. Panics.
///
/// # Caller contract (preconditions)
///
/// 1. **Valid match interval.** The `mlen - 1` interior angle pairs
///    along the match must satisfy the revcomp relation
///    `a[(ns + i) % n] == -b[(ne + m - i) % m]` for every
///    `i ∈ 1..mlen`. The canonical way to obtain such a tuple is
///    [`crate::geom::rat::Rat::get_match`] or
///    [`crate::geom::patch::forward_match_length`] /
///    [`crate::geom::patch::GrowingPatch::get_all_matches`].
///    A `debug_assert!` enforces this in test builds; in release,
///    `glue_raw_angles` is unchecked and may produce a geometrically
///    meaningless boundary if fed a bogus interval.
///
/// 2. **No `±hturn` junctions.** Callers must reject results whose
///    `a_yx` or `a_xy` equals `±hturn` (degenerate pinched junction).
///    `glue_raw_angles` only rejects this internally in the
///    keystone-glue branch (where the two junctions collapse to one
///    merged angle, and that single angle is checked).
///
/// # Non-local validity (post-glue self-intersection)
///
/// `glue_raw_angles` only builds the new angle sequence. Whether the
/// resulting boundary is a valid simple polygon (no self-intersection,
/// hole-free) is the caller's responsibility, enforced by one of:
///
/// - [`crate::geom::snake::Snake::try_from`] on the returned `angles`
///   (used by `Rat::try_glue` and `GrowingPatch::init_from_first_add`);
/// - incremental grid-segment checks (`check_segment_clear`) used by
///   `GrowingPatch::add_tile_growing`;
/// - upstream canonical enumeration (`forward_match_length` chains in
///   `flower_petal_glue`) which guarantees the result is a valid simple
///   polygon when fed back into `GrowingPatch::from_parts`.
///
/// # Returns
///
/// `Some(GlueResult)` on success; `None` only in the keystone-glue case
/// if the merged junction angle would be the degenerate `±hturn`.
/// Other unreachable configurations panic rather than return `None`,
/// so a `None` here unambiguously means "valid keystone match, but the
/// merged angle pinches".
pub fn glue_raw_angles<T: IsRing>(
    self_angles: &[i8],
    other_angles: &[i8],
    start_a: usize,
    mlen: usize,
    start_b: usize,
) -> Option<GlueResult> {
    let n = self_angles.len();
    let m = other_angles.len();

    debug_assert!(
        mlen == 0
            || mlen > n
            || mlen > m
            || (1..mlen)
                .all(|i| self_angles[(start_a + i) % n] == -other_angles[(start_b + m - i) % m]),
        "glue_raw_angles: claimed match interval (start_a={start_a}, mlen={mlen}, \
         start_b={start_b}) violates the revcomp relation at some interior offset. \
         The caller must pass a real match — use Rat::get_match, \
         forward_match_length, or GrowingPatch::get_all_matches to obtain one."
    );

    // Surviving-edge counts on each side.
    let surv_a = n - mlen;
    let surv_b = m - mlen;

    // Both sides fully consumed -- see function docstring. Panic
    // rather than silently returning `None`.
    if surv_a == 0 && surv_b == 0 {
        unreachable!(
            "glue_raw_angles: both sides fully consumed (mlen={mlen} == n == m) -- \
             this would require B to be A's mirror image glued along its entire \
             perimeter, which has no planar polygon realisation. start_a={start_a}, \
             start_b={start_b}."
        );
    }

    // Build the surviving boundary in two passes: A's surviving edges
    // first, then B's. The two angles at the seam (positions 0 and
    // surv_a, when both sides contribute) get overwritten with the
    // junction angles below.
    let mut angles = Vec::with_capacity(surv_a + surv_b);
    for i in 0..surv_a {
        angles.push(self_angles[(start_a + mlen + i) % n]);
    }
    for i in 0..surv_b {
        angles.push(other_angles[(start_b + i) % m]);
    }

    if surv_a > 0 && surv_b > 0 {
        // Normal glue: both sides contribute. Two distinct new junctions.
        //
        //   ccw (a_yx) at position 0:        between B-survivors and A-survivors
        //   cw  (a_xy) at position surv_a:   between A-survivors and B-survivors
        let x_first = self_angles[(start_a + mlen) % n]; // first surviving A angle
        let x_last = self_angles[start_a % n]; // = self_angles[(start_a+n) % n]; last cyclic A angle
        let y_first = other_angles[start_b % m]; // first surviving B angle
        let y_last = other_angles[(start_b + surv_b) % m]; // last cyclic B angle

        let a_yx = junction_angle::<T>(x_first, y_last);
        let a_xy = junction_angle::<T>(y_first, x_last);

        angles[0] = a_yx;
        angles[surv_a] = a_xy;

        Some(GlueResult {
            angles,
            cw_junction_idx: Some(surv_a),
            ccw_junction_idx: Some(0),
            a_yx: Some(a_yx),
            a_xy: Some(a_xy),
        })
    } else if surv_b == 0 {
        // B-keystone (hole-filling) glue. See function docstring for
        // the geometric meaning -- in short, B fills a tile-shaped
        // cavity in a weakly-simple patch A, healing the pinch point
        // into a strictly simple polygon. The CW and CCW junctions
        // collapse into a single merged seam vertex.
        let x_first = self_angles[(start_a + mlen) % n];
        let x_last = self_angles[start_a % n];
        let y_first = other_angles[start_b % m];
        let merged = normalize_angle::<T>(x_first + x_last + y_first - T::turn());
        if merged.abs() == T::hturn() {
            return None; // degenerate ±hturn junction
        }
        angles[0] = merged;

        Some(GlueResult {
            angles,
            cw_junction_idx: None,
            ccw_junction_idx: Some(0),
            a_yx: Some(merged),
            a_xy: Some(merged),
        })
    } else {
        // A-keystone -- mirror of the B-keystone above with A and B
        // swapped. No production caller exercises this; see function
        // docstring. Panic rather than emit a malformed result.
        unreachable!(
            "glue_raw_angles: A-keystone case (surv_a=0, surv_b={surv_b}) is not implemented \
             -- no caller is expected to exercise it. start_a={start_a}, mlen={mlen}, \
             start_b={start_b}, n={n}, m={m}"
        );
    }
}
