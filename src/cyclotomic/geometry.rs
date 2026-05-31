//! Utility functions to use cyclotomic rings for 2D geometry
use std::cmp::Ordering;

use super::linalg::{dot_sign, is_between, is_ccw, wedge_sign};
use super::traits::{HasZZ4, IntersectUnitSegments, IsRing, OneImag};

/// Where a query point sits relative to a simple polygon boundary.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PointLocation {
    Inside,
    Outside,
    Boundary,
}

/// Classify where `q` sits relative to the polygon given by `polygon`.
///
/// `polygon` is a slice of `n` distinct vertices in CCW order with **no**
/// closing duplicate; edges are `(polygon[i], polygon[(i+1) % n])`.
///
/// Returns [`Boundary`] if `q` coincides with a polygon vertex or lies
/// strictly between the endpoints of any edge. Otherwise the winding number
/// (horizontal-ray test, all sign arithmetic via `im_sign` + `wedge_sign`)
/// decides [`Inside`] vs [`Outside`].
///
/// All geometry is exact; no f64 is used at any point.
///
/// # Panics
///
/// Panics if `polygon` has fewer than 3 vertices.
pub fn point_in_polygon<ZZ: IsRing>(q: &ZZ, polygon: &[ZZ]) -> PointLocation {
    let n = polygon.len();
    assert!(n >= 3, "point_in_polygon: need at least 3 vertices");

    for v in polygon {
        if q == v {
            return PointLocation::Boundary;
        }
    }

    for i in 0..n {
        let a = &polygon[i];
        let b = &polygon[(i + 1) % n];
        let edge = *b - *a;
        let to_q = *q - *a;
        if wedge_sign(&edge, &to_q) == 0 && is_between(q, (a, b)) {
            return PointLocation::Boundary;
        }
    }

    let mut wn: i32 = 0;
    for i in 0..n {
        let a = &polygon[i];
        let b = &polygon[(i + 1) % n];
        let da = (*a - *q).im_sign();
        let db = (*b - *q).im_sign();
        if da <= 0 && db > 0 {
            let edge = *b - *a;
            let to_q = *q - *a;
            if wedge_sign(&edge, &to_q) > 0 {
                wn += 1;
            }
        } else if da > 0 && db <= 0 {
            let edge = *b - *a;
            let to_q = *q - *a;
            if wedge_sign(&edge, &to_q) < 0 {
                wn -= 1;
            }
        }
    }

    if wn != 0 {
        PointLocation::Inside
    } else {
        PointLocation::Outside
    }
}

/// Lexicographic comparison of two ring elements by `(Re, Im)`.
///
/// Exact: compares `Re(a) - Re(b)` and `Im(a) - Im(b)` via `re_sign`/`im_sign`
/// on the difference, with no f64 intermediate.
pub fn cmp_xy<ZZ: IsRing>(a: &ZZ, b: &ZZ) -> Ordering {
    let d = *a - *b;
    match d.re_sign().cmp(&0) {
        Ordering::Equal => d.im_sign().cmp(&0),
        ord => ord,
    }
}

/// Strict convex hull via Andrew's monotone chain.
///
/// Returns indices into `points` (in CCW order, starting from the
/// lexicographically smallest point) of the vertices that lie on the
/// strict convex hull. Colinear intermediate points are **excluded**
/// (when `wedge_sign == 0` the middle point is popped).
///
/// Every comparison and orientation test stays in exact ring arithmetic;
/// no f64 is used at any point.
///
/// # Panics
///
/// Panics if `points` is empty.
pub fn convex_hull<ZZ: IsRing>(points: &[ZZ]) -> Vec<usize> {
    assert!(!points.is_empty(), "convex_hull: empty input");
    let n = points.len();
    if n <= 2 {
        return (0..n).collect();
    }

    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&i, &j| cmp_xy(&points[i], &points[j]));

    let mut hull: Vec<usize> = Vec::with_capacity(n);

    // Lower hull: scan left to right.
    for &pi in &idx {
        while hull.len() >= 2 {
            let a = hull[hull.len() - 2];
            let b = hull[hull.len() - 1];
            if wedge_sign(&(points[b] - points[a]), &(points[pi] - points[b])) >= 0 {
                hull.pop();
            } else {
                break;
            }
        }
        hull.push(pi);
    }

    // Upper hull: scan right to left.
    let lower_len = hull.len() + 1;
    for &pi in idx.iter().rev() {
        while hull.len() >= lower_len {
            let a = hull[hull.len() - 2];
            let b = hull[hull.len() - 1];
            if wedge_sign(&(points[b] - points[a]), &(points[pi] - points[b])) >= 0 {
                hull.pop();
            } else {
                break;
            }
        }
        hull.push(pi);
    }

    hull.pop();

    hull
}

/// Label each point as "on strict convex hull" or not.
///
/// Returns `Vec<bool>` of the same length as `points`. A `true` at index `i`
/// means `points[i]` lies on the strict convex hull, with colinear
/// intermediates excluded.
///
/// Thin wrapper around [`convex_hull`]; all geometry stays in exact ring
/// arithmetic.
pub fn hull_labels<ZZ: IsRing>(points: &[ZZ]) -> Vec<bool> {
    let hull = convex_hull(points);
    let mut labels = vec![false; points.len()];
    for idx in hull {
        labels[idx] = true;
    }
    labels
}

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
pub fn lines_parallel<ZZ: IsRing>((a, b): (&ZZ, &ZZ), (c, d): (&ZZ, &ZZ)) -> bool {
    wedge_sign(&(*b - *a), &(*d - *c)) == 0
}

/// Return whether the line through `(a, b)` is perpendicular to the line through `(c, d)`.
///
/// Uses the identity `dot(b - a, d - c) == 0`.
pub fn lines_perp<ZZ: IsRing>((a, b): (&ZZ, &ZZ), (c, d): (&ZZ, &ZZ)) -> bool {
    dot_sign(&(*b - *a), &(*d - *c)) == 0
}

/// Return whether line segments AB and CD intersect.
///
/// Sharing an endpoint (`a == c`, etc.) is NOT counted as intersection.
/// Proper interior crossings, colinear overlaps, and T-touches (one
/// segment's endpoint strictly interior to the other segment) all
/// return true.
///
/// The fast-path [`intersect_unit_segments`] handles the same predicate
/// when both segments are unit-length; this generic version is the
/// reference for arbitrary-length segments. Both must agree on every
/// shared input, which is what the trait tests in
/// `integral_basis::tests::intersect_unit_segments_matches_generic`
/// enforce.
pub fn intersect<ZZ: IsRing + PartialEq>(&(a, b): &(ZZ, ZZ), &(c, d): &(ZZ, ZZ)) -> bool {
    if a == c || a == d || b == c || b == d {
        return false;
    }
    let c_on_ab = is_colinear(&c, &a, &b);
    let d_on_ab = is_colinear(&d, &a, &b);
    if c_on_ab && d_on_ab {
        // Colinear: any interior point of one segment falling on the
        // other counts as overlap.
        return is_between(&c, (&a, &b)) || is_between(&d, (&a, &b));
    }
    // T-touch cases: one endpoint lies strictly between the other
    // segment's endpoints. `is_between` is strict (endpoints already
    // excluded above), so this catches exactly the "endpoint on
    // interior" case that the strict-CCW branch below misses (wedge==0
    // gets read as wedge<0 by `is_ccw`).
    if c_on_ab && is_between(&c, (&a, &b)) {
        return true;
    }
    if d_on_ab && is_between(&d, (&a, &b)) {
        return true;
    }
    if is_colinear(&a, &c, &d) && is_between(&a, (&c, &d)) {
        return true;
    }
    if is_colinear(&b, &c, &d) && is_between(&b, (&c, &d)) {
        return true;
    }
    // Proper interior crossing: c, d strictly on opposite sides of line
    // ab AND a, b strictly on opposite sides of line cd.
    is_ccw(&a, (&c, &d)) != is_ccw(&b, (&c, &d)) && is_ccw(&a, (&b, &c)) != is_ccw(&a, (&b, &d))
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
pub fn intersect_unit_segments<ZZ: IntersectUnitSegments>(s1: &(ZZ, ZZ), s2: &(ZZ, ZZ)) -> bool {
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

/// Return `p` modulo an axis-aligned rectangle of integer dimensions
/// `(width, height)` anchored at `anchor`, by repeated lattice shifts.
///
/// `anchor` is the lower-left corner; `width` shifts along the real axis,
/// `height` along the imaginary axis. Result lands in
/// `[anchor.re, anchor.re + width] x [anchor.im, anchor.im + height]`
/// (closed upper bound, matching the prior `mod_bound` semantics).
///
/// Restricted to `HasZZ4` rings because the lattice shifts use
/// `ZZ::one_i()` (imaginary unit must be in the ring).
pub fn point_mod_rect<ZZ: HasZZ4>(p: &ZZ, anchor: &ZZ, (width, height): (i64, i64)) -> ZZ {
    let w_step: ZZ = ZZ::from(width);
    let h_step: ZZ = ZZ::one_i() * ZZ::from(height);
    let upper: ZZ = *anchor + w_step + h_step;

    let mut ret = *p;

    // Real axis: bring ret.re into [anchor.re, anchor.re + width].
    let mut was_less = false;
    while (ret - *anchor).re_sign() < 0 {
        was_less = true;
        ret = ret + w_step;
    }
    if !was_less {
        while (upper - ret).re_sign() < 0 {
            ret = ret - w_step;
        }
    }

    // Imaginary axis: same shape, using h_step and im_sign.
    let mut was_less = false;
    while (ret - *anchor).im_sign() < 0 {
        was_less = true;
        ret = ret + h_step;
    }
    if !was_less {
        while (upper - ret).im_sign() < 0 {
            ret = ret - h_step;
        }
    }

    ret
}

/// Return whether `p` lies in the half-open integer-bounded box
/// `[x_min, x_max) x [y_min, y_max)`.
///
/// Works for every ring (no `HasZZ4` requirement). The point's cell in
/// the unit square grid is computed via the ring's `CellFloor` impl
/// (exact for blessed rings), then the bounding-box test is a plain
/// integer compare. Half-open semantics match `cell_floor`'s natural
/// `[cx, cx+1) x [cy, cy+1)` cell shape.
pub fn point_in_int_bbox<ZZ: IsRing>(
    p: &ZZ,
    (x_min, y_min): (i64, i64),
    (x_max, y_max): (i64, i64),
) -> bool {
    let (cx, cy) = p.cell_floor();
    x_min <= cx && cx < x_max && y_min <= cy && cy < y_max
}

#[cfg(test)]
mod tests {
    use num_traits::{One, Zero};

    use super::super::rings::ZZ12;
    use super::super::traits::{Ccw, OneImag, Units};
    use super::*;

    type ZZi = ZZ12;

    // ---- intersect: T-touch bug (vertex on segment interior) ----

    /// Two unit-length ZZ12 edges where one segment's endpoint lies
    /// exactly on the strict interior of the other segment's line.
    /// In Z[omega_12] (rank-4 lattice) this is a real configuration:
    /// (1 - omega + omega^2) lies on the line from -2 to (-2 + omega)
    /// at parameter t = sqrt(3) - 1.
    ///
    /// Both `intersect` and `intersect_unit_segments` should return
    /// `true` for this pair -- it is a genuine T-touch self-intersection.
    /// This is the bug behind the OEIS A316192 mismatch at ZZ12 n>=9.
    #[test]
    fn test_intersect_zz12_t_touch_endpoint_on_segment() {
        let omega = ZZ12::ccw();
        // V_1 = -omega, V_2 = -1 - omega + omega^2.
        let v1 = -omega;
        let v2 = -ZZ12::one() - omega + omega * omega;
        // V_5 = -2, V_6 = -2 + omega.
        let v5 = -ZZ12::from(2);
        let v6 = v5 + omega;

        // V_2 is on the line through V_5..V_6: wedge((V_6 - V_5), (V_2 - V_5)) == 0.
        assert_eq!(
            wedge_sign(&(v6 - v5), &(v2 - v5)),
            0,
            "V_2 should be colinear with V_5-V_6"
        );
        // ...and strictly between V_5 and V_6 (not at either endpoint).
        assert!(
            is_between(&v2, (&v5, &v6)),
            "V_2 should be strictly interior to V_5-V_6"
        );

        // The two segments share no endpoint.
        assert_ne!(v1, v5);
        assert_ne!(v1, v6);
        assert_ne!(v2, v5);
        assert_ne!(v2, v6);

        // So (V_1, V_2) and (V_5, V_6) genuinely touch (at V_2 = interior
        // point of V_5-V_6). Both checks should report intersection.
        assert!(
            intersect(&(v1, v2), &(v5, v6)),
            "generic intersect missed T-touch (V_2 lies on V_5-V_6 interior)"
        );
        assert!(
            intersect_unit_segments(&(v1, v2), &(v5, v6)),
            "intersect_unit_segments missed T-touch (V_2 lies on V_5-V_6 interior)"
        );
    }

    // ---- cmp_xy ----

    #[test]
    fn test_cmp_xy() {
        let a: ZZi = ZZi::zero();
        let b: ZZi = ZZi::one();
        let e: ZZi = ZZi::one_i();

        assert_eq!(cmp_xy(&a, &a), Ordering::Equal);
        assert_eq!(cmp_xy(&a, &b), Ordering::Less);
        assert_eq!(cmp_xy(&b, &a), Ordering::Greater);
        assert_eq!(cmp_xy(&a, &e), Ordering::Less);
        assert_eq!(cmp_xy(&e, &b), Ordering::Less);
    }

    // ---- convex_hull basic ----

    #[test]
    fn test_convex_hull_triangle() {
        let pts: Vec<ZZi> = vec![ZZi::zero(), ZZi::one(), ZZi::one_i()];
        let hull = convex_hull(&pts);
        assert_eq!(hull.len(), 3);
    }

    #[test]
    fn test_convex_hull_square() {
        let pts: Vec<ZZi> = vec![
            ZZi::zero(),
            ZZi::one(),
            ZZi::one() + ZZi::one_i(),
            ZZi::one_i(),
        ];
        let hull = convex_hull(&pts);
        assert_eq!(hull.len(), 4);
    }

    #[test]
    fn test_convex_hull_colinear_excluded() {
        // (0,0), (1,0), (2,0), (2,1), (2,2), (1,2), (0,2), (0,1)
        let pts: Vec<ZZi> = vec![
            ZZi::zero(),                                // 0: (0,0)
            ZZi::one(),                                 // 1: (1,0)
            ZZi::from(2),                               // 2: (2,0)
            ZZi::from(2) + ZZi::one_i(),                // 3: (2,1)
            ZZi::from(2) + ZZi::from(2) * ZZi::one_i(), // 4: (2,2)
            ZZi::one() + ZZi::from(2) * ZZi::one_i(),   // 5: (1,2)
            ZZi::from(2) * ZZi::one_i(),                // 6: (0,2)
            ZZi::one_i(),                               // 7: (0,1)
        ];
        let hull = convex_hull(&pts);
        assert_eq!(hull.len(), 4);
        let hull_set: std::collections::HashSet<usize> = hull.into_iter().collect();
        assert!(hull_set.contains(&0));
        assert!(hull_set.contains(&2));
        assert!(hull_set.contains(&4));
        assert!(hull_set.contains(&6));
        assert!(!hull_set.contains(&1));
        assert!(!hull_set.contains(&3));
        assert!(!hull_set.contains(&5));
        assert!(!hull_set.contains(&7));
    }

    // ---- convex_hull via Snake ----

    use crate::geom::snake::Snake;
    use crate::geom::tiles;

    fn snake_hull_labels<T: IsRing>(snake: &Snake<T>) -> Vec<bool> {
        assert!(snake.is_closed());
        let n = snake.len();
        hull_labels(&snake.representative()[..n])
    }

    fn validate_hull_labels<T: IsRing>(snake: &Snake<T>, labels: &[bool]) {
        let n = snake.len();
        assert_eq!(labels.len(), n);

        let hull_indices: Vec<usize> = (0..n).filter(|&i| labels[i]).collect();
        assert!(!hull_indices.is_empty(), "hull must be non-empty");

        let pts = &snake.representative()[..n];

        // Strict convexity: every triple of consecutive hull vertices is CCW.
        let h = &hull_indices;
        for k in 0..h.len() {
            let a = h[k];
            let b = h[(k + 1) % h.len()];
            let c = h[(k + 2) % h.len()];
            let w = wedge_sign(&(pts[b] - pts[a]), &(pts[c] - pts[b]));
            assert!(
                w > 0,
                "hull not strictly convex at hull vertices [{a},{b},{c}], wedge={w}"
            );
        }

        // Non-hull vertices are inside or on the hull boundary.
        for i in 0..n {
            if labels[i] {
                continue;
            }
            let v = pts[i];
            for k in 0..h.len() {
                let a = pts[h[k]];
                let b = pts[h[(k + 1) % h.len()]];
                let w = wedge_sign(&(b - a), &(v - a));
                assert!(
                    w >= 0,
                    "non-hull vertex {i} is outside hull edge [{}..{}], wedge={w}",
                    h[k],
                    h[(k + 1) % h.len()]
                );
            }
        }
    }

    #[test]
    fn test_hull_labels_hexagon_all_on_hull() {
        let hex: Snake<ZZ12> = tiles::hexagon();
        let labels = snake_hull_labels(&hex);
        assert_eq!(labels.len(), 6);
        assert!(
            labels.iter().all(|&b| b),
            "all hex vertices should be on hull"
        );
        validate_hull_labels(&hex, &labels);
    }

    #[test]
    fn test_hull_labels_triangle_all_on_hull() {
        let tri: Snake<ZZ12> = tiles::triangle();
        let labels = snake_hull_labels(&tri);
        assert_eq!(labels.len(), 3);
        assert!(labels.iter().all(|&b| b));
        validate_hull_labels(&tri, &labels);
    }

    #[test]
    fn test_hull_labels_square_all_on_hull() {
        let sq: Snake<ZZ12> = tiles::square();
        let labels = snake_hull_labels(&sq);
        assert_eq!(labels.len(), 4);
        assert!(labels.iter().all(|&b| b));
        validate_hull_labels(&sq, &labels);
    }

    #[test]
    fn test_hull_labels_rect_colinear_intermediates() {
        // 2x2 rect with extra vertices on each edge:
        // (0,0),(1,0),(2,0),(2,1),(2,2),(1,2),(0,2),(0,1)
        let rect: Snake<ZZ12> = Snake::try_from(&[0i8, 0, 3, 0, 3, 0, 3, 0]).unwrap();
        let labels = snake_hull_labels(&rect);
        assert_eq!(labels.len(), 8);
        // Corners on hull, intermediates not.
        assert!(labels[0], "(0,0) should be on hull");
        assert!(!labels[1], "(1,0) colinear intermediate");
        assert!(labels[2], "(2,0) should be on hull");
        assert!(!labels[3], "(2,1) colinear intermediate");
        assert!(labels[4], "(2,2) should be on hull");
        assert!(!labels[5], "(1,2) colinear intermediate");
        assert!(labels[6], "(0,2) should be on hull");
        assert!(!labels[7], "(0,1) colinear intermediate");
        validate_hull_labels(&rect, &labels);
    }

    #[test]
    fn test_hull_labels_l_shape() {
        // L-shape boundary, 8 unit segments, CCW:
        //   (0,0)→(1,0)→(2,0)→(2,1)→(1,1)→(1,2)→(0,2)→(0,1)→(0,0)
        // Vertex 4=(1,1) is the concave corner. (1,0) and (0,1) are colinear
        // intermediates on hull edges.
        let l: Snake<ZZ12> = Snake::try_from(&[0i8, 0, 3, 3, -3, 3, 3, 0]).unwrap();
        assert!(l.is_closed());
        let labels = snake_hull_labels(&l);
        assert_eq!(labels.len(), 8);
        assert!(labels[0], "(0,0) on hull");
        assert!(!labels[1], "(1,0) colinear on bottom edge");
        assert!(labels[2], "(2,0) on hull");
        assert!(labels[3], "(2,1) on hull");
        assert!(!labels[4], "(1,1) concave corner");
        assert!(labels[5], "(1,2) on hull");
        assert!(labels[6], "(0,2) on hull");
        assert!(!labels[7], "(0,1) colinear on left edge");
        validate_hull_labels(&l, &labels);
    }

    #[test]
    fn test_hull_labels_longer_rect_multiple_colinear() {
        // 3x2 rect: (0,0),(1,0),(2,0),(3,0),(3,1),(3,2),(2,2),(1,2),(0,2),(0,1)
        let rect: Snake<ZZ12> = Snake::try_from(&[0i8, 0, 0, 3, 0, 3, 0, 0, 3, 0]).unwrap();
        assert!(rect.is_closed());
        let labels = snake_hull_labels(&rect);
        assert_eq!(labels.len(), 10);
        // Hull corners: (0,0)=0, (3,0)=3, (3,2)=5, (0,2)=8
        assert!(labels[0]);
        assert!(!labels[1], "(1,0) colinear");
        assert!(!labels[2], "(2,0) colinear");
        assert!(labels[3], "(3,0) corner");
        assert!(!labels[4], "(3,1) colinear");
        assert!(labels[5], "(3,2) corner");
        assert!(!labels[6], "(2,2) colinear");
        assert!(!labels[7], "(1,2) colinear");
        assert!(labels[8], "(0,2) corner");
        assert!(!labels[9], "(0,1) colinear");
        validate_hull_labels(&rect, &labels);
    }

    #[test]
    fn test_hull_labels_spectre() {
        let spec: Snake<ZZ12> = tiles::spectre();
        let labels = snake_hull_labels(&spec);
        assert_eq!(labels.len(), 14);
        let hull_count = labels.iter().filter(|&&b| b).count();
        assert!(
            hull_count < 14,
            "spectre is non-convex, some vertices must be non-hull"
        );
        assert!(hull_count >= 4, "hull should have at least 4 vertices");
        validate_hull_labels(&spec, &labels);
    }

    #[test]
    fn test_hull_labels_cyclic_invariant() {
        use crate::geom::rat::Rat;

        let hex: Snake<ZZ12> = tiles::hexagon();
        let labels_base = snake_hull_labels(&hex);

        let rat: Rat<ZZ12> = Rat::from_unchecked(&hex);
        for shift in 0..6 {
            let shifted = rat.cycled(shift);
            let shifted_snake: Snake<ZZ12> = Snake::from_slice_unsafe(shifted.seq());
            let labels_shifted = snake_hull_labels(&shifted_snake);
            let hull_count_base = labels_base.iter().filter(|&&b| b).count();
            let hull_count_shifted = labels_shifted.iter().filter(|&&b| b).count();
            assert_eq!(
                hull_count_base, hull_count_shifted,
                "hull vertex count changed at shift {shift}"
            );
            validate_hull_labels(&shifted_snake, &labels_shifted);
        }
    }

    #[test]
    fn test_hull_labels_cross_ring() {
        use crate::cyclotomic::ZZ24;
        use crate::geom::tiles as tiles24;

        let hex12: Snake<ZZ12> = tiles::hexagon();
        let hex24: Snake<ZZ24> = tiles24::hexagon();

        let labels12 = snake_hull_labels(&hex12);
        let labels24 = snake_hull_labels(&hex24);

        assert_eq!(
            labels12, labels24,
            "same shape should give same hull labels"
        );
        validate_hull_labels(&hex12, &labels12);
        validate_hull_labels(&hex24, &labels24);
    }

    #[test]
    fn test_hull_labels_u_shape() {
        // 3x2 bounding box with 1x1 notch carved from top-center.
        //   (0,2) ___ (1,2) (2,2) (3,2)
        //        |   |           |
        //        |   +---+       |      (1,1),(2,1) = notch bottom
        //        +---------------+
        //      (0,0)           (3,0)
        //
        // 12 vertices, CCW from (0,0):
        //   0:(0,0) 1:(1,0) 2:(2,0) 3:(3,0)
        //   4:(3,1) 5:(3,2) 6:(2,2) 7:(2,1)
        //   8:(1,1) 9:(1,2) 10:(0,2) 11:(0,1)
        let u: Snake<ZZ12> = Snake::try_from(&[0i8, 0, 0, 3, 0, 3, 3, -3, -3, 3, 3, 0]).unwrap();
        assert!(u.is_closed());
        let labels = snake_hull_labels(&u);
        assert_eq!(labels.len(), 12);

        assert!(labels[0], "(0,0) BL corner");
        assert!(!labels[1], "(1,0) colinear bottom");
        assert!(!labels[2], "(2,0) colinear bottom");
        assert!(labels[3], "(3,0) BR corner");
        assert!(!labels[4], "(3,1) colinear right");
        assert!(labels[5], "(3,2) TR corner");
        assert!(!labels[6], "(2,2) colinear top");
        assert!(!labels[7], "(2,1) concave inner");
        assert!(!labels[8], "(1,1) concave inner");
        assert!(!labels[9], "(1,2) colinear top");
        assert!(labels[10], "(0,2) TL corner");
        assert!(!labels[11], "(0,1) colinear left");

        let hull_count = labels.iter().filter(|&&b| b).count();
        assert_eq!(hull_count, 4, "U-shape hull should have exactly 4 vertices");
        validate_hull_labels(&u, &labels);
    }

    #[test]
    fn test_hull_labels_u_shape_cyclic_invariant() {
        use crate::geom::rat::Rat;

        let u: Snake<ZZ12> = Snake::try_from(&[0i8, 0, 0, 3, 0, 3, 3, -3, -3, 3, 3, 0]).unwrap();
        let labels_base = snake_hull_labels(&u);
        let hull_count_base = labels_base.iter().filter(|&&b| b).count();
        assert_eq!(hull_count_base, 4);

        let rat: Rat<ZZ12> = Rat::from_unchecked(&u);
        for shift in 0..12 {
            let shifted = rat.cycled(shift);
            let shifted_snake: Snake<ZZ12> = Snake::from_slice_unsafe(shifted.seq());
            let labels_shifted = snake_hull_labels(&shifted_snake);
            let hull_count_shifted = labels_shifted.iter().filter(|&&b| b).count();
            assert_eq!(
                hull_count_base, hull_count_shifted,
                "hull vertex count changed at shift {shift}"
            );
            validate_hull_labels(&shifted_snake, &labels_shifted);
        }
    }

    // ---- point_in_polygon ----

    fn snake_points<T: IsRing>(snake: &Snake<T>) -> Vec<T> {
        let n = snake.len();
        snake.representative()[..n].to_vec()
    }

    #[test]
    fn test_pip_rect_interior() {
        // 2x2 rect: (0,0),(1,0),(2,0),(2,1),(2,2),(1,2),(0,2),(0,1)
        // (1,1) is strictly interior (not a vertex, not on any edge).
        let rect: Snake<ZZ12> = Snake::try_from(&[0i8, 0, 3, 0, 3, 0, 3, 0]).unwrap();
        let pts = snake_points(&rect);
        assert_eq!(
            point_in_polygon(&(ZZ12::one() + ZZ12::one_i()), &pts),
            PointLocation::Inside,
            "(1,1) should be Inside the 2x2 rect"
        );
    }

    #[test]
    fn test_pip_rect_exterior() {
        let rect: Snake<ZZ12> = Snake::try_from(&[0i8, 0, 3, 0, 3, 0, 3, 0]).unwrap();
        let pts = snake_points(&rect);
        let far: ZZ12 = ZZ12::from(5) + ZZ12::from(5) * ZZ12::one_i();
        assert_eq!(
            point_in_polygon(&far, &pts),
            PointLocation::Outside,
            "(5,5) should be Outside"
        );
        let below: ZZ12 = -ZZ12::one_i();
        assert_eq!(
            point_in_polygon(&below, &pts),
            PointLocation::Outside,
            "(0,-1) should be Outside"
        );
    }

    #[test]
    fn test_pip_square_boundary_vertex() {
        let sq: Snake<ZZ12> = tiles::square();
        let pts = snake_points(&sq);
        for i in 0..pts.len() {
            assert_eq!(
                point_in_polygon(&pts[i], &pts),
                PointLocation::Boundary,
                "square vertex {i} should be Boundary"
            );
        }
    }

    #[test]
    fn test_pip_hex_boundary_vertex() {
        let hex: Snake<ZZ12> = tiles::hexagon();
        let pts = snake_points(&hex);
        for i in 0..pts.len() {
            assert_eq!(
                point_in_polygon(&pts[i], &pts),
                PointLocation::Boundary,
                "hex vertex {i} should be Boundary"
            );
        }
    }

    #[test]
    fn test_pip_triangle_boundary_vertex() {
        let tri: Snake<ZZ12> = tiles::triangle();
        let pts = snake_points(&tri);
        for i in 0..pts.len() {
            assert_eq!(
                point_in_polygon(&pts[i], &pts),
                PointLocation::Boundary,
                "triangle vertex {i} should be Boundary"
            );
        }
    }

    #[test]
    fn test_pip_rect_colinear_vertex_boundary() {
        // 2x2 rect with colinear intermediates.
        // (1,0) is vertex 1 but also colinear on bottom edge — still Boundary.
        let rect: Snake<ZZ12> = Snake::try_from(&[0i8, 0, 3, 0, 3, 0, 3, 0]).unwrap();
        let pts = snake_points(&rect);
        assert_eq!(
            point_in_polygon(&ZZ12::one(), &pts),
            PointLocation::Boundary,
            "(1,0) is a polygon vertex, must be Boundary"
        );
    }

    #[test]
    fn test_pip_u_shape_notch_vertices_are_boundary() {
        // U-shape notch corners (1,1) and (2,1) are polygon vertices.
        let u: Snake<ZZ12> = Snake::try_from(&[0i8, 0, 0, 3, 0, 3, 3, -3, -3, 3, 3, 0]).unwrap();
        let pts = snake_points(&u);
        for i in 0..pts.len() {
            assert_eq!(
                point_in_polygon(&pts[i], &pts),
                PointLocation::Boundary,
                "U vertex {i} should be Boundary"
            );
        }
    }

    #[test]
    fn test_pip_u_shape_exterior() {
        let u: Snake<ZZ12> = Snake::try_from(&[0i8, 0, 0, 3, 0, 3, 3, -3, -3, 3, 3, 0]).unwrap();
        let pts = snake_points(&u);
        let far: ZZ12 = ZZ12::from(5) + ZZ12::from(5) * ZZ12::one_i();
        assert_eq!(
            point_in_polygon(&far, &pts),
            PointLocation::Outside,
            "(5,5) should be Outside the U"
        );
    }

    #[test]
    fn test_pip_l_shape_vertices_boundary() {
        let l: Snake<ZZ12> = Snake::try_from(&[0i8, 0, 3, 3, -3, 3, 3, 0]).unwrap();
        let pts = snake_points(&l);
        for i in 0..pts.len() {
            assert_eq!(
                point_in_polygon(&pts[i], &pts),
                PointLocation::Boundary,
                "L vertex {i} should be Boundary"
            );
        }
    }

    #[test]
    fn test_pip_cyclic_invariant() {
        use crate::geom::rat::Rat;

        let rect: Snake<ZZ12> = Snake::try_from(&[0i8, 0, 3, 0, 3, 0, 3, 0]).unwrap();
        let rat: Rat<ZZ12> = Rat::from_unchecked(&rect);

        for shift in 0..8 {
            let shifted = rat.cycled(shift);
            let shifted_snake: Snake<ZZ12> = Snake::from_slice_unsafe(shifted.seq());
            let shifted_pts = snake_points(&shifted_snake);
            for j in 0..shifted_pts.len() {
                assert_eq!(
                    point_in_polygon(&shifted_pts[j], &shifted_pts),
                    PointLocation::Boundary,
                    "shift {shift}, vertex {j}: must be Boundary"
                );
            }
        }
    }

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
