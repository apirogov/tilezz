use crate::zzutil::intersect;

use super::gaussint::GaussInt;
use super::zzbase::ZZNum;
use super::zzutil::{cell_of, indices_from_cells, normalize_angle, seg_neighborhood_of};
use num_complex::Complex;
use num_traits::{ToPrimitive, Zero};
use std::collections::HashMap;
use std::fmt::Debug;
use std::fmt::Display;

/// Representation of a turtle (i.e. an oriented point).
pub struct Turtle<T: ZZNum> {
    /// Position in the complex integer plane.
    pub position: T,

    /// Facing direction (interpreted modulo full turn).
    pub direction: i8,
}

impl<T: ZZNum> Default for Turtle<T> {
    /// Return a canonical turtle, i.e. located at the origin
    /// and looking in the direction of the positive real axis.
    fn default() -> Self {
        Turtle {
            position: T::zero(),
            direction: 0,
        }
    }
}

/// Blueprint of a polyline or polygon that consists of
/// unit-length segments in a chosen complex integer ring.
#[derive(Debug, Clone)]
pub struct Snake<T: ZZNum> {
    /// Abstract sequence of unit segments that point into
    /// different directions (i.e. turtle movement instructions).
    angles: Vec<i8>,

    /// Sum of outer angles (i.e. sum of the angle sequence).
    ang_sum: i64,

    /// Representative polyline vertices
    /// (always non-empty and starting at origin).
    points: Vec<T>,

    /// Structure to collect indices of points that fall
    /// into the same unit square grid cell.
    /// Used to minimize number of segments to be compared
    /// when checking for self-intersections
    /// (including touching segment endpoints).
    ///
    /// Each point with index i is part of at most two segments,
    /// with indices (i-1, i) and (i, i+1), if 0 < i < num_points
    grid: HashMap<GaussInt<i64>, Vec<usize>>,
}

impl<I: ToPrimitive, T: ZZNum> From<&[I]> for Snake<T> {
    /// Create a snake from an angle sequence.
    /// Equivalent to adding all angles sequentially.
    fn from(angles: &[I]) -> Self {
        let mut result = Self::new();
        for angle in angles {
            if !result.add(angle.to_i8().unwrap()) {
                panic!("Self-intersecting angle sequence!")
            }
        }
        result
    }
}

impl<const N: usize, I: ToPrimitive, T: ZZNum> From<&[I; N]> for Snake<T> {
    fn from(angles: &[I; N]) -> Self {
        Self::from(angles.as_slice())
    }
}

impl<T: ZZNum> Snake<T> {
    /// Return a new empty snake.
    pub fn new() -> Self {
        Self {
            points: vec![T::zero(); 1],
            angles: Vec::new(),
            ang_sum: 0,
            grid: HashMap::from([(GaussInt::zero(), vec![0])]),
        }
    }

    // ----

    /// Return the number of unit segments of this snake.
    pub fn len(&self) -> usize {
        self.angles.len()
    }

    /// Returns whether the snake is empty,
    /// i.e. that it has no segments.
    pub fn is_empty(&self) -> bool {
        self.angles.is_empty()
    }

    /// Return angle sequence this snake is defined by.
    pub fn angles(&self) -> &[i8] {
        self.angles.as_slice()
    }

    /// Return the sum of angles. If the snake is closed,
    /// corresponds to the sum of outer angles.
    pub fn angle_sum(&self) -> i64 {
        self.ang_sum
    }

    /// Return canonical representative polyline vertex sequence.
    pub fn representative(&self) -> &[T] {
        self.points.as_slice()
    }

    /// Returns true if the represented path is closed, i.e.
    /// not empty + first and last point is the same (and equal to 0).
    ///
    /// As we always ensure that the sequence is not self-intersecting,
    /// this also equals to being a simple polygon.
    pub fn is_closed(&self) -> bool {
        !self.is_empty() && *self.points.last().unwrap() == T::zero()
    }

    // ----

    /// Returns the offset from the origin, i.e. effective
    /// direct Euclidean distance travelled by the snake.
    pub fn offset(&self) -> T {
        *self.points.last().unwrap()
    }

    /// Return the current facing direction of the snake.
    pub fn direction(&self) -> i8 {
        (self.ang_sum % (T::turn() as i64)) as i8
    }

    /// The head of the snake is the final turtle state
    /// after tracing out all the steps in the given snake.
    /// The tail of a snake is always fixed at the origin.
    pub fn head(&self) -> Turtle<T> {
        Turtle {
            position: self.offset(),
            direction: self.direction(),
        }
    }

    // ----

    /// Return the two end points of the next segment,
    /// obtained from adding a new directed step.
    fn next_seg(&self, angle: i8) -> (T, T) {
        let old_direction = self.direction();
        let last = *self.points.last().unwrap();
        let new_pt = last + (T::unit(old_direction) * T::unit(angle));
        (last, new_pt)
    }

    /// Add a segment to the snake without checking for
    /// denormalization, degeneracy or self-intersection.
    fn add_unsafe(&mut self, angle: i8) {
        // compute next representative point (uses current orientation!)
        let (_, new_pt) = self.next_seg(angle);

        // append segment to symbolic angle sequence
        self.angles.push(angle);
        self.ang_sum += angle as i64;

        // register point in the grid
        self.grid
            .entry(cell_of(new_pt))
            .or_default()
            .push(self.points.len());

        // add point to representative polyline
        self.points.push(new_pt);
    }

    /// Return a polyline by tracing out snake with given turtle.
    pub fn to_polyline(&self, turtle: &Turtle<T>) -> Vec<T> {
        let mut result = Self::new();
        *result.points.last_mut().unwrap() = turtle.position;
        result.ang_sum = turtle.direction as i64;

        for angle in self.angles.iter() {
            result.add_unsafe(*angle);
        }

        result.points
    }

    pub fn to_polyline_f64(&self, turtle: &Turtle<T>) -> Vec<(f64, f64)> {
        self.to_polyline(&turtle)
            .iter()
            .map(|p| {
                let Complex { re: x, im: y } = p.complex();
                (x, y)
            })
            .collect()
    }

    /// Return all representative segments of the current snakes
    /// that intersect with at least one of the given cells.
    fn cell_segs(&self, cells: &[GaussInt<i64>]) -> Vec<(T, T)> {
        let mut seg_pt_indices: Vec<(usize, usize)> = Vec::new();
        for pt_idx in indices_from_cells(&self.grid, cells) {
            // each point is part of at least one and at most 2 segments
            if pt_idx != 0 {
                seg_pt_indices.push((pt_idx - 1, pt_idx));
            }
            // NOTE: self.len() = self.points.len()-1 = last pt idx
            if pt_idx != self.len() {
                seg_pt_indices.push((pt_idx, pt_idx + 1));
            }
        }
        seg_pt_indices.sort();
        seg_pt_indices.dedup();
        seg_pt_indices
            .iter()
            .map(|(ix1, ix2)| (self.points[*ix1], self.points[*ix2]))
            .collect()
    }

    /// Check if the new segment can be added without
    /// causing a self-intersection.
    ///
    /// Note that the correctness of the optimized check strongly relies
    /// on the assumption that all segments have unit length.
    fn can_add(&self, angle: i8) -> bool {
        if self.is_closed() {
            return false; // already closed -> adding anything makes it non-simple
        }

        // end points of new candidate segment
        let new_seg @ (prev_pt, new_pt) = self.next_seg(angle);
        // unit square lattice cell neighborhood of new segment
        let neighborhood = seg_neighborhood_of(prev_pt, new_pt);
        // all candidates for segment intersection
        let neighbor_segs = self.cell_segs(&neighborhood);

        let new_pt_nz = !new_pt.is_zero();
        for s @ (x, y) in neighbor_segs {
            // if the new point is NOT the origin (i.e. closes the snake),
            // explicitly check endpoints for intersection.
            // (the segment intersection ignores this edge case)
            if new_pt_nz && (new_pt == x || new_pt == y) {
                return false;
            }
            if intersect(&new_seg, &s) {
                return false;
            }
        }
        return true;
    }

    /// Add a new segment to the snake.
    /// Returns true on success or false if the new segment
    /// would cause a self-intersection (point is rejected).
    pub fn add(&mut self, angle: i8) -> bool {
        // normalize angle to be in (-halfturn, ..., halfturn)
        let a = normalize_angle::<T>(angle);
        // check that angle is not degenerate
        // for a simple polyline or polygon
        // (a half-turn step means walking the
        // last segment backwards, i.e. self-intersects)
        assert!(a.abs() != T::hturn());

        // check for induced self-intersections
        if !self.can_add(a) {
            return false;
        }

        // everything is ok -> add new segment
        self.add_unsafe(a);

        // if the snake is closed, we need to fix up the start angle
        if self.is_closed() {
            // the missing angle must complete one full turn
            // clockwise or counter-clockwise (simple polygon property)
            let target = (T::turn() as i64) * self.ang_sum.signum();
            let missing = target - (self.ang_sum - self.angles[0] as i64);
            self.angles[0] = missing as i8;
            self.ang_sum = target;
        }
        return true;
    }

    /// Concatenate given snake to this one in-place.
    pub fn extend(&mut self, other: &Snake<T>) {
        for angle in other.angles.iter() {
            // NOTE: cannot use unsafe_add here,
            // because there is no guarantee that
            // concatenating two snakes is valid.
            self.add(*angle);
        }
    }

    /// Return concatenation of two snakes.
    pub fn concat(&self, other: &Snake<T>) -> Self {
        let mut result = Self::new();
        result.extend(self);
        result.extend(other);
        result
    }
}

// For comparisons (Eq, Ord) and presentation (Display)
// we only care about the angles, all other data is derivative.
impl<T: ZZNum> PartialEq for Snake<T> {
    fn eq(&self, other: &Self) -> bool {
        self.angles == other.angles
    }
}
impl<T: ZZNum> Eq for Snake<T> {}
impl<T: ZZNum> PartialOrd for Snake<T> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.angles.partial_cmp(&other.angles)
    }
}
impl<T: ZZNum> Ord for Snake<T> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.angles.cmp(&other.angles)
    }
}
impl<T: ZZNum> Display for Snake<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.angles.fmt(f)
    }
}

pub mod constants {
    use super::*;

    /// Return sequence of the spectre tile over a compatible ring (divisible by 12).
    pub fn spectre<T: ZZNum>() -> Snake<T> {
        let ring = T::zz_params().full_turn_steps;
        assert_eq!(ring % 12, 0);

        let seq_zz12: Vec<i8> = vec![3, 2, 0, 2, -3, 2, 3, 2, -3, 2, 3, -2, 3, -2];
        let scale = ring / 12;
        let seq: Vec<i8> = seq_zz12.iter().map(|x| x * scale).collect();

        Snake::from(seq.as_slice())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::traits::Ccw;
    use crate::zz::{ZZ12, ZZ24};
    use crate::zzbase::ZZBase;
    use constants::spectre;
    use num_traits::{One, Zero};

    #[test]
    fn test_spectre() {
        // test adaptive scaling dependent on chosen ring
        let seq_zz24: Vec<i8> = vec![6, 4, 0, 4, -6, 4, 6, 4, -6, 4, 6, -4, 6, -4];
        let s: Snake<ZZ24> = spectre();
        assert_eq!(s.angles(), seq_zz24);
    }

    #[test]
    fn test_cell_of() {
        assert_eq!(cell_of(ZZ12::zero()), GaussInt::new(0, 0));
        assert_eq!(cell_of(ZZ12::unit(0)), GaussInt::new(1, 0));
        assert_eq!(cell_of(ZZ12::unit(1)), GaussInt::new(1, 1));
        assert_eq!(cell_of(ZZ12::unit(2)), GaussInt::new(1, 1));
        assert_eq!(cell_of(ZZ12::unit(3)), GaussInt::new(0, 1));
        assert_eq!(cell_of(ZZ12::unit(4)), GaussInt::new(-1, 1));
        assert_eq!(cell_of(ZZ12::unit(5)), GaussInt::new(-1, 1));
        assert_eq!(cell_of(ZZ12::unit(6)), GaussInt::new(-1, 0));
        assert_eq!(cell_of(ZZ12::unit(7)), GaussInt::new(-1, -1));
        assert_eq!(cell_of(ZZ12::unit(8)), GaussInt::new(-1, -1));
        assert_eq!(cell_of(ZZ12::unit(9)), GaussInt::new(0, -1));
        assert_eq!(cell_of(ZZ12::unit(10)), GaussInt::new(1, -1));
        assert_eq!(cell_of(ZZ12::unit(11)), GaussInt::new(1, -1));
    }

    #[test]
    fn test_seg_neighborhood_of() {
        let gint_vec = |pts: &[(i64, i64)]| {
            pts.iter()
                .map(|(x, y)| GaussInt::new(*x, *y))
                .collect::<Vec<GaussInt<i64>>>()
        };

        // segment fully contained in one cell
        let tmp1 = (ZZ12::one().scale(2) - ZZ12::unit(2) - ZZ12::unit(-1).scale(2)) * ZZ12::unit(1);
        let tmp2 = (ZZ12::one().scale(3) - ZZ12::unit(2) - ZZ12::unit(-1).scale(3)) * ZZ12::unit(2);
        let p1 = (tmp1 - tmp2) * ZZ12::unit(-2);
        let p2 = p1 + ZZ12::unit(2);
        let n0 = seg_neighborhood_of(p1, p2);
        let n0_exp = gint_vec(&[(-1, 0), (0, -1), (0, 0), (0, 1), (1, 0)]);
        assert_eq!(n0, n0_exp);

        // segment touches 2 cells (horizontal)
        let n1 = seg_neighborhood_of(ZZ12::zero(), ZZ12::unit(0));
        let n1_exp = gint_vec(&[
            (-1, 0),
            (0, -1),
            (0, 0),
            (0, 1),
            (1, -1),
            (1, 0),
            (1, 1),
            (2, 0),
        ]);
        assert_eq!(n1, n1_exp);

        // segment touches 2 cells (vertical)
        let n2 = seg_neighborhood_of(ZZ12::zero(), ZZ12::unit(3));
        let n2_exp = gint_vec(&[
            (-1, 0),
            (-1, 1),
            (0, -1),
            (0, 0),
            (0, 1),
            (0, 2),
            (1, 0),
            (1, 1),
        ]);
        assert_eq!(n2, n2_exp);

        // segment touches 2 cells (diagonal)
        let n3 = seg_neighborhood_of(ZZ12::zero(), ZZ12::unit(1));
        let n3_exp = gint_vec(&[
            (-1, 0),
            (0, -1),
            (0, 0),
            (0, 1),
            (1, 0),
            (1, 1),
            (1, 2),
            (2, 1),
        ]);
        assert_eq!(n3, n3_exp);
    }

    #[test]
    #[should_panic]
    fn test_add_invalid_angle() {
        let mut s: Snake<ZZ12> = Snake::new();
        s.add(6);
    }

    #[test]
    #[should_panic]
    fn test_add_invalid_angle_neg() {
        let mut s: Snake<ZZ12> = Snake::new();
        s.add(-6);
    }

    #[test]
    #[should_panic]
    fn test_add_invalid_angle_mod() {
        let mut s: Snake<ZZ12> = Snake::new();
        s.add(18);
    }

    #[test]
    fn test_basic() {
        let mut s: Snake<ZZ12> = Snake::new();
        assert!(s.is_empty());
        assert!(!s.is_closed());
        assert_eq!(s.len(), 0);
        assert_eq!(s.angle_sum(), 0);
        assert_eq!(s.direction(), 0);
        assert_eq!(s.offset(), ZZ12::zero());

        s.add(-9); // same as 3
        assert!(!s.is_empty());
        assert!(!s.is_closed());
        assert_eq!(*s.angles().last().unwrap(), 3);
        assert_eq!(s.direction(), 3);
        assert_eq!(s.offset(), ZZ12::unit(3));

        s.add(2);
        s.add(2);
        s.add(5);
        assert_eq!(s.len(), 4);

        // check behavior of extend and concat
        let mut s2: Snake<ZZ12> = Snake::new();
        s2.add(-2);
        s2.add(0);
        s2.add(5);
        assert_eq!(s2.len(), 3);

        let s3 = s.concat(&s2);
        s.extend(&s2);
        assert_eq!(s, s3);

        // check that state is as expected
        assert_eq!(s.len(), 7);
        assert_eq!(s.angle_sum(), 15);
        assert_eq!(s.direction(), 3);
        assert_eq!(s.head().direction, s.direction());
        assert_eq!(s.head().position, s.offset());

        // rep. polyline = instantiated with default turtle
        let l_exp: Vec<ZZ12> = s.representative().to_vec();
        let l = s.to_polyline(&Turtle::default());
        assert_eq!(l, l_exp);

        // check instantiation with a non-trivial turtle
        let t = Turtle {
            position: ZZ12::ccw().scale(7),
            direction: -2,
        };
        let l2_exp: Vec<ZZ12> = l_exp
            .iter()
            .map(|pt| (*pt * ZZ12::unit(t.direction)) + t.position)
            .collect();
        let l2 = s.to_polyline(&t);
        assert_eq!(l2, l2_exp);
    }

    #[test]
    fn test_closed_snake() {
        // test that once a snake is closed, the initial angle is fixed up
        // (the first angle is relative to the initially undefined predecessor)
        let sq1: Snake<ZZ12> = Snake::from(&[0, 3, 3, 3]);
        assert!(sq1.is_closed());
        assert_eq!(sq1.angle_sum(), ZZ12::turn() as i64);
        let sq2: Snake<ZZ12> = Snake::from(&[0, -3, -3, -3]);
        assert!(sq2.is_closed());
        assert_eq!(sq2.angle_sum(), -ZZ12::turn() as i64);

        let mut square: Snake<ZZ12> = Snake::new();
        assert!(!square.is_closed()); // empty

        for _ in 0..3 {
            square.add(3);
        }
        assert!(!square.is_closed()); // not closed

        square.add(3);
        assert!(square.is_closed()); // valid
    }

    #[test]
    fn test_cell_segs() {
        let mut s: Snake<ZZ12> = Snake::new();
        assert!(s.cell_segs(&[GaussInt::new(0, 0)]).is_empty());

        s.add(0);
        assert_eq!(
            s.cell_segs(&[GaussInt::new(0, 0)]),
            &[(ZZ12::zero(), ZZ12::one())]
        );

        s.add(0);
        assert_eq!(
            s.cell_segs(&[GaussInt::new(0, 0)]),
            &[(ZZ12::zero(), ZZ12::one())]
        );
        assert_eq!(
            s.cell_segs(&[GaussInt::new(1, 0)]),
            &[
                (ZZ12::zero(), ZZ12::one()),
                (ZZ12::one(), ZZ12::one().scale(2))
            ]
        );
        assert_eq!(
            s.cell_segs(&[GaussInt::new(2, 0)]),
            &[(ZZ12::one(), ZZ12::one().scale(2))]
        );
    }

    #[test]
    fn test_can_add() {
        let mut s: Snake<ZZ12> = Snake::from(&[3, 3, 3]);
        assert!(s.add(3)); // can add, closes shape
        assert!(!s.add(0)); // cannot add, shape is already closed

        let mut s2: Snake<ZZ12> = Snake::from(&[0, 3, 3, 3]);
        assert!(!s2.add(3)); // cannot add, closes shape not in origin

        let mut s3: Snake<ZZ12> = Snake::from(&[0, 4]);
        assert!(!s3.add(5)); // cannot add, crosses segment

        // try a more complicated example
        let mut s4: Snake<ZZ12> = Snake::new();
        assert!(s4.add(0));
        assert!(s4.add(5));

        // take sharpest turns possible
        assert!(!s4.add(4));
        assert!(s4.add(3));

        assert!(!s4.add(5));
        assert!(s4.add(4));

        assert!(!s4.add(2));
        assert!(s4.add(1));

        assert!(!s4.add(5));
        assert!(s4.add(4));

        // create a dead-end where we cannot continue the snake
        s4.extend(&Snake::from(&[1, 1, 3, 5]));
        assert_eq!(s4.len(), 10);
        for i in (-ZZ12::hturn() + 1)..(ZZ12::hturn() - 1) {
            assert!(!s4.add(i));
        }
    }

    #[test]
    fn test_eq_ord() {
        let s1: Snake<ZZ12> = Snake::from(&[1, 2, 3]);
        let mut s2: Snake<ZZ12> = Snake::from(s1.angles.as_slice());
        assert_eq!(s1, s2);

        s2.add(4);
        assert!(s1 != s2);
        assert!(s1 < s2);

        let s3: Snake<ZZ12> = Snake::from(&[-1, 4, 3]);
        assert!(s3 < s1);
        assert!(s3 < s2);

        let s4: Snake<ZZ12> = Snake::new();
        assert!(s4 < s1);
        assert!(s4 < s2);
        assert!(s4 < s3);
    }

    #[test]
    fn test_display() {
        let mut s: Snake<ZZ12> = Snake::new();
        for i in -2..2 {
            s.add(i);
        }
        s.add(7);
        assert_eq!(format!("{s}"), "[-2, -1, 0, 1, -5]");
    }
}
