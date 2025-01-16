use crate::grid::UnitSquareGrid;
use crate::zz::{HasZZ12, HasZZ4, HasZZ6};
use crate::zzbase::ZZNum;
use crate::zzutil::intersect;
use crate::zzutil::{normalize_angle, upscale_angles};
use num_complex::Complex;
use num_traits::ToPrimitive;
use std::fmt::Debug;
use std::fmt::Display;

/// Representation of a turtle (i.e. an oriented point).
pub struct Turtle<T: ZZNum> {
    /// Position in the complex integer plane.
    pub pos: T,

    /// Facing direction (interpreted modulo full turn).
    pub dir: i8,
}

impl<T: ZZNum> Turtle<T> {
    /// Create a new turtle with given location and orientation.
    pub fn new(pos: T, dir: i8) -> Self {
        Self { pos, dir }
    }
}

impl<T: ZZNum> Default for Turtle<T> {
    /// Return a canonical turtle, i.e. located at the origin
    /// and looking in the direction of the positive real axis.
    fn default() -> Self {
        Turtle {
            pos: T::zero(),
            dir: 0,
        }
    }
}

/// Blueprint of a polyline or polygon that consists of
/// unit-length segments in a chosen complex integer ring.
///
/// It can be seen as (a subset of) the free monoid of angle sequences,
/// that is interpreted as an action on a point in the plane
/// and can be quotiented by equal total dislocation.
/// If we allow self intersecting sequences, we actually have
/// the full free monoid.
///
/// From this perspective, a closed snake corresponds to the quotient class that
/// contains the identity action / monoid element (i.e. no movement) as a representative.
///
/// Note that this structure does not support a useful inversion,
/// because we cannot include "turning around" at the end of the chain
/// into the formalism. The snake is in so far asymmetric, until it is closed.
/// A closed snake traces out a path that can be inverted (see the `Rat` type).
#[derive(Debug, Clone)]
pub struct Snake<T: ZZNum> {
    /// Abstract sequence of unit segments that point into different directions
    /// (i.e. turtle movement instructions).
    /// Each number corresponds to the corresponding rotated unit in the
    /// underlying cyclotomic ring.
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
    grid: UnitSquareGrid,

    /// If set to true, the snake will NOT prevent self-intersections.
    allow_intersections: bool,
}

impl<I: ToPrimitive, T: ZZNum> TryFrom<&[I]> for Snake<T> {
    type Error = &'static str;

    /// Create a snake from an angle sequence.
    /// Equivalent to adding all angles sequentially.
    fn try_from(angles: &[I]) -> Result<Self, Self::Error> {
        let mut result = Self::new();
        match result.extend_from_slice(angles) {
            Ok(()) => Ok(result),
            Err(_) => Err("Self-intersecting angle sequence!"),
        }
    }
}
impl<const N: usize, I: ToPrimitive, T: ZZNum> TryFrom<&[I; N]> for Snake<T> {
    type Error = &'static str;

    fn try_from(angles: &[I; N]) -> Result<Self, Self::Error> {
        Self::try_from(angles.as_slice())
    }
}

impl<T: ZZNum> Snake<T> {
    /// Return a new empty snake that is guaranteed to be free of self-intersections.
    pub fn new() -> Self {
        let mut grid = UnitSquareGrid::new();
        grid.add((0, 0), 0);
        Self {
            points: vec![T::zero(); 1],
            angles: Vec::new(),
            ang_sum: 0,
            grid: grid,
            allow_intersections: false,
        }
    }

    /// Return a new empty snake that allows self-intersecting segment chains.
    pub fn new_unchecked() -> Self {
        Self {
            allow_intersections: true,
            ..Self::new()
        }
    }

    /// Return a new empty snake that allows self-intersecting segment chains
    /// from a given sequence of angles.
    pub fn from_slice_unchecked<I: ToPrimitive>(angles: &[I]) -> Self {
        let mut result = Self::new_unchecked();
        result.extend_from_slice(angles).unwrap(); // Err will not happen with unchecked
        result
    }

    // ----

    /// Returns whether the snake is unchecked (i.e. created with `new_unchecked`).
    pub fn is_unchecked(&self) -> bool {
        self.allow_intersections
    }

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
            pos: self.offset(),
            dir: self.direction(),
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
            .add(UnitSquareGrid::cell_of(new_pt), self.points.len());
        // add point to representative polyline
        self.points.push(new_pt);
    }

    /// Return a polyline by tracing out snake with given turtle.
    pub fn to_polyline(&self, turtle: Turtle<T>) -> Vec<T> {
        let mut result = Self::new();
        *result.points.last_mut().unwrap() = turtle.pos;
        result.ang_sum = turtle.dir as i64;

        for angle in self.angles.iter() {
            result.add_unsafe(*angle);
        }

        result.points
    }

    pub fn to_polyline_f64(&self, turtle: Turtle<T>) -> Vec<(f64, f64)> {
        self.to_polyline(turtle)
            .iter()
            .map(|p| {
                let Complex { re: x, im: y } = p.complex();
                (x, y)
            })
            .collect()
    }

    /// Return all representative segments of the current snakes
    /// that intersect with at least one of the given cells.
    fn cell_segs(&self, cells: &[(i64, i64)]) -> Vec<(T, T)> {
        let mut seg_pt_indices: Vec<(usize, usize)> = Vec::new();
        for pt_idx in self.grid.get_cells(cells) {
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
        if self.allow_intersections {
            return true; // anything goes! WOOHOOO!!!
        }

        if self.is_closed() {
            return false; // already closed -> adding anything makes it non-simple
        }

        // end points of new candidate segment
        let new_seg @ (prev_pt, new_pt) = self.next_seg(angle);
        // unit square lattice cell neighborhood of new segment
        let neighborhood = UnitSquareGrid::seg_neighborhood_of(prev_pt, new_pt);
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

    /// Extend a snake from a sequence of angles in-place.
    pub fn extend_from_slice<I: ToPrimitive>(&mut self, angles: &[I]) -> Result<(), &str> {
        for angle in angles {
            if !self.add(angle.to_i8().unwrap()) {
                return Err("Cannot extend, self-intersecting sequence!");
            }
        }
        Ok(())
    }

    /// Concatenate given snake to this one in-place.
    pub fn extend(&mut self, other: &Snake<T>) -> Result<(), &str> {
        self.extend_from_slice(other.angles.as_slice())
    }

    /// Return concatenation of two snakes.
    pub fn concat(&self, other: &Snake<T>) -> Result<Self, &str> {
        let errmsg = "Cannot concatenate, self-intersecting sequence!";
        let mut result = Self::new();
        // FIXME: why can't I return the error message from the inner result?! (borrow issues)
        match result.extend(self) {
            Ok(()) => {}
            Err(_) => return Err(errmsg),
        }
        match result.extend(other) {
            Ok(()) => {}
            Err(_) => return Err(errmsg),
        }
        return Ok(result);
    }

    /// Return twice the area of the represented polygon, computed using the shoelace formula.
    /// If the snake is not closed, will implicitly add the segment from the last to first point.
    /// See: https://en.wikipedia.org/wiki/Shoelace_formula
    pub fn double_area(&self) -> T {
        let mut result = T::zero();
        for i in 1..self.len() {
            // println!("{} x {} = {}", self.points[i - 1].wedge(&self.points[i]);
            result = result + self.points[i - 1].wedge(&self.points[i]);
        }
        if !self.is_closed() {
            result = result + self.points[self.len() - 1].wedge(&self.points[0]);
        }

        return result;
    }

    /// Point-in-polygon check using raycasting.
    pub fn is_point_inside(&self) -> bool {
        if !self.is_closed() {
            return false;
        }

        return true;
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

    /// Return sequence of a square tile over a compatible ring (divisible by 4).
    pub fn square<T: ZZNum + HasZZ4>() -> Snake<T> {
        Snake::try_from(upscale_angles::<T>(4, &[1, 1, 1, 1]).as_slice()).unwrap()
    }

    /// Return sequence of a equilateral triangle tile over a compatible ring (divisible by 6).
    pub fn triangle<T: ZZNum + HasZZ6>() -> Snake<T> {
        Snake::try_from(upscale_angles::<T>(6, &[2, 2, 2]).as_slice()).unwrap()
    }

    /// Return sequence of a hexagon tile over a compatible ring (divisible by 6).
    pub fn hexagon<T: ZZNum + HasZZ6>() -> Snake<T> {
        Snake::try_from(upscale_angles::<T>(6, &[1, 1, 1, 1, 1, 1]).as_slice()).unwrap()
    }

    /// Return sequence of the spectre tile over a compatible ring (divisible by 12).
    pub fn spectre<T: ZZNum + HasZZ12>() -> Snake<T> {
        Snake::try_from(
            upscale_angles::<T>(12, &[3, 2, 0, 2, -3, 2, 3, 2, -3, 2, 3, -2, 3, -2]).as_slice(),
        )
        .unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::traits::Ccw;
    use crate::zz::{ZZ12, ZZ24};
    use crate::zzbase::ZZBase;
    use constants::{hexagon, spectre, square, triangle};
    use num_rational::Ratio;
    use num_traits::{One, Zero};

    #[test]
    fn test_spectre() {
        // test adaptive scaling dependent on chosen ring
        let seq_zz24: Vec<i8> = vec![6, 4, 0, 4, -6, 4, 6, 4, -6, 4, 6, -4, 6, -4];
        let s: Snake<ZZ24> = spectre();
        assert_eq!(s.angles(), seq_zz24);
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

        let s3 = s.concat(&s2).unwrap();
        s.extend(&s2).unwrap();
        assert_eq!(s, s3);

        // check that state is as expected
        assert_eq!(s.len(), 7);
        assert_eq!(s.angle_sum(), 15);
        assert_eq!(s.direction(), 3);
        assert_eq!(s.head().dir, s.direction());
        assert_eq!(s.head().pos, s.offset());

        // rep. polyline = instantiated with default turtle
        let l_exp: Vec<ZZ12> = s.representative().to_vec();
        let l = s.to_polyline(Turtle::default());
        assert_eq!(l, l_exp);

        // check instantiation with a non-trivial turtle
        let t = Turtle {
            pos: ZZ12::ccw().scale(7),
            dir: -2,
        };
        let l2_exp: Vec<ZZ12> = l_exp
            .iter()
            .map(|pt| (*pt * ZZ12::unit(t.dir)) + t.pos)
            .collect();
        let l2 = s.to_polyline(t);
        assert_eq!(l2, l2_exp);
    }

    #[test]
    fn test_closed_snake() {
        // test that once a snake is closed, the initial angle is fixed up
        // (the first angle is relative to the initially undefined predecessor)
        let sq1: Snake<ZZ12> = Snake::try_from(&[0, 3, 3, 3]).unwrap();
        assert!(sq1.is_closed());
        assert_eq!(sq1.angle_sum(), ZZ12::turn() as i64);
        let sq2: Snake<ZZ12> = Snake::try_from(&[0, -3, -3, -3]).unwrap();
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
        assert!(s.cell_segs(&[(0, 0)]).is_empty());

        s.add(0);
        assert_eq!(s.cell_segs(&[(0, 0)]), &[(ZZ12::zero(), ZZ12::one())]);

        s.add(0);
        assert_eq!(s.cell_segs(&[(0, 0)]), &[(ZZ12::zero(), ZZ12::one())]);
        assert_eq!(
            s.cell_segs(&[(1, 0)]),
            &[
                (ZZ12::zero(), ZZ12::one()),
                (ZZ12::one(), ZZ12::one().scale(2))
            ]
        );
        assert_eq!(
            s.cell_segs(&[(2, 0)]),
            &[(ZZ12::one(), ZZ12::one().scale(2))]
        );
    }

    #[test]
    fn test_can_add() {
        let mut s: Snake<ZZ12> = Snake::try_from(&[3, 3, 3]).unwrap();
        assert!(s.add(3)); // can add, closes shape
        assert!(!s.add(0)); // cannot add, shape is already closed

        let mut s2: Snake<ZZ12> = Snake::try_from(&[0, 3, 3, 3]).unwrap();
        assert!(!s2.add(3)); // cannot add, closes shape not in origin

        let mut s3: Snake<ZZ12> = Snake::try_from(&[0, 4]).unwrap();
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
        s4.extend(&Snake::try_from(&[1, 1, 3, 5]).unwrap()).unwrap();
        assert_eq!(s4.len(), 10);
        for i in (-ZZ12::hturn() + 1)..(ZZ12::hturn() - 1) {
            assert!(!s4.add(i));
        }
    }

    #[test]
    fn test_area() {
        use crate::gaussint::GaussInt;

        // area of an equilateral triangle with side length 1 is sqrt(3)/4
        let tri_area_sq3 = GaussInt::new(Ratio::<i64>::new_raw(1, 2), Ratio::<i64>::zero());
        assert_eq!(
            triangle::<ZZ12>().double_area(),
            ZZ12::new(&[0.into(), tri_area_sq3])
        );

        assert_eq!(
            square::<ZZ12>().double_area(),
            ZZ12::new(&[2.into(), 0.into()])
        );

        // area of a regular hexagon with side length 1 is 3*sqrt(3)/2
        assert_eq!(
            hexagon::<ZZ12>().double_area(),
            ZZ12::new(&[0.into(), 3.into()])
        );
    }

    #[test]
    fn test_eq_ord() {
        let s1: Snake<ZZ12> = Snake::try_from(&[1, 2, 3]).unwrap();
        let mut s2: Snake<ZZ12> = Snake::try_from(s1.angles.as_slice()).unwrap();
        assert_eq!(s1, s2);

        s2.add(4);
        assert!(s1 != s2);
        assert!(s1 < s2);

        let s3: Snake<ZZ12> = Snake::try_from(&[-1, 4, 3]).unwrap();
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
