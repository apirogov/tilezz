use super::zz::ZZNum;
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
#[derive(Debug, Clone, Hash)]
pub struct Snake<T: ZZNum> {
    /// Abstract sequence of unit segments that point into
    /// different directions (i.e. turtle movement instructions).
    angles: Vec<i8>,

    /// Sum of outer angles (i.e. sum of the angle sequence).
    ang_sum: i64,

    /// Representative polyline vertices
    /// (always non-empty and starting at origin).
    points: Vec<T>,
}

impl<T: ZZNum> Snake<T> {
    /// Return a new empty snake.
    pub fn new() -> Self {
        Self {
            points: vec![T::zero(); 1],
            angles: Vec::new(),
            ang_sum: 0,
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

    /// Returns true if the represented path is closed,
    /// i.e. the first and last point is the same (and equal to 0).
    pub fn is_closed(&self) -> bool {
        *self.points.last().unwrap() == T::zero()
    }

    // ----

    pub fn offset(&self) -> T {
        *self.points.last().unwrap()
    }

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

    /// Add a segment to the snake without checking for
    /// denormalization, degeneracy or self-intersection.
    fn add_unsafe(&mut self, angle: i8) {
        let old_direction = self.direction();

        self.angles.push(angle);
        self.ang_sum += angle as i64;

        let last = self.points.last().unwrap();
        self.points
            .push(*last + (T::unit(old_direction) * T::unit(angle)));
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

    /// Check that an angle does not represent a half-turn
    /// (it is invalid because it is self-intersecting).
    fn is_angle_valid(angle: i8) -> bool {
        (angle % T::turn()).abs() != T::hturn()
    }

    /// Normalize an angle to (-half-turn, half-turn). This is
    /// used to have a unique symbolic snake representation.
    fn normalize_angle(angle: i8) -> i8 {
        let a = angle % T::turn();
        if a.abs() <= T::hturn() {
            a
        } else {
            -(a.signum() * T::turn() - a)
        }
    }

    /// Check if the new segment can be added without
    /// causing a self-intersection.
    fn can_add(&self, _angle: i8) -> bool {
        // TODO: implement, ensure that we NEVER self-intersect
        true
    }

    /// Add a new segment to the snake.
    /// Returns true on success or false if the new segment
    /// would cause a self-intersection (point is rejected).
    pub fn add(&mut self, angle: i8) -> bool {
        // check that angle is not degenerate
        // for a simple polyline or polygon
        assert!(Self::is_angle_valid(angle));

        // normalize angle to be in (-halfturn, ..., halfturn)
        let a = Self::normalize_angle(angle);

        // check for induced self-intersections
        if !self.can_add(a) {
            return false;
        }

        // everything is ok -> add new segment
        self.add_unsafe(a);
        return true;
    }

    /// Concatenate two snakes in-place.
    pub fn extend(&mut self, other: &Snake<T>) {
        for angle in other.angles.iter() {
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

    /// Check that a snake is representing a simple polygon
    /// (and thus can be converted into a rat).
    pub fn is_rat(&self) -> bool {
        // NOTE: no need to check that the polygon is simple,
        // as we ensure this continuously (prevent self-intersections).
        !self.is_empty() && self.is_closed()
    }

    pub fn slice(&self, _: usize, _: usize) -> Self {
        Self::new()
        // TODO: implement slicing
        // (do we need a ISnake trait for Snake, SnakeView and Rat ?)
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

// TODO: implement Rat<T>
// (snake must be non-empty + closed + have ccw angle sum)

#[cfg(test)]

mod tests {
    use super::*;
    use crate::traits::Ccw;
    use crate::zz::ZZ12;
    use crate::zzbase::ZZBase;
    use num_traits::Zero;

    #[test]
    fn test_valid_angles() {
        assert!(Snake::<ZZ12>::is_angle_valid(0));
        assert!(Snake::<ZZ12>::is_angle_valid(1));
        assert!(Snake::<ZZ12>::is_angle_valid(7));
        assert!(Snake::<ZZ12>::is_angle_valid(12));
        assert!(!Snake::<ZZ12>::is_angle_valid(6));
        assert!(!Snake::<ZZ12>::is_angle_valid(-6));
        assert!(!Snake::<ZZ12>::is_angle_valid(18));
    }

    #[test]
    fn test_normalize_angles() {
        assert_eq!(Snake::<ZZ12>::normalize_angle(0), 0);
        assert_eq!(Snake::<ZZ12>::normalize_angle(1), 1);
        assert_eq!(Snake::<ZZ12>::normalize_angle(-1), -1);
        assert_eq!(Snake::<ZZ12>::normalize_angle(5), 5);
        assert_eq!(Snake::<ZZ12>::normalize_angle(-5), -5);
        assert_eq!(Snake::<ZZ12>::normalize_angle(6), 6);
        assert_eq!(Snake::<ZZ12>::normalize_angle(-6), -6);
        assert_eq!(Snake::<ZZ12>::normalize_angle(7), -5);
        assert_eq!(Snake::<ZZ12>::normalize_angle(-7), 5);
        assert_eq!(Snake::<ZZ12>::normalize_angle(11), -1);
        assert_eq!(Snake::<ZZ12>::normalize_angle(-11), 1);
        assert_eq!(Snake::<ZZ12>::normalize_angle(12), 0);
        assert_eq!(Snake::<ZZ12>::normalize_angle(-12), 0);
        assert_eq!(Snake::<ZZ12>::normalize_angle(13), 1);
        assert_eq!(Snake::<ZZ12>::normalize_angle(-13), -1);
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

    // TODO: test adding invalid angle (value, or self-intersecting)

    #[test]
    fn test_basic() {
        let mut s: Snake<ZZ12> = Snake::new();
        assert!(s.is_empty());
        assert!(s.is_closed());
        assert!(!s.is_rat());
        assert_eq!(s.len(), 0);
        assert_eq!(s.angle_sum(), 0);
        assert_eq!(s.direction(), 0);
        assert_eq!(s.offset(), ZZ12::zero());

        s.add(-9);
        assert!(!s.is_empty());
        assert!(!s.is_closed());
        assert!(!s.is_rat());
        assert_eq!(*s.angles().last().unwrap(), 3);
        assert_eq!(s.direction(), 3);
        assert_eq!(s.offset(), ZZ12::unit(3));

        s.add(4);
        s.add(5);

        let mut s2: Snake<ZZ12> = Snake::new();
        s2.add(0);
        s2.add(-2);
        s2.add(5);

        // check behavior of extend and concat
        let s3 = s.concat(&s2);
        s.extend(&s2);
        assert_eq!(s, s3);

        // check that state is as expected
        assert_eq!(s.len(), 6);
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

        let mut square: Snake<ZZ12> = Snake::new();
        for _ in 0..4 {
            square.add(3);
        }
        assert!(square.is_rat());
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

    #[test]
    fn test_comparisons() {}
}
