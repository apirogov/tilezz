use crate::snake::{Snake, Turtle};
use crate::zzbase::ZZNum;
use crate::zzutil::revcomp;
use std::fmt::{Debug, Display};
use std::marker::PhantomData;

/// Booth's lexicographically minimal string rotation algorithm.
/// Returns index of the beginning of the lex. minimal rotation of given sequence in O(n).
///
/// See https://en.wikipedia.org/wiki/Lexicographically_minimal_string_rotation
///
fn lex_min_rot<T: Eq + Ord>(s: &[T]) -> usize {
    let n = s.len() as isize;
    let mut f: Vec<isize> = vec![-1; 2 * s.len()];
    let mut k: isize = 0;
    let comp_idx = |k, i| ((k + i + 1) % n) as usize;
    for j in 1..(2 * n) {
        let j_mod_n = (j % n) as usize;
        let mut i = f[(j - k - 1) as usize];
        while i != -1 && s[j_mod_n] != s[comp_idx(k, i)] {
            if s[j_mod_n] < s[comp_idx(k, i)] {
                k = j - i - 1;
            }
            i = f[i as usize];
        }
        if i == -1 && s[j_mod_n] != s[comp_idx(k, i)] {
            if s[j_mod_n] < s[comp_idx(k, i)] {
                k = j;
            }
            f[(j - k) as usize] = -1;
        } else {
            f[(j - k) as usize] = i + 1;
        }
    }
    return k as usize;
}

/// Given an angle sequence, computes its lexicographically
/// minimal rotation and returns two copies of that,
/// with an offset of the beginning of the sequence wrt.
/// the original input sequence.
fn prepare_seq(angles: &[i8]) -> (Vec<i8>, usize) {
    // normalize (lex. minimal rotation of cyclic sequence)
    let mut canonical = angles.to_vec();
    let offset = lex_min_rot(canonical.as_slice());
    canonical.rotate_left(offset);

    // add a second copy (convenient for slicing and algorithms)
    canonical.extend(&canonical.clone());

    (canonical, offset)
}

/// Given two slices of sequences, return length of longest match
/// (assuming that we start to match both at index 0)
/// and possibly a offset to be subtracted from the left side
/// (if the match can be extended to the left, assuming the sequences are cyclic)
fn match_length(x: &[i8], y: &[i8]) -> (usize, usize) {
    let min_len = x.len().min(y.len());
    if min_len < 2 {
        // at least one sequence has
        // no nodes or only single node -> no segment
        return (0, 0);
    }

    // we start checking at the 2nd point and look for the right side end
    let mut len = 1;
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
    let mut offset = 0;
    for i in 1..(min_len - len) {
        if x[x.len() - i] != y[y.len() - i] {
            // first non-complementary angle -> left end of match
            offset = i;
            len += i;
            break;
        }
    }
    return (len, offset);
}

#[derive(Debug, Clone)]
pub struct Rat<T: ZZNum> {
    /// Even though we don't need to use the underlying complex integer ring,
    /// we parametrize by that type for improved type safety (to not mix incompatible polygons).
    phantom: PhantomData<T>,

    /// Canonical lexicographically minimal angle sequence (used for Eq/Ord),
    /// twice (for easier use of cyclic structure with string algorithms)
    angles: Vec<i8>,

    /// Angle sum (can be used to get chirality and also infer the complex integer ring)
    angle_sum: i64,

    /// Indexing offset (used e.g. for slicing).
    cyc: usize,
}

impl<T: ZZNum> TryFrom<&Snake<T>> for Rat<T> {
    type Error = &'static str;
    /// Construct a new rat from a snake.
    fn try_from(snake: &Snake<T>) -> Result<Self, Self::Error> {
        if snake.is_closed() {
            Ok(Self::from_unchecked(snake))
        } else {
            Err("Given snake is not closed!")
        }
    }
}

impl<T: ZZNum> Rat<T> {
    /// Create a rat from an angle sequence.
    /// Assumes that the sequence describes a valid simple polygon.
    pub fn from_slice_unchecked(angles: &[i8]) -> Self {
        let angle_sum: i64 = angles.iter().map(|x| *x as i64).sum();
        assert_eq!(angle_sum.abs(), T::turn() as i64);
        let (seq, offset) = prepare_seq(angles);

        Self {
            phantom: PhantomData,
            angles: seq,
            angle_sum: angle_sum,
            cyc: (angles.len() - offset) % angles.len(),
        }
    }

    /// Create a rat from a snake.
    /// Assumes that the sequence describes a valid simple polygon.
    pub fn from_unchecked(snake: &Snake<T>) -> Self {
        Self::from_slice_unchecked(snake.angles())
    }

    fn revcomp_seq(&self, offset: i64) -> Vec<i8> {
        // NOTE: needed off-by-one offset to ensure that
        // the first angle refers to the same node as before
        revcomp(self.slice_from(offset + 1, self.len()))
    }

    /// Return chirally reflected description of the same shape
    /// (its sequence is equal to the reverse complement)
    pub fn reflect(&self) -> Self {
        Self::from_slice_unchecked(self.revcomp_seq(0).as_slice())
    }

    // ----

    /// Returns the number of line segments (= number of vertices).
    pub fn len(&self) -> usize {
        // NOTE: we store the sequence twice!
        self.angles.len() / 2
    }

    /// Returns the chirality (1 = ccw, -1 = cw).
    pub fn chirality(&self) -> i8 {
        self.angle_sum.signum() as i8
    }

    // ----

    /// Set the start node to be the beginning of the canonical sequence.
    pub fn canonical(mut self) -> Self {
        self.cyc = 0;
        self
    }

    /// Shift the starting node of the tile (does not affect equivalence).
    pub fn cycle(mut self, offset: i64) -> Self {
        self.cyc = ((self.cyc as i64 + offset) % (self.len() as i64)) as usize;
        self
    }

    /// Return a copy with shifted starting node.
    pub fn cycled(&self, offset: i64) -> Self {
        self.clone().cycle(offset)
    }

    // ----

    /// Return a slice of the angle sequence with length n, starting
    /// at given offset index with respect to the canonical sequence.
    fn slice_from_canonical(&self, start: usize, len: usize) -> &[i8] {
        assert!(len <= self.len());
        &self.angles[start..(start + len)]
    }

    /// Return a slice of the angle sequence with length n.
    pub fn slice_from(&self, start: i64, len: usize) -> &[i8] {
        self.slice_from_canonical(
            ((self.cyc as i64 + start) % self.len() as i64) as usize,
            len,
        )
    }

    /// Return a slice of the angle sequence with length n.
    pub fn slice(&self, len: usize) -> &[i8] {
        self.slice_from(0, len)
    }

    /// Return the outer angle sequence.
    pub fn seq(&self) -> &[i8] {
        self.slice(self.len())
    }

    // ----

    /// Return the maximal length of a match described by the given pair of indices.
    /// Note that this still does NOT mean that the tiles can be legally glued along the match,
    /// because outside of the match the sequence could have self-intersections.
    pub fn match_length(&self, (self_start, other_end): (i64, i64), other: &Self) -> usize {
        let x_seq = self.slice_from(self_start, self.len());
        let y_seq = other.revcomp_seq(other_end);
        match_length(x_seq, y_seq.as_slice()).0
    }

    /// Return the maximal match that touches the given pair of indices.
    /// The returned indices are normalized to be non-negative.
    /// Note that this still does NOT mean that the tiles can be legally glued along the match.
    pub fn get_match(
        &self,
        (self_start, other_end): (i64, i64),
        other: &Self,
    ) -> (i64, usize, i64) {
        let x_seq = self.slice_from(self_start, self.len());
        let y_seq = other.revcomp_seq(other_end);
        let (len, offset) = match_length(x_seq, y_seq.as_slice());

        let ext_start: i64 = (self_start - offset as i64).rem_euclid(self.len() as i64);
        let ext_end: i64 = (other_end + offset as i64).rem_euclid(self.len() as i64);
        (ext_start, len, ext_end)
    }

    /// Return a new rat resulting from glueing two rats along a common boundary
    /// that is given by a pair of two glued indices (the match is extended in both directions).
    pub fn try_glue(&self, matched_indices: (i64, i64), other: &Self) -> Result<Self, &str> {
        let (norm_start, mlen, norm_end) = self.get_match(matched_indices, other);

        // get boundaries without the match (but keep the endpoints)
        let x = self.slice_from(norm_start + mlen as i64, self.len() - mlen + 1);
        let y = self.slice_from(norm_end, other.len() - mlen + 1);

        // concatenate the sequences (without duplicating the endpoints)
        let mut glued_seq: Vec<i8> = x[..x.len() - 1].to_vec();
        glued_seq.extend_from_slice(&y[..y.len() - 1]);

        // fix glued vertex angles at the transition of the boundaries
        let a_yx = x[0] + y[y.len() - 1] - T::hturn();
        let a_xy = y[0] + x[x.len() - 1] - T::hturn();
        glued_seq[0] = a_yx;
        glued_seq[x.len() - 1] = a_xy;

        Snake::try_from(glued_seq.as_slice()).and_then(|s| Ok(Self::from_unchecked(&s)))
    }

    pub fn glue(&self, matched_indices: (i64, i64), other: &Self) -> Self {
        self.try_glue(matched_indices, other).unwrap()
    }

    // ----

    /// Return sequence of coordinates of the line segment chain describing a
    /// realized polygonal tile boundary. The first and last point will be equal.
    pub fn to_polyline_f64(&self, turtle: Turtle<T>) -> Vec<(f64, f64)> {
        Snake::from_slice_unchecked(self.seq()).to_polyline_f64(&turtle)
    }
}

impl<T: ZZNum> PartialEq for Rat<T> {
    fn eq(&self, other: &Self) -> bool {
        self.angles == other.angles
    }
}
impl<T: ZZNum> Eq for Rat<T> {}
impl<T: ZZNum> PartialOrd for Rat<T> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.angles.partial_cmp(&other.angles)
    }
}
impl<T: ZZNum> Ord for Rat<T> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.angles.cmp(&other.angles)
    }
}
impl<T: ZZNum> Display for Rat<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.seq().fmt(f)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::snake::constants::{hexagon, spectre};
    use crate::zz::ZZ12;

    #[test]
    fn test_lex_min_rot() {
        assert_eq!(lex_min_rot(&['h', 'e', 'l', 'l', 'o']), 1);
        assert_eq!(lex_min_rot(&[1, 3, 2, 1, 2, 1, 2, 0]), 7);
        assert_eq!(lex_min_rot(&[1, 3, 2, 1, -2, 1, 2, 0]), 4);
        assert_eq!(lex_min_rot(&[1, 1, 2, 1, 1, 1, 2, 1, 1, 2, 1, 1]), 10);
    }

    #[test]
    fn test_match_length() {
        // degenerate cases
        assert_eq!(match_length(&[1, 2, 3], &[]), (0, 0));
        assert_eq!(match_length(&[1, 2, 3], &[1]), (0, 0));

        // match cannot be longer than shortest sequence
        assert_eq!(match_length(&[1, 2, 3], &[1, 2]), (1, 0));
        // the end points do not matter
        assert_eq!(match_length(&[1, 2, 3], &[4, 5]), (1, 0));
        // the inside points must be equal
        assert_eq!(match_length(&[1, 2, 3, 4], &[5, 2, -3, 6]), (2, 0));
        assert_eq!(match_length(&[1, 2, 3, 4], &[5, 2, 3, 6]), (3, 0));
        // match can be extended to the left (input is starting in the middle)
        assert_eq!(match_length(&[3, 4, 1, 2], &[3, 6, 5, 2]), (3, 2));
    }

    #[test]
    fn test_len_display() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&spectre());
        // should show with start offset as used for initialization
        // (even though its normalized internally)
        let exp = "[3, 2, 0, 2, -3, 2, 3, 2, -3, 2, 3, -2, 3, -2]";
        assert_eq!(format!("{r}"), exp);
    }

    #[test]
    fn test_new_len_chirality_revcomp() {
        let spectre: Snake<ZZ12> = spectre();
        let spectre_angles = spectre.angles();

        let mut canonical = spectre_angles.to_vec();
        let offset = lex_min_rot(canonical.as_slice());
        canonical.rotate_left(offset);

        let r: Rat<ZZ12> = Rat::from_unchecked(&spectre);
        assert_eq!(
            r,
            Rat::from_unchecked(&Snake::from_slice_unchecked(canonical.as_slice()))
        );

        // we store two copies of the angle sequence, but it is reported for one
        assert_eq!(r.len(), spectre_angles.len());
        assert_eq!(r.seq(), spectre_angles);

        let c: Rat<ZZ12> = r.reflect();
        assert!(r != c);

        // check chirality flipping due to revcomp
        assert_eq!(r.chirality(), 1);
        assert_eq!(c.chirality(), -1);

        // twice revcomp = original
        let cc: Rat<ZZ12> = c.reflect();
        assert_eq!(r, cc);
    }

    #[test]
    fn test_reflect() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&spectre());
        assert!(r != r.reflect());

        // twice reflect -> restores original sequence + offset
        assert_eq!(r.reflect().reflect().seq(), r.seq());

        // reflecting the sequence keeps the current node constant
        // (needs off-by-one offset in implementation)
        let exp: &[i8] = &[-3, 2, -3, 2, -3, -2, 3, -2, -3, -2, 3, -2, 0, -2];
        assert_eq!(r.reflect().seq(), exp);
    }

    #[test]
    fn test_slice() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&spectre()).canonical();

        assert_eq!(r.slice(0).len(), 0);
        assert_eq!(r.slice(1).len(), 1);
        assert_eq!(r.slice(r.len()), r.seq());
        assert_eq!(r.slice_from(2, 3), &r.seq()[2..5]);
    }

    #[test]
    fn test_cycle() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&spectre()).canonical();

        assert!(r.seq() != r.cycled(1).seq());
        assert_eq!(r.seq(), r.cycled(1).cycled(-3).cycled(2).seq());
        assert_eq!(r.cycled(7).slice_from(-2, 3), r.slice_from(5, 3));

        // offset cycling does not affect equivalence, it is just convenience
        assert_eq!(r, r.cycled(2));

        // revcomp works respecting the current start offset
        assert_eq!(r.cycled(2), r.cycled(2).reflect().reflect());
    }

    #[test]
    fn test_rat_match_length() {
        // test manually verified match combinations
        let r: Rat<ZZ12> = Rat::from_unchecked(&spectre());
        assert_eq!(r.match_length((0, 0), &r), 2);
        assert_eq!(r.match_length((2, 0), &r), 4);
        assert_eq!(r.match_length((2, 6), &r), 1);
    }

    #[test]
    fn test_glue() {
        let h: Rat<ZZ12> = Rat::from_unchecked(&hexagon());
        let t1 = h.glue((1, 5), &h);
        assert_eq!(t1.seq(), &[-2, 2, 2, 2, 2, -2, 2, 2, 2, 2]);

        let r: Rat<ZZ12> = Rat::from_unchecked(&spectre());
        let t2 = r.glue((2, 0), &r);
        let mystic = &[
            0, 2, -3, 2, 3, -2, 3, -2, 3, 2, -3, 2, 0, 2, -3, 2, 3, 2, -3, 2,
        ];
        assert_eq!(t2.seq(), mystic);
    }
}
