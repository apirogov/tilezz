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

impl<T: ZZNum> Rat<T> {
    /// Create a rat from an angle sequence.
    /// Assumes that the sequence describes a valid simple polygon.
    fn from_seq_unsafe(angles: &[i8]) -> Self {
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

    /// Construct a new rat from a snake.
    pub fn new(snake: &Snake<T>) -> Self {
        assert!(snake.is_closed());
        Self::from_seq_unsafe(snake.angles())
    }

    /// Return chirally reflected description of the same shape
    /// (its sequence is equal to the reverse complement)
    pub fn reflect(&self) -> Self {
        Self::from_seq_unsafe(revcomp(self.seq()).as_slice())
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
    pub fn slice_from_canonical(&self, start: usize, len: usize) -> &[i8] {
        assert!(len <= self.len());
        &self.angles[start..(start + len)]
    }

    /// Return a slice of the angle sequence with length n.
    pub fn slice(&self, len: usize) -> &[i8] {
        self.slice_from_canonical(self.cyc, len)
    }

    /// Return the outer angle sequence.
    pub fn seq(&self) -> &[i8] {
        self.slice(self.len())
    }

    // ----

    /// Return sequence of coordinates of the line segment chain describing a
    /// realized polygonal tile boundary. The first and last point will be equal.
    pub fn to_polyline_f64(&self, turtle: Turtle<T>) -> Vec<(f64, f64)> {
        Snake::from(self.seq()).to_polyline_f64(&turtle)
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
    use crate::snake::constants::spectre;
    use crate::zz::{ZZ12, ZZ24};

    #[test]
    fn test_lex_min_rot() {
        assert_eq!(lex_min_rot(&['h', 'e', 'l', 'l', 'o']), 1);
        assert_eq!(lex_min_rot(&[1, 3, 2, 1, 2, 1, 2, 0]), 7);
        assert_eq!(lex_min_rot(&[1, 3, 2, 1, -2, 1, 2, 0]), 4);
        assert_eq!(lex_min_rot(&[1, 1, 2, 1, 1, 1, 2, 1, 1, 2, 1, 1]), 10);
    }

    #[test]
    fn test_len_display() {
        let r: Rat<ZZ12> = Rat::new(&spectre());
        // should show with start offset as used for initialization
        // (even though its normalized internally)
        let exp = "[3, 2, 0, 2, -3, 2, 3, 2, -3, 2, 3, -2, 3, -2]";
        assert_eq!(format!("{r}"), exp);
    }

    #[test]
    fn test_new_len_chirality_revcomp() {
        let s: Snake<ZZ12> = spectre();
        let spectre = s.angles();

        let mut canonical = spectre.to_vec();
        let offset = lex_min_rot(canonical.as_slice());
        canonical.rotate_left(offset);

        let r: Rat<ZZ12> = Rat::new(&Snake::from(spectre));
        assert_eq!(r, Rat::new(&Snake::from(canonical.as_slice())));

        // we store two copies of the angle sequence, but it is reported for one
        assert_eq!(r.len(), spectre.len());
        assert_eq!(r.seq(), spectre);

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
    fn test_cycle() {
        let r: Rat<ZZ24> = Rat::new(&spectre()).canonical();

        assert!(r.seq() != r.cycled(1).seq());
        assert_eq!(r.seq(), r.cycled(1).cycled(-3).cycled(2).seq());
        assert_eq!(r.slice_from_canonical(5, 3), r.cycled(5).slice(3));

        // offset cycling does not affect equivalence, it is just convenience
        assert_eq!(r, r.cycled(2));

        // revcomp works respecting the current start offset
        assert_eq!(r.cycled(2), r.cycled(2).reflect().reflect());
    }
}
