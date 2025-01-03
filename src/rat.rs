use crate::snake::Snake;
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
}

impl<T: ZZNum> Rat<T> {
    /// Returns the number of line segments (= number of vertices).
    pub fn len(&self) -> usize {
        // NOTE: we store the sequence twice!
        self.angles.len() / 2
    }

    /// Returns the chirality (1 = ccw, -1 = cw).
    pub fn chirality(&self) -> i8 {
        self.angle_sum.signum() as i8
    }

    pub fn new(snake: &Snake<T>) -> Self {
        assert!(snake.is_closed());

        // normalize (lex. minimal rotation of cyclic sequence)
        let mut canonical = snake.angles().to_vec();
        let offset = lex_min_rot(canonical.as_slice());
        canonical.rotate_left(offset);

        // add a second copy (convenient for slicing and algorithms)
        canonical.extend(&canonical.clone());

        Self {
            phantom: PhantomData,
            angles: canonical,
            angle_sum: snake.angle_sum(),
        }
    }

    /// Return a slice of the angle sequence with length n, starting
    /// at given offset index with respect to the canonical sequence.
    pub fn slice(&self, start: usize, len: usize) -> &[i8] {
        assert!(len <= self.len());
        &self.angles[start..(start + len)]
    }

    /// Return the canonical outer angle sequence.
    pub fn seq(&self) -> &[i8] {
        self.slice(0, self.len())
    }

    /// Return chirally reflected description of the same shape
    /// (its sequence is equal to the reverse complement)
    pub fn revcomp(&self) -> Self {
        println!("{:?} {:?}", &self.seq(), &revcomp(self.seq()));
        // TODO: unsafe from slice, bypassing Snake (we know the revcomp is also a valid polygon)
        Self::new(&Snake::from(revcomp(self.seq()).as_slice()))
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
        self.angles.as_slice()[0..self.len()].fmt(f)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::snake::constants::spectre;
    use crate::zz::ZZ12;

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

        // should be sequence of the minimally rotated string
        let exp = "[-3, 2, 3, -2, 3, -2, 3, 2, 0, 2, -3, 2, 3, 2]";
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
        assert_eq!(r.len(), canonical.len());
        assert_eq!(r.seq(), canonical.as_slice());

        let c: Rat<ZZ12> = r.revcomp();
        assert!(r != c);

        // check chirality flipping due to revcomp
        assert_eq!(r.chirality(), 1);
        assert_eq!(c.chirality(), -1);

        // twice revcomp = original
        let cc: Rat<ZZ12> = c.revcomp();
        assert_eq!(r, cc);
    }
}
