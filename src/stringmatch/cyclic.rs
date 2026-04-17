use crate::stringmatch::MatchIndex;

/// A maximal reverse-complementary match between two cyclic angle sequences.
///
/// Represents a shared boundary segment between two polygonal tiles:
/// tile A's subsequence starting at `pos_a` matches tile B's subsequence
/// traced in reverse (CW direction) starting at `pos_b`, for `len` edges.
///
/// The match is maximal: it cannot be extended in either direction cyclically.
pub struct CyclicMatch {
    pub tile_a: usize,
    pub pos_a: usize,
    pub tile_b: usize,
    pub pos_b: usize,
    pub len: usize,
}

impl std::fmt::Debug for CyclicMatch {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "CyclicMatch(tile_a={}, pos_a={}, tile_b={}, pos_b={}, len={})",
            self.tile_a, self.pos_a, self.tile_b, self.pos_b, self.len
        )
    }
}

/// Index for finding all maximal reverse-complementary matches across
/// a collection of cyclic angle sequences (polygonal tile boundaries).
///
/// Each input sequence represents the exterior angles of a closed polygon.
/// The index finds all maximal-length boundary segments shared between any
/// pair of tiles, where one tile's segment is the reverse complement of the other's
/// (representing the same geometric boundary traced in opposite directions).
///
/// # Construction
///
/// Internally, each sequence is doubled (for cyclic wraparound) and its reverse
/// complement is also indexed. A single `MatchIndex` is built over all 2n sequences.
///
/// # Usage
///
/// ```ignore
/// use tilezz::stringmatch::CyclicMatchIndex;
///
/// let tiles: Vec<Vec<i8>> = vec![
///     vec![3, 2, 0, 2, -3, 2, 3, 2, -3, 2, 3, -2, 3, -2],
///     vec![1, 2, 3, 4],
/// ];
/// let idx = CyclicMatchIndex::new(&tiles);
/// let matches = idx.maximal_rc_matches(0, 1);
/// for m in &matches {
///     println!("tile {} pos {} <-> tile {} pos {}, len {}",
///              m.tile_a, m.pos_a, m.tile_b, m.pos_b, m.len);
/// }
/// ```
pub struct CyclicMatchIndex {
    originals: Vec<Vec<i8>>,
    n: usize,
    index: MatchIndex,
}

impl CyclicMatchIndex {
    /// Build a cyclic match index over the given angle sequences.
    ///
    /// # Panics
    ///
    /// Panics if `strings` is empty.
    pub fn new(strings: &[Vec<i8>]) -> Self {
        assert!(!strings.is_empty(), "need at least one string");

        let n = strings.len();

        let doubled: Vec<Vec<i8>> = strings
            .iter()
            .map(|s| s.iter().copied().chain(s.iter().copied()).collect())
            .collect();

        let doubled_rc: Vec<Vec<i8>> = strings
            .iter()
            .map(|s| {
                let rc: Vec<i8> = s.iter().rev().map(|&a| -a).collect();
                rc.iter().copied().chain(rc.iter().copied()).collect()
            })
            .collect();

        let mut all: Vec<Vec<i8>> = Vec::with_capacity(2 * n);
        all.extend(doubled);
        all.extend(doubled_rc);

        CyclicMatchIndex {
            originals: strings.to_vec(),
            n,
            index: MatchIndex::new(&all),
        }
    }

    /// Number of tiles in the index.
    pub fn num_tiles(&self) -> usize {
        self.n
    }

    /// Length (number of angles) of tile `i`.
    pub fn tile_len(&self, i: usize) -> usize {
        self.originals[i].len()
    }

    /// Find all maximal reverse-complementary matches between tiles `i` and `j`.
    ///
    /// Each `CyclicMatch` describes a shared boundary segment where tile A's
    /// angles at positions `pos_a..pos_a+len` (cyclically) match the reverse
    /// complement of tile B's angles at positions `pos_b-len+1..=pos_b` (cyclically).
    pub fn maximal_rc_matches(&self, i: usize, j: usize) -> Vec<CyclicMatch> {
        let raw = self.index.all_positive_matches(i, self.n + j);

        let n_a = self.originals[i].len();
        let n_b = self.originals[j].len();
        let max_len = n_a.min(n_b);

        if max_len == 0 {
            return vec![];
        }

        let a = &self.originals[i];
        let b = &self.originals[j];

        let mut result: Vec<CyclicMatch> = raw
            .iter()
            .filter_map(|m| {
                let len = m.a.len.min(max_len);
                if len == 0 {
                    return None;
                }
                let pos_a = m.a.start % n_a;
                let pos_b = n_b - 1 - m.b.start % n_b;

                if !is_cyclic_left_maximal(a, b, pos_a, pos_b, len) {
                    return None;
                }

                Some(CyclicMatch {
                    tile_a: i,
                    pos_a,
                    tile_b: j,
                    pos_b,
                    len,
                })
            })
            .collect();

        dedup_cyclic(&mut result);
        result
    }
}

fn is_cyclic_left_maximal(a: &[i8], b: &[i8], pos_a: usize, pos_b: usize, len: usize) -> bool {
    if len >= a.len() || len >= b.len() {
        return true;
    }
    let prev_a = a[(pos_a + a.len() - 1) % a.len()];
    let prev_b = -b[(pos_b + 1) % b.len()];
    prev_a != prev_b
}

fn dedup_cyclic(matches: &mut Vec<CyclicMatch>) {
    matches.sort_by(|a, b| {
        a.pos_a
            .cmp(&b.pos_a)
            .then_with(|| a.pos_b.cmp(&b.pos_b))
            .then_with(|| b.len.cmp(&a.len))
    });
    matches.dedup_by(|a, b| a.pos_a == b.pos_a && a.pos_b == b.pos_b);
}

#[cfg(test)]
mod tests {
    use super::*;

    fn naive_cyclic_rc_matches(a: &[i8], b: &[i8]) -> Vec<(usize, usize, usize)> {
        let n_a = a.len();
        let n_b = b.len();
        if n_a == 0 || n_b == 0 {
            return vec![];
        }
        let max_len = n_a.min(n_b);
        let mut result = Vec::new();

        for ia in 0..n_a {
            for ib_rc in 0..n_b {
                let pos_b = n_b - 1 - ib_rc;
                let mut l = 0;
                while l < max_len {
                    let pa = (ia + l) % n_a;
                    let prcb = (ib_rc + l) % n_b;
                    let rc_val = -b[n_b - 1 - prcb];
                    if a[pa] != rc_val {
                        break;
                    }
                    l += 1;
                }
                if l == 0 {
                    continue;
                }
                if !is_cyclic_left_maximal(a, b, ia, pos_b, l) {
                    continue;
                }
                result.push((ia, pos_b, l));
            }
        }

        result.sort();
        result.dedup();
        result
    }

    fn check_rc(strings: &[Vec<i8>], i: usize, j: usize) {
        let idx = CyclicMatchIndex::new(strings);
        let matches = idx.maximal_rc_matches(i, j);
        let expected = naive_cyclic_rc_matches(&strings[i], &strings[j]);

        let mut got: Vec<(usize, usize, usize)> =
            matches.iter().map(|m| (m.pos_a, m.pos_b, m.len)).collect();
        got.sort();
        got.dedup();

        assert_eq!(
            got, expected,
            "RC mismatch: tile {i} vs tile {j}\n  got={got:?}\n  exp={expected:?}"
        );
    }

    #[test]
    fn identical_tiles_full_match() {
        let s = vec![1i8, 2, 3];
        check_rc(&[s.clone(), s.clone()], 0, 1);
    }

    #[test]
    fn revcomp_tile_full_match() {
        let a = vec![1i8, 2, 3];
        let b: Vec<i8> = a.iter().rev().map(|&x| -x).collect();
        check_rc(&[a, b], 0, 1);
    }

    #[test]
    fn no_match() {
        let a = vec![1i8, 2, 3];
        let b = vec![4i8, 5, 6];
        check_rc(&[a, b], 0, 1);
    }

    #[test]
    fn partial_rc_match() {
        let a = vec![1i8, 2, 3, 4];
        let b = vec![5i8, 6, -3, -4];
        check_rc(&[a, b], 0, 1);
    }

    #[test]
    fn cyclic_wraparound() {
        let a = vec![1i8, 2, 3];
        let b = vec![-1i8, 5, -3];
        check_rc(&[a, b], 0, 1);
    }

    #[test]
    fn self_match() {
        let a = vec![1i8, 2, 1];
        check_rc(&[a], 0, 0);
    }

    #[test]
    fn square_tiles() {
        let sq = vec![2i8, 2, 2, -6];
        check_rc(&[sq.clone(), sq.clone()], 0, 1);
    }

    #[test]
    fn spectre_like_angles() {
        let a: Vec<i8> = vec![3, 2, 0, 2, -3, 2, 3, 2, -3, 2, 3, -2, 3, -2];
        check_rc(&[a.clone(), a.clone()], 0, 1);
    }

    #[test]
    fn empty_tile() {
        check_rc(&[vec![], vec![1i8, 2]], 0, 1);
        check_rc(&[vec![1i8, 2], vec![]], 0, 1);
    }

    #[test]
    fn single_angle_tiles() {
        check_rc(&[vec![2i8], vec![-2i8]], 0, 1);
        check_rc(&[vec![2i8], vec![2i8]], 0, 1);
    }

    #[test]
    fn three_tiles_pairwise() {
        let tiles: Vec<Vec<i8>> = vec![vec![1i8, 2, 3], vec![2i8, 3, 4], vec![-3i8, -2, 5]];
        for i in 0..tiles.len() {
            for j in 0..tiles.len() {
                check_rc(&tiles, i, j);
            }
        }
    }

    #[test]
    fn repeated_angles() {
        let a = vec![1i8, 1, 1];
        check_rc(&[a.clone(), a.clone()], 0, 1);
    }

    #[test]
    fn negative_angles() {
        let a = vec![-1i8, -2, -3];
        let b = vec![3i8, 2, 1];
        check_rc(&[a, b], 0, 1);
    }

    #[test]
    fn vs_naive_random() {
        let mut seed: u64 = 54321;
        for _ in 0..40 {
            let tiles: Vec<Vec<i8>> = (0..4)
                .map(|_| {
                    (0..6)
                        .map(|_| {
                            seed = seed
                                .wrapping_mul(6364136223846793005)
                                .wrapping_add(1442695040888963407);
                            ((seed % 7) as i8) - 3
                        })
                        .collect()
                })
                .collect();
            for i in 0..tiles.len() {
                for j in 0..tiles.len() {
                    check_rc(&tiles, i, j);
                }
            }
        }
    }

    fn verify_rc_content(a: &[i8], b: &[i8], pos_a: usize, pos_b: usize, len: usize) {
        let n_a = a.len();
        let n_b = b.len();
        assert!(len > 0, "match length must be positive");
        assert!(len <= n_a && len <= n_b, "match length exceeds tile size");
        for i in 0..len {
            let pa = (pos_a + i) % n_a;
            let pb = (n_b + pos_b - i) % n_b;
            assert_eq!(
                a[pa], -b[pb],
                "RC content mismatch at offset {i}: a[{pa}]={} vs -b[{pb}]={}",
                a[pa], -b[pb],
            );
        }
    }

    fn verify_all_matches(idx: &CyclicMatchIndex, i: usize, j: usize) {
        let matches = idx.maximal_rc_matches(i, j);
        let a = &idx.originals[i];
        let b = &idx.originals[j];
        for m in &matches {
            assert_eq!(m.tile_a, i);
            assert_eq!(m.tile_b, j);
            assert!(m.pos_a < a.len());
            assert!(m.pos_b < b.len());
            assert!(m.len > 0);
            verify_rc_content(a, b, m.pos_a, m.pos_b, m.len);
        }
    }

    #[test]
    fn e2e_spectre_self_match() {
        let angles: Vec<i8> = vec![3, 2, 0, 2, -3, 2, 3, 2, -3, 2, 3, -2, 3, -2];
        let tiles = vec![angles.clone(), angles.clone()];
        let idx = CyclicMatchIndex::new(&tiles);

        check_rc(&tiles, 0, 1);
        verify_all_matches(&idx, 0, 1);

        let matches = idx.maximal_rc_matches(0, 1);
        assert!(!matches.is_empty(), "spectre should have self-matches");

        let max_len = matches.iter().map(|m| m.len).max().unwrap();
        assert_eq!(max_len, 3, "spectre max RC self-match should be 3");
    }

    #[test]
    fn e2e_spectre_pairwise() {
        let a: Vec<i8> = vec![3, 2, 0, 2, -3, 2, 3, 2, -3, 2, 3, -2, 3, -2];
        let b: Vec<i8> = vec![0, 2, -3, 2, 3, 2, -3, 2, 3, -2, 3, -2, 3, 2];
        let tiles = vec![a, b];
        let idx = CyclicMatchIndex::new(&tiles);

        for i in 0..2 {
            for j in 0..2 {
                check_rc(&tiles, i, j);
                verify_all_matches(&idx, i, j);
            }
        }
    }

    #[test]
    fn e2e_tetrominos() {
        let t_o: Vec<i8> = vec![0, 1, 0, 1, 0, 1, 0, 1];
        let t_i: Vec<i8> = vec![0, 0, 0, 1, 1, 0, 0, 0, 1, 1];
        let t_t: Vec<i8> = vec![-1, 1, 1, -1, 1, 1, 0, 0, 1, 1];
        let tiles = vec![t_o, t_i, t_t];

        let idx = CyclicMatchIndex::new(&tiles);
        for i in 0..tiles.len() {
            for j in 0..tiles.len() {
                check_rc(&tiles, i, j);
                verify_all_matches(&idx, i, j);
            }
        }
    }

    #[test]
    fn e2e_tetromino_self_matches() {
        let t_s: Vec<i8> = vec![-1, 1, 1, 0, 1, -1, 1, 1, 0, 1];
        let tiles = vec![t_s];
        let idx = CyclicMatchIndex::new(&tiles);

        let matches = idx.maximal_rc_matches(0, 0);
        verify_all_matches(&idx, 0, 0);
        assert!(!matches.is_empty(), "S-tetromino should have self-matches");
    }

    #[test]
    fn e2e_hexagons_no_rc_match() {
        let hex: Vec<i8> = vec![1, 1, 1, 1, 1, 1];
        let tiles = vec![hex.clone(), hex.clone()];
        let idx = CyclicMatchIndex::new(&tiles);

        let matches = idx.maximal_rc_matches(0, 1);
        assert!(
            matches.is_empty(),
            "identical hexagons have no RC matches (all positive vs all negative)"
        );
    }

    #[test]
    fn e2e_hexagon_vs_reversed_hexagon() {
        let hex: Vec<i8> = vec![1, 1, 1, 1, 1, 1];
        let hex_rev: Vec<i8> = vec![-1, -1, -1, -1, -1, -1];
        let tiles = vec![hex, hex_rev];
        let idx = CyclicMatchIndex::new(&tiles);

        let matches = idx.maximal_rc_matches(0, 1);
        verify_all_matches(&idx, 0, 1);

        let max_len = matches.iter().map(|m| m.len).max().unwrap_or(0);
        assert_eq!(
            max_len, 6,
            "hexagon and reversed hexagon should fully match"
        );
    }

    #[test]
    fn e2e_penrose_rhombs() {
        let wide: Vec<i8> = vec![
            2, 0, -1, 2, -1, 0, 0, 3, 0, 0, 1, -2, 1, 0, 0, 2, 0, 0, -1, 2, -1, 0, 0, 3, 0, 1, -2,
            1, 0,
        ];
        let narrow: Vec<i8> = vec![
            1, 0, -1, 2, -1, 0, 0, 4, 0, 0, -1, 2, -1, 0, 0, 1, 0, 1, -2, 1, 0, 0, 4, 0, 0, 1, -2,
            1, 0,
        ];
        let tiles = vec![wide, narrow];
        let idx = CyclicMatchIndex::new(&tiles);

        for i in 0..2 {
            for j in 0..2 {
                check_rc(&tiles, i, j);
                verify_all_matches(&idx, i, j);
            }
        }
    }

    #[test]
    fn e2e_mixed_shapes() {
        let square: Vec<i8> = vec![1, 1, 1, 1];
        let triangle: Vec<i8> = vec![2, 2, 2];
        let rect: Vec<i8> = vec![0, 0, 1, 1, 0, 0, 1, 1];
        let tiles = vec![square, triangle, rect];

        let idx = CyclicMatchIndex::new(&tiles);
        for i in 0..tiles.len() {
            for j in 0..tiles.len() {
                check_rc(&tiles, i, j);
                verify_all_matches(&idx, i, j);
            }
        }
    }

    #[test]
    fn e2e_many_tiles_exhaustive() {
        let mut seed: u64 = 99999;
        let mut next_angle = || {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((seed % 7) as i8) - 3
        };

        let tiles: Vec<Vec<i8>> = (0..8)
            .map(|_| (0..10).map(|_| next_angle()).collect())
            .collect();

        let idx = CyclicMatchIndex::new(&tiles);
        for i in 0..tiles.len() {
            for j in 0..tiles.len() {
                check_rc(&tiles, i, j);
                verify_all_matches(&idx, i, j);
            }
        }
    }
}
