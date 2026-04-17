use crate::stringmatch::MatchIndex;

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

pub struct CyclicMatchIndex {
    originals: Vec<Vec<i8>>,
    n: usize,
    index: MatchIndex,
}

impl CyclicMatchIndex {
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

    pub fn num_tiles(&self) -> usize {
        self.n
    }

    pub fn tile_len(&self, i: usize) -> usize {
        self.originals[i].len()
    }

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
}
