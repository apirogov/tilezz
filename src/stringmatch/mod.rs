mod cyclic;
mod lcp;
mod maximal;
mod rmq;
mod sais;

use std::collections::HashMap;

pub struct Interval {
    pub string: usize,
    pub start: usize,
    pub len: usize,
}

pub struct Match {
    pub a: Interval,
    pub b: Interval,
}

impl Match {
    fn new(string_i: usize, string_j: usize, start_i: usize, start_j: usize, len: usize) -> Self {
        if string_i == string_j && start_i > start_j {
            Match {
                a: Interval {
                    string: string_j,
                    start: start_j,
                    len,
                },
                b: Interval {
                    string: string_i,
                    start: start_i,
                    len,
                },
            }
        } else {
            Match {
                a: Interval {
                    string: string_i,
                    start: start_i,
                    len,
                },
                b: Interval {
                    string: string_j,
                    start: start_j,
                    len,
                },
            }
        }
    }
}

impl PartialEq for Match {
    fn eq(&self, other: &Self) -> bool {
        self.a.string == other.a.string
            && self.a.start == other.a.start
            && self.a.len == other.a.len
            && self.b.string == other.b.string
            && self.b.start == other.b.start
            && self.b.len == other.b.len
    }
}

impl Eq for Match {}

impl PartialOrd for Match {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Match {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.a
            .string
            .cmp(&other.a.string)
            .then_with(|| self.a.start.cmp(&other.a.start))
            .then_with(|| self.b.string.cmp(&other.b.string))
            .then_with(|| self.b.start.cmp(&other.b.start))
            .then_with(|| self.a.len.cmp(&other.a.len))
    }
}

pub struct MatchIndex {
    concat: Vec<u32>,
    sa: Vec<usize>,
    lcp: Vec<usize>,
    rank: Vec<usize>,
    doc: Vec<usize>,
    doc_pos: Vec<usize>,
    boundaries: Vec<usize>,
    string_lens: Vec<usize>,
    rmq: rmq::SparseTable,
    num_strings: usize,
}

impl MatchIndex {
    pub fn new(strings: &[Vec<i8>]) -> Self {
        assert!(!strings.is_empty(), "need at least one string");

        let num_strings = strings.len();
        let (concat, boundaries, string_lens) = Self::concatenate(strings);

        let sa = sais::build_suffix_array(&concat);
        let lcp = lcp::build_lcp_array(&concat, &sa);
        let rank = Self::build_rank(&sa);
        let (doc, doc_pos) = Self::build_doc_arrays(concat.len(), &boundaries, &string_lens);
        let rmq = rmq::SparseTable::new(&lcp);

        MatchIndex {
            concat,
            sa,
            lcp,
            rank,
            doc,
            doc_pos,
            boundaries,
            string_lens,
            rmq,
            num_strings,
        }
    }

    pub fn num_strings(&self) -> usize {
        self.num_strings
    }

    pub fn string_len(&self, i: usize) -> usize {
        self.string_lens[i]
    }

    pub fn maximal_matches(&self, i: usize, j: usize) -> Vec<Match> {
        maximal::find_maximal_matches(self, i, j)
    }

    pub(crate) fn all_positive_matches(&self, i: usize, j: usize) -> Vec<Match> {
        maximal::find_all_positive_matches(self, i, j)
    }

    pub fn lcp_between(&self, a: usize, b: usize) -> usize {
        if a == b {
            return self.concat.len() - a;
        }
        let ra = self.rank[a];
        let rb = self.rank[b];
        let (lo, hi) = if ra < rb { (ra + 1, rb) } else { (rb + 1, ra) };
        self.rmq.query(lo, hi)
    }

    fn concatenate(strings: &[Vec<i8>]) -> (Vec<u32>, Vec<usize>, Vec<usize>) {
        let n = strings.len();

        let mut chars: Vec<i8> = strings.iter().flat_map(|s| s.iter().copied()).collect();
        chars.sort();
        chars.dedup();

        let char_offset = n as u32;
        let char_map: HashMap<i8, u32> = chars
            .into_iter()
            .enumerate()
            .map(|(i, c)| (c, char_offset + i as u32))
            .collect();

        let total_len: usize = strings.iter().map(|s| s.len()).sum::<usize>() + n;
        let mut concat = Vec::with_capacity(total_len);
        let mut boundaries = Vec::with_capacity(n);
        let mut string_lens = Vec::with_capacity(n);

        for (i, s) in strings.iter().enumerate() {
            boundaries.push(concat.len());
            string_lens.push(s.len());
            for &c in s {
                concat.push(char_map[&c]);
            }
            concat.push((n - 1 - i) as u32);
        }

        (concat, boundaries, string_lens)
    }

    fn build_rank(sa: &[usize]) -> Vec<usize> {
        let mut rank = vec![0usize; sa.len()];
        for (i, &pos) in sa.iter().enumerate() {
            rank[pos] = i;
        }
        rank
    }

    fn build_doc_arrays(
        n: usize,
        boundaries: &[usize],
        string_lens: &[usize],
    ) -> (Vec<usize>, Vec<usize>) {
        let mut doc = vec![0usize; n];
        let mut doc_pos = vec![0usize; n];
        for (sid, (&start, &len)) in boundaries.iter().zip(string_lens.iter()).enumerate() {
            for i in 0..=len {
                doc[start + i] = sid;
                doc_pos[start + i] = i;
            }
        }
        (doc, doc_pos)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_string() {
        let idx = MatchIndex::new(&[vec![1i8, 2, 3]]);
        assert_eq!(idx.num_strings(), 1);
        assert_eq!(idx.string_len(0), 3);
    }

    #[test]
    fn two_strings() {
        let idx = MatchIndex::new(&[vec![1i8, 2], vec![3i8, 4]]);
        assert_eq!(idx.num_strings(), 2);
        assert_eq!(idx.string_len(0), 2);
        assert_eq!(idx.string_len(1), 2);
    }

    #[test]
    fn empty_strings() {
        let idx = MatchIndex::new(&[vec![], vec![1i8], vec![]]);
        assert_eq!(idx.num_strings(), 3);
        assert_eq!(idx.string_len(0), 0);
        assert_eq!(idx.string_len(1), 1);
        assert_eq!(idx.string_len(2), 0);
    }

    #[test]
    fn concat_last_is_sentinel_zero() {
        let idx = MatchIndex::new(&[vec![1i8, 2], vec![3i8]]);
        assert_eq!(*idx.concat.last().unwrap(), 0);
    }

    #[test]
    fn concat_sentinels_unique_and_below_chars() {
        let strings: Vec<Vec<i8>> = vec![vec![5i8, 6], vec![7i8], vec![8i8, 9, 10]];
        let idx = MatchIndex::new(&strings);
        let n = strings.len();
        let sentinel_positions: Vec<usize> = strings
            .iter()
            .enumerate()
            .map(|(i, s)| idx.boundaries[i] + s.len())
            .collect();
        let sentinels: Vec<u32> = sentinel_positions.iter().map(|&p| idx.concat[p]).collect();
        let char_offset = n as u32;
        for &s in &sentinels {
            assert!(s < char_offset, "sentinel {s} >= char_offset {char_offset}");
        }
        let char_values: Vec<u32> = idx
            .concat
            .iter()
            .enumerate()
            .filter(|(i, _)| !sentinel_positions.contains(i))
            .map(|(_, &v)| v)
            .collect();
        for &c in &char_values {
            assert!(c >= char_offset, "char {c} < char_offset {char_offset}");
        }
    }

    #[test]
    fn rank_is_inverse_of_sa() {
        let idx = MatchIndex::new(&[vec![1i8, 2, 1], vec![2i8, 1, 2]]);
        for (i, &pos) in idx.sa.iter().enumerate() {
            assert_eq!(idx.rank[pos], i);
        }
        for i in 0..idx.concat.len() {
            assert_eq!(idx.sa[idx.rank[i]], i);
        }
    }

    #[test]
    fn doc_arrays_correct() {
        let strings: Vec<Vec<i8>> = vec![vec![1i8, 2], vec![3i8]];
        let idx = MatchIndex::new(&strings);
        assert_eq!(idx.doc[idx.boundaries[0]], 0);
        assert_eq!(idx.doc[idx.boundaries[0] + 1], 0);
        assert_eq!(idx.doc[idx.boundaries[1]], 1);
        assert_eq!(idx.doc_pos[idx.boundaries[0]], 0);
        assert_eq!(idx.doc_pos[idx.boundaries[0] + 1], 1);
        assert_eq!(idx.doc_pos[idx.boundaries[1]], 0);
    }

    #[test]
    fn lcp_between_correct() {
        let strings: Vec<Vec<i8>> = vec![vec![1i8, 2, 3], vec![2i8, 3, 4]];
        let idx = MatchIndex::new(&strings);
        let pos_23_in_0 = idx.boundaries[0] + 1;
        let pos_23_in_1 = idx.boundaries[1];
        assert_eq!(idx.lcp_between(pos_23_in_0, pos_23_in_1), 2);
    }
}
