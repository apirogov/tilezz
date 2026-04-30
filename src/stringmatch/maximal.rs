use crate::stringmatch::{Match, MatchIndex};

pub fn find_maximal_matches(index: &MatchIndex, i: usize, j: usize) -> Vec<Match> {
    find_matches_with_filter(index, i, j, |index, pos_a, local_a, pos_b, local_b| {
        is_left_maximal(index, pos_a, local_a, pos_b, local_b)
    })
}

pub fn find_positive_matches_at_positions(
    index: &MatchIndex,
    i: usize,
    j: usize,
    positions: &[usize],
) -> Vec<Match> {
    find_matches_at_positions_with_filter(index, i, j, positions, |_, _, _, _, _| true)
}

fn find_matches_with_filter(
    index: &MatchIndex,
    i: usize,
    j: usize,
    accept: impl Fn(&MatchIndex, usize, usize, usize, usize) -> bool,
) -> Vec<Match> {
    let len_i = index.string_lens[i];
    let positions: Vec<usize> = (0..len_i).collect();
    find_matches_at_positions_with_filter(index, i, j, &positions, accept)
}

fn find_matches_at_positions_with_filter(
    index: &MatchIndex,
    i: usize,
    j: usize,
    positions: &[usize],
    accept: impl Fn(&MatchIndex, usize, usize, usize, usize) -> bool,
) -> Vec<Match> {
    let mut matches = Vec::new();

    let start_i = index.boundaries[i];

    for &local_a in positions {
        let pos_a = start_i + local_a;
        let r = index.rank[pos_a];

        scan_side(
            index,
            i,
            j,
            pos_a,
            local_a,
            (0..r).rev(),
            |k| index.lcp[k + 1],
            &accept,
            &mut matches,
        );

        scan_side(
            index,
            i,
            j,
            pos_a,
            local_a,
            (r + 1)..index.concat.len(),
            |k| index.lcp[k],
            &accept,
            &mut matches,
        );
    }

    matches.sort();
    matches.dedup();
    matches
}

#[allow(clippy::too_many_arguments)]
fn scan_side(
    index: &MatchIndex,
    string_i: usize,
    string_j: usize,
    pos_a: usize,
    local_a: usize,
    range: impl Iterator<Item = usize>,
    lcp_at: impl Fn(usize) -> usize,
    accept: &dyn Fn(&MatchIndex, usize, usize, usize, usize) -> bool,
    matches: &mut Vec<Match>,
) {
    let mut running_min = usize::MAX;

    for k in range {
        running_min = running_min.min(lcp_at(k));
        if running_min == 0 {
            break;
        }
        if index.doc[index.sa[k]] != string_j {
            continue;
        }
        let pos_b = index.sa[k];
        let local_b = index.doc_pos[pos_b];
        if string_i == string_j && local_a == local_b {
            continue;
        }
        if accept(index, pos_a, local_a, pos_b, local_b) {
            matches.push(Match::new(
                string_i,
                string_j,
                local_a,
                local_b,
                running_min,
            ));
        }
    }
}

fn is_left_maximal(
    index: &MatchIndex,
    pos_a: usize,
    local_a: usize,
    pos_b: usize,
    local_b: usize,
) -> bool {
    if local_a == 0 || local_b == 0 {
        return true;
    }
    index.concat[pos_a - 1] != index.concat[pos_b - 1]
}

#[cfg(test)]
mod tests {
    use crate::stringmatch::MatchIndex;

    fn naive_maximal_matches(a: &[i8], b: &[i8]) -> Vec<(usize, usize, usize)> {
        let mut result = Vec::new();
        for ia in 0..a.len() {
            for ib in 0..b.len() {
                let mut l = 0;
                while ia + l < a.len() && ib + l < b.len() && a[ia + l] == b[ib + l] {
                    l += 1;
                }
                if l == 0 {
                    continue;
                }
                let left_ok = ia == 0 || ib == 0 || a[ia - 1] != b[ib - 1];
                if left_ok {
                    result.push((ia, ib, l));
                }
            }
        }
        result.sort();
        result.dedup();
        result
    }

    fn normalize_matches(tuples: &mut Vec<(usize, usize, usize)>, same_string: bool) {
        if same_string {
            tuples.retain(|&(a, b, _)| a != b);
            for t in tuples.iter_mut() {
                let (a, b, l) = *t;
                *t = (a.min(b), a.max(b), l);
            }
        }
        tuples.sort();
        tuples.dedup();
    }

    fn check_matches(strings: &[Vec<i8>], i: usize, j: usize) {
        let idx = MatchIndex::new(strings);
        let matches = idx.maximal_matches(i, j);

        let mut got: Vec<(usize, usize, usize)> = matches
            .iter()
            .map(|m| (m.a.start, m.b.start, m.a.len))
            .collect();
        let mut expected = naive_maximal_matches(&strings[i], &strings[j]);
        let same = i == j;
        normalize_matches(&mut got, same);
        normalize_matches(&mut expected, same);

        assert_eq!(got, expected, "mismatch for strings[{i}] vs strings[{j}]");
    }

    #[test]
    fn no_matches() {
        check_matches(&[vec![1i8, 2], vec![3i8, 4]], 0, 1);
    }

    #[test]
    fn full_match() {
        check_matches(&[vec![1i8, 2, 3], vec![1i8, 2, 3]], 0, 1);
    }

    #[test]
    fn partial_overlap() {
        check_matches(&[vec![1i8, 2, 3, 4], vec![3i8, 4, 5]], 0, 1);
    }

    #[test]
    fn partial_overlap_reversed() {
        check_matches(&[vec![3i8, 4, 5], vec![1i8, 2, 3, 4]], 0, 1);
    }

    #[test]
    fn multiple_overlaps() {
        check_matches(&[vec![1i8, 2, 1, 2], vec![1i8, 2, 3, 1]], 0, 1);
    }

    #[test]
    fn same_string_self_matches() {
        check_matches(&[vec![1i8, 2, 1, 2]], 0, 0);
    }

    #[test]
    fn single_char_strings() {
        check_matches(&[vec![5i8], vec![5i8]], 0, 1);
        check_matches(&[vec![5i8], vec![6i8]], 0, 1);
    }

    #[test]
    fn empty_string_involved() {
        check_matches(&[vec![], vec![1i8, 2]], 0, 1);
        check_matches(&[vec![1i8, 2], vec![]], 0, 1);
    }

    #[test]
    fn all_same_string() {
        check_matches(&[vec![1i8, 1, 1]], 0, 0);
    }

    #[test]
    fn three_strings_pairwise() {
        let strings: Vec<Vec<i8>> = vec![vec![1i8, 2, 3], vec![2i8, 3, 4], vec![5i8, 1, 2]];
        check_matches(&strings, 0, 1);
        check_matches(&strings, 0, 2);
        check_matches(&strings, 1, 2);
    }

    #[test]
    fn long_repeated_pattern() {
        let s: Vec<i8> = [1i8, 2].iter().cycle().take(20).copied().collect();
        check_matches(&[s.clone(), s.clone()], 0, 1);
    }

    #[test]
    fn vs_naive_random() {
        let mut seed: u64 = 98765;
        for _ in 0..30 {
            let strings = (0..3)
                .map(|_| {
                    (0..8)
                        .map(|_| {
                            seed = seed
                                .wrapping_mul(6364136223846793005)
                                .wrapping_add(1442695040888963407);
                            ((seed % 5) as i8) - 2
                        })
                        .collect()
                })
                .collect::<Vec<Vec<i8>>>();
            for i in 0..strings.len() {
                for j in i..strings.len() {
                    check_matches(&strings, i, j);
                }
            }
        }
    }

    #[test]
    fn negative_values() {
        check_matches(&[vec![-1i8, -2, -3], vec![-2i8, -3, -4]], 0, 1);
    }
}
