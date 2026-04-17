pub fn build_lcp_array(text: &[u32], sa: &[usize]) -> Vec<usize> {
    let n = text.len();
    assert_eq!(sa.len(), n, "SA length must equal text length");

    if n == 0 {
        return vec![];
    }

    let rank = build_rank(sa);
    let mut lcp = vec![0usize; n];
    let mut k: usize = 0;

    for i in 0..n {
        if rank[i] == 0 {
            k = 0;
            continue;
        }
        let j = sa[rank[i] - 1];
        while i + k < n && j + k < n && text[i + k] == text[j + k] {
            k += 1;
        }
        lcp[rank[i]] = k;
        k = k.saturating_sub(1);
    }

    lcp
}

fn build_rank(sa: &[usize]) -> Vec<usize> {
    let n = sa.len();
    let mut rank = vec![0usize; n];
    for (i, &pos) in sa.iter().enumerate() {
        rank[pos] = i;
    }
    rank
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stringmatch::sais::build_suffix_array;

    fn naive_lcp(text: &[u32], sa: &[usize]) -> Vec<usize> {
        let n = text.len();
        let mut lcp = vec![0usize; n];
        for i in 1..n {
            let a = sa[i - 1];
            let b = sa[i];
            let mut k = 0;
            while a + k < n && b + k < n && text[a + k] == text[b + k] {
                k += 1;
            }
            lcp[i] = k;
        }
        lcp
    }

    fn build_and_check(text: &[u32]) -> Vec<usize> {
        let sa = build_suffix_array(text);
        let lcp = build_lcp_array(text, &sa);
        let expected = naive_lcp(text, &sa);
        assert_eq!(lcp, expected, "text={text:?}");
        lcp
    }

    #[test]
    fn single_sentinel() {
        let text = vec![0u32];
        let sa = build_suffix_array(&text);
        let lcp = build_lcp_array(&text, &sa);
        assert_eq!(lcp, vec![0]);
    }

    #[test]
    fn two_chars() {
        build_and_check(&[1, 0]);
    }

    #[test]
    fn banana() {
        let text: Vec<u32> = vec![2, 1, 3, 1, 3, 1, 0];
        let sa = build_suffix_array(&text);
        let lcp = build_lcp_array(&text, &sa);
        assert_eq!(sa, vec![6, 5, 3, 1, 0, 4, 2]);
        assert_eq!(lcp, vec![0, 0, 1, 3, 0, 0, 2]);
    }

    #[test]
    fn all_same() {
        build_and_check(&[1, 1, 1, 0]);
    }

    #[test]
    fn already_sorted() {
        build_and_check(&[1, 2, 3, 0]);
    }

    #[test]
    fn reverse_sorted() {
        build_and_check(&[3, 2, 1, 0]);
    }

    #[test]
    fn ababab() {
        build_and_check(&[1, 2, 1, 2, 1, 2, 0]);
    }

    #[test]
    fn two_equal_blocks() {
        build_and_check(&[1, 2, 3, 1, 2, 3, 0]);
    }

    #[test]
    fn long_repeated_pattern() {
        let mut text = vec![];
        for _ in 0..20 {
            text.push(1u32);
            text.push(2u32);
        }
        text.push(0);
        build_and_check(&text);
    }

    #[test]
    fn vs_naive_random() {
        let mut s: u64 = 77777;
        for _ in 0..50 {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let n = (s as usize % 80) + 5;
            let text: Vec<u32> = (0..n)
                .map(|i| {
                    if i == n - 1 {
                        0
                    } else {
                        s = s
                            .wrapping_mul(6364136223846793005)
                            .wrapping_add(1442695040888963407);
                        (s % 10) as u32 + 1
                    }
                })
                .collect();
            build_and_check(&text);
        }
    }

    #[test]
    fn rank_is_inverse_of_sa() {
        let text: Vec<u32> = vec![2, 1, 3, 1, 3, 1, 0];
        let sa = build_suffix_array(&text);
        let rank = build_rank(&sa);
        for (i, &pos) in sa.iter().enumerate() {
            assert_eq!(rank[pos], i, "sa[{i}]={pos} but rank[{pos}]={}", rank[pos]);
        }
        for i in 0..text.len() {
            assert_eq!(sa[rank[i]], i, "sa[rank[{i}]] != {i}");
        }
    }
}
