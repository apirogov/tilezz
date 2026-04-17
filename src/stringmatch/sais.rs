//! SA-IS (Suffix Array Induced Sorting) for integer alphabets.
//!
//! Reference: Nong, Zhang, Chan -
//! "Two Efficient Algorithms for Linear Time Suffix Array Construction" (2011)
//!
//! Input: a sequence of `u32` values in `[0, K)` where `K` is the alphabet size.
//! The last element **must** be `0` (the sentinel, which is the unique minimum).

const EMPTY: usize = usize::MAX;

pub fn build_suffix_array(text: &[u32]) -> Vec<usize> {
    assert!(!text.is_empty(), "text must not be empty");
    assert_eq!(
        *text.last().unwrap(),
        0,
        "last element must be sentinel (0)"
    );

    let n = text.len();
    if n == 1 {
        return vec![0];
    }

    let alphabet_size = (*text.iter().max().unwrap()) as usize + 1;
    let mut sa = vec![EMPTY; n];
    sais(text, &mut sa, alphabet_size);
    sa
}

fn sais(text: &[u32], sa: &mut [usize], alphabet_size: usize) {
    let n = text.len();

    if n == 1 {
        sa[0] = 0;
        return;
    }

    let stypes = classify_types(text);

    sa.fill(EMPTY);
    place_lms_right_to_left(text, sa, &stypes, alphabet_size);
    induce_l(text, sa, &stypes, alphabet_size);
    induce_s(text, sa, &stypes, alphabet_size);

    let n1 = compact_sorted_lms(sa, &stypes);
    let sorted_lms = &sa[..n1];

    let (names, num_names) = name_lms_substrings(text, &stypes, sorted_lms);

    let mut pos_to_name = vec![0u32; n];
    for j in 0..n1 {
        pos_to_name[sorted_lms[j]] = names[j];
    }

    let lms_in_order: Vec<usize> = (0..n).filter(|&i| is_lms(&stypes, i)).collect();
    let s1: Vec<u32> = lms_in_order.iter().map(|&p| pos_to_name[p]).collect();

    let sorted_lms_final = if num_names as usize == n1 {
        sorted_lms.to_vec()
    } else {
        let k1 = num_names as usize;
        let mut sa1 = vec![EMPTY; n1];
        sais(&s1, &mut sa1, k1);
        sa1.iter().map(|&rank| lms_in_order[rank]).collect()
    };

    sa.fill(EMPTY);
    place_lms_sorted(text, sa, alphabet_size, &sorted_lms_final);
    induce_l(text, sa, &stypes, alphabet_size);
    induce_s(text, sa, &stypes, alphabet_size);
}

fn classify_types(text: &[u32]) -> Vec<bool> {
    let n = text.len();
    let mut stypes = vec![false; n];
    stypes[n - 1] = true;
    for i in (0..n - 1).rev() {
        if text[i] < text[i + 1] {
            stypes[i] = true;
        } else if text[i] == text[i + 1] {
            stypes[i] = stypes[i + 1];
        }
    }
    stypes
}

fn is_lms(stypes: &[bool], i: usize) -> bool {
    i > 0 && stypes[i] && !stypes[i - 1]
}

fn place_lms_right_to_left(text: &[u32], sa: &mut [usize], stypes: &[bool], alphabet_size: usize) {
    let mut bk = Buckets::from_text(text, alphabet_size);
    for i in (1..text.len()).rev() {
        if is_lms(stypes, i) {
            let c = text[i] as usize;
            sa[bk.tail(c)] = i;
            bk.decr_tail(c);
        }
    }
}

fn place_lms_sorted(text: &[u32], sa: &mut [usize], alphabet_size: usize, sorted: &[usize]) {
    let mut bk = Buckets::from_text(text, alphabet_size);
    for &pos in sorted.iter().rev() {
        let c = text[pos] as usize;
        sa[bk.tail(c)] = pos;
        bk.decr_tail(c);
    }
}

fn induce_l(text: &[u32], sa: &mut [usize], stypes: &[bool], alphabet_size: usize) {
    let n = sa.len();
    let mut bk = Buckets::from_text(text, alphabet_size);
    for i in 0..n {
        if sa[i] == EMPTY {
            continue;
        }
        let j = sa[i];
        if j == 0 {
            continue;
        }
        if !stypes[j - 1] {
            let c = text[j - 1] as usize;
            sa[bk.head(c)] = j - 1;
            bk.incr_head(c);
        }
    }
}

fn induce_s(text: &[u32], sa: &mut [usize], stypes: &[bool], alphabet_size: usize) {
    let n = sa.len();
    let mut bk = Buckets::from_text(text, alphabet_size);
    for i in (0..n).rev() {
        if sa[i] == EMPTY {
            continue;
        }
        let j = sa[i];
        if j == 0 {
            continue;
        }
        if stypes[j - 1] {
            let c = text[j - 1] as usize;
            sa[bk.tail(c)] = j - 1;
            bk.decr_tail(c);
        }
    }
}

fn compact_sorted_lms(sa: &mut [usize], stypes: &[bool]) -> usize {
    let n = sa.len();
    let mut n1 = 0;
    for i in 0..n {
        if sa[i] != EMPTY && is_lms(stypes, sa[i]) {
            sa[n1] = sa[i];
            n1 += 1;
        }
    }
    for slot in sa.iter_mut().take(n).skip(n1) {
        *slot = EMPTY;
    }
    n1
}

fn name_lms_substrings(text: &[u32], stypes: &[bool], sorted_lms: &[usize]) -> (Vec<u32>, u32) {
    let n1 = sorted_lms.len();
    let mut names = vec![0u32; n1];
    let mut num_names = 1u32;
    for j in 1..n1 {
        if !lms_substrings_equal(text, stypes, sorted_lms[j - 1], sorted_lms[j]) {
            num_names += 1;
        }
        names[j] = num_names - 1;
    }
    (names, num_names)
}

fn lms_substrings_equal(text: &[u32], stypes: &[bool], p1: usize, p2: usize) -> bool {
    let n = text.len();
    let mut d = 0;
    loop {
        let a = p1 + d;
        let b = p2 + d;
        if a >= n || b >= n {
            return a >= n && b >= n;
        }
        if text[a] != text[b] || stypes[a] != stypes[b] {
            return false;
        }
        if d > 0 && is_lms(stypes, a) {
            return is_lms(stypes, b);
        }
        d += 1;
    }
}

struct Buckets {
    head: Vec<usize>,
    tail: Vec<usize>,
}

impl Buckets {
    fn from_text(text: &[u32], alphabet_size: usize) -> Self {
        let mut count = vec![0usize; alphabet_size];
        for &c in text {
            count[c as usize] += 1;
        }
        let mut head = vec![0usize; alphabet_size];
        let mut tail = vec![0usize; alphabet_size];
        let mut sum = 0;
        for i in 0..alphabet_size {
            head[i] = sum;
            sum += count[i];
            tail[i] = sum.saturating_sub(1);
        }
        Buckets { head, tail }
    }

    fn head(&self, c: usize) -> usize {
        self.head[c]
    }

    fn tail(&self, c: usize) -> usize {
        self.tail[c]
    }

    fn incr_head(&mut self, c: usize) {
        self.head[c] += 1;
    }

    fn decr_tail(&mut self, c: usize) {
        self.tail[c] = self.tail[c].saturating_sub(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn naive_sa(text: &[u32]) -> Vec<usize> {
        let mut sa: Vec<usize> = (0..text.len()).collect();
        sa.sort_by(|&a, &b| text[a..].cmp(&text[b..]));
        sa
    }

    #[test]
    fn single_sentinel() {
        assert_eq!(build_suffix_array(&[0]), vec![0]);
    }

    #[test]
    fn two_chars() {
        assert_eq!(build_suffix_array(&[1, 0]), vec![1, 0]);
    }

    #[test]
    fn banana() {
        let text: Vec<u32> = vec![2, 1, 3, 1, 3, 1, 0];
        assert_eq!(build_suffix_array(&text), vec![6, 5, 3, 1, 0, 4, 2]);
    }

    #[test]
    fn all_same() {
        assert_eq!(build_suffix_array(&[1, 1, 1, 0]), vec![3, 2, 1, 0]);
    }

    #[test]
    fn already_sorted() {
        let text: Vec<u32> = vec![1, 2, 3, 0];
        let sa = build_suffix_array(&text);
        assert_eq!(sa, naive_sa(&text));
    }

    #[test]
    fn reverse_sorted() {
        let text: Vec<u32> = vec![3, 2, 1, 0];
        let sa = build_suffix_array(&text);
        assert_eq!(sa, naive_sa(&text));
    }

    #[test]
    fn ababab() {
        let text: Vec<u32> = vec![1, 2, 1, 2, 1, 2, 0];
        let sa = build_suffix_array(&text);
        assert_eq!(sa, naive_sa(&text));
    }

    #[test]
    fn descending_stairs() {
        let text: Vec<u32> = vec![4, 3, 2, 1, 0];
        assert_eq!(build_suffix_array(&text), vec![4, 3, 2, 1, 0]);
    }

    #[test]
    fn ascending_stairs() {
        let text: Vec<u32> = vec![1, 2, 3, 4, 0];
        assert_eq!(build_suffix_array(&text), naive_sa(&text));
    }

    #[test]
    fn two_equal_blocks() {
        let text: Vec<u32> = vec![1, 2, 3, 1, 2, 3, 0];
        assert_eq!(build_suffix_array(&text), naive_sa(&text));
    }

    #[test]
    fn single_char_before_sentinel() {
        for c in [1u32, 2, 5, 10] {
            let text = vec![c, 0];
            assert_eq!(build_suffix_array(&text), vec![1, 0], "c={c}");
        }
    }

    #[test]
    fn long_repeated_pattern() {
        let mut text = vec![];
        for _ in 0..20 {
            text.push(1u32);
            text.push(2u32);
        }
        text.push(0);
        assert_eq!(build_suffix_array(&text), naive_sa(&text));
    }

    #[test]
    fn every_char_different() {
        let text: Vec<u32> = vec![5, 3, 1, 4, 2, 0];
        assert_eq!(build_suffix_array(&text), naive_sa(&text));
    }

    #[test]
    fn vs_naive_random_small() {
        let seeds = [42u64, 137, 999, 0, u64::MAX];
        for seed in seeds {
            let mut s = seed;
            for n in [2, 3, 5, 7, 13, 21, 37] {
                let text: Vec<u32> = (0..n)
                    .map(|i| {
                        if i == n - 1 {
                            0
                        } else {
                            s = s
                                .wrapping_mul(6364136223846793005)
                                .wrapping_add(1442695040888963407);
                            (s % 5) as u32 + 1
                        }
                    })
                    .collect();
                let sais_sa = build_suffix_array(&text);
                let naive = naive_sa(&text);
                assert_eq!(sais_sa, naive, "seed={seed}, n={n}");
            }
        }
    }

    #[test]
    fn vs_naive_random_medium() {
        let mut s: u64 = 12345;
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
            let sais_sa = build_suffix_array(&text);
            let naive = naive_sa(&text);
            assert_eq!(sais_sa, naive);
        }
    }

    #[test]
    fn classify_types_known() {
        let text: Vec<u32> = vec![2, 1, 3, 1, 3, 1, 0];
        let stypes = classify_types(&text);
        assert_eq!(stypes, vec![false, true, false, true, false, false, true]);
    }

    #[test]
    fn lms_positions_banana() {
        let text: Vec<u32> = vec![2, 1, 3, 1, 3, 1, 0];
        let stypes = classify_types(&text);
        let lms: Vec<usize> = (0..text.len()).filter(|&i| is_lms(&stypes, i)).collect();
        assert_eq!(lms, vec![1, 3, 6]);
    }

    #[test]
    fn buckets_correct() {
        let text: Vec<u32> = vec![2, 1, 3, 1, 3, 1, 0];
        let bk = Buckets::from_text(&text, 4);
        assert_eq!(bk.head(0), 0);
        assert_eq!(bk.tail(0), 0);
        assert_eq!(bk.head(1), 1);
        assert_eq!(bk.tail(1), 3);
        assert_eq!(bk.head(2), 4);
        assert_eq!(bk.tail(2), 4);
        assert_eq!(bk.head(3), 5);
        assert_eq!(bk.tail(3), 6);
    }
}
