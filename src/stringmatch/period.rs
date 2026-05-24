pub fn repetition_factor<T: Eq>(s: &[T]) -> usize {
    let n = s.len();
    if n <= 1 {
        return 1;
    }
    let pi = prefix_function(s);
    let period = n - pi[n - 1];
    if n.is_multiple_of(period) {
        n / period
    } else {
        1
    }
}

pub(crate) fn prefix_function<T: Eq>(s: &[T]) -> Vec<usize> {
    let n = s.len();
    let mut pi = vec![0usize; n];
    for i in 1..n {
        let mut j = pi[i - 1];
        while j > 0 && s[i] != s[j] {
            j = pi[j - 1];
        }
        if s[i] == s[j] {
            j += 1;
        }
        pi[i] = j;
    }
    pi
}

#[cfg(test)]
mod tests {
    use super::*;

    fn naive_repetition_factor<T: Eq>(s: &[T]) -> usize {
        let n = s.len();
        if n <= 1 {
            return 1;
        }
        for period in 1..n {
            if !n.is_multiple_of(period) {
                continue;
            }
            let pattern = &s[..period];
            if (0..n).all(|i| s[i] == pattern[i % period]) {
                return n / period;
            }
        }
        1
    }

    fn check<T: Eq + std::fmt::Debug>(s: &[T], expected: usize) {
        let got = repetition_factor(s);
        assert_eq!(got, expected, "input={s:?}");
        let naive = naive_repetition_factor(s);
        assert_eq!(got, naive, "input={s:?}");
    }

    #[test]
    fn empty_string() {
        check::<u8>(&[], 1);
    }

    #[test]
    fn single_char() {
        check(&[1u32], 1);
    }

    #[test]
    fn two_same() {
        check(&[1u32, 1], 2);
    }

    #[test]
    fn two_different() {
        check(&[1u32, 2], 1);
    }

    #[test]
    fn ababab() {
        check(&[1u32, 2, 1, 2, 1, 2], 3);
    }

    #[test]
    fn ababac() {
        check(&[1u32, 2, 1, 2, 1, 3], 1);
    }

    #[test]
    fn abcabc() {
        check(&[1u32, 2, 3, 1, 2, 3], 2);
    }

    #[test]
    fn all_same_short() {
        check(&[5u32, 5, 5, 5], 4);
    }

    #[test]
    fn all_same_long() {
        let s: Vec<u32> = vec![7u32; 100];
        check(&s, 100);
    }

    #[test]
    fn single_repetition() {
        check(&[1u32, 2, 3, 4, 5], 1);
    }

    #[test]
    fn four_repetitions() {
        check(&[1u32, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3], 4);
    }

    #[test]
    fn not_quite_periodic() {
        check(&[1u32, 2, 1, 2, 1, 3], 1);
    }

    #[test]
    fn prefix_function_correctness() {
        let pi = prefix_function(&[1u32, 2, 1, 2, 1, 2]);
        assert_eq!(pi, vec![0, 0, 1, 2, 3, 4]);
    }

    #[test]
    fn prefix_function_no_overlap() {
        let pi = prefix_function(&[1u32, 2, 3, 4, 5]);
        assert_eq!(pi, vec![0, 0, 0, 0, 0]);
    }

    #[test]
    fn prefix_function_full_overlap() {
        let pi = prefix_function(&[7u32, 7, 7, 7]);
        assert_eq!(pi, vec![0, 1, 2, 3]);
    }

    #[test]
    fn vs_naive_random() {
        let mut s: u64 = 42;
        for _ in 0..200 {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let n = (s as usize % 40) + 2;
            let text: Vec<u32> = (0..n)
                .map(|_| {
                    s = s
                        .wrapping_mul(6364136223846793005)
                        .wrapping_add(1442695040888963407);
                    (s % 4) as u32
                })
                .collect();
            let got = repetition_factor(&text);
            let expected = naive_repetition_factor(&text);
            assert_eq!(got, expected, "text={text:?}");
        }
    }

    #[test]
    fn works_with_i8() {
        assert_eq!(repetition_factor(&[1i8, 2, 1, 2, 1, 2]), 3);
        assert_eq!(repetition_factor(&[1i8, 2, 1, 2, 1, 3]), 1);
    }

    #[test]
    fn works_with_u8_bytes() {
        assert_eq!(repetition_factor(b"abcabc"), 2);
        assert_eq!(repetition_factor(b"ababab"), 3);
        assert_eq!(repetition_factor(b"aaaa"), 4);
        assert_eq!(repetition_factor(b"abcd"), 1);
    }
}
