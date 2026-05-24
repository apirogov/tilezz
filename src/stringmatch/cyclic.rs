//! Linear-time cyclic substring containment.

use super::period::prefix_function;

/// Returns `true` iff `needle` appears as a contiguous **cyclic**
/// substring of `haystack` — i.e. there exists some starting offset
/// `s` such that `needle[i] == haystack[(s + i) % n]` for all
/// `0 <= i < needle.len()`.
///
/// Empty needles return `false`; needles longer than `haystack` return
/// `false` (they can't fit cyclically).
///
/// O(`haystack.len()` + `needle.len()`) time via KMP run against
/// `haystack` walked cyclically for `n + m - 1` steps — that range is
/// sufficient to surface every cyclic alignment of `needle`.
pub fn cyclic_contains<T: Eq>(haystack: &[T], needle: &[T]) -> bool {
    let n = haystack.len();
    let m = needle.len();
    if m == 0 || m > n {
        return false;
    }

    let pi = prefix_function(needle);

    let total = n + m - 1;
    let mut k = 0usize;
    for i in 0..total {
        let c = &haystack[i % n];
        while k > 0 && needle[k] != *c {
            k = pi[k - 1];
        }
        if needle[k] == *c {
            k += 1;
        }
        if k == m {
            return true;
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    // Brute-force O(n*m) reference for cross-checking.
    fn naive<T: Eq>(haystack: &[T], needle: &[T]) -> bool {
        let n = haystack.len();
        let m = needle.len();
        if m == 0 || m > n {
            return false;
        }
        (0..n).any(|s| (0..m).all(|i| haystack[(s + i) % n] == needle[i]))
    }

    #[test]
    fn empty_needle_is_false() {
        assert!(!cyclic_contains::<i8>(&[1, 2, 3], &[]));
    }

    #[test]
    fn empty_haystack_is_false() {
        assert!(!cyclic_contains::<i8>(&[], &[1]));
    }

    #[test]
    fn needle_longer_than_haystack_is_false() {
        assert!(!cyclic_contains(&[1i8, 2], &[1i8, 2, 3]));
    }

    #[test]
    fn linear_substring() {
        assert!(cyclic_contains(&[1i8, 2, 3, 4, 5], &[2i8, 3, 4]));
    }

    #[test]
    fn wrapping_substring() {
        // [1,2,3,4,5] cyclically contains [4,5,1,2]
        assert!(cyclic_contains(&[1i8, 2, 3, 4, 5], &[4i8, 5, 1, 2]));
    }

    #[test]
    fn full_length_rotation() {
        // Any rotation of haystack is itself a cyclic substring at full length.
        assert!(cyclic_contains(&[1i8, 2, 3], &[3i8, 1, 2]));
        assert!(cyclic_contains(&[1i8, 2, 3], &[2i8, 3, 1]));
        assert!(cyclic_contains(&[1i8, 2, 3], &[1i8, 2, 3]));
    }

    #[test]
    fn full_length_non_rotation() {
        assert!(!cyclic_contains(&[1i8, 2, 3], &[1i8, 3, 2]));
    }

    #[test]
    fn missing_substring() {
        assert!(!cyclic_contains(&[1i8, 2, 3, 4, 5], &[2i8, 4]));
    }

    #[test]
    fn single_element_needle() {
        assert!(cyclic_contains(&[1i8, 2, 3], &[2i8]));
        assert!(!cyclic_contains(&[1i8, 2, 3], &[4i8]));
    }

    #[test]
    fn repeating_haystack() {
        // haystack is "aaaaaa", needle "aaa" — must match anywhere.
        assert!(cyclic_contains(&[1u8; 6], &[1u8; 3]));
        // haystack all-same, needle includes a different element — no.
        assert!(!cyclic_contains(&[1u8; 6], &[1u8, 1, 2]));
    }

    #[test]
    fn does_not_double_count_wrap() {
        // haystack length n, needle length n-1: only n distinct cyclic
        // starting positions are possible. None of these spuriously
        // matches a non-rotation needle.
        assert!(!cyclic_contains(&[1i8, 2, 3, 4], &[1i8, 3, 4]));
    }

    #[test]
    fn vs_naive_random() {
        let mut s: u64 = 0xc0ffee;
        let next = |state: &mut u64| {
            *state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            *state
        };
        for _ in 0..500 {
            let n = (next(&mut s) as usize % 30) + 1;
            let m = (next(&mut s) as usize % (n + 2)) + 1;
            let haystack: Vec<u8> = (0..n).map(|_| (next(&mut s) % 3) as u8).collect();
            let needle: Vec<u8> = (0..m).map(|_| (next(&mut s) % 3) as u8).collect();
            assert_eq!(
                cyclic_contains(&haystack, &needle),
                naive(&haystack, &needle),
                "haystack={haystack:?} needle={needle:?}"
            );
        }
    }
}
