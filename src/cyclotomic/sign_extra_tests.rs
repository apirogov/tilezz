// Extra tests for cyclotomic signum helpers.
// Kept in a separate module to avoid cluttering sign.rs.

#[cfg(test)]
mod extra_tests {
    use crate::cyclotomic::sign::{
        signum_sum_sqrt_expr_2, signum_sum_sqrt_expr_4, signum_sum_sqrt_expr_4_pentagonal,
    };

    // Slow reference using f64, only for tests.
    fn sign_f64(x: f64) -> i64 {
        if x == 0.0 {
            0
        } else if x > 0.0 {
            1
        } else {
            -1
        }
    }

    #[test]
    fn test_signum_sum_sqrt_expr_2_matches_f64_on_small_range() {
        // cover sign changes, zeros, and asymmetry
        for a in -10..=10 {
            for b in -10..=10 {
                for (m, n) in [(2f64, 3f64), (2f64, 5f64), (3f64, 7f64)] {
                    let got = signum_sum_sqrt_expr_2(a, m as i64, b, n as i64);
                    let exp = sign_f64((a as f64) * m.sqrt() + (b as f64) * n.sqrt());
                    assert_eq!(got, exp, "a={a} b={b} m={m} n={n}");
                }
            }
        }
    }

    #[test]
    fn test_signum_sum_sqrt_expr_4_matches_f64_on_small_range() {
        // parameters modeled after ZZ24: sqrt(1), sqrt(2), sqrt(3), sqrt(6)
        for a in -6..=6 {
            for b in -6..=6 {
                for c in -6..=6 {
                    for d in -6..=6 {
                        let got = signum_sum_sqrt_expr_4(a, 1, b, 2, c, 3, d, 6);
                        let val = (a as f64)
                            + (b as f64) * 2f64.sqrt()
                            + (c as f64) * 3f64.sqrt()
                            + (d as f64) * 6f64.sqrt();
                        let exp = sign_f64(val);
                        assert_eq!(got, exp, "a={a} b={b} c={c} d={d}");
                    }
                }
            }
        }
    }

    #[test]
    fn test_signum_sum_sqrt_expr_4_pentagonal_matches_f64_on_small_range() {
        let sqrt5: f64 = 2.236_067_977_499_79;
        let y: f64 = 10.0 - 2.0 * sqrt5;
        for a in -6..=6 {
            for b in -6..=6 {
                for c in -6..=6 {
                    for d in -6..=6 {
                        let got = signum_sum_sqrt_expr_4_pentagonal(a, b, c, d);
                        let val = (a as f64)
                            + (b as f64) * sqrt5
                            + (c as f64) * y.sqrt()
                            + (d as f64) * (5.0 * y).sqrt();
                        let exp = sign_f64(val);
                        assert_eq!(got, exp, "a={a} b={b} c={c} d={d}");
                    }
                }
            }
        }
    }
}
