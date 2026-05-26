//! Floating-point-free sign extraction for sums of square-root expressions.
//!
//! Every helper operates on bare integer coefficients (typed `T: IntRing +
//! ZSigned`) and returns the sign of the algebraic expression. Each
//! cyclotomic ring picks the helper that matches the algebraic structure
//! of its real subring; concretely:
//!
//! * `signum_sum_sqrt_expr_2`         -- ZZ8 (`Q(sqrt(2))`)
//! * `signum_sum_sqrt_expr_4`         -- ZZ24 (`Q(sqrt(2), sqrt(3))`)
//! * `signum_sum_sqrt_expr_4_pentagonal` -- ZZ10, ZZ20 (`Q(sqrt(5), sqrt(10-2sqrt(5)))`)
//! * `signum_sum_sqrt_expr_4_zz16`    -- ZZ16 (`Q(sqrt(2), sqrt(2+sqrt(2)))`)
//! * `signum_sum_sqrt_expr_8_zz32`    -- ZZ32 (`Q(sqrt(2), sqrt(2+sqrt(2)), sqrt(2+sqrt(2+sqrt(2))))`)
//! * `signum_sum_sqrt_expr_8_zz60`    -- ZZ60 (`Q(sqrt(3), sqrt(5), sqrt(10-2sqrt(5)))`)
//!
//! ZZ4 and ZZ12 inline their own sign extraction (trivial integer sign and
//! `sign_m_plus_n_sqrt3` respectively).

use num_traits::FromPrimitive;

use super::numtraits::{IntRing, ZSigned};

/// Floating-point-free solution to get sign of an expression
/// a*sqrt(n) + b*sqrt(m)
/// where a,b,m,n are all integers
/// and m,n are squarefree coprime constants > 1.
pub fn signum_sum_sqrt_expr_2<T: IntRing + ZSigned>(a: T, m: T, b: T, n: T) -> T {
    // if a and b are both positive or negative -> trivial,
    let sgn_a = a.signum();
    let sgn_b = b.signum();
    if sgn_a == sgn_b {
        return sgn_a;
    }
    if a.is_zero() {
        return sgn_b;
    }
    if b.is_zero() {
        return sgn_a;
    }
    // if a>0, b<0: z > 0 <=> n*a^2 > m*b^2 <=> n*a^2 - m*b^2 > 0
    // if a<0, b>0: z > 0 <=> n*a^2 < m*b^2 <=> m*b^2 - n*a^2 > 0
    // (and symmetrically for z < 0), which can be summarized
    // by using the sign of a and b to avoid case splitting.
    (sgn_a * a * a * m + sgn_b * b * b * n).signum()
}

/// Floating-point-free solution to get sign of an expression
/// a + b*sqrt(m) + c*sqrt(n) + d*sqrt(m*n)
/// where a,b,c,m,n are all integers
/// and m,n are squarefree coprime constants > 1.
#[allow(clippy::too_many_arguments)]
pub fn signum_sum_sqrt_expr_4<T: IntRing + ZSigned + FromPrimitive>(
    a: T,
    k: T,
    b: T,
    m: T,
    c: T,
    n: T,
    d: T,
    l: T,
) -> T {
    // reduce to 2x 2 roots case:
    // a*sqrt(k) + b*sqrt(m) + c*sqrt(n) + d*sqrt(l) > 0
    // <=> b*sqrt(m) + c*sqrt(n) > -(a*sqrt(k) + d*sqrt(l))

    // sign(a*sqrt(k) + d*sqrt(l))
    let sgn_ad_terms = signum_sum_sqrt_expr_2(a, k, d, l);
    // sign(b*sqrt(m) + c*sqrt(n))
    let sgn_bc_terms = signum_sum_sqrt_expr_2(b, m, c, n);
    // both half-expressions have same sign -> trivial
    if sgn_bc_terms == sgn_ad_terms {
        return sgn_ad_terms;
    }

    // (at least) one half-expression is zero -> trivial
    //
    // NOTE: https://qchu.wordpress.com/2009/07/02/square-roots-have-no-unexpected-linear-relationships/
    // or this question: https://math.stackexchange.com/a/30695
    // so we have: expression is zero <=> all linear coeffs are zero
    // (as we have distinct square-free roots)
    if sgn_ad_terms.is_zero() {
        return sgn_bc_terms;
    }
    if sgn_bc_terms.is_zero() {
        return sgn_ad_terms;
    }

    // now w.l.o.g. assume b/c term is pos. and a/d term neg.
    // (i.e. in the inequality both LHS and RHS positive),
    // to account for other case -> multiply result by sgn_bc_terms.

    // assume k = 1 and l = m*n
    // (typical structure of the rings we have)
    if !(k.is_one() && l == m * n) {
        panic!("Unhandled general case!");
    }

    // => use rewritten simplified expression (more efficient):
    // b*sqrt(m) + c*sqrt(n) > -(a + d*sqrt(mn))
    // <=> b^2*m + c^2*n + 2bc*sqrt(mn)
    //   > a^2 + d^2*mn + 2ad*sqrt(mn)
    // <=> b^2*m + c^2*n - d^2*mn - a^2
    //   > 2 * (ad - bc) * sqrt(mn)
    // <=> (b^2*m + c^2*n - d^2*mn - a^2)^2
    //   > 4 * mn * (ad - bc)^2
    // <=> (b^2*m + c^2*n - d^2*mn - a^2)^2 - 4mn*(ad-bc)^2 > 0
    let four = T::from_i8(4).unwrap();
    let mn = l;
    let lhs = (b * b * m) + (c * c * n) - (d * d * mn) - (a * a);
    let sq_lhs = lhs.signum() * lhs * lhs;
    let ad_m_bc = (a * d) - (b * c);
    let sq_rhs = four * mn * ad_m_bc.signum() * ad_m_bc * ad_m_bc;

    sgn_bc_terms.signum() * (sq_lhs - sq_rhs).signum()
}

/// Floating-point-free sign of a + b*sqrt(2) + c*sqrt(2+sqrt(2)) + d*sqrt(2*(2+sqrt(2)))
/// via recursive reduction from Q(sqrt(2), sqrt(2+sqrt(2))) to Q(sqrt(2)).
///
/// This is the ZZ16 real subring. Group as z = P + Q*sqrt(y) where
/// P = a + b*sqrt(2), Q = c + d*sqrt(2), y = 2 + sqrt(2). Then sign(z) is
/// determined by sign(P), sign(Q), and sign(P^2 - Q^2*y), with
/// P^2 - Q^2*y in Q(sqrt(2)). Same closed-form shape as the pentagonal
/// helper below; differs only in `y`.
pub fn signum_sum_sqrt_expr_4_zz16<T: IntRing + ZSigned + FromPrimitive>(
    a: T,
    b: T,
    c: T,
    d: T,
) -> T {
    let sp = signum_sum_sqrt_expr_2(a, T::one(), b, T::from_i8(2).unwrap());
    let sq = signum_sum_sqrt_expr_2(c, T::one(), d, T::from_i8(2).unwrap());

    if sp == sq {
        return sp;
    }
    if sq.is_zero() {
        return sp;
    }
    if sp.is_zero() {
        return sq;
    }

    let int2 = T::from_i8(2).unwrap();
    let int4 = T::from_i8(4).unwrap();

    let aa = a * a;
    let bb = b * b;
    let cc = c * c;
    let dd = d * d;
    let cd = c * d;
    let ab = a * b;

    // P^2 - Q^2 * (2 + sqrt(2)) in Z[sqrt(2)]:
    //   P^2     = a^2 + 2*b^2 + 2*a*b*sqrt(2)
    //   Q^2     = c^2 + 2*d^2 + 2*c*d*sqrt(2)
    //   Q^2 * y = (2*c^2 + 4*d^2 + 4*c*d) + (c^2 + 2*d^2 + 4*c*d) * sqrt(2)
    let alpha = aa + int2 * bb - int2 * cc - int4 * dd - int4 * cd;
    let beta = int2 * ab - cc - int2 * dd - int4 * cd;

    let spq = signum_sum_sqrt_expr_2(alpha, T::one(), beta, int2);

    -sq * spq
}

/// Floating-point-free sign of a + b*sqrt(5) + c*sqrt(10-2*sqrt(5)) + d*sqrt(5*(10-2*sqrt(5)))
/// via recursive reduction from Q(sqrt(5), sqrt(10-2*sqrt(5))) to Q(sqrt(5)).
///
/// Group as z = P + Q*sqrt(y) where P = a + b*sqrt(5), Q = c + d*sqrt(5), y = 10 - 2*sqrt(5).
/// Then sign(z) is determined by sign(P), sign(Q), and sign(P^2 - Q^2*y), with P^2 - Q^2*y in Q(sqrt(5)).
pub fn signum_sum_sqrt_expr_4_pentagonal<T: IntRing + ZSigned + FromPrimitive>(
    a: T,
    b: T,
    c: T,
    d: T,
) -> T {
    let sp = signum_sum_sqrt_expr_2(a, T::one(), b, T::from_i8(5).unwrap());
    let sq = signum_sum_sqrt_expr_2(c, T::one(), d, T::from_i8(5).unwrap());

    if sp == sq {
        return sp;
    }
    if sq.is_zero() {
        return sp;
    }

    let int2 = T::from_i8(2).unwrap();
    let int5 = T::from_i8(5).unwrap();
    let int10 = T::from_i8(10).unwrap();
    let int20 = T::from_i8(20).unwrap();
    let int50 = T::from_i8(50).unwrap();

    let aa = a * a;
    let bb = b * b;
    let cc = c * c;
    let dd = d * d;

    let alpha = aa + int5 * bb - int10 * cc + int20 * c * d - int50 * dd;
    let beta = int2 * a * b + int2 * cc - int20 * c * d + int10 * dd;

    let spq = signum_sum_sqrt_expr_2(alpha, T::one(), beta, int5);

    -sq * spq
}

/// Floating-point-free sign of
///   a0 + a1*sqrt(2) + a2*sqrt(2+sqrt(2)) + a3*sqrt(2+sqrt(2+sqrt(2)))
///     + a4*sqrt(2(2+sqrt(2))) + a5*sqrt(2(2+sqrt(2+sqrt(2))))
///     + a6*sqrt((2+sqrt(2))(2+sqrt(2+sqrt(2))))
///     + a7*sqrt(2(2+sqrt(2))(2+sqrt(2+sqrt(2))))
/// via two levels of recursive reduction.
///
/// This is the ZZ32 real subring. Group as `z = P + Q*sqrt(y)` with
///   P = a0 + a1*sqrt(2) + a2*sqrt(2+sqrt(2)) + a4*sqrt(2(2+sqrt(2)))
///   Q = a3 + a5*sqrt(2) + a6*sqrt(2+sqrt(2)) + a7*sqrt(2(2+sqrt(2)))
///   y = 2 + sqrt(2+sqrt(2))
/// Both `P, Q` live in the ZZ16 real subring `Q(sqrt(2), sqrt(2+sqrt(2)))`,
/// so `sign(P), sign(Q)` use `signum_sum_sqrt_expr_4_zz16`. The
/// inner reduction `P^2 - Q^2 * y` also lives in the ZZ16 real subring, so
/// that sign call closes the recursion (which itself reduces to
/// `Q(sqrt(2))`).
#[allow(clippy::too_many_arguments)]
pub fn signum_sum_sqrt_expr_8_zz32<T: IntRing + ZSigned + FromPrimitive>(
    a0: T,
    a1: T,
    a2: T,
    a3: T,
    a4: T,
    a5: T,
    a6: T,
    a7: T,
) -> T {
    let int2 = T::from_i8(2).unwrap();
    let int4 = T::from_i8(4).unwrap();

    let sp = signum_sum_sqrt_expr_4_zz16(a0, a1, a2, a4);
    let sq = signum_sum_sqrt_expr_4_zz16(a3, a5, a6, a7);

    if sp == sq {
        return sp;
    }
    if sq.is_zero() {
        return sp;
    }
    if sp.is_zero() {
        return sq;
    }

    // P^2 in basis {1, sqrt(2), sqrt(2+sqrt(2)), sqrt(2(2+sqrt(2)))}:
    //   const:    a0^2 + 2*a1^2 + 2*a2^2 + 4*a4^2 + 4*a2*a4
    //   sqrt(2):  a2^2 + 2*a4^2 + 2*a0*a1 + 4*a2*a4
    //   sqrt(2+sqrt(2)):       2*a0*a2 + 4*a1*a4
    //   sqrt(2(2+sqrt(2))):    2*a0*a4 + 2*a1*a2
    let p0 = a0 * a0 + int2 * a1 * a1 + int2 * a2 * a2 + int4 * a4 * a4 + int4 * a2 * a4;
    let p1 = a2 * a2 + int2 * a4 * a4 + int2 * a0 * a1 + int4 * a2 * a4;
    let p2 = int2 * a0 * a2 + int4 * a1 * a4;
    let p3 = int2 * a0 * a4 + int2 * a1 * a2;

    // Q^2 in same basis (substituting a3, a5, a6, a7 for a0, a1, a2, a4).
    let q0 = a3 * a3 + int2 * a5 * a5 + int2 * a6 * a6 + int4 * a7 * a7 + int4 * a6 * a7;
    let q1 = a6 * a6 + int2 * a7 * a7 + int2 * a3 * a5 + int4 * a6 * a7;
    let q2 = int2 * a3 * a6 + int4 * a5 * a7;
    let q3 = int2 * a3 * a7 + int2 * a5 * a6;

    // Q^2 * y where y = 2 + sqrt(2+sqrt(2)) = 2 + b2 in basis above:
    //   const:    2*q0 + 2*q2 + 2*q3   (q3*b1^2 = 2*q3 from b1*b2*b2)
    //   sqrt(2):  2*q1 + q2 + 2*q3
    //   sqrt(2+sqrt(2)):       2*q2 + q0
    //   sqrt(2(2+sqrt(2))):    2*q3 + q1
    let qy_0 = int2 * q0 + int2 * q2 + int2 * q3;
    let qy_1 = int2 * q1 + q2 + int2 * q3;
    let qy_2 = int2 * q2 + q0;
    let qy_3 = int2 * q3 + q1;

    let alpha = p0 - qy_0;
    let beta = p1 - qy_1;
    let gamma = p2 - qy_2;
    let delta = p3 - qy_3;

    let spq = signum_sum_sqrt_expr_4_zz16(alpha, beta, gamma, delta);

    -sq * spq
}

/// Floating-point-free sign of
///   a + b*sqrt(3) + c*sqrt(5) + d*sqrt(10-2*sqrt(5))
///     + e*sqrt(15) + f*sqrt(3(10-2*sqrt(5)))
///     + g*sqrt(5(10-2*sqrt(5))) + h*sqrt(15(10-2*sqrt(5)))
/// via recursive reduction from `Q(sqrt(3), sqrt(5), sqrt(10-2*sqrt(5)))` to
/// `Q(sqrt(3), sqrt(5))` (where the inner closed-form
/// `signum_sum_sqrt_expr_4` applies, since `m = 3, n = 5, l = 15 = m*n`).
///
/// This is the ZZ60 real subring. Group as `z = P + Q*sqrt(y)` with
///   P = a + b*sqrt(3) + c*sqrt(5) + e*sqrt(15)
///   Q = d + f*sqrt(3) + g*sqrt(5) + h*sqrt(15)
///   y = 10 - 2*sqrt(5)
/// Both `P, Q` live in the biquadratic field `Q(sqrt(3), sqrt(5))`. Same
/// closed-form shape as the pentagonal helper, with one extra
/// `Q(sqrt(5))` -> `Q(sqrt(3), sqrt(5))` lifting on the inner sign helper.
#[allow(clippy::too_many_arguments)]
pub fn signum_sum_sqrt_expr_8_zz60<T: IntRing + ZSigned + FromPrimitive>(
    a: T,
    b: T,
    c: T,
    d: T,
    e: T,
    f: T,
    g: T,
    h: T,
) -> T {
    let int1 = T::one();
    let int2 = T::from_i8(2).unwrap();
    let int3 = T::from_i8(3).unwrap();
    let int5 = T::from_i8(5).unwrap();
    let int6 = T::from_i8(6).unwrap();
    let int10 = T::from_i8(10).unwrap();
    let int15 = T::from_i8(15).unwrap();

    // sign(P) and sign(Q), each via the exact biquadratic
    // `signum_sum_sqrt_expr_4(_, 1, _, 3, _, 5, _, 15)`.
    let sp = signum_sum_sqrt_expr_4(a, int1, b, int3, c, int5, e, int15);
    let sq = signum_sum_sqrt_expr_4(d, int1, f, int3, g, int5, h, int15);

    if sp == sq {
        return sp;
    }
    if sq.is_zero() {
        return sp;
    }
    if sp.is_zero() {
        return sq;
    }

    // P^2 in basis {1, sqrt(3), sqrt(5), sqrt(15)}:
    //   const:    a^2 + 3*b^2 + 5*c^2 + 15*e^2
    //   sqrt(3):  2*a*b + 10*c*e        (sqrt(5)*sqrt(15) = 5*sqrt(3))
    //   sqrt(5):  2*a*c + 6*b*e         (sqrt(3)*sqrt(15) = 3*sqrt(5))
    //   sqrt(15): 2*a*e + 2*b*c         (sqrt(3)*sqrt(5) = sqrt(15))
    let p2_0 = a * a + int3 * b * b + int5 * c * c + int15 * e * e;
    let p2_1 = int2 * a * b + int10 * c * e;
    let p2_2 = int2 * a * c + int6 * b * e;
    let p2_3 = int2 * a * e + int2 * b * c;

    // Q^2 in the same basis:
    let q2_0 = d * d + int3 * f * f + int5 * g * g + int15 * h * h;
    let q2_1 = int2 * d * f + int10 * g * h;
    let q2_2 = int2 * d * g + int6 * f * h;
    let q2_3 = int2 * d * h + int2 * f * g;

    // Q^2 * (10 - 2*sqrt(5)) in {1, sqrt(3), sqrt(5), sqrt(15)}:
    //   const:    10*q2_0 - 10*q2_2     (2*sqrt(5)*sqrt(5) = 10)
    //   sqrt(3):  10*q2_1 - 10*q2_3
    //   sqrt(5):  10*q2_2 - 2*q2_0
    //   sqrt(15): 10*q2_3 - 2*q2_1
    let qy_0 = int10 * q2_0 - int10 * q2_2;
    let qy_1 = int10 * q2_1 - int10 * q2_3;
    let qy_2 = int10 * q2_2 - int2 * q2_0;
    let qy_3 = int10 * q2_3 - int2 * q2_1;

    let alpha = p2_0 - qy_0;
    let beta = p2_1 - qy_1;
    let gamma = p2_2 - qy_2;
    let delta = p2_3 - qy_3;

    let spq = signum_sum_sqrt_expr_4(alpha, int1, beta, int3, gamma, int5, delta, int15);

    -sq * spq
}

// ----------------
// Tests

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sum_root_expr_sign_2() {
        assert_eq!(signum_sum_sqrt_expr_2(0, 2, 0, 3), 0);
        assert_eq!(signum_sum_sqrt_expr_2(1, 2, 0, 3), 1);
        assert_eq!(signum_sum_sqrt_expr_2(0, 2, -1, 3), -1);
        assert_eq!(signum_sum_sqrt_expr_2(2, 2, -1, 3), 1);
        assert_eq!(signum_sum_sqrt_expr_2(-5, 2, 4, 3), -1);
        assert_eq!(signum_sum_sqrt_expr_2(-5, 2, 5, 3), 1);
    }

    #[test]
    fn test_sum_root_expr_sign_4() {
        let sign_zz24 = |a, b, c, d| signum_sum_sqrt_expr_4(a, 1, b, 2, c, 3, d, 6);

        // trivial sanity-checks
        assert_eq!(sign_zz24(0, 0, 0, 0), 0);
        assert_eq!(sign_zz24(1, 1, 1, 1), 1);
        assert_eq!(sign_zz24(-1, -1, -1, -1), -1);
        assert_eq!(sign_zz24(1, 0, 0, 0), 1);
        assert_eq!(sign_zz24(0, -1, 0, 0), -1);
        assert_eq!(sign_zz24(0, 0, 1, 0), 1);
        assert_eq!(sign_zz24(0, 0, 0, -1), -1);
        // non-trivial tests
        assert_eq!(sign_zz24(5, 7, 11, -13), 1);
        assert_eq!(sign_zz24(5, 7, 11, -14), -1);
        assert_eq!(sign_zz24(17, -11, 9, -7), -1);
        assert_eq!(sign_zz24(18, -11, 9, -7), 1);
        assert_eq!(sign_zz24(18, -11, 8, -7), -1);
        assert_eq!(sign_zz24(18, -11, 8, -6), 1);

        // try with parameters where terms are all really close
        {
            let (a, b, c, d) = (130, 92, 75, 53);
            assert_eq!(sign_zz24(-a, -b, c, d), -1);
            assert_eq!(sign_zz24(-a, b, -c, d), 1);
            assert_eq!(sign_zz24(-a, b, c, -d), 1);
            assert_eq!(sign_zz24(a, -b, -c, d), -1);
            assert_eq!(sign_zz24(a, -b, c, -d), -1);
            assert_eq!(sign_zz24(a, b, -c, -d), 1);
        }
        {
            let (a, b, c, d) = (485, 343, 280, 198);
            assert_eq!(sign_zz24(-a, -b, c, d), -1);
            assert_eq!(sign_zz24(-a, b, -c, d), 1);
            assert_eq!(sign_zz24(-a, b, c, -d), 1);
            assert_eq!(sign_zz24(a, -b, -c, d), -1);
            assert_eq!(sign_zz24(a, -b, c, -d), -1);
            assert_eq!(sign_zz24(a, b, -c, -d), 1);
        }
    }

    #[test]
    fn test_sum_root_expr_sign_4_pentagonal() {
        let sign_zz10 = signum_sum_sqrt_expr_4_pentagonal::<i64>;

        assert_eq!(sign_zz10(0, 0, 0, 0), 0);
        assert_eq!(sign_zz10(1, 0, 0, 0), 1);
        assert_eq!(sign_zz10(-1, 0, 0, 0), -1);
        assert_eq!(sign_zz10(0, 1, 0, 0), 1);
        assert_eq!(sign_zz10(0, -1, 0, 0), -1);
        assert_eq!(sign_zz10(0, 0, 1, 0), 1);
        assert_eq!(sign_zz10(0, 0, -1, 0), -1);
        assert_eq!(sign_zz10(0, 0, 0, 1), 1);
        assert_eq!(sign_zz10(0, 0, 0, -1), -1);

        assert_eq!(sign_zz10(1, 1, 0, 0), 1);
        assert_eq!(sign_zz10(-1, -1, 0, 0), -1);
        assert_eq!(sign_zz10(0, 0, 1, 1), 1);
        assert_eq!(sign_zz10(0, 0, -1, -1), -1);

        assert_eq!(sign_zz10(-3, 0, 1, 0), -1);
        assert_eq!(sign_zz10(3, 0, -1, 0), 1);
        assert_eq!(sign_zz10(1, 0, 0, -2), -1);
        assert_eq!(sign_zz10(-1, 0, 0, 2), 1);
    }

    #[test]
    fn test_sum_root_expr_sign_4_zz16() {
        // sign of a + b*sqrt(2) + c*sqrt(2+sqrt(2)) + d*sqrt(2*(2+sqrt(2)))
        let sign_zz16 = signum_sum_sqrt_expr_4_zz16::<i64>;

        // Trivial axes.
        assert_eq!(sign_zz16(0, 0, 0, 0), 0);
        assert_eq!(sign_zz16(1, 0, 0, 0), 1);
        assert_eq!(sign_zz16(-1, 0, 0, 0), -1);
        assert_eq!(sign_zz16(0, 1, 0, 0), 1);
        assert_eq!(sign_zz16(0, -1, 0, 0), -1);
        assert_eq!(sign_zz16(0, 0, 1, 0), 1);
        assert_eq!(sign_zz16(0, 0, -1, 0), -1);
        assert_eq!(sign_zz16(0, 0, 0, 1), 1);
        assert_eq!(sign_zz16(0, 0, 0, -1), -1);

        // Same-sign sums.
        assert_eq!(sign_zz16(1, 1, 0, 0), 1);
        assert_eq!(sign_zz16(-1, -1, 0, 0), -1);
        assert_eq!(sign_zz16(0, 0, 1, 1), 1);
        assert_eq!(sign_zz16(0, 0, -1, -1), -1);

        // Mixed-sign with a definite winner (numerical sanity, exact via
        // closed-form): with sqrt(2) ~ 1.414, sqrt(2+sqrt(2)) ~ 1.848,
        // sqrt(2(2+sqrt(2))) ~ 2.613.
        //   1 + (-1)*sqrt(2) = ~ -0.414       -> -1
        assert_eq!(sign_zz16(1, -1, 0, 0), -1);
        //   2 + (-1)*sqrt(2) = ~ 0.586        -> 1
        assert_eq!(sign_zz16(2, -1, 0, 0), 1);
        //   1 + (-1)*sqrt(2+sqrt(2)) ~ -0.848 -> -1
        assert_eq!(sign_zz16(1, 0, -1, 0), -1);
        //   2 + (-1)*sqrt(2+sqrt(2)) ~ 0.152  -> 1
        assert_eq!(sign_zz16(2, 0, -1, 0), 1);
        //   sqrt(2(2+sqrt(2))) ~ 2.613, sqrt(2) ~ 1.414, so 2*sqrt(2) ~ 2.828 > 2.613
        //   2*sqrt(2) - sqrt(2(2+sqrt(2))) ~ 0.215 > 0
        assert_eq!(sign_zz16(0, 2, 0, -1), 1);
        //   sqrt(2) - sqrt(2(2+sqrt(2))) ~ -1.199 -> -1
        assert_eq!(sign_zz16(0, 1, 0, -1), -1);

        // P ~ 0 cases (sign of Q wins).
        //   sqrt(2+sqrt(2)) + sqrt(2(2+sqrt(2))) = sqrt(2+sqrt(2)) * (1 + sqrt(2)) > 0
        assert_eq!(sign_zz16(0, 0, 1, 1), 1);
        // Q ~ 0 cases (sign of P wins).
        //   1 + sqrt(2) > 0
        assert_eq!(sign_zz16(1, 1, 0, 0), 1);

        // P > 0, Q > 0 -> 1; P < 0, Q < 0 -> -1.
        assert_eq!(sign_zz16(1, 1, 1, 1), 1);
        assert_eq!(sign_zz16(-1, -1, -1, -1), -1);
    }

    #[test]
    fn test_sum_root_expr_sign_8_zz60() {
        // sign of a + b*sqrt(3) + c*sqrt(5) + d*sqrt(10-2*sqrt(5))
        //        + e*sqrt(15) + f*sqrt(3(10-2*sqrt(5)))
        //        + g*sqrt(5(10-2*sqrt(5))) + h*sqrt(15(10-2*sqrt(5)))
        let s = signum_sum_sqrt_expr_8_zz60::<i64>;

        // Trivial axes (each basis element is positive).
        assert_eq!(s(0, 0, 0, 0, 0, 0, 0, 0), 0);
        assert_eq!(s(1, 0, 0, 0, 0, 0, 0, 0), 1);
        assert_eq!(s(-1, 0, 0, 0, 0, 0, 0, 0), -1);
        assert_eq!(s(0, 1, 0, 0, 0, 0, 0, 0), 1);
        assert_eq!(s(0, -1, 0, 0, 0, 0, 0, 0), -1);
        assert_eq!(s(0, 0, 1, 0, 0, 0, 0, 0), 1);
        assert_eq!(s(0, 0, -1, 0, 0, 0, 0, 0), -1);
        assert_eq!(s(0, 0, 0, 1, 0, 0, 0, 0), 1);
        assert_eq!(s(0, 0, 0, -1, 0, 0, 0, 0), -1);
        assert_eq!(s(0, 0, 0, 0, 1, 0, 0, 0), 1);
        assert_eq!(s(0, 0, 0, 0, 0, 1, 0, 0), 1);
        assert_eq!(s(0, 0, 0, 0, 0, 0, 1, 0), 1);
        assert_eq!(s(0, 0, 0, 0, 0, 0, 0, 1), 1);

        // All-1 / all-(-1) are unambiguous.
        assert_eq!(s(1, 1, 1, 1, 1, 1, 1, 1), 1);
        assert_eq!(s(-1, -1, -1, -1, -1, -1, -1, -1), -1);

        // Mixed-sign sanity (numerical): sqrt(3) ~ 1.732, sqrt(5) ~ 2.236,
        //   sqrt(15) ~ 3.873, sqrt(10-2*sqrt(5)) ~ 2.351,
        //   sqrt(3*(10-2*sqrt(5))) ~ 4.072, sqrt(5*(10-2*sqrt(5))) ~ 5.257,
        //   sqrt(15*(10-2*sqrt(5))) ~ 9.106.
        //
        //   1 + (-1)*sqrt(3) ~ -0.732 -> -1
        assert_eq!(s(1, -1, 0, 0, 0, 0, 0, 0), -1);
        //   2 + (-1)*sqrt(3) ~ 0.268 -> 1
        assert_eq!(s(2, -1, 0, 0, 0, 0, 0, 0), 1);
        //   sqrt(15) - sqrt(5) - sqrt(3) ~ -0.095 -> -1
        assert_eq!(s(0, -1, -1, 0, 1, 0, 0, 0), -1);
        //   sqrt(15) - sqrt(5) ~ 1.637 -> 1
        assert_eq!(s(0, 0, -1, 0, 1, 0, 0, 0), 1);

        // P near zero, Q strictly positive: result should be sign(Q).
        //   d*sqrt(y) -- d > 0 -> sign = +1.
        assert_eq!(s(0, 0, 0, 1, 0, 0, 0, 0), 1);
        assert_eq!(s(0, 0, 0, -1, 0, 0, 0, 0), -1);

        // Mixed P, Q opposite signs: nontrivial reduction kicks in.
        //   P = -1 (so P^2 = 1), Q = 1 (so Q^2 = 1), y ~ 5.528.
        //   z = -1 + sqrt(y) ~ -1 + 2.351 = 1.351 > 0 -> 1.
        assert_eq!(s(-1, 0, 0, 1, 0, 0, 0, 0), 1);
        //   P = -3, Q = 1, y ~ 5.528: z = -3 + 2.351 < 0 -> -1.
        assert_eq!(s(-3, 0, 0, 1, 0, 0, 0, 0), -1);
        //   P = 3, Q = -1: z = 3 - 2.351 > 0 -> 1.
        assert_eq!(s(3, 0, 0, -1, 0, 0, 0, 0), 1);
    }

    #[test]
    fn test_sum_root_expr_sign_8_zz32() {
        // sign of a0 + a1*sqrt(2) + a2*sqrt(2+sqrt(2)) + a3*sqrt(2+sqrt(2+sqrt(2)))
        //        + a4*sqrt(2(2+sqrt(2))) + a5*sqrt(2(2+sqrt(2+sqrt(2))))
        //        + a6*sqrt((2+sqrt(2))(2+sqrt(2+sqrt(2))))
        //        + a7*sqrt(2(2+sqrt(2))(2+sqrt(2+sqrt(2))))
        let s = signum_sum_sqrt_expr_8_zz32::<i64>;

        // Numerical reference values:
        //   sqrt(2)                                   ~ 1.4142
        //   sqrt(2+sqrt(2))                           ~ 1.8478
        //   sqrt(2+sqrt(2+sqrt(2)))                   ~ 1.9616
        //   sqrt(2(2+sqrt(2)))                        ~ 2.6131
        //   sqrt(2(2+sqrt(2+sqrt(2))))                ~ 2.7741
        //   sqrt((2+sqrt(2))(2+sqrt(2+sqrt(2))))      ~ 3.6249
        //   sqrt(2(2+sqrt(2))(2+sqrt(2+sqrt(2))))     ~ 5.1262

        // Trivial axes (each basis element is positive).
        assert_eq!(s(0, 0, 0, 0, 0, 0, 0, 0), 0);
        for i in 0..8 {
            let mut v = [0i64; 8];
            v[i] = 1;
            assert_eq!(
                s(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7]),
                1,
                "axis {i} positive should be +1"
            );
            v[i] = -1;
            assert_eq!(
                s(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7]),
                -1,
                "axis {i} negative should be -1"
            );
        }

        // All-1 / all-(-1).
        assert_eq!(s(1, 1, 1, 1, 1, 1, 1, 1), 1);
        assert_eq!(s(-1, -1, -1, -1, -1, -1, -1, -1), -1);

        // Mixed-sign with definite winners.
        //   1 + (-1)*sqrt(2) ~ -0.414       -> -1
        assert_eq!(s(1, -1, 0, 0, 0, 0, 0, 0), -1);
        //   2 + (-1)*sqrt(2) ~ 0.586        -> 1
        assert_eq!(s(2, -1, 0, 0, 0, 0, 0, 0), 1);
        //   sqrt(2(2+sqrt(2))) - sqrt(2+sqrt(2+sqrt(2))) ~ 2.613 - 1.962 ~ 0.651 -> 1
        assert_eq!(s(0, 0, 0, -1, 1, 0, 0, 0), 1);
        //   sqrt(2+sqrt(2+sqrt(2))) - sqrt(2(2+sqrt(2))) ~ -0.651 -> -1
        assert_eq!(s(0, 0, 0, 1, -1, 0, 0, 0), -1);

        // Recursion-triggering: P and Q nonzero with opposite signs.
        //   P = -1, Q = 1 -> z = -1 + sqrt(y) where y = 2+sqrt(2+sqrt(2)) ~ 3.848
        //   sqrt(y) ~ 1.962, so z ~ 0.962 > 0 -> 1.
        assert_eq!(s(-1, 0, 0, 1, 0, 0, 0, 0), 1);
        //   P = -2, Q = 1: z ~ -0.038 < 0 -> -1.
        assert_eq!(s(-2, 0, 0, 1, 0, 0, 0, 0), -1);
        //   P = -3, Q = 1: z ~ -1.038 -> -1.
        assert_eq!(s(-3, 0, 0, 1, 0, 0, 0, 0), -1);
        //   P = 3, Q = -1: z = 3 - 1.962 > 0 -> 1.
        assert_eq!(s(3, 0, 0, -1, 0, 0, 0, 0), 1);
    }
}
