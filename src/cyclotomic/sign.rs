use num_traits::{FromPrimitive, Zero};

use super::numtraits::{IntRing, ZSigned};
use super::traits::IsReal;

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

    // println!("ad: {sgn_ad_terms:?} bc: {sgn_bc_terms:?}");
    // println!("{a:?} {k:?} {b:?} {m:?} {c:?} {n:?} {d:?} {l:?}");

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

    // println!("{sgn_bc_terms:?} | {lhs:?} | {sq_lhs:?} || {ad_m_bc:?} | {sq_rhs:?}");

    sgn_bc_terms.signum() * (sq_lhs - sq_rhs).signum()
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

// --------

/// This implementation uses f64 (works in all ZZ rings).
pub fn zz_partial_signum_fallback<Z: IsReal>(val: &Z) -> Z {
    let val = val.complex64().re;
    if val.is_zero() {
        Z::zero()
    } else {
        let mut result = Z::zero();
        result.zz_coeffs_mut()[0] = Z::Scalar::from_f64(val.signum()).unwrap();
        result
    }
}

/// Trivial if there is just one term in the expression,
/// regardless of its symbolic root (which usually will be 1).
pub fn zz_partial_signum_1_sym<Z: IsReal>(val: &Z) -> Z {
    let cs = val.zz_coeffs();

    let mut result = Z::zero();
    result.zz_coeffs_mut()[0] = cs[0].signum();
    result
}

pub fn zz_partial_signum_2_sym<Z: IsReal>(val: &Z) -> Z {
    let ftoi = |f| Z::Scalar::from_f64(f).unwrap();
    let cs = val.zz_coeffs();
    let rs = Z::zz_params().sym_roots_sqs;

    let (a, b) = (cs[0], cs[1]);
    let (m, n) = (ftoi(rs[0]), ftoi(rs[1]));
    let sgn = signum_sum_sqrt_expr_2(a, m, b, n);

    let mut result = Z::zero();
    result.zz_coeffs_mut()[0] = sgn;
    result
}

pub fn zz_partial_signum_4_sym<Z: IsReal>(val: &Z) -> Z {
    // Avoid heap allocations: we know this is only used for rings with exactly 4 symbolic roots.
    let cs = val.zz_coeffs();
    debug_assert!(cs.len() == 4);

    // NOTE: sym_roots_sqs are f64 parameters, but they are compile-time constants.
    // Convert them directly without allocating a Vec.
    let rs = Z::zz_params().sym_roots_sqs;
    debug_assert!(rs.len() == 4);

    let (a, b, c, d) = (cs[0], cs[1], cs[2], cs[3]);
    let (k, m, n, l) = (
        Z::Scalar::from_f64(rs[0]).unwrap(),
        Z::Scalar::from_f64(rs[1]).unwrap(),
        Z::Scalar::from_f64(rs[2]).unwrap(),
        Z::Scalar::from_f64(rs[3]).unwrap(),
    );

    let sgn = signum_sum_sqrt_expr_4(a, k, b, m, c, n, d, l);

    let mut result = Z::zero();
    result.zz_coeffs_mut()[0] = sgn;
    result
}

pub fn zz_partial_signum_4_pentagonal<Z: IsReal>(val: &Z) -> Z {
    let cs = val.zz_coeffs();
    debug_assert!(cs.len() == 4);

    let rs = Z::zz_params().sym_roots_sqs;
    debug_assert!(rs.len() == 4);

    let sgn = signum_sum_sqrt_expr_4_pentagonal(cs[0], cs[1], cs[2], cs[3]);

    let mut result = Z::zero();
    result.zz_coeffs_mut()[0] = sgn;
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    use super::super::constants::*;
    use super::super::symnum::ZZComplex;
    use super::super::types::{Z10, Z24, ZZ10, ZZ24};
    use super::super::units::Units;
    use num_traits::One;

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
    fn test_signum() {
        let sq2: ZZ24 = sqrt2();
        let sq3: ZZ24 = sqrt3();
        let sq6: ZZ24 = sqrt6();

        let z = Z24::zero();
        let p = Z24::one();
        let m = -p;

        // use same test as above
        let sign_zz24 = |a, b, c, d| {
            (ZZ24::from(a) + ZZ24::from(b) * sq2 + ZZ24::from(c) * sq3 + ZZ24::from(d) * sq6)
                .re()
                .signum()
        };

        let (a, b, c, d) = (485, 343, 280, 198);
        assert_eq!(sign_zz24(0, 0, 0, 0), z);
        assert_eq!(sign_zz24(-a, -b, c, d), m);
        assert_eq!(sign_zz24(-a, b, -c, d), p);
        assert_eq!(sign_zz24(-a, b, c, -d), p);
        assert_eq!(sign_zz24(a, -b, -c, d), m);
        assert_eq!(sign_zz24(a, -b, c, -d), m);
        assert_eq!(sign_zz24(a, b, -c, -d), p);
    }

    #[test]
    fn test_sum_root_expr_sign_4_pentagonal() {
        let sign_zz10 = |a, b, c, d| signum_sum_sqrt_expr_4_pentagonal(a, b, c, d);

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
    fn test_signum_zz10() {
        let sq5: ZZ10 = sqrt5();
        let sqy: ZZ10 = <ZZ10 as Units>::unit(1).im().into();
        let sq5y: ZZ10 = sq5 * sqy;

        let z = Z10::zero();
        let p = Z10::one();
        let m = -p;

        let sign_zz10 = |a, b, c, d| {
            (ZZ10::from(a) + ZZ10::from(b) * sq5 + ZZ10::from(c) * sqy + ZZ10::from(d) * sq5y)
                .re()
                .signum()
        };

        let (a, b, c, d) = (485, 343, 280, 198);
        assert_eq!(sign_zz10(0, 0, 0, 0), z);
        assert_eq!(sign_zz10(-a, -b, c, d), m);
        assert_eq!(sign_zz10(-a, b, -c, d), p);
        assert_eq!(sign_zz10(-a, b, c, -d), p);
        assert_eq!(sign_zz10(a, -b, -c, d), m);
        assert_eq!(sign_zz10(a, -b, c, -d), m);
        assert_eq!(sign_zz10(a, b, -c, -d), p);
    }
}
