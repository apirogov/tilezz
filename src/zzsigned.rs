use num_traits::{FromPrimitive, Signed, Zero};

use crate::traits::IntRing;
use crate::zzbase::ZNum;

/// ZSigned is like num_traits::Signed, but without the Num-like constraints.
pub trait ZSigned: IntRing {
    fn signum(&self) -> Self;

    fn abs(&self) -> Self {
        *self * self.signum()
    }
    fn is_positive(&self) -> bool {
        self.signum() == Self::one()
    }
    fn is_negative(&self) -> bool {
        self.signum() == -Self::one()
    }
    fn abs_sub(&self, other: &Self) -> Self {
        <Self as ZSigned>::abs(&(*self - *other))
    }
}

impl<T: IntRing + Signed> ZSigned for T {
    fn signum(&self) -> Self {
        self.signum()
    }
}




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

// --------

/// This implementation uses f64 (works in all ZZ rings).
pub fn zz_partial_signum_fallback<Z: ZNum>(val: &Z) -> Z {
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
pub fn zz_partial_signum_1_sym<Z: ZNum>(val: &Z) -> Z {
    let cs = val.zz_coeffs();

    let mut result = Z::zero();
    result.zz_coeffs_mut()[0] = Z::Scalar::from(cs[0].signum());
    result
}

pub fn zz_partial_signum_2_sym<Z: ZNum>(val: &Z) -> Z {
    let ftoi = |f| Z::Scalar::from_f64(f).unwrap();
    let cs = val.zz_coeffs();
    let rs = Z::zz_params().sym_roots_sqs;

    let (a, b) = (cs[0], cs[1]);
    let (m, n) = (ftoi(rs[0]), ftoi(rs[1]));
    let sgn = signum_sum_sqrt_expr_2(a, m, b, n);

    let mut result = Z::zero();
    result.zz_coeffs_mut()[0] = Z::Scalar::from(sgn);
    result
}

pub fn zz_partial_signum_4_sym<Z: ZNum>(val: &Z) -> Z {
    let cs: Vec<Z::Scalar> = val.zz_coeffs().to_vec();
    let rs: Vec<Z::Scalar> = Z::zz_params()
        .sym_roots_sqs
        .iter()
        .map(|r| Z::Scalar::from_f64(*r).unwrap())
        .collect();

    let (a, b, c, d) = (cs[0], cs[1], cs[2], cs[3]);
    let (k, m, n, l) = (rs[0], rs[1], rs[2], rs[3]);

    let sgn = signum_sum_sqrt_expr_4(a, k, b, m, c, n, d, l);

    let mut result = Z::zero();
    result.zz_coeffs_mut()[0] = Z::Scalar::from(sgn);
    result
}

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
}
