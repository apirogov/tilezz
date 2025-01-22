use super::gaussint::GaussInt;
use super::traits::{Ccw, CycIntRing, InnerIntType, IntField, IntRing, RealSigned};
use num_complex::Complex64;
use num_integer::Integer;
use num_rational::Ratio;
use num_traits::{FromPrimitive, One, PrimInt, Zero};
use std::fmt::Debug;
use std::fmt::Display;

// configure a signed primitive integer size here
pub type MyInt = i64;

// internally, ZZn use a Gaussian integer where all coefficients are
// a ratio with a fixed denominator (modulo simplification)
pub type Frac = Ratio<MyInt>;
pub type GInt = GaussInt<Frac>;

#[derive(Debug)]
pub struct ZZParams<'a> {
    // pub phantom: PhantomData<&'a T>,
    /// Squares of symbolic roots. IMPORTANT: we assume that the first one is 1
    pub sym_roots_sqs: &'a [f64],
    /// Labels of symbolic roots.
    pub sym_roots_lbls: &'a [&'a str],
    /// Number of symbolic roots in array sym_roots_sqs.
    pub sym_roots_num: usize,
    /// Number of steps in this complex integer ring that makes a full rotation.
    pub full_turn_steps: i8,
    /// Scaling factor 1/l common to all terms in the symbolic root number sum.
    pub scaling_fac: MyInt,
    /// Unit of rotation in coefficients of this ZZ type
    pub ccw_unit_coeffs: &'a [[MyInt; 2]],
}

impl ZZParams<'static> {
    // helper function to lift the param coeffs into a GaussInt of suitable type
    pub fn ccw_unit<T: PrimInt + Integer + IntRing>(&self) -> Vec<GaussInt<Ratio<T>>> {
        let sc = T::from(self.scaling_fac).unwrap();
        self.ccw_unit_coeffs
            .into_iter()
            .map(|x| {
                GaussInt::new(
                    Ratio::<T>::new_raw(T::from(x[0]).unwrap(), sc),
                    Ratio::<T>::new_raw(T::from(x[1]).unwrap(), sc),
                )
            })
            .collect()
    }
}

// --------------------------------

/// Floating-point-free solution to get sign of an expression
/// a*sqrt(n) + b*sqrt(m)
/// where a,b,m,n are all integers
/// and m,n are squarefree coprime constants > 1.
pub fn signum_sum_sqrt_expr_2<T: IntRing + RealSigned>(a: T, m: T, b: T, n: T) -> T {
    // if a and b are both positive or negative -> trivial,
    let sgn_a = a.re_signum();
    let sgn_b = b.re_signum();
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
    (sgn_a * a * a * m + sgn_b * b * b * n).re_signum()
}

/// Floating-point-free solution to get sign of an expression
/// a + b*sqrt(m) + c*sqrt(n) + d*sqrt(m*n)
/// where a,b,c,m,n are all integers
/// and m,n are squarefree coprime constants > 1.
pub fn signum_sum_sqrt_expr_4<T: IntRing + RealSigned + FromPrimitive>(
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
    let sq_lhs = lhs.re_signum() * lhs * lhs;
    let ad_m_bc = (a * d) - (b * c);
    let sq_rhs = four * mn * ad_m_bc.re_signum() * ad_m_bc * ad_m_bc;

    // println!("{sgn_bc_terms:?} | {lhs:?} | {sq_lhs:?} || {ad_m_bc:?} | {sq_rhs:?}");

    sgn_bc_terms.re_signum() * (sq_lhs - sq_rhs).re_signum()
}

// --------------------------------

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
    result.zz_coeffs_mut()[0] = Z::Scalar::from(cs[0].re_signum());
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

// --------------------------------

/// Compute inverse by repeated rationalizing of the denominator.
///
/// IMPORTANT: This only works correctly if there are no nested roots in
/// the representation of the cyclotomic field, i.e. it does NOT work
/// correctly for fields that contain ZZ5 or ZZ8 as subfields.
#[inline]
pub fn zz_inv<Z: ZNum>(val: &Z) -> Z {
    // for x/y where y = a + b, b being a single (scaled) square root,
    // we compute y' = a - b and produce x*y'/(a+b)(a-b) = x*y'/a^2-b^2
    // where a^2 is a simpler term and b is rational.
    // we repeat this for all square roots in the denominator.

    let num_terms = Z::zz_params().sym_roots_num;

    let mut numer = Z::one();
    let mut denom = val.clone();

    let mut root_ix = 0;
    let mut non_root_ix = 0;
    let mut root_found = true;
    while root_found {
        root_found = false;
        for i in 0..num_terms {
            if Z::zz_params().sym_roots_sqs[i].is_one() {
                non_root_ix = i; // we need the non-irrational part later
                continue; // non-irrational term
            }
            let c = denom.zz_coeffs()[i];
            if !c.is_zero() {
                root_found = true;
                root_ix = i;
                break;
            }
        }

        if !root_found {
            break; // rational denominator
        }

        // "conjugate" one square root to use the binomial trick
        let mut conjugated = denom.clone();
        let curr_root_coeff = conjugated.zz_coeffs()[root_ix];
        conjugated.zz_coeffs_mut()[root_ix] = -curr_root_coeff;

        // compute b^2
        let mut curr_root_term_sq = Z::zero();
        curr_root_term_sq.zz_coeffs_mut()[root_ix] = curr_root_coeff; // = b
        curr_root_term_sq = curr_root_term_sq * curr_root_term_sq; // = b^2

        // update numerator (= x * y')
        numer = numer * conjugated;

        // update denominator (= (a + b)(a - b)= a^2 - b^2)
        denom.zz_coeffs_mut()[root_ix] = Z::Scalar::zero(); // = a
        denom = denom * denom - curr_root_term_sq; // = a^2 - b^2
    }

    // now we have a rational denominator (i.e. no square root terms)
    // so we can just flip it and multiply with the numerator.
    let mut inv_denom = Z::zero();
    inv_denom.zz_coeffs_mut()[non_root_ix] = Z::Scalar::one() / denom.zz_coeffs()[non_root_ix];

    return numer * inv_denom;
}

// --------------------------------

pub trait ZZBase {
    type Scalar: PartialOrd + IntField + InnerIntType + FromPrimitive + Debug + RealSigned + 'static;
    type Real: ZNum;

    /// Return angle representing one full turn.
    #[inline]
    fn turn() -> i8 {
        Self::zz_params().full_turn_steps
    }

    /// Return angle representing half a turn.
    #[inline]
    fn hturn() -> i8 {
        Self::turn() / 2
    }

    /// Return whether this ring supports a quarter turn,
    /// i.e. can represent the imaginary unit i.
    #[inline]
    fn has_qturn() -> bool {
        Self::turn() % 4 == 0
    }

    /// Return angle representing a quarter turn (if ring supports it).
    #[inline]
    fn qturn() -> i8 {
        assert_eq!(Self::turn() % 4, 0);
        Self::turn() / 4
    }

    /// Return unit length vector pointing in direction of given angle.
    fn unit(i: i8) -> Self
    where
        Self: ZZNum + One + Ccw,
    {
        let one = Self::one();
        let ccw = Self::ccw();
        let j = i.rem_euclid(<Self as ZZBase>::turn());
        one * <Self as ZZBase>::pow(&ccw, j)
    }

    /// Raise to an integer power.
    // NOTE: using i8 instead of u8 for convenience (angles use i8)
    fn pow(&self, i: i8) -> Self
    where
        Self: ZNum,
    {
        assert!(i >= 0, "Negative powers are not supported!");
        let mut x = Self::one();
        for _ in 0..i {
            x = x * (*self);
        }
        return x;
    }

    /// Return imaginary unit (if the ring supports it).
    #[inline]
    fn one_i() -> Self
    where
        Self: ZZNum,
    {
        let qt = <Self as ZZBase>::qturn();
        <Self as ZZBase>::unit(qt)
    }

    /// Scalar multiplication.
    #[inline]
    fn scale(&self, scalar: i64) -> Self
    where
        Self: Sized,
    {
        let cs: Vec<Self::Scalar> = Self::zz_mul_scalar(self.zz_coeffs(), scalar);
        Self::new(&cs)
    }

    /// Convert to a complex floating point number.
    fn complex64(&self) -> Complex64;

    // functions that can be implemented via zz_base_impl!
    // --------
    fn new(coeffs: &[Self::Scalar]) -> Self;

    fn conj(&self) -> Self;

    fn zz_coeffs(&self) -> &[Self::Scalar];
    fn zz_coeffs_mut(&mut self) -> &mut [Self::Scalar];

    fn zz_params() -> &'static ZZParams<'static>;
    // fn zz_ccw_vec() -> Vec<Self::Scalar>;
    fn zz_mul_arrays(x: &[Self::Scalar], y: &[Self::Scalar]) -> Vec<Self::Scalar>;
    fn zz_mul_scalar(arr: &[Self::Scalar], scalar: i64) -> Vec<Self::Scalar>;

    // implementations for implementing other traits using zz_ops_impl!
    // --------
    #[inline]
    fn zz_zero_vec() -> Vec<Self::Scalar> {
        vec![Self::Scalar::zero(); Self::zz_params().sym_roots_num]
    }

    fn zz_one_vec() -> Vec<Self::Scalar> {
        let mut ret = vec![Self::Scalar::zero(); Self::zz_params().sym_roots_num];
        ret[0] = Self::Scalar::one();
        ret
    }

    fn zz_add(&self, other: &Self) -> Vec<Self::Scalar> {
        let mut ret = Self::zz_zero_vec();
        for (i, (aval, bval)) in self.zz_coeffs().iter().zip(other.zz_coeffs()).enumerate() {
            ret[i] = *aval + *bval;
        }
        ret
    }

    fn zz_sub(&self, other: &Self) -> Vec<Self::Scalar> {
        let mut ret = Self::zz_zero_vec();
        for (i, (aval, bval)) in self.zz_coeffs().iter().zip(other.zz_coeffs()).enumerate() {
            ret[i] = *aval - *bval;
        }
        ret
    }

    fn zz_neg(&self) -> Vec<Self::Scalar> {
        let mut ret = Self::zz_zero_vec();
        for (i, val) in self.zz_coeffs().iter().enumerate() {
            ret[i] = -(*val);
        }
        ret
    }

    #[inline]
    fn zz_mul(&self, other: &Self) -> Vec<Self::Scalar> {
        Self::zz_mul_arrays(self.zz_coeffs(), other.zz_coeffs())
    }
}

pub trait ZZComplex {
    /// Return true if the value is purely real.
    fn is_real(&self) -> bool;

    /// Return true if the value is purely imaginary.
    fn is_imag(&self) -> bool;

    /// Return true if the value is mixed real and imaginary.
    fn is_complex(&self) -> bool {
        !self.is_real() && !self.is_imag()
    }

    /// Return the real part of the value,
    /// i.e. the value (z + z.conj()) / 2
    fn re(&self) -> <Self as ZZBase>::Real
    where
        Self: ZZBase;

    /// Return the imaginary part of the value (rotated onto the real axis),
    /// i.e. the value (z - z.conj()) / 2i
    fn im(&self) -> <Self as ZZBase>::Real
    where
        Self: ZZBase;

    /// Split the value into its real and imaginary contributions.
    /// Note that the imaginary component is converted to real.
    ///
    // NOTE: z = dot(1, z) + i*wedge(1, z), as dot(1, z) = Re(z), wedge(1, z) = Im(z)
    fn re_im(&self) -> (<Self as ZZBase>::Real, <Self as ZZBase>::Real)
    where
        Self: ZZBase,
    {
        (self.re(), self.im())
    }
}

/// Quadratic real extension corresponding to a cyclotomic ring
/// (used to e.g. split a cyclotomic value into separate real and imaginary parts)
pub trait ZNum: ZZBase + InnerIntType + RealSigned + IntRing + Display {}

/// A cyclotomic ring. You probably want to parametrize generic code over this trait.
pub trait ZZNum: ZNum + ZZComplex + CycIntRing {}

// --------------------------------

#[macro_export]
macro_rules! zz_base_impl {
    ($name:ident, $name_real:ident, $params:ident, $mul_func:ident, $re_signum_func:ident) => {
        #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
        pub struct $name_real {
            coeffs: [Frac; $params.sym_roots_num],
        }
        #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
        pub struct $name {
            coeffs: [GInt; $params.sym_roots_num],
        }

        impl InnerIntType for $name_real {
            type IntType = <Frac as InnerIntType>::IntType;
        }
        impl InnerIntType for $name {
            type IntType = <GInt as InnerIntType>::IntType;
        }

        impl ZZBase for $name_real {
            type Scalar = Frac;
            type Real = $name_real;

            fn conj(&self) -> Self {
                self.clone()
            }

            #[inline]
            fn zz_coeffs(&self) -> &[Self::Scalar] {
                &self.coeffs
            }

            #[inline]
            fn zz_coeffs_mut(&mut self) -> &mut [Self::Scalar] {
                &mut self.coeffs
            }

            #[inline]
            fn zz_params() -> &'static ZZParams<'static> {
                &$params
            }

            #[inline]
            fn zz_mul_arrays(x: &[Self::Scalar], y: &[Self::Scalar]) -> Vec<Self::Scalar> {
                $mul_func(x, y)
            }

            #[inline]
            fn zz_mul_scalar(x: &[Self::Scalar], scalar: i64) -> Vec<Self::Scalar> {
                let sc: Self::Scalar = Self::Scalar::from(scalar);
                x.into_iter().map(|c| *c * sc).collect()
            }

            fn new(coeffs: &[Self::Scalar]) -> Self {
                let mut ret = Self {
                    coeffs: [Self::Scalar::zero(); $params.sym_roots_num],
                };
                ret.coeffs.clone_from_slice(coeffs);
                ret
            }

            fn complex64(&self) -> Complex64 {
                let nums: Vec<Complex64> = self
                    .zz_coeffs()
                    .into_iter()
                    .map(|x| {
                        let re = x.to_f64().unwrap();
                        Complex64::new(re, 0.0)
                    })
                    .collect();
                let units: Vec<Complex64> = Self::zz_params()
                    .sym_roots_sqs
                    .into_iter()
                    .map(|x| Complex64::new(x.sqrt(), 0.0))
                    .collect();
                let mut ret = Complex64::zero();
                for (n, u) in nums.iter().zip(units.iter()) {
                    ret += n * u;
                }
                ret
            }
        }
        impl ZZBase for $name {
            type Scalar = GaussInt<Frac>;
            type Real = $name_real;

            fn conj(&self) -> Self {
                let cs: Vec<GInt> = self.zz_coeffs().iter().map(|c| c.conj()).collect();
                Self::new(&cs)
            }

            #[inline]
            fn zz_coeffs(&self) -> &[Self::Scalar] {
                &self.coeffs
            }

            #[inline]
            fn zz_coeffs_mut(&mut self) -> &mut [Self::Scalar] {
                &mut self.coeffs
            }

            #[inline]
            fn zz_params() -> &'static ZZParams<'static> {
                &$params
            }

            #[inline]
            fn zz_mul_arrays(x: &[Self::Scalar], y: &[Self::Scalar]) -> Vec<Self::Scalar> {
                $mul_func(x, y)
            }

            #[inline]
            fn zz_mul_scalar(x: &[Self::Scalar], scalar: i64) -> Vec<Self::Scalar> {
                let sc: Self::Scalar = Self::Scalar::from(scalar);
                x.into_iter().map(|c| *c * sc).collect()
            }

            fn new(coeffs: &[Self::Scalar]) -> Self {
                let mut ret = Self {
                    coeffs: [Self::Scalar::zero(); $params.sym_roots_num],
                };
                ret.coeffs.clone_from_slice(coeffs);
                ret
            }

            fn complex64(&self) -> Complex64 {
                let nums: Vec<Complex64> = self
                    .zz_coeffs()
                    .into_iter()
                    .map(|x| {
                        let re = x.real.to_f64().unwrap();
                        let im = x.imag.to_f64().unwrap();
                        Complex64::new(re, im)
                    })
                    .collect();
                let units: Vec<Complex64> = Self::zz_params()
                    .sym_roots_sqs
                    .into_iter()
                    .map(|x| Complex64::new(x.sqrt(), 0.0))
                    .collect();
                let mut ret = Complex64::zero();
                for (n, u) in nums.iter().zip(units.iter()) {
                    ret += n * u;
                }
                ret
            }
        }

        impl From<<$name_real as InnerIntType>::IntType> for $name_real {
            fn from(value: <GInt as InnerIntType>::IntType) -> Self {
                Self::one().scale(value)
            }
        }

        impl From<<$name as InnerIntType>::IntType> for $name {
            fn from(value: <GInt as InnerIntType>::IntType) -> Self {
                Self::one().scale(value)
            }
        }

        impl From<$name_real> for $name {
            /// Lift real-valued ring value into the corresponding cyclomatic ring
            fn from(value: $name_real) -> Self {
                let cs: Vec<<$name as ZZBase>::Scalar> =
                    value.zz_coeffs().iter().map(|z| (*z).into()).collect();
                Self::new(cs.as_slice())
            }
        }

        impl<T: Integer + ToPrimitive> From<(T, T)> for $name {
            /// Convert a Gaussian integer into a cyclotomic integer (only for fields including them)
            fn from((re, im): (T, T)) -> Self {
                let r = <$name as InnerIntType>::IntType::from_i64(re.to_i64().unwrap()).unwrap();
                let i = <$name as InnerIntType>::IntType::from_i64(im.to_i64().unwrap()).unwrap();
                Self::one().scale(r) + Self::one_i().scale(i)
            }
        }

        impl ZNum for $name_real {}
        impl ZNum for $name {}

        impl IntRing for $name_real {}
        impl IntRing for $name {}

        impl RealSigned for $name_real {
            fn re_signum(&self) -> Self {
                $re_signum_func(self)
            }
        }
        impl RealSigned for $name {
            fn re_signum(&self) -> Self {
                $re_signum_func(self)
            }
        }

        impl ZZComplex for $name {
            fn is_real(&self) -> bool {
                self.zz_coeffs().iter().all(|c| c.imag.is_zero())
            }

            fn is_imag(&self) -> bool {
                self.zz_coeffs().iter().all(|c| c.real.is_zero())
            }

            fn re(&self) -> <Self as ZZBase>::Real {
                let cs: Vec<Frac> = self.zz_coeffs().iter().map(|c| c.real).collect();
                $name_real::new(cs.as_slice())
            }

            fn im(&self) -> <Self as ZZBase>::Real {
                let cs: Vec<Frac> = self.zz_coeffs().iter().map(|c| c.imag).collect();
                $name_real::new(cs.as_slice())
            }
        }

        impl Ccw for $name {
            fn ccw() -> Self {
                Self::new(Self::zz_params().ccw_unit().as_slice())
            }
            fn is_ccw(&self) -> bool {
                *self == Self::ccw()
            }
        }
        impl CycIntRing for $name {}
        impl ZZNum for $name {}
    };
}

#[macro_export]
macro_rules! zz_ops_impl {
    ($($t:ty)*) => ($(
        impl Add<$t> for $t {
            type Output = Self;
            fn add(self, other: Self) -> Self {
                Self::new(&Self::zz_add(&self, &other))
            }
        }
        impl Sub<$t> for $t {
            type Output = Self;
            fn sub(self, other: Self) -> Self {
                Self::new(&Self::zz_sub(&self, &other))
            }
        }
        impl Neg for $t {
            type Output = Self;
            fn neg(self) -> Self {
                Self::new(&Self::zz_neg(&self))
            }
        }
        impl Mul<$t> for $t {
            type Output = Self;
            fn mul(self, other: Self) -> Self {
                Self::new(&Self::zz_mul(&self, &other))
            }
        }

        // TODO: refactor and implement this for QQn types only!
        // QQn should be defined as a newtype of ZZn, with extra ability to divide.
        impl Div<$t> for $t {
            type Output = Self;

            fn div(self, other: Self) -> Self {
                Self::mul(self, zz_inv(&other))
            }
        }

        impl Zero for $t {
            fn zero() -> Self {
                Self::new(&Self::zz_zero_vec())
            }
            fn is_zero(&self) -> bool {
                self.coeffs.to_vec() == Self::zz_zero_vec()
            }
        }
        impl One for $t {
            fn one() -> Self {
                Self::new(&Self::zz_one_vec())
            }
            fn is_one(&self) -> bool {
                self.coeffs.to_vec() == Self::zz_one_vec()
            }
        }

        impl Display for $t {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                let nums: Vec<String> = self.coeffs.into_iter().map(|x| format!("{x}")).collect();
                let units: Vec<String> = <$t>::zz_params().sym_roots_lbls.into_iter().map(|x| format!("sqrt({x})")).collect();
                let parts: Vec<String> = nums.iter().zip(units.iter()).filter(|(x, _)| x != &"0").map(|(x, y)| {
                    let is_real_unit = y == "sqrt(1)";
                    if (x == "1") {
                        if (is_real_unit) { "1".to_string() } else { y.to_string() }
                    } else if (is_real_unit) {
                        format!("{x}")
                    } else {
                        format!("({x})*{y}")
                    }
                }).collect();
                let joined = parts.join(" + ");
                let result = if (joined.is_empty()){ "0".to_string() } else { joined };
                return write!(f, "{result}");
            }
        }
    )*)
}
