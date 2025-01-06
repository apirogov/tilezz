use super::gaussint::GaussInt;
use super::traits::{Ccw, ComplexIntRing, InnerIntType, IntRing};
use num_complex::Complex64;
use num_integer::Integer;
use num_rational::Ratio;
use num_traits::{FromPrimitive, One, PrimInt, Signed, ToPrimitive, Zero};
use std::fmt::Debug;
use std::fmt::Display;
use std::marker::PhantomData;

// configure a signed primitive integer size here
pub type MyInt = i64;

// internally, ZZn use a Gaussian integer where all coefficients are
// a ratio with a fixed denominator (modulo simplification)
pub type Frac = Ratio<MyInt>;
pub type GInt = GaussInt<Frac>;

#[derive(Debug)]
pub struct ZZParams<'a, T> {
    pub phantom: PhantomData<&'a T>,

    /// Squares of symbolic roots.
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

impl<T> ZZParams<'static, Ratio<T>>
where
    T: PrimInt + Integer + Signed + IntRing,
{
    // helper function to lift the param coeffs into a GaussInt of suitable type
    pub fn ccw_unit(self) -> Vec<GaussInt<Ratio<T>>> {
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

/// Floating-point-free solution to get sign of an expression
/// a*sqrt(n) + b*sqrt(m)
/// where a,b,m,n are all integers
/// and m,n are squarefree coprime constants > 1.
pub fn signum_sum_sqrt_expr_2<T: IntRing + Signed + Copy + FromPrimitive + Debug>(
    a: T,
    m: T,
    b: T,
    n: T,
) -> T {
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
pub fn signum_sum_sqrt_expr_4<T: IntRing + Signed + Copy + FromPrimitive + Debug>(
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

/// This implementation uses f64 (works in all ZZ rings).
pub fn zz_partial_signum_fallback<ZZ: ZZNum>(val: &ZZ) -> Frac {
    let val = val.complex().re;
    if val.is_zero() {
        Frac::zero()
    } else {
        Frac::from_f64(val.signum()).unwrap()
    }
}

/// Trivial if there is just one term in the expression,
/// regardless of its symbolic root (which usually will be 1).
pub fn zz_partial_signum_1_sym<ZZ: ZZNum>(val: &ZZ) -> Frac {
    let cs = val.zz_coeffs();

    cs[0].real.signum()
}

pub fn zz_partial_signum_2_sym<ZZ: ZZNum>(val: &ZZ) -> Frac {
    let ftoi = |f| Frac::from_f64(f).unwrap();
    let cs = val.zz_coeffs();
    let rs = ZZ::zz_params().sym_roots_sqs;

    let (a, b) = (cs[0].real, cs[1].real);
    let (m, n) = (ftoi(rs[0]), ftoi(rs[1]));
    signum_sum_sqrt_expr_2(a, m, b, n)
}

pub fn zz_partial_signum_4_sym<ZZ: ZZNum>(val: &ZZ) -> Frac {
    let cs: Vec<Frac> = val.zz_coeffs().iter().map(|x| x.real).collect();
    let rs: Vec<Frac> = ZZ::zz_params()
        .sym_roots_sqs
        .iter()
        .map(|r| Frac::from_f64(*r).unwrap())
        .collect();

    let (a, b, c, d) = (cs[0], cs[1], cs[2], cs[3]);
    let (k, m, n, l) = (rs[0], rs[1], rs[2], rs[3]);
    signum_sum_sqrt_expr_4(a, k, b, m, c, n, d, l)
}

pub trait ZZBase<
    T: Sized
        + Signed
        + PartialOrd
        + IntRing
        + InnerIntType
        + ToPrimitive
        + FromPrimitive
        + Debug
        + 'static,
>
{
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
        let j = i.rem_euclid(<Self as ZZBase<T>>::turn());
        one * <Self as ZZBase<T>>::powi(&ccw, j)
    }

    /// Raise to an integer power.
    fn powi(&self, i: i8) -> Self
    where
        Self: ZZNum,
    {
        if i < 0 {
            panic!("Negative powers are not supported!");
        }
        let mut x = Self::one();
        for _ in 0..i {
            x = x * (*self);
        }
        return x;
    }

    /// Return imaginary unit (if ring supports it).
    #[inline]
    fn one_i() -> Self
    where
        Self: ZZNum,
    {
        let qt = <Self as ZZBase<T>>::qturn();
        <Self as ZZBase<T>>::unit(qt)
    }

    /// Complex conjugation.
    #[inline]
    fn conj(&self) -> Self
    where
        Self: Sized,
    {
        let cs: Vec<GaussInt<T>> = self.zz_coeffs().iter().map(|c| c.conj()).collect();
        Self::new(&cs)
    }

    /// Scalar multiplication.
    #[inline]
    fn scale(&self, scalar: <GaussInt<T> as InnerIntType>::IntType) -> Self
    where
        Self: Sized,
    {
        let cs: Vec<GaussInt<T>> = Self::zz_mul_scalar(self.zz_coeffs(), scalar);
        Self::new(&cs)
    }

    /// Convert to a complex floating point number.
    fn complex(&self) -> Complex64 {
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

    /// Return true if the value is purely real.
    fn is_real(&self) -> bool {
        self.zz_coeffs().iter().all(|c| c.imag.is_zero())
    }

    /// Return true if the value is purely imaginary.
    fn is_imag(&self) -> bool {
        self.zz_coeffs().iter().all(|c| c.real.is_zero())
    }

    fn is_complex(&self) -> bool {
        !self.is_real() && !self.is_imag()
    }

    /// Split the value into its real and imaginary contributions.
    /// Note that the imaginary component is converted to real.
    fn xy(&self) -> (Self, Self)
    where
        Self: ZZNum,
    {
        let cs_x: Vec<GaussInt<T>> = self
            .zz_coeffs()
            .iter()
            .map(|c| GaussInt::new(c.real, T::zero()))
            .collect();
        let cs_y: Vec<GaussInt<T>> = self
            .zz_coeffs()
            .iter()
            .map(|c| GaussInt::new(c.imag, T::zero()))
            .collect();
        (Self::new(&cs_x), Self::new(&cs_y))
    }

    /// Return the sign of the real part of the value.
    /// Note that ZZ cannot not support Signed trait in general.
    ///
    /// Implementation depends on the underlying ring.
    fn re_signum(&self) -> T
    where
        Self: ZZNum;

    /// Get (a, b, c) for line ax + by + c = 0 based on two points.
    fn line_through(&self, other: &Self) -> (Self, Self, Self)
    where
        Self: ZZNum,
    {
        let (p1x, p1y) = <Self as ZZBase<T>>::xy(self);
        let (p2x, p2y) = <Self as ZZBase<T>>::xy(other);
        (p1y - p2y, p2x - p1x, p1x * p2y - p2x * p1y)
    }

    /// Return whether the point is on the given line.
    fn is_colinear(&self, (a, b, c): &(Self, Self, Self)) -> bool
    where
        Self: ZZNum,
    {
        let (px, py) = <Self as ZZBase<T>>::xy(self);
        (*a * px + *b * py + *c).is_zero()
    }

    /// Dot product between two values.
    ///
    /// Result is always real-valued.
    fn dot(&self, other: &Self) -> Self
    where
        Self: ZZNum,
    {
        let (p1x, p1y) = <Self as ZZBase<T>>::xy(self);
        let (p2x, p2y) = <Self as ZZBase<T>>::xy(other);
        (p1x * p2x) + (p1y * p2y)
    }

    /// Squared norm, i.e. dot product of value with itself.
    fn norm_sq(&self) -> Self
    where
        Self: ZZNum,
    {
        <Self as ZZBase<T>>::dot(self, self)
    }

    /// 2D "cross" product (wedge is the general name).
    ///
    /// Result is always real-valued.
    fn wedge(&self, other: &Self) -> Self
    where
        Self: ZZNum,
    {
        let (p1x, p1y) = <Self as ZZBase<T>>::xy(self);
        let (p2x, p2y) = <Self as ZZBase<T>>::xy(other);
        (p1x * p2y) - (p1y * p2x)
    }

    /// Return whether this point is strictly between the other two.
    /// Note that we assume all three involved points are colinear.
    fn is_between(&self, p1: &Self, p2: &Self) -> bool
    where
        Self: ZZNum,
    {
        // Alternative approach (from Wilderberger, "Divine Proportions"):
        /*
        let l = <Self as ZZBase<T>>::line_through(&p1, &p2);
        if !<Self as ZZBase<T>>::is_colinear(&self, &l) {
            return false;
        }

        // NOTE: self = p3, i.e. check that p1 ---- p3 ---- p2
        let q1 = <Self as ZZBase<T>>::q(&(*self - *p1));
        let q2 = <Self as ZZBase<T>>::q(&(*self - *p2));
        let q3 = <Self as ZZBase<T>>::q(&(*p1 - *p2));

        // if d > 0 -> p3 on segment (p1, p2)
        // if d = 0 -> p3 is equal to an endpoint
        // if d < 0 -> p3 colinear, but outside
        let d = q3 - q1 - q2; // = +- 2*sqrt(q1*q2)
        let sgn = <Self as ZZBase<T>>::zz_partial_signum(&d);
        return sgn.is_one();
        */

        // classical approach:
        // ----
        let v = *p1 - *self;
        let w = *self - *p2;
        let wed = <Self as ZZBase<T>>::wedge(&v, &w);
        let dot = <Self as ZZBase<T>>::dot(&v, &w);
        let dot_sign = <Self as ZZBase<T>>::re_signum(&dot);
        wed.is_zero() && dot_sign.is_one()
    }

    // functions that can be implemented via zz_base_impl!
    // --------
    fn new(coeffs: &[GaussInt<T>]) -> Self;

    fn zz_coeffs(&self) -> &[GaussInt<T>];
    fn zz_coeffs_mut(&mut self) -> &mut [GaussInt<T>];

    fn zz_params() -> &'static ZZParams<'static, T>;
    fn zz_ccw_vec() -> Vec<GaussInt<T>>;
    fn zz_mul_arrays(x: &[GaussInt<T>], y: &[GaussInt<T>]) -> Vec<GaussInt<T>>;
    fn zz_mul_scalar(
        arr: &[GaussInt<T>],
        scalar: <GaussInt<T> as InnerIntType>::IntType,
    ) -> Vec<GaussInt<T>>;

    // implementations for implementing other traits using zz_ops_impl!
    // --------
    #[inline]
    fn zz_zero_vec() -> Vec<GaussInt<T>> {
        vec![GaussInt::zero(); Self::zz_params().sym_roots_num]
    }

    fn zz_one_vec() -> Vec<GaussInt<T>> {
        let mut ret = vec![GaussInt::zero(); Self::zz_params().sym_roots_num];
        ret[0] = GaussInt::one();
        ret
    }

    fn zz_add(&self, other: &Self) -> Vec<GaussInt<T>> {
        let mut ret = Self::zz_zero_vec();
        for (i, (aval, bval)) in self.zz_coeffs().iter().zip(other.zz_coeffs()).enumerate() {
            ret[i] = *aval + *bval;
        }
        ret
    }

    fn zz_sub(&self, other: &Self) -> Vec<GaussInt<T>> {
        let mut ret = Self::zz_zero_vec();
        for (i, (aval, bval)) in self.zz_coeffs().iter().zip(other.zz_coeffs()).enumerate() {
            ret[i] = *aval - *bval;
        }
        ret
    }

    fn zz_neg(&self) -> Vec<GaussInt<T>> {
        let mut ret = Self::zz_zero_vec();
        for (i, val) in self.zz_coeffs().iter().enumerate() {
            ret[i] = -(*val);
        }
        ret
    }

    #[inline]
    fn zz_mul(&self, other: &Self) -> Vec<GaussInt<T>> {
        Self::zz_mul_arrays(self.zz_coeffs(), other.zz_coeffs())
    }
}

/// The trait that other structures should parametrize on
/// to work with any complex integer ring.
pub trait ZZNum: ZZBase<Frac> + InnerIntType + ComplexIntRing + Display {}

#[macro_export]
macro_rules! zz_base_impl {
    ($name:ident, $params:ident, $mul_func:ident, $re_signum_func:ident) => {
        #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
        pub struct $name {
            coeffs: [GInt; $params.sym_roots_num],
        }

        impl InnerIntType for $name {
            type IntType = <GInt as InnerIntType>::IntType;
        }

        impl From<<GInt as InnerIntType>::IntType> for $name {
            fn from(value: <GInt as InnerIntType>::IntType) -> Self {
                Self::one().scale(value)
            }
        }

        impl<T: Integer + ToPrimitive> From<(T, T)> for $name {
            fn from((re, im): (T, T)) -> Self {
                let r = <GInt as InnerIntType>::IntType::from_i64(re.to_i64().unwrap()).unwrap();
                let i = <GInt as InnerIntType>::IntType::from_i64(im.to_i64().unwrap()).unwrap();
                Self::one().scale(r) + Self::one_i().scale(i)
            }
        }

        impl ZZBase<Frac> for $name {
            #[inline]
            fn zz_coeffs(&self) -> &[GInt] {
                &self.coeffs
            }

            #[inline]
            fn zz_coeffs_mut(&mut self) -> &mut [GInt] {
                &mut self.coeffs
            }

            #[inline]
            fn zz_params() -> &'static ZZParams<'static, Frac> {
                &$params
            }

            #[inline]
            fn zz_ccw_vec() -> Vec<GInt> {
                $params.ccw_unit()
            }

            #[inline]
            fn zz_mul_arrays(x: &[GInt], y: &[GInt]) -> Vec<GInt> {
                $mul_func(x, y)
            }

            #[inline]
            fn re_signum(&self) -> Frac
            where
                Self: ZZNum,
            {
                $re_signum_func(self)
            }

            #[inline]
            fn zz_mul_scalar(x: &[GInt], scalar: i64) -> Vec<GInt> {
                let sc: GInt = GInt::from(scalar);
                x.into_iter().map(|c| *c * sc).collect()
            }

            fn new(coeffs: &[GInt]) -> Self {
                let mut ret = Self {
                    coeffs: [GInt::zero(); $params.sym_roots_num],
                };
                ret.coeffs.clone_from_slice(coeffs);
                ret
            }
        }

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
        impl Ccw for $t {
            fn ccw() -> Self {
                Self::new(&Self::zz_ccw_vec())
            }
            fn is_ccw(&self) -> bool {
                self.coeffs.to_vec() == Self::zz_ccw_vec()
            }
        }

        impl IntRing for $t {}
        impl ComplexIntRing for $t {}

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
