use super::gaussint::GaussInt;
use super::traits::{ComplexIntRing, InnerIntType, IntRing};
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
    pub fn ccw_unit(self, scaling_fac: T) -> Vec<GaussInt<Ratio<T>>> {
        self.ccw_unit_coeffs
            .into_iter()
            .map(|x| {
                GaussInt::new(
                    Ratio::<T>::new_raw(T::from(x[0]).unwrap(), scaling_fac.clone()),
                    Ratio::<T>::new_raw(T::from(x[1]).unwrap(), scaling_fac.clone()),
                )
            })
            .collect()
    }
}

pub trait ZZBase<
    T: Signed + PartialOrd + IntRing + InnerIntType + ToPrimitive + FromPrimitive + Debug + 'static,
>
{
    #[inline]
    fn turn() -> i8 {
        Self::zz_params().full_turn_steps
    }
    #[inline]
    fn hturn() -> i8 {
        Self::turn() / 2
    }

    /// Return quarter turn angle (if ring supports it).
    #[inline]
    fn qturn() -> i8 {
        assert_eq!(Self::turn() % 4, 0);
        Self::turn() / 4
    }

    /// Return imaginary unit (if ring supports it).
    #[inline]
    fn one_i() -> Self
    where
        Self: Sized,
    {
        Self::unit(Self::qturn())
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

    /// Construct a new ZZ number by lifting an integer.
    #[inline]
    fn from_int(scalar: <GaussInt<T> as InnerIntType>::IntType) -> Self
    where
        Self: Sized + One,
    {
        Self::one().scale(scalar)
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

    // --------

    /// This implementation uses f64 (works in all ZZ rings).
    fn zz_partial_signum_fallback(&self) -> T
    where
        Self: ZZNum,
    {
        let val = <Self as ZZBase<T>>::complex(self).re;
        if val.is_zero() {
            T::zero()
        } else {
            T::from_f64(val.signum()).unwrap()
        }
    }

    /// Trivial implementation if there is just one component,
    /// regardless of its symbolic root (which usually will be 1).
    fn zz_partial_signum_1_sym(&self) -> T
    where
        Self: ZZNum,
    {
        let cs = <Self as ZZBase<T>>::zz_coeffs(self);
        cs[0].real.signum()
    }

    /// Solution without using floating point operations.
    ///
    /// z = a*sqrt(n) + b*sqrt(m) > 0 can be checked based
    /// on the squares, assuming that n,m  are integers > 0:
    /// if a and b are both positive or negative -> trivial,
    /// if a>0, b<0: z > 0 <=> n*a^2 > m*b^2 <=> n*a^2 - m*b^2 > 0
    /// if a<0, b>0: z > 0 <=> n*a^2 < m*b^2 <=> m*b^2 - n*a^2 > 0
    /// (and symmetrically for z < 0), which can be summarized
    /// by using the sign of a and b to avoid case splitting.
    fn zz_partial_signum_2_sym(&self) -> T
    where
        Self: ZZNum,
    {
        let p: &ZZParams<'_, T> = Self::zz_params();
        let n = T::from_f64(p.sym_roots_sqs[0]).unwrap();
        let m = T::from_f64(p.sym_roots_sqs[1]).unwrap();

        let cs = <Self as ZZBase<T>>::zz_coeffs(self);
        let a = cs[0].real;
        let b = cs[1].real;

        let signed_n_a_sq = n * a * a * a.signum();
        let signed_m_b_sq = m * b * b * b.signum();
        (signed_n_a_sq + signed_m_b_sq).signum()
    }

    /// Return the sign of the value (assuming it is real-valued).
    /// Note that ZZ cannot not support Signed trait in general.
    fn partial_signum(&self) -> T
    where
        Self: ZZNum,
    {
        assert!(<Self as ZZBase<T>>::is_real(self));

        // TODO: implement special case for 4 (to get ZZ24 float-free)
        match <Self as ZZBase<T>>::zz_params().sym_roots_num {
            1 => <Self as ZZBase<T>>::zz_partial_signum_1_sym(&self),
            2 => <Self as ZZBase<T>>::zz_partial_signum_2_sym(&self),
            _ => <Self as ZZBase<T>>::zz_partial_signum_fallback(&self),
        }
    }

    // --------

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
        let sgn = <Self as ZZBase<T>>::partial_signum(&d);
        return sgn.is_one();
        */

        // classical approach:
        // ----
        let v = *p1 - *self;
        let w = *self - *p2;
        let wed = <Self as ZZBase<T>>::wedge(&v, &w);
        let dot = <Self as ZZBase<T>>::dot(&v, &w);
        let dot_sign = <Self as ZZBase<T>>::partial_signum(&dot);
        wed.is_zero() && dot_sign.is_one()
    }

    // functions that can be implemented via zz_base_impl!
    // --------
    fn new(coeffs: &[GaussInt<T>]) -> Self;
    fn unit(i: i8) -> Self;
    fn powi(&self, i: i8) -> Self;

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

#[macro_export]
macro_rules! zz_base_impl {
    ($name:ident, $params:ident, $mul_func:ident) => {
        #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
        pub struct $name {
            coeffs: [GInt; $params.sym_roots_num],
        }

        impl InnerIntType for $name {
            type IntType = <GInt as InnerIntType>::IntType;
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
                $params.ccw_unit($params.scaling_fac)
            }

            #[inline]
            fn zz_mul_arrays(x: &[GInt], y: &[GInt]) -> Vec<GInt> {
                $mul_func(x, y)
            }

            #[inline]
            fn zz_mul_scalar(x: &[GInt], scalar: i64) -> Vec<GInt> {
                let sc: GInt = to_gint(scalar);
                x.into_iter().map(|c| *c * sc).collect()
            }

            fn new(coeffs: &[GInt]) -> Self {
                let mut ret = Self {
                    coeffs: [GInt::zero(); $params.sym_roots_num],
                };
                ret.coeffs.clone_from_slice(coeffs);
                ret
            }

            fn unit(i: i8) -> Self {
                return Self::one() * Self::ccw().powi(i.rem_euclid(Self::turn()));
            }

            fn powi(&self, i: i8) -> Self {
                if (i < 0) {
                    panic!("Negative powers are not supported!");
                }
                let mut x = Self::one();
                for _ in 0..i {
                    x = x * (*self);
                }
                return x;
            }
        }
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
                let units: Vec<String> = <$t>::zz_params().sym_roots_sqs.into_iter().map(|x| format!("sqrt({x})")).collect();
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

/// The trait that other structures should parametrize on
/// to work with any complex integer ring.
pub trait ZZNum: ZZBase<Frac> + InnerIntType + ComplexIntRing + Display {}
