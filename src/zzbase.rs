//! In this module the core traits for all rings and fields are defined.

use std::fmt::{Debug, Display};

use num_complex::Complex64;
use num_integer::Integer;
use num_rational::Ratio;
use num_traits::{FromPrimitive, One, PrimInt, ToPrimitive, Zero};

use crate::gaussint::GaussInt;
use crate::traits::{Ccw, Conj, InnerIntType, IntField, IntRing};
use crate::zzsigned::ZSigned;

/// We fix a general-purpose signed primitive integer size here. i64 is a natural choice.
// NOTE: If the need arises, everything could be refactored to parametrize over this type.
pub type MyInt = i64;

/// Internally, ZZn use a Gaussian integer where all coefficients are
/// a ratio with a fixed denominator (modulo simplification).
pub type Frac = Ratio<MyInt>;

/// Internally, the ring is constructed in a slightly non-standard way.
/// For implementation convenience, instead of having each root twice (real and imaginary),
/// this symmetry is pulled out by using Gaussian integers as coefficients for real-valued roots.
pub type GInt = GaussInt<Frac>;

/// Core parameters governing the behavior of a ring.
#[derive(Debug)]
pub struct ZZParams<'a> {
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
    /// Helper function to lift the param coeffs into a GaussInt of suitable type.
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

/// Assumptions about a scalar type to be used for coefficients representing ring elements.
pub trait ZScalar:
    IntField + InnerIntType + FromPrimitive + ZSigned + Conj + Debug + Display
{
}
impl<T: IntField + InnerIntType + FromPrimitive + ZSigned + Conj + Debug + Display> ZScalar for T {}

/// Common trait for various things required to implement cyclotomic rings.
pub trait ZZBase {
    type Scalar: ZScalar;
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
    ///
    /// Useful in the rare situation when this information is needed at runtime.
    #[inline]
    fn has_qturn() -> bool {
        Self::turn() % 4 == 0
    }

    /// Return angle representing a quarter turn (if ring supports it).
    /// **NOTE:** In unsuitable rings this value will be incorrect (see `has_turn()`).
    #[inline]
    fn qturn() -> i8 {
        Self::turn() / 4
    }

    /// Return angle representing a quarter turn (if ring supports it).
    #[inline]
    fn opt_qturn() -> Option<i8> {
        if Self::has_qturn() {
            Some(Self::qturn())
        } else {
            None
        }
    }

    /// Return unit length vector pointing in direction of given angle.
    fn unit(angle: i8) -> Self
    where
        Self: ZZNum,
    {
        Self::one() * Self::pow(&Self::ccw(), angle.rem_euclid(Self::turn()))
    }

    /// Raise to an integer power.
    // NOTE: using i8 instead of u8 for convenience (angles use i8)
    fn pow(&self, i: i8) -> Self
    where
        Self: ZCommon,
    {
        assert!(i >= 0, "Negative powers are not supported!");
        let mut x = Self::one();
        for _ in 0..i {
            x = x * (*self);
        }
        return x;
    }

    /// Scalar multiplication.
    #[inline]
    fn scale<I: Integer + ToPrimitive>(&self, scalar: I) -> Self
    where
        Self: Sized,
    {
        let cs: Vec<Self::Scalar> = Self::zz_mul_scalar(self.zz_coeffs(), scalar.to_i64().unwrap());
        Self::new(&cs)
    }

    /// Convert to a complex floating point number.
    fn complex64(&self) -> Complex64;

    // functions that can be implemented via zz_base_impl!
    // --------
    fn new(coeffs: &[Self::Scalar]) -> Self;

    fn zz_coeffs(&self) -> &[Self::Scalar];
    fn zz_coeffs_mut(&mut self) -> &mut [Self::Scalar];

    fn zz_params() -> &'static ZZParams<'static>;
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

    /// Return true if the value is mixed real **and** imaginary.
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

/// Cyclotomic ring or its real part.
pub trait ZCommon: ZZBase + InnerIntType + IntRing + Conj + Display {}

/// Quadratic real extension corresponding to a cyclotomic ring
/// (used to e.g. split a cyclotomic value into separate real and imaginary parts)
pub trait ZNum: ZCommon + ZSigned {}

/// A cyclotomic ring. You probably want to parametrize generic code over this trait.
pub trait ZZNum: ZCommon + ZZComplex + Ccw {}

// --------------------------------

#[macro_export]
macro_rules! zz_triv_impl {
    ($trait_name: ident, $($t:ty)*) => ($(
        impl $trait_name for $t {}
    )*)
}

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

        impl ZZBase for $name_real {
            type Scalar = Frac;
            type Real = $name_real;

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

        impl From<$name_real> for $name {
            /// Lift real-valued ring value into the corresponding cyclomatic ring
            fn from(value: $name_real) -> Self {
                let cs: Vec<<$name as ZZBase>::Scalar> =
                    value.zz_coeffs().iter().map(|z| (*z).into()).collect();
                Self::new(cs.as_slice())
            }
        }

        impl Conj for $name_real {
            fn conj(&self) -> Self {
                self.clone()
            }
            fn co_conj(&self) -> Self {
                self.neg()
            }
        }
        impl Conj for $name {
            fn conj(&self) -> Self {
                let cs: Vec<GInt> = self.zz_coeffs().iter().map(|c| c.conj()).collect();
                Self::new(&cs)
            }
            fn co_conj(&self) -> Self {
                let cs: Vec<GInt> = self.zz_coeffs().iter().map(|c| c.co_conj()).collect();
                Self::new(&cs)
            }
        }

        impl ZSigned for $name_real {
            fn signum(&self) -> Self {
                $re_signum_func(self)
            }
        }

        impl ZNum for $name_real {}

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

        impl InnerIntType for $t {
            type IntType = <<Self as ZZBase>::Scalar as InnerIntType>::IntType;
        }

        impl From<<$t as InnerIntType>::IntType> for $t {
            fn from(value: <<Self as ZZBase>::Scalar as InnerIntType>::IntType) -> Self {
                Self::one().scale(value)
            }
        }

        impl IntRing for $t {}
        impl ZCommon for $t {}
    )*)
}
