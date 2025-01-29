//! In this module the core traits for all rings and fields are defined.

use std::fmt::{Debug, Display};
use std::hash::Hash;

use num_complex::Complex64;
use num_integer::Integer;
use num_rational::Ratio;
use num_traits::{FromPrimitive, One, PrimInt, ToPrimitive, Zero};

use super::gaussint::GaussInt;
use super::numtraits::{IntField, IntRing, ZSigned};
use super::traits::{IsComplex, IsRingOrField};

/// We fix a general-purpose signed primitive integer size here. i64 is a natural choice.
// NOTE: If the need arises, everything could be refactored to parametrize over this type.
pub type IntT = i64;

/// Internally, ZZn use a Gaussian integer where all coefficients are
/// a ratio with a fixed denominator (modulo simplification).
pub type RatioT = Ratio<IntT>;

/// Internally, the ring is constructed in a slightly non-standard way.
/// For implementation convenience, instead of having each root twice (real and imaginary),
/// this symmetry is pulled out by using Gaussian integers as coefficients for real-valued roots.
pub type GIntT = GaussInt<RatioT>;

/// Core parameters governing the behavior of a cyclotomic ring.
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
    pub scaling_fac: IntT,
    /// Unit of rotation in coefficients of this ZZ type
    pub ccw_unit_coeffs: &'a [[IntT; 2]],
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
pub trait SymScalar: IntField + FromPrimitive + ZSigned + Debug + Display {}
impl<T: IntField + FromPrimitive + ZSigned + Debug + Display> SymScalar for T {}

/// Common trait for various things required to implement cyclotomic rings.
pub trait SymNum: Clone + Copy + PartialEq + Eq + Hash + Debug + Display {
    /// Type used for coefficients for the symbolic vectors.
    type Scalar: SymScalar;

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
        Self: IsComplex + IsRingOrField,
    {
        Self::one() * Self::pow(&Self::ccw(), angle.rem_euclid(Self::turn()))
    }

    /// Raise to an integer power.
    // NOTE: using i8 instead of u8 for convenience (angles use i8)
    fn pow(&self, i: i8) -> Self
    where
        Self: IsRingOrField,
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

    /// Convert to two coordinates.
    fn xy(&self) -> (f64, f64) {
        let c = self.complex64();
        (c.re, c.im)
    }

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

// ----------------

#[macro_export]
macro_rules! zz_triv_impl {
    ($trait_name: ident, $($t:ty)*) => ($(
        impl $trait_name for $t {}
    )*)
}

macro_rules! impl_primint_traits {
    ($t:ident) => {
        impl InnerIntType for $t {
            type IntType = <<Self as SymNum>::Scalar as InnerIntType>::IntType;
        }

        impl From<<$t as InnerIntType>::IntType> for $t {
            fn from(value: <Self as InnerIntType>::IntType) -> Self {
                Self::one().scale(value)
            }
        }
    };
}

macro_rules! impl_symnum_q_from_z {
    ($q:ident, $z:ident) => {

        impl Display for $z {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                let nums: Vec<String> = self.coeffs.into_iter().map(|x| format!("{x}")).collect();
                let units: Vec<String> = <$z>::zz_params()
                    .sym_roots_lbls
                    .into_iter()
                    .map(|x| format!("sqrt({x})"))
                    .collect();
                let parts: Vec<String> = nums
                    .iter()
                    .zip(units.iter())
                    .filter(|(x, _)| x != &"0")
                    .map(|(x, y)| {
                        let is_real_unit = y == "sqrt(1)";
                        if (x == "1") {
                            if (is_real_unit) {
                                "1".to_string()
                            } else {
                                y.to_string()
                            }
                        } else if (is_real_unit) {
                            format!("{x}")
                        } else {
                            format!("({x})*{y}")
                        }
                    })
                    .collect();
                let joined = parts.join(" + ");
                let result = if (joined.is_empty()) {
                    "0".to_string()
                } else {
                    joined
                };
                return write!(f, "{result}");
            }
        }

        impl Display for $q {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                self.0.fmt(f)
            }
        }

        impl SymNum for $q {
            type Scalar = <$z as SymNum>::Scalar;

            #[inline]
            fn new(coeffs: &[Self::Scalar]) -> Self {
                Self($z::new(coeffs))
            }
            #[inline]
            fn complex64(&self) -> Complex64 {
                self.0.complex64()
            }
            #[inline]
            fn zz_coeffs(&self) -> &[Self::Scalar] {
                self.0.zz_coeffs()
            }
            #[inline]
            fn zz_coeffs_mut(&mut self) -> &mut [Self::Scalar] {
                self.0.zz_coeffs_mut()
            }
            #[inline]
            fn zz_params() -> &'static ZZParams<'static> {
                $z::zz_params()
            }
            #[inline]
            fn zz_mul_arrays(x: &[Self::Scalar], y: &[Self::Scalar]) -> Vec<Self::Scalar> {
                $z::zz_mul_arrays(x, y)
            }
            #[inline]
            fn zz_mul_scalar(x: &[Self::Scalar], scalar: i64) -> Vec<Self::Scalar> {
                $z::zz_mul_scalar(x, scalar)
            }
        }
    };
}

#[macro_export]
macro_rules! impl_symnum {
    ($name:ident, $name_real: ident, $qname: ident, $qname_real: ident, $params:ident, $mul_func:ident) => {
        #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
        pub struct $name_real {
            coeffs: [RatioT; $params.sym_roots_num],
        }
        #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
        pub struct $name {
            coeffs: [GIntT; $params.sym_roots_num],
        }
        #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
        pub struct $qname_real(pub $name_real);
        #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
        pub struct $qname(pub $name);

        impl SymNum for $name_real {
            type Scalar = RatioT;

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

        impl SymNum for $name {
            type Scalar = GaussInt<RatioT>;

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

        impl_symnum_q_from_z!($qname, $name);
        impl_symnum_q_from_z!($qname_real, $name_real);
    };
}

// ----------------

#[macro_export]
macro_rules! impl_real_traits {
    ($t:ident, $signum_func:ident) => {
        impl ZSigned for $t {
            fn signum(&self) -> Self {
                $signum_func(self)
            }
        }

        impl Ord for $t {
            fn cmp(&self, other: &Self) -> Ordering {
                let diff = *self - *other;
                if diff.is_zero() {
                    return Ordering::Equal;
                }
                if diff.is_positive() {
                    Ordering::Greater
                } else {
                    Ordering::Less
                }
            }
        }

        impl PartialOrd for $t {
            fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
                Some(self.cmp(other))
            }
        }

        impl RealTraits for $t {}
    };
}

// ----------------

pub trait ZZComplex: IsRingOrField {
    /// Return true if the value is purely real.
    fn is_real(&self) -> bool;

    /// Return true if the value is purely imaginary.
    fn is_imag(&self) -> bool;

    /// Return true if the value is mixed real **and** imaginary.
    fn is_complex(&self) -> bool {
        !self.is_real() && !self.is_imag()
    }

    // FIXME: projection is NOT supported in general by the ring, needs field!
    // --------
    /// Return the real part of the value,
    /// i.e. the value (z + z.conj()) / 2
    fn re(&self) -> <Self as IsRingOrField>::Real
    where
        Self: SymNum;

    /// Return the imaginary part of the value (rotated onto the real axis),
    /// i.e. the value (z - z.conj()) / 2i
    fn im(&self) -> <Self as IsRingOrField>::Real
    where
        Self: SymNum;
    // --------

    /// Split the value into its real and imaginary contributions.
    /// Note that the imaginary component is converted to real.
    ///
    // NOTE: z = dot(1, z) + i*wedge(1, z), as dot(1, z) = Re(z), wedge(1, z) = Im(z)
    fn re_im(&self) -> (<Self as IsRingOrField>::Real, <Self as IsRingOrField>::Real)
    where
        Self: SymNum,
    {
        (self.re(), self.im())
    }
}

#[macro_export]
macro_rules! impl_complex_traits {
    ($t:ident) => {
        impl Ccw for $t {
            fn ccw() -> Self {
                Self::new(Self::zz_params().ccw_unit().as_slice())
            }
            fn is_ccw(&self) -> bool {
                *self == Self::ccw()
            }
        }

        impl ZZComplex for $t {
            fn is_real(&self) -> bool {
                self.zz_coeffs().iter().all(|c| c.imag.is_zero())
            }

            fn is_imag(&self) -> bool {
                self.zz_coeffs().iter().all(|c| c.real.is_zero())
            }

            fn re(&self) -> <Self as IsRingOrField>::Real {
                let cs: Vec<RatioT> = self.zz_coeffs().iter().map(|c| c.real).collect();
                <Self as IsRingOrField>::Real::new(cs.as_slice())
            }

            fn im(&self) -> <Self as IsRingOrField>::Real {
                let cs: Vec<RatioT> = self.zz_coeffs().iter().map(|c| c.imag).collect();
                <Self as IsRingOrField>::Real::new(cs.as_slice())
            }
        }

        impl ComplexTraits for $t {}
    };
}

// ----------------

macro_rules! impl_conj {
    ($name:ident, $name_real: ident) => {
        impl Conj for $name_real {
            fn conj(&self) -> Self {
                self.clone()
            }
        }
        impl Conj for $name {
            fn conj(&self) -> Self {
                let cs: Vec<GIntT> = self.zz_coeffs().iter().map(|c| c.conj()).collect();
                Self::new(&cs)
            }
        }
    };
}

#[macro_export]
macro_rules! impl_intring_traits {
    ($t:ty) => {
        impl Neg for $t {
            type Output = Self;
            fn neg(self) -> Self {
                Self::new(&Self::zz_neg(&self))
            }
        }
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
                self.zz_coeffs().to_vec() == Self::zz_zero_vec()
            }
        }
        impl One for $t {
            fn one() -> Self {
                Self::new(&Self::zz_one_vec())
            }
            fn is_one(&self) -> bool {
                self.zz_coeffs().to_vec() == Self::zz_one_vec()
            }
        }

        impl IntRing for $t {}
    };
}

#[macro_export]
macro_rules! impl_ring_traits {
    // ($name:ident, $name_real: ident, $qname: ident, $qname_real: ident, $params:ident, $mul_func:ident, $re_signum_func:ident) => {
    ($name:ident, $name_real: ident, $qname: ident, $qname_real: ident) => {
        impl_primint_traits!($name_real);
        impl_primint_traits!($qname_real);
        impl_primint_traits!($name);
        impl_primint_traits!($qname);

        impl_intring_traits!($name_real);
        impl_intring_traits!($qname_real);
        impl_intring_traits!($name);
        impl_intring_traits!($qname);

        impl_conj!($name, $name_real);
        impl_conj!($qname, $qname_real);

        impl RingTraits for $name_real {}
        impl RingTraits for $qname_real {}
        impl RingTraits for $name {}
        impl RingTraits for $qname {}
    };
}

#[macro_export]
macro_rules! impl_field_traits {
    ($q:ty) => {
        impl Div<$q> for $q {
            type Output = Self;
            fn div(self, other: Self) -> Self {
                Self::from(self.0.mul(zz_inv(&other.0)))
            }
        }

        impl IntField for $q {}

        impl FieldTraits for $q {
            fn coerce_ring(self) -> <Self as IsRealOrComplex>::Ring {
                self.0
            }
        }
    };
}

// --------
#[macro_export]
macro_rules! impl_functional_traits {
    ($name:ident, $name_real: ident, $qname: ident, $qname_real: ident, $signum_func:ident) => {
        impl_ring_traits!($name, $name_real, $qname, $qname_real);

        impl_field_traits!($qname);
        impl_field_traits!($qname_real);

        impl_real_traits!($name_real, $signum_func);
        impl_real_traits!($qname_real, $signum_func);

        impl_complex_traits!($name);
        impl_complex_traits!($qname);

        impl IsRingOrField for $name {
            type Real = $name_real;
            type Complex = $name;
        }
        impl IsRingOrField for $name_real {
            type Real = $name_real;
            type Complex = $name;
        }
        impl IsRingOrField for $qname {
            type Real = $qname_real;
            type Complex = $qname;
        }
        impl IsRingOrField for $qname_real {
            type Real = $qname_real;
            type Complex = $qname;
        }

        impl IsRealOrComplex for $name {
            type Ring = $name;
            type Field = $qname;
        }
        impl IsRealOrComplex for $qname {
            type Ring = $name;
            type Field = $qname;
        }
        impl IsRealOrComplex for $name_real {
            type Ring = $name_real;
            type Field = $qname_real;
        }
        impl IsRealOrComplex for $qname_real {
            type Ring = $name_real;
            type Field = $qname_real;
        }

        impl IsRing for $name {
            type Real = $name_real;
            type Complex = $name;
        }
        impl IsRing for $name_real {
            type Real = $name_real;
            type Complex = $name;
        }

        impl IsField for $qname {
            type Real = $qname_real;
            type Complex = $qname;
        }
        impl IsField for $qname_real {
            type Real = $qname_real;
            type Complex = $qname;
        }

        impl IsReal for $name_real {
            type Ring = $name_real;
            type Field = $qname_real;
        }
        impl IsReal for $qname_real {
            type Ring = $name_real;
            type Field = $qname_real;
        }

        impl IsComplex for $name {
            type Ring = $name;
            type Field = $qname;
        }
        impl IsComplex for $qname {
            type Ring = $name;
            type Field = $qname;
        }

        impl ZType for $name_real {
            type Complex = $name;
            type Field = $qname_real;
        }
        impl ZZType for $name {
            type Real = $name_real;
            type Field = $qname;
        }
        impl QType for $qname_real {
            type Complex = $qname;
            type Ring = $name_real;
        }
        impl QQType for $qname {
            type Real = $qname_real;
            type Ring = $name;
        }
    };
}

#[macro_export]
macro_rules! impl_conversions {
    ($zz:ident, $z:ident, $qq:ident, $q:ident) => {
        impl From<$z> for $zz {
            /// Lift real-valued ring value into the corresponding cyclomatic ring
            fn from(value: $z) -> Self {
                let cs: Vec<<$zz as SymNum>::Scalar> =
                    value.zz_coeffs().iter().map(|z| (*z).into()).collect();
                Self::new(cs.as_slice())
            }
        }
        impl From<$z> for $qq {
            /// Lift real-valued ring value into the corresponding cyclomatic ring
            fn from(value: $z) -> Self {
                let cs: Vec<<$zz as SymNum>::Scalar> =
                    value.zz_coeffs().iter().map(|z| (*z).into()).collect();
                Self::new(cs.as_slice())
            }
        }
        impl From<$q> for $qq {
            /// Lift real-valued ring value into the corresponding cyclomatic ring
            fn from(value: $q) -> Self {
                let cs: Vec<<$qq as SymNum>::Scalar> =
                    value.zz_coeffs().iter().map(|z| (*z).into()).collect();
                Self::new(cs.as_slice())
            }
        }

        impl From<$z> for $q {
            /// Lift ring value into the corresponding cyclomatic field
            fn from(value: $z) -> Self {
                Self(value)
            }
        }
        impl From<$zz> for $qq {
            /// Lift ring value into the corresponding cyclomatic field
            fn from(value: $zz) -> Self {
                Self(value)
            }
        }
    };
}

// --------

#[macro_export]
macro_rules! zz_test {
    ($name:ident, $type:ty) => {
        mod $name {
            use super::*;

            type ZZi = $type;

            #[test]
            fn test_units_sum_is_complex() {
                let p: ZZi = zz_units_sum();
                assert!(p.is_complex());
            }

            #[test]
            fn test_basic() {
                // # of full turn steps is an even natural (so a half turn is possible)
                assert!(ZZi::turn() > 0);
                assert!(ZZi::turn() % 2 == 0);

                // check vector sizes
                let roots_num = ZZi::zz_params().sym_roots_num;
                assert_eq!(ZZi::zz_params().sym_roots_sqs.len(), roots_num);
                assert_eq!(ZZi::zz_params().ccw_unit_coeffs.len(), roots_num);

                // check zero-vector
                let z = ZZi::zero();
                let cs = z.zz_coeffs();
                assert_eq!(cs.len(), roots_num);
                for i in 0..roots_num {
                    assert_eq!(cs[i], GIntT::zero());
                }

                // check one-vector
                let o = ZZi::one();
                let cs = o.zz_coeffs();
                assert_eq!(cs.len(), roots_num);
                assert_eq!(cs[0], GIntT::one());
                for i in 1..roots_num {
                    assert_eq!(cs[i], GIntT::zero());
                }
            }

            #[test]
            fn test_add_sub() {
                // test addition / subtraction
                assert_eq!(ZZi::zero() + ZZi::zero(), ZZi::zero());
                assert_eq!(ZZi::one() + ZZi::zero(), ZZi::one());
                assert_eq!(ZZi::zero() + ZZi::one(), ZZi::one());
                assert_eq!(-ZZi::one() + ZZi::one(), ZZi::zero());
                assert_eq!(ZZi::one() - ZZi::one(), ZZi::zero());
            }

            #[test]
            fn test_mul() {
                // test scalar multiplication
                assert_eq!(ZZi::zero().scale(2), ZZi::zero());
                assert_eq!(ZZi::ccw().scale(3), ZZi::ccw() + ZZi::ccw() + ZZi::ccw());
                assert_eq!(ZZi::one().scale(-1), -ZZi::one());
                assert_eq!(ZZi::one().scale(-42), ZZi::from(-42));

                // test multiplication
                assert_eq!(ZZi::zero() * ZZi::zero(), ZZi::zero());
                assert_eq!(ZZi::one() * ZZi::zero(), ZZi::zero());
                assert_eq!(ZZi::zero() * ZZi::one(), ZZi::zero());
                assert_eq!(ZZi::one() * ZZi::one(), ZZi::one());
                assert_eq!(-ZZi::one() * ZZi::one(), -ZZi::one());
                assert_eq!(ZZi::one() * (-ZZi::one()), -ZZi::one());
                assert_eq!((-ZZi::one()) * (-ZZi::one()), ZZi::one());

                // check result numerically (simplifies debugging of multiplication)
                let r = ZZi::zz_params().sym_roots_num;
                for i in 0..r {
                    println!("------------");
                    for j in 0..r {
                        println!("----");
                        println!("l={i} r={j}");
                        let mut vec1: Vec<GIntT> = vec![0.into(); r];
                        let mut vec2: Vec<GIntT> = vec![0.into(); r];
                        vec1[i] = 1.into();
                        vec2[j] = 1.into();
                        let v1 = ZZi::new(vec1.as_slice());
                        let v2 = ZZi::new(vec2.as_slice());

                        let exp = v1.complex64().re * v2.complex64().re;
                        let prod = (v1 * v2).complex64().re;
                        println!("x={}\ny={}\nx*y={}", v1, v2, v1 * v2);
                        assert!((exp - prod).abs() < 1e-8);
                    }
                }
            }

            #[test]
            fn test_rotations() {
                // test ccw()
                assert_eq!(ZZi::ccw() * ZZi::ccw().conj(), ZZi::one());
                assert_eq!(-(-(ZZi::one()) * ZZi::ccw()), ZZi::ccw());

                // test going around the unit circle step by step
                let mut x = ZZi::one();
                for _ in 0..ZZi::turn() {
                    x = x * ZZi::ccw();
                }
                assert_eq!(x, ZZi::one());

                // test unit()
                assert_eq!(ZZi::unit(0), ZZi::one());
                assert_eq!(ZZi::unit(-1), ZZi::unit(ZZi::turn() - 1));
                assert_eq!(ZZi::unit(1), ZZi::unit(ZZi::turn() + 1));
                assert_eq!(ZZi::unit(-ZZi::hturn()), ZZi::unit(ZZi::hturn()));
                assert_eq!(ZZi::unit(ZZi::hturn()), -ZZi::one());
                if ZZi::has_qturn() {
                    assert_eq!(ZZi::unit(ZZi::qturn()).zz_coeffs()[0], GIntT::from((0, 1)));
                }

                // test powi()
                assert_eq!(ZZi::ccw().pow(ZZi::hturn()), -ZZi::one());
                assert_eq!(ZZi::ccw().pow(ZZi::turn()), ZZi::one());
                assert_eq!(ZZi::ccw().pow(ZZi::hturn()).pow(2), ZZi::one());
            }

            #[test]
            #[should_panic]
            fn test_neg_powi() {
                ZZi::one().pow(-1);
            }

            #[test]
            fn test_scaling_fac() {
                // test scaling fac is correct by checking denom. of coeffs of all units
                // (that the denom. always can be expressed as multple of scaling factor)
                // and that the chosen constant factor is indeed minimal
                let sc_fac = ZZi::zz_params().scaling_fac;
                let mut max_fac: i64 = 0;
                for i in 0..ZZi::turn() {
                    let x = ZZi::unit(i);
                    println!("{x}");
                    for c in x.zz_coeffs() {
                        assert_eq!(sc_fac % c.real.denom(), 0);
                        assert_eq!(sc_fac % c.imag.denom(), 0);
                        max_fac = max_fac.max(*c.real.denom());
                        max_fac = max_fac.max(*c.imag.denom());
                    }
                }
                assert_eq!(sc_fac, max_fac);
            }

            #[test]
            fn test_xy() {
                if !ZZi::has_qturn() {
                    return;
                }

                // test correctness of splitting and reconstruction
                for a in 0..ZZi::hturn() {
                    let p = ZZi::unit(a);
                    let (x, y) = p.re_im();
                    assert_eq!(p, ZZi::from(x) + ZZi::from(y) * ZZi::unit(ZZi::qturn()));
                }
            }

            #[test]
            fn test_is_real_imag_complex() {
                assert!(ZZi::zero().is_real());
                assert!(ZZi::zero().is_imag());

                assert!(ZZi::one().is_real());
                assert!((-ZZi::one()).is_real());
                assert!(!ZZi::one().is_imag());
                assert!(!ZZi::one().is_complex());

                if ZZi::has_qturn() && ZZi::qturn() > 1 {
                    let p = ZZi::ccw();
                    assert!(!p.is_real());
                    assert!(!p.is_imag());
                    assert!(p.is_complex());
                }

                if ZZi::has_qturn() {
                    let imag_unit = ZZi::unit(ZZi::qturn());
                    assert!(!imag_unit.is_real());
                    assert!(imag_unit.is_imag());
                    assert!((-imag_unit).is_imag());
                    assert!(!imag_unit.is_complex());
                }
            }

            #[test]
            fn test_complex() {
                use num_complex::Complex64;

                let x = ZZi::zero();
                assert_eq!(x.complex64(), Complex64::zero());
                let x = ZZi::one();
                assert_eq!(x.complex64(), Complex64::one());
                let x = -ZZi::one();
                assert_eq!(x.complex64(), -Complex64::one());
                let x = ZZi::one() + ZZi::one();
                assert_eq!(x.complex64(), Complex64::new(2.0, 0.0));

                let x = ZZi::ccw();
                let c = x.complex64();
                println!("{c} = {x}");
            }

            #[test]
            fn test_hashable() {
                use std::collections::HashSet;

                let mut s: HashSet<ZZi> = HashSet::new();
                s.insert(ZZi::zero());
                s.insert(ZZi::one());
                s.insert(ZZi::ccw());
                assert!(s.contains(&ZZi::ccw()));
                assert!(s.contains(&(ZZi::ccw() + ZZi::zero())));
                assert!(!s.contains(&(ZZi::ccw() + ZZi::one())));
            }
        }
    };
}
