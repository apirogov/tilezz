use super::gaussint::GaussInt;
use super::traits::IntRing;
use num_integer::Integer;
use num_rational::Ratio;
use num_complex::Complex64;
use num_traits::{One, PrimInt, Zero, ToPrimitive};
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
    T: PrimInt + Integer + IntRing,
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

pub trait ZZBase<T: IntRing + ToPrimitive + 'static>  { 
    #[inline]
    fn turn() -> i8 {
        Self::zz_params().full_turn_steps
    }
    #[inline]
    fn hturn() -> i8 {
        Self::turn()/2
    }

    fn complex(&self) -> Complex64 {
        let nums: Vec<Complex64> = self.zz_coeffs().into_iter().map(|x| {
            let re = x.real.to_f64().unwrap();
            let im = x.imag.to_f64().unwrap();
            Complex64::new(re, im)
        }).collect();
        let units: Vec<Complex64> = Self::zz_params().sym_roots_sqs.into_iter().map(|x| {
            Complex64::new(x.sqrt(), 0.0)
        }).collect();
        let mut ret = Complex64::zero();
        for (n, u) in nums.iter().zip(units.iter()) {
            ret += n * u;
        }
        ret
    }

    // functions that can be implemented via zz_base_impl!
    fn new(coeffs: &[GaussInt<T>]) -> Self;
    fn unit(i: i8) -> Self;
    fn powi(&self, i: i8) -> Self;

    fn zz_coeffs(&self) -> &[GaussInt<T>];
    fn zz_coeffs_mut(&mut self) -> &mut [GaussInt<T>];

    fn zz_params() -> &'static ZZParams<'static, T>;
    fn zz_ccw_vec() -> Vec<GaussInt<T>>;
    fn zz_mul_arrays(x: &[GaussInt<T>], y: &[GaussInt<T>]) -> Vec<GaussInt<T>>;

    // implementations for implementing other traits using zz_ops_impl!
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

    #[inline]
    fn conj(&self) -> Self
    where
        Self: Sized,
    {
        let cs: Vec<GaussInt<T>> = self.zz_coeffs().iter().map(|c| c.conj()).collect();
        Self::new(&cs)
    }
}

#[macro_export]
macro_rules! zz_base_impl {
    ($name:ident, $params:ident, $mul_func:ident) => {
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Hash)]
        pub struct $name {
            coeffs: [GInt; $params.sym_roots_num],
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
                return $mul_func(x, y);
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
                if (i<0) {
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
