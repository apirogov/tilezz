//! Manual ZZ12 implementation using a pure-`i64` integral cyclotomic basis.
//!
//! Storage is `[i64; 4]` representing `a + b*zeta + c*zeta^2 + d*zeta^3`
//! where `zeta = e^(i*pi/6) = sqrt(3)/2 + i/2`. The minimal polynomial of
//! `zeta` over `Z` is `Phi_12(x) = x^4 - x^2 + 1`, so
//! `zeta^4 = zeta^2 - 1`, `zeta^5 = zeta^3 - zeta`, `zeta^6 = -1`, etc.
//!
//! This replaces the macro-generated `ZZ12` that stored
//! `[GaussInt<Ratio<i64>>; 2]` over the symbolic `{1, sqrt(3)}` basis.
//! Eliminating `num_rational::Ratio` from the multiplication path is a
//! significant perf win (rat_enum benchmark).
//!
//! Cartesian projection:
//!
//! ```text
//!   Re(zeta^k):    k=0 -> 1,           k=1 -> sqrt(3)/2,   k=2 -> 1/2,        k=3 -> 0
//!   Im(zeta^k):    k=0 -> 0,           k=1 -> 1/2,         k=2 -> sqrt(3)/2,  k=3 -> 1
//! ```
//!
//! So for `(a, b, c, d)`:
//!
//! ```text
//!   Re = (2a + c)/2 + (b/2) sqrt(3)
//!   Im = (b + 2d)/2 + (c/2) sqrt(3)
//! ```
//!
//! Equivalently, `Z12` (with scaling_fac=2, basis `{sqrt(1), sqrt(3)}`)
//! storing numerators-over-2 `(p, q)` denotes `p/2 + (q/2) sqrt(3)`. The
//! re/im conversions therefore feed Z12 raw numerators `(2a+c, b)` and
//! `(b+2d, c)` respectively (each with implicit denominator 2 contributed
//! by Z12's scaling factor).

use std::fmt;
use std::fmt::Display;
use std::hash::Hash;
use std::ops::{Add, Mul, Neg, Sub};

use num_complex::Complex64;
use num_rational::Ratio;
use num_traits::{One, Pow, Zero};

use super::numtraits::{Ccw, Conj, InnerIntType, IntRing, ReImSign, WithinRadius};
use super::params::ZZ12_PARAMS;
use super::symnum::{SymNum, ZZComplex, ZZParams};
use super::traits::{
    ComplexTraits, HasZZ12Impl, HasZZ4Impl, HasZZ6Impl, IsComplex, IsReal, IsRing, IsRingOrField,
    RingTraits, ZType, ZZType,
};
use super::types::Z12;
use super::units::Units;

/// ZZ12 stored in the integral cyclotomic basis `{1, zeta, zeta^2, zeta^3}`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct ZZ12 {
    /// `coeffs[0] + coeffs[1]*zeta + coeffs[2]*zeta^2 + coeffs[3]*zeta^3`
    coeffs: [i64; 4],
}

impl ZZ12 {
    #[inline]
    pub const fn from_int_coeffs(coeffs: [i64; 4]) -> Self {
        Self { coeffs }
    }

    /// Raw access to integral-basis coefficients (for debugging / tests).
    #[inline]
    pub const fn int_coeffs(&self) -> [i64; 4] {
        self.coeffs
    }
}

// ----------------
// SymNum

impl SymNum for ZZ12 {
    type Scalar = i64;

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
        // We reuse the existing `ZZ12_PARAMS` only for its `full_turn_steps = 12`
        // value, which is the only field consumers outside this module read
        // (e.g. `Self::turn()`). The `sym_roots_*` and `ccw_unit_coeffs` fields
        // are NOT meaningful in the integral-basis storage and must not be
        // used to interpret `zz_coeffs()` of this type.
        &ZZ12_PARAMS
    }

    #[inline]
    fn zz_mul_arrays(_x: &[Self::Scalar], _y: &[Self::Scalar]) -> Vec<Self::Scalar> {
        // We provide a direct `Mul<ZZ12>` impl that operates on the fixed-size
        // `[i64; 4]` storage. This vec-shaped trampoline is only required by
        // the `SymNum` trait surface; nothing on the hot path reaches it.
        unreachable!(
            "ZZ12 multiplication goes through the direct `Mul` impl, not `zz_mul_arrays`"
        )
    }

    #[inline]
    fn zz_mul_scalar(x: &[Self::Scalar], scalar: i64) -> Vec<Self::Scalar> {
        x.iter().map(|c| *c * scalar).collect()
    }

    fn new(coeffs: &[Self::Scalar]) -> Self {
        assert_eq!(coeffs.len(), 4, "ZZ12::new expects 4 integral-basis coefficients");
        Self {
            coeffs: [coeffs[0], coeffs[1], coeffs[2], coeffs[3]],
        }
    }

    fn complex64(&self) -> Complex64 {
        // zeta = sqrt(3)/2 + i/2
        // Re = a + b * sqrt(3)/2 + c * 1/2 + d * 0
        // Im = 0 + b * 1/2       + c * sqrt(3)/2 + d
        const HALF_SQRT_3: f64 = 0.866_025_403_784_438_6_f64; // sqrt(3)/2
        let [a, b, c, d] = self.coeffs;
        let (a, b, c, d) = (a as f64, b as f64, c as f64, d as f64);
        let re = a + b * HALF_SQRT_3 + 0.5 * c;
        let im = 0.5 * b + c * HALF_SQRT_3 + d;
        Complex64::new(re, im)
    }
}

// ----------------
// Display

impl Display for ZZ12 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Display via the projected Z12 representation so output matches
        // the previous symbolic-basis storage's format.
        let z: Zz12DisplayBridge = (*self).into();
        z.fmt(f)
    }
}

// Bridge type that exists purely so we can reuse Z12's existing `Display`
// implementation (which prints `(p/2) + (q/2)*sqrt(3)` style). For ZZ12 we
// show `Re + i*Im` using the existing GaussInt<Z12> display logic.
struct Zz12DisplayBridge {
    re: Z12,
    im: Z12,
}

impl From<ZZ12> for Zz12DisplayBridge {
    fn from(z: ZZ12) -> Self {
        let (re, im) = z.re_im();
        Self { re, im }
    }
}

impl Display for Zz12DisplayBridge {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Emulate `GaussInt<Z12>` display so existing tests keep matching.
        let real_zero = self.re.is_zero();
        let imag_zero = self.im.is_zero();
        if real_zero && imag_zero {
            return write!(f, "0");
        }
        if imag_zero {
            return write!(f, "{}", self.re);
        }
        if real_zero {
            return write!(f, "{}i", self.im);
        }
        // both nonzero
        write!(f, "{}+{}i", self.re, self.im)
    }
}

// ----------------
// Ring ops (Neg / Add / Sub / Mul)

impl Neg for ZZ12 {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Self {
            coeffs: [
                -self.coeffs[0],
                -self.coeffs[1],
                -self.coeffs[2],
                -self.coeffs[3],
            ],
        }
    }
}

impl Add<ZZ12> for ZZ12 {
    type Output = Self;
    #[inline]
    fn add(self, other: Self) -> Self {
        Self {
            coeffs: [
                self.coeffs[0] + other.coeffs[0],
                self.coeffs[1] + other.coeffs[1],
                self.coeffs[2] + other.coeffs[2],
                self.coeffs[3] + other.coeffs[3],
            ],
        }
    }
}

impl Sub<ZZ12> for ZZ12 {
    type Output = Self;
    #[inline]
    fn sub(self, other: Self) -> Self {
        Self {
            coeffs: [
                self.coeffs[0] - other.coeffs[0],
                self.coeffs[1] - other.coeffs[1],
                self.coeffs[2] - other.coeffs[2],
                self.coeffs[3] - other.coeffs[3],
            ],
        }
    }
}

impl Mul<ZZ12> for ZZ12 {
    type Output = Self;
    /// Multiplication in the integral cyclotomic basis `{1, zeta, zeta^2, zeta^3}`.
    ///
    /// Derivation: expand `x * y` polynomially in `zeta`, then reduce powers
    /// `zeta^4 = zeta^2 - 1`, `zeta^5 = zeta^3 - zeta`, `zeta^6 = -1`.
    ///
    /// For `x = (a, b, c, d)` and `y = (e, f, g, h)`:
    ///
    /// ```text
    ///   result_0 = ae - bh - cg - df - dh
    ///   result_1 = af + be - ch - dg
    ///   result_2 = ag + bf + ce + bh + cg + df
    ///   result_3 = ah + bg + cf + de + ch + dg
    /// ```
    #[inline]
    fn mul(self, other: Self) -> Self {
        let [a, b, c, d] = self.coeffs;
        let [e, f, g, h] = other.coeffs;
        let r0 = a * e - b * h - c * g - d * f - d * h;
        let r1 = a * f + b * e - c * h - d * g;
        let r2 = a * g + b * f + c * e + b * h + c * g + d * f;
        let r3 = a * h + b * g + c * f + d * e + c * h + d * g;
        Self {
            coeffs: [r0, r1, r2, r3],
        }
    }
}

impl Pow<u8> for ZZ12 {
    type Output = Self;
    fn pow(self, other: u8) -> Self {
        self.zz_pow(other)
    }
}
impl Pow<i8> for ZZ12 {
    type Output = Self;
    fn pow(self, other: i8) -> Self {
        assert!(other >= 0, "Negative powers are not supported!");
        self.zz_pow(other as u8)
    }
}

// ----------------
// Zero / One

impl Zero for ZZ12 {
    #[inline]
    fn zero() -> Self {
        Self { coeffs: [0; 4] }
    }
    #[inline]
    fn is_zero(&self) -> bool {
        self.coeffs == [0; 4]
    }
}

impl One for ZZ12 {
    #[inline]
    fn one() -> Self {
        Self {
            coeffs: [1, 0, 0, 0],
        }
    }
    #[inline]
    fn is_one(&self) -> bool {
        self.coeffs == [1, 0, 0, 0]
    }
}

impl IntRing for ZZ12 {}

// ----------------
// Conjugation

impl Conj for ZZ12 {
    /// Conjugation in `{1, zeta, ..., zeta^3}`:
    ///
    /// `conj(zeta) = zeta^11 = zeta - zeta^3`  ->  `(0,  1, 0, -1)`
    /// `conj(zeta^2) = zeta^10 = 1 - zeta^2`   ->  `(1,  0, -1, 0)`
    /// `conj(zeta^3) = zeta^9 = -zeta^3`        ->  `(0, 0, 0, -1)`
    ///
    /// So `conj((a, b, c, d)) = (a + c, b, -c, -b - d)`.
    #[inline]
    fn conj(&self) -> Self {
        let [a, b, c, d] = self.coeffs;
        Self {
            coeffs: [a + c, b, -c, -b - d],
        }
    }
}

// ----------------
// InnerIntType + From<i64>

impl InnerIntType for ZZ12 {
    type IntType = i64;
}

impl From<i64> for ZZ12 {
    #[inline]
    fn from(value: i64) -> Self {
        Self {
            coeffs: [value, 0, 0, 0],
        }
    }
}

impl From<(i64, i64)> for ZZ12 {
    #[inline]
    fn from((re, im): (i64, i64)) -> Self {
        // a + b*i, where i = zeta^3, so im contributes to coeffs[3].
        Self {
            coeffs: [re, 0, 0, im],
        }
    }
}

// Conversion: lift Z12 (real, basis {sqrt(1), sqrt(3)}, scaling_fac=2) into ZZ12.
//
// Z12 stores numerators-over-2: `(p, q)` -> `p/2 + (q/2)*sqrt(3)`.
//
// In ZZ12's integral basis with `zeta^2 = 1/2 + sqrt(3)/2 * i ...` actually
// `sqrt(3) = 2*Re(zeta) = zeta + conj(zeta) = zeta + (zeta - zeta^3) = 2*zeta - zeta^3`.
// And `1 = (1, 0, 0, 0)`.
//
// So lifting `p/2 + (q/2)*sqrt(3)`:
//   (p/2) * (1, 0, 0, 0) + (q/2) * (0, 2, 0, -1) * 1/?
//
// Hmm, this requires fractions. Let's instead use the formula that `2 * sqrt(3) = 2*(2*zeta - zeta^3) = 4*zeta - 2*zeta^3`,
// so `(q/2)*sqrt(3) = q*zeta - (q/2)*zeta^3` ... still fractions for odd q.
//
// Actually `sqrt(3) = 2*Re(zeta)`. But we want integral ZZ12 values, so we
// need `p` and `q` to satisfy divisibility. For the cases that arise in the
// codebase (sums of unit vectors, norms, dots) the Z12 values DO end up
// integral in ZZ12, but in general a Z12 can be a half-integer combination.
//
// To handle all valid Z12 inputs, we exploit that Z12 stores `(p, q)` as
// `Ratio<i64>` with potentially-non-trivial denominators. The safe path is
// to read the rational coefficients directly.
impl From<Z12> for ZZ12 {
    fn from(value: Z12) -> Self {
        // Z12.coeffs = [Ratio<i64>; 2] over basis [sqrt(1) = 1, sqrt(3)].
        // value = c0 + c1 * sqrt(3), where c0, c1 are Ratio<i64>.
        // Express in integral basis using sqrt(3) = 2*zeta - zeta^3:
        //   value = c0 * (1,0,0,0) + c1 * (0, 2, 0, -1)
        //         = (c0, 2*c1, 0, -c1).
        //
        // For Z12 values arising from re_im() of an integral ZZ12, all
        // coefficients are integers (the Z12 denominators are 1 or 2 but
        // those halves disappear because Z12 stores `(2a+c, b)` and
        // `(b+2d, c)` over its scaling_fac=2; once interpreted as actual
        // rationals they are integers because `c0 = (2a+c)/2` may be a
        // half-integer when `c` is odd, but then `c1 = b/2` is paired such
        // that the result is integral in ZZ12).
        //
        // To handle the general case we round to integers, panicking if the
        // resulting ZZ12 would have non-integral coefficients (which is a
        // logic bug rather than a runtime error).
        let c0: Ratio<i64> = value.zz_coeffs()[0];
        let c1: Ratio<i64> = value.zz_coeffs()[1];
        // (c0, 2*c1, 0, -c1) must be integral.
        let to_int = |r: Ratio<i64>| -> i64 {
            assert!(
                r.denom() == &1,
                "Z12 -> ZZ12 conversion produced a non-integral coefficient: {r:?}"
            );
            *r.numer()
        };
        let two = Ratio::<i64>::from_integer(2);
        Self {
            coeffs: [to_int(c0), to_int(two * c1), 0, to_int(-c1)],
        }
    }
}

// ----------------
// RingTraits (auto from InnerIntType + IntRing + Conj + From<IntT> + From<IntType>)

impl RingTraits for ZZ12 {}

// ----------------
// Ccw + ZZComplex

impl Ccw for ZZ12 {
    #[inline]
    fn ccw() -> Self {
        // zeta = (0, 1, 0, 0)
        Self {
            coeffs: [0, 1, 0, 0],
        }
    }
    #[inline]
    fn is_ccw(&self) -> bool {
        self.coeffs == [0, 1, 0, 0]
    }
}

impl ZZComplex for ZZ12 {
    fn is_real(&self) -> bool {
        // Im((a,b,c,d)) = (b + 2d) + c*sqrt(3) all halved.
        // For Im == 0 we need (b + 2d) == 0 AND c == 0 (sqrt(3) is irrational).
        let [_, b, c, d] = self.coeffs;
        c == 0 && (b + 2 * d) == 0
    }

    fn is_imag(&self) -> bool {
        // Re((a,b,c,d)) = (2a + c) + b*sqrt(3) all halved.
        // For Re == 0 we need (2a + c) == 0 AND b == 0.
        let [a, b, c, _] = self.coeffs;
        b == 0 && (2 * a + c) == 0
    }

    fn re(&self) -> Z12 {
        // Re = (2a + c)/2 + (b/2) * sqrt(3)
        // Z12's storage basis is [sqrt(1), sqrt(3)] with `Ratio<i64>` coeffs
        // (no per-ring scaling factor applied at this layer -- the
        // `scaling_fac=2` only governs how unit vectors are constructed).
        let [a, b, c, _] = self.coeffs;
        let half = Ratio::<i64>::new_raw(1, 2);
        let c0 = Ratio::<i64>::from_integer(2 * a + c) * half;
        let c1 = Ratio::<i64>::from_integer(b) * half;
        Z12::new(&[c0, c1])
    }

    fn im(&self) -> Z12 {
        // Im = (b + 2d)/2 + (c/2) * sqrt(3)
        let [_, b, c, d] = self.coeffs;
        let half = Ratio::<i64>::new_raw(1, 2);
        let c0 = Ratio::<i64>::from_integer(b + 2 * d) * half;
        let c1 = Ratio::<i64>::from_integer(c) * half;
        Z12::new(&[c0, c1])
    }
}

impl ComplexTraits for ZZ12 {}

// ----------------
// IsRingOrField + IsRing + IsComplex + ZZType

impl IsRingOrField for ZZ12 {
    type Real = Z12;
    type Complex = ZZ12;
}
impl IsRingOrField for Z12 {
    type Real = Z12;
    type Complex = ZZ12;
}
impl IsRing for ZZ12 {
    type Real = Z12;
    type Complex = ZZ12;
}
impl IsRing for Z12 {
    type Real = Z12;
    type Complex = ZZ12;
}
impl IsReal for Z12 {
    type Ring = Z12;
}
impl IsComplex for ZZ12 {
    type Ring = ZZ12;
}
impl ZType for Z12 {
    type Complex = ZZ12;
}
impl ZZType for ZZ12 {
    type Real = Z12;
}

// ----------------
// HasZZ4Impl, HasZZ6Impl, HasZZ12Impl

impl HasZZ4Impl for ZZ12 {}
impl HasZZ6Impl for ZZ12 {}
impl HasZZ12Impl for ZZ12 {}

// ----------------
// Units::unit

impl Units for ZZ12 {
    #[inline]
    fn unit(angle: i8) -> Self {
        // zeta^k lookup for k in 0..12. Negative k normalizes via mod 12.
        // Inline because the table is tiny and avoids static + OnceLock cost.
        static UNIT_TABLE: [[i64; 4]; 12] = [
            [1, 0, 0, 0],   // 1
            [0, 1, 0, 0],   // zeta
            [0, 0, 1, 0],   // zeta^2
            [0, 0, 0, 1],   // zeta^3 = i
            [-1, 0, 1, 0],  // zeta^4 = zeta^2 - 1
            [0, -1, 0, 1],  // zeta^5 = zeta^3 - zeta
            [-1, 0, 0, 0],  // zeta^6 = -1
            [0, -1, 0, 0],  // zeta^7
            [0, 0, -1, 0],  // zeta^8
            [0, 0, 0, -1],  // zeta^9
            [1, 0, -1, 0],  // zeta^10
            [0, 1, 0, -1],  // zeta^11
        ];
        let idx = angle.rem_euclid(12) as usize;
        Self {
            coeffs: UNIT_TABLE[idx],
        }
    }
}

// ----------------
// ReImSign override: this is the perf-critical hot path.

impl ReImSign for ZZ12 {
    /// Sign of `Re(z) = ((2a + c) + b*sqrt(3)) / 2`.
    ///
    /// Denominator 2 is positive, so this is `sign((2a + c) + b*sqrt(3))`.
    #[inline]
    fn re_sign(&self) -> i8 {
        let [a, b, c, _] = self.coeffs;
        sign_m_plus_n_sqrt3(2 * a + c, b)
    }

    /// Sign of `Im(z) = ((b + 2d) + c*sqrt(3)) / 2`.
    #[inline]
    fn im_sign(&self) -> i8 {
        let [_, b, c, d] = self.coeffs;
        sign_m_plus_n_sqrt3(b + 2 * d, c)
    }
}

/// Sign of `m + n*sqrt(3)` where m, n are i64.
///
/// Returns -1, 0, or 1.
///
/// Strategy: cheap cases (one zero, both same sign) resolve without
/// multiplications; mixed-sign case compares `m^2` vs `3 * n^2`.
#[inline]
fn sign_m_plus_n_sqrt3(m: i64, n: i64) -> i8 {
    if m == 0 && n == 0 {
        return 0;
    }
    if m >= 0 && n >= 0 {
        return 1;
    }
    if m <= 0 && n <= 0 {
        return -1;
    }
    // Mixed signs.
    let m_sq = m
        .checked_mul(m)
        .expect("ZZ12 re_sign/im_sign overflow: m*m");
    let n_sq = n
        .checked_mul(n)
        .expect("ZZ12 re_sign/im_sign overflow: n*n");
    let three_n_sq = 3i64
        .checked_mul(n_sq)
        .expect("ZZ12 re_sign/im_sign overflow: 3*n*n");
    if m > 0 {
        // n < 0: positive iff m^2 > 3*n^2
        match m_sq.cmp(&three_n_sq) {
            std::cmp::Ordering::Less => -1,
            std::cmp::Ordering::Equal => 0,
            std::cmp::Ordering::Greater => 1,
        }
    } else {
        // m < 0, n > 0: positive iff 3*n^2 > m^2
        match three_n_sq.cmp(&m_sq) {
            std::cmp::Ordering::Less => -1,
            std::cmp::Ordering::Equal => 0,
            std::cmp::Ordering::Greater => 1,
        }
    }
}

// ----------------
// WithinRadius override: pure-i64 squared-norm comparison.
//
// For z = (a, b, c, d) in the integral basis {1, zeta, zeta^2, zeta^3}:
//   Re(z) = ((2a + c) + b*sqrt(3)) / 2
//   Im(z) = ((b + 2d) + c*sqrt(3)) / 2
//   |z|^2 = Re(z)^2 + Im(z)^2 = (M + N*sqrt(3)) / 4
// where
//   M = (2a + c)^2 + 3*b^2 + (b + 2d)^2 + 3*c^2
//   N = 2 * ((2a + c)*b + (b + 2d)*c)
//
// |z|^2 <= radius^2 iff (M + N*sqrt(3)) / 4 <= radius^2
//                  iff (M - 4*radius^2) + N*sqrt(3) <= 0
// which is just `sign_m_plus_n_sqrt3(M - 4 * radius^2, N) <= 0`.

impl WithinRadius for ZZ12 {
    #[inline]
    fn within_radius(&self, radius: i64) -> bool {
        let [a, b, c, d] = self.coeffs;
        let s = 2 * a + c;
        let t = b + 2 * d;
        let m = s * s + 3 * b * b + t * t + 3 * c * c;
        let n = 2 * (s * b + t * c);
        let four_r_sq = 4 * radius * radius;
        sign_m_plus_n_sqrt3(m - four_r_sq, n) <= 0
    }
}

// ----------------
// OneImag (via HasZZ4 + Units blanket impl in `traits.rs`). Nothing to do.

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic::constants::zz_units_sum;
    use crate::cyclotomic::numtraits::OneImag;
    use crate::cyclotomic::symnum::SymNum;

    type ZZi = ZZ12;

    #[test]
    fn test_units_periodic() {
        for k in 0i8..12 {
            assert_eq!(<ZZi as Units>::unit(k), <ZZi as Units>::unit(k + 12));
            assert_eq!(<ZZi as Units>::unit(k), <ZZi as Units>::unit(k - 12));
        }
    }

    #[test]
    fn test_units_group() {
        for a in 0..12 {
            for b in 0..12 {
                let lhs = <ZZi as Units>::unit(a) * <ZZi as Units>::unit(b);
                let rhs = <ZZi as Units>::unit(((a as i16 + b as i16) % 12) as i8);
                assert_eq!(lhs, rhs, "unit({a})*unit({b})");
            }
        }
    }

    #[test]
    fn test_zeta_squared() {
        let z = ZZi::ccw();
        assert_eq!(z * z, <ZZi as Units>::unit(2));
    }

    #[test]
    fn test_zeta_cubed_is_i() {
        let z = ZZi::ccw();
        let i = z * z * z;
        assert_eq!(i, ZZ12::one_i());
        assert_eq!(i * i, -ZZ12::one());
    }

    #[test]
    fn test_complex64_roundtrip() {
        for k in 0..12 {
            let u = <ZZi as Units>::unit(k);
            let c = u.complex64();
            let exp_angle = (k as f64) * std::f64::consts::PI / 6.0;
            assert!((c.re - exp_angle.cos()).abs() < 1e-12);
            assert!((c.im - exp_angle.sin()).abs() < 1e-12);
        }
    }

    #[test]
    fn test_conj_round_trip() {
        for k in 0..12 {
            let u = <ZZi as Units>::unit(k);
            assert_eq!(u.conj().conj(), u);
            // |u|^2 == 1 -> u * conj(u) == 1
            assert_eq!(u * u.conj(), ZZi::one(), "k={k}");
        }
    }

    #[test]
    fn test_re_im_signs() {
        for k in 0..12 {
            let u = <ZZi as Units>::unit(k);
            let c = u.complex64();
            let exp_re_sign = if c.re > 1e-9 {
                1
            } else if c.re < -1e-9 {
                -1
            } else {
                0
            };
            let exp_im_sign = if c.im > 1e-9 {
                1
            } else if c.im < -1e-9 {
                -1
            } else {
                0
            };
            assert_eq!(u.re_sign(), exp_re_sign, "k={k}");
            assert_eq!(u.im_sign(), exp_im_sign, "k={k}");
        }
    }

    #[test]
    fn test_basic_arith() {
        assert_eq!(ZZi::zero() + ZZi::zero(), ZZi::zero());
        assert_eq!(ZZi::one() + ZZi::zero(), ZZi::one());
        assert_eq!(ZZi::one() * ZZi::one(), ZZi::one());
        assert_eq!(-ZZi::one() * -ZZi::one(), ZZi::one());
        assert_eq!(ZZi::ccw().pow(12u8), ZZi::one());
        assert_eq!(ZZi::ccw().pow(6u8), -ZZi::one());
    }

    #[test]
    fn test_units_sum_complex() {
        let p: ZZi = zz_units_sum();
        assert!(p.is_complex());
    }
}
