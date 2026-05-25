//! General traits that are needed or convenient for technical reasons to simplify implementation.

use num_integer::Integer;
use num_rational::Ratio;
use num_traits::{One, Signed, Zero};
use std::marker::Copy;
use std::ops::{Add, Div, Mul, Neg, Sub};

/// Trait with all traits a reasonable integer ring should provide.
pub trait IntRing:
    Copy
    + Clone
    + Zero
    + One
    + Neg<Output = Self>
    + Add<Self, Output = Self>
    + Sub<Self, Output = Self>
    + Mul<Self, Output = Self>
    // + Pow<i8, Output = Self>
    // + Pow<u8, Output = Self>
    + PartialEq
    + Eq
{
}

/// A field is a ring that supports division.
pub trait IntField: IntRing + Div<Self, Output = Self> {}

impl IntRing for i32 {}
impl IntField for i32 {}

impl IntRing for i64 {}
impl IntField for i64 {}

impl<T: Integer + IntRing> IntRing for Ratio<T> {}
impl<T: Integer + IntField> IntField for Ratio<T> {}

// NOTE: this only makes sense for the complex fields.
pub trait Ccw {
    /// Return unit of (counterclockwise) rotation for a complex integer ring,
    /// i.e. unit vector rotated by the smallest possible rotation step.
    ///
    /// Technically, this is the actual single generator of the cyclotomic ring,
    /// as we can recover the normal unit by taking a power of the rotation unit.
    fn ccw() -> Self;

    /// Return true if the value is equal to the rotation unit.
    fn is_ccw(&self) -> bool;
}

// NOTE: Cyclotomic rings are closed under complex conjugation, so this is ok to do.
pub trait Conj {
    /// For representations of complex numbers, return conjugated complex number.
    /// For one-dimensional numbers, return same number.
    fn conj(&self) -> Self;
}

impl<T: Integer + Clone + std::ops::Neg<Output = T>> Conj for Ratio<T> {
    fn conj(&self) -> Self {
        self.clone()
    }
}

/// Trait to be implemented by rings that support the imaginary unit i directly.
// NOTE: Only to be implemented for types that are a superset of the Gaussian integers.
pub trait OneImag {
    /// Return imaginary unit.
    fn one_i() -> Self;

    /// Return true if the value is equal to the imaginary unit.
    fn is_one_i(&self) -> bool;
}

/// Utility trait to access the integer type on top of which
/// more complex numeric types are built.
pub trait InnerIntType {
    type IntType;
}

impl InnerIntType for i32 {
    type IntType = i32;
}
impl InnerIntType for i64 {
    type IntType = i64;
}

impl<T: Integer + IntRing + InnerIntType> InnerIntType for Ratio<T> {
    type IntType = T::IntType;
}

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

/// Sign extraction for the real and imaginary parts of a (complex) value.
///
/// This trait exists so that geometric predicates can be expressed purely in
/// terms of `ZZ` arithmetic plus a per-ring sign-of-component query, without
/// the predicates themselves ever materializing a `ZZ::Real` value. That
/// allows future ring storage migrations (e.g. ZZ12 stored in a pure-i64
/// integral basis) to override sign extraction without touching the geometry
/// code that consumes it.
///
/// `re_sign` / `im_sign` must return one of `-1`, `0`, or `1`.
pub trait ReImSign {
    /// Sign of the real component: -1, 0, or 1.
    fn re_sign(&self) -> i8;
    /// Sign of the imaginary component: -1, 0, or 1.
    fn im_sign(&self) -> i8;
}

/// Default `ReImSign` implementation: any complex-valued ring with
/// `ZSigned` real components can extract signs via `re_im()` + `ZSigned`
/// predicates.
///
/// This is semantically a no-op vs the previous direct usage of `is_positive`
/// / `is_negative` / `is_zero` on the projected real parts. Per-ring
/// overrides can specialize this for storage formats that admit a cheaper
/// sign query (notably ZZ12 with its `[i64; 4]` integral-basis storage).
///
/// This is **not** a blanket impl on every cyclotomic ring: that would
/// block per-ring overrides on stable Rust (no specialization). Instead we
/// expose a macro that generates the impl for each ring type, and any type
/// that wants a custom impl simply opts out by not invoking the macro.
///
/// The macro takes the ring type plus a per-ring `coeffs_real_sign_*`
/// function from `cyclotomic::sign` whose shape is
/// `fn(&[Ratio<i64>], &[f64]) -> i8`. Real-sign feeds it the `.real` halves
/// of the GaussInt coefficient array, imag-sign the `.imag` halves.
#[macro_export]
macro_rules! impl_re_im_sign_via_proj {
    ($t:ty, $sign_func:path) => {
        impl $crate::cyclotomic::ReImSign for $t {
            fn re_sign(&self) -> i8 {
                use $crate::cyclotomic::SymNum;
                let cs = <$t as SymNum>::zz_coeffs(self);
                let n = cs.len();
                // Stack-friendly: ZZ rings have at most 8 symbolic roots.
                let mut buf: [num_rational::Ratio<i64>; 8] =
                    [num_rational::Ratio::<i64>::new_raw(0, 1); 8];
                for (i, c) in cs.iter().enumerate() {
                    buf[i] = c.real;
                }
                let roots = <$t as SymNum>::zz_params().sym_roots_sqs;
                $sign_func(&buf[..n], roots)
            }
            fn im_sign(&self) -> i8 {
                use $crate::cyclotomic::SymNum;
                let cs = <$t as SymNum>::zz_coeffs(self);
                let n = cs.len();
                let mut buf: [num_rational::Ratio<i64>; 8] =
                    [num_rational::Ratio::<i64>::new_raw(0, 1); 8];
                for (i, c) in cs.iter().enumerate() {
                    buf[i] = c.imag;
                }
                let roots = <$t as SymNum>::zz_params().sym_roots_sqs;
                $sign_func(&buf[..n], roots)
            }
        }
    };
}

/// Per-ring specialization hook for the unit-length segment intersection test.
///
/// `intersect_unit_segments(s1, s2)` is semantically equivalent to the generic
/// `cyclotomic::geometry::intersect(s1, s2)`, but **assumes** both segments
/// have unit length (i.e. `s1.1 - s1.0` and `s2.1 - s2.0` are unit vectors of
/// the ring's CCW unit group).
///
/// The hot-path implementation lives in ZZ12 and uses a 3-multiplication
/// pure-i64 fast path that exploits the unit-length structure. Every other
/// ring defers to the generic `intersect` via `impl_intersect_unit_segments_via_general!`.
pub trait IntersectUnitSegments: Sized {
    /// Return whether the unit-length segments `s1` and `s2` intersect.
    /// Touching only in endpoints does **not** count as intersection.
    fn intersect_unit_segments(s1: &(Self, Self), s2: &(Self, Self)) -> bool;
}

/// Default `IntersectUnitSegments` impl: just call the generic `intersect`.
///
/// Per-ring overrides (notably ZZ12) opt out by not invoking this macro and
/// providing a specialized impl directly.
#[macro_export]
macro_rules! impl_intersect_unit_segments_via_general {
    ($t:ty) => {
        impl $crate::cyclotomic::IntersectUnitSegments for $t {
            #[inline]
            fn intersect_unit_segments(
                s1: &($t, $t),
                s2: &($t, $t),
            ) -> bool {
                $crate::cyclotomic::geometry::intersect::<$t>(s1, s2)
            }
        }
    };
}

/// Predicate "this point is within Euclidean radius `r` of the origin".
///
/// Equivalent to `|self|^2 <= r * r`, with `r >= 0`. Exists so that the DFS
/// reachability heuristic in `rat_enum` does not have to materialize a
/// real-subring intermediate value -- per-ring overrides can stay in
/// pure-i64 arithmetic.
///
/// Per-ring overrides should compute the squared norm in pure i64 arithmetic
/// where possible. The default-via-complex64 macro is provided for rings
/// that don't yet have a direct fast path.
pub trait WithinRadius {
    /// `true` iff `|self|^2 <= radius * radius`. `radius` must be `>= 0`.
    fn within_radius(&self, radius: i64) -> bool;
}

/// Default `WithinRadius` impl via `complex64()` magnitude with a small
/// epsilon to tolerate float rounding around the boundary.
///
/// This is only the slow fallback path for rings without a pure-i64
/// specialization (currently every ring except ZZ12 uses this). The DFS
/// pruning in `rat_enum` is monotone -- false-positives only widen the
/// search, never compromise correctness -- so a small epsilon bias toward
/// "inside" is safe.
#[macro_export]
macro_rules! impl_within_radius_via_norm_sq {
    ($t:ty) => {
        impl $crate::cyclotomic::WithinRadius for $t {
            fn within_radius(&self, radius: i64) -> bool {
                use $crate::cyclotomic::SymNum;
                let c = <$t as SymNum>::complex64(self);
                let norm_sq = c.re * c.re + c.im * c.im;
                let r_sq = (radius as f64) * (radius as f64);
                // 1e-9 tolerance: well below the resolution of any
                // representable lattice point in any ring of practical
                // interest, but enough to absorb f64 rounding for very
                // dense rings like ZZ60.
                norm_sq <= r_sq + 1e-9
            }
        }
    };
}
