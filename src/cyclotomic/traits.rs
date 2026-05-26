//! Cyclotomic-ring trait surface.
//!
//! Every `ZZn` ring (`n in {4, 8, 10, 12, 16, 20, 24, 32, 60}`) is generated
//! by `define_integral_zz!` in `integral_basis.rs` and lives in `rings.rs`.
//! All trait declarations that those rings implement -- and the marker
//! traits downstream code uses to constrain "is this a ring with feature
//! X" -- collect in this file.
//!
//! The split with `numtraits.rs` is "pure-numeric foundation lives there,
//! everything specific to a complex-valued cyclotomic ring lives here".

use std::fmt::{Debug, Display};
use std::hash::Hash;
use std::ops::Neg;

use num_complex::Complex64;
use num_integer::Integer;
use num_rational::Ratio;
use num_traits::{FromPrimitive, One, ToPrimitive};

use super::numtraits::{InnerIntType, IntField, IntRing, ZSigned};

// ----------------
// Type aliases.

/// The integer scalar used throughout: every ring stores `[i64; phi(n)]`.
pub type IntT = i64;

// ----------------
// Per-element abstract operations.

/// "Counterclockwise unit": the primitive generator of the ring's CCW
/// rotation group, i.e. `e^(2*pi*i/n)` for `ZZn`.
pub trait Ccw {
    fn ccw() -> Self;
    fn is_ccw(&self) -> bool;
}

/// Complex conjugation. For 1-dim values (`i32`, `i64`, `Ratio<T>`) this
/// is the identity.
pub trait Conj {
    fn conj(&self) -> Self;
}

impl<T: Integer + Clone + Neg<Output = T>> Conj for Ratio<T> {
    fn conj(&self) -> Self {
        self.clone()
    }
}

/// Rings that contain the imaginary unit `i` (equivalently `4 | n`).
pub trait OneImag {
    fn one_i() -> Self;
    fn is_one_i(&self) -> bool;
}

/// Sign extraction for the real and imaginary parts of a (complex) value.
/// Each helper returns `-1`, `0`, or `1`.
///
/// Used by geometric predicates (`geometry::intersect`) to avoid
/// materializing a real-subring intermediate value; each ring routes the
/// query through its exact-sign helper from `sign.rs`.
pub trait ReImSign {
    fn re_sign(&self) -> i8;
    fn im_sign(&self) -> i8;
}

/// Per-ring specialization hook for the unit-length segment intersection
/// test.
///
/// `intersect_unit_segments(s1, s2)` is equivalent to
/// `geometry::intersect(s1, s2)` *given* that both segments have unit
/// length. ZZ12 uses a hand-rolled 3-multiplication pure-i64 fast path;
/// every other ring routes through `intersect_unit_segments_basis` via
/// the `impl_integral_intersect_unit_segments_via_basis!` macro.
pub trait IntersectUnitSegments: Sized {
    fn intersect_unit_segments(s1: &(Self, Self), s2: &(Self, Self)) -> bool;
}

/// Predicate "this point is within Euclidean radius `r` of the origin",
/// i.e. `|self|^2 <= r * r` with `r >= 0`.
///
/// Used by the DFS reachability pruning in `rat_enum`; each ring routes
/// through `impl_integral_within_radius_via_complex64!` (f64 magnitude
/// with a tiny epsilon -- safe because DFS pruning is monotone in
/// inclusions) or a hand-rolled pure-i64 override (ZZ12).
pub trait WithinRadius {
    fn within_radius(&self, radius: i64) -> bool;
}

// ----------------
// Top-level shared trait every ring implements.

/// Bound for the scalar coefficient type of a `SymNum` ring. In the
/// integer-basis storage `Self::Scalar = i64` for all rings, but kept as
/// an associated type so the trait can absorb a future field-of-fractions
/// variant.
pub trait SymScalar: IntField + FromPrimitive + ZSigned + Debug + Display {}
impl<T: IntField + FromPrimitive + ZSigned + Debug + Display> SymScalar for T {}

/// Shared cyclotomic-ring trait. The `define_integral_zz!` macro emits
/// the per-ring impl.
pub trait SymNum: Clone + Copy + PartialEq + Eq + Hash + Debug + Display {
    /// Coefficient type of the ring's storage vector.
    type Scalar: SymScalar;

    /// Number of distinct angles in one full turn (`n` for `ZZn`).
    /// Implemented per-ring as a literal.
    fn turn() -> i8;

    /// Half-turn angle: `unit(hturn()) = -1`.
    #[inline]
    fn hturn() -> i8 {
        Self::turn() / 2
    }

    /// `true` if the ring contains `i` (equivalently `4 | n`).
    #[inline]
    fn has_qturn() -> bool {
        Self::turn() % 4 == 0
    }

    /// Quarter-turn angle. Only meaningful when `has_qturn()`.
    #[inline]
    fn qturn() -> i8 {
        Self::turn() / 4
    }

    /// `Some(qturn)` when the ring contains `i`, `None` otherwise.
    #[inline]
    fn opt_qturn() -> Option<i8> {
        if Self::has_qturn() {
            Some(Self::qturn())
        } else {
            None
        }
    }

    /// Scalar multiplication by an integer.
    #[inline]
    fn scale<I: Integer + ToPrimitive>(&self, scalar: I) -> Self
    where
        Self: Sized + Copy,
    {
        let sc = Self::Scalar::from_i64(scalar.to_i64().unwrap()).unwrap();
        let mut ret = *self;
        for c in ret.zz_coeffs_mut().iter_mut() {
            *c = *c * sc;
        }
        ret
    }

    /// Project to a `Complex64`.
    fn complex64(&self) -> Complex64;

    /// `(re, im)` of `complex64()`.
    #[inline]
    fn xy(&self) -> (f64, f64) {
        let c = self.complex64();
        (c.re, c.im)
    }

    fn zz_coeffs(&self) -> &[Self::Scalar];
    fn zz_coeffs_mut(&mut self) -> &mut [Self::Scalar];

    /// Exponentiation by squaring. The macro emits `Pow<u8>` / `Pow<i8>`
    /// impls that call this.
    fn zz_pow(&self, i: u8) -> Self
    where
        Self: IsRingOrField,
    {
        let mut base = *self;
        let mut exp = i;
        let mut acc = Self::one();
        while exp != 0 {
            if (exp & 1) == 1 {
                acc = acc * base;
            }
            exp >>= 1;
            if exp != 0 {
                base = base * base;
            }
        }
        acc
    }
}

/// Predicates that classify a cyclotomic ring element as purely real,
/// purely imaginary, or strictly complex.
pub trait ZZComplex {
    fn is_real(&self) -> bool;
    fn is_imag(&self) -> bool;

    /// Default: strictly complex iff neither real nor imaginary.
    fn is_complex(&self) -> bool {
        !self.is_real() && !self.is_imag()
    }
}

// ----------------
// Composite ring traits and `Units`.

/// Sub-bag of trait obligations: the trio composed by `IsRingOrField`.
pub trait ComplexTraits: Ccw + ZZComplex {}

pub trait RingTraits:
    InnerIntType + IntRing + Conj + From<IntT> + From<<Self as InnerIntType>::IntType>
{
}

pub trait IsSymNum: SymNum + RingTraits {}
impl<T: SymNum + RingTraits> IsSymNum for T {}

/// Catch-all "this is a cyclotomic ring or field" bound. Today every
/// implementing type is a ring; a future `IsField` would also implement
/// `IsRingOrField` and add a `Div<Self, Output = Self>` bound.
pub trait IsRingOrField: IsSymNum + ComplexTraits + ReImSign + IntersectUnitSegments {}

/// Integer ring (no division). Every current `ZZn` implements this.
pub trait IsRing: IsRingOrField {}

// Future: pub trait IsField: IsRingOrField + Div<Self, Output = Self> {}

/// Marker for the standalone `Z[zeta_n]` rings (today the `ZZn` types).
pub trait ZZType: IsRing {}

/// Cached lookup of roots of unity for a concrete ring.
///
/// Implementations are emitted by `impl_integral_units_via_basis!` (a
/// `OnceLock`-cached unit table built from `derive_units_lookup`) or by a
/// hand-rolled per-ring `static UNIT_TABLE: [[i64; PHI]; N] = [...]` for
/// hot-path rings (ZZ12).
pub trait Units: Copy + One + Ccw + IsRingOrField {
    /// Unit-length vector pointing in direction of `angle`.
    fn unit(angle: i8) -> Self;

    /// Build the full unit table of size `turn()` by walking `ccw()^k`.
    /// Used by helpers that need to iterate over all roots of unity.
    #[inline]
    fn build_unit_table() -> Vec<Self> {
        let n = Self::turn();
        (0..n).map(|k| Self::ccw().zz_pow(k as u8)).collect()
    }
}

// ----------------
// Subring-containment markers.

/// Ring is exactly `ZZ4` (all edges are axis-aligned unit segments).
/// Unlike `HasZZ4` (which includes `ZZ8`, `ZZ12`, etc.), only `ZZ4` itself
/// implements this.
pub trait IsZZ4Impl {}

/// Rings containing `ZZ4`.
pub trait HasZZ4Impl {}
pub trait HasZZ4: IsRingOrField + HasZZ4Impl + From<(IntT, IntT)> {}
impl<T: IsRingOrField + HasZZ4Impl + From<(IntT, IntT)>> HasZZ4 for T {}

pub trait IsZZ4: IsRingOrField + IsZZ4Impl {}
impl<T: IsRingOrField + IsZZ4Impl> IsZZ4 for T {}

/// Rings containing `ZZ6`.
pub trait HasZZ6Impl {}
pub trait HasZZ6: IsRingOrField + HasZZ6Impl {}
impl<T: IsRingOrField + HasZZ6Impl> HasZZ6 for T {}

/// Rings containing `ZZ8`.
pub trait HasZZ8Impl {}
pub trait HasZZ8: IsRingOrField + HasZZ8Impl {}
impl<T: IsRingOrField + HasZZ8Impl> HasZZ8 for T {}

/// Rings containing `ZZ10`.
pub trait HasZZ10Impl {}
pub trait HasZZ10: IsRingOrField + HasZZ10Impl {}
impl<T: IsRingOrField + HasZZ10Impl> HasZZ10 for T {}

/// Rings containing `ZZ12`.
pub trait HasZZ12Impl {}
pub trait HasZZ12: IsRingOrField + HasZZ12Impl {}
impl<T: IsRingOrField + HasZZ12Impl> HasZZ12 for T {}

impl<T: HasZZ4 + Units> OneImag for T {
    fn one_i() -> Self {
        <Self as Units>::unit(Self::qturn())
    }

    fn is_one_i(&self) -> bool {
        *self == Self::one_i()
    }
}
