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

use num_complex::Complex64;
use num_integer::Integer;
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

/// Complex conjugation. Each cyclotomic ring implements its own;
/// real-only scalars (i32/i64) would be the identity but the trait is
/// only used by complex ring elements.
pub trait Conj {
    fn conj(&self) -> Self;
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

/// Integer cell coordinates `(floor(Re(self)), floor(Im(self)))` for
/// spatial bucketing in `intgeom::grid::UnitSquareGrid`.
///
/// Two entry points:
///
/// * [`cell_floor`](Self::cell_floor) -- the fast path. Returns the
///   f64-only floor of `complex64()`. In debug builds it cross-checks
///   against `cell_floor_exact` via `debug_assert_eq!`. This is what
///   the grid uses in the hot path (spatial bucketing only needs
///   deterministic answers, not bit-exact ones at boundary points).
///
/// * [`cell_floor_exact`](Self::cell_floor_exact) -- the slow path,
///   ~20% of `patch_enum`/`rat_enum` runtime when called per insertion.
///   Per-ring impls use the f64 floor as a hint, then verify each axis
///   via integer ring arithmetic + `re_sign`/`im_sign` checks. f64 is
///   only a hint; output is bit-exact.
///
/// Per-ring exact implementation strategy:
///
/// * For rings containing `i` (`HasZZ4`-style: ZZ4, ZZ8, ZZ12, ZZ16,
///   ZZ20, ZZ24, ZZ32, ZZ60): cell-corner construction goes through
///   `From<(i64, i64)>` and verification through the `rect_signs`
///   primitive.
/// * For rings without `i` (today: ZZ10): per-axis sign verification
///   using ring-specific shifted-sign primitives (no `i` needed).
pub trait CellFloor: SymNum {
    /// Exact `(floor(Re), floor(Im))`. Bit-exact at every point; can
    /// disagree with the fast `cell_floor` only on the measure-zero
    /// set of points exactly on a half-integer grid line.
    fn cell_floor_exact(&self) -> (i64, i64);

    /// Fast `(floor(Re), floor(Im))` via f64. In debug builds,
    /// `debug_assert_eq!`s against `cell_floor_exact` so any
    /// boundary-case divergence shows up in tests.
    ///
    /// Canonical ring elements with small integer coefficients
    /// essentially never land exactly on a half-integer line, so the
    /// fast path matches the exact path in practice.
    #[inline]
    fn cell_floor(&self) -> (i64, i64) {
        let c = self.complex64();
        let fast = (c.re.floor() as i64, c.im.floor() as i64);
        debug_assert_eq!(
            fast,
            self.cell_floor_exact(),
            "fast cell_floor disagrees with cell_floor_exact",
        );
        fast
    }
}

/// Predicate "this point is within Euclidean radius `r` of the origin",
/// i.e. `|self|^2 <= r * r` with `r >= 0`.
///
/// Used by the DFS reachability pruning in `rat_enum`. Each ring routes
/// through `impl_integral_within_radius_via_norm_sq!` (exact via
/// `(r^2 - z*conj(z)).re_sign() >= 0`, pure ring arithmetic) or a
/// hand-rolled pure-i64 fast path (ZZ12).
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
        Self: IsRing,
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
// Catch-all ring bound and structural marker traits.

/// Cached lookup of roots of unity for a concrete ring.
///
/// Implementations are emitted by `impl_integral_units_via_basis!` (a
/// `OnceLock`-cached unit table built from `derive_units_lookup`) or by a
/// hand-rolled per-ring `static UNIT_TABLE: [[i64; PHI]; N] = [...]` for
/// hot-path rings (ZZ12).
pub trait Units: Copy + One + Ccw {
    /// Unit-length vector pointing in direction of `angle`.
    fn unit(angle: i8) -> Self;
}

/// **The cyclotomic-ring bound.** Every concrete `ZZn` ring implements
/// this catch-all collection, so downstream callers only need to write
/// `T: IsRing` (plus any subring-containment marker like `HasZZ4`) rather
/// than spelling out every individual capability.
///
/// In math, a field is a ring with division. The hierarchy here mirrors
/// that: a future `IsField` would be `IsRing + Div<Self, Output = Self>`,
/// so a function written against `T: IsRing` automatically accepts a
/// field too. Functions that actually require division ask for
/// `T: IsField` (when that exists).
pub trait IsRing:
    SymNum
    + Ccw
    + Conj
    + ZZComplex
    + ReImSign
    + IntersectUnitSegments
    + WithinRadius
    + CellFloor
    + Units
    + InnerIntType
    + IntRing
    + From<IntT>
    + From<<Self as InnerIntType>::IntType>
{
}

// Future: pub trait IsField: IsRing + Div<Self, Output = Self> {}

// ----------------
// Subring-containment markers.

/// Ring is exactly `ZZ4` (all edges are axis-aligned unit segments).
/// Unlike `HasZZ4` (which includes `ZZ8`, `ZZ12`, etc.), only `ZZ4` itself
/// implements this.
pub trait IsZZ4Impl {}

/// Rings containing `ZZ4`. Implies `IsRing` plus the lattice-point
/// constructor `From<(i64, i64)>`, since `i = unit(qturn)` is in the ring.
pub trait HasZZ4Impl {}
pub trait HasZZ4: IsRing + HasZZ4Impl + From<(IntT, IntT)> {}
impl<T: IsRing + HasZZ4Impl + From<(IntT, IntT)>> HasZZ4 for T {}

pub trait IsZZ4: IsRing + IsZZ4Impl {}
impl<T: IsRing + IsZZ4Impl> IsZZ4 for T {}

/// Rings containing `ZZ6`.
pub trait HasZZ6Impl {}
pub trait HasZZ6: IsRing + HasZZ6Impl {}
impl<T: IsRing + HasZZ6Impl> HasZZ6 for T {}

/// Rings containing `ZZ8`.
pub trait HasZZ8Impl {}
pub trait HasZZ8: IsRing + HasZZ8Impl {}
impl<T: IsRing + HasZZ8Impl> HasZZ8 for T {}

/// Rings containing `ZZ10`.
pub trait HasZZ10Impl {}
pub trait HasZZ10: IsRing + HasZZ10Impl {}
impl<T: IsRing + HasZZ10Impl> HasZZ10 for T {}

/// Rings containing `ZZ12`.
pub trait HasZZ12Impl {}
pub trait HasZZ12: IsRing + HasZZ12Impl {}
impl<T: IsRing + HasZZ12Impl> HasZZ12 for T {}

impl<T: HasZZ4> OneImag for T {
    fn one_i() -> Self {
        <Self as Units>::unit(Self::qturn())
    }

    fn is_one_i(&self) -> bool {
        *self == Self::one_i()
    }
}
