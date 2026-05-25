use super::numtraits::{
    Ccw, Conj, InnerIntType, IntRing, IntersectUnitSegments, OneImag, ReImSign,
};
use super::symnum::SymNum;
use super::symnum::{IntT, ZZComplex};
use super::units::Units;

// ----------------
// Building blocks composed by `IsRingOrField`.

pub trait ComplexTraits: Ccw + ZZComplex {}

pub trait RingTraits:
    InnerIntType + IntRing + Conj + From<IntT> + From<<Self as InnerIntType>::IntType>
{
}

pub trait IsSymNum: SymNum + RingTraits {}
impl<T: SymNum + RingTraits> IsSymNum for T {}

// ----------------
// Top of the hierarchy.
//
// `IsRingOrField` is the catch-all bound for a cyclotomic ring or field.
// Today every implementing type is a ring; a future `IsField` would also
// implement `IsRingOrField` and add a `Div<Self, Output = Self>` bound.
// `IsRing` keeps the ring-vs-field structural split alive so adding
// `IsField` later is a pure extension, not a rearrangement.

pub trait IsRingOrField:
    IsSymNum + ComplexTraits + ReImSign + IntersectUnitSegments
{
}

/// Integer ring (no division). Every current ZZ* type implements this.
pub trait IsRing: IsRingOrField {}

// Future: pub trait IsField: IsRingOrField + Div<Self, Output = Self> {}

/// Marker for the standalone Z[zeta_n] integer rings (today the ZZ* types).
///
/// `ZZType` is structural: a future `QQType` for cyclotomic fields would
/// sit alongside it on top of `IsField`.
pub trait ZZType: IsRing {}

// ----------------
// Subring-containment markers.
//
// `HasZZk<T>` says: T contains the k-th cyclotomic subring; `IsZZ4` says T
// *is* exactly ZZ4. These markers are unchanged in shape -- only the bound
// is loosened from the gone `IsComplex` to `IsRingOrField`.

/// Ring is exactly ZZ4 (all edges are axis-aligned unit segments).
/// Unlike HasZZ4 (which includes ZZ8, ZZ12, etc.), only ZZ4 itself implements this.
pub trait IsZZ4Impl {}

/// rings containing ZZ4
pub trait HasZZ4Impl {}
pub trait HasZZ4: IsRingOrField + HasZZ4Impl + From<(IntT, IntT)> {}
impl<T: IsRingOrField + HasZZ4Impl + From<(IntT, IntT)>> HasZZ4 for T {}

pub trait IsZZ4: IsRingOrField + IsZZ4Impl {}
impl<T: IsRingOrField + IsZZ4Impl> IsZZ4 for T {}

/// rings containing ZZ6
pub trait HasZZ6Impl {}
pub trait HasZZ6: IsRingOrField + HasZZ6Impl {}
impl<T: IsRingOrField + HasZZ6Impl> HasZZ6 for T {}
/// rings containing ZZ8
pub trait HasZZ8Impl {}
pub trait HasZZ8: IsRingOrField + HasZZ8Impl {}
impl<T: IsRingOrField + HasZZ8Impl> HasZZ8 for T {}
/// rings containing ZZ10
pub trait HasZZ10Impl {}
pub trait HasZZ10: IsRingOrField + HasZZ10Impl {}
impl<T: IsRingOrField + HasZZ10Impl> HasZZ10 for T {}
/// rings containing ZZ12
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
