use super::numtraits::{Ccw, Conj, InnerIntType, IntRing, OneImag, ZSigned};
use super::symnum::SymNum;
use super::symnum::{IntT, ZZComplex};
use super::units::Units;

pub trait RealTraits: ZSigned + Ord {}
pub trait ComplexTraits: Ccw + ZZComplex {}
pub trait RingTraits:
    InnerIntType + IntRing + Conj + From<IntT> + From<<Self as InnerIntType>::IntType>
{
}

pub trait IsSymNum: SymNum + RingTraits {}
impl<T: SymNum + RingTraits> IsSymNum for T {}

pub trait IsRingOrField: IsSymNum {
    type Real: IsRingOrField<Real = Self::Real, Complex = Self::Complex> + IsReal;
    type Complex: IsRingOrField<Real = Self::Real, Complex = Self::Complex> + IsComplex;
}

pub trait IsRing: IsSymNum {
    type Real: IsReal<Ring = Self::Real>
        + IsRing<Real = Self::Real, Complex = Self::Complex>
        + IsRingOrField<Real = Self::Real, Complex = Self::Complex>;
    type Complex: IsComplex<Ring = Self::Complex>
        + IsRing<Real = Self::Real, Complex = Self::Complex>
        + IsRingOrField<Real = Self::Real, Complex = Self::Complex>;
}

pub trait IsReal: IsSymNum + RealTraits {
    type Ring: IsRing<Real = Self::Ring>
        + IsRingOrField<Real = Self::Ring>
        + IsReal<Ring = Self::Ring>;
}
pub trait IsComplex: IsSymNum + ComplexTraits {
    type Ring: IsRing<Complex = Self::Ring>
        + IsRingOrField<Complex = Self::Ring>
        + IsComplex<Ring = Self::Ring>;
}

pub trait ZType: IsRing<Real = Self> + IsReal<Ring = Self> {
    type Complex: ZZType<Real = Self>;
}
pub trait ZZType: IsRing<Complex = Self> + IsComplex<Ring = Self> {
    type Real: ZType<Complex = Self>;
}

/// Ring is exactly ZZ4 (all edges are axis-aligned unit segments).
/// Unlike HasZZ4 (which includes ZZ8, ZZ12, etc.), only ZZ4 itself implements this.
pub trait IsZZ4Impl {}

/// rings containing ZZ4
pub trait HasZZ4Impl {}
pub trait HasZZ4: IsComplex + HasZZ4Impl + From<(IntT, IntT)> {}
impl<T: IsComplex + HasZZ4Impl + From<(IntT, IntT)>> HasZZ4 for T {}

pub trait IsZZ4: IsComplex + IsZZ4Impl {}
impl<T: IsComplex + IsZZ4Impl> IsZZ4 for T {}

/// rings containing ZZ6
pub trait HasZZ6Impl {}
pub trait HasZZ6: IsComplex + HasZZ6Impl {}
impl<T: IsComplex + HasZZ6Impl> HasZZ6 for T {}
/// rings containing ZZ8
pub trait HasZZ8Impl {}
pub trait HasZZ8: IsComplex + HasZZ8Impl {}
impl<T: IsComplex + HasZZ8Impl> HasZZ8 for T {}
/// rings containing ZZ10
pub trait HasZZ10Impl {}
pub trait HasZZ10: IsComplex + HasZZ10Impl {}
impl<T: IsComplex + HasZZ10Impl> HasZZ10 for T {}
/// rings containing ZZ12
pub trait HasZZ12Impl {}
pub trait HasZZ12: IsComplex + HasZZ12Impl {}
impl<T: IsComplex + HasZZ12Impl> HasZZ12 for T {}

impl<T: HasZZ4 + Units> OneImag for T {
    fn one_i() -> Self {
        <Self as Units>::unit(Self::qturn())
    }

    fn is_one_i(&self) -> bool {
        *self == Self::one_i()
    }
}
