use super::numtraits::{Ccw, Conj, InnerIntType, IntField, IntRing, OneImag, ZSigned};
use super::symnum::SymNum;
use super::symnum::{IntT, ZZComplex};

pub trait RealTraits: ZSigned + Ord {}
pub trait ComplexTraits: Ccw + ZZComplex {}
pub trait RingTraits:
    InnerIntType + IntRing + Conj + From<IntT> + From<<Self as InnerIntType>::IntType>
{
}

pub trait IsSymNum: SymNum + RingTraits {}
impl<T: SymNum + RingTraits> IsSymNum for T {}

pub trait FieldTraits: RingTraits + IntField + IsRealOrComplex {
    fn coerce_ring(self) -> <Self as IsRealOrComplex>::Ring
    where
        Self: Sized;
}

pub trait IsRingOrField: IsSymNum {
    type Real: IsRingOrField<Real = Self::Real, Complex = Self::Complex> + IsReal;
    type Complex: IsRingOrField<Real = Self::Real, Complex = Self::Complex> + IsComplex;
}
pub trait IsRealOrComplex: IsSymNum {
    type Ring: IsRealOrComplex<Ring = Self::Ring, Field = Self::Field> + IsRing;
    type Field: IsRealOrComplex<Ring = Self::Ring, Field = Self::Field> + IsField;
}

pub trait IsRing: IsSymNum {
    type Real: IsReal<Ring = Self::Real>
        + IsRealOrComplex<Ring = Self::Real>
        + IsRing<Real = Self::Real, Complex = Self::Complex>
        + IsRingOrField<Real = Self::Real, Complex = Self::Complex>;
    type Complex: IsComplex<Ring = Self::Complex>
        + IsRealOrComplex<Ring = Self::Complex>
        + IsRing<Real = Self::Real, Complex = Self::Complex>
        + IsRingOrField<Real = Self::Real, Complex = Self::Complex>;
}
pub trait IsField: IsSymNum + FieldTraits {
    type Real: IsReal<Field = Self::Real>
        + IsRealOrComplex<Field = Self::Real>
        + IsField<Real = Self::Real, Complex = Self::Complex>
        + IsRingOrField<Real = Self::Real, Complex = Self::Complex>;
    type Complex: IsComplex<Field = Self::Complex>
        + IsRealOrComplex<Field = Self::Complex>
        + IsField<Real = Self::Real, Complex = Self::Complex>
        + IsRingOrField<Real = Self::Real, Complex = Self::Complex>;
}

pub trait IsReal: IsSymNum + RealTraits {
    type Ring: IsRing<Real = Self::Ring>
        + IsRingOrField<Real = Self::Ring>
        + IsReal<Ring = Self::Ring, Field = Self::Field>
        + IsRealOrComplex<Ring = Self::Ring, Field = Self::Field>;
    type Field: IsField<Real = Self::Field>
        + IsRingOrField<Real = Self::Field>
        + IsReal<Ring = Self::Ring, Field = Self::Field>
        + IsRealOrComplex<Ring = Self::Ring, Field = Self::Field>;
}
pub trait IsComplex: IsSymNum + ComplexTraits {
    type Ring: IsRing<Complex = Self::Ring>
        + IsRingOrField<Complex = Self::Ring>
        + IsComplex<Ring = Self::Ring, Field = Self::Field>
        + IsRealOrComplex<Ring = Self::Ring, Field = Self::Field>;
    type Field: IsField<Complex = Self::Field>
        + IsRingOrField<Complex = Self::Field>
        + IsComplex<Ring = Self::Ring, Field = Self::Field>
        + IsRealOrComplex<Ring = Self::Ring, Field = Self::Field>;
}

pub trait ZType: IsRing<Real = Self> + IsReal<Ring = Self> {
    type Complex: ZZType<Real = Self>;
    type Field: QType<Ring = Self>;
}
pub trait ZZType: IsRing<Complex = Self> + IsComplex<Ring = Self> {
    type Real: ZType<Complex = Self>;
    type Field: QQType<Ring = Self>;
}
pub trait QType: IsField<Real = Self> + IsReal<Field = Self> {
    type Complex: QQType<Real = Self>;
    type Ring: ZType<Field = Self>;
}
pub trait QQType: IsField<Complex = Self> + IsComplex<Field = Self> {
    type Real: QType<Complex = Self>;
    type Ring: ZZType<Field = Self>;
}

/// rings containing ZZ4
pub trait HasZZ4Impl {}
pub trait HasZZ4: IsComplex + HasZZ4Impl + From<(IntT, IntT)> {}
impl<T: IsComplex + HasZZ4Impl + From<(IntT, IntT)>> HasZZ4 for T {}

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

impl<T: HasZZ4> OneImag for T {
    fn one_i() -> Self {
        Self::unit(Self::qturn())
    }

    fn is_one_i(&self) -> bool {
        *self == Self::one_i()
    }
}

impl<T: IsRealOrComplex + IsField> HasZZ4Impl for T where T::Ring: HasZZ4Impl {}
impl<T: IsRealOrComplex + IsField> HasZZ6Impl for T where T::Ring: HasZZ6Impl {}
impl<T: IsRealOrComplex + IsField> HasZZ8Impl for T where T::Ring: HasZZ8Impl {}
impl<T: IsRealOrComplex + IsField> HasZZ10Impl for T where T::Ring: HasZZ10Impl {}
impl<T: IsRealOrComplex + IsField> HasZZ12Impl for T where T::Ring: HasZZ12Impl {}
