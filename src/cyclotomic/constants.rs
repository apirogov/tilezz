//! Useful constants available in the cyclotomic rings and fields.

use num_traits::One;

use super::symnum::SymNum;
use super::traits::{HasZZ10, HasZZ12, HasZZ8, IsComplex, IsRingOrField};
use super::types::{ZZ10, ZZ20, ZZ6};

// Returns the sum of all units of a complex integer ring.
pub fn zz_units_sum<T: IsRingOrField + IsComplex>() -> T {
    let mut p = T::zero();
    for i in 0..T::turn() {
        p = p + T::unit(i).scale(i as i64);
    }
    p
}

// NOTE: as we can get the real-valued square roots represented,
// it means that we can represent any linear combination
// in a ring that supports quarter turn rotation (i.e. ZZDiv4).

pub fn sqrt2<T: IsComplex + HasZZ8>() -> T {
    let sc = T::zz_params().full_turn_steps / 8;
    T::unit(sc) + T::unit(-sc)
}

pub fn sqrt3<T: IsComplex + HasZZ12>() -> T {
    let sc = T::zz_params().full_turn_steps / 12;
    T::unit(sc) + T::unit(-sc)
}

pub fn sqrt5<T: IsComplex + HasZZ10>() -> T {
    let sc = T::zz_params().full_turn_steps / 10;
    (T::unit(sc) + T::unit(-sc)) * T::one().scale(2) - T::one()
}

pub fn sqrt6<T: IsComplex + HasZZ8 + HasZZ12>() -> T {
    let sc = T::zz_params().full_turn_steps / 24;
    (T::unit(sc) + T::unit(-sc)) * T::one().scale(2) - sqrt2::<T>()
}

// misc. irregular
pub fn isqrt3() -> ZZ6 {
    ZZ6::unit(1) + ZZ6::unit(2)
}
pub fn zz10_isqrt_penta() -> ZZ10 {
    ZZ10::unit(1) * ZZ10::from(4) - ZZ10::one() - sqrt5()
}
pub fn zz20_half_sqrt_penta() -> ZZ20 {
    ZZ20::unit(3) + ZZ20::unit(-3)
}

#[cfg(test)]
mod tests {
    use super::super::params::ZZ10_Y;
    use super::super::types::{ZZ10, ZZ12, ZZ16, ZZ20, ZZ24, ZZ30, ZZ32, ZZ60, ZZ8};
    use super::*;

    #[test]
    fn test_constants() {
        use std::f64::consts::SQRT_2;

        let sq2 = SQRT_2;
        let sq3 = 3.0_f64.sqrt();
        let sq_penta = ZZ10_Y.sqrt();
        let hsq_penta = 0.5 * ZZ10_Y.sqrt();
        let sq5 = 5.0_f64.sqrt();
        let sq6 = 6.0_f64.sqrt();

        assert_eq!(sqrt2::<ZZ8>().complex64().re, sq2);
        assert_eq!(sqrt2::<ZZ16>().complex64().re, sq2);
        assert_eq!(sqrt2::<ZZ24>().complex64().re, sq2);
        assert_eq!(sqrt2::<ZZ32>().complex64().re, sq2);

        assert_eq!(sqrt3::<ZZ12>().complex64().re, sq3);
        assert_eq!(sqrt3::<ZZ24>().complex64().re, sq3);
        assert_eq!(sqrt3::<ZZ60>().complex64().re, sq3);

        assert_eq!(sqrt5::<ZZ10>().complex64().re, sq5);
        assert_eq!(sqrt5::<ZZ20>().complex64().re, sq5);
        assert_eq!(sqrt5::<ZZ30>().complex64().re, sq5);
        assert_eq!(sqrt5::<ZZ60>().complex64().re, sq5);

        assert_eq!(sqrt6::<ZZ24>().complex64().re, sq6);
        // assert_eq!(sqrt6::<ZZ120>().complex().re, sq6);
        // assert_eq!(sqrt6::<ZZ240>().complex().re, sq6);

        assert_eq!(isqrt3().complex64().im, sq3);
        assert_eq!(zz10_isqrt_penta().complex64().im, sq_penta);
        assert_eq!(zz20_half_sqrt_penta().complex64().re, hsq_penta);
    }
}
