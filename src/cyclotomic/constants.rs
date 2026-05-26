//! Useful constants available in the cyclotomic rings and fields.

use num_traits::One;

use super::traits::{HasZZ10, HasZZ12, HasZZ8, IsRingOrField};
use super::rings::{ZZ10, ZZ20};
use super::traits::Units;

// Returns the sum of all units of a complex integer ring.
pub fn zz_units_sum<T: IsRingOrField + Units>() -> T {
    let mut p = T::zero();
    for i in 0..T::turn() {
        p = p + <T as Units>::unit(i).scale(i as i64);
    }
    p
}

// NOTE: as we can get the real-valued square roots represented,
// it means that we can represent any linear combination
// in a ring that supports quarter turn rotation (i.e. ZZDiv4).

pub fn sqrt2<T: IsRingOrField + HasZZ8 + Units>() -> T {
    let sc = T::turn() / 8;
    <T as Units>::unit(sc) + <T as Units>::unit(-sc)
}

pub fn sqrt3<T: IsRingOrField + HasZZ12 + Units>() -> T {
    let sc = T::turn() / 12;
    <T as Units>::unit(sc) + <T as Units>::unit(-sc)
}

pub fn sqrt5<T: IsRingOrField + HasZZ10 + Units>() -> T {
    let sc = T::turn() / 10;
    (<T as Units>::unit(sc) + <T as Units>::unit(-sc)) * T::one().scale(2) - T::one()
}

pub fn sqrt6<T: IsRingOrField + HasZZ8 + HasZZ12 + Units>() -> T {
    let sc = T::turn() / 24;
    (<T as Units>::unit(sc) + <T as Units>::unit(-sc)) * T::one().scale(2) - sqrt2::<T>()
}

// misc. irregular
pub fn zz10_isqrt_penta() -> ZZ10 {
    <ZZ10 as Units>::unit(1) * ZZ10::from(4) - ZZ10::one() - sqrt5()
}
pub fn zz20_half_sqrt_penta() -> ZZ20 {
    <ZZ20 as Units>::unit(3) + <ZZ20 as Units>::unit(-3)
}

#[cfg(test)]
mod tests {
    use super::super::rings::{ZZ10, ZZ12, ZZ16, ZZ20, ZZ24, ZZ32, ZZ60, ZZ8, ZZ10_Y};
    use super::super::traits::SymNum;
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
        assert_eq!(sqrt5::<ZZ60>().complex64().re, sq5);

        assert_eq!(sqrt6::<ZZ24>().complex64().re, sq6);
        // assert_eq!(sqrt6::<ZZ120>().complex().re, sq6);
        // assert_eq!(sqrt6::<ZZ240>().complex().re, sq6);

        assert_eq!(zz10_isqrt_penta().complex64().im, sq_penta);
        assert_eq!(zz20_half_sqrt_penta().complex64().re, hsq_penta);
    }
}
