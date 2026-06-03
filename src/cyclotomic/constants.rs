//! Distinguished ring elements expressed in `ZZn` arithmetic.
//!
//! The `sqrt2` / `sqrt3` / `sqrt5` / `sqrt6` family builds the named real
//! square root as a ring element via `unit(k) + unit(-k) = 2*cos(2*pi*k/n)`
//! at a suitably chosen angle `k`. The result is a *real* ring element
//! whose `complex64().re` is bit-exact with the corresponding f64
//! `sqrt(N)`. Each function is gated on a `HasZZk` marker that guarantees
//! the chosen angle exists in the ring.

use num_traits::One;

use super::rings::{ZZ10, ZZ20};
use super::traits::{HasZZ8, HasZZ10, HasZZ12, IsRing, Units};

/// Sum of all unit roots scaled by their angle index: `sum_k k * unit(k)`.
/// Lands at a generic, "interesting" ring point useful for property tests.
pub fn zz_units_sum<T: IsRing>() -> T {
    let mut p = T::zero();
    for i in 0..T::turn() {
        p = p + <T as Units>::unit(i).scale(i as i64);
    }
    p
}

/// `sqrt(2)` as a ring element of any ring containing ZZ8.
pub fn sqrt2<T: HasZZ8>() -> T {
    let sc = T::turn() / 8;
    <T as Units>::unit(sc) + <T as Units>::unit(-sc)
}

/// `sqrt(3)` as a ring element of any ring containing ZZ12.
pub fn sqrt3<T: HasZZ12>() -> T {
    let sc = T::turn() / 12;
    <T as Units>::unit(sc) + <T as Units>::unit(-sc)
}

/// `sqrt(5)` as a ring element of any ring containing ZZ10. Identity:
/// `2*(unit(turn/10) + unit(-turn/10)) - 1 = 2*(2*cos(36)) - 1 = sqrt(5)`.
pub fn sqrt5<T: HasZZ10>() -> T {
    let sc = T::turn() / 10;
    (<T as Units>::unit(sc) + <T as Units>::unit(-sc)) * T::one().scale(2) - T::one()
}

/// `sqrt(6)` as a ring element of any ring containing both ZZ8 and ZZ12
/// (i.e. ZZ24, the smallest ring with both). Identity:
/// `2*(unit(turn/24) + unit(-turn/24)) - sqrt(2) = 2*(2*cos(15)) - sqrt(2) = sqrt(6)`.
pub fn sqrt6<T: HasZZ8 + HasZZ12>() -> T {
    let sc = T::turn() / 24;
    (<T as Units>::unit(sc) + <T as Units>::unit(-sc)) * T::one().scale(2) - sqrt2::<T>()
}

// ----------------
// Misc ring-specific square roots that don't fit the HasZZk pattern.

/// The pentagonal radical `i * sqrt(10 - 2*sqrt(5))` as a ZZ10 element.
pub fn zz10_isqrt_penta() -> ZZ10 {
    <ZZ10 as Units>::unit(1) * ZZ10::from(4) - ZZ10::one() - sqrt5()
}

/// `sqrt(10 - 2*sqrt(5)) / 2` as a ZZ20 element.
pub fn zz20_half_sqrt_penta() -> ZZ20 {
    <ZZ20 as Units>::unit(3) + <ZZ20 as Units>::unit(-3)
}

#[cfg(test)]
mod tests {
    use super::super::rings::{ZZ8, ZZ10, ZZ10_Y, ZZ12, ZZ16, ZZ20, ZZ24, ZZ32, ZZ60};
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
