//! Generic implementation of division for cyclotomic rings (to get fields).

use num_traits::{One, Zero};

use super::symnum::SymNum;
use super::traits::IsRing;

/// Compute inverse by repeated rationalizing of the denominator.
///
/// IMPORTANT: This only works correctly if there are no nested roots in
/// the representation of the cyclotomic field
// FIXME: resolve this and unlock the other fields, if possible.
pub fn zz_inv<Z: SymNum + IsRing + One>(val: &Z) -> Z {
    // for x/y where y = a + b, b being a single (scaled) square root,
    // we compute y' = a - b and produce x*y'/(a+b)(a-b) = x*y'/a^2-b^2
    // where a^2 is a simpler term and b is rational.
    // we repeat this for all square roots in the denominator.

    // println!("------------------------\n{val}");

    let num_terms = Z::zz_params().sym_roots_num;

    let mut numer = Z::one();
    let mut denom = val.clone();

    // first ensure that the denominator is real
    let denom_conj = denom.conj();
    numer = numer * denom_conj;
    denom = denom * denom_conj;

    let mut root_ix = 0;
    let mut non_root_ix = 0;
    let mut root_found = true;
    while root_found {
        // println!("\n{numer}\n----\n{denom}\n");

        root_found = false;
        for i in 0..num_terms {
            if Z::zz_params().sym_roots_sqs[i].is_one() {
                non_root_ix = i; // we need the non-irrational part later
                continue; // non-irrational term
            }
            let c = denom.zz_coeffs()[i];
            if !c.is_zero() {
                root_found = true;
                root_ix = i;
                break;
            }
        }

        if !root_found {
            // println!("NO MORE ROOTS");
            break; // rational denominator
        }
        // println!(
        //     "ROOT FOUND: (#{root_ix}), coeff: {}",
        //     denom.zz_coeffs()[root_ix]
        // );

        // compute y' = a - b
        let denom_conj = {
            let curr_root_coeff = denom.zz_coeffs()[root_ix];
            let mut yy = denom.clone();
            yy.zz_coeffs_mut()[root_ix] = -curr_root_coeff;
            yy
        };
        // println!("CONJUGATED DENOM: {}", denom_conj);

        // update numerator (= x * y')
        numer = numer * denom_conj;
        // update denominator (= (a + b)(a - b) = a^2 - b^2)
        denom = denom * denom_conj;
    }

    // now we have a rational denominator (i.e. no square root terms)

    // just flip it and multiply with the numerator to get the final result
    let mut inv_denom = Z::zero();
    inv_denom.zz_coeffs_mut()[non_root_ix] = Z::Scalar::one() / denom.zz_coeffs()[non_root_ix];
    let inv_denom = inv_denom;

    return numer * inv_denom;
}
