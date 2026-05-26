//! One-shot codegen for ZZ32's `re_decomp` and `im_decomp` tables.
//! See `codegen_zz60.rs` for the same recipe.

use tilezz::cyclotomic::{SymNum, Units, ZZ32};

fn main() {
    let phi = 16usize; // phi(32)
    let k = 8usize; // sym_roots_num
    let scaling_fac = ZZ32::zz_params().scaling_fac;
    eprintln!("ZZ32: phi = {phi}, K = {k}, scaling_fac = {scaling_fac}");

    let scale = |r: &num_rational::Ratio<i64>| -> i64 {
        let scaled = r * num_rational::Ratio::<i64>::from_integer(scaling_fac);
        assert_eq!(
            *scaled.denom(),
            1,
            "post-scaling denom != 1: {}/{}",
            scaled.numer(),
            scaled.denom()
        );
        *scaled.numer()
    };

    println!("    re_decomp: [");
    for kk in 0..phi {
        let z = <ZZ32 as Units>::unit(kk as i8);
        let coeffs = z.zz_coeffs();
        let row: Vec<i64> = (0..k).map(|j| scale(&coeffs[j].real)).collect();
        println!(
            "        [{}],",
            row.iter().map(i64::to_string).collect::<Vec<_>>().join(", ")
        );
    }
    println!("    ],");

    println!("    im_decomp: [");
    for kk in 0..phi {
        let z = <ZZ32 as Units>::unit(kk as i8);
        let coeffs = z.zz_coeffs();
        let row: Vec<i64> = (0..k).map(|j| scale(&coeffs[j].imag)).collect();
        println!(
            "        [{}],",
            row.iter().map(i64::to_string).collect::<Vec<_>>().join(", ")
        );
    }
    println!("    ],");

    eprintln!("// done.");
}
