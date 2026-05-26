//! One-shot codegen for ZZ60's `re_decomp` and `im_decomp` tables.
//!
//! Queries the legacy `ZZ60` (still present in `types.rs` at the time
//! this script runs) for `unit(k)` over `k = 0..16` and extracts the
//! `Re`/`Im` coefficient K-vectors in the symbolic basis
//! `{sqrt(roots[j])}` with the implicit `/scaling_fac` denominator. Emits
//! the resulting tables as Rust source ready to paste into `rings.rs`.
//!
//! Run via:
//!   cargo run --bin codegen_zz60
//!
//! Then paste the printed `re_decomp:` and `im_decomp:` blocks into the
//! corresponding `define_integral_zz!(name: ZZ60, ...)` invocation.
//!
//! This binary is intentionally a one-shot tool -- delete after use.

use tilezz::cyclotomic::{SymNum, Units, ZZ60};

fn main() {
    let phi = 16usize; // phi(60)
    let k = 8usize; // sym_roots_num
    let scaling_fac = ZZ60::zz_params().scaling_fac;
    eprintln!("ZZ60: phi = {phi}, K = {k}, scaling_fac = {scaling_fac}");

    println!("    // Auto-generated from legacy ZZ60 (codegen_zz60.rs).");
    println!("    // Re(zeta^k) in basis {{1, sqrt(3), sqrt(5), sqrt(10-2sqrt(5)),");
    println!("    //                        sqrt(15), sqrt(3(10-2sqrt(5))),");
    println!("    //                        sqrt(5(10-2sqrt(5))), sqrt(15(10-2sqrt(5)))}}");
    println!("    // with implicit /{scaling_fac}.");
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
        let z = <ZZ60 as Units>::unit(kk as i8);
        let coeffs = z.zz_coeffs();
        let row: Vec<i64> = (0..k).map(|j| scale(&coeffs[j].real)).collect();
        let entries = row
            .iter()
            .map(|v| format!("{v}"))
            .collect::<Vec<_>>()
            .join(", ");
        println!("        [{entries}],");
    }
    println!("    ],");

    println!("    im_decomp: [");
    for kk in 0..phi {
        let z = <ZZ60 as Units>::unit(kk as i8);
        let coeffs = z.zz_coeffs();
        let row: Vec<i64> = (0..k).map(|j| scale(&coeffs[j].imag)).collect();
        let entries = row
            .iter()
            .map(|v| format!("{v}"))
            .collect::<Vec<_>>()
            .join(", ");
        println!("        [{entries}],");
    }
    println!("    ],");

    eprintln!("// done.");
}
