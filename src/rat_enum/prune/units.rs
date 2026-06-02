//! Per-ring extraction of unit-vector coefficient arrays.

use crate::cyclotomic::{IsRing, Units, ZZ10, ZZ12, ZZ16, ZZ20, ZZ24, ZZ32, ZZ4, ZZ60, ZZ6, ZZ8};

/// Returns ALL `n` unit vectors `unit(0)..unit(n-1)` for the chosen
/// ring (NOT just the DFS's `n-1`-direction window). A partial rat's
/// cumulative facing can land on any absolute direction over enough
/// turns, so closure feasibility must be assessed against the full
/// set of possible unit vectors. Used by [`super::ModularPrune`] to
/// feed its BFS over modular displacements.
pub fn unit_vectors_for_ring(ring: u8) -> (Vec<Vec<i64>>, usize) {
    fn extract<ZZ: IsRing, const PHI: usize, F>(coeffs_fn: F) -> (Vec<Vec<i64>>, usize)
    where
        F: Fn(&ZZ) -> [i64; PHI],
    {
        let mut out = Vec::new();
        for d in 0..ZZ::turn() {
            let u: ZZ = <ZZ as Units>::unit(d);
            out.push(coeffs_fn(&u).to_vec());
        }
        (out, PHI)
    }
    match ring {
        4 => extract::<ZZ4, 2, _>(|x| x.int_coeffs()),
        6 => extract::<ZZ6, 2, _>(|x| x.int_coeffs()),
        8 => extract::<ZZ8, 4, _>(|x| x.int_coeffs()),
        10 => extract::<ZZ10, 4, _>(|x| x.int_coeffs()),
        12 => extract::<ZZ12, 4, _>(|x| x.int_coeffs()),
        16 => extract::<ZZ16, 8, _>(|x| x.int_coeffs()),
        20 => extract::<ZZ20, 8, _>(|x| x.int_coeffs()),
        24 => extract::<ZZ24, 8, _>(|x| x.int_coeffs()),
        32 => extract::<ZZ32, 16, _>(|x| x.int_coeffs()),
        60 => extract::<ZZ60, 16, _>(|x| x.int_coeffs()),
        _ => panic!("invalid ring selected"),
    }
}
