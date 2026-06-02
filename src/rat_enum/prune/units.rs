//! Per-ring extraction of unit-vector coefficient arrays.

use crate::cyclotomic::{IsRing, Units};

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
    crate::dispatch_ring!(ring, extract::<ZZ, PHI, _>(|x| x.int_coeffs()))
}
