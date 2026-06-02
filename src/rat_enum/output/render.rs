//! Polyline conversion for the `--mode render` GIF path.
//!
//! `polygons<ZZ>` converts canonical angle sequences to floating-
//! point polylines (one per rat) for downstream rendering via
//! `tilezz::vis`. `run_rat_enum_polylines` is the runtime-ring
//! wrapper used by the CLI.

use crate::cyclotomic::{IsRing, ZZ10, ZZ12, ZZ16, ZZ20, ZZ24, ZZ32, ZZ4, ZZ60, ZZ8};
use crate::geom::rat::Rat;
use crate::geom::snake::Turtle;
use crate::rat_enum::enumerate_dispatch;
use crate::vis::plotutils::P64;

pub fn polygons<ZZ: IsRing>(rats: Vec<Vec<i8>>) -> Vec<Vec<P64>> {
    rats.into_iter()
        .map(|seq| Rat::<ZZ>::from_slice_unchecked(&seq).to_polyline_f64(Turtle::default()))
        .collect()
}

pub fn run_rat_enum_polylines(
    ring: u8,
    max_steps: usize,
    step: i8,
    n_threads: usize,
    dihedral: bool,
    paranoid: bool,
) -> Vec<Vec<P64>> {
    match ring {
        4 => polygons::<ZZ4>(
            enumerate_dispatch::<ZZ4>(max_steps, step, n_threads, dihedral, paranoid).0,
        ),
        8 => polygons::<ZZ8>(
            enumerate_dispatch::<ZZ8>(max_steps, step, n_threads, dihedral, paranoid).0,
        ),
        10 => polygons::<ZZ10>(
            enumerate_dispatch::<ZZ10>(max_steps, step, n_threads, dihedral, paranoid).0,
        ),
        12 => polygons::<ZZ12>(
            enumerate_dispatch::<ZZ12>(max_steps, step, n_threads, dihedral, paranoid).0,
        ),
        16 => polygons::<ZZ16>(
            enumerate_dispatch::<ZZ16>(max_steps, step, n_threads, dihedral, paranoid).0,
        ),
        20 => polygons::<ZZ20>(
            enumerate_dispatch::<ZZ20>(max_steps, step, n_threads, dihedral, paranoid).0,
        ),
        24 => polygons::<ZZ24>(
            enumerate_dispatch::<ZZ24>(max_steps, step, n_threads, dihedral, paranoid).0,
        ),
        32 => polygons::<ZZ32>(
            enumerate_dispatch::<ZZ32>(max_steps, step, n_threads, dihedral, paranoid).0,
        ),
        60 => polygons::<ZZ60>(
            enumerate_dispatch::<ZZ60>(max_steps, step, n_threads, dihedral, paranoid).0,
        ),
        _ => panic!("invalid ring selected"),
    }
}
