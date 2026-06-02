//! Polyline conversion for the `--mode render` GIF path.
//!
//! `polygons<ZZ>` converts canonical angle sequences to floating-
//! point polylines (one per rat) for downstream rendering via
//! `tilezz::vis`. `run_rat_enum_polylines` is the runtime-ring
//! wrapper used by the CLI.

use crate::cyclotomic::IsRing;
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
    free: bool,
    paranoid: bool,
) -> Vec<Vec<P64>> {
    crate::dispatch_ring!(
        ring,
        polygons::<ZZ>(enumerate_dispatch::<ZZ>(max_steps, step, n_threads, free, paranoid).0)
    )
}
