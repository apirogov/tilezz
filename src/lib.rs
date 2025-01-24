#[macro_use]
extern crate arrayref;

pub mod traits;

pub mod gaussint;

#[macro_use]
pub mod zzbase;

pub mod zzsigned;

pub mod zzparams;

#[macro_use]
pub mod zz;

pub mod qq;

pub mod zzgeom;

pub mod angles;

pub mod grid;

pub mod snake;

pub mod rat;

pub mod cyclotomic;

pub mod prelude;

pub mod plotutils;

#[cfg(feature = "plotters")]
pub mod plotters;
