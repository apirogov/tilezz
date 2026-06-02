//! Seed-based orchestration: the multi-threaded path + the
//! external `--list-seeds` / `--seed <prefix>` flow.
//!
//! See [`parallel`] for the in-process worker dispatch, and
//! [`external`] for the seed-listing + seed-resume helpers that
//! external orchestrators (xargs, GNU parallel, Slurm, ...) drive.
//!
//! [`super::EnumResult`] is the common `(rats, stats)` return type
//! used by every dispatch entry point.

pub mod external;
pub mod parallel;

pub use external::{
    collect_seed_prefixes, dispatch_collect_seed_prefixes, dispatch_enumerate_from_seed,
    enumerate_from_seed,
};
pub use parallel::{parallel_drain_seeds, rat_enum_parallel, splitting_depth};
