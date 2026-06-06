/// Exact source commit this build was compiled from, captured by
/// `build.rs` via `git rev-parse HEAD`. A `-dirty` suffix means the
/// working tree had uncommitted changes at build time (so the bare
/// commit does NOT fully describe the binary); `"unknown"` means it was
/// built outside a git checkout. Recorded in dataset RO-Crate metadata
/// and surfaced by each CLI's `--version`.
pub const GIT_COMMIT: &str = env!("TILEZZ_GIT_COMMIT");

/// Long version string for CLI `--version`: the crate semver plus the
/// exact build commit, e.g. `0.1.0 (dd4fd14...)` or, for a non-pristine
/// build, `0.1.0 (dd4fd14...-dirty)`. Lets a human or a reproduction
/// script confirm precisely which binary produced a result.
pub const VERSION: &str = concat!(
    env!("CARGO_PKG_VERSION"),
    " (",
    env!("TILEZZ_GIT_COMMIT"),
    ")"
);

pub mod cyclotomic;

pub mod geom;

pub mod analysis;

pub mod util;

pub mod stringmatch;

pub mod vis;

#[cfg(feature = "cli")]
pub mod rat_enum;

#[cfg(feature = "rat_explorer")]
pub mod web;
