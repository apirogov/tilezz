//! Rat enumeration: depth-first walk over a cyclotomic ring's
//! unit-direction graph that finds every simple closed polygon
//! (`Rat`) up to a chosen perimeter.
//!
//! The enumeration is the engine behind the `rat_enum` binary; this
//! module exposes it as a library so the binary stays thin and so
//! tests (correctness cross-validation, OEIS pinning) can drive the
//! same code paths from `#[cfg(test)]` without spawning subprocesses.
//!
//! # Optional optimizations
//!
//! All on top of the baseline DFS, all gated by [`prune::Prunes`]:
//!
//! - **Modular reachability prune** ([`prune::modular`]): for a set
//!   of small moduli `m`, precomputes the set of mod-`m`
//!   displacements reachable by sums of `<= r` unit vectors; rejects
//!   any candidate whose post-step displacement isn't in the table
//!   for the remaining-steps slot. ~10-376× wall-time win on rings
//!   with non-trivial mod-2/3/4/6 structure.
//!
//! - **Closure-table prune** ([`prune::closure_table`]): pre-enumerates
//!   every simple open snake up to length `L` and stores their
//!   (endpoint, facing) "closure keys". When the DFS's remaining
//!   step budget is `<= L`, rejects any candidate whose required
//!   suffix key isn't in the table. Strictly stronger than any
//!   modular projection (uses exact lattice + facing info). Adds
//!   2-5× on top of the modular prune.
//!
//! # Output modes
//!
//! Driven by the binary's `--mode` flag (see `Mode` in
//! `src/bin/rat_enum.rs`); the library exposes the underlying
//! functions independent of any CLI framing.
//!
//! - In-memory: [`enumerate_dispatch`] / [`run_rat_enum_seqs`]
//!   return the full set as `Vec<Vec<i8>>`. Fine up to ~5 GB of
//!   peak RSS.
//! - Streaming (memory-bounded, for large `n`): [`stream`] runs the
//!   same DFS but routes closures through `FnMut(&[i8])` callbacks
//!   into per-thread sort-buffer run files, then a k-way merge
//!   stage produces a deduplicated sorted artifact plus a BLAKE3
//!   certificate. See [`stream::stream_enum_dispatch`] and
//!   [`stream::merge_runs`].
//!
//! See submodule docs for the specific contracts.
//!
//! # Correctness
//!
//! Every optimization has dedicated cross-validation tests against
//! the baseline DFS, plus an OEIS A316192 anchor through n=10. See
//! `opt_correctness_tests` in `src/bin/rat_enum.rs`.

pub mod canonical;
pub mod dfs;
pub mod output;
pub mod prune;
pub mod seed;
pub mod stats;
pub mod stream;

use crate::cyclotomic::IsRing;
use crate::rat_enum::canonical::make_ops;
use crate::rat_enum::dfs::rat_enum_with;
use crate::rat_enum::prune::snapshot_prunes;
use crate::rat_enum::seed::parallel::rat_enum_parallel;
use crate::rat_enum::stats::DfsStats;

/// Output of a single-ring enumeration: the canonical sequences
/// (sorted by length), plus the DFS stats counters.
pub type EnumResult = (Vec<Vec<i8>>, DfsStats);

/// Generic ring-specific dispatcher: picks single-threaded vs parallel
/// based on `n_threads`, builds the canonical-ops pair, snapshots the
/// global prune state, and runs the enumeration. Returns the canonical
/// sequences (sorted by length) plus DFS stats.
pub fn enumerate_dispatch<ZZ: IsRing>(
    max_steps: usize,
    step: i8,
    n_threads: usize,
    free: bool,
    paranoid: bool,
) -> EnumResult {
    let ops = make_ops(free);
    let label = if free {
        "free enumeration"
    } else {
        "enumeration"
    };
    let prefix = if free { "free " } else { "" };
    let prunes = snapshot_prunes();
    if n_threads <= 1 {
        rat_enum_with::<ZZ>(max_steps, step, ops, label, prefix, paranoid, &prunes)
    } else {
        rat_enum_parallel::<ZZ>(
            max_steps, step, n_threads, ops, label, prefix, paranoid, &prunes,
        )
    }
}

/// Runtime-ring dispatcher: like [`enumerate_dispatch`] but takes a
/// `ring: u8` and selects the underlying type by match. Used by the
/// CLI which only knows the ring number at runtime.
pub fn run_rat_enum_seqs(
    ring: u8,
    max_steps: usize,
    step: i8,
    n_threads: usize,
    free: bool,
    paranoid: bool,
) -> EnumResult {
    crate::dispatch_ring!(
        ring,
        enumerate_dispatch::<ZZ>(max_steps, step, n_threads, free, paranoid)
    )
}
