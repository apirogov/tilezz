//! External-orchestration entry points.
//!
//! `--mode list-seeds` walks the DFS down to `split_depth` and
//! prints (1) `SEED <comma,sep,angles>` lines for alive prefixes
//! and (2) `RAT [...]` lines for polygons that closed before
//! reaching `split_depth`. `--seed <prefix>` then resumes the DFS
//! from a given prefix and prints `RAT` lines for what it finds.
//! The two outputs (one `--list-seeds` plus per-seed jobs) are
//! mergeable via `sort -u` on the `RAT` lines to recover the
//! one-shot result.

use std::collections::HashSet;

use crate::cyclotomic::IsRing;
use crate::geom::snake::Snake;
use crate::rat_enum::canonical::make_ops;
use crate::rat_enum::dfs::{collect_seeds, hashset_recorder, rat_enum_step, SeedGather};
use crate::rat_enum::prune::snapshot_prunes;
use crate::rat_enum::seed::parallel::{parallel_drain_seeds, splitting_depth};
use crate::rat_enum::stats::DfsStats;

type SeedListing = (HashSet<Vec<i8>>, Vec<Vec<i8>>);

/// Walk the DFS down to `split_depth` and return
/// `(closed, alive_prefixes)`. `closed` contains polygons that
/// closed *during the seed walk* (perimeter <= split_depth). The
/// caller is responsible for handling both: closed ones go straight
/// into the final result; alive prefixes are the work units to be
/// dispatched to per-seed jobs.
pub fn collect_seed_prefixes<ZZ: IsRing>(
    max_steps: usize,
    step: i8,
    split_depth: usize,
    dihedral: bool,
) -> (HashSet<Vec<i8>>, Vec<Vec<i8>>) {
    let ops = make_ops(dihedral);
    let prunes = snapshot_prunes();
    let mut snake: Snake<ZZ> = Snake::new();
    let mut seeds: Vec<Vec<i8>> = Vec::new();
    let mut closed: HashSet<Vec<i8>> = HashSet::new();
    let mut stats = DfsStats::default();
    {
        let mut record_closed = hashset_recorder(&mut closed);
        let mut gather = SeedGather {
            seeds: &mut seeds,
            record_closed: &mut record_closed,
            stats: &mut stats,
        };
        collect_seeds::<ZZ>(
            &mut snake,
            max_steps,
            step,
            split_depth,
            &mut gather,
            ops,
            false,
            &prunes,
        );
    }
    (closed, seeds)
}

/// Resume the DFS from `seed` (a walked angle prefix) and return
/// the unique canonical sequences found below it.
///
/// * `n_threads == 1` (recommended): single-threaded fast path. Runs
///   `rat_enum_step` directly from the seed's snake state. Holds at
///   most one `HashSet` (this seed's share of the final set).
/// * `n_threads > 1`: sub-splits the seed's subtree at depth
///   `seed.len() + splitting_depth(n_threads, branching)` (clamped
///   to `max_steps`) to produce ~10*n_threads work units, then
///   dispatches them across worker threads via `parallel_drain_seeds`.
///
/// # Memory profile and why separate processes are strictly better
///
/// **Single-threaded:** ~1x the final set size at peak (one HashSet,
/// growing with resize transients to ~3x during a bucket-array
/// doubling).
///
/// **Multi-threaded:** ~5x the final set size at peak. Each of
/// `n_threads` workers holds a thread-local HashSet that lives for
/// the whole DFS phase, plus the merge double-buffers as the main
/// `closed_main` absorbs each local. This is the same memory profile
/// as the whole-tree parallel mode.
///
/// **If memory matters and you have multiple seeds to process,
/// running them as separate single-threaded processes (e.g.
/// `xargs -P 16 ... --seed {} --threads 1`) is strictly better:**
/// each process holds only its own share of the final set, and the
/// OS reclaims everything (allocator arenas included) when each
/// process exits. Use `n_threads > 1` only when an orchestrator
/// constrains you to one process at a time.
pub fn enumerate_from_seed<ZZ: IsRing>(
    max_steps: usize,
    step: i8,
    seed: &[i8],
    n_threads: usize,
    dihedral: bool,
    paranoid: bool,
) -> HashSet<Vec<i8>> {
    let ops = make_ops(dihedral);
    let prunes = snapshot_prunes();

    if n_threads <= 1 {
        // Single-threaded fast path. No sub-split, no worker spawn:
        // just walk the DFS from this seed's snake state directly.
        let mut snake: Snake<ZZ> = Snake::from_slice_unsafe(seed);
        let mut local: HashSet<Vec<i8>> = HashSet::new();
        let mut stats = DfsStats::default();
        {
            let mut record = hashset_recorder(&mut local);
            rat_enum_step::<ZZ>(
                &mut snake,
                max_steps,
                step,
                &mut record,
                &mut stats,
                ops,
                paranoid,
                &prunes,
            );
        }
        return local;
    }

    // Multi-threaded: sub-split the seed's subtree to ~10*n_threads
    // work units (capped at `max_steps`), then dispatch.
    let hm1 = (ZZ::hturn() as usize).saturating_sub(1);
    let branching = 2 * (hm1 / step.max(1) as usize) + 1;
    let sub_split = splitting_depth(n_threads, branching);
    let split_depth = (seed.len() + sub_split).min(max_steps);

    let mut closed_main: HashSet<Vec<i8>> = HashSet::new();
    let mut sub_seeds: Vec<Vec<i8>> = Vec::new();
    let mut seed_stats = DfsStats::default();
    {
        let mut snake: Snake<ZZ> = Snake::from_slice_unsafe(seed);
        let mut record_closed = hashset_recorder(&mut closed_main);
        let mut gather = SeedGather {
            seeds: &mut sub_seeds,
            record_closed: &mut record_closed,
            stats: &mut seed_stats,
        };
        collect_seeds::<ZZ>(
            &mut snake,
            max_steps,
            step,
            split_depth,
            &mut gather,
            ops,
            paranoid,
            &prunes,
        );
    }

    let (merged, _) = parallel_drain_seeds::<ZZ>(
        &sub_seeds,
        closed_main,
        seed_stats,
        max_steps,
        step,
        n_threads,
        ops,
        paranoid,
        &prunes,
    );
    merged
}

/// Runtime-ring dispatch wrapper for [`collect_seed_prefixes`].
pub fn dispatch_collect_seed_prefixes(
    ring: u8,
    max_steps: usize,
    step: i8,
    split_depth: usize,
    dihedral: bool,
) -> SeedListing {
    crate::dispatch_ring!(
        ring,
        collect_seed_prefixes::<ZZ>(max_steps, step, split_depth, dihedral)
    )
}

/// Runtime-ring dispatch wrapper for [`enumerate_from_seed`].
pub fn dispatch_enumerate_from_seed(
    ring: u8,
    max_steps: usize,
    step: i8,
    seed: &[i8],
    n_threads: usize,
    dihedral: bool,
    paranoid: bool,
) -> HashSet<Vec<i8>> {
    crate::dispatch_ring!(
        ring,
        enumerate_from_seed::<ZZ>(max_steps, step, seed, n_threads, dihedral, paranoid)
    )
}
