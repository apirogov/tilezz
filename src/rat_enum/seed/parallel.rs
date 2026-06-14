//! In-process multi-threaded DFS.
//!
//! `rat_enum_parallel` is the entry: walk down to `split_depth` to
//! produce alive prefixes (the "seeds"), then dispatch them across
//! `n_threads` workers via an atomic counter. Each worker keeps its
//! own `HashSet` of canonical sequences; the main thread merges them
//! at the end. Memory profile: `~5x` the final set size at peak (per-
//! thread sets + merge double-buffer). See the file-level docs on
//! the binary for why running independent single-threaded processes
//! is strictly better when memory matters.

use std::collections::HashSet;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;

use crate::cyclotomic::IsRing;
use crate::geom::snake::Snake;
use crate::rat_enum::canonical::CanonicalOps;
use crate::rat_enum::dfs::{SeedGather, collect_seeds, hashset_recorder, rat_enum_step};
use crate::rat_enum::prune::Prunes;
use crate::rat_enum::stats::DfsStats;

/// Per-level DFS branching factor: the number of step-valid turn
/// directions the loop `(-hturn+1)..hturn` walks, i.e.
/// `2 * ((hturn - 1) / step) + 1`. E.g. ZZ4 step1 -> 3, ZZ12 step1 -> 11,
/// ZZ14 step2 (the ZZ7 subring) -> 7. Single source for the seed-split
/// sizing (`splitting_depth`) and the cylinder odometer
/// (`super::super::stream::progress::odometer_fraction`).
pub fn branch_factor(hturn: i8, step: i8) -> usize {
    let hm1 = (hturn.max(1) - 1) as usize;
    2 * (hm1 / step.max(1) as usize) + 1
}

/// Pick a DFS splitting depth such that the seed-walk produces
/// roughly `10 * n_threads` work units. With a branching factor of
/// `b` candidate directions per level, depth `d` enumerates at most
/// `b^d` seeds (the canonical-rotation + intersect + reachability
/// prunes knock that down further, but the raw count is the right
/// upper bound for sizing). We invert: `d = ceil(log_b(10 * threads))`.
///
/// `branching` is the per-level branching factor of the DFS:
/// `2 * hturn - 1` for a ZZ ring (the loop walks `(-hturn+1)..hturn`).
/// E.g. ZZ4 -> 3, ZZ12 -> 11, ZZ24 -> 23.
pub fn splitting_depth(n_threads: usize, branching: usize) -> usize {
    if n_threads <= 1 || branching <= 1 {
        return 0;
    }
    let target = (10 * n_threads) as f64;
    let depth = (target.ln() / (branching as f64).ln()).ceil() as usize;
    depth.max(1)
}

/// Parallel variant of single-threaded `rat_enum_with`: splits the
/// DFS at `split_depth` (selected via [`splitting_depth`]), then
/// hands the resulting alive prefixes out to `n_threads` worker
/// threads via a shared atomic counter. Each worker keeps its own
/// `HashSet` and the main thread merges the per-worker sets at the
/// end.
#[allow(clippy::too_many_arguments)]
pub fn rat_enum_parallel<ZZ: IsRing>(
    max_steps: usize,
    step: i8,
    n_threads: usize,
    ops: CanonicalOps,
    label: &str,
    prefix: &str,
    paranoid: bool,
    prunes: &Prunes,
) -> (Vec<Vec<i8>>, DfsStats) {
    let branching = branch_factor(ZZ::hturn(), step);
    let split_depth = splitting_depth(n_threads, branching);

    println!("-------- {label} started --------");
    if paranoid {
        println!("paranoid: per-step fresh-snake cross-check enabled");
    }
    println!("parallel: n_threads={n_threads} branching={branching} split_depth={split_depth}");

    let mut closed_main: HashSet<Vec<i8>> = HashSet::new();
    let mut seeds: Vec<Vec<i8>> = Vec::new();
    let mut seed_stats = DfsStats::default();
    {
        let mut snake: Snake<ZZ> = Snake::new();
        let mut record_closed = hashset_recorder(&mut closed_main);
        let mut gather = SeedGather {
            seeds: &mut seeds,
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
            prunes,
        );
    }
    println!("parallel: {} seed states collected", seeds.len());

    let (merged, worker_stats) = parallel_drain_seeds::<ZZ>(
        &seeds,
        closed_main,
        seed_stats,
        max_steps,
        step,
        n_threads,
        ops,
        paranoid,
        prunes,
    );

    println!(
        "-------- {label} completed --------\n{prefix}{} rats found",
        merged.len()
    );

    let mut result: Vec<Vec<i8>> = merged.into_iter().collect();
    result.sort_by_key(|x| x.len());
    (result, worker_stats)
}

/// Dispatch a collected list of seed prefixes across `n_threads`
/// worker threads via an atomic counter. Each worker takes seeds one
/// at a time, runs `rat_enum_step` from the seed's snake state, and
/// accumulates canonical sequences into a thread-local HashSet. At
/// the end the locals are folded into `closed_main` (which already
/// contains any polygons that closed during seed collection).
///
/// Used by both [`rat_enum_parallel`] (whole-tree enumeration, seeds
/// collected from the root) and `enumerate_from_seed` with threads>1
/// (single-seed sub-tree enumeration, sub-seeds collected from a
/// given prefix).
#[allow(clippy::too_many_arguments)]
pub fn parallel_drain_seeds<ZZ: IsRing>(
    seeds: &[Vec<i8>],
    closed_main: HashSet<Vec<i8>>,
    seed_stats: DfsStats,
    max_steps: usize,
    step: i8,
    n_threads: usize,
    ops: CanonicalOps,
    paranoid: bool,
    prunes: &Prunes,
) -> (HashSet<Vec<i8>>, DfsStats) {
    let next_idx = AtomicUsize::new(0);
    let next_ref = &next_idx;

    thread::scope(|s| {
        let mut handles = Vec::with_capacity(n_threads);
        for _ in 0..n_threads {
            handles.push(s.spawn(move || -> (HashSet<Vec<i8>>, DfsStats) {
                let mut local: HashSet<Vec<i8>> = HashSet::new();
                let mut stats = DfsStats::default();
                loop {
                    let i = next_ref.fetch_add(1, Ordering::Relaxed);
                    if i >= seeds.len() {
                        break;
                    }
                    let mut snake: Snake<ZZ> = Snake::from_slice_unsafe(&seeds[i]);
                    let mut record = hashset_recorder(&mut local);
                    rat_enum_step::<ZZ>(
                        &mut snake,
                        max_steps,
                        step,
                        &mut record,
                        &mut stats,
                        ops,
                        paranoid,
                        prunes,
                        None,
                    );
                }
                (local, stats)
            }));
        }
        let mut merged = closed_main;
        let mut total_stats = seed_stats;
        for h in handles {
            let (local, wstats) = h.join().expect("worker panic");
            merged.extend(local);
            total_stats.merge(&wstats);
        }
        (merged, total_stats)
    })
}
