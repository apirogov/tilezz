//! Depth-first walk over the cyclotomic unit-direction graph.
//!
//! Three layered entry points:
//!
//! - [`rat_enum_with`]: single-threaded top-level. Constructs the
//!   initial snake, recursive entry into [`rat_enum_step`], collects
//!   the closed canonical sequences into a `HashSet`.
//!
//! - [`rat_enum_step`]: the recursive worker. Holds the snake as
//!   `&mut`, tries every direction (subject to `step`), applies all
//!   prunes in order (canonical -> reachability -> modular ->
//!   closure-key -> intersect), recurses or records.
//!
//! - [`collect_seeds`]: the parallel variant's seed-collection
//!   prepass. Walks down to `split_depth` and accumulates alive
//!   prefixes into [`SeedGather::seeds`]; polygons that close before
//!   reaching the split are reported via [`SeedGather::record_closed`].
//!
//! Both [`rat_enum_step`] and [`SeedGather`] take a `&mut dyn
//! FnMut(&[i8])` to receive closed canonical sequences, so the same
//! engine drives the in-memory ([`hashset_recorder`]) and streaming
//! ([`super::stream`]) sinks.
//!
//! All three respect the optional prunes through the [`prune::Prunes`]
//! parameter; see [`super::prune`] for what each prune does.

use std::collections::HashSet;
use std::sync::atomic::{AtomicBool, Ordering};

use crate::cyclotomic::{IsRing, Units};
use crate::geom::rat::Rat;
use crate::geom::snake::Snake;
use crate::rat_enum::canonical::CanonicalOps;
use crate::rat_enum::prune::Prunes;
use crate::rat_enum::stats::DfsStats;

/// When `true`, the DFS streams each newly-discovered rat to stdout
/// as `RAT [...]` on the fly. This is the protocol used by `--seed`
/// and `--mode bench`; modes like `--mode dafsa` that write a single
/// binary artifact set this to `false` before enumerating to avoid
/// millions of useless stdout lines.
pub static STREAM_RAT_LINES: AtomicBool = AtomicBool::new(true);

/// Single-threaded enumeration entry point. Initializes an empty
/// snake, runs the recursive DFS, returns the `(sorted-by-length)
/// list, stats)` pair.
///
/// `label` is the header to print at start/end (e.g. `"free
/// enumeration"`); `prefix` is the prefix for the closing summary
/// line (e.g. `"free "`).
pub fn rat_enum_with<ZZ: IsRing>(
    max_steps: usize,
    step: i8,
    ops: CanonicalOps,
    label: &str,
    prefix: &str,
    paranoid: bool,
    prunes: &Prunes,
) -> (Vec<Vec<i8>>, DfsStats) {
    let mut result: HashSet<Vec<i8>> = HashSet::new();
    let mut snake: Snake<ZZ> = Snake::new();
    let mut stats = DfsStats::default();

    println!("-------- {label} started --------");
    if paranoid {
        println!("paranoid: per-step fresh-snake cross-check enabled");
    }
    {
        let mut record = hashset_recorder(&mut result);
        rat_enum_step::<ZZ>(
            &mut snake,
            max_steps,
            step,
            &mut record,
            &mut stats,
            ops,
            paranoid,
            prunes,
        );
    }
    println!(
        "-------- {label} completed --------\n{prefix}{} rats found",
        result.len()
    );

    let mut result: Vec<Vec<i8>> = result.into_iter().collect();
    result.sort_by_key(|x| x.len());
    (result, stats)
}

/// Build the canonical "record into a `HashSet`" callback used by
/// the in-memory enumeration paths. Inserts each new sequence and,
/// if [`STREAM_RAT_LINES`] is set, emits a `RAT [...]` line to
/// stdout. Lifted here so every in-memory caller -- single-thread,
/// parallel workers, seed-walk closures -- uses the same recipe.
pub fn hashset_recorder<'a>(set: &'a mut HashSet<Vec<i8>>) -> impl FnMut(&[i8]) + 'a {
    move |seq: &[i8]| {
        if set.insert(seq.to_vec()) && STREAM_RAT_LINES.load(Ordering::Relaxed) {
            println!("RAT {seq:?}");
        }
    }
}

/// Recursive DFS worker. Mutates `snake` in place across recursion
/// (via `add` / `pop`). On every closure the canonical sequence is
/// passed to `record`; the callback decides where it goes (e.g.
/// `HashSet::insert` for the in-memory path, or a sort-buffer push
/// for the streaming path -- see [`super::stream`]). The callback
/// returns nothing; `STREAM_RAT_LINES` printing now lives inside the
/// callbacks the in-memory callers wrap their `HashSet` with.
#[allow(clippy::too_many_arguments)]
pub fn rat_enum_step<ZZ: IsRing>(
    snake: &mut Snake<ZZ>,
    max_steps: usize,
    step: i8,
    record: &mut dyn FnMut(&[i8]),
    stats: &mut DfsStats,
    ops: CanonicalOps,
    paranoid: bool,
    prunes: &Prunes,
) {
    let depth = snake.angles().len();
    if depth >= max_steps {
        return;
    }
    let remaining = (max_steps - depth) as i64;

    for direction in ((-ZZ::hturn() + 1)..ZZ::hturn()).rev() {
        if direction.rem_euclid(step) != 0 {
            continue;
        }
        // Canonical prune before the radius and geometry checks so
        // that rejected branches pay nothing.
        if !(ops.is_canonical)(snake.angles(), direction) {
            stats.canonical_skip += 1;
            continue;
        }

        // Early reachability prune: compute the next head position
        // before paying the cyclotomic intersect cost of Snake::add.
        // If the new point is too far from the origin to ever close,
        // skip this direction entirely.
        let new_pt = snake.offset()
            + <ZZ as Units>::unit(snake.direction()) * <ZZ as Units>::unit(direction);
        if !new_pt.is_zero() && !new_pt.within_radius(remaining) {
            stats.too_far += 1;
            continue;
        }

        // Optional modular reachability prune (set via `--mod-prune`).
        // After taking this direction, the snake will be at `new_pt`
        // with `remaining - 1` directions still to add for closure.
        // If no sum of (remaining-1) unit vectors equals `-new_pt`
        // (modulo any active modulus), closure is impossible.
        let remaining_after = (remaining as usize).saturating_sub(1);
        if let Some(mp) = prunes.mod_prune.as_deref() {
            if !mp.allows_closure(new_pt.int_coeffs_slice(), remaining_after) {
                stats.mod_skip += 1;
                continue;
            }
        }
        // Optional closure-key prune (set via `--closure-key-prune`).
        // Only fires when `remaining_after <= max_l` (otherwise the
        // tabulated suffix lengths aren't enough to cover the
        // closing range, and pruning would be unsound).
        if let Some(ck) = prunes.closure_key_prune.as_deref() {
            if remaining_after <= ck.max_l {
                let turn = ZZ::turn();
                let new_facing = (snake.direction() + direction).rem_euclid(turn);
                let neg_facing = (-new_facing).rem_euclid(turn);
                // target suffix endpoint = -unit(-new_facing) * new_pt.
                let target: ZZ = -(<ZZ as Units>::unit(neg_facing) * new_pt);
                let key = (target.int_coeffs_slice().to_vec(), neg_facing);
                if !ck.keys.contains(&key) {
                    stats.closure_key_skip += 1;
                    continue;
                }
            }
        }
        if !snake.add(direction) {
            stats.intersected += 1;
            continue;
        }
        if paranoid {
            // After each successful add(), replay the entire current
            // angle prefix in a fresh Snake. Any disagreement means
            // the stateful incremental check accepts a prefix the
            // from-scratch check rejects -- a Snake bug.
            let angles = snake.angles().to_vec();
            let mut fresh: Snake<ZZ> = Snake::new();
            for (i, &a) in angles.iter().enumerate() {
                assert!(
                    fresh.add(a),
                    "stateful snake accepted full prefix {:?} but fresh snake \
                     rejected angle {} at step {}",
                    &angles,
                    a,
                    i
                );
            }
            assert_eq!(
                fresh.is_closed(),
                snake.is_closed(),
                "fresh snake disagrees on is_closed for {:?}",
                &angles
            );
        }
        if snake.is_closed() {
            stats.closed += 1;
            let r = {
                let tmp = Rat::from_unchecked(snake);
                if tmp.chirality() > 0 {
                    tmp
                } else {
                    tmp.reversed()
                }
                .canonical()
            };
            let seq = (ops.canonicalize)(r.seq());
            record(&seq);
        } else {
            stats.recursed += 1;
            rat_enum_step::<ZZ>(snake, max_steps, step, record, stats, ops, paranoid, prunes);
        }
        snake.pop();
    }
}

/// Output accumulators for [`collect_seeds`]: incomplete-walk
/// prefixes (seeds for worker threads), a callback to receive
/// already-closed polygons that closed during the seed walk, and DFS
/// stats. The callback mirrors `rat_enum_step`'s `record` argument
/// so streaming sinks can plug in for both layers.
pub struct SeedGather<'a> {
    pub seeds: &'a mut Vec<Vec<i8>>,
    pub record_closed: &'a mut dyn FnMut(&[i8]),
    pub stats: &'a mut DfsStats,
}

/// Walk the existing DFS only down to `split_depth` and collect every
/// alive snake state (angle prefix) that survives the canonical-rotation
/// prune, the `Snake::add` self-intersect check, and the reachability
/// heuristic. These are the seeds the worker threads will pick up.
///
/// Polygons that already close at or above `split_depth` are recorded
/// directly into `gather.closed` -- they have no remaining work to delegate.
#[allow(clippy::too_many_arguments)]
pub fn collect_seeds<ZZ: IsRing>(
    snake: &mut Snake<ZZ>,
    max_steps: usize,
    step: i8,
    split_depth: usize,
    gather: &mut SeedGather<'_>,
    ops: CanonicalOps,
    paranoid: bool,
    prunes: &Prunes,
) {
    let depth = snake.angles().len();
    if depth >= max_steps {
        return;
    }
    let remaining = (max_steps - depth) as i64;

    for direction in ((-ZZ::hturn() + 1)..ZZ::hturn()).rev() {
        if direction.rem_euclid(step) != 0 {
            continue;
        }
        if !(ops.is_canonical)(snake.angles(), direction) {
            gather.stats.canonical_skip += 1;
            continue;
        }

        // Early reachability prune (see `rat_enum_step`).
        let new_pt = snake.offset()
            + <ZZ as Units>::unit(snake.direction()) * <ZZ as Units>::unit(direction);
        if !new_pt.is_zero() && !new_pt.within_radius(remaining) {
            gather.stats.too_far += 1;
            continue;
        }

        // Modular reachability prune (see `rat_enum_step`).
        let remaining_after = (remaining as usize).saturating_sub(1);
        if let Some(mp) = prunes.mod_prune.as_deref() {
            if !mp.allows_closure(new_pt.int_coeffs_slice(), remaining_after) {
                gather.stats.mod_skip += 1;
                continue;
            }
        }
        if let Some(ck) = prunes.closure_key_prune.as_deref() {
            if remaining_after <= ck.max_l {
                let turn = ZZ::turn();
                let new_facing = (snake.direction() + direction).rem_euclid(turn);
                let neg_facing = (-new_facing).rem_euclid(turn);
                let target: ZZ = -(<ZZ as Units>::unit(neg_facing) * new_pt);
                let key = (target.int_coeffs_slice().to_vec(), neg_facing);
                if !ck.keys.contains(&key) {
                    gather.stats.closure_key_skip += 1;
                    continue;
                }
            }
        }
        if !snake.add(direction) {
            gather.stats.intersected += 1;
            continue;
        }
        if paranoid {
            let angles = snake.angles().to_vec();
            let mut fresh: Snake<ZZ> = Snake::new();
            for (i, &a) in angles.iter().enumerate() {
                assert!(
                    fresh.add(a),
                    "stateful snake accepted full prefix {:?} but fresh snake \
                     rejected angle {} at step {}",
                    &angles,
                    a,
                    i
                );
            }
            assert_eq!(
                fresh.is_closed(),
                snake.is_closed(),
                "fresh snake disagrees on is_closed for {:?}",
                &angles
            );
        }
        if snake.is_closed() {
            gather.stats.closed += 1;
            let r = {
                let tmp = Rat::from_unchecked(snake);
                if tmp.chirality() > 0 {
                    tmp
                } else {
                    tmp.reversed()
                }
                .canonical()
            };
            let seq = (ops.canonicalize)(r.seq());
            (gather.record_closed)(&seq);
        } else {
            gather.stats.recursed += 1;
            if snake.angles().len() >= split_depth {
                gather.seeds.push(snake.angles().to_vec());
            } else {
                collect_seeds::<ZZ>(
                    snake,
                    max_steps,
                    step,
                    split_depth,
                    gather,
                    ops,
                    paranoid,
                    prunes,
                );
            }
        }
        snake.pop();
    }
}
