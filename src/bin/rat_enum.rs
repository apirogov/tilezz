use std::collections::HashSet;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::thread;
use std::time::Instant;

use clap::{Parser, ValueEnum};

use tilezz::cyclotomic::*;
use tilezz::geom::rat::{lex_min_rot, Rat};
use tilezz::geom::snake::{Snake, Turtle};
use tilezz::stringmatch::repetition_factor;
use tilezz::vis::animation::render_gif;
use tilezz::vis::draw::{MarkerStyle, TileStyle};
use tilezz::vis::plotutils::P64;
use tilezz::vis::scene::{Color, Fill, Scene, Stroke, TextStyle, Viewport};

static VERBOSE: Mutex<bool> = Mutex::new(false);

// -------- DFS stats --------

#[derive(Default)]
struct DfsStats {
    closed: u64,
    intersected: u64,
    too_far: u64,
    recursed: u64,
    canonical_skip: u64,
}

impl DfsStats {
    fn total(&self) -> u64 {
        self.closed + self.intersected + self.too_far + self.recursed + self.canonical_skip
    }

    fn merge(&mut self, other: &DfsStats) {
        self.closed += other.closed;
        self.intersected += other.intersected;
        self.too_far += other.too_far;
        self.recursed += other.recursed;
        self.canonical_skip += other.canonical_skip;
    }
}

impl std::fmt::Display for DfsStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let t = self.total();
        writeln!(f, "DFS stats ({} total direction attempts):", t)?;
        writeln!(
            f,
            "  canonical_skip:  {:>10} ({:>5.1}%)",
            self.canonical_skip,
            100.0 * self.canonical_skip as f64 / t as f64
        )?;
        writeln!(
            f,
            "  intersected:    {:>10} ({:>5.1}%)",
            self.intersected,
            100.0 * self.intersected as f64 / t as f64
        )?;
        writeln!(
            f,
            "  closed:         {:>10} ({:>5.1}%)",
            self.closed,
            100.0 * self.closed as f64 / t as f64
        )?;
        writeln!(
            f,
            "  recursed:       {:>10} ({:>5.1}%)",
            self.recursed,
            100.0 * self.recursed as f64 / t as f64
        )?;
        write!(
            f,
            "  too_far:        {:>10} ({:>5.1}%)",
            self.too_far,
            100.0 * self.too_far as f64 / t as f64
        )
    }
}

// --------

/// Pair of canonical-check and output-mapping functions that parameterise
/// the DFS.  Both rotation-canonical and dihedral-canonical enumeration
/// share the same core walk; they differ only in which pair of functions
/// they supply.
///
/// * `is_canonical` — prefix prune applied *before* `Snake::add`.
///   Returns false when the extended walk cannot be the lex-min rotation
///   (or dihedral image) of any closure it could grow into.
/// * `canonicalize` — applied to the chirality-normalised canonical
///   rotation at closure, producing the key inserted into the result
///   `HashSet`.  For rotation-canonical this is the identity; for
///   dihedral-canonical it picks the lex-min over rotations and
///   reversed-rotations.
#[derive(Clone, Copy)]
struct CanonicalOps {
    is_canonical: fn(&[i8], i8) -> bool,
    canonicalize: fn(&[i8]) -> Vec<i8>,
}

/// Canonical-rotation prune for the open walk `prefix` virtually
/// extended by one more angle `new`. Lets the caller skip walks
/// that cannot be the lex-min cyclic rotation of any closure they
/// could grow into, *before* paying the cyclotomic intersect cost
/// of `Snake::add`.
///
/// Returns `false` when there's a rotation index `k > 0` such that
/// the rotation `[prefix++new][k..]` is strictly lex-less than
/// `[prefix++new][0..]` within the wrap-free comparable region; that
/// decision is permanent (no future angles can flip it) so the
/// eventual closed polygon's canonical rotation cannot start at
/// position 0 -- pruning is safe.
fn is_canonical_extended(prefix: &[i8], new: i8) -> bool {
    let base = prefix.len();
    let d = base + 1;
    let get = |i: usize| -> i8 {
        if i < base {
            prefix[i]
        } else {
            new
        }
    };
    for k in 1..d {
        for i in 0..(d - k) {
            match get(k + i).cmp(&get(i)) {
                std::cmp::Ordering::Less => return false,
                std::cmp::Ordering::Greater => break,
                std::cmp::Ordering::Equal => continue,
            }
        }
    }
    true
}

fn canonical_identity(seq: &[i8]) -> Vec<i8> {
    seq.to_vec()
}

/// Enumerate every simple polygon with boundary length up to
/// `max_steps` over the cyclotomic ring `ZZ`, in canonical-CCW form.
///
/// Depth-first walk: a single `Snake` is mutated in place via
/// [`Snake::add`] / [`Snake::pop`] across the recursion, so the inner
/// loop allocates only when a polygon actually closes (the canonical
/// [`Rat`] needs a fresh vec for hashing). Each polygon's `n` cyclic
/// walks collapse to a single canonical walk via the lex-min
/// rotation prune in [`is_canonical_extended`].
///
/// `step` restricts the DFS to directions that are multiples of `step`.
/// `step = 1` (default) walks every direction. `step = 2` on ZZ20 walks
/// only even-indexed directions, enumerating the ZZ10-equivalent
/// subset; `step = 2` on ZZ24 enumerates the ZZ12 subset; etc.
#[cfg(test)]
fn rat_enum<ZZ: IsRing>(max_steps: usize, step: i8) -> (Vec<Vec<i8>>, DfsStats) {
    rat_enum_with::<ZZ>(
        max_steps,
        step,
        CanonicalOps {
            is_canonical: is_canonical_extended,
            canonicalize: canonical_identity,
        },
        "enumeration",
        "",
    )
}

/// Dihedral-canonical variant of [`rat_enum`].
///
/// Two-stage design (see module-level comment above for rationale):
///
/// 1. [`is_dihedral_canonical_extended`] prunes walk prefixes whose
///    complement rotation is lex-smaller, reducing the search tree
///    (~2.3x speedup at ZZ12 n=10).
///
/// 2. [`dihedral_canonical`] at closure maps the chirality-normalized
///    canonical rotation to the lex-min dihedral form, so both
///    surviving members of each chiral pair hash to the same key.
///
/// Returns the same set as running [`rat_enum`] and quotienting by
/// dihedral equivalence.
#[cfg(test)]
fn rat_enum_dihedral<ZZ: IsRing>(max_steps: usize, step: i8) -> (Vec<Vec<i8>>, DfsStats) {
    rat_enum_with::<ZZ>(
        max_steps,
        step,
        CanonicalOps {
            is_canonical: is_dihedral_canonical_extended,
            canonicalize: dihedral_canonical,
        },
        "dihedral enumeration",
        "dihedral ",
    )
}

fn rat_enum_with<ZZ: IsRing>(
    max_steps: usize,
    step: i8,
    ops: CanonicalOps,
    label: &str,
    prefix: &str,
) -> (Vec<Vec<i8>>, DfsStats) {
    let mut result: HashSet<Vec<i8>> = HashSet::new();
    let mut snake: Snake<ZZ> = Snake::new();
    let mut stats = DfsStats::default();

    println!("-------- {label} started --------");
    rat_enum_step::<ZZ>(&mut snake, max_steps, step, &mut result, &mut stats, ops);
    println!(
        "-------- {label} completed --------\n{prefix}{} rats found",
        result.len()
    );

    let mut result: Vec<Vec<i8>> = result.into_iter().collect();
    result.sort_by_key(|x| x.len());
    (result, stats)
}

fn rat_enum_step<ZZ: IsRing>(
    snake: &mut Snake<ZZ>,
    max_steps: usize,
    step: i8,
    result: &mut HashSet<Vec<i8>>,
    stats: &mut DfsStats,
    ops: CanonicalOps,
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

        if !snake.add(direction) {
            stats.intersected += 1;
            continue;
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
            if result.insert(seq.clone()) {
                println!("RAT {seq:?}");
            }
        } else {
            stats.recursed += 1;
            rat_enum_step::<ZZ>(snake, max_steps, step, result, stats, ops);
        }
        snake.pop();
    }
}

// -------- dihedral-canonical enumeration --------
//
// Two-stage design: prefix pruning for speed + output dedup for
// correctness.  The prefix prune alone is NOT sufficient to produce
// exactly one representative per dihedral class.  Here is why.
//
// A polygon's dihedral images include all n rotations and all n
// reflections.  For a chiral pair {P, mirror(P)}, both have distinct
// rotation-canonical forms.  To produce dihedral classes, we must
// keep exactly one.
//
// Stage 1 — is_dihedral_canonical_extended:
//
//   Compares the walk prefix against complement rotations (negated
//   angle sequences).  At depth 0 this prunes all positive first
//   angles (-a_0 < a_0), cutting the root branching factor roughly in
//   half.  At deeper depths it prunes additional branches where a
//   complement rotation is lex-smaller.
//
//   This is a HEURISTIC prefix prune, not a dihedral guarantee.  The
//   check guards the raw walk prefix, but the DFS output at closure
//   goes through chirality normalization (.reversed() if chirality <
//   0) then canonical rotation (.canonical()).  These transformations
//   produce a different sequence than what the prefix check was
//   comparing, so both members of some chiral pairs still reach
//   closure as distinct canonical-rotation outputs.
//
// Stage 2 — dihedral_canonical at closure:
//
//   Maps the chirality-normalized canonical rotation to the lex-min
//   over all rotations AND reversed-rotations.  Both canonical
//   rotations of a chiral pair map to the same dihedral-canonical
//   form, so the HashSet deduplicates them.  Without this, the output
//   would contain duplicate chiral pairs — valid polygons, but not
//   dihedral-canonical.
//
// The speedup (~2.3x at ZZ12 n=10) comes from stage 1 reducing the
// search tree.  Stage 2 is cheap (O(n^2) per closed polygon) and
// handles the correctness gap left by stage 1.

/// Heuristic dihedral prefix prune.  Like [`is_canonical_extended`]
/// but also compares the walk prefix against complement rotations
/// (negated angle sequences).  Returns `false` when any complement
/// rotation of the eventual closure has a lex-smaller prefix than the
/// identity walk, indicating the walk is suboptimal under the
/// dihedral group.
///
/// This is sound (never prunes a walk that could produce the
/// dihedral-min output) but incomplete (does not eliminate all
/// chiral duplicates — see module-level comment above).  At depth 0
/// it prunes all positive first angles (`-a_0 < a_0`), roughly
/// halving the root branching factor.
///
/// Complement rotation `k` maps position `i` to `-a_{(k-i) mod n}`.
/// For `k <= d` (where `d = len(prefix) + 1`), positions `i = 0..k`
/// have both sides known.
fn is_dihedral_canonical_extended(prefix: &[i8], new: i8) -> bool {
    let base = prefix.len();
    let d = base + 1;
    let get = |i: usize| -> i8 {
        if i < base {
            prefix[i]
        } else {
            new
        }
    };

    // Rotation check (unchanged from is_canonical_extended).
    for k in 1..d {
        for i in 0..(d - k) {
            match get(k + i).cmp(&get(i)) {
                std::cmp::Ordering::Less => return false,
                std::cmp::Ordering::Greater => break,
                std::cmp::Ordering::Equal => continue,
            }
        }
    }

    // Complement-reflection check: for k = 0..d, compare -a_{k-i}
    // with a_i at positions i = 0..k.
    for k in 0..=d {
        for i in 0..=k {
            match (-get(k - i)).cmp(&get(i)) {
                std::cmp::Ordering::Less => return false,
                std::cmp::Ordering::Greater => break,
                std::cmp::Ordering::Equal => continue,
            }
        }
    }

    true
}

/// Correctness layer for the dihedral DFS.  Computes the lex-minimum
/// over all cyclic rotations and reversed rotations of `seq`.  Both
/// canonical rotations of a chiral pair map to the same value, so the
/// HashSet deduplicates them.  Without this step, the output would
/// contain both members of some chiral pairs
/// (see module-level comment above).
fn dihedral_canonical(seq: &[i8]) -> Vec<i8> {
    let n = seq.len();
    let mut best: Vec<i8> = seq.to_vec();
    let mut rot: Vec<i8> = seq.to_vec();
    for _ in 0..n {
        if rot < best {
            best = rot.clone();
        }
        let mut rev = rot.clone();
        rev.reverse();
        if rev < best {
            best = rev.clone();
        }
        rot.rotate_left(1);
    }
    best
}

// -------- multi-threaded enumeration --------

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
fn splitting_depth(n_threads: usize, branching: usize) -> usize {
    if n_threads <= 1 || branching <= 1 {
        return 0;
    }
    let target = (10 * n_threads) as f64;
    let depth = (target.ln() / (branching as f64).ln()).ceil() as usize;
    depth.max(1)
}

/// Output accumulators for [`collect_seeds`]: incomplete-walk prefixes
/// (seeds for worker threads), already-closed polygons, and DFS stats.
struct SeedGather<'a> {
    seeds: &'a mut Vec<Vec<i8>>,
    closed: &'a mut HashSet<Vec<i8>>,
    stats: &'a mut DfsStats,
}

/// Walk the existing DFS only down to `split_depth` and collect every
/// alive snake state (angle prefix) that survives the canonical-rotation
/// prune, the `Snake::add` self-intersect check, and the reachability
/// heuristic. These are the seeds the worker threads will pick up.
///
/// Polygons that already close at or above `split_depth` are recorded
/// directly into `gather.closed` -- they have no remaining work to delegate.
fn collect_seeds<ZZ: IsRing>(
    snake: &mut Snake<ZZ>,
    max_steps: usize,
    step: i8,
    split_depth: usize,
    gather: &mut SeedGather<'_>,
    ops: CanonicalOps,
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

        if !snake.add(direction) {
            gather.stats.intersected += 1;
            continue;
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
            if gather.closed.insert(seq.clone()) {
                println!("RAT {seq:?}");
            }
        } else {
            gather.stats.recursed += 1;
            if snake.angles().len() >= split_depth {
                gather.seeds.push(snake.angles().to_vec());
            } else {
                collect_seeds::<ZZ>(snake, max_steps, step, split_depth, gather, ops);
            }
        }
        snake.pop();
    }
}

/// Parallel variant of [`rat_enum`]: splits the DFS at `split_depth`
/// (selected via [`splitting_depth`]), then hands the resulting alive
/// prefixes out to `n_threads` worker threads via a shared atomic
/// counter. Each worker keeps its own `HashSet` and the main thread
/// merges the per-worker sets at the end.
fn rat_enum_parallel<ZZ: IsRing>(
    max_steps: usize,
    step: i8,
    n_threads: usize,
    ops: CanonicalOps,
    label: &str,
    prefix: &str,
) -> (Vec<Vec<i8>>, DfsStats) {
    let hm1 = (ZZ::hturn() as usize).saturating_sub(1);
    let branching = 2 * (hm1 / step.max(1) as usize) + 1;
    let split_depth = splitting_depth(n_threads, branching);

    println!("-------- {label} started --------");
    println!("parallel: n_threads={n_threads} branching={branching} split_depth={split_depth}");

    let mut closed_main: HashSet<Vec<i8>> = HashSet::new();
    let mut seeds: Vec<Vec<i8>> = Vec::new();
    let mut seed_stats = DfsStats::default();
    {
        let mut snake: Snake<ZZ> = Snake::new();
        let mut gather = SeedGather {
            seeds: &mut seeds,
            closed: &mut closed_main,
            stats: &mut seed_stats,
        };
        collect_seeds::<ZZ>(&mut snake, max_steps, step, split_depth, &mut gather, ops);
    }
    println!("parallel: {} seed states collected", seeds.len());

    let next_idx = AtomicUsize::new(0);
    let seeds_ref: &[Vec<i8>] = &seeds;
    let next_ref = &next_idx;

    let (merged, worker_stats): (HashSet<Vec<i8>>, DfsStats) = thread::scope(|s| {
        let mut handles = Vec::with_capacity(n_threads);
        for _ in 0..n_threads {
            handles.push(s.spawn(move || -> (HashSet<Vec<i8>>, DfsStats) {
                let mut local: HashSet<Vec<i8>> = HashSet::new();
                let mut stats = DfsStats::default();
                loop {
                    let i = next_ref.fetch_add(1, Ordering::Relaxed);
                    if i >= seeds_ref.len() {
                        break;
                    }
                    let mut snake: Snake<ZZ> = Snake::from_slice_unsafe(&seeds_ref[i]);
                    rat_enum_step::<ZZ>(&mut snake, max_steps, step, &mut local, &mut stats, ops);
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
    });

    println!(
        "-------- {label} completed --------\n{prefix}{} rats found",
        merged.len()
    );

    let mut result: Vec<Vec<i8>> = merged.into_iter().collect();
    result.sort_by_key(|x| x.len());
    (result, worker_stats)
}

fn polygons<ZZ: IsRing>(rats: Vec<Vec<i8>>) -> Vec<Vec<P64>> {
    rats.into_iter()
        .map(|seq| Rat::<ZZ>::from_slice_unchecked(&seq).to_polyline_f64(Turtle::default()))
        .collect()
}

type EnumResult = (Vec<Vec<i8>>, DfsStats);

fn enumerate_dispatch<ZZ: IsRing>(
    max_steps: usize,
    step: i8,
    n_threads: usize,
    dihedral: bool,
) -> EnumResult {
    let ops = if dihedral {
        CanonicalOps {
            is_canonical: is_dihedral_canonical_extended,
            canonicalize: dihedral_canonical,
        }
    } else {
        CanonicalOps {
            is_canonical: is_canonical_extended,
            canonicalize: canonical_identity,
        }
    };

    let label = if dihedral {
        "dihedral enumeration"
    } else {
        "enumeration"
    };
    let prefix = if dihedral { "dihedral " } else { "" };

    if n_threads <= 1 {
        rat_enum_with::<ZZ>(max_steps, step, ops, label, prefix)
    } else {
        rat_enum_parallel::<ZZ>(max_steps, step, n_threads, ops, label, prefix)
    }
}

fn run_rat_enum_polylines(
    ring: u8,
    max_steps: usize,
    step: i8,
    n_threads: usize,
    dihedral: bool,
) -> Vec<Vec<P64>> {
    match ring {
        4 => polygons::<ZZ4>(enumerate_dispatch::<ZZ4>(max_steps, step, n_threads, dihedral).0),
        8 => polygons::<ZZ8>(enumerate_dispatch::<ZZ8>(max_steps, step, n_threads, dihedral).0),
        10 => polygons::<ZZ10>(enumerate_dispatch::<ZZ10>(max_steps, step, n_threads, dihedral).0),
        12 => polygons::<ZZ12>(enumerate_dispatch::<ZZ12>(max_steps, step, n_threads, dihedral).0),
        16 => polygons::<ZZ16>(enumerate_dispatch::<ZZ16>(max_steps, step, n_threads, dihedral).0),
        20 => polygons::<ZZ20>(enumerate_dispatch::<ZZ20>(max_steps, step, n_threads, dihedral).0),
        24 => polygons::<ZZ24>(enumerate_dispatch::<ZZ24>(max_steps, step, n_threads, dihedral).0),
        32 => polygons::<ZZ32>(enumerate_dispatch::<ZZ32>(max_steps, step, n_threads, dihedral).0),
        60 => polygons::<ZZ60>(enumerate_dispatch::<ZZ60>(max_steps, step, n_threads, dihedral).0),
        _ => panic!("invalid ring selected"),
    }
}

/// True iff this canonical sequence equals the canonical form of its
/// mirror image. The mirror polygon's CCW angle sequence is just
/// `reverse(seq)` -- traversing the mirrored vertices in CCW order (in
/// original coordinates) hits the original vertices in reverse order,
/// and each turn keeps its sign because the mirror's orientation flip
/// and the CCW-vs-CW flip cancel. Achiral polygons survive a chirality
/// quotient as a single class.
fn is_achiral(canonical: &[i8]) -> bool {
    let mut reflected: Vec<i8> = canonical.iter().rev().copied().collect();
    let offset = lex_min_rot(&reflected);
    reflected.rotate_left(offset);
    reflected == canonical
}

/// Print per-boundary-length statistics over the enumerated canonical
/// sequences: total count, achiral count, and a histogram of cyclic
/// rotational symmetry orders.
fn print_stats(rats: &[Vec<i8>]) {
    use std::collections::BTreeMap;
    // length -> (total, achiral, rep_factor -> count)
    let mut per_len: BTreeMap<usize, (usize, usize, BTreeMap<usize, usize>)> = BTreeMap::new();
    for seq in rats {
        let n = seq.len();
        let entry = per_len.entry(n).or_default();
        entry.0 += 1;
        if is_achiral(seq) {
            entry.1 += 1;
        }
        let rf = repetition_factor(seq);
        *entry.2.entry(rf).or_insert(0) += 1;
    }
    println!("statistics by boundary length:");
    println!("  length |  total | achiral (%) | rotational symmetry histogram (rep_factor: count)");
    println!("  -------+--------+-------------+--------------------------------------------------");
    for (n, (total, achiral, hist)) in &per_len {
        let pct = if *total > 0 {
            100.0 * (*achiral as f64) / (*total as f64)
        } else {
            0.0
        };
        let mut hist_str = String::new();
        for (rf, cnt) in hist {
            if !hist_str.is_empty() {
                hist_str.push_str("  ");
            }
            hist_str.push_str(&format!("x{rf}={cnt}"));
        }
        println!("  {n:>6} | {total:>6} | {achiral:>5} ({pct:>5.1}%) | {hist_str}");
    }
}

fn run_rat_enum_seqs(
    ring: u8,
    max_steps: usize,
    step: i8,
    n_threads: usize,
    dihedral: bool,
) -> EnumResult {
    let f: fn(usize, i8, usize, bool) -> EnumResult = match ring {
        4 => enumerate_dispatch::<ZZ4>,
        8 => enumerate_dispatch::<ZZ8>,
        10 => enumerate_dispatch::<ZZ10>,
        12 => enumerate_dispatch::<ZZ12>,
        16 => enumerate_dispatch::<ZZ16>,
        20 => enumerate_dispatch::<ZZ20>,
        24 => enumerate_dispatch::<ZZ24>,
        32 => enumerate_dispatch::<ZZ32>,
        60 => enumerate_dispatch::<ZZ60>,
        _ => panic!("invalid ring selected"),
    };
    f(max_steps, step, n_threads, dihedral)
}

// --------

#[derive(Copy, Clone, Debug, ValueEnum)]
enum Mode {
    /// Enumerate and render output (GIF)
    Render,
    /// Enumerate only and report elapsed time
    Bench,
}

#[derive(Parser, Debug)]
#[command(version, about = "Compute all simple polygons over a ring with a maximal boundary length", long_about = None)]
struct Cli {
    #[arg(short = 'r', long)]
    ring: u8,

    #[arg(short = 'n', long)]
    max_steps: usize,

    #[arg(long, value_enum, default_value_t = Mode::Render)]
    mode: Mode,

    #[arg(short = 'o', long, help = "Output GIF filename (render mode only)")]
    filename: Option<String>,

    #[arg(short, long)]
    verbose: bool,

    /// In Bench mode, write a flamegraph SVG to this path (requires
    /// the `pprof` cargo feature).
    #[arg(long)]
    profile: Option<String>,

    /// In Bench mode, print a per-boundary-length breakdown of total /
    /// achiral / rotational-symmetry histogram after enumeration.
    #[arg(long)]
    stats: bool,

    /// Number of worker threads for the DFS.
    ///
    /// `1` (default) runs the original single-threaded path
    /// unchanged. `0` resolves to `num_cpus::get()`. Any other
    /// positive value is used as-is, capped at `num_cpus::get()`.
    #[arg(long, default_value_t = 1)]
    threads: usize,

    /// Restrict the DFS to directions that are multiples of `step`.
    ///
    /// `1` (default) walks every direction in `(-hturn+1)..hturn`.
    /// `step = 2` on ZZ20 walks only even-indexed directions, which
    /// enumerates exactly the ZZ10-equivalent polygons (b=9 vs 19).
    /// Similarly `step = 2` on ZZ24 -> ZZ12 subset, `step = 3` on
    /// ZZ24 -> ZZ8 subset, etc.
    #[arg(long, default_value_t = 1)]
    step: i8,

    /// Use dihedral-canonical enumeration: outputs one representative
    /// per chiral pair (lex-min over rotations and reflections).
    /// Faster than the default rotation-canonical DFS because
    /// complement-reflection pruning halves the search space at
    /// the root.
    #[arg(long)]
    dihedral: bool,
}

/// Resolve the `--threads` CLI value to an actual worker count.
/// `0` -> `num_cpus::get()`; otherwise cap at `num_cpus::get()`.
fn resolve_n_threads(requested: usize) -> usize {
    let max = num_cpus::get().max(1);
    if requested == 0 {
        max
    } else {
        requested.min(max)
    }
}

fn main() {
    let cli = Cli::parse();
    if cli.verbose {
        let mut verbose = VERBOSE.lock().unwrap();
        *verbose = true;
    }

    let n_threads = resolve_n_threads(cli.threads);

    match cli.mode {
        Mode::Bench => {
            let profile = tilezz::util::profile::ProfileGuard::start(cli.profile.as_deref());

            let t0 = Instant::now();
            let (rats, stats) =
                run_rat_enum_seqs(cli.ring, cli.max_steps, cli.step, n_threads, cli.dihedral);
            let dt = t0.elapsed();

            // Use the result so it can't be trivially optimized away.
            let total_boundary_len: usize = rats.iter().map(|s| s.len()).sum();

            println!(
                "benchmark: ring={} step={} max_steps={} -> {} unique rats (total boundary len={}) in {:?}",
                cli.ring,
                cli.step,
                cli.max_steps,
                rats.len(),
                total_boundary_len,
                dt
            );

            println!("{stats}");

            profile.finish();

            if cli.stats {
                print_stats(&rats);
            }
        }
        Mode::Render => {
            let rats: Vec<Vec<P64>> =
                run_rat_enum_polylines(cli.ring, cli.max_steps, cli.step, n_threads, cli.dihedral);

            let Some(filename) = cli.filename else {
                return;
            };

            // Symmetric square viewport sized to fit every polygon
            // centered at the origin. Each polygon goes through
            // Rat::to_polyline_f64 with the default Turtle (origin,
            // direction 0), so we just need the max |x|/|y| across
            // all polylines.
            let mut max_abs: f64 = 1.0;
            for poly in &rats {
                for (x, y) in poly {
                    max_abs = max_abs.max(x.abs()).max(y.abs());
                }
            }
            let pad = 0.15 * max_abs;
            let r = max_abs + pad;
            let bounds = ((-r, -r), (r, r));

            // Style: filled yellow tile with thin black border,
            // red vertex markers labelled with their index.
            let style = TileStyle::filled(
                Fill::solid(Color::YELLOW.with_alpha(80)),
                Stroke::solid(Color::BLACK, 0.01 * r),
            )
            .with_vertex_marker(MarkerStyle::filled_circle(0.06 * r, Color::RED))
            .with_vertex_labels(TextStyle::new(0.04 * r, Color::WHITE).bold());

            let frames: Vec<Scene> = rats
                .iter()
                .enumerate()
                .map(|(ix, tile)| {
                    println!("rendering frame {}/{}", ix + 1, rats.len());
                    let mut scene = Scene::new().with_background(Color::WHITE);
                    scene.draw_tile(tile, &style);
                    scene
                })
                .collect();

            let w = 500u32;
            let vp = Viewport::square_for(w, bounds, 8);
            let gif_bytes = render_gif(&frames, &vp, 500).expect("render GIF");
            std::fs::write(&filename, gif_bytes).expect("write GIF");
            println!("wrote {filename}");
        }
    }
}

#[cfg(test)]
mod dihedral_tests {
    use super::*;
    use tilezz::cyclotomic::ZZ12;

    /// Cross-validate the dihedral DFS variant against the original
    /// rotation-canonical DFS: the original output quotiented by
    /// dihedral_canonical must equal the dihedral variant's output.
    #[test]
    fn test_dihedral_enum_matches_dfs_quotient() {
        let (all_rats, _stats) = super::rat_enum::<ZZ12>(9, 1);
        let mut expected: HashSet<Vec<i8>> = HashSet::new();
        for seq in &all_rats {
            expected.insert(super::dihedral_canonical(seq));
        }

        let (dihedral_rats, _stats) = super::rat_enum_dihedral::<ZZ12>(9, 1);
        let actual: HashSet<Vec<i8>> = dihedral_rats.into_iter().collect();

        if expected != actual {
            let missing: Vec<_> = expected.difference(&actual).collect();
            let extra: Vec<_> = actual.difference(&expected).collect();
            if !missing.is_empty() {
                eprintln!("MISSING from dihedral DFS:");
                for s in &missing {
                    eprintln!("  {s:?}");
                }
            }
            if !extra.is_empty() {
                eprintln!("EXTRA in dihedral DFS:");
                for s in &extra {
                    eprintln!("  {s:?}");
                }
            }
        }
        assert_eq!(
            expected.len(),
            actual.len(),
            "DFS quotient has {}, dihedral DFS found {}",
            expected.len(),
            actual.len()
        );
        assert_eq!(
            expected, actual,
            "dihedral DFS output differs from DFS quotient"
        );
        eprintln!(
            "OK — {} dihedral classes match (from {} rotation-canonical rats)",
            actual.len(),
            all_rats.len()
        );
    }
}
