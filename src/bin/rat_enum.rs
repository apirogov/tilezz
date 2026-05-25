use std::collections::HashSet;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::thread;
use std::time::Instant;

use clap::{Parser, ValueEnum};

use tilezz::cyclotomic::*;
use tilezz::intgeom::rat::{lex_min_rot, Rat};
use tilezz::intgeom::snake::{Snake, Turtle};
use tilezz::stringmatch::repetition_factor;
use tilezz::vis::animation::render_gif;
use tilezz::vis::draw::{MarkerStyle, TileStyle};
use tilezz::vis::plotutils::P64;
use tilezz::vis::scene::{Color, Fill, Scene, Stroke, TextStyle, Viewport};

static VERBOSE: Mutex<bool> = Mutex::new(false);

// --------

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
        if i < base { prefix[i] } else { new }
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

/// Enumerate every simple polygon with boundary length up to
/// `max_steps` over the cyclotomic ring `ZZ`, in canonical-CCW form.
///
/// Depth-first walk: a single `Snake` is mutated in place via
/// [`Snake::add`] / [`Snake::pop`] across the recursion, so the inner
/// loop allocates only when a polygon actually closes (the canonical
/// [`Rat`] needs a fresh vec for hashing). Each polygon's `n` cyclic
/// walks collapse to a single canonical walk via the lex-min
/// rotation prune in [`is_canonical_extended`].
pub fn rat_enum<ZZ: ZZType + Units + WithinRadius>(max_steps: usize) -> Vec<Vec<i8>> {
    let mut result: HashSet<Vec<i8>> = HashSet::new();
    let mut snake: Snake<ZZ> = Snake::new();

    println!("-------- enumeration started --------");
    rat_enum_step::<ZZ>(&mut snake, max_steps, &mut result);
    println!(
        "-------- enumeration completed --------\n{} rats found",
        result.len()
    );

    let mut result: Vec<Vec<i8>> = result.into_iter().collect();
    result.sort_by_key(|x| x.len());
    result
}

fn rat_enum_step<ZZ: ZZType + Units + WithinRadius>(
    snake: &mut Snake<ZZ>,
    max_steps: usize,
    result: &mut HashSet<Vec<i8>>,
) {
    let depth = snake.angles().len();
    if depth >= max_steps {
        return;
    }
    // Heuristic: the snake's current head must lie within `remaining`
    // unit steps of the origin, else the path can never close in time.
    let remaining = (max_steps - depth) as i64;

    for direction in ((-ZZ::hturn() + 1)..ZZ::hturn()).rev() {
        // Canonical-rotation prune (see `is_canonical_extended`):
        // pre-check *before* calling `Snake::add` so that
        // canonical-rejected branches don't pay the cyclotomic
        // intersect cost.
        if !is_canonical_extended(snake.angles(), direction) {
            continue;
        }
        if !snake.add(direction) {
            continue;
        }
        if snake.is_closed() {
            // polygon completed -> record canonical ccw description.
            let r = {
                let tmp = Rat::from_unchecked(snake);
                if tmp.chirality() > 0 {
                    tmp
                } else {
                    tmp.reversed()
                }
                .canonical()
            };
            let seq = r.seq().to_vec();
            if result.insert(seq.clone()) {
                println!("RAT {seq:?}");
            }
        } else if snake.head().pos.within_radius(remaining) {
            // Still reachable -- recurse to extend further.
            rat_enum_step::<ZZ>(snake, max_steps, result);
        }
        // Backtrack: restore the snake to the pre-add state.
        snake.pop();
    }
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

/// Walk the existing DFS only down to `split_depth` and collect every
/// alive snake state (angle prefix) that survives the canonical-rotation
/// prune, the `Snake::add` self-intersect check, and the reachability
/// heuristic. These are the seeds the worker threads will pick up.
///
/// Polygons that already close at or above `split_depth` are recorded
/// directly into `closed` -- they have no remaining work to delegate.
fn collect_seeds<ZZ: ZZType + Units + WithinRadius>(
    snake: &mut Snake<ZZ>,
    max_steps: usize,
    split_depth: usize,
    seeds: &mut Vec<Vec<i8>>,
    closed: &mut HashSet<Vec<i8>>,
) {
    let depth = snake.angles().len();
    if depth >= max_steps {
        return;
    }
    let remaining = (max_steps - depth) as i64;

    for direction in ((-ZZ::hturn() + 1)..ZZ::hturn()).rev() {
        if !is_canonical_extended(snake.angles(), direction) {
            continue;
        }
        if !snake.add(direction) {
            continue;
        }
        if snake.is_closed() {
            let r = {
                let tmp = Rat::from_unchecked(snake);
                if tmp.chirality() > 0 {
                    tmp
                } else {
                    tmp.reversed()
                }
                .canonical()
            };
            let seq = r.seq().to_vec();
            if closed.insert(seq.clone()) {
                println!("RAT {seq:?}");
            }
        } else if snake.head().pos.within_radius(remaining) {
            if snake.angles().len() >= split_depth {
                // Reached splitting horizon -- emit this prefix as a
                // seed for a worker thread to expand.
                seeds.push(snake.angles().to_vec());
            } else {
                collect_seeds::<ZZ>(snake, max_steps, split_depth, seeds, closed);
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
pub fn rat_enum_parallel<ZZ: ZZType + Units + WithinRadius>(
    max_steps: usize,
    n_threads: usize,
) -> Vec<Vec<i8>> {
    // Branching = number of candidate directions per DFS level,
    // i.e. `(-hturn + 1)..hturn` which has length `2*hturn - 1`.
    let branching = (2 * ZZ::hturn() as usize).saturating_sub(1);
    let split_depth = splitting_depth(n_threads, branching);

    println!("-------- enumeration started --------");
    println!(
        "parallel: n_threads={n_threads} branching={branching} split_depth={split_depth}"
    );

    let mut closed_main: HashSet<Vec<i8>> = HashSet::new();
    let mut seeds: Vec<Vec<i8>> = Vec::new();
    {
        let mut snake: Snake<ZZ> = Snake::new();
        collect_seeds::<ZZ>(
            &mut snake,
            max_steps,
            split_depth,
            &mut seeds,
            &mut closed_main,
        );
    }
    println!("parallel: {} seed states collected", seeds.len());

    let next_idx = AtomicUsize::new(0);
    let seeds_ref: &[Vec<i8>] = &seeds;
    let next_ref = &next_idx;

    let merged: HashSet<Vec<i8>> = thread::scope(|s| {
        let mut handles = Vec::with_capacity(n_threads);
        for _ in 0..n_threads {
            handles.push(s.spawn(move || -> HashSet<Vec<i8>> {
                let mut local: HashSet<Vec<i8>> = HashSet::new();
                loop {
                    let i = next_ref.fetch_add(1, Ordering::Relaxed);
                    if i >= seeds_ref.len() {
                        break;
                    }
                    // Rebuild snake from this trusted (already self-
                    // intersect-checked) prefix and resume DFS.
                    let mut snake: Snake<ZZ> = Snake::from_slice_unsafe(&seeds_ref[i]);
                    rat_enum_step::<ZZ>(&mut snake, max_steps, &mut local);
                }
                local
            }));
        }
        let mut merged = closed_main;
        for h in handles {
            let local = h.join().expect("worker panic");
            merged.extend(local);
        }
        merged
    });

    println!(
        "-------- enumeration completed --------\n{} rats found",
        merged.len()
    );

    let mut result: Vec<Vec<i8>> = merged.into_iter().collect();
    result.sort_by_key(|x| x.len());
    result
}

fn polygons<ZZ: ZZType + Units>(rats: Vec<Vec<i8>>) -> Vec<Vec<P64>> {
    rats.into_iter()
        .map(|seq| Rat::<ZZ>::from_slice_unchecked(&seq).to_polyline_f64(Turtle::default()))
        .collect()
}

fn enumerate_dispatch<ZZ: ZZType + Units + WithinRadius>(
    max_steps: usize,
    n_threads: usize,
) -> Vec<Vec<i8>> {
    if n_threads <= 1 {
        rat_enum::<ZZ>(max_steps)
    } else {
        rat_enum_parallel::<ZZ>(max_steps, n_threads)
    }
}

fn run_rat_enum_polylines(ring: u8, max_steps: usize, n_threads: usize) -> Vec<Vec<P64>> {
    match ring {
        4 => polygons::<ZZ4>(enumerate_dispatch::<ZZ4>(max_steps, n_threads)),
        8 => polygons::<ZZ8>(enumerate_dispatch::<ZZ8>(max_steps, n_threads)),
        12 => polygons::<ZZ12>(enumerate_dispatch::<ZZ12>(max_steps, n_threads)),
        16 => polygons::<ZZ16>(enumerate_dispatch::<ZZ16>(max_steps, n_threads)),
        20 => polygons::<ZZ20>(enumerate_dispatch::<ZZ20>(max_steps, n_threads)),
        24 => polygons::<ZZ24>(enumerate_dispatch::<ZZ24>(max_steps, n_threads)),
        60 => polygons::<ZZ60>(enumerate_dispatch::<ZZ60>(max_steps, n_threads)),
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

fn run_rat_enum_seqs(ring: u8, max_steps: usize, n_threads: usize) -> Vec<Vec<i8>> {
    let f: fn(usize, usize) -> Vec<Vec<i8>> = match ring {
        4 => enumerate_dispatch::<ZZ4>,
        8 => enumerate_dispatch::<ZZ8>,
        12 => enumerate_dispatch::<ZZ12>,
        16 => enumerate_dispatch::<ZZ16>,
        20 => enumerate_dispatch::<ZZ20>,
        24 => enumerate_dispatch::<ZZ24>,
        60 => enumerate_dispatch::<ZZ60>,
        _ => panic!("invalid ring selected"),
    };
    f(max_steps, n_threads)
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
            #[cfg(feature = "pprof")]
            let _guard = cli.profile.as_ref().map(|_| {
                pprof::ProfilerGuardBuilder::default()
                    .frequency(1000)
                    .blocklist(&["libc", "libgcc", "pthread", "vdso"])
                    .build()
                    .expect("profiler start")
            });

            let t0 = Instant::now();
            let rats: Vec<Vec<i8>> = run_rat_enum_seqs(cli.ring, cli.max_steps, n_threads);
            let dt = t0.elapsed();

            // Use the result so it can't be trivially optimized away.
            let total_boundary_len: usize = rats.iter().map(|s| s.len()).sum();

            println!(
                "benchmark: ring={} max_steps={} -> {} unique rats (total boundary len={}) in {:?}",
                cli.ring,
                cli.max_steps,
                rats.len(),
                total_boundary_len,
                dt
            );

            #[cfg(feature = "pprof")]
            if let (Some(path), Some(guard)) = (cli.profile.as_ref(), _guard.as_ref()) {
                let report = guard.report().build().expect("profiler report");
                let svg = std::fs::File::create(path).expect("create flamegraph");
                report.flamegraph(svg).expect("write flamegraph");
                println!("wrote flamegraph: {path}");
            }
            #[cfg(not(feature = "pprof"))]
            if cli.profile.is_some() {
                eprintln!("rat_enum: --profile requires the `pprof` cargo feature");
            }

            if cli.stats {
                print_stats(&rats);
            }
        }
        Mode::Render => {
            let rats: Vec<Vec<P64>> = run_rat_enum_polylines(cli.ring, cli.max_steps, n_threads);

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
            .with_vertex_marker(MarkerStyle::filled_circle(
                0.06 * r,
                Color::RED,
            ))
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
