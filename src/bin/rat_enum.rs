//! `rat_enum`: enumerate every simple polygon on a cyclotomic-ring
//! lattice with boundary length up to `n`. Output mode (Render/Bench/
//! ListSeeds) is controlled by `--mode`.
//!
//! # Memory profile and how to scale to large `n`
//!
//! The result set is a `HashSet<Vec<i8>>` of canonical (or
//! dihedral-canonical) angle sequences. The *raw* sequence data is
//! tiny (`n` bytes per polygon), but several effects multiply the
//! peak memory of a single process by ~50-100x:
//!
//! 1. **HashSet bucket overhead.** Hashbrown stores each entry with
//!    a control byte, a cached hash, and the `Vec<i8>` header (24
//!    bytes pointer/len/cap). The `Vec`'s heap allocation is
//!    malloc-rounded (~16-32 bytes for an n=13 payload). With load
//!    factor ~0.875, one fully-populated set is ~80 bytes per
//!    polygon -- already 6x the raw payload.
//!
//! 2. **Hashbrown grows by doubling.** Inside a resize from 2^k to
//!    2^(k+1) buckets, *both* arrays are alive while entries are
//!    copied. A set that just resized briefly holds ~3x its final
//!    structure size. When 16 worker threads are all near the same
//!    fill level (atomic-counter dispatch keeps work balanced),
//!    several can resize concurrently.
//!
//! 3. **Per-thread duplication during parallel DFS.** In
//!    `rat_enum_parallel`, each of `n_threads` workers maintains its
//!    own thread-local `HashSet`. These sets are live for the entire
//!    DFS phase and only consolidated at the end. Peak during DFS:
//!    `n_threads * per-thread-final-size * resize-factor`.
//!
//! 4. **Merge double-buffer.** Once workers complete, the main
//!    thread folds each local into a `merged` set via `extend`.
//!    The local being drained is freed only at end of the
//!    for-iteration; the others sit in their `JoinHandle`s until
//!    consumed. Right before the last drain: 1 unconsumed local,
//!    plus `merged` at ~final size, plus a transient resize buffer
//!    inside `merged` -- another ~2x at peak.
//!
//! 5. **glibc malloc arena retention.** Each worker thread's malloc
//!    arena keeps its high-water-mark pages mapped, even after the
//!    transient `Vec<i8>` allocations from inside the close-event
//!    closure are freed. This isn't huge per thread, but it does
//!    add up at scale.
//!
//! Empirically at ZZ12 dihedral n=13 (4.08M polygons, ~50 MB raw),
//! single-process 16-thread peak is **~5.5 GB** -- a ~100x blowup.
//!
//! # Recommendation: separate processes
//!
//! **Running N independent single-threaded processes is strictly
//! better than running one process with `--threads N`** for memory:
//!
//! * Each process holds only its share of the polygons (~1x final
//!   set size for that share, plus the resize transient: peak ~3x).
//! * No per-thread duplication: there's only one HashSet per
//!   process.
//! * No cross-thread merge double-buffer: each process's final set
//!   *is* the final set for that share.
//! * On process exit, the OS reclaims **every page** (arenas
//!   included). Sequential or wall-clock-parallel batches never
//!   share memory across batches.
//!
//! Mechanics: use `--mode list-seeds` to emit a list of work-unit
//! prefixes, then run one `--seed <prefix> --threads 1` process per
//! prefix. An external orchestrator (`xargs -P N`, GNU `parallel`,
//! Slurm, etc.) picks how many run concurrently -- and that becomes
//! the only knob that affects total wall time. Per-process peak
//! memory is *unchanged* by orchestrator parallelism, because each
//! process holds only its own share.
//!
//! Example for n=14 at ZZ12 with 16 cores:
//!
//! ```sh
//! ./target/release/rat_enum --ring 12 -n 14 --dihedral \
//!     --mode list-seeds > listing.txt
//! grep '^SEED ' listing.txt | sed 's/SEED //' | \
//!   xargs -I{} -P 16 \
//!     ./target/release/rat_enum --ring 12 -n 14 --dihedral \
//!     --seed "{}" --threads 1 \
//!   > all_rats.txt
//! { grep '^RAT ' listing.txt; cat all_rats.txt; } \
//!   | sort -u | grep -c '^RAT '
//! ```
//!
//! Peak memory: ~16 processes * ~50-100 MB each ≈ ~1 GB total,
//! versus ~30-40 GB for the equivalent single-process run.
//!
//! `--seed <prefix> --threads N>1` exists (it routes through
//! `enumerate_from_seed`'s sub-split + parallel-drain path) for
//! environments where an orchestrator gives one big process at a
//! time, but it has the same per-process memory profile as
//! `--mode bench --threads N` -- it does NOT recover the
//! single-threaded-per-process savings.

use std::fs::File;
use std::io::BufWriter;
use std::sync::Mutex;
use std::time::Instant;

use clap::{Parser, ValueEnum};

use tilezz::rat_enum::dfs::STREAM_RAT_LINES;
use tilezz::rat_enum::output::{print_stats, run_rat_enum_polylines};
use tilezz::rat_enum::prune::{install_closure_key_prune, install_mod_prune};
use tilezz::rat_enum::run_rat_enum_seqs;
use tilezz::rat_enum::seed::{dispatch_collect_seed_prefixes, dispatch_enumerate_from_seed};
use tilezz::stringmatch::RatDafsa;

// Test-only imports. The two private test helpers below (`rat_enum`,
// `rat_enum_dihedral`) and the `dihedral_tests` / `opt_correctness_tests`
// modules pull these via `super::*`; in production code none of these
// is referenced from the binary.
#[cfg(test)]
use std::collections::HashSet;
#[cfg(test)]
use std::sync::Arc;
#[cfg(test)]
use tilezz::cyclotomic::{IsRing, ZZ10, ZZ12, ZZ16, ZZ20, ZZ24, ZZ32, ZZ4, ZZ60, ZZ8};
#[cfg(test)]
use tilezz::rat_enum::canonical::{dihedral_canonical, make_ops};
#[cfg(test)]
use tilezz::rat_enum::dfs::rat_enum_with;
#[cfg(test)]
use tilezz::rat_enum::prune::{
    closure_key::{collect_closure_keys, ClosureKeyPrune},
    modular::ModularPrune,
    units::unit_vectors_for_ring,
    Prunes,
};
#[cfg(test)]
use tilezz::rat_enum::seed::rat_enum_parallel;
#[cfg(test)]
use tilezz::rat_enum::stats::DfsStats;
use tilezz::vis::animation::render_gif;
use tilezz::vis::draw::{MarkerStyle, TileStyle};
use tilezz::vis::plotutils::P64;
use tilezz::vis::scene::{Color, Fill, Scene, Stroke, TextStyle, Viewport};

static VERBOSE: Mutex<bool> = Mutex::new(false);


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
        make_ops(false),
        "enumeration",
        "",
        false,
        &Prunes::default(),
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
        make_ops(true),
        "dihedral enumeration",
        "dihedral ",
        false,
        &Prunes::default(),
    )
}


// --------

#[derive(Copy, Clone, Debug, ValueEnum)]
enum Mode {
    /// Enumerate and render output (GIF)
    Render,
    /// Enumerate only and report elapsed time
    Bench,
    /// Walk the DFS down to `split_depth` (default 3) and print alive
    /// prefixes one per line as "SEED a,b,c", plus any polygon that
    /// closed before reaching that depth as "RAT [...]". Use with
    /// `--seed` to dispatch per-prefix jobs externally and merge.
    ListSeeds,
    /// Enumerate, then write a gzipped `RatDafsa` to the path given
    /// by `-o`. The on-disk format is `tilezz-rat-dafsa`; assigned
    /// external indices follow `(length asc, lex asc)` order. The
    /// length-prefix encoding used internally is described in
    /// `RatDafsa::JSON_SCHEMA_DOC`.
    Dafsa,
}

#[derive(Parser, Debug)]
#[command(
    version,
    about = "Enumerate simple cyclotomic matchstick polygons (rats) up to a given perimeter",
    long_about = "\
Enumerate every simple polygon (closed self-avoiding boundary) of\n\
unit-length edges on a cyclotomic ring `ZZn`, with perimeter up to `-n`.\n\
\n\
Output is selected via `--mode`:\n\
  * bench         enumerate, time, report counts (default workflow)\n\
  * render        enumerate and render the polygons as a GIF (`-o file.gif`)\n\
  * dafsa         enumerate and write a gzipped DAFSA (`-o file.bin.gz`)\n\
  * dafsa-blocks  same, but as a directory of lazy-loadable blocks\n\
  * list-seeds    walk DFS to `--split-depth` and print SEED + RAT lines\n\
                  (for external orchestration; see `--seed`)\n\
\n\
Performance opts (combine for maximum speedup):\n\
  --mod-prune          modular reachability prune; up to 376x on some rings\n\
  --closure-key-prune  exact lattice closure-key prune; another 2-5x on top\n\
  --threads N          parallel DFS; sub-linear scaling due to set merging\n\
  --dihedral           dihedral-canonical output (one rep per chiral pair);\n\
                       also accelerates DFS via complement-reflection prune\n\
\n\
For memory at large n: prefer `--mode list-seeds` + per-seed processes\n\
(via xargs / GNU parallel) over `--threads N>1`. See the file-level\n\
docstring in src/bin/rat_enum.rs for the full memory accounting.\n"
)]
struct Cli {
    /// Cyclotomic ring index n in `ZZn`. Supported:
    /// 4, 8, 10, 12, 16, 20, 24, 32, 60.
    #[arg(short = 'r', long)]
    ring: u8,

    /// Maximum boundary perimeter to enumerate up to (inclusive).
    /// Memory and time scale exponentially in `n`; expect each
    /// `+1` to multiply runtime by `~3-10` depending on ring.
    #[arg(short = 'n', long)]
    max_steps: usize,

    /// Output mode. See the top-level `--help` for the menu.
    #[arg(long, value_enum, default_value_t = Mode::Render)]
    mode: Mode,

    /// Output filename. Meaning depends on `--mode`:
    ///   * render        path/to/file.gif
    ///   * dafsa         path/to/file.bin.gz
    ///   * dafsa-blocks  path/to/output_dir/ (will be created)
    #[arg(short = 'o', long)]
    filename: Option<String>,

    /// Print extra progress / diagnostic chatter.
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

    /// After every successful `Snake::add`, replay the entire current
    /// angle prefix in a fresh `Snake::new()` and assert that fresh
    /// snake accepts the same prefix and reaches the same
    /// `is_closed` state. Panics on any disagreement -- traps cases
    /// where the stateful incremental check (post-pop grid state)
    /// admits prefixes that a from-scratch check would reject.
    /// Roughly O(n) extra work per step; expect ~2x slowdown.
    #[arg(long)]
    paranoid: bool,

    /// Enable the modular reachability prune. For each modulus
    /// `m in {2..=16}` whose state space `m^phi` fits the cell
    /// budget, precomputes the cumulative reachable mod-m
    /// displacement set and rejects any candidate direction whose
    /// post-extension displacement can't reach 0 mod m in any
    /// number of remaining steps. Cheap per-node (one packed hash
    /// lookup per modulus). See `--mode probe-modular` for the
    /// per-ring saturation curves and `--mod-prune-moduli` for
    /// A/B testing different modulus sets.
    #[arg(long)]
    mod_prune: bool,

    /// Comma-separated list of moduli to use with `--mod-prune`
    /// (default: 2..=16, filtered by ring's state-space budget).
    /// Useful for A/B testing: e.g. `--mod-prune-moduli 2,3,4,6`
    /// reproduces the original hardcoded set.
    #[arg(long, value_delimiter = ',', num_args = 1..)]
    mod_prune_moduli: Option<Vec<i64>>,

    /// Enable the closure-key prune: pre-enumerate every simple
    /// open snake up to length L (= `--closure-key-depth`) and
    /// store their (endpoint, facing) keys. During DFS, when the
    /// remaining-after-this-direction count is <= L, the candidate
    /// is pruned iff its required suffix's closure key isn't in
    /// the precomputed set. Strictly stronger than any modular
    /// projection (uses exact lattice + facing info). More memory
    /// and pre-pass cost; per-node cost is one ring mul + one
    /// hash lookup.
    #[arg(long)]
    closure_key_prune: bool,

    /// Suffix-length cap for `--closure-key-prune`. Memory/pre-pass
    /// cost scales as |S_L| ~ c^L (c ~ 5-10 for typical rings);
    /// the prune only fires when `remaining_after <= L`. L=4 was
    /// empirically optimal across ZZ8/12/16 at n=10..15: larger L
    /// finds more skips but the table-build cost dominates the
    /// per-node savings.
    #[arg(long, default_value_t = 4)]
    closure_key_depth: usize,

    /// Resume the DFS from a specific seed prefix (comma-separated
    /// angles like `-5,3,-4`). When set, the DFS skips seed
    /// collection and resumes from this prefix; output is the same
    /// RAT lines as Bench mode. Honours `--threads`. Used together
    /// with `--mode list-seeds` to dispatch per-prefix jobs
    /// externally.
    #[arg(long, value_delimiter = ',', allow_hyphen_values = true, num_args = 1)]
    seed: Option<Vec<i8>>,

    /// DFS splitting depth for `--mode list-seeds` and for the
    /// parallel (`--threads > 1`) seed-collection step. Defaults to
    /// 3 in list-seeds mode (~50-300 seeds for typical n), or to
    /// the auto-picked value otherwise.
    #[arg(long)]
    split_depth: Option<usize>,
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

    if cli.mod_prune {
        install_mod_prune(
            cli.ring,
            cli.max_steps,
            cli.mod_prune_moduli.as_deref(),
        );
    }
    if cli.closure_key_prune {
        install_closure_key_prune(cli.ring, cli.closure_key_depth);
    }
    // `--mode list-seeds`: walk the DFS down to split_depth, print
    // alive prefixes as `SEED a,b,c` lines and any already-closed
    // polygons as `RAT [...]` lines, then exit. The split_depth
    // defaults to 3 (typically ~50-300 seeds at usable n).
    if matches!(cli.mode, Mode::ListSeeds) {
        let split_depth = cli.split_depth.unwrap_or(3);
        let (closed, seeds) = dispatch_collect_seed_prefixes(
            cli.ring,
            cli.max_steps,
            cli.step,
            split_depth,
            cli.dihedral,
        );
        for prefix in &seeds {
            let s = prefix
                .iter()
                .map(|a| a.to_string())
                .collect::<Vec<_>>()
                .join(",");
            println!("SEED {s}");
        }
        // Polygons that closed before reaching split_depth (typically
        // perimeter-3 only, e.g. the triangle at ZZ12). These would
        // be missed if the orchestrator only ran the SEED jobs.
        for seq in &closed {
            println!("RAT {seq:?}");
        }
        eprintln!(
            "list-seeds: {} alive prefixes at depth {}, {} polygons closed before split",
            seeds.len(),
            split_depth,
            closed.len(),
        );
        return;
    }

    // `--seed <prefix>`: skip seed collection, resume DFS from the
    // given prefix, print RAT lines. With `--threads 1` (default and
    // recommended), this is one work unit; the orchestrator runs
    // multiple `--seed` jobs as separate processes to use more cores.
    // With `--threads N>1`, the seed's subtree is itself sub-split
    // for in-process parallelism -- useful when an orchestrator
    // gives one big process at a time, but per-process memory is
    // ~N x higher. See the file-level docstring for why separate
    // processes are strictly better for memory.
    if let Some(seed) = cli.seed.as_deref() {
        let t0 = Instant::now();
        let rats = dispatch_enumerate_from_seed(
            cli.ring,
            cli.max_steps,
            cli.step,
            seed,
            n_threads,
            cli.dihedral,
            cli.paranoid,
        );
        let dt = t0.elapsed();
        for seq in &rats {
            println!("RAT {seq:?}");
        }
        eprintln!("seed {:?}: {} unique rats in {dt:?}", seed, rats.len());
        return;
    }

    match cli.mode {
        Mode::ListSeeds => unreachable!("handled above"),
        Mode::Dafsa => {
            let Some(filename) = cli.filename.as_deref() else {
                eprintln!("--mode dafsa requires -o <output path>");
                std::process::exit(2);
            };

            // The DFS streams each newly-found rat to stdout by
            // default (for `--seed` / `--mode bench`); in dafsa mode
            // the artifact is the binary file at `-o`, so silence the
            // streaming lines before enumerating.
            STREAM_RAT_LINES.store(false, std::sync::atomic::Ordering::Relaxed);

            let t0 = Instant::now();
            let (rats, _stats) = run_rat_enum_seqs(
                cli.ring,
                cli.max_steps,
                cli.step,
                n_threads,
                cli.dihedral,
                cli.paranoid,
            );
            eprintln!("enumerated {} rats in {:?}", rats.len(), t0.elapsed());

            // `RatDafsa::from_rats` handles the (length, lex) sort,
            // dedup, and length-prefix encoding internally; the
            // resulting external index of each rat is its position in
            // the (length, lex)-ordered sequence.
            let t1 = Instant::now();
            let dafsa = RatDafsa::from_rats(rats.iter().map(|r| r.as_slice()));
            eprintln!(
                "built RatDafsa ({} entries) in {:?}",
                dafsa.len(),
                t1.elapsed()
            );

            let t2 = Instant::now();
            let file = File::create(filename).expect("create output file");
            dafsa
                .write_json_gz(BufWriter::new(file))
                .expect("write gzipped RatDafsa");
            let bytes = std::fs::metadata(filename).map(|m| m.len()).unwrap_or(0);
            eprintln!("wrote {filename} ({bytes} bytes) in {:?}", t2.elapsed());

            if cli.stats {
                print_stats(&rats);
            }
        }
        Mode::Bench => {
            let profile = tilezz::util::profile::ProfileGuard::start(cli.profile.as_deref());

            let t0 = Instant::now();
            let (rats, stats) = run_rat_enum_seqs(
                cli.ring,
                cli.max_steps,
                cli.step,
                n_threads,
                cli.dihedral,
                cli.paranoid,
            );
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
            let rats: Vec<Vec<P64>> = run_rat_enum_polylines(
                cli.ring,
                cli.max_steps,
                cli.step,
                n_threads,
                cli.dihedral,
                cli.paranoid,
            );

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
    use std::collections::HashMap;
    use tilezz::cyclotomic::geometry::intersect;
    use tilezz::cyclotomic::{IsRing, Units, ZZ12, ZZ4, ZZ8};

    /// Independent validator: reconstruct the polygon from `angles` in
    /// exact ring arithmetic and check three properties without any
    /// dependence on `Snake::can_add`:
    ///   (a) the walk closes at the origin;
    ///   (b) the n vertices visited (origin + n-1 intermediates) are
    ///       all pairwise distinct;
    ///   (c) no two non-adjacent boundary edges share any point
    ///       (proper crossing or shared interior).
    fn validate_simple_polygon<ZZ: IsRing>(angles: &[i8]) -> Result<(), String> {
        let n = angles.len();
        if n < 3 {
            return Err(format!("perimeter {n} < 3"));
        }
        // Reconstruct vertices.
        let mut pts: Vec<ZZ> = Vec::with_capacity(n + 1);
        pts.push(ZZ::zero());
        let mut dir: i64 = 0;
        for &a in angles {
            dir = (dir + a as i64).rem_euclid(ZZ::turn() as i64);
            let next = *pts.last().unwrap() + <ZZ as Units>::unit(dir as i8);
            pts.push(next);
        }
        // (a) closes.
        if !pts.last().unwrap().is_zero() {
            return Err(format!("does not close: last={:?}", pts.last().unwrap()));
        }
        pts.pop(); // drop the closing duplicate; now pts has the n cycle vertices

        // (b) distinct vertices. Use a HashMap to find duplicates and
        // report both indices.
        let mut seen: HashMap<ZZ, usize> = HashMap::new();
        for (i, &p) in pts.iter().enumerate() {
            if let Some(&j) = seen.get(&p) {
                return Err(format!("duplicate vertex: pts[{j}] == pts[{i}]"));
            }
            seen.insert(p, i);
        }

        // (c) no two non-adjacent edges share any point. Brute O(n^2)
        // check via `intersect`. Adjacent edges share an endpoint by
        // construction, so we skip cyclically-adjacent pairs.
        for i in 0..n {
            let a1 = pts[i];
            let a2 = pts[(i + 1) % n];
            for j in (i + 2)..n {
                // Skip the cyclic-adjacent pair (last edge vs first edge).
                if i == 0 && j == n - 1 {
                    continue;
                }
                let b1 = pts[j];
                let b2 = pts[(j + 1) % n];
                if intersect(&(a1, a2), &(b1, b2)) {
                    return Err(format!(
                        "edges {i} ({a1:?}->{a2:?}) and {j} ({b1:?}->{b2:?}) intersect/touch"
                    ));
                }
            }
        }

        Ok(())
    }

    /// Per-boundary-length counts the ZZ12 dihedral enumeration must
    /// match: OEIS A316192 "Number of self-avoiding polygons with
    /// perimeter n and sides = 1 ... not counting rotations and
    /// reflections as distinct" (the boundary path may nowhere touch
    /// or intersect itself).
    ///
    /// Reference values, indexed by perimeter `n`:
    ///   n=1..10 -> [0, 0, 1, 3, 4, 22, 69, 418, 2210, 14024]
    ///
    /// Historical note: this test used to fail at n=9 (2217 vs 2210)
    /// and n=10 (14124 vs 14024). The extras were polygons whose
    /// boundary had a T-touch -- a vertex of one edge sitting on the
    /// strict interior of another edge -- which the segment
    /// intersection predicate missed because the orientation tests
    /// collapsed `wedge == 0` ("on the line") into the same branch
    /// as `wedge < 0` ("right side"). Fixed in `b2df8f4` by adding
    /// explicit endpoint-on-other-segment checks for the colinear
    /// case. T-touches first occur at n=9 on ZZ12, so this test now
    /// covers the bug-sensitive range exactly. See the focused unit
    /// test `test_intersect_zz12_t_touch_endpoint_on_segment` in
    /// `cyclotomic::geometry::tests`.
    #[test]
    fn test_oeis_a316192_zz12() {
        const OEIS: &[(usize, usize)] = &[
            (3, 1),
            (4, 3),
            (5, 4),
            (6, 22),
            (7, 69),
            (8, 418),
            (9, 2210),
            (10, 14024),
        ];
        // Enumerate up to the largest n in the reference. The DFS
        // walks the dihedral-canonical variant; output is bucketed
        // per perimeter length.
        let max_n = OEIS.iter().map(|&(n, _)| n).max().unwrap();
        let (rats, _) = super::rat_enum_dihedral::<ZZ12>(max_n, 1);

        // Bucket per perimeter length.
        let mut by_len: std::collections::BTreeMap<usize, usize> =
            std::collections::BTreeMap::new();
        for seq in &rats {
            *by_len.entry(seq.len()).or_insert(0) += 1;
        }

        let mut mismatches: Vec<(usize, usize, usize)> = Vec::new();
        for &(n, expected) in OEIS {
            let got = by_len.get(&n).copied().unwrap_or(0);
            if got != expected {
                mismatches.push((n, got, expected));
            }
        }
        if !mismatches.is_empty() {
            for (n, got, expected) in &mismatches {
                eprintln!(
                    "n={n}: got {got}, expected (OEIS A316192) {expected}, diff {:+}",
                    *got as i64 - *expected as i64
                );
            }
            panic!(
                "dihedral ZZ12 enumeration differs from OEIS A316192 at {} perimeter length(s)",
                mismatches.len()
            );
        }
    }

    #[test]
    fn test_dihedral_output_polygons_are_simple_and_unique() {
        let (dihedral_rats, _) = super::rat_enum_dihedral::<ZZ12>(9, 1);

        // (i) Encodings are pairwise distinct (the HashSet via which
        // they were collected enforces this; assert defensively).
        let unique: std::collections::HashSet<Vec<i8>> = dihedral_rats.iter().cloned().collect();
        assert_eq!(
            unique.len(),
            dihedral_rats.len(),
            "duplicate sequences in dihedral output"
        );

        // (ii) Every polygon is geometrically simple per the
        // independent validator above. Collect failures so a single
        // run reports every case.
        let mut failures: Vec<(Vec<i8>, String)> = Vec::new();
        for seq in &dihedral_rats {
            if let Err(why) = validate_simple_polygon::<ZZ12>(seq) {
                failures.push((seq.clone(), why));
            }
        }
        if !failures.is_empty() {
            eprintln!("{} polygons failed independent validation:", failures.len());
            for (seq, why) in failures.iter().take(20) {
                eprintln!("  {seq:?} -- {why}");
            }
            panic!("non-simple polygons in dihedral enumeration output");
        }
    }

    /// For a single `(ring, max_steps)` pair, check that the
    /// dihedral DFS output equals the rotation DFS output quotiented
    /// by `dihedral_canonical`. Stage-1 over-pruning shows up as
    /// MISSING dihedral classes; stage-2 under-deduplication shows
    /// up as EXTRA. The body prints both lists before asserting so
    /// any future regression is diagnosable from the test log.
    fn check_quotient_match<ZZ: IsRing>(ring_label: &str, max_steps: usize) {
        let (all_rats, _) = super::rat_enum::<ZZ>(max_steps, 1);
        let mut expected: HashSet<Vec<i8>> = HashSet::new();
        for seq in &all_rats {
            expected.insert(super::dihedral_canonical(seq));
        }
        let (dihedral_rats, _) = super::rat_enum_dihedral::<ZZ>(max_steps, 1);
        let actual: HashSet<Vec<i8>> = dihedral_rats.into_iter().collect();

        if expected != actual {
            for s in expected.difference(&actual) {
                eprintln!("[{ring_label} n={max_steps}] MISSING: {s:?}");
            }
            for s in actual.difference(&expected) {
                eprintln!("[{ring_label} n={max_steps}] EXTRA:   {s:?}");
            }
        }
        assert_eq!(
            expected, actual,
            "{ring_label} n={max_steps}: dihedral DFS != rotation DFS / dihedral",
        );
        eprintln!(
            "[{ring_label} n={max_steps}] OK -- {} dihedral classes from {} rotation-canonical rats",
            actual.len(),
            all_rats.len(),
        );
    }

    /// Sweep multiple `(ring, max_steps)` pairs. ZZ4 and ZZ8 are
    /// included to catch ring-specific over-pruning: with only 3 or
    /// 7 candidate directions per step, the per-rotation comparisons
    /// are more constrained and an off-by-one in the prefix prune is
    /// easier to spot than at ZZ12.
    #[test]
    fn test_dihedral_enum_matches_dfs_quotient() {
        for n in [4, 5, 6, 7, 8, 9] {
            check_quotient_match::<ZZ12>("ZZ12", n);
        }
        for n in [4, 6, 8, 10, 12] {
            check_quotient_match::<ZZ8>("ZZ8", n);
        }
        for n in [4, 6, 8, 10, 12, 14] {
            check_quotient_match::<ZZ4>("ZZ4", n);
        }
    }

    /// Mechanics of the `--mode list-seeds` + per-seed dispatch:
    /// running the DFS in one shot must yield the same set of
    /// canonical sequences as collecting seeds at some `split_depth`,
    /// then enumerating from each seed and union-ing with the
    /// already-closed polygons. Tested across:
    ///   * small `max_steps` (5-7) for tractable runtime;
    ///   * dihedral on AND off (each uses a different `CanonicalOps`
    ///     pair and a different dedup hash);
    ///   * `split_depth` in 1-3 (1: every root-direction branch is a
    ///     seed; 3: small polygons close before split, exercising
    ///     the `closed` half of the merge).
    ///
    /// This is the regression target for the seed-partition API.
    /// Catches any future drift where collect_seeds and the per-seed
    /// runner stop agreeing on what "the same enumeration" means.
    #[test]
    fn test_seed_partitioning_matches_one_shot() {
        for &n in &[5usize, 6, 7] {
            for &dihedral in &[false, true] {
                for &split_depth in &[1usize, 2, 3] {
                    if split_depth >= n {
                        continue;
                    }
                    let (one_shot_seqs, _) = if dihedral {
                        super::rat_enum_dihedral::<ZZ12>(n, 1)
                    } else {
                        super::rat_enum::<ZZ12>(n, 1)
                    };
                    let one_shot: HashSet<Vec<i8>> = one_shot_seqs.into_iter().collect();

                    let (mut seeded, prefixes) =
                        tilezz::rat_enum::seed::collect_seed_prefixes::<ZZ12>(n, 1, split_depth, dihedral);
                    for prefix in &prefixes {
                        // Sweep both branches of the n_threads dispatch
                        // (single-threaded direct path AND sub-split +
                        // parallel_drain_seeds) within the same test --
                        // they must produce the same set.
                        for nthreads in &[1usize, 4] {
                            seeded.extend(tilezz::rat_enum::seed::enumerate_from_seed::<ZZ12>(
                                n, 1, prefix, *nthreads, dihedral, false,
                            ));
                        }
                    }

                    assert_eq!(
                        one_shot,
                        seeded,
                        "n={n} dihedral={dihedral} split_depth={split_depth}: \
                         seed-partitioned != one-shot \
                         ({} prefixes + {} pre-closed)",
                        prefixes.len(),
                        seeded.len() - prefixes.iter().map(|_| 0).sum::<usize>(),
                    );
                }
            }
        }
    }
}

/// Correctness tests for the optional DFS prunes.
///
/// For each ring x mode (rotation/dihedral) x thread-count, the test
/// matrix runs the enumeration under every subset of the three
/// optional prunes (`mod`, `coord-proj`, `closure-key`) and asserts
/// the resulting canonical-rat set is bit-identical to the baseline
/// (no prunes). A regression in any single prune -- or any
/// composition -- shows up as a set-equality mismatch on these tests.
///
/// Caps: ring n is chosen so the full sweep stays under ~30 s on a
/// release build. ZZ12 stops at n=8 (no T-touch bug); ZZ8 / ZZ4 go
/// further within their compute budgets.
#[cfg(test)]
mod opt_correctness_tests {
    use super::*;
    use std::sync::Arc;

    /// Dispatch helper: collect closure keys for `ring` up to length
    /// `max_l`. Mirrors the per-ring match in `install_closure_key_prune`.
    fn closure_keys_for(
        ring: u8,
        max_l: usize,
    ) -> rustc_hash::FxHashSet<(Vec<i64>, i8)> {
        match ring {
            4 => collect_closure_keys::<ZZ4>(max_l),
            8 => collect_closure_keys::<ZZ8>(max_l),
            10 => collect_closure_keys::<ZZ10>(max_l),
            12 => collect_closure_keys::<ZZ12>(max_l),
            16 => collect_closure_keys::<ZZ16>(max_l),
            20 => collect_closure_keys::<ZZ20>(max_l),
            24 => collect_closure_keys::<ZZ24>(max_l),
            32 => collect_closure_keys::<ZZ32>(max_l),
            60 => collect_closure_keys::<ZZ60>(max_l),
            _ => panic!("unknown ring {ring}"),
        }
    }

    /// Build a `Prunes` for testing per the chosen opt subset.
    /// Both prunes are constructed against `(ring, max_steps)`
    /// regardless of whether they're enabled, so each test case is
    /// self-contained and independent of any global state.
    fn build_prunes(
        ring: u8,
        max_steps: usize,
        with_mod: bool,
        with_ck: bool,
    ) -> Prunes {
        let (units, phi) = unit_vectors_for_ring(ring);
        let mut prunes = Prunes::default();
        if with_mod {
            let mp = ModularPrune::build(&units, phi, max_steps, None);
            prunes.mod_prune = Some(Arc::new(mp));
        }
        if with_ck {
            let max_l = 4;
            let keys = closure_keys_for(ring, max_l);
            prunes.closure_key_prune = Some(Arc::new(ClosureKeyPrune { max_l, keys }));
        }
        prunes
    }

    /// Run the DFS with explicit prunes on ring `ZZ` and return the
    /// resulting canonical-rat set. Bypasses the global PRUNES static
    /// so different test cases don't fight over it.
    fn run<ZZ: IsRing>(
        max_steps: usize,
        dihedral: bool,
        n_threads: usize,
        prunes: &Prunes,
    ) -> std::collections::HashSet<Vec<i8>> {
        let ops = make_ops(dihedral);
        let (rats, _) = if n_threads <= 1 {
            rat_enum_with::<ZZ>(max_steps, 1, ops, "test", "", false, prunes)
        } else {
            rat_enum_parallel::<ZZ>(max_steps, 1, n_threads, ops, "test", "", false, prunes)
        };
        rats.into_iter().collect()
    }

    /// Enumerate the 4 subsets of {mod, closure-key}.
    fn opt_subsets() -> impl Iterator<Item = (bool, bool)> {
        (0..4).map(|mask| (mask & 1 != 0, mask & 2 != 0))
    }

    /// Cross-validation matrix: every opt subset x every thread count
    /// must produce the same canonical-rat set as the baseline.
    fn check_ring<ZZ: IsRing>(ring: u8, max_steps: usize) {
        for &dihedral in &[false, true] {
            let baseline = run::<ZZ>(max_steps, dihedral, 1, &Prunes::default());
            for (with_mod, with_ck) in opt_subsets() {
                let prunes = build_prunes(ring, max_steps, with_mod, with_ck);
                for &n_threads in &[1usize, 4] {
                    let got = run::<ZZ>(max_steps, dihedral, n_threads, &prunes);
                    assert_eq!(
                        got,
                        baseline,
                        "ZZ{ring} n={max_steps} dihedral={dihedral} \
                         mod={with_mod} ck={with_ck} threads={n_threads}: \
                         result set differs from baseline ({} vs {} rats)",
                        got.len(),
                        baseline.len(),
                    );
                }
            }
        }
    }

    #[test]
    fn cross_validate_zz4() {
        check_ring::<ZZ4>(4, 10);
    }

    #[test]
    fn cross_validate_zz8() {
        check_ring::<ZZ8>(8, 10);
    }

    #[test]
    fn cross_validate_zz12() {
        // Capped at n=8 for runtime budget -- the 32-combo matrix
        // at higher n explodes runtime quickly. The OEIS test below
        // covers ZZ12 dihedral up to n=10 for every opt combo,
        // which is the bug-sensitive range (T-touches first occur
        // at n=9 -- see the historical note in test_oeis_a316192_zz12).
        check_ring::<ZZ12>(12, 8);
    }

    /// External anchor: OEIS A316192 per-length counts for ZZ12
    /// dihedral matchstick polygons, n=3..=10. The range
    /// deliberately includes n=9, 10 -- the T-touch-sensitive
    /// lengths that the now-fixed segment-intersection bug used to
    /// miscount. Every opt combination must reproduce the OEIS
    /// reference there: a regression in any prune that re-introduces
    /// the missed T-touch (or any other geometry weakness) shows up
    /// as a count mismatch.
    ///
    /// Slow test (~60-90 s release): n=10 baseline alone is ~22s,
    /// repeated across opt subsets. Worth it -- this is the one
    /// test that ties our enumeration to a published external
    /// reference at the bug-sensitive lengths.
    #[test]
    fn oeis_a316192_each_opt_combo() {
        const OEIS: &[(usize, usize)] = &[
            (3, 1),
            (4, 3),
            (5, 4),
            (6, 22),
            (7, 69),
            (8, 418),
            (9, 2210),  // T-touches first appear here; bug-fix-sensitive
            (10, 14024), // bug-fix-sensitive
        ];
        let max_n = OEIS.iter().map(|&(n, _)| n).max().unwrap();

        for (with_mod, with_ck) in opt_subsets() {
            let prunes = build_prunes(12, max_n, with_mod, with_ck);
            for &n_threads in &[1usize, 4] {
                let rats = run::<ZZ12>(max_n, true, n_threads, &prunes);
                let mut by_len: std::collections::BTreeMap<usize, usize> =
                    std::collections::BTreeMap::new();
                for seq in &rats {
                    *by_len.entry(seq.len()).or_insert(0) += 1;
                }

                let mut mismatches: Vec<(usize, usize, usize)> = Vec::new();
                for &(n, expected) in OEIS {
                    let got = by_len.get(&n).copied().unwrap_or(0);
                    if got != expected {
                        mismatches.push((n, got, expected));
                    }
                }
                if !mismatches.is_empty() {
                    for (n, got, expected) in &mismatches {
                        eprintln!(
                            "n={n}: got {got}, expected (OEIS A316192) {expected}, diff {:+}",
                            *got as i64 - *expected as i64
                        );
                    }
                    panic!(
                        "OEIS mismatch with mod={with_mod} ck={with_ck} threads={n_threads}: \
                         {} length(s) differ",
                        mismatches.len()
                    );
                }
            }
        }
    }
}
