//! `rat_enum`: enumerate every simple polygon on a cyclotomic-ring
//! lattice with boundary length up to `n`. Output is controlled by
//! `--mode` -- see the `Mode` enum for the full menu.
//!
//! # Scaling to large `n`
//!
//! The in-memory modes (`--mode bench`, `--mode render`, `--mode
//! dafsa`, `--mode dafsa-blocks`) hold the entire result set as a
//! `HashSet<Vec<i8>>` of canonical (or free, i.e. dihedral-reduced) angle
//! sequences while the DFS runs. The *raw* payload is tiny (~`n`
//! bytes per polygon), but several effects multiply peak RSS by
//! ~50-100x:
//!
//! 1. **HashSet bucket overhead.** Per-entry control byte + cached
//!    hash + `Vec<i8>` header (24 B) + the malloc-rounded heap slab
//!    for the angles. ~80 bytes per polygon at typical loads.
//! 2. **Hashbrown doubling resize.** Mid-resize, both old and new
//!    bucket arrays are live: ~3x transient blowup per set.
//! 3. **Per-thread duplication.** Each parallel worker keeps a
//!    private HashSet for the whole DFS; final merge fans them in.
//! 4. **Final-merge double-buffer.** Worker locals sit in their
//!    `JoinHandle`s while the main set grows by `extend`; peak
//!    adds another ~2x.
//! 5. **Glibc arena retention.** Per-thread arenas don't release
//!    pages back to the OS once the high-water mark is reached.
//!
//! Empirically: ZZ12 free n=13 (4.08M polygons, ~50 MB raw)
//! peaks at **~5.5 GB** in a 16-thread `--mode bench`. n=14 (~30M
//! polygons) is out of reach for any single in-memory run on a
//! commodity workstation.
//!
//! # Recommendation: the streaming pipeline (`--mode stream` then
//! `--mode merge`)
//!
//! Workers stream closures to per-thread, bounded sort-buffer run
//! files instead of accumulating in HashSets; a separate merge
//! pass k-way dedupes the runs into a single sorted `unique.bin`
//! plus a `certificate.json` with the BLAKE3 hash. Memory is
//! bounded by `n_threads * buffer_size` (~16 MB/worker default),
//! regardless of the final rat count. Single host, single
//! invocation, no orchestrator:
//!
//! ```sh
//! ./target/release/rat_enum --ring 12 -n 14 --free \
//!     --mode stream -o out/ --threads 16
//! ./target/release/rat_enum --ring 12 -n 14 --free \
//!     --mode merge  -o out/
//! ./target/release/rat_enum --ring 12 -n 14 --free \
//!     --mode build  -o out/ --target-block-bytes 1048576
//! # out/dafsa/ now holds the blocked RatDafsa, readable by
//! # LazyRatDafsa / LazyRatDafsaAsync. out/certificate.json
//! # carries the BLAKE3 of unique.bin so the build is auditable.
//! ```
//!
//! Stage 3 (`--mode build`) reads `unique.bin` and streams each
//! record straight through `RatDafsa::from_sorted_unique_rats`
//! into a DAFSA builder, so peak RSS scales with the compressed
//! automaton (~12-32 bytes per state) rather than with the input
//! set (~80 bytes per HashSet entry).
//!
//! # Alternative: multi-host distribution (`--mode list-seeds`)
//!
//! Still relevant when you want to fan one enumeration out across
//! multiple machines (each runs `--seed <prefix> --threads 1` on
//! its share of the seed list). Per-process memory stays bounded
//! by that process's share -- the same multipliers above apply,
//! just to a smaller subset. For a single host though, prefer the
//! streaming pipeline: tighter memory bound, no result-set
//! concatenation step needed.
//!
//! `--seed <prefix> --threads N>1` exists for environments where
//! the orchestrator hands you one big slot at a time, but it has
//! the same per-process memory profile as `--mode bench --threads
//! N` -- it does NOT recover the single-threaded-per-process
//! savings.

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

// Test-only imports needed by the file-top helpers (`rat_enum`,
// `rat_enum_free`) which the test modules call via `super::`.
// Imports used *only* inside a specific test module live inside that
// module's `use super::*;` block rather than at file top.
#[cfg(test)]
use std::collections::HashSet;
#[cfg(test)]
use tilezz::cyclotomic::IsRing;
#[cfg(test)]
use tilezz::rat_enum::canonical::{free_canonical, make_ops};
#[cfg(test)]
use tilezz::rat_enum::dfs::rat_enum_with;
#[cfg(test)]
use tilezz::rat_enum::prune::Prunes;
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

/// Free variant of [`rat_enum`].
///
/// Two-stage design (see module-level comment above for rationale):
///
/// 1. [`is_free_canonical_extended`] prunes walk prefixes whose
///    complement rotation is lex-smaller, reducing the search tree
///    (~2.3x speedup at ZZ12 n=10).
///
/// 2. [`free_canonical`] at closure maps the chirality-normalized
///    canonical rotation to the free (lex-min dihedral) form, so both
///    surviving members of each chiral pair hash to the same key.
///
/// Returns the same set as running [`rat_enum`] and quotienting by
/// free equivalence.
#[cfg(test)]
fn rat_enum_free<ZZ: IsRing>(max_steps: usize, step: i8) -> (Vec<Vec<i8>>, DfsStats) {
    rat_enum_with::<ZZ>(
        max_steps,
        step,
        make_ops(true),
        "free enumeration",
        "free ",
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
    /// Enumerate, then write a blocked (lazy-loadable) `RatDafsa`
    /// asset directory at the path given by `-o`. Produces
    /// `block_index.json` + `blocks/<sha256>.bin` gzipped files
    /// readable by [`tilezz::stringmatch::LazyRatDafsa`] (or an
    /// async equivalent). Pair with `--target-block-bytes` if the
    /// default is too coarse or too fine for your set size.
    DafsaBlocks,
    /// Stage 1 of the streaming pipeline. Runs the parallel DFS but
    /// writes each worker's closures to per-thread sort-buffer runs
    /// under `-o <dir>/runs/run_tNN_rMM.bin`, instead of accumulating
    /// in an in-memory HashSet. Memory is bounded by the per-thread
    /// buffer (~16 MB/worker default), so very large enumerations
    /// (e.g. ZZ12 n=15+) fit on a commodity workstation. Pair with
    /// `--mode merge` to fold the runs into a single sorted
    /// `unique.bin` + `certificate.json`, then `--mode build` to
    /// produce the blocked `RatDafsa` asset under `<-o dir>/dafsa/`
    /// without materialising the full set in memory.
    Stream,
    /// Stage 2 of the streaming pipeline. K-way merges `<-o
    /// dir>/runs/*.bin` (produced by `--mode stream`) into
    /// `<-o dir>/unique.bin` (deduped, sorted) and writes
    /// `<-o dir>/certificate.json` with the BLAKE3 hash + headline
    /// counts. Doesn't re-run the DFS -- just folds existing run
    /// files. Requires the same `--ring`, `-n`, `--step`,
    /// `--free` flags as the producing `--mode stream` call (they
    /// are recorded in the certificate verbatim).
    Merge,
    /// Stage 3 of the streaming pipeline. Reads `<-o dir>/unique.bin`
    /// (produced by `--mode merge`), streams the records through
    /// `RatDafsa::from_sorted_unique_rats`, and writes the blocked
    /// asset to `<-o dir>/dafsa/` via `RatDafsa::write_blocks`. Pair
    /// with `--target-block-bytes`. Does not run the DFS or
    /// re-merge; it's the streaming-friendly equivalent of `--mode
    /// dafsa-blocks` for inputs too large to fit in memory as a
    /// Vec<Vec<i8>>.
    Build,
}

#[derive(Parser, Debug)]
#[command(
    version,
    about = "Enumerate simple cyclotomic matchstick polygons (rats) up to a given perimeter",
    long_about = "\
Enumerate every simple polygon (closed self-avoiding boundary) of\n\
unit-length edges on a cyclotomic ring `ZZn`, with perimeter up to `-n`.\n\
\n\
Output is selected via `--mode` (default: render):\n\
  * render        enumerate and render the polygons as a GIF (`-o file.gif`)\n\
  * bench         enumerate, time, report counts (no -o needed)\n\
  * dafsa         enumerate and write a gzipped DAFSA (`-o file.bin.gz`)\n\
  * dafsa-blocks  same, but as a directory of lazy-loadable blocks\n\
  * stream        stage 1 of the streaming pipeline: write per-thread runs\n\
                  to `<-o dir>/runs/run_tNN_rMM.bin` (bounded memory)\n\
  * merge         stage 2 of the streaming pipeline: k-way merge runs into\n\
                  `<-o dir>/unique.bin` + `certificate.json`\n\
  * build         stage 3 of the streaming pipeline: stream `unique.bin`\n\
                  through a DAFSA builder and write blocks to `<-o dir>/dafsa/`\n\
  * list-seeds    walk DFS to `--split-depth` and print SEED + RAT lines\n\
                  (for multi-host orchestration; see `--seed`)\n\
\n\
Performance opts (combine for maximum speedup):\n\
  --mod-prune          modular reachability prune; up to 376x on some rings\n\
  --closure-key-prune  exact lattice closure-key prune; another 2-5x on top\n\
  --threads N          parallel DFS; near-linear in streaming modes,\n\
                       sub-linear in HashSet modes (set-merge overhead)\n\
  --free               free (= full dihedral reduction) output:\n\
                       one rep per chiral pair (lex-min over rotations and\n\
                       reflections); also accelerates DFS via\n\
                       complement-reflection prune\n\
\n\
For memory at large n: prefer the streaming pipeline (`--mode stream`\n\
then `--mode merge`) over `--mode bench --threads N`. See the file-level\n\
docstring in src/bin/rat_enum.rs for the full memory accounting.\n"
)]
struct Cli {
    /// Cyclotomic ring index n in `ZZn`. Supported:
    /// 4, 6, 8, 10, 12, 16, 20, 24, 32, 60.
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
    ///   * stream        path/to/output_dir/ (will be created; runs go in `<dir>/runs/`)
    ///   * merge         path/to/output_dir/ (same one used by `--mode stream`)
    ///   * build         path/to/output_dir/ (same; blocks go in `<dir>/dafsa/`)
    #[arg(short = 'o', long)]
    filename: Option<String>,

    /// Print extra progress / diagnostic chatter.
    #[arg(short, long)]
    verbose: bool,

    /// In Bench mode, write a flamegraph SVG to this path (requires
    /// the `debug` cargo feature).
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

    /// Free (= full dihedral symmetry reduction) enumeration: outputs
    /// one representative per chiral pair (lex-min over rotations
    /// AND reflections). Faster than the default rotation-canonical
    /// DFS because complement-reflection pruning halves the search
    /// space at the root.
    #[arg(long)]
    free: bool,

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
    /// lookup per modulus). Use `--mod-prune-moduli` to A/B test
    /// different modulus sets.
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

    /// Target uncompressed bytes per block file when writing the
    /// `--mode dafsa-blocks` asset. The writer closes a block once
    /// its serialised size crosses this threshold (the final block
    /// may be smaller). 1 MiB (1048576) is the sweet spot for HTTP-
    /// served assets; reduce for tiny example assets so the test
    /// fans out across multiple blocks.
    #[arg(long, default_value_t = 1 << 20)]
    target_block_bytes: u32,

    /// Skip emission of `ro-crate-metadata.json` and the `schemas/`
    /// directory when writing a `--mode dafsa-blocks` asset. Default
    /// is to write a self-describing RO-Crate 1.2 bundle alongside
    /// the manifest and block files; this flag is the escape hatch
    /// for callers that want only the bare wire format (e.g. an
    /// internal test fixture).
    #[arg(long)]
    no_rocrate: bool,

    /// Optional OEIS A-number this enumeration realises (e.g.
    /// `A316192` for ZZ12 free dihedral). Surfaced into the
    /// `ro-crate-metadata.json` as a `subjectOf` contextual entity
    /// so an archivist can pivot from the asset to the OEIS
    /// sequence. Honoured by the asset-writing modes (`--mode
    /// dafsa-blocks` and the streaming pipeline's `--mode build`).
    #[arg(long)]
    oeis_a_number: Option<String>,

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
        install_mod_prune(cli.ring, cli.max_steps, cli.mod_prune_moduli.as_deref());
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
            cli.free,
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
    // ~N x higher because this path still uses the HashSet-backed
    // `enumerate_from_seed` (i.e. it's bench-shaped, not stream-
    // shaped). For memory-bounded large-n enumeration, prefer
    // `--mode stream` instead, where multi-threading within one
    // process and multi-process orchestration have ~the same
    // memory profile (no per-thread result accumulation to
    // multiply).
    if let Some(seed) = cli.seed.as_deref() {
        let t0 = Instant::now();
        let rats = dispatch_enumerate_from_seed(
            cli.ring,
            cli.max_steps,
            cli.step,
            seed,
            n_threads,
            cli.free,
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
        Mode::Build => {
            let Some(out_dir) = cli.filename.as_deref() else {
                eprintln!("--mode build requires -o <output directory>");
                std::process::exit(2);
            };
            let out_dir = std::path::Path::new(out_dir);
            let unique_path = out_dir.join(tilezz::rat_enum::stream::UNIQUE_FILENAME);
            if !unique_path.exists() {
                eprintln!(
                    "--mode build: {} not found; run `--mode merge` first",
                    unique_path.display()
                );
                std::process::exit(2);
            }

            let t0 = Instant::now();
            // Stream records straight into the streaming DAFSA
            // builder. Peak RSS scales with the in-progress
            // automaton, not with the input set.
            let records = tilezz::rat_enum::stream::read_unique_records(&unique_path)
                .expect("open unique.bin")
                .map(|r| r.expect("read unique record"));
            let dafsa = RatDafsa::from_sorted_unique_rats(records);
            let t_build = t0.elapsed();
            eprintln!(
                "build: streamed {} rats into RatDafsa in {:?}",
                dafsa.len(),
                t_build
            );

            let blocks_dir = out_dir.join("dafsa");
            std::fs::create_dir_all(&blocks_dir).expect("create dafsa/");
            let t1 = Instant::now();
            dafsa
                .write_blocks(&blocks_dir, cli.target_block_bytes)
                .expect("write_blocks");

            // Same self-describing bundle as Mode::DafsaBlocks, but the
            // reproduction recipe records the streaming pipeline's
            // three-stage shape instead of a single --mode dafsa-blocks
            // invocation.
            if !cli.no_rocrate {
                use tilezz::stringmatch::dafsa::{
                    AssetParams, ProducedVia, write_archival_extras, write_ro_crate,
                };
                let params = AssetParams {
                    ring: cli.ring,
                    max_steps: cli.max_steps,
                    step: cli.step,
                    free: cli.free,
                    target_block_bytes: cli.target_block_bytes,
                    n_sequences: dafsa.len() as u64,
                    oeis_a_number: cli.oeis_a_number.as_deref(),
                    produced_via: ProducedVia::StreamingPipeline,
                };
                write_archival_extras(&blocks_dir, &params).expect("write archival extras");
                write_ro_crate(&blocks_dir, &params).expect("write ro-crate-metadata.json");
            }

            fn dir_size_recursive(p: &std::path::Path) -> u64 {
                let mut total = 0u64;
                if let Ok(rd) = std::fs::read_dir(p) {
                    for entry in rd.flatten() {
                        let path = entry.path();
                        if path.is_dir() {
                            total += dir_size_recursive(&path);
                        } else if let Ok(m) = entry.metadata() {
                            total += m.len();
                        }
                    }
                }
                total
            }
            let total_bytes = dir_size_recursive(&blocks_dir);
            eprintln!(
                "build: wrote {} ({} bytes) in {:?}",
                blocks_dir.display(),
                total_bytes,
                t1.elapsed()
            );
        }
        Mode::Merge => {
            let Some(out_dir) = cli.filename.as_deref() else {
                eprintln!("--mode merge requires -o <output directory>");
                std::process::exit(2);
            };
            let t0 = Instant::now();
            let cert = tilezz::rat_enum::stream::merge_runs(
                std::path::Path::new(out_dir),
                cli.ring,
                cli.max_steps,
                cli.step,
                cli.free,
            )
            .expect("merge_runs");
            let dt = t0.elapsed();
            println!(
                "merge: ring={} step={} max_steps={} -> {} unique rats in {:?}",
                cli.ring, cli.step, cli.max_steps, cert.unique_records, dt
            );
        }
        Mode::Stream => {
            let Some(out_dir) = cli.filename.as_deref() else {
                eprintln!("--mode stream requires -o <output directory>");
                std::process::exit(2);
            };
            // Stage 1 is artifact-producing; don't double up by also
            // dumping RAT lines on stdout.
            STREAM_RAT_LINES.store(false, std::sync::atomic::Ordering::Relaxed);

            let t0 = Instant::now();
            let stats = tilezz::rat_enum::stream::stream_enum_dispatch(
                cli.ring,
                cli.max_steps,
                cli.step,
                n_threads,
                cli.free,
                cli.paranoid,
                std::path::Path::new(out_dir),
            )
            .expect("stream_enum_dispatch");
            let dt = t0.elapsed();
            println!(
                "stream: ring={} step={} max_steps={} -> wrote runs to {} in {:?}",
                cli.ring, cli.step, cli.max_steps, out_dir, dt
            );
            println!("{stats}");
        }
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
                cli.free,
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
        Mode::DafsaBlocks => {
            let Some(dir) = cli.filename.as_deref() else {
                eprintln!("--mode dafsa-blocks requires -o <output directory>");
                std::process::exit(2);
            };

            STREAM_RAT_LINES.store(false, std::sync::atomic::Ordering::Relaxed);

            let t0 = Instant::now();
            let (rats, _stats) = run_rat_enum_seqs(
                cli.ring,
                cli.max_steps,
                cli.step,
                n_threads,
                cli.free,
                cli.paranoid,
            );
            eprintln!("enumerated {} rats in {:?}", rats.len(), t0.elapsed());

            let t1 = Instant::now();
            let dafsa = RatDafsa::from_rats(rats.iter().map(|r| r.as_slice()));
            eprintln!(
                "built RatDafsa ({} entries) in {:?}",
                dafsa.len(),
                t1.elapsed()
            );

            let t2 = Instant::now();
            let path = std::path::Path::new(dir);
            std::fs::create_dir_all(path).expect("create output dir");
            dafsa
                .write_blocks(path, cli.target_block_bytes)
                .expect("write blocked RatDafsa");

            // Self-describing RO-Crate bundle: schemas/ + tools/decode.py
            // + REPRODUCE.md + ro-crate-metadata.json. Skipped under
            // --no-rocrate (e.g. internal test fixtures that want only
            // the bare wire format).
            if !cli.no_rocrate {
                use tilezz::stringmatch::dafsa::{
                    AssetParams, ProducedVia, write_archival_extras, write_ro_crate,
                };
                let params = AssetParams {
                    ring: cli.ring,
                    max_steps: cli.max_steps,
                    step: cli.step,
                    free: cli.free,
                    target_block_bytes: cli.target_block_bytes,
                    n_sequences: dafsa.len() as u64,
                    oeis_a_number: cli.oeis_a_number.as_deref(),
                    produced_via: ProducedVia::InMemory,
                };
                write_archival_extras(path, &params).expect("write archival extras");
                write_ro_crate(path, &params).expect("write ro-crate-metadata.json");
            }

            // Walk the directory tree so the total includes blocks/
            // and schemas/ subdir contents, not just top-level files.
            fn dir_size(p: &std::path::Path) -> u64 {
                let mut total = 0u64;
                if let Ok(rd) = std::fs::read_dir(p) {
                    for entry in rd.flatten() {
                        let path = entry.path();
                        if path.is_dir() {
                            total += dir_size(&path);
                        } else if let Ok(m) = entry.metadata() {
                            total += m.len();
                        }
                    }
                }
                total
            }
            let total_bytes = dir_size(path);
            let extras = if cli.no_rocrate {
                ""
            } else {
                " + schemas + ro-crate-metadata"
            };
            eprintln!(
                "wrote {dir}/ ({total_bytes} bytes across manifest + blocks{extras}) in {:?}",
                t2.elapsed()
            );

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
                cli.free,
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
                cli.free,
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
mod free_tests {
    use super::*;
    use std::collections::HashMap;
    use tilezz::cyclotomic::geometry::intersect;
    use tilezz::cyclotomic::{IsRing, Units, ZZ4, ZZ8, ZZ12};

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

    /// Per-boundary-length counts the ZZ12 free enumeration must
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
    #[cfg_attr(
        debug_assertions,
        ignore = "release-only: full A316192 pin to n>=10 takes ~12 min debug / ~30 s release"
    )]
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
        // walks the free variant; output is bucketed
        // per perimeter length.
        let max_n = OEIS.iter().map(|&(n, _)| n).max().unwrap();
        let (rats, _) = super::rat_enum_free::<ZZ12>(max_n, 1);

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
                "free ZZ12 enumeration differs from OEIS A316192 at {} perimeter length(s)",
                mismatches.len()
            );
        }
    }

    #[test]
    fn test_free_output_polygons_are_simple_and_unique() {
        let (free_rats, _) = super::rat_enum_free::<ZZ12>(9, 1);

        // (i) Encodings are pairwise distinct (the HashSet via which
        // they were collected enforces this; assert defensively).
        let unique: std::collections::HashSet<Vec<i8>> = free_rats.iter().cloned().collect();
        assert_eq!(
            unique.len(),
            free_rats.len(),
            "duplicate sequences in free output"
        );

        // (ii) Every polygon is geometrically simple per the
        // independent validator above. Collect failures so a single
        // run reports every case.
        let mut failures: Vec<(Vec<i8>, String)> = Vec::new();
        for seq in &free_rats {
            if let Err(why) = validate_simple_polygon::<ZZ12>(seq) {
                failures.push((seq.clone(), why));
            }
        }
        if !failures.is_empty() {
            eprintln!("{} polygons failed independent validation:", failures.len());
            for (seq, why) in failures.iter().take(20) {
                eprintln!("  {seq:?} -- {why}");
            }
            panic!("non-simple polygons in free enumeration output");
        }
    }

    /// For a single `(ring, max_steps)` pair, check that the
    /// free DFS output equals the rotation DFS output quotiented
    /// by `free_canonical`. Stage-1 over-pruning shows up as
    /// MISSING free classes; stage-2 under-deduplication shows
    /// up as EXTRA. The body prints both lists before asserting so
    /// any future regression is diagnosable from the test log.
    fn check_quotient_match<ZZ: IsRing>(ring_label: &str, max_steps: usize) {
        let (all_rats, _) = super::rat_enum::<ZZ>(max_steps, 1);
        let mut expected: HashSet<Vec<i8>> = HashSet::new();
        for seq in &all_rats {
            expected.insert(super::free_canonical(seq));
        }
        let (free_rats, _) = super::rat_enum_free::<ZZ>(max_steps, 1);
        let actual: HashSet<Vec<i8>> = free_rats.into_iter().collect();

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
            "{ring_label} n={max_steps}: free DFS != rotation DFS / free",
        );
        eprintln!(
            "[{ring_label} n={max_steps}] OK -- {} free classes from {} rotation-canonical rats",
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
    #[cfg_attr(
        debug_assertions,
        ignore = "release-only: ZZ12 n<=9 + ZZ8 n<=12 dihedral cross-check is ~9 min debug / ~25 s release"
    )]
    fn test_free_enum_matches_dfs_quotient() {
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
    ///   * free on AND off (each uses a different `CanonicalOps`
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
            for &free in &[false, true] {
                for &split_depth in &[1usize, 2, 3] {
                    if split_depth >= n {
                        continue;
                    }
                    let (one_shot_seqs, _) = if free {
                        super::rat_enum_free::<ZZ12>(n, 1)
                    } else {
                        super::rat_enum::<ZZ12>(n, 1)
                    };
                    let one_shot: HashSet<Vec<i8>> = one_shot_seqs.into_iter().collect();

                    let (mut seeded, prefixes) = tilezz::rat_enum::seed::collect_seed_prefixes::<
                        ZZ12,
                    >(n, 1, split_depth, free);
                    for prefix in &prefixes {
                        // Sweep both branches of the n_threads dispatch
                        // (single-threaded direct path AND sub-split +
                        // parallel_drain_seeds) within the same test --
                        // they must produce the same set.
                        for nthreads in &[1usize, 4] {
                            seeded.extend(tilezz::rat_enum::seed::enumerate_from_seed::<ZZ12>(
                                n, 1, prefix, *nthreads, free, false,
                            ));
                        }
                    }

                    assert_eq!(
                        one_shot,
                        seeded,
                        "n={n} free={free} split_depth={split_depth}: \
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
/// For each ring x mode (rotation/free) x thread-count, the test
/// matrix runs the enumeration under every subset of the optional
/// prunes (`mod`, `closure-key`) and asserts the resulting
/// canonical-rat set is bit-identical to the baseline (no prunes). A
/// regression in any single prune -- or any composition -- shows up
/// as a set-equality mismatch on these tests.
///
/// Caps: ring n is chosen so the full sweep stays under ~30 s on a
/// release build. ZZ12 stops at n=8 (no T-touch bug); ZZ8 / ZZ4 go
/// further within their compute budgets.
#[cfg(test)]
mod opt_correctness_tests {
    use super::*;
    use std::sync::Arc;
    use tilezz::cyclotomic::{
        ZZ4, ZZ6, ZZ8, ZZ10, ZZ12, ZZ14, ZZ16, ZZ18, ZZ20, ZZ24, ZZ32, ZZ60,
    };
    use tilezz::rat_enum::prune::closure_key::{ClosureKeyPrune, collect_closure_keys};
    use tilezz::rat_enum::prune::modular::ModularPrune;
    use tilezz::rat_enum::prune::units::unit_vectors_for_ring;
    use tilezz::rat_enum::seed::rat_enum_parallel;

    /// Dispatch helper: collect closure keys for `ring` up to length
    /// `max_l`. Mirrors the per-ring match in `install_closure_key_prune`.
    fn closure_keys_for(ring: u8, max_l: usize) -> rustc_hash::FxHashSet<(Vec<i64>, i8)> {
        tilezz::dispatch_ring!(ring, collect_closure_keys::<ZZ>(max_l))
    }

    /// Build a `Prunes` for testing per the chosen opt subset.
    /// Both prunes are constructed against `(ring, max_steps)`
    /// regardless of whether they're enabled, so each test case is
    /// self-contained and independent of any global state.
    fn build_prunes(ring: u8, max_steps: usize, with_mod: bool, with_ck: bool) -> Prunes {
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
        free: bool,
        n_threads: usize,
        prunes: &Prunes,
    ) -> std::collections::HashSet<Vec<i8>> {
        let ops = make_ops(free);
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
        for &free in &[false, true] {
            let baseline = run::<ZZ>(max_steps, free, 1, &Prunes::default());
            for (with_mod, with_ck) in opt_subsets() {
                let prunes = build_prunes(ring, max_steps, with_mod, with_ck);
                for &n_threads in &[1usize, 4] {
                    let got = run::<ZZ>(max_steps, free, n_threads, &prunes);
                    assert_eq!(
                        got,
                        baseline,
                        "ZZ{ring} n={max_steps} free={free} \
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
        // covers ZZ12 free up to n=10 for every opt combo,
        // which is the bug-sensitive range (T-touches first occur
        // at n=9 -- see the historical note in test_oeis_a316192_zz12).
        check_ring::<ZZ12>(12, 8);
    }

    #[test]
    #[cfg_attr(
        debug_assertions,
        ignore = "release-only: ZZ14 n=7 32-prune-combo cross-check is ~5 min debug / ~15 s release"
    )]
    fn cross_validate_zz14() {
        // Capped at n=7 for runtime budget. ZZ14 has no OEIS oracle
        // for the full-ring counts (A316197 is the bipartite subset),
        // so cross-prune correctness is the strongest independent
        // check we have. Even at n=7, the 32-combo matrix verifies
        // that the cubic-root-based real_sign + multivariate-sign
        // CellFloor produce a consistent rat set across all prune
        // subsets.
        check_ring::<ZZ14>(14, 7);
    }

    #[test]
    #[cfg_attr(
        debug_assertions,
        ignore = "release-only: ZZ18 n=7 32-prune-combo cross-check is ~13 min debug / ~40 s release"
    )]
    fn cross_validate_zz18() {
        // ZZ18 shares the cubic-root sign infrastructure with ZZ14
        // (different minpoly, same algorithm). Same runtime-budget
        // cap; no OEIS oracle for full-ring counts (A316199 is the
        // bipartite subset).
        check_ring::<ZZ18>(18, 7);
    }

    /// **External verification of ZZ18's sign infrastructure** via a
    /// completely independent code path. ZZ18 contains ZZ6 (since
    /// `6 | 18`), and `--step 3` on ZZ18 restricts the DFS to turns
    /// that are multiples of `pi/3` -- which is exactly ZZ6's unit
    /// turn angle. So ZZ18-step-3 enumerates the same polygons as
    /// ZZ6-step-1, but through ZZ18's algebraic machinery (cubic-
    /// root sign helper, multivariate sign helper, custom CellFloor).
    ///
    /// ZZ6 is OEIS-pinned to A284869 (12 published terms verified by
    /// `oeis_a284869_zz6_pin`). If ZZ18-step-3 matches ZZ6, that
    /// transitively verifies ZZ18's sign machinery against external
    /// data -- otherwise the only oracle for ZZ18 is internal cross-
    /// prune consistency, which uses the same sign code throughout
    /// and so couldn't catch a sign-helper bug.
    ///
    /// Runs to perim 8 (~10 ms each) covering the bug-sensitive
    /// range where our sign infrastructure differs structurally
    /// from the simpler ZZ6 (which uses `sign_m_plus_n_sqrt3`).
    #[test]
    fn cross_validate_zz18_step3_matches_zz6() {
        // ZZ6 reference enumeration (uses `sign_m_plus_n_sqrt3`,
        // which is independently OEIS-verified).
        let ops_zz6 = make_ops(true);
        let (zz6_rats, _) =
            rat_enum_with::<ZZ6>(8, 1, ops_zz6, "zz6_ref", "", false, &Prunes::default());
        let mut zz6_by_len: std::collections::BTreeMap<usize, usize> =
            std::collections::BTreeMap::new();
        for r in &zz6_rats {
            *zz6_by_len.entry(r.len()).or_insert(0) += 1;
        }

        // ZZ18 step=3 enumeration (uses our new cubic-root sign +
        // multivariate sign + custom CellFloor for the imag axis).
        let ops_zz18 = make_ops(true);
        let (zz18_rats, _) =
            rat_enum_with::<ZZ18>(8, 3, ops_zz18, "zz18_step3", "", false, &Prunes::default());
        let mut zz18_by_len: std::collections::BTreeMap<usize, usize> =
            std::collections::BTreeMap::new();
        for r in &zz18_rats {
            *zz18_by_len.entry(r.len()).or_insert(0) += 1;
        }

        // Compare per-length counts.
        let mut all_lens: std::collections::BTreeSet<usize> = zz6_by_len.keys().copied().collect();
        all_lens.extend(zz18_by_len.keys().copied());
        let mut mismatches: Vec<(usize, usize, usize)> = Vec::new();
        for &n in &all_lens {
            let z6 = zz6_by_len.get(&n).copied().unwrap_or(0);
            let z18 = zz18_by_len.get(&n).copied().unwrap_or(0);
            if z6 != z18 {
                mismatches.push((n, z6, z18));
            }
        }
        if !mismatches.is_empty() {
            for (n, z6, z18) in &mismatches {
                eprintln!(
                    "perim={n}: ZZ6={z6}, ZZ18-step-3={z18}, diff={:+}",
                    *z18 as i64 - *z6 as i64
                );
            }
            panic!(
                "ZZ18-step-3 disagrees with ZZ6 at {} perimeter(s) -- sign helper or \
                 cell_floor regression in ZZ18",
                mismatches.len()
            );
        }
    }

    /// Generic step-subset cross-check. `ZZ_big --step k` walks only
    /// turns that are multiples of `k`; when `k * small = big` that
    /// enumerates exactly the `ZZ_small`-equivalent polygons -- a
    /// bijection on turn sequences, so the per-perimeter counts must
    /// match. This drives the BIG ring's own sign + cell_floor
    /// machinery (the i128 nested-sqrt helpers, for the deep rings)
    /// yet compares against a ring that is OEIS-pinned, transitively
    /// anchoring the big ring to external data with no brute-force
    /// oracle. `--step k` collapses the effective branching to ~n/k
    /// directions, so at the small `max_steps` used here these stay
    /// in the fast default suite (each well under a second).
    fn assert_step_subset<Big: IsRing, Small: IsRing>(max_steps: usize, step: i8, label: &str) {
        let count_by_len = |rats: &[Vec<i8>]| {
            let mut m = std::collections::BTreeMap::<usize, usize>::new();
            for r in rats {
                *m.entry(r.len()).or_insert(0) += 1;
            }
            m
        };
        let (big, _) = rat_enum_with::<Big>(
            max_steps,
            step,
            make_ops(true),
            label,
            "",
            false,
            &Prunes::default(),
        );
        let (small, _) = rat_enum_with::<Small>(
            max_steps,
            1,
            make_ops(true),
            "ref",
            "",
            false,
            &Prunes::default(),
        );
        assert_eq!(
            count_by_len(&big),
            count_by_len(&small),
            "{label}: step-{step} subset counts disagree with the reference ring \
             -- a sign-helper or cell_floor bug in the bigger ring",
        );
    }

    /// Step-subset anchors for the rings WITHOUT their own OEIS
    /// sequence: each reduces (via `--step`) to an OEIS-pinned ring,
    /// transitively verifying its exact-geometry machinery. ZZ14
    /// (= 2*7) is the sole ring with no such reduction and is checked
    /// separately. Small n keeps the suite fast; a deeper sweep lives
    /// in the `#[ignore]` release tests.
    #[test]
    fn zz16_step2_matches_zz8() {
        assert_step_subset::<ZZ16, ZZ8>(8, 2, "zz16_step2");
    }
    #[test]
    fn zz20_step2_matches_zz10() {
        assert_step_subset::<ZZ20, ZZ10>(7, 2, "zz20_step2");
    }
    #[test]
    fn zz24_step2_matches_zz12() {
        assert_step_subset::<ZZ24, ZZ12>(7, 2, "zz24_step2");
    }
    #[test]
    fn zz32_step4_matches_zz8() {
        assert_step_subset::<ZZ32, ZZ8>(6, 4, "zz32_step4");
    }
    #[test]
    fn zz60_step5_matches_zz12() {
        assert_step_subset::<ZZ60, ZZ12>(6, 5, "zz60_step5");
    }

    /// External anchor: OEIS A316192 per-length counts for ZZ12
    /// free matchstick polygons, n=3..=10. The range
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
    #[cfg_attr(
        debug_assertions,
        ignore = "release-only: A316192 pin x 32 prune combos is ~17 min debug / ~60-90 s release"
    )]
    fn oeis_a316192_each_opt_combo() {
        const OEIS: &[(usize, usize)] = &[
            (3, 1),
            (4, 3),
            (5, 4),
            (6, 22),
            (7, 69),
            (8, 418),
            (9, 2210),   // T-touches first appear here; bug-fix-sensitive
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

    /// Frontier regression guard for the A316192 EXTENSION terms
    /// (ZZ12 free, n >= 11) -- the ones OEIS does not publish and that
    /// we would submit, so they have no external oracle. Two oracle-
    /// free checks:
    ///   * prune-invariance at n=11: the optional prunes (mod +
    ///     closure-key) must not change the free count vs the
    ///     baseline (free/canonical/reachability prunes only). An
    ///     over-aggressive prune that dropped a valid rat would
    ///     silently undercount -- this catches it one step past the
    ///     OEIS-pinned range.
    ///   * value pins for a(11)=89075 and a(12)=597581, locking the
    ///     first extension terms against any future silent drift.
    /// n<=12 only: the unpruned baseline at n>=13 is too slow even
    /// for a release test; deeper terms rely on the prunes'
    /// analytical (n-independent) soundness, the verification gate,
    /// and an independent recompute at submission time (see
    /// docs/oeis-A316200-correction for the recompute methodology).
    #[test]
    #[ignore = "opt-in (cargo test -- --include-ignored): ~4 min; pre-submission / thorough CI \
                guard, not part of the default debug or release suite"]
    fn zz12_extension_frontier_guard() {
        // Prune-invariance at n=11 (optional prunes on vs off).
        let none = build_prunes(12, 11, false, false);
        let all11 = build_prunes(12, 11, true, true);
        let base = by_length(&run::<ZZ12>(11, true, 0, &none));
        let pruned = by_length(&run::<ZZ12>(11, true, 0, &all11));
        assert_eq!(
            base, pruned,
            "ZZ12 free n=11: optional prunes changed the count -- over-pruning at the frontier",
        );
        // Value pins for the first two extension terms.
        let all12 = build_prunes(12, 12, true, true);
        let c12 = by_length(&run::<ZZ12>(12, true, 0, &all12));
        assert_eq!(c12.get(&11).copied(), Some(89075), "ZZ12 a(11) regression");
        assert_eq!(c12.get(&12).copied(), Some(597581), "ZZ12 a(12) regression");
    }

    /// Group an enumerated `Vec<Vec<i8>>` into `len -> count`.
    fn by_length(
        rats: &std::collections::HashSet<Vec<i8>>,
    ) -> std::collections::BTreeMap<usize, usize> {
        let mut by_len = std::collections::BTreeMap::new();
        for seq in rats {
            *by_len.entry(seq.len()).or_insert(0) += 1;
        }
        by_len
    }

    /// Shared pin assertion: enumerate ring `ZZ` in free mode up to
    /// `max_n` with all prunes enabled (single-threaded for
    /// determinism), then assert each `(perim, expected)` pair from
    /// `oeis`. Used by the cheap external-anchor pins below.
    fn assert_oeis_pins<ZZ: IsRing>(ring: u8, oeis_name: &str, oeis: &[(usize, usize)]) {
        let max_n = oeis.iter().map(|&(n, _)| n).max().unwrap();
        let prunes = build_prunes(ring, max_n, true, true);
        let rats = run::<ZZ>(max_n, true, 1, &prunes);
        let by_len = by_length(&rats);
        let mut mismatches: Vec<(usize, usize, usize)> = Vec::new();
        for &(n, expected) in oeis {
            let got = by_len.get(&n).copied().unwrap_or(0);
            if got != expected {
                mismatches.push((n, got, expected));
            }
        }
        if !mismatches.is_empty() {
            for (n, got, expected) in &mismatches {
                eprintln!(
                    "{oeis_name} ZZ{ring} perim={n}: got {got}, OEIS says {expected}, diff {:+}",
                    *got as i64 - *expected as i64
                );
            }
            panic!("{oeis_name} ZZ{ring}: {} term(s) differ", mismatches.len());
        }
    }

    /// External anchor: OEIS A266549 per-length counts for ZZ4
    /// free matchstick polygons on the square lattice. Indexed by
    /// perim 2n (square-lattice polygons only close at even perim).
    /// Independently confirmed by our enumeration through perim 14;
    /// see `docs/oeis.md` for the full overview.
    #[test]
    fn oeis_a266549_zz4_pin() {
        const OEIS: &[(usize, usize)] = &[
            (4, 1),   // a(2)
            (6, 1),   // a(3)
            (8, 3),   // a(4)
            (10, 6),  // a(5)
            (12, 25), // a(6)
            (14, 86), // a(7)
        ];
        assert_oeis_pins::<ZZ4>(4, "A266549", OEIS);
    }

    /// External anchor: OEIS A316198 per-length counts for ZZ8
    /// free matchstick polygons. Indexed by perim 2n. Hugo's
    /// data covers perim 4..12; our perim 14 = 240549 is a new term
    /// we can contribute back to OEIS (see `docs/oeis.md`).
    #[test]
    fn oeis_a316198_zz8_pin() {
        const OEIS: &[(usize, usize)] = &[
            (4, 2),      // a(2)
            (6, 6),      // a(3)
            (8, 59),     // a(4)
            (10, 695),   // a(5)
            (12, 12198), // a(6)
        ];
        assert_oeis_pins::<ZZ8>(8, "A316198", OEIS);
    }

    /// External anchor: OEIS A316200 per-length counts for ZZ10
    /// free matchstick polygons. Pins perim 4..10 ONLY: at
    /// perim 11 our enumeration disagrees with OEIS (we say 9883,
    /// OEIS says 19405) and we believe the OEIS term is in error.
    /// Five independent enumeration paths (ring 10 step 1, ring 20
    /// step 2, single/multi-threaded, with/without prunes, plus the
    /// full streaming pipeline) all return 9883. See `docs/oeis.md`
    /// for the discrepancy analysis.
    #[test]
    fn oeis_a316200_zz10_pin() {
        const OEIS: &[(usize, usize)] = &[
            (4, 2),   // a(4)
            (5, 2),   // a(5)
            (6, 10),  // a(6)
            (7, 15),  // a(7)
            (8, 124), // a(8)
            (9, 352), // a(9)
            (10, 2378), // a(10)
                      // a(11): OEIS says 19405, we say 9883. Pin omitted
                      // until the discrepancy is resolved with OEIS.
        ];
        assert_oeis_pins::<ZZ10>(10, "A316200", OEIS);
    }

    /// External anchor: OEIS A284869 per-length counts for ZZ6
    /// free matchstick polygons on the triangular lattice. This
    /// is the deepest external oracle in our coverage -- A284869
    /// publishes a(3..24), giving us 22 verification anchors against
    /// an independent source. We pin perim 3..14 here (matches first
    /// 12 of the 22 terms exactly) since enumeration past n=14 grows
    /// by ~3-4x per step and becomes test-runtime expensive; deeper
    /// terms are spot-checked manually rather than under CI.
    #[test]
    fn oeis_a284869_zz6_pin() {
        const OEIS: &[(usize, usize)] = &[
            (3, 1),      // a(3)
            (4, 1),      // a(4)
            (5, 1),      // a(5)
            (6, 4),      // a(6)
            (7, 5),      // a(7)
            (8, 16),     // a(8)
            (9, 37),     // a(9)
            (10, 120),   // a(10)
            (11, 344),   // a(11)
            (12, 1175),  // a(12)
            (13, 3807),  // a(13)
            (14, 13224), // a(14)
        ];
        assert_oeis_pins::<ZZ6>(6, "A284869", OEIS);
    }
}
