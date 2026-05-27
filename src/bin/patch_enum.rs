use std::collections::BTreeMap;
use std::sync::Arc;
use std::time::Instant;

use clap::Parser;
use rustc_hash::FxHashSet;

use tilezz::analysis::patch_enum::enum_patches;
use tilezz::cyclotomic::{IsRing, ZZ12, ZZ4};
use tilezz::geom::rat::Rat;
use tilezz::geom::tileset::{self, TileSet};

#[derive(Parser)]
#[command(
    name = "patch_enum",
    about = "Enumerate tile patches by incremental growth"
)]
struct Args {
    /// Tileset to grow patches from.
    #[arg(long, value_enum, default_value = "hex")]
    tileset: TileSetKind,

    /// Maximum patch size (number of tile copies).
    #[arg(long, default_value_t = 6)]
    max_size: usize,

    /// Cyclotomic ring. `zz4` only supports square / tetrominoes.
    #[arg(long, value_enum, default_value = "zz12")]
    ring: RingChoice,

    /// Cross-check against a brute-force `try_glue` enumeration.
    /// Slow for `max_size > 4`.
    #[arg(long)]
    validate: bool,

    /// Flamegraph output path (requires --features pprof).
    #[arg(long)]
    pprof: Option<String>,
}

#[derive(Clone, Debug, clap::ValueEnum)]
enum TileSetKind {
    Hex,
    Square,
    Mixed,
    Tetris,
    Spectre,
}

#[derive(Clone, Debug, clap::ValueEnum)]
enum RingChoice {
    ZZ4,
    ZZ12,
}

fn ts_zz4(kind: &TileSetKind) -> Arc<TileSet<ZZ4>> {
    match kind {
        TileSetKind::Square => tileset::square::<ZZ4>(),
        TileSetKind::Tetris => tileset::tetrominoes::<ZZ4>(),
        other => {
            eprintln!("error: {other:?} not supported with --ring zz4 (need ring with ZZ4 only)");
            std::process::exit(1);
        }
    }
}

fn ts_zz12(kind: &TileSetKind) -> Arc<TileSet<ZZ12>> {
    match kind {
        TileSetKind::Hex => tileset::hex::<ZZ12>(),
        TileSetKind::Square => tileset::square::<ZZ12>(),
        TileSetKind::Mixed => tileset::mixed::<ZZ12>(),
        TileSetKind::Tetris => tileset::tetrominoes::<ZZ12>(),
        TileSetKind::Spectre => tileset::spectre::<ZZ12>(),
    }
}

/// Brute-force reference enumeration: glue every
/// `(patch, tile, ia, ib)` combination via `try_glue`. No
/// canonicalisation shortcuts, just the most direct path possible.
///
/// Outer-loop parallelism only: each patch's gluing attempts run on a
/// worker thread, then the per-worker `FxHashSet`s are merged. The
/// nested `(tile_idx, ia, ib)` loops and the `try_glue` call stay
/// pointwise identical to the single-threaded brute force, so the set
/// of attempted glue operations is unchanged.
fn brute_force_grow<T>(
    tileset: &Arc<TileSet<T>>,
    max_size: usize,
) -> BTreeMap<usize, FxHashSet<Rat<T>>>
where
    T: IsRing + Send + Sync,
{
    let mut results: BTreeMap<usize, FxHashSet<Rat<T>>> = BTreeMap::new();
    if max_size == 0 || tileset.num_tiles() == 0 {
        return results;
    }
    results.insert(1, tileset.rats().iter().cloned().collect());

    let n_threads = num_cpus::get().max(1);

    for k in 2..=max_size {
        let prev: Vec<Rat<T>> = results[&(k - 1)].iter().cloned().collect();
        let chunk_size = prev.len().div_ceil(n_threads).max(1);

        let next: FxHashSet<Rat<T>> = std::thread::scope(|s| {
            let mut handles = Vec::with_capacity(n_threads);
            for chunk in prev.chunks(chunk_size) {
                handles.push(s.spawn(move || {
                    let mut local: FxHashSet<Rat<T>> = FxHashSet::default();
                    for patch in chunk {
                        for tile_idx in 0..tileset.num_tiles() {
                            let tile = tileset.rat(tile_idx);
                            for ia in 0..patch.len() {
                                for ib in 0..tile.len() {
                                    if let Ok(glued) = patch.try_glue((ia as i64, ib as i64), tile)
                                    {
                                        local.insert(glued);
                                    }
                                }
                            }
                        }
                    }
                    local
                }));
            }
            let mut merged: FxHashSet<Rat<T>> = FxHashSet::default();
            for h in handles {
                merged.extend(h.join().expect("brute-force worker panic"));
            }
            merged
        });

        results.insert(k, next);
    }
    results
}

fn run<T>(tileset: Arc<TileSet<T>>, max_size: usize, validate: bool, label: &str)
where
    T: IsRing + Send + Sync,
{
    eprintln!(
        "tileset: {label} ({} tile{}), max_size: {max_size}, validate: {validate}",
        tileset.num_tiles(),
        if tileset.num_tiles() == 1 { "" } else { "s" },
    );

    let t0 = Instant::now();
    let patches = enum_patches(Arc::clone(&tileset), max_size);
    let elapsed = t0.elapsed();

    if validate {
        let t_bf = Instant::now();
        let bf = brute_force_grow(&tileset, max_size);
        let dt_bf = t_bf.elapsed();
        let mut ok = true;
        for k in 1..=max_size {
            let fast = patches.get(&k).cloned().unwrap_or_default();
            let brute = bf.get(&k).cloned().unwrap_or_default();
            if fast != brute {
                eprintln!(
                    "  size {k}: MISMATCH (enum_patches={}, brute_force={})",
                    fast.len(),
                    brute.len()
                );
                for r in fast.difference(&brute).take(5) {
                    eprintln!("    only in enum_patches: {:?}", r.seq());
                }
                for r in brute.difference(&fast).take(5) {
                    eprintln!("    only in brute force: {:?}", r.seq());
                }
                ok = false;
            }
        }
        if !ok {
            std::process::exit(1);
        }
        let speedup = dt_bf.as_secs_f64() / elapsed.as_secs_f64().max(1e-9);
        eprintln!(
            "  validate: OK | enum_patches: {elapsed:.2?} | brute_force: {dt_bf:.2?} | {speedup:.1}x speedup"
        );
    } else {
        eprintln!("  enum_patches: {elapsed:.2?}");
    }

    eprintln!("\n--- Summary ---");
    let mut total = 0;
    for k in 1..=max_size {
        let count = patches.get(&k).map(|s| s.len()).unwrap_or(0);
        eprintln!("  size {k:>3}: {count:>7} patches");
        total += count;
    }
    eprintln!("  total:    {total:>7} patches");
}

fn main() {
    let args = Args::parse();
    if args.max_size < 1 {
        eprintln!("max_size must be at least 1");
        std::process::exit(1);
    }
    if args.validate && args.max_size > 4 {
        eprintln!("note: --validate with max_size > 4 may be very slow");
    }

    let profile = tilezz::util::profile::ProfileGuard::start(args.pprof.as_deref());

    let label = format!("{:?}", args.tileset);
    match args.ring {
        RingChoice::ZZ4 => run(ts_zz4(&args.tileset), args.max_size, args.validate, &label),
        RingChoice::ZZ12 => run(ts_zz12(&args.tileset), args.max_size, args.validate, &label),
    }

    profile.finish();
}
