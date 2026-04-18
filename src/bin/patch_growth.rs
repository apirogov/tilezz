use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

use clap::Parser;

use tilezz::cyclotomic::{IsComplex, IsRingOrField, Units, ZZ12, ZZ4};
use tilezz::intgeom::rat::Rat;
use tilezz::intgeom::snake::Snake;
use tilezz::intgeom::tiles;
use tilezz::intgeom::tileset::{GlueStats, TileSet};

#[derive(Parser)]
#[command(
    name = "patch_growth",
    about = "Enumerate tile patches by incremental growth"
)]
struct Args {
    /// Seed tile shape
    #[arg(long, value_enum, default_value = "hexagon")]
    tile: TileShape,

    /// Maximum patch size (number of tiles)
    #[arg(long, default_value_t = 6)]
    max_size: usize,

    /// Cyclotomic ring (ZZ4 only supports square tiles)
    #[arg(long, value_enum, default_value = "zz12")]
    ring: RingChoice,

    /// Print per-pair glue statistics
    #[arg(long)]
    verbose: bool,

    /// Validate results against brute force (caps max_size at 4)
    #[arg(long)]
    validate: bool,
}

#[derive(Clone, Debug, clap::ValueEnum)]
enum TileShape {
    Hexagon,
    Square,
    Triangle,
    Spectre,
}

#[derive(Clone, Debug, clap::ValueEnum)]
enum RingChoice {
    ZZ4,
    ZZ12,
}

fn seed_tile_zz4(shape: &TileShape) -> Rat<ZZ4> {
    match shape {
        TileShape::Square => Rat::from_unchecked(&tiles::square::<ZZ4>()),
        _ => {
            eprintln!(
                "error: {:?} tile not supported with --ring zz4 (only square)",
                shape
            );
            std::process::exit(1);
        }
    }
}

fn seed_tile_zz12(shape: &TileShape) -> Rat<ZZ12> {
    let snake: Snake<ZZ12> = match shape {
        TileShape::Hexagon => tiles::hexagon(),
        TileShape::Square => tiles::square(),
        TileShape::Triangle => tiles::triangle(),
        TileShape::Spectre => tiles::spectre(),
    };
    Rat::from_unchecked(&snake)
}

// --- TileSet-based matching ---

fn tileset_match<T: IsComplex + IsRingOrField + Units>(
    patches: &[Rat<T>],
    seed: &[Rat<T>],
    verbose: bool,
    validate: bool,
) -> (BTreeSet<Rat<T>>, GlueStats) {
    if patches.is_empty() || seed.is_empty() {
        return (BTreeSet::new(), GlueStats::default());
    }
    let count_a = patches.len();
    let count_b = seed.len();

    let mut all_tiles: Vec<Rat<T>> = patches.to_vec();
    all_tiles.extend(seed.iter().cloned());

    let ts = TileSet::new(all_tiles);

    let pairs: Vec<(usize, usize)> = (0..count_a)
        .flat_map(|i| (count_a..count_a + count_b).map(move |j| (i, j)))
        .collect();

    let (results, stats) = if validate {
        let t_fast = Instant::now();
        let (fast, fast_stats) = ts.valid_rats_for_pairs(&pairs);
        let dt_fast = t_fast.elapsed();

        let mut bf_results = BTreeSet::new();
        let mut bf_ops = 0;
        let t_bf = Instant::now();
        for a in patches {
            for b in seed {
                for ia in 0..a.len() {
                    for ib in 0..b.len() {
                        bf_ops += 1;
                        if let Ok(glued) = a.try_glue((ia as i64, ib as i64), b) {
                            bf_results.insert(glued);
                        }
                    }
                }
            }
        }
        let dt_bf = t_bf.elapsed();

        if fast != bf_results {
            eprintln!("MISMATCH!");
            for r in fast.difference(&bf_results).take(5) {
                eprintln!("  tileset extra: {:?}", r.seq());
            }
            for r in bf_results.difference(&fast).take(5) {
                eprintln!("  brute force extra: {:?}", r.seq());
            }
            std::process::exit(1);
        }
        let speedup = dt_bf.as_secs_f64() / dt_fast.as_secs_f64().max(1e-9);
        eprintln!(
            "    validate: OK | tileset: {:.2?} ({} seqs) | brute_force: {:.2?} ({} ops) | {:.1}x speedup",
            dt_fast, fast_stats.unique_sequences, dt_bf, bf_ops, speedup,
        );
        (fast, fast_stats)
    } else {
        ts.valid_rats_for_pairs(&pairs)
    };

    if verbose {
        eprintln!(
            "    tileset: {} x {} tiles, {} distinct",
            count_a,
            count_b,
            results.len(),
        );
        eprintln!("    {stats}");
    }

    (results, stats)
}

fn run<T: IsComplex + IsRingOrField + Units>(
    seed: Rat<T>,
    max_size: usize,
    verbose: bool,
    validate: bool,
    ring_label: &str,
    tile_label: &str,
) {
    eprintln!(
        "Seed: {}, edges: {}, ring: {}, max_size: {}, validate: {}",
        tile_label,
        seed.len(),
        ring_label,
        max_size,
        validate,
    );

    let mut patches: BTreeMap<usize, Vec<Rat<T>>> = BTreeMap::new();
    let mut all_stats: BTreeMap<usize, GlueStats> = BTreeMap::new();
    let seed_vec = vec![seed];
    patches.insert(1, seed_vec.clone());

    for k in 2..=max_size {
        let t0 = Instant::now();
        let (results, stats) = tileset_match::<T>(&patches[&(k - 1)], &seed_vec, verbose, validate);
        let elapsed = t0.elapsed();

        eprintln!(
            "  size {:>3} = {:>3} + {:>3}: {} distinct patches [{:.2?}]",
            k,
            k - 1,
            1,
            results.len(),
            elapsed,
        );
        eprintln!("    {stats}");

        all_stats.insert(k, stats);
        patches.insert(k, results.into_iter().collect());
    }

    eprintln!("\n--- Summary ---");
    let mut total = 0usize;
    for (&size, pats) in &patches {
        eprintln!("  size {:>3}: {:>6} patches", size, pats.len());
        total += pats.len();
    }
    eprintln!("  total:   {total:>6} patches");

    let mut aggregate = GlueStats::default();
    for s in all_stats.values() {
        aggregate += s.clone();
    }
    eprintln!("\n--- Aggregate Stats ---");
    eprintln!("  {aggregate}");
}

fn main() {
    let args = Args::parse();

    if args.max_size < 1 {
        eprintln!("max_size must be at least 1");
        std::process::exit(1);
    }
    if args.validate && args.max_size > 4 {
        eprintln!("note: --validate with max_size > 4 may be slow (brute force scales poorly)");
    }

    let tile_label = format!("{:?}", args.tile);
    let ring_label = format!("{:?}", args.ring).to_lowercase();

    match args.ring {
        RingChoice::ZZ4 => {
            let seed = seed_tile_zz4(&args.tile);
            run::<ZZ4>(
                seed,
                args.max_size,
                args.verbose,
                args.validate,
                &ring_label,
                &tile_label,
            );
        }
        RingChoice::ZZ12 => {
            let seed = seed_tile_zz12(&args.tile);
            run::<ZZ12>(
                seed,
                args.max_size,
                args.verbose,
                args.validate,
                &ring_label,
                &tile_label,
            );
        }
    }
}
