use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

use clap::Parser;
use tilezz::cyclotomic::{IsComplex, IsRingOrField, Units, ZZ12};
use tilezz::intgeom::growing::grow_redelmeier;
use tilezz::intgeom::rat::Rat;
use tilezz::intgeom::snake::Snake;
use tilezz::intgeom::tiles;
use tilezz::intgeom::tileset::TileSet;

#[derive(Parser)]
#[command(name = "patch_bench", about = "Benchmark patch enumeration approaches")]
struct Args {
    /// Maximum patch size
    #[arg(long, default_value_t = 8)]
    max_size: usize,

    /// Approach to use
    #[arg(long, value_enum, default_value = "both")]
    approach: Approach,

    /// Tile shape
    #[arg(long, value_enum, default_value = "hexagon")]
    tile: TileShape,
}

#[derive(Clone, Debug, clap::ValueEnum)]
enum Approach {
    Old,
    New,
    Both,
}

#[derive(Clone, Debug, clap::ValueEnum)]
enum TileShape {
    Hexagon,
    Square,
    Triangle,
    Spectre,
}

fn make_seed(shape: &TileShape) -> Rat<ZZ12> {
    let snake: Snake<ZZ12> = match shape {
        TileShape::Hexagon => tiles::hexagon(),
        TileShape::Square => tiles::square(),
        TileShape::Triangle => tiles::triangle(),
        TileShape::Spectre => tiles::spectre(),
    };
    Rat::from_unchecked(&snake)
}

fn run_old<T>(seed: Rat<T>, max_size: usize) -> BTreeMap<usize, BTreeSet<Rat<T>>>
where
    T: IsComplex + IsRingOrField + Units,
{
    let mut patches: BTreeMap<usize, BTreeSet<Rat<T>>> = BTreeMap::new();
    patches.insert(1, std::iter::once(seed.clone()).collect());

    for k in 2..=max_size {
        let prev: Vec<Rat<T>> = patches[&(k - 1)].iter().cloned().collect();
        let count_a = prev.len();
        let mut all_tiles = prev;
        all_tiles.push(seed.clone());
        let ts = TileSet::new(all_tiles);
        let pairs: Vec<(usize, usize)> = (0..count_a).map(|i| (i, count_a)).collect();
        let t = Instant::now();
        let (results, _) = ts.valid_rats_for_pairs(&pairs);
        let dt = t.elapsed();
        eprintln!(
            "  [old] size {:>3}: {:>6} patches [{:.2?}]",
            k,
            results.len(),
            dt,
        );
        patches.insert(k, results);
    }

    patches
}

fn run_new<T>(seed: Rat<T>, max_size: usize) -> BTreeMap<usize, BTreeSet<Rat<T>>>
where
    T: IsComplex + IsRingOrField + Units,
{
    grow_redelmeier(&seed, max_size)
}

fn print_results(
    label: &str,
    results: &BTreeMap<usize, BTreeSet<Rat<ZZ12>>>,
    elapsed: std::time::Duration,
) {
    println!("=== {label} ===");
    let mut total = 0;
    for (&size, patches) in results {
        println!("  size {size:>3}: {:>6} patches", patches.len());
        total += patches.len();
    }
    println!("  total:   {total:>6} patches");
    println!("  time:    {elapsed:.2?}");
}

fn main() {
    let args = Args::parse();

    if args.max_size < 1 {
        eprintln!("max_size must be at least 1");
        std::process::exit(1);
    }

    let seed = make_seed(&args.tile);
    let tile_label = format!("{:?}", args.tile);

    match args.approach {
        Approach::Old => {
            let t0 = Instant::now();
            let results = run_old(seed, args.max_size);
            let elapsed = t0.elapsed();
            print_results(&format!("Old (TileSet) - {tile_label}"), &results, elapsed);
        }
        Approach::New => {
            let t0 = Instant::now();
            let results = run_new(seed, args.max_size);
            let elapsed = t0.elapsed();
            print_results(
                &format!("New (Redelmeier) - {tile_label}"),
                &results,
                elapsed,
            );
        }
        Approach::Both => {
            let t0 = Instant::now();
            let old_results = run_old(seed.clone(), args.max_size);
            let old_elapsed = t0.elapsed();

            let t0 = Instant::now();
            let new_results = run_new(seed, args.max_size);
            let new_elapsed = t0.elapsed();

            println!();
            print_results(
                &format!("Old (TileSet) - {tile_label}"),
                &old_results,
                old_elapsed,
            );
            println!();
            print_results(
                &format!("New (Redelmeier) - {tile_label}"),
                &new_results,
                new_elapsed,
            );

            println!("\n=== Comparison ===");
            let matches = old_results == new_results;
            println!("Results match: {}", if matches { "YES" } else { "NO" });
            if !matches {
                for k in 1..=args.max_size {
                    let old_count = old_results.get(&k).map(|s| s.len()).unwrap_or(0);
                    let new_count = new_results.get(&k).map(|s| s.len()).unwrap_or(0);
                    if old_count != new_count {
                        eprintln!("  size {k}: old={old_count} new={new_count}");
                    }
                }
            }
            let speedup = old_elapsed.as_secs_f64() / new_elapsed.as_secs_f64().max(1e-9);
            println!("Speedup: {speedup:.2}x");
        }
    }
}
