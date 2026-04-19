use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

use clap::Parser;
use tilezz::cyclotomic::{ZZ12, IsComplex, IsRingOrField, Units};
use tilezz::intgeom::growing::{grow_redelmeier, make_free};
use tilezz::intgeom::rat::Rat;
use tilezz::intgeom::snake::Snake;
use tilezz::intgeom::tiles;
use tilezz::intgeom::tileset::TileSet;

#[derive(Parser)]
#[command(name = "patch_bench", about = "Benchmark patch enumeration approaches - ZZ12 only")]
struct Args {
    /// Maximum patch size
    #[arg(long, default_value_t = 7)]
    max_size: usize,

    /// Approach to use
    #[arg(long, value_enum, default_value = "both")]
    approach: Approach,

    /// Tile shape
    #[arg(long, value_enum, default_value = "hexagon")]
    tile: TileShape,

    /// Mod out reflections (report free instead of one-sided counts)
    #[arg(long, default_value_t = false)]
    free: bool,
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

fn run_old(seed: Rat<ZZ12>, max_size: usize) -> BTreeMap<usize, BTreeSet<Rat<ZZ12>>> {
    let mut patches: BTreeMap<usize, BTreeSet<Rat<ZZ12>>> = BTreeMap::new();
    patches.insert(1, std::iter::once(seed.clone()).collect());

    for k in 2..=max_size {
        let prev: Vec<Rat<ZZ12>> = patches[&(k - 1)].iter().cloned().collect();
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

fn run_new(seed: Rat<ZZ12>, max_size: usize) -> BTreeMap<usize, BTreeSet<Rat<ZZ12>>> {
    grow_redelmeier(&seed, max_size)
}

fn to_free<T: IsComplex + IsRingOrField + Units>(
    onesided: &BTreeMap<usize, BTreeSet<Rat<T>>>,
) -> BTreeMap<usize, BTreeSet<Rat<T>>>
where
    T: Clone,
{
    onesided
        .iter()
        .map(|(&k, set)| (k, make_free(set)))
        .collect()
}

fn main() {
    let args = Args::parse();

    let seed = make_seed(&args.tile);

    let full_label = format!("{:?} [ZZ12]", args.tile);
    let label_free = if args.free { " (free)" } else { "" };

    match args.approach {
        Approach::Old => {
            let results = run_old(seed, args.max_size);
            let results = if args.free { to_free(&results) } else { results };
            println!("=== Old - {}{} ===", full_label, label_free);
            let mut total = 0usize;
            for k in 1..=args.max_size {
                let count = results.get(&k).map(|s| s.len()).unwrap_or(0);
                println!("  size {:>3}: {:>6} patches", k, count);
                total += count;
            }
            println!("  total: {:>6} patches", total);
        }
        Approach::New => {
            let results = run_new(seed, args.max_size);
            let results = if args.free { to_free(&results) } else { results };
            println!("=== New (Redelmeier) - {}{} ===", full_label, label_free);
            let mut total = 0usize;
            for k in 1..=args.max_size {
                let count = results.get(&k).map(|s| s.len()).unwrap_or(0);
                println!("  size {:>3}: {:>6} patches", k, count);
                total += count;
            }
            println!("  total: {:>6} patches", total);
        }
        Approach::Both => {
            let t_old = Instant::now();
            let old_results = run_old(seed.clone(), args.max_size);
            let old_time = t_old.elapsed();

            let t_new = Instant::now();
            let new_results = run_new(seed, args.max_size);
            let new_time = t_new.elapsed();

            let old_final = if args.free { to_free(&old_results) } else { old_results };
            let new_final = if args.free { to_free(&new_results) } else { new_results };

            println!("=== Old - {}{} ===", full_label, label_free);
            let mut old_total = 0usize;
            for k in 1..=args.max_size {
                let count = old_final.get(&k).map(|s| s.len()).unwrap_or(0);
                println!("  size {:>3}: {:>6} patches", k, count);
                old_total += count;
            }
            println!("  total: {:>6} patches", old_total);
            println!("  time:  {:.2?}", old_time);

            println!();
            println!("=== New (Redelmeier) - {}{} ===", full_label, label_free);
            let mut new_total = 0usize;
            for k in 1..=args.max_size {
                let count = new_final.get(&k).map(|s| s.len()).unwrap_or(0);
                println!("  size {:>3}: {:>6} patches", k, count);
                new_total += count;
            }
            println!("  total: {:>6} patches", new_total);
            println!("  time:  {:.2?}", new_time);

            println!();
            println!("=== Comparison ===");
            let _match_old = old_final.len();
            let _match_new = new_final.len();
            let results_match = old_total == new_total;
            println!("Results match: {}", if results_match { "YES" } else { "NO" });

            let speedup = old_time.as_secs_f64() / new_time.as_secs_f64();
            println!("Speedup: {:.2}x", speedup);
        }
    }
}
