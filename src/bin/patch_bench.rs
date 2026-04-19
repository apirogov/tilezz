use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

use clap::Parser;
use tilezz::cyclotomic::{IsComplex, IsRingOrField, Units, ZZ12};
use tilezz::intgeom::growing::{grow_redelmeier, make_free};
use tilezz::intgeom::rat::Rat;
use tilezz::intgeom::snake::Snake;
use tilezz::intgeom::tiles;

#[derive(Parser)]
#[command(name = "patch_bench", about = "Benchmark Redelmeier patch enumeration - ZZ12 only")]
struct Args {
    #[arg(long, default_value_t = 9)]
    max_size: usize,

    #[arg(long, value_enum, default_value = "hexagon")]
    tile: TileShape,

    #[arg(long, default_value_t = false)]
    free: bool,
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

fn run_redelmeier(seed: Rat<ZZ12>, max_size: usize) -> BTreeMap<usize, BTreeSet<Rat<ZZ12>>> {
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

    let t0 = Instant::now();
    let results = run_redelmeier(seed, args.max_size);
    let elapsed = t0.elapsed();

    let final_results = if args.free { to_free(&results) } else { results };

    println!("=== Redelmeier - {}{} ===", full_label, label_free);
    let mut total = 0usize;
    for k in 1..=args.max_size {
        let count = final_results.get(&k).map(|s| s.len()).unwrap_or(0);
        println!("  size {:>3}: {:>6} patches", k, count);
        total += count;
    }
    println!("  total: {:>6} patches", total);
    println!("  time:  {:.2?}", elapsed);
}