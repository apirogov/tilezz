use std::collections::BTreeMap;
use std::time::Instant;

use clap::Parser;
use rustc_hash::FxHashSet;
use tilezz::cyclotomic::ZZ12;
use tilezz::intgeom::growing::{grow_redelmeier, grow_redelmeier_free};
use tilezz::intgeom::rat::Rat;
use tilezz::intgeom::snake::Snake;
use tilezz::intgeom::tiles;

#[cfg(feature = "pprof")]
use pprof::ProfilerGuardBuilder;

#[derive(Parser)]
#[command(
    name = "patch_bench",
    about = "Benchmark Redelmeier patch enumeration - ZZ12 only"
)]
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

fn run_redelmeier(
    seed: Rat<ZZ12>,
    max_size: usize,
    free: bool,
) -> BTreeMap<usize, FxHashSet<Rat<ZZ12>>> {
    if free {
        grow_redelmeier_free(&seed, max_size)
    } else {
        grow_redelmeier(&seed, max_size)
    }
}

fn main() {
    let args = Args::parse();

    #[cfg(feature = "pprof")]
    let guard = ProfilerGuardBuilder::default()
        .frequency(1000)
        .blocklist(&["libc", "libgcc", "pthread", "vdso"])
        .build()
        .unwrap();

    let seed = make_seed(&args.tile);

    let full_label = format!("{:?} [ZZ12]", args.tile);
    let label_free = if args.free { " (free)" } else { "" };

    let t0 = Instant::now();
    let final_results = run_redelmeier(seed, args.max_size, args.free);
    let elapsed = t0.elapsed();

    println!("=== Redelmeier - {}{} ===", full_label, label_free);
    let mut total = 0usize;
    for k in 1..=args.max_size {
        let count = final_results.get(&k).map(|s| s.len()).unwrap_or(0);
        println!("  size {:>3}: {:>6} patches", k, count);
        total += count;
    }
    println!("  total: {:>6} patches", total);
    println!("  time:  {:.2?}", elapsed);

    #[cfg(feature = "pprof")]
    {
        if let Ok(report) = guard.report().build() {
            let flamegraph_file = "flamegraph.svg";
            let mut file = std::fs::File::create(flamegraph_file).unwrap();
            report.flamegraph(&mut file).unwrap();
            eprintln!("\nFlame graph written to {}", flamegraph_file);
        }
    }
}
