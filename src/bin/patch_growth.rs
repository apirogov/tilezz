use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

use clap::Parser;

use tilezz::cyclotomic::{IsComplex, IsRingOrField, Units, ZZ12};
use tilezz::intgeom::rat::Rat;
use tilezz::intgeom::snake::Snake;
use tilezz::intgeom::tiles;
use tilezz::intgeom::tileset::TileSet;

#[derive(Parser)]
#[command(
    name = "patch_growth",
    about = "Enumerate tile patches via binary decomposition"
)]
struct Args {
    /// Seed tile shape
    #[arg(long, value_enum, default_value = "hexagon")]
    tile: TileShape,

    /// Maximum patch size to compute
    #[arg(long, default_value_t = 6)]
    max_size: usize,

    /// Print per-pair glue statistics
    #[arg(long)]
    verbose: bool,
}

#[derive(Clone, Debug, clap::ValueEnum)]
enum TileShape {
    Hexagon,
    Square,
    Triangle,
    Spectre,
}

fn seed_tile(shape: &TileShape) -> Rat<ZZ12> {
    let snake: Snake<ZZ12> = match shape {
        TileShape::Hexagon => tiles::hexagon(),
        TileShape::Square => tiles::square(),
        TileShape::Triangle => tiles::triangle(),
        TileShape::Spectre => tiles::spectre(),
    };
    Rat::from_unchecked(&snake)
}

fn largest_power_of_2_leq(n: usize) -> usize {
    assert!(n > 0);
    1 << (usize::BITS - 1 - n.leading_zeros())
}

fn is_power_of_2(n: usize) -> bool {
    n > 0 && (n & (n - 1)) == 0
}

fn self_match<T: IsComplex + IsRingOrField + Units>(
    patches: &[Rat<T>],
    verbose: bool,
) -> (BTreeSet<Rat<T>>, usize, usize) {
    if patches.is_empty() {
        return (BTreeSet::new(), 0, 0);
    }
    let ts = TileSet::new(patches.to_vec());
    let glues = ts.all_valid_glues();
    let attempts = glues.len();
    let mut results = BTreeSet::new();
    for g in &glues {
        results.insert(g.result.clone());
    }
    if verbose {
        let mut by_len = BTreeMap::new();
        for g in &glues {
            by_len.entry(g.match_len).or_insert_with(Vec::new).push(g);
        }
        for (len, gs) in &by_len {
            eprintln!("    match_len={}: {} glues", len, gs.len(),);
        }
    }
    (results, attempts, ts.num_tiles())
}

fn cross_match<T: IsComplex + IsRingOrField + Units>(
    patches_a: &[Rat<T>],
    patches_b: &[Rat<T>],
    verbose: bool,
) -> (BTreeSet<Rat<T>>, usize) {
    if patches_a.is_empty() || patches_b.is_empty() {
        return (BTreeSet::new(), 0);
    }

    let count_a = patches_a.len();
    let count_b = patches_b.len();

    let mut all_tiles: Vec<Rat<T>> = patches_a.to_vec();
    all_tiles.extend(patches_b.iter().cloned());

    let ts = TileSet::new(all_tiles);

    let mut results = BTreeSet::new();
    let mut attempts = 0;

    for i in 0..count_a {
        for j in count_a..(count_a + count_b) {
            let glues = ts.valid_glues(i, j);
            attempts += glues.len();
            for g in &glues {
                results.insert(g.result.clone());
            }
        }
    }

    if verbose {
        eprintln!(
            "    cross_match: {} x {} tiles, {} glue ops, {} distinct results",
            count_a,
            count_b,
            attempts,
            results.len(),
        );
    }

    (results, attempts)
}

fn main() {
    let args = Args::parse();

    if args.max_size < 1 {
        eprintln!("max_size must be at least 1");
        std::process::exit(1);
    }

    let seed = seed_tile(&args.tile);
    eprintln!(
        "Seed: {:?}, edges: {}, max_size: {}",
        args.tile,
        seed.len(),
        args.max_size,
    );

    let mut patches: BTreeMap<usize, Vec<Rat<ZZ12>>> = BTreeMap::new();
    patches.insert(1, vec![seed]);

    eprintln!("\n--- Phase 1: Powers of 2 ---");
    let max_power = largest_power_of_2_leq(args.max_size);
    let mut power = 1usize;
    while power < max_power {
        let next = power * 2;
        let t0 = Instant::now();
        let (results, attempts, num_tiles) = self_match(&patches[&power], args.verbose);
        let elapsed = t0.elapsed();
        eprintln!(
            "  size {:>3} = {:>3} + {:>3}: {} distinct patches ({} glue ops on {} tiles) [{:.2?}]",
            next,
            power,
            power,
            results.len(),
            attempts,
            num_tiles,
            elapsed,
        );
        patches.insert(next, results.into_iter().collect());
        power = next;
    }

    if args.max_size > 1 {
        eprintln!("\n--- Phase 2: Intermediate sizes ---");
    }
    for k in 2..=args.max_size {
        if is_power_of_2(k) {
            continue;
        }
        let a = largest_power_of_2_leq(k);
        let b = k - a;
        assert!(
            patches.contains_key(&a),
            "size {a} not yet computed (needed for size {k})",
        );
        assert!(
            patches.contains_key(&b),
            "size {b} not yet computed (needed for size {k})",
        );

        let t0 = Instant::now();
        let (results, attempts) = cross_match(&patches[&a], &patches[&b], args.verbose);
        let elapsed = t0.elapsed();
        eprintln!(
            "  size {:>3} = {:>3} + {:>3}: {} distinct patches ({} glue ops) [{:.2?}]",
            k,
            a,
            b,
            results.len(),
            attempts,
            elapsed,
        );
        patches.insert(k, results.into_iter().collect());
    }

    eprintln!("\n--- Summary ---");
    let mut total = 0usize;
    for (&size, pats) in &patches {
        eprintln!("  size {:>3}: {:>6} patches", size, pats.len());
        total += pats.len();
    }
    eprintln!("  total:   {total:>6} patches");
}
