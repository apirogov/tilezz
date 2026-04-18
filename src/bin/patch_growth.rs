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
    about = "Enumerate tile patches via binary decomposition"
)]
struct Args {
    /// Seed tile shape
    #[arg(long, value_enum, default_value = "hexagon")]
    tile: TileShape,

    /// Maximum patch size to compute
    #[arg(long, default_value_t = 6)]
    max_size: usize,

    /// Cyclotomic ring (ZZ4 only supports square tiles)
    #[arg(long, value_enum, default_value = "zz12")]
    ring: RingChoice,

    /// Growth strategy: naive (one-by-one), binary, or both (compare)
    #[arg(long, value_enum, default_value = "naive")]
    strategy: Strategy,

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

#[derive(Clone, Debug, clap::ValueEnum)]
enum Strategy {
    Naive,
    Binary,
    Both,
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

// --- TileSet-based matching (CMI + heuristic) ---

fn tileset_match<T: IsComplex + IsRingOrField + Units>(
    patches_a: &[Rat<T>],
    patches_b: &[Rat<T>],
    verbose: bool,
) -> (BTreeSet<Rat<T>>, GlueStats) {
    if patches_a.is_empty() || patches_b.is_empty() {
        return (BTreeSet::new(), GlueStats::default());
    }
    let count_a = patches_a.len();
    let count_b = patches_b.len();

    let mut all_tiles: Vec<Rat<T>> = patches_a.to_vec();
    all_tiles.extend(patches_b.iter().cloned());

    let ts = TileSet::new(all_tiles);

    let self_match = std::ptr::eq(patches_a.as_ptr(), patches_b.as_ptr());
    let pairs: Vec<(usize, usize)> = if self_match {
        (0..count_a)
            .flat_map(|i| (0..count_a).map(move |j| (i, j)))
            .collect()
    } else {
        (0..count_a)
            .flat_map(|i| (count_a..count_a + count_b).map(move |j| (i, j)))
            .collect()
    };

    let (results, stats) = ts.valid_rats_for_pairs(&pairs);

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

// --- Brute-force matching (raw try_glue, no CMI) ---

fn brute_force_match<T: IsComplex + IsRingOrField + Units>(
    patches_a: &[Rat<T>],
    patches_b: &[Rat<T>],
) -> (BTreeSet<Rat<T>>, usize) {
    if patches_a.is_empty() || patches_b.is_empty() {
        return (BTreeSet::new(), 0);
    }

    let mut results = BTreeSet::new();
    let mut attempts = 0;

    for a in patches_a {
        for b in patches_b {
            for ia in 0..a.len() {
                for ib in 0..b.len() {
                    attempts += 1;
                    if let Ok(glued) = a.try_glue((ia as i64, ib as i64), b) {
                        results.insert(glued);
                    }
                }
            }
        }
    }

    (results, attempts)
}

// --- Growth strategies ---

struct RoundResult<T: IsComplex> {
    patches: Vec<Rat<T>>,
    elapsed: std::time::Duration,
    num_tiles_a: usize,
    num_tiles_b: usize,
    stats: GlueStats,
}

fn compute_round<T: IsComplex + IsRingOrField + Units>(
    pa: &[Rat<T>],
    pb: &[Rat<T>],
    verbose: bool,
    mode: &Mode,
) -> RoundResult<T> {
    let t0 = Instant::now();
    let (results, stats) = match mode {
        Mode::Fast => tileset_match(pa, pb, verbose),
        Mode::Validate => {
            let t_fast = Instant::now();
            let (fast, fast_stats) = tileset_match(pa, pb, verbose);
            let dt_fast = t_fast.elapsed();

            let t_bf = Instant::now();
            let (bf, bf_ops) = brute_force_match(pa, pb);
            let dt_bf = t_bf.elapsed();

            if fast != bf {
                let fast_only: Vec<_> = fast.difference(&bf).collect();
                let bf_only: Vec<_> = bf.difference(&fast).collect();
                eprintln!("MISMATCH!");
                if !fast_only.is_empty() {
                    eprintln!(
                        "  tileset found {} extra patches (false positives):",
                        fast_only.len()
                    );
                    for r in fast_only.iter().take(5) {
                        eprintln!("    {:?}", r.seq());
                    }
                }
                if !bf_only.is_empty() {
                    eprintln!(
                        "  brute force found {} extra patches (missed by tileset):",
                        bf_only.len()
                    );
                    for r in bf_only.iter().take(5) {
                        eprintln!("    {:?}", r.seq());
                    }
                }
                std::process::exit(1);
            }
            let speedup = dt_bf.as_secs_f64() / dt_fast.as_secs_f64().max(1e-9);
            eprintln!(
                "    validate: OK | tileset: {:.2?} ({} seqs) | brute_force: {:.2?} ({} ops) | {:.1}x speedup | {} patches",
                dt_fast, fast_stats.unique_sequences, dt_bf, bf_ops, speedup, fast.len(),
            );
            (fast, fast_stats)
        }
    };
    RoundResult {
        patches: results.into_iter().collect(),
        elapsed: t0.elapsed(),
        num_tiles_a: pa.len(),
        num_tiles_b: pb.len(),
        stats,
    }
}

enum Mode {
    Fast,
    Validate,
}

fn largest_power_of_2_leq(n: usize) -> usize {
    assert!(n > 0);
    1 << (usize::BITS - 1 - n.leading_zeros())
}

fn is_power_of_2(n: usize) -> bool {
    n > 0 && (n & (n - 1)) == 0
}

fn print_summary<T: IsComplex>(
    patches: &BTreeMap<usize, Vec<Rat<T>>>,
    all_stats: &BTreeMap<usize, GlueStats>,
) {
    eprintln!("\n--- Summary ---");
    let mut total = 0usize;
    for (&size, pats) in patches {
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

fn run_naive<T: IsComplex + IsRingOrField + Units>(
    seed: Rat<T>,
    max_size: usize,
    verbose: bool,
    mode: &Mode,
) -> BTreeMap<usize, Vec<Rat<T>>> {
    let mut patches: BTreeMap<usize, Vec<Rat<T>>> = BTreeMap::new();
    let mut all_stats: BTreeMap<usize, GlueStats> = BTreeMap::new();
    let seed_vec = vec![seed];
    patches.insert(1, seed_vec.clone());

    for k in 2..=max_size {
        let r = compute_round(&patches[&(k - 1)], &seed_vec, verbose, mode);
        eprintln!(
            "  size {:>3} = {:>3} + {:>3}: {} distinct patches ({} x {} tiles) [{:.2?}]",
            k,
            k - 1,
            1,
            r.patches.len(),
            r.num_tiles_a,
            r.num_tiles_b,
            r.elapsed,
        );
        eprintln!("    {}", r.stats);
        all_stats.insert(k, r.stats);
        patches.insert(k, r.patches);
    }

    print_summary(&patches, &all_stats);
    patches
}

fn run_binary<T: IsComplex + IsRingOrField + Units>(
    seed: Rat<T>,
    max_size: usize,
    verbose: bool,
    mode: &Mode,
) -> BTreeMap<usize, Vec<Rat<T>>> {
    let mut patches: BTreeMap<usize, Vec<Rat<T>>> = BTreeMap::new();
    let mut all_stats: BTreeMap<usize, GlueStats> = BTreeMap::new();
    patches.insert(1, vec![seed]);

    eprintln!("\n--- Phase 1: Powers of 2 ---");
    let max_power = largest_power_of_2_leq(max_size);
    let mut power = 1usize;
    while power < max_power {
        let next = power * 2;
        let r = compute_round(&patches[&power], &patches[&power], verbose, mode);
        eprintln!(
            "  size {:>3} = {:>3} + {:>3}: {} distinct patches ({} x {} tiles) [{:.2?}]",
            next,
            power,
            power,
            r.patches.len(),
            r.num_tiles_a,
            r.num_tiles_b,
            r.elapsed,
        );
        eprintln!("    {}", r.stats);
        all_stats.insert(next, r.stats);
        patches.insert(next, r.patches);
        power = next;
    }

    if max_size > 1 {
        eprintln!("\n--- Phase 2: Intermediate sizes ---");
    }
    for k in 2..=max_size {
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

        let r = compute_round(&patches[&a], &patches[&b], verbose, mode);
        eprintln!(
            "  size {:>3} = {:>3} + {:>3}: {} distinct patches ({} x {} tiles) [{:.2?}]",
            k,
            a,
            b,
            r.patches.len(),
            r.num_tiles_a,
            r.num_tiles_b,
            r.elapsed,
        );
        eprintln!("    {}", r.stats);
        all_stats.insert(k, r.stats);
        patches.insert(k, r.patches);
    }

    print_summary(&patches, &all_stats);
    patches
}

fn run<T: IsComplex + IsRingOrField + Units>(
    seed: Rat<T>,
    args: &Args,
    mode: &Mode,
    ring_label: &str,
    tile_label: &str,
) {
    eprintln!(
        "Seed: {}, edges: {}, ring: {}, max_size: {}, strategy: {}, mode: {}",
        tile_label,
        seed.len(),
        ring_label,
        args.max_size,
        match args.strategy {
            Strategy::Naive => "naive",
            Strategy::Binary => "binary",
            Strategy::Both => "both",
        },
        match mode {
            Mode::Fast => "fast",
            Mode::Validate => "validate",
        },
    );

    match args.strategy {
        Strategy::Naive => {
            run_naive::<T>(seed, args.max_size, args.verbose, mode);
        }
        Strategy::Binary => {
            run_binary::<T>(seed, args.max_size, args.verbose, mode);
        }
        Strategy::Both => {
            eprintln!("\n=== Binary decomposition ===");
            let binary = run_binary::<T>(seed.clone(), args.max_size, args.verbose, mode);
            eprintln!("\n=== Naive (one-by-one) ===");
            let naive = run_naive::<T>(seed, args.max_size, args.verbose, mode);

            eprintln!("\n=== Comparison ===");
            for k in 2..=args.max_size {
                let nb = binary.get(&k).map(|v| v.len()).unwrap_or(0);
                let nn = naive.get(&k).map(|v| v.len()).unwrap_or(0);
                let missed = nn.saturating_sub(nb);
                let extra = nb.saturating_sub(nn);
                let tag = if missed > 0 { " <<< MISSED" } else { "" };
                eprintln!(
                    "  size {:>3}: binary={:>6}  naive={:>6}  missed={:>5}  extra={:>5}{}",
                    k, nb, nn, missed, extra, tag,
                );
            }
        }
    }
}

fn main() {
    let args = Args::parse();

    let max_size = args.max_size;
    if args.validate && max_size > 4 {
        eprintln!("note: --validate with max_size > 4 may be slow (brute force scales poorly)");
    }
    if max_size < 1 {
        eprintln!("max_size must be at least 1");
        std::process::exit(1);
    }

    let mode = if args.validate {
        Mode::Validate
    } else {
        Mode::Fast
    };

    let tile_label = format!("{:?}", args.tile);
    let ring_label = format!("{:?}", args.ring).to_lowercase();

    match args.ring {
        RingChoice::ZZ4 => {
            let seed = seed_tile_zz4(&args.tile);
            run::<ZZ4>(seed, &args, &mode, &ring_label, &tile_label);
        }
        RingChoice::ZZ12 => {
            let seed = seed_tile_zz12(&args.tile);
            run::<ZZ12>(seed, &args, &mode, &ring_label, &tile_label);
        }
    }
}
