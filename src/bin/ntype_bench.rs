use std::sync::Arc;
use std::time::Instant;

use clap::Parser;
use tilezz::cyclotomic::{IsComplex, IsRingOrField, Units, ZZ10, ZZ12};
use tilezz::intgeom::neighborhood::{NeighborhoodIndex, NtEntry, NtKind};
use tilezz::intgeom::tileset::{self, TileSet};

#[cfg(feature = "pprof")]
use pprof::ProfilerGuardBuilder;

#[derive(Parser)]
#[command(
    name = "ntype_bench",
    about = "Benchmark neighborhood type enumeration"
)]
struct Args {
    #[arg(long, value_enum, default_value = "hex")]
    tile: TileSetKind,
}

#[derive(Clone, Debug, clap::ValueEnum)]
enum TileSetKind {
    Square,
    Hex,
    Mixed,
    Tetris,
    Spectre,
    Penrose,
}

fn run_bench<T: IsComplex + IsRingOrField + Units>(label: &str, ts: Arc<TileSet<T>>) {
    eprintln!("=== Neighborhood types: {} ===", label);

    #[cfg(feature = "pprof")]
    let guard = ProfilerGuardBuilder::default()
        .frequency(1000)
        .blocklist(&["libc", "libgcc", "pthread", "vdso"])
        .build()
        .unwrap();

    let t0 = Instant::now();
    let idx = NeighborhoodIndex::new(Arc::clone(&ts));
    let elapsed = t0.elapsed();

    let t1 = Instant::now();
    let kinds = idx.classify_all();
    let classify_time = t1.elapsed();

    let transitions = idx.transitions();
    // Closed phase-2 SurroundedTile entries are the terminal Blessed
    // states (= the "closed corona" leaves). No more synthetic
    // closing transitions in the BFS.
    let closed_phase2 = idx
        .entries()
        .iter()
        .filter(|e| matches!(e, NtEntry::Phase2(st) if st.is_closed))
        .count();

    let dead = kinds.iter().filter(|&&k| k == NtKind::Dead).count();
    let undead = kinds.iter().filter(|&&k| k == NtKind::Undead).count();
    let blessed = kinds.iter().filter(|&&k| k == NtKind::Blessed).count();
    let free = kinds.iter().filter(|&&k| k == NtKind::Free).count();

    println!("=== {} ===", label);
    println!("  total: {} types", idx.num_types());
    println!("  time:  {:.2?} (classify: {:.2?})", elapsed, classify_time);
    println!();
    println!("  Classification:");
    println!(
        "    dead={} undead={} blessed={} free={}",
        dead, undead, blessed, free
    );
    println!();
    println!(
        "  Transitions: {} (closed phase-2 entries: {})",
        transitions.len(),
        closed_phase2
    );

    #[cfg(feature = "pprof")]
    {
        if let Ok(report) = guard.report().build() {
            let flamegraph_file = "flamegraph_ntype.svg";
            let mut file = std::fs::File::create(flamegraph_file).unwrap();
            report.flamegraph(&mut file).unwrap();
            eprintln!("\nFlame graph written to {}", flamegraph_file);
        }
    }
}

fn main() {
    let args = Args::parse();

    match args.tile {
        TileSetKind::Square => run_bench("Square", tileset::square::<ZZ12>()),
        TileSetKind::Hex => run_bench("Hex", tileset::hex::<ZZ12>()),
        TileSetKind::Mixed => run_bench("Mixed", tileset::mixed::<ZZ12>()),
        TileSetKind::Tetris => run_bench("Tetris", tileset::tetrominoes::<ZZ12>()),
        TileSetKind::Spectre => run_bench("Spectre", tileset::spectre::<ZZ12>()),
        TileSetKind::Penrose => run_bench("Penrose", tileset::penrose::<ZZ10>()),
    }
}
