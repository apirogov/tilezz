use std::collections::BTreeMap;
use std::sync::Arc;
use std::time::Instant;

use clap::Parser;
use tilezz::cyclotomic::{IsComplex, IsRingOrField, Units, ZZ10, ZZ12};
use tilezz::intgeom::neighborhood::{NeighborhoodIndex, NtKind, NT_CLOSED_ID};
use tilezz::intgeom::rat::Rat;

use tilezz::intgeom::tiles;
use tilezz::intgeom::tileset::TileSet;

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
    let closed_transitions = transitions
        .iter()
        .filter(|t| t.dst_id == NT_CLOSED_ID)
        .count();
    let close_branch = transitions.len() - closed_transitions;

    let mut by_gap: BTreeMap<usize, usize> = BTreeMap::new();
    let mut by_gap_kind: BTreeMap<(usize, NtKind), usize> = BTreeMap::new();
    for (i, nhood) in idx.entries().iter().enumerate() {
        *by_gap.entry(nhood.gap_len()).or_insert(0) += 1;
        *by_gap_kind.entry((nhood.gap_len(), kinds[i])).or_insert(0) += 1;
    }

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
        "  Transitions: {} (close-branch={}, CLOSED={})",
        transitions.len(),
        close_branch,
        closed_transitions
    );
    println!();
    println!("  By gap_len:");
    for (gl, cnt) in &by_gap {
        let d = by_gap_kind.get(&(*gl, NtKind::Dead)).copied().unwrap_or(0);
        let u = by_gap_kind
            .get(&(*gl, NtKind::Undead))
            .copied()
            .unwrap_or(0);
        let b = by_gap_kind
            .get(&(*gl, NtKind::Blessed))
            .copied()
            .unwrap_or(0);
        let f = by_gap_kind.get(&(*gl, NtKind::Free)).copied().unwrap_or(0);
        println!(
            "    gap_len={}: {} (d={} u={} b={} f={})",
            gl, cnt, d, u, b, f
        );
    }

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
        TileSetKind::Square => {
            let rat = Rat::try_from(&tiles::square::<ZZ12>()).unwrap();
            run_bench("Square", Arc::new(TileSet::new(vec![rat])));
        }
        TileSetKind::Hex => {
            let rat = Rat::try_from(&tiles::hexagon::<ZZ12>()).unwrap();
            run_bench("Hex", Arc::new(TileSet::new(vec![rat])));
        }
        TileSetKind::Mixed => {
            let sq = Rat::try_from(&tiles::square::<ZZ12>()).unwrap();
            let hex = Rat::try_from(&tiles::hexagon::<ZZ12>()).unwrap();
            run_bench("Mixed", Arc::new(TileSet::new(vec![sq, hex])));
        }
        TileSetKind::Spectre => {
            let rat = Rat::try_from(&tiles::spectre::<ZZ12>()).unwrap();
            run_bench("Spectre", Arc::new(TileSet::new(vec![rat])));
        }
        TileSetKind::Penrose => {
            let n = Rat::try_from(&tiles::penrose_p3_narrow::<ZZ10>()).unwrap();
            let w = Rat::try_from(&tiles::penrose_p3_wide::<ZZ10>()).unwrap();
            run_bench("Penrose", Arc::new(TileSet::new(vec![n, w])));
        }
    }
}
