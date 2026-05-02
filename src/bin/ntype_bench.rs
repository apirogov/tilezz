use std::collections::BTreeMap;
use std::sync::Arc;
use std::time::Instant;

use clap::Parser;
use tilezz::cyclotomic::{IsComplex, IsRingOrField, Units, ZZ10, ZZ12};
use tilezz::intgeom::neighborhood::OpenNeighborhoodIndex;
use tilezz::intgeom::rat::Rat;
use tilezz::intgeom::snake::Snake;
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
    let idx = OpenNeighborhoodIndex::new(Arc::clone(&ts));
    let elapsed = t0.elapsed();

    let mut by_gap: BTreeMap<usize, usize> = BTreeMap::new();
    for nhood in idx.entries() {
        *by_gap.entry(nhood.gap_len()).or_insert(0) += 1;
    }

    println!("=== {} ===", label);
    for (gl, cnt) in &by_gap {
        println!("  gap_len={}: {} types", gl, cnt);
    }
    println!("  total: {} types", idx.num_types());
    println!("  time:  {:.2?}", elapsed);

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
