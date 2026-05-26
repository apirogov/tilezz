use std::sync::Arc;
use std::time::Instant;

use clap::{Parser, Subcommand};
use tilezz::cyclotomic::{IsRingOrField, ZZ10, ZZ12};
use tilezz::intgeom::rat::Rat;
use tilezz::intgeom::tileset::{self, TileSet};
use tilezz::intgeom::vertextypes::{Collection, OpenVertexTypeIndex};

#[cfg(feature = "pprof")]
use pprof::ProfilerGuardBuilder;

#[derive(Parser)]
#[command(
    name = "vtype_enum",
    about = "Vertex type enumeration: collect + validate"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    Collect {
        #[arg(long, value_enum, default_value = "hex")]
        tile: TileSetKind,
        #[arg(long)]
        output: String,
    },
    Validate {
        #[arg(long)]
        input: String,
    },
}

#[derive(Clone, Debug, clap::ValueEnum)]
enum TileSetKind {
    Hex,
    Square,
    Mixed,
    Tetris,
    Spectre,
    Penrose,
}

fn make_ts_12(kind: &TileSetKind) -> Arc<TileSet<ZZ12>> {
    match kind {
        TileSetKind::Hex => tileset::hex::<ZZ12>(),
        TileSetKind::Square => tileset::square::<ZZ12>(),
        TileSetKind::Spectre => tileset::spectre::<ZZ12>(),
        TileSetKind::Mixed => tileset::mixed::<ZZ12>(),
        TileSetKind::Tetris => tileset::tetrominoes::<ZZ12>(),
        TileSetKind::Penrose => panic!("penrose requires ZZ10"),
    }
}

fn make_ts_10() -> Arc<TileSet<ZZ10>> {
    tileset::penrose::<ZZ10>()
}

fn collect_generic<T: IsRingOrField>(
    ts: Arc<TileSet<T>>,
    ring_name: &str,
    output: &str,
    label: &str,
) {
    eprintln!("[{}] Running BFS...", label);
    let t0 = Instant::now();
    let idx = OpenVertexTypeIndex::new(ts);
    let elapsed = t0.elapsed();

    let n_total = idx.num_types();
    let n_initial = idx.entries().iter().filter(|e| e.is_initial()).count();
    let n_free = idx.entries().iter().filter(|e| e.is_free()).count();
    let n_blessed = idx.entries().iter().filter(|e| e.is_blessed()).count();
    let n_undead = idx.entries().iter().filter(|e| e.is_undead()).count();
    let n_dead = idx.entries().iter().filter(|e| e.is_dead()).count();
    let n_closed = idx.transitions().iter().filter(|t| t.is_closed()).count();
    let n_open = idx.transitions().len() - n_closed;
    eprintln!(
        "[{}] types={} (initial={}, free={}, blessed={}, undead={}, dead={}) transitions={} ({} open, {} closed) time={:.2?}",
        label,
        n_total,
        n_initial,
        n_free,
        n_blessed,
        n_undead,
        n_dead,
        idx.transitions().len(),
        n_open,
        n_closed,
        elapsed,
    );

    let collection = Collection::from_index(&idx, ring_name);
    let t_write = Instant::now();
    let file = std::fs::File::create(output).unwrap_or_else(|e| panic!("create {output}: {e}"));
    serde_json::to_writer(std::io::BufWriter::new(file), &collection)
        .unwrap_or_else(|e| panic!("serialize {output}: {e}"));
    eprintln!(
        "  Wrote {} vtypes, {} transitions to {} in {:.2?}",
        collection.vtypes.len(),
        collection.transitions.len(),
        output,
        t_write.elapsed(),
    );
}

fn typed_tileset<T, F>(collection: &Collection, make_rat: F) -> Arc<TileSet<T>>
where
    T: IsRingOrField,
    F: Fn(&[i8]) -> Rat<T>,
{
    let rats: Vec<Rat<T>> = collection.tile_angles.iter().map(|s| make_rat(s)).collect();
    Arc::new(TileSet::new(rats))
}

fn validate_typed<T: IsRingOrField>(
    collection: &Collection,
    tile_ts: &Arc<TileSet<T>>,
) -> Result<(), String> {
    let t0 = Instant::now();
    eprintln!(
        "  Parsed: ring={}, tiles={}, vtypes={}, transitions={}",
        collection.ring,
        collection.tile_angles.len(),
        collection.vtypes.len(),
        collection.transitions.len(),
    );

    eprintln!("  Phase 1: Reconstructing witnesses...");
    let t1 = Instant::now();
    let witnesses = collection.reconstruct_witnesses(tile_ts)?;
    eprintln!(
        "  Reconstructed {} witnesses in {:.2?}",
        witnesses.len(),
        t1.elapsed(),
    );

    eprintln!("  Phase 2: Verifying vertex types...");
    let t2 = Instant::now();
    let vt_errors = collection.vtype_errors(&witnesses);
    for (id, msg) in vt_errors.iter().take(5) {
        eprintln!("  ERROR: VTYPE {}: {}", id, msg);
    }
    if !vt_errors.is_empty() {
        return Err(format!("{} vertex type mismatches", vt_errors.len()));
    }
    eprintln!(
        "  All {} vertex types verified in {:.2?}",
        collection.vtypes.len(),
        t2.elapsed(),
    );

    eprintln!("  Phase 3: Verifying transitions...");
    let t3 = Instant::now();
    let tr_errors = collection.transition_errors(tile_ts, &witnesses);
    for ((src, dst), msg) in tr_errors.iter().take(5) {
        eprintln!("  ERROR: TRANS {} -> {}: {}", src, dst, msg);
    }
    if !tr_errors.is_empty() {
        return Err(format!("{} transition errors", tr_errors.len()));
    }
    eprintln!(
        "  All {} transitions verified in {:.2?}",
        collection.transitions.len(),
        t3.elapsed(),
    );

    eprintln!("  Phase 4: Completeness check...");
    let t4 = Instant::now();
    let report = collection.completeness_errors(&witnesses);
    for id in report.missing.iter().take(10) {
        eprintln!("  MISSING: VTYPE {} produces unknown junction type", id);
    }
    if !report.is_complete() {
        return Err(format!(
            "Completeness FAILED: {} vertex types produce unknown junction types",
            report.missing.len(),
        ));
    }
    eprintln!(
        "  Completeness: {} match checks passed in {:.2?}",
        report.matches_checked,
        t4.elapsed(),
    );

    eprintln!("  Validation PASSED in {:.2?}", t0.elapsed());
    Ok(())
}

fn main() {
    #[cfg(feature = "pprof")]
    let guard = ProfilerGuardBuilder::default()
        .frequency(1000)
        .blocklist(&["libc", "libgcc", "pthread", "vdso"])
        .build()
        .unwrap();

    let cli = Cli::parse();

    match &cli.command {
        Commands::Collect { tile, output } => match tile {
            TileSetKind::Penrose => {
                let ts = make_ts_10();
                collect_generic(ts, "ZZ10", output, "penrose");
            }
            _ => {
                let ts = make_ts_12(tile);
                let label = format!("{:?}", tile).to_lowercase();
                collect_generic(ts, "ZZ12", output, &label);
            }
        },
        Commands::Validate { input } => {
            eprintln!("=== Validating: {} ===", input);
            let file = std::fs::File::open(input).unwrap_or_else(|e| {
                eprintln!("Open error: {}", e);
                std::process::exit(1);
            });
            let collection: Collection = serde_json::from_reader(std::io::BufReader::new(file))
                .unwrap_or_else(|e| {
                    eprintln!("Parse error: {}", e);
                    std::process::exit(1);
                });

            let result = match collection.ring.as_str() {
                "ZZ12" => {
                    let ts = typed_tileset(&collection, Rat::<ZZ12>::from_slice_unchecked);
                    validate_typed(&collection, &ts)
                }
                "ZZ10" => {
                    let ts = typed_tileset(&collection, Rat::<ZZ10>::from_slice_unchecked);
                    validate_typed(&collection, &ts)
                }
                other => {
                    eprintln!("Unsupported ring: {}", other);
                    std::process::exit(1);
                }
            };
            if let Err(e) = result {
                eprintln!("FAIL: {}", e);
                std::process::exit(1);
            }
            println!("OK");
        }
    }

    #[cfg(feature = "pprof")]
    {
        if let Ok(report) = guard.report().build() {
            let path = "flamegraph_vtype.svg";
            let mut file = std::fs::File::create(path).unwrap();
            report.flamegraph(&mut file).unwrap();
            eprintln!("Flame graph written to {}", path);
        }
    }
}
