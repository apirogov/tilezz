use std::sync::Arc;
use std::time::Instant;

use clap::{Parser, Subcommand};
use tilezz::cyclotomic::{ZZ10, ZZ12};
use tilezz::intgeom::rat::Rat;
use tilezz::intgeom::tiles;
use tilezz::intgeom::tileset::TileSet;
use tilezz::misc::seq_explorer::{check_fixed_point, Collection, SeqExplorer};

#[derive(Parser)]
#[command(name = "seq_collect", about = "SeqExplorer: collect + validate")]
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
    Spectre,
    Penrose,
}

fn make_ts_12(kind: &TileSetKind) -> Arc<TileSet<ZZ12>> {
    match kind {
        TileSetKind::Hex => {
            let rat = Rat::try_from(&tiles::hexagon::<ZZ12>()).unwrap();
            Arc::new(TileSet::new(vec![rat]))
        }
        TileSetKind::Square => {
            let rat = Rat::try_from(&tiles::square::<ZZ12>()).unwrap();
            Arc::new(TileSet::new(vec![rat]))
        }
        TileSetKind::Mixed => {
            let sq = Rat::try_from(&tiles::square::<ZZ12>()).unwrap();
            let hex = Rat::try_from(&tiles::hexagon::<ZZ12>()).unwrap();
            Arc::new(TileSet::new(vec![sq, hex]))
        }
        TileSetKind::Spectre => {
            let rat = Rat::try_from(&tiles::spectre::<ZZ12>()).unwrap();
            Arc::new(TileSet::new(vec![rat]))
        }
        TileSetKind::Penrose => panic!("penrose requires ZZ10"),
    }
}

fn make_ts_10() -> Arc<TileSet<ZZ10>> {
    let n = Rat::try_from(&tiles::penrose_p3_narrow()).unwrap();
    let w = Rat::try_from(&tiles::penrose_p3_wide()).unwrap();
    Arc::new(TileSet::new(vec![n, w]))
}

/// Rebuild a typed `TileSet<T>` from a `Collection`'s raw tile angle
/// sequences. The closure converts angle slices into the ring's `Rat`.
fn typed_tileset<T, F>(collection: &Collection, make_rat: F) -> Arc<TileSet<T>>
where
    T: tilezz::cyclotomic::IsComplex + tilezz::cyclotomic::IsRingOrField + tilezz::cyclotomic::Units,
    F: Fn(&[i8]) -> Rat<T>,
{
    let rats: Vec<Rat<T>> = collection
        .tile_angles
        .iter()
        .map(|s| make_rat(s))
        .collect();
    Arc::new(TileSet::new(rats))
}

fn validate_typed<T>(collection: &Collection, tile_ts: &Arc<TileSet<T>>) -> Result<(), String>
where
    T: tilezz::cyclotomic::IsComplex + tilezz::cyclotomic::IsRingOrField + tilezz::cyclotomic::Units,
{
    let t0 = Instant::now();
    let k = tile_ts.rats().iter().map(|r| r.len()).max().unwrap_or(0);

    eprintln!(
        "  Parsed: ring={}, tiles={}, rats={}, subseqs={}, k={}",
        collection.ring,
        collection.tile_angles.len(),
        collection.provenances.len(),
        collection.subseqs.len(),
        k,
    );

    eprintln!("  Replaying glues...");
    let t1 = Instant::now();
    let rats = collection.replay_rats(tile_ts)?;
    eprintln!("  Replayed {} rats in {:.2?}", rats.len(), t1.elapsed());

    eprintln!("  Checking subseq presence...");
    let t2 = Instant::now();
    let presence_errors = collection.presence_errors(&rats);
    for (rat_id, seq) in presence_errors.iter().take(5) {
        eprintln!(
            "  ERROR: subseq {:?} not found in rat {} (len={})",
            seq,
            rat_id,
            rats[*rat_id].len()
        );
    }
    if !presence_errors.is_empty() {
        return Err(format!("{} subseq presence errors", presence_errors.len()));
    }
    eprintln!(
        "  All {} subseqs verified present in {:.2?}",
        collection.subseqs.len(),
        t2.elapsed(),
    );

    eprintln!("  Checking completeness (witnesses × tiles)...");
    let t3 = Instant::now();
    let witness_ids: Vec<usize> = {
        let mut ids: Vec<usize> = collection.subseqs.iter().map(|(id, _)| *id).collect();
        ids.sort_unstable();
        ids.dedup();
        ids
    };
    let known: std::collections::BTreeSet<Vec<i8>> =
        collection.subseqs.iter().map(|(_, s)| s.clone()).collect();

    let report = check_fixed_point(tile_ts, &rats, &witness_ids, &known, k);
    for sub in report.missing.iter().take(10) {
        eprintln!("  MISSING: subseq {:?} not in collection", sub);
    }
    eprintln!(
        "  Completeness: {} matches checked, {} new subseqs found in {:.2?}",
        report.matches_checked,
        report.missing.len(),
        t3.elapsed(),
    );
    if !report.is_complete() {
        return Err(format!(
            "Completeness check FAILED: {} new subseqs not in collection",
            report.missing.len()
        ));
    }

    eprintln!("  Validation PASSED in {:.2?}", t0.elapsed());
    Ok(())
}

fn main() {
    #[cfg(feature = "pprof")]
    let guard = pprof::ProfilerGuardBuilder::default()
        .frequency(1000)
        .blocklist(&["libc", "libgcc", "pthread", "vdso"])
        .build()
        .unwrap();

    let cli = Cli::parse();

    match &cli.command {
        Commands::Collect { tile, output } => {
            let collection = match tile {
                TileSetKind::Penrose => {
                    let ts = make_ts_10();
                    let t0 = Instant::now();
                    let explorer = SeqExplorer::new(ts);
                    let elapsed = t0.elapsed();
                    eprintln!(
                        "[penrose] subseqs={} witnesses={} rats={} k={} time={:.2?}",
                        explorer.num_subseqs(),
                        explorer.num_contributing_rats(),
                        explorer.num_known_rats(),
                        explorer.max_subseq_len(),
                        elapsed,
                    );
                    Collection::from_explorer(&explorer, "ZZ10")
                }
                _ => {
                    let ts = make_ts_12(tile);
                    let t0 = Instant::now();
                    let explorer = SeqExplorer::new(ts);
                    let elapsed = t0.elapsed();
                    eprintln!(
                        "[{:?}] subseqs={} witnesses={} rats={} k={} time={:.2?}",
                        tile,
                        explorer.num_subseqs(),
                        explorer.num_contributing_rats(),
                        explorer.num_known_rats(),
                        explorer.max_subseq_len(),
                        elapsed,
                    );
                    Collection::from_explorer(&explorer, "ZZ12")
                }
            };

            let t_write = Instant::now();
            let file = std::fs::File::create(output)
                .unwrap_or_else(|e| panic!("create {output}: {e}"));
            serde_json::to_writer(std::io::BufWriter::new(file), &collection)
                .unwrap_or_else(|e| panic!("serialize {output}: {e}"));
            eprintln!(
                "  Wrote {} subseqs to {} in {:.2?}",
                collection.subseqs.len(),
                output,
                t_write.elapsed(),
            );
        }
        Commands::Validate { input } => {
            eprintln!("=== Validating: {} ===", input);
            let file = std::fs::File::open(input).unwrap_or_else(|e| {
                eprintln!("Open error: {}", e);
                std::process::exit(1);
            });
            let collection: Collection =
                serde_json::from_reader(std::io::BufReader::new(file)).unwrap_or_else(|e| {
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
            let path = "flamegraph_seq.svg";
            let mut file = std::fs::File::create(path).unwrap();
            report.flamegraph(&mut file).unwrap();
            eprintln!("Flame graph written to {}", path);

            let mut text_report = Vec::new();
            use std::io::Write;
            report.text_detailed(&mut text_report).unwrap();
            let text_path = "pprof_seq.txt";
            std::fs::write(text_path, &text_report).unwrap();
            eprintln!("Text report written to {}", text_path);
        }
    }
}
