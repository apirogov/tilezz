//! Load a saved `SegSeqBFS` snapshot and run the rule analysis
//! (terminals + cursed + nonlocal classification).
//!
//! Usage:
//!   cargo run --release --bin segseqbfs_analyze -- --tile spectre --snapshot snapshots/spectre.bin

use std::path::PathBuf;
use std::sync::Arc;
use std::time::Instant;

use clap::Parser;
use tilezz::cyclotomic::{IsComplex, IsRingOrField, Units, ZZ10, ZZ12};
use tilezz::intgeom::rat::Rat;
use tilezz::intgeom::segrewrite::{RuleAnalysis, SegRewriteTable};
use tilezz::intgeom::segseqbfs::SegSeqBFS;
use tilezz::intgeom::tiles;
use tilezz::intgeom::tileset::TileSet;

#[derive(Parser)]
#[command(
    name = "segseqbfs_analyze",
    about = "Load a SegSeqBFS snapshot and run terminal/cursed/nonlocal analysis"
)]
struct Cli {
    #[arg(long, value_enum)]
    tile: TileSetKind,
    #[arg(long)]
    snapshot: PathBuf,
}

#[derive(Clone, Debug, clap::ValueEnum)]
enum TileSetKind {
    Hex,
    Square,
    Mixed,
    Spectre,
    Penrose,
}

fn analyze<T: IsComplex + IsRingOrField + Units>(
    tileset: Arc<TileSet<T>>,
    label: &str,
    snapshot: &PathBuf,
) {
    eprintln!("[{}] loading snapshot from {}", label, snapshot.display());
    let t0 = Instant::now();
    let bfs = match SegSeqBFS::load_from(tileset, snapshot) {
        Ok(b) => b,
        Err(e) => {
            eprintln!("[{}] load failed: {}", label, e);
            std::process::exit(1);
        }
    };
    eprintln!(
        "[{}] loaded in {:.1}s: sequences={}, rules={}, witnesses={}",
        label,
        t0.elapsed().as_secs_f64(),
        bfs.num_sequences(),
        bfs.num_rules(),
        bfs.num_witness_patches(),
    );

    let t1 = Instant::now();
    let table = SegRewriteTable::build_with_stats(&bfs);
    eprintln!(
        "[{}] built rewrite table in {:.1}s",
        label,
        t1.elapsed().as_secs_f64()
    );

    let t2 = Instant::now();
    let analysis = RuleAnalysis::build(&bfs, &table);
    eprintln!(
        "[{}] analysis in {:.1}s:",
        label,
        t2.elapsed().as_secs_f64()
    );
    eprintln!(
        "  terminals     = {}",
        analysis.terminals.len()
    );
    eprintln!(
        "  dead_rhs      = {}",
        analysis.dead_rhs.len()
    );
    eprintln!(
        "  ctx_sens_rhs  = {}",
        analysis.context_sensitive_rhs.len()
    );
    // Diagnostic: how many dead RHSes also contain a terminal? If all
    // of them do, the dead-rhs check is subsumed by terminal-detection
    // for this tileset.
    let dead_with_terminal = analysis
        .dead_rhs
        .iter()
        .filter(|rhs| rhs.iter().any(|s| analysis.terminals.contains(s)))
        .count();
    eprintln!(
        "  dead_rhs containing a terminal seg: {}/{}",
        dead_with_terminal,
        analysis.dead_rhs.len()
    );
    eprintln!(
        "  cursed_rules  = {}/{} ({:.1}%)",
        analysis.cursed_rules.len(),
        table.num_rules(),
        100.0 * analysis.cursed_rules.len() as f64 / table.num_rules() as f64,
    );
    eprintln!(
        "  cursed_lhs    = {}/{} ({:.1}%)",
        analysis.cursed_lhs.len(),
        table.num_lhs(),
        100.0 * analysis.cursed_lhs.len() as f64 / table.num_lhs() as f64,
    );
    eprintln!(
        "  nonlocal_rules = {}/{} ({:.1}%)",
        analysis.nonlocal_rules.len(),
        table.num_rules(),
        100.0 * analysis.nonlocal_rules.len() as f64 / table.num_rules() as f64,
    );
}

fn main() {
    let cli = Cli::parse();
    match cli.tile {
        TileSetKind::Hex => {
            let rat = Rat::<ZZ12>::try_from(&tiles::hexagon()).unwrap();
            let ts = Arc::new(TileSet::new(vec![rat]));
            analyze(ts, "hex", &cli.snapshot);
        }
        TileSetKind::Square => {
            let rat = Rat::<ZZ12>::try_from(&tiles::square()).unwrap();
            let ts = Arc::new(TileSet::new(vec![rat]));
            analyze(ts, "square", &cli.snapshot);
        }
        TileSetKind::Mixed => {
            let sq = Rat::<ZZ12>::try_from(&tiles::square()).unwrap();
            let hex = Rat::<ZZ12>::try_from(&tiles::hexagon()).unwrap();
            let ts = Arc::new(TileSet::new(vec![sq, hex]));
            analyze(ts, "mixed", &cli.snapshot);
        }
        TileSetKind::Spectre => {
            let rat = Rat::<ZZ12>::try_from(&tiles::spectre()).unwrap();
            let ts = Arc::new(TileSet::new(vec![rat]));
            analyze(ts, "spectre", &cli.snapshot);
        }
        TileSetKind::Penrose => {
            let n = Rat::<ZZ10>::try_from(&tiles::penrose_p3_narrow()).unwrap();
            let w = Rat::<ZZ10>::try_from(&tiles::penrose_p3_wide()).unwrap();
            let ts = Arc::new(TileSet::new(vec![n, w]));
            analyze(ts, "penrose", &cli.snapshot);
        }
    }
}
