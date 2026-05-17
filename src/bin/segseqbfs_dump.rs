//! Run a `SegSeqBFS` for a chosen tileset and dump the result to a
//! file so analyses (terminal/cursed classification, witness
//! inspection, etc.) can load it in seconds rather than re-running
//! the 10-15 minute BFS each time.
//!
//! Progress lines are emitted by the BFS itself via `eprintln!`.
//!
//! Usage:
//!   cargo run --release --bin segseqbfs_dump -- --tile spectre --out /tmp/spectre.bin

use std::path::PathBuf;
use std::sync::Arc;
use std::time::Instant;

use clap::Parser;
use tilezz::cyclotomic::{IsComplex, IsRingOrField, Units, ZZ10, ZZ12};
use tilezz::intgeom::rat::Rat;
use tilezz::intgeom::segseqbfs::SegSeqBFS;
use tilezz::intgeom::tiles;
use tilezz::intgeom::tileset::TileSet;

#[derive(Parser)]
#[command(
    name = "segseqbfs_dump",
    about = "Run SegSeqBFS for a tileset and serialize the result for later analysis"
)]
struct Cli {
    /// Tileset to enumerate.
    #[arg(long, value_enum)]
    tile: TileSetKind,
    /// Output file (binary, bincode-encoded snapshot).
    #[arg(long)]
    out: PathBuf,
    /// Cap on distinct sequences (safety bound).
    #[arg(long, default_value_t = 100_000_000)]
    seq_cap: usize,
    /// Cap on distinct enqueued patches (safety bound).
    #[arg(long, default_value_t = 100_000_000)]
    patch_cap: usize,
}

#[derive(Clone, Debug, clap::ValueEnum)]
enum TileSetKind {
    Hex,
    Square,
    Mixed,
    Spectre,
    Penrose,
}

fn run_and_save<T: IsComplex + IsRingOrField + Units>(
    tileset: Arc<TileSet<T>>,
    label: &str,
    out: &PathBuf,
    seq_cap: usize,
    patch_cap: usize,
) {
    eprintln!("[{}] running SegSeqBFS (seq_cap={}, patch_cap={})", label, seq_cap, patch_cap);
    let t0 = Instant::now();
    let bfs = match SegSeqBFS::run(tileset, seq_cap, patch_cap) {
        Ok(b) => b,
        Err(e) => {
            eprintln!("[{}] BFS failed: {:?}", label, e);
            std::process::exit(1);
        }
    };
    let elapsed = t0.elapsed();
    eprintln!(
        "[{}] done in {:.1}s: sequences={}, rules={}, witnesses={}, patches_seen={}",
        label,
        elapsed.as_secs_f64(),
        bfs.num_sequences(),
        bfs.num_rules(),
        bfs.num_witness_patches(),
        bfs.num_patches_seen(),
    );

    eprintln!("[{}] saving to {}", label, out.display());
    let t1 = Instant::now();
    if let Err(e) = bfs.save_to(out) {
        eprintln!("[{}] save failed: {}", label, e);
        std::process::exit(1);
    }
    let bytes = std::fs::metadata(out).map(|m| m.len()).unwrap_or(0);
    eprintln!(
        "[{}] wrote {} bytes ({:.1} MB) in {:.1}s",
        label,
        bytes,
        bytes as f64 / 1_048_576.0,
        t1.elapsed().as_secs_f64(),
    );
}

fn main() {
    let cli = Cli::parse();
    match cli.tile {
        TileSetKind::Hex => {
            let rat = Rat::<ZZ12>::try_from(&tiles::hexagon()).unwrap();
            let ts = Arc::new(TileSet::new(vec![rat]));
            run_and_save(ts, "hex", &cli.out, cli.seq_cap, cli.patch_cap);
        }
        TileSetKind::Square => {
            let rat = Rat::<ZZ12>::try_from(&tiles::square()).unwrap();
            let ts = Arc::new(TileSet::new(vec![rat]));
            run_and_save(ts, "square", &cli.out, cli.seq_cap, cli.patch_cap);
        }
        TileSetKind::Mixed => {
            let sq = Rat::<ZZ12>::try_from(&tiles::square()).unwrap();
            let hex = Rat::<ZZ12>::try_from(&tiles::hexagon()).unwrap();
            let ts = Arc::new(TileSet::new(vec![sq, hex]));
            run_and_save(ts, "mixed", &cli.out, cli.seq_cap, cli.patch_cap);
        }
        TileSetKind::Spectre => {
            let rat = Rat::<ZZ12>::try_from(&tiles::spectre()).unwrap();
            let ts = Arc::new(TileSet::new(vec![rat]));
            run_and_save(ts, "spectre", &cli.out, cli.seq_cap, cli.patch_cap);
        }
        TileSetKind::Penrose => {
            let n = Rat::<ZZ10>::try_from(&tiles::penrose_p3_narrow()).unwrap();
            let w = Rat::<ZZ10>::try_from(&tiles::penrose_p3_wide()).unwrap();
            let ts = Arc::new(TileSet::new(vec![n, w]));
            run_and_save(ts, "penrose", &cli.out, cli.seq_cap, cli.patch_cap);
        }
    }
}
