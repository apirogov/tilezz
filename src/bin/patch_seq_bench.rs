use std::sync::Arc;
use std::time::Instant;

use clap::{Parser, Subcommand};
use tilezz::cyclotomic::{IsComplex, IsRingOrField, Units, ZZ10, ZZ12};
use tilezz::intgeom::patch::GrowingPatch;
use tilezz::intgeom::rat::Rat;
use tilezz::intgeom::tiles;
use tilezz::intgeom::tileset::TileSet;
use tilezz::intgeom::vtseq_explorer::PatchSeqExplorer;

#[derive(Parser)]
#[command(name = "patch_seq_bench", about = "Benchmark PatchSeqExplorer")]
struct Cli {
    #[command(subcommand)]
    command: BenchCommands,
}

#[derive(Subcommand)]
enum BenchCommands {
    Run {
        #[arg(long, value_enum)]
        tile: TileSetKind,
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

fn run_generic<T: IsComplex + IsRingOrField + Units>(ts: Arc<TileSet<T>>, label: &str) {
    let t0 = Instant::now();
    let explorer = PatchSeqExplorer::new(Arc::clone(&ts));
    let elapsed = t0.elapsed();
    eprintln!(
        "[{}] subseqs={} witnesses={} patches={} k={} time={:.2?}",
        label,
        explorer.num_subseqs(),
        explorer.num_witnesses(),
        explorer.num_patches(),
        explorer.max_subseq_len(),
        elapsed,
    );
}

fn main() {
    let cli = Cli::parse();
    match &cli.command {
        BenchCommands::Run { tile } => match tile {
            TileSetKind::Hex => {
                let rat = Rat::<ZZ12>::try_from(&tiles::hexagon()).unwrap();
                let ts = Arc::new(TileSet::new(vec![rat]));
                run_generic(ts, "hex");
            }
            TileSetKind::Square => {
                let rat = Rat::<ZZ12>::try_from(&tiles::square()).unwrap();
                let ts = Arc::new(TileSet::new(vec![rat]));
                run_generic(ts, "square");
            }
            TileSetKind::Mixed => {
                let sq = Rat::<ZZ12>::try_from(&tiles::square()).unwrap();
                let hex = Rat::<ZZ12>::try_from(&tiles::hexagon()).unwrap();
                let ts = Arc::new(TileSet::new(vec![sq, hex]));
                run_generic(ts, "mixed");
            }
            TileSetKind::Spectre => {
                let rat = Rat::<ZZ12>::try_from(&tiles::spectre()).unwrap();
                let ts = Arc::new(TileSet::new(vec![rat]));
                run_generic(ts, "spectre");
            }
            TileSetKind::Penrose => {
                let n = Rat::<ZZ10>::try_from(&tiles::penrose_p3_narrow()).unwrap();
                let w = Rat::<ZZ10>::try_from(&tiles::penrose_p3_wide()).unwrap();
                let ts = Arc::new(TileSet::new(vec![n, w]));
                run_generic(ts, "penrose");
            }
        },
    }
}
