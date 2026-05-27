//! `tileset_collect` -- unified collect-and-validate binary across the
//! fixed-point classifications we run on tile sets.
//!
//! Three kinds, one CLI:
//!   * `nbhd`  -- `NeighborhoodIndex` (corona / phase-2 classification).
//!   * `vtype` -- `OpenVertexTypeIndex` (vertex-type BFS).
//!   * `seq`   -- `SeqExplorer` (subseq fixed-point enumeration).
//!
//! Each kind already exposes a serializable `Collection` struct plus a
//! validator that rebuilds from saved JSON and cross-checks against a
//! freshly reconstructed tileset. This binary just wraps them in a common
//! CLI envelope (`{kind, payload}`) so `validate` can dispatch by tag.

use std::sync::Arc;
use std::time::Instant;

use clap::{Parser, Subcommand, ValueEnum};
use serde::{Deserialize, Serialize};

use tilezz::analysis::neighborhood::{self, NeighborhoodIndex};
use tilezz::analysis::seq_explorer::{self, check_fixed_point, SeqExplorer};
use tilezz::analysis::vertextypes::{self, OpenVertexTypeIndex};
use tilezz::cyclotomic::{IsRing, ZZ10, ZZ12};
use tilezz::geom::rat::Rat;
use tilezz::geom::tileset::{self, TileSet};
use tilezz::util::profile::ProfileGuard;

#[derive(Parser)]
#[command(
    name = "tileset_collect",
    about = "Collect + validate fixed-point classifications over a tileset"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Run the chosen classification on the chosen tileset. If `--output`
    /// is given, the collected `Collection` is serialized to JSON (wrapped
    /// in a kind-tagged envelope so `validate` can dispatch later).
    Collect {
        #[arg(long, value_enum, default_value = "hex")]
        tileset: TileSetKind,
        #[arg(long, value_enum)]
        kind: CollectKind,
        #[arg(long)]
        output: Option<String>,
        #[arg(long, help = "Flamegraph output path (requires --features pprof)")]
        pprof: Option<String>,
    },
    /// Replay a saved collection file: reconstruct its tileset from the
    /// embedded angle sequences, then run the kind-specific cross-checks.
    Validate {
        #[arg(long)]
        input: String,
    },
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum TileSetKind {
    Hex,
    Square,
    Mixed,
    Tetris,
    Spectre,
    Penrose,
}

impl TileSetKind {
    fn label(&self) -> &'static str {
        match self {
            TileSetKind::Hex => "hex",
            TileSetKind::Square => "square",
            TileSetKind::Mixed => "mixed",
            TileSetKind::Tetris => "tetris",
            TileSetKind::Spectre => "spectre",
            TileSetKind::Penrose => "penrose",
        }
    }
}

#[derive(Clone, Copy, Debug, ValueEnum, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
enum CollectKind {
    Nbhd,
    Vtype,
    Seq,
}

/// Envelope wrapping the kind-specific `Collection` payload, so the
/// `validate` subcommand can dispatch by tag without peeking at the
/// payload schema first.
#[derive(Serialize, Deserialize)]
struct Envelope {
    kind: CollectKind,
    payload: serde_json::Value,
}

// ---------------- tileset construction ----------------

fn ts_zz12(kind: TileSetKind) -> Arc<TileSet<ZZ12>> {
    match kind {
        TileSetKind::Hex => tileset::hex::<ZZ12>(),
        TileSetKind::Square => tileset::square::<ZZ12>(),
        TileSetKind::Mixed => tileset::mixed::<ZZ12>(),
        TileSetKind::Tetris => tileset::tetrominoes::<ZZ12>(),
        TileSetKind::Spectre => tileset::spectre::<ZZ12>(),
        TileSetKind::Penrose => panic!("penrose requires ZZ10"),
    }
}

fn ts_zz10(kind: TileSetKind) -> Arc<TileSet<ZZ10>> {
    match kind {
        TileSetKind::Penrose => tileset::penrose::<ZZ10>(),
        _ => panic!("only penrose uses ZZ10"),
    }
}

fn rebuild_tileset_zz12(tile_angles: &[Vec<i8>]) -> Arc<TileSet<ZZ12>> {
    let rats: Vec<Rat<ZZ12>> = tile_angles
        .iter()
        .map(|s| Rat::<ZZ12>::from_slice_unchecked(s))
        .collect();
    Arc::new(TileSet::new(rats))
}

fn rebuild_tileset_zz10(tile_angles: &[Vec<i8>]) -> Arc<TileSet<ZZ10>> {
    let rats: Vec<Rat<ZZ10>> = tile_angles
        .iter()
        .map(|s| Rat::<ZZ10>::from_slice_unchecked(s))
        .collect();
    Arc::new(TileSet::new(rats))
}

// ---------------- per-kind collectors ----------------

fn collect_nbhd<T: IsRing>(
    ts: Arc<TileSet<T>>,
    ring: &str,
    label: &str,
) -> neighborhood::Collection {
    eprintln!("[{label}] Running neighborhood-type BFS...");
    let t0 = Instant::now();
    let idx = NeighborhoodIndex::new(Arc::clone(&ts));
    let elapsed = t0.elapsed();

    let kinds = idx.classify_all();
    let dead = kinds
        .iter()
        .filter(|&&k| k == neighborhood::NtKind::Dead)
        .count();
    let undead = kinds
        .iter()
        .filter(|&&k| k == neighborhood::NtKind::Undead)
        .count();
    let blessed = kinds
        .iter()
        .filter(|&&k| k == neighborhood::NtKind::Blessed)
        .count();
    let free = kinds
        .iter()
        .filter(|&&k| k == neighborhood::NtKind::Free)
        .count();
    eprintln!(
        "[{label}] types={} (dead={} undead={} blessed={} free={}) transitions={} time={:.2?}",
        idx.num_types(),
        dead,
        undead,
        blessed,
        free,
        idx.transitions().len(),
        elapsed,
    );

    idx.to_collection(ring)
}

fn collect_vtype<T: IsRing>(
    ts: Arc<TileSet<T>>,
    ring: &str,
    label: &str,
) -> vertextypes::Collection {
    eprintln!("[{label}] Running vertex-type BFS...");
    let t0 = Instant::now();
    let idx = OpenVertexTypeIndex::new(Arc::clone(&ts));
    let elapsed = t0.elapsed();

    let entries = idx.entries();
    let n_initial = entries.iter().filter(|e| e.is_initial()).count();
    let n_free = entries.iter().filter(|e| e.is_free()).count();
    let n_blessed = entries.iter().filter(|e| e.is_blessed()).count();
    let n_undead = entries.iter().filter(|e| e.is_undead()).count();
    let n_dead = entries.iter().filter(|e| e.is_dead()).count();
    let n_closed = idx.transitions().iter().filter(|t| t.is_closed()).count();
    let n_open = idx.transitions().len() - n_closed;
    eprintln!(
        "[{label}] types={} (initial={} free={} blessed={} undead={} dead={}) transitions={} ({} open, {} closed) time={:.2?}",
        idx.num_types(),
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

    vertextypes::Collection::from_index(&idx, ring)
}

fn collect_seq<T: IsRing>(
    ts: Arc<TileSet<T>>,
    ring: &str,
    label: &str,
) -> seq_explorer::Collection {
    eprintln!("[{label}] Running subseq fixed-point BFS...");
    let t0 = Instant::now();
    let explorer = SeqExplorer::new(Arc::clone(&ts));
    let elapsed = t0.elapsed();

    eprintln!(
        "[{label}] subseqs={} rats={} k={} time={:.2?}",
        explorer.num_subseqs(),
        explorer.num_rats(),
        explorer.max_subseq_len(),
        elapsed,
    );

    seq_explorer::Collection::from_explorer(&explorer, ring)
}

// ---------------- per-kind validators ----------------

fn validate_nbhd<T: IsRing>(
    coll: neighborhood::Collection,
    ts: Arc<TileSet<T>>,
) -> Result<(), String> {
    eprintln!(
        "  Parsed: ring={}, tiles={}, entries={}, transitions={}, kinds={}",
        coll.ring,
        coll.tile_angles.len(),
        coll.entries.len(),
        coll.transitions.len(),
        coll.kinds.len(),
    );

    let t0 = Instant::now();
    // `from_collection` cross-checks ids, classifies, and rejects on
    // kind mismatches; success means the saved data round-trips against
    // a freshly built index.
    let idx = NeighborhoodIndex::from_collection(ts, coll)?;
    let invalid = idx.validate();
    if !invalid.is_empty() {
        return Err(format!(
            "{} entries failed NeighborhoodIndex::validate (ids: {:?}...)",
            invalid.len(),
            &invalid[..invalid.len().min(5)],
        ));
    }
    eprintln!(
        "  Rebuilt + validated {} entries in {:.2?}",
        idx.num_types(),
        t0.elapsed(),
    );
    Ok(())
}

fn validate_vtype<T: IsRing>(
    coll: vertextypes::Collection,
    ts: Arc<TileSet<T>>,
) -> Result<(), String> {
    let t0 = Instant::now();
    eprintln!(
        "  Parsed: ring={}, tiles={}, vtypes={}, transitions={}",
        coll.ring,
        coll.tile_angles.len(),
        coll.vtypes.len(),
        coll.transitions.len(),
    );

    eprintln!("  Phase 1: Reconstructing witnesses...");
    let t1 = Instant::now();
    let witnesses = coll.reconstruct_witnesses(&ts)?;
    eprintln!(
        "  Reconstructed {} witnesses in {:.2?}",
        witnesses.len(),
        t1.elapsed(),
    );

    eprintln!("  Phase 2: Verifying vertex types...");
    let t2 = Instant::now();
    let vt_errors = coll.vtype_errors(&witnesses);
    for (id, msg) in vt_errors.iter().take(5) {
        eprintln!("  ERROR: VTYPE {}: {}", id, msg);
    }
    if !vt_errors.is_empty() {
        return Err(format!("{} vertex type mismatches", vt_errors.len()));
    }
    eprintln!(
        "  All {} vertex types verified in {:.2?}",
        coll.vtypes.len(),
        t2.elapsed(),
    );

    eprintln!("  Phase 3: Verifying transitions...");
    let t3 = Instant::now();
    let tr_errors = coll.transition_errors(&ts, &witnesses);
    for ((src, dst), msg) in tr_errors.iter().take(5) {
        eprintln!("  ERROR: TRANS {} -> {}: {}", src, dst, msg);
    }
    if !tr_errors.is_empty() {
        return Err(format!("{} transition errors", tr_errors.len()));
    }
    eprintln!(
        "  All {} transitions verified in {:.2?}",
        coll.transitions.len(),
        t3.elapsed(),
    );

    eprintln!("  Phase 4: Completeness check...");
    let t4 = Instant::now();
    let report = coll.completeness_errors(&witnesses);
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

fn validate_seq<T: IsRing>(
    coll: seq_explorer::Collection,
    ts: Arc<TileSet<T>>,
) -> Result<(), String> {
    let t0 = Instant::now();
    let k = ts.rats().iter().map(|r| r.len()).max().unwrap_or(0);

    eprintln!(
        "  Parsed: ring={}, tiles={}, rats={}, subseqs={}, k={}",
        coll.ring,
        coll.tile_angles.len(),
        coll.provenances.len(),
        coll.subseqs.len(),
        k,
    );

    eprintln!("  Replaying glues...");
    let t1 = Instant::now();
    let rats = coll.replay_rats(&ts)?;
    eprintln!("  Replayed {} rats in {:.2?}", rats.len(), t1.elapsed());

    eprintln!("  Checking subseq presence...");
    let t2 = Instant::now();
    let presence_errors = coll.presence_errors(&rats);
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
        coll.subseqs.len(),
        t2.elapsed(),
    );

    eprintln!("  Checking completeness (witnesses x tiles)...");
    let t3 = Instant::now();
    let witness_ids: Vec<usize> = {
        let mut ids: Vec<usize> = coll.subseqs.iter().map(|(id, _)| *id).collect();
        ids.sort_unstable();
        ids.dedup();
        ids
    };
    let known: std::collections::BTreeSet<Vec<i8>> =
        coll.subseqs.iter().map(|(_, s)| s.clone()).collect();

    let report = check_fixed_point(&ts, &rats, &witness_ids, &known, k);
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

// ---------------- collect / validate dispatch ----------------

fn run_collect_zz12(ts: Arc<TileSet<ZZ12>>, kind: CollectKind, label: &str) -> serde_json::Value {
    match kind {
        CollectKind::Nbhd => serde_json::to_value(collect_nbhd(ts, "ZZ12", label)).unwrap(),
        CollectKind::Vtype => serde_json::to_value(collect_vtype(ts, "ZZ12", label)).unwrap(),
        CollectKind::Seq => serde_json::to_value(collect_seq(ts, "ZZ12", label)).unwrap(),
    }
}

fn run_collect_zz10(ts: Arc<TileSet<ZZ10>>, kind: CollectKind, label: &str) -> serde_json::Value {
    match kind {
        CollectKind::Nbhd => serde_json::to_value(collect_nbhd(ts, "ZZ10", label)).unwrap(),
        CollectKind::Vtype => serde_json::to_value(collect_vtype(ts, "ZZ10", label)).unwrap(),
        CollectKind::Seq => serde_json::to_value(collect_seq(ts, "ZZ10", label)).unwrap(),
    }
}

fn run_validate(env: Envelope) -> Result<(), String> {
    // All three Collection types embed `ring: String`. Peek at it before
    // committing to a typed deserialize so we know which ring to rebuild.
    let ring = env
        .payload
        .get("ring")
        .and_then(|v| v.as_str())
        .ok_or("payload missing `ring` field")?
        .to_string();

    let tile_angles: Vec<Vec<i8>> = env
        .payload
        .get("tile_angles")
        .ok_or("payload missing `tile_angles`")?
        .clone()
        .pipe(serde_json::from_value)
        .map_err(|e| e.to_string())?;

    match ring.as_str() {
        "ZZ12" => {
            let ts = rebuild_tileset_zz12(&tile_angles);
            dispatch_validate_zz12(env.kind, env.payload, ts)
        }
        "ZZ10" => {
            let ts = rebuild_tileset_zz10(&tile_angles);
            dispatch_validate_zz10(env.kind, env.payload, ts)
        }
        other => Err(format!("unsupported ring: {other}")),
    }
}

fn dispatch_validate_zz12(
    kind: CollectKind,
    payload: serde_json::Value,
    ts: Arc<TileSet<ZZ12>>,
) -> Result<(), String> {
    match kind {
        CollectKind::Nbhd => {
            let coll: neighborhood::Collection =
                serde_json::from_value(payload).map_err(|e| e.to_string())?;
            validate_nbhd(coll, ts)
        }
        CollectKind::Vtype => {
            let coll: vertextypes::Collection =
                serde_json::from_value(payload).map_err(|e| e.to_string())?;
            validate_vtype(coll, ts)
        }
        CollectKind::Seq => {
            let coll: seq_explorer::Collection =
                serde_json::from_value(payload).map_err(|e| e.to_string())?;
            validate_seq(coll, ts)
        }
    }
}

fn dispatch_validate_zz10(
    kind: CollectKind,
    payload: serde_json::Value,
    ts: Arc<TileSet<ZZ10>>,
) -> Result<(), String> {
    match kind {
        CollectKind::Nbhd => {
            let coll: neighborhood::Collection =
                serde_json::from_value(payload).map_err(|e| e.to_string())?;
            validate_nbhd(coll, ts)
        }
        CollectKind::Vtype => {
            let coll: vertextypes::Collection =
                serde_json::from_value(payload).map_err(|e| e.to_string())?;
            validate_vtype(coll, ts)
        }
        CollectKind::Seq => {
            let coll: seq_explorer::Collection =
                serde_json::from_value(payload).map_err(|e| e.to_string())?;
            validate_seq(coll, ts)
        }
    }
}

/// Tiny extension helper so we can chain `.pipe(f)` on owned values
/// (avoids repeating `let foo = serde_json::from_value(json).map_err(...);`).
trait Pipe: Sized {
    fn pipe<R>(self, f: impl FnOnce(Self) -> R) -> R {
        f(self)
    }
}
impl<T> Pipe for T {}

// ---------------- main ----------------

fn main() {
    let cli = Cli::parse();

    match cli.command {
        Commands::Collect {
            tileset,
            kind,
            output,
            pprof,
        } => {
            let profile = ProfileGuard::start(pprof.as_deref());

            let payload = match tileset {
                TileSetKind::Penrose => run_collect_zz10(ts_zz10(tileset), kind, tileset.label()),
                _ => run_collect_zz12(ts_zz12(tileset), kind, tileset.label()),
            };

            if let Some(path) = output {
                let env = Envelope { kind, payload };
                let t_write = Instant::now();
                let file =
                    std::fs::File::create(&path).unwrap_or_else(|e| panic!("create {path}: {e}"));
                serde_json::to_writer(std::io::BufWriter::new(file), &env)
                    .unwrap_or_else(|e| panic!("serialize {path}: {e}"));
                eprintln!("  Wrote {} in {:.2?}", path, t_write.elapsed());
            }

            profile.finish();
        }
        Commands::Validate { input } => {
            eprintln!("=== Validating: {} ===", input);
            let file = std::fs::File::open(&input).unwrap_or_else(|e| {
                eprintln!("Open error: {e}");
                std::process::exit(1);
            });
            let env: Envelope = serde_json::from_reader(std::io::BufReader::new(file))
                .unwrap_or_else(|e| {
                    eprintln!("Parse error: {e}");
                    std::process::exit(1);
                });

            match run_validate(env) {
                Ok(()) => println!("OK"),
                Err(e) => {
                    eprintln!("FAIL: {e}");
                    std::process::exit(1);
                }
            }
        }
    }
}
