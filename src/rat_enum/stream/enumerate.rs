//! Stage 1 of the streaming pipeline: parallel DFS where each
//! worker pushes closures into a per-thread sort buffer that flushes
//! to `out_dir/runs/run_tNN_rMM.bin` instead of accumulating them in
//! an in-memory HashSet.
//!
//! Memory profile: bounded by `n_threads × buffer_size`, regardless
//! of the final rat count. ZZ12 n=15 (~230M rats) becomes feasible
//! on a 16 GB workstation; the same workload via `--mode bench`
//! needs ~30 GB.
//!
//! Caveat: the buffer-local dedup is best-effort. Records may
//! reappear across different runs (worker A finds a rat that
//! worker B also finds). Stage 2's k-way merge does the final
//! global dedup.

use std::path::Path;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::thread;
use std::time::{Duration, Instant};

use crate::cyclotomic::IsRing;
use crate::geom::snake::Snake;
use crate::rat_enum::canonical::make_ops;
use crate::rat_enum::dfs::{SeedGather, collect_seeds, rat_enum_step};
use crate::rat_enum::prune::{Prunes, snapshot_prunes};
use crate::rat_enum::seed::parallel::{branch_factor, splitting_depth};
use crate::rat_enum::stats::DfsStats;
use crate::rat_enum::stream::progress::{
    STATE_DONE, STATE_RUNNING, SeedCost, SeedCostSummary, WorkerCell, fmt_dur, run_monitor,
};
use crate::rat_enum::stream::runs::RunWriter;

/// Filesystem layout: `out_dir/runs/run_tNN_rMM.bin`. Created on first flush.
pub const RUNS_SUBDIR: &str = "runs";

/// Parallel streaming enumeration over ring `ZZ`. Each worker owns a
/// [`RunWriter`] pointed at `out_dir/runs/`. The seed walk's
/// early-closure stream goes into a thread-id=0 writer (so all run
/// files share the same shape).
///
/// Returns the combined `DfsStats` once every worker has finished
/// and flushed its buffer.
#[allow(clippy::too_many_arguments)]
pub fn stream_enum_parallel<ZZ: IsRing>(
    max_steps: usize,
    step: i8,
    n_threads: usize,
    free: bool,
    paranoid: bool,
    prunes: &Prunes,
    out_dir: &Path,
    heartbeat: Option<Duration>,
) -> std::io::Result<DfsStats> {
    let ops = make_ops(free);
    let runs_dir = out_dir.join(RUNS_SUBDIR);
    std::fs::create_dir_all(&runs_dir)?;

    let label = if free {
        "free stream"
    } else {
        "rotation stream"
    };
    println!("-------- {label} (out_dir={}) --------", out_dir.display());
    if paranoid {
        println!("paranoid: per-step fresh-snake cross-check enabled");
    }

    let branching = branch_factor(ZZ::hturn(), step);
    let split_depth = splitting_depth(n_threads.max(1), branching);
    println!("stream: n_threads={n_threads} branching={branching} split_depth={split_depth}");

    // Seed walk -- alive prefixes for workers, plus a per-thread-0
    // writer for any polygons that close before reaching split_depth.
    let mut seeds: Vec<Vec<i8>> = Vec::new();
    let mut seed_stats = DfsStats::default();
    {
        let mut seed_writer = RunWriter::new(&runs_dir, 0);
        let mut snake: Snake<ZZ> = Snake::new();
        let mut record_closed = |seq: &[i8]| seed_writer.record(seq);
        let mut gather = SeedGather {
            seeds: &mut seeds,
            record_closed: &mut record_closed,
            stats: &mut seed_stats,
        };
        collect_seeds::<ZZ>(
            &mut snake,
            max_steps,
            step,
            split_depth,
            &mut gather,
            ops,
            paranoid,
            prunes,
        );
        // seed_writer drops here, flushing its buffer.
    }
    println!("stream: {} seed states collected", seeds.len());

    // Parallel workers + an optional live monitor share these. Workers
    // consume seeds via the shared atomic counter, own a `RunWriter`
    // keyed by thread index, publish live telemetry into their own
    // `WorkerCell`, and record per-seed cost as each seed retires.
    let next_idx = AtomicUsize::new(0);
    let next_ref = &next_idx;
    let completed = AtomicUsize::new(0);
    let completed_ref = &completed;
    let runs_dir_ref = &runs_dir;
    let seeds_ref: &[Vec<i8>] = &seeds;
    let n_workers = n_threads.max(1);
    let seeds_total = seeds.len();

    let board: Vec<WorkerCell> = (0..n_workers).map(|_| WorkerCell::default()).collect();
    let board_ref = &board;
    let monitor_done = AtomicBool::new(false);
    let monitor_done_ref = &monitor_done;
    let started = Instant::now();

    let (worker_stats, seed_costs): (Vec<DfsStats>, Vec<SeedCost>) = thread::scope(|s| {
        // Live heartbeat (opt-in via --heartbeat); runs until the
        // workers finish and flip `monitor_done`.
        if let Some(interval) = heartbeat {
            s.spawn(move || {
                run_monitor(
                    board_ref,
                    next_ref,
                    completed_ref,
                    seeds_total,
                    runs_dir_ref,
                    started,
                    interval,
                    monitor_done_ref,
                );
            });
        }

        let mut handles = Vec::with_capacity(n_workers);
        for (worker_id, cell) in board_ref.iter().enumerate() {
            // Thread ids start at 1 -- thread 0 is reserved for the
            // seed-walk early closures above.
            let tid = worker_id + 1;
            handles.push(s.spawn(move || -> (DfsStats, Vec<SeedCost>) {
                let mut local_stats = DfsStats::default();
                let mut costs: Vec<SeedCost> = Vec::new();
                let mut writer = RunWriter::new(runs_dir_ref, tid);
                loop {
                    let i = next_ref.fetch_add(1, Ordering::Relaxed);
                    if i >= seeds_ref.len() {
                        break;
                    }
                    cell.seed_idx.store(i as u32, Ordering::Relaxed);
                    cell.seed_len
                        .store(seeds_ref[i].len() as u32, Ordering::Relaxed);
                    cell.seed_start_ms
                        .store(started.elapsed().as_millis() as u64, Ordering::Relaxed);
                    cell.progress_ppm.store(0, Ordering::Relaxed);
                    cell.state.store(STATE_RUNNING, Ordering::Relaxed);
                    let closed_before = local_stats.closed;
                    let t0 = Instant::now();
                    let mut snake: Snake<ZZ> = Snake::from_slice_unsafe(&seeds_ref[i]);
                    let mut record = |seq: &[i8]| writer.record(seq);
                    rat_enum_step::<ZZ>(
                        &mut snake,
                        max_steps,
                        step,
                        &mut record,
                        &mut local_stats,
                        ops,
                        paranoid,
                        prunes,
                        Some(cell),
                    );
                    costs.push(SeedCost {
                        elapsed_ns: t0.elapsed().as_nanos() as u64,
                        closures: local_stats.closed - closed_before,
                    });
                    completed_ref.fetch_add(1, Ordering::Relaxed);
                }
                drop(writer); // explicit; flushes remaining buffer
                cell.state.store(STATE_DONE, Ordering::Relaxed);
                (local_stats, costs)
            }));
        }
        // Join every worker FIRST (collecting Results, not unwrapping), then
        // signal the monitor, then propagate any panic. Doing this before the
        // unwrap matters: `thread::scope` joins the monitor on exit, so if a
        // worker panicked and we unwound before flipping `monitor_done`, the
        // monitor's loop would never terminate and the scope would deadlock.
        let joined: Vec<std::thread::Result<(DfsStats, Vec<SeedCost>)>> =
            handles.into_iter().map(|h| h.join()).collect();
        monitor_done_ref.store(true, Ordering::Relaxed);
        let mut stats_acc: Vec<DfsStats> = Vec::with_capacity(n_workers);
        let mut costs_acc: Vec<SeedCost> = Vec::new();
        for r in joined {
            let (st, mut costs) = r.expect("worker panic");
            stats_acc.push(st);
            costs_acc.append(&mut costs);
        }
        (stats_acc, costs_acc)
    });

    let mut total_stats = seed_stats;
    for ws in &worker_stats {
        total_stats.merge(ws);
    }

    // Per-seed cost distribution: the empirical answer to "are the seeds
    // roughly equal work, or heavy-tailed?" (skew = max / median).
    if let Some(sum) = SeedCostSummary::from_costs(&seed_costs) {
        println!(
            "stream: per-seed cost over {} seeds -- elapsed min/median/p90/max = \
             {}/{}/{}/{}, skew(max/median)={:.1}x, closures total={}",
            sum.n,
            fmt_dur(Duration::from_nanos(sum.min_ns)),
            fmt_dur(Duration::from_nanos(sum.median_ns)),
            fmt_dur(Duration::from_nanos(sum.p90_ns)),
            fmt_dur(Duration::from_nanos(sum.max_ns)),
            sum.skew,
            sum.total_closures,
        );
    }

    // Inventory the runs we produced. Useful for users / Stage 2.
    let run_files = crate::rat_enum::stream::runs::list_run_files(&runs_dir)?;
    let total_bytes: u64 = run_files
        .iter()
        .filter_map(|p| std::fs::metadata(p).ok())
        .map(|m| m.len())
        .sum();
    println!(
        "stream: wrote {} run file(s), {} bytes total",
        run_files.len(),
        total_bytes
    );

    Ok(total_stats)
}

/// Runtime-ring dispatcher for [`stream_enum_parallel`]. Snapshots
/// the global prune state, picks the typed ring impl, runs the
/// stream.
#[allow(clippy::too_many_arguments)]
pub fn stream_enum_dispatch(
    ring: u8,
    max_steps: usize,
    step: i8,
    n_threads: usize,
    free: bool,
    paranoid: bool,
    out_dir: &Path,
    heartbeat_secs: u64,
) -> std::io::Result<DfsStats> {
    let prunes = snapshot_prunes();
    let n = n_threads.max(1);
    let heartbeat = (heartbeat_secs > 0).then(|| Duration::from_secs(heartbeat_secs));
    crate::dispatch_ring!(
        ring,
        stream_enum_parallel::<ZZ>(
            max_steps, step, n, free, paranoid, &prunes, out_dir, heartbeat
        )
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic::{ZZ8, ZZ12};
    use crate::rat_enum::canonical::make_ops as canonical_make_ops;
    use crate::rat_enum::dfs::rat_enum_with;
    use crate::rat_enum::prune::Prunes;
    use crate::rat_enum::stream::merge::{UNIQUE_FILENAME, merge_runs, read_unique_records};
    use std::path::PathBuf;
    use std::sync::atomic::{AtomicUsize, Ordering as AOrd};

    fn tempdir() -> PathBuf {
        static C: AtomicUsize = AtomicUsize::new(0);
        let n = C.fetch_add(1, AOrd::Relaxed);
        let pid = std::process::id();
        let path = std::env::temp_dir().join(format!("rat_enum_stream_e2e_{pid}_{n}"));
        std::fs::create_dir_all(&path).unwrap();
        path
    }

    /// Sort a baseline `Vec<Vec<i8>>` into the (length asc, lex asc)
    /// order that `unique.bin` is in.
    fn sort_by_len_then_lex(mut v: Vec<Vec<i8>>) -> Vec<Vec<i8>> {
        v.sort_by(|a, b| a.len().cmp(&b.len()).then_with(|| a.cmp(b)));
        v
    }

    /// Drive stream + merge for the given ring, then compare the
    /// recovered set against the baseline DFS. The streaming pipeline
    /// must produce exactly the same set of canonical rats as
    /// `rat_enum_with`, in the same (length, lex) order.
    fn check_stream_matches_baseline<ZZ: crate::cyclotomic::IsRing>(
        ring: u8,
        max_steps: usize,
        free: bool,
    ) {
        let dir = tempdir();
        let prunes = Prunes::default();

        let stats = stream_enum_parallel::<ZZ>(max_steps, 1, 4, free, false, &prunes, &dir, None)
            .expect("stream_enum_parallel");
        assert!(stats.closed > 0, "no closures recorded -- did Stage 1 run?");

        let cert = merge_runs(&dir, ring, max_steps, 1, free).expect("merge_runs");
        assert_eq!(
            cert.ring, ring,
            "certificate.ring does not match the request"
        );

        let from_stream: Vec<Vec<i8>> = read_unique_records(&dir.join(UNIQUE_FILENAME))
            .unwrap()
            .map(|r| r.unwrap())
            .collect();
        assert_eq!(
            from_stream.len(),
            cert.unique_records as usize,
            "read_unique_records count diverges from certificate"
        );

        let ops = canonical_make_ops(free);
        let (baseline, _) = rat_enum_with::<ZZ>(max_steps, 1, ops, "baseline", "", false, &prunes);
        let expected = sort_by_len_then_lex(baseline);

        assert_eq!(
            from_stream.len(),
            expected.len(),
            "stream/baseline cardinality mismatch (ZZ{ring} n={max_steps} free={free}): \
             {} vs {}",
            from_stream.len(),
            expected.len(),
        );
        assert_eq!(
            from_stream, expected,
            "stream/baseline content mismatch (ZZ{ring} n={max_steps} free={free})"
        );

        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn stream_matches_baseline_zz8_n10_rotation() {
        check_stream_matches_baseline::<ZZ8>(8, 10, false);
    }

    #[test]
    fn stream_matches_baseline_zz8_n10_free() {
        check_stream_matches_baseline::<ZZ8>(8, 10, true);
    }

    #[test]
    fn stream_matches_baseline_zz12_n8_rotation() {
        check_stream_matches_baseline::<ZZ12>(12, 8, false);
    }

    #[test]
    fn stream_matches_baseline_zz12_n8_free() {
        check_stream_matches_baseline::<ZZ12>(12, 8, true);
    }

    /// Stage 3 end-to-end: stream -> merge -> build a streaming
    /// RatDafsa via `from_sorted_unique_rats`, and check it's
    /// observationally identical to a buffering `from_rats` built
    /// from the baseline DFS. Guards against any drift between the
    /// two RatDafsa constructors and against the streaming pipeline
    /// silently producing rats in the wrong order (which the
    /// `from_sorted_unique_rats` debug-assert would catch first).
    fn check_stream_build_matches_baseline<ZZ: crate::cyclotomic::IsRing>(
        ring: u8,
        max_steps: usize,
        free: bool,
    ) {
        use crate::stringmatch::RatDafsa;

        let dir = tempdir();
        let prunes = Prunes::default();

        stream_enum_parallel::<ZZ>(max_steps, 1, 4, free, false, &prunes, &dir, None)
            .expect("stream_enum_parallel");
        merge_runs(&dir, ring, max_steps, 1, free).expect("merge_runs");

        // Streaming build: feed unique.bin's records directly into
        // `from_sorted_unique_rats`. No Vec<Vec<i8>> in the middle.
        let records = read_unique_records(&dir.join(UNIQUE_FILENAME))
            .unwrap()
            .map(|r| r.unwrap());
        let streamed_dafsa = RatDafsa::from_sorted_unique_rats(records);

        // Reference: baseline DFS -> buffering `from_rats`.
        let ops = canonical_make_ops(free);
        let (baseline, _) = rat_enum_with::<ZZ>(max_steps, 1, ops, "baseline", "", false, &prunes);
        let buffered_dafsa = RatDafsa::from_rats(baseline.iter().map(|v| v.as_slice()));

        // Headline checks: same count, same (length, lex) iteration,
        // same index_of for every rat.
        assert_eq!(
            streamed_dafsa.len(),
            buffered_dafsa.len(),
            "stream-build/baseline cardinality mismatch (ZZ{ring} n={max_steps} free={free})"
        );
        let streamed_iter: Vec<Vec<i8>> = streamed_dafsa.iter().collect();
        let buffered_iter: Vec<Vec<i8>> = buffered_dafsa.iter().collect();
        assert_eq!(
            streamed_iter, buffered_iter,
            "stream-build/baseline iter mismatch (ZZ{ring} n={max_steps} free={free})"
        );
        for rat in &streamed_iter {
            assert_eq!(
                streamed_dafsa.index_of(rat.as_slice()),
                buffered_dafsa.index_of(rat.as_slice()),
                "index_of mismatch for {:?}",
                rat
            );
        }

        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn stream_build_matches_baseline_zz8_n10_free() {
        check_stream_build_matches_baseline::<ZZ8>(8, 10, true);
    }

    #[test]
    fn stream_build_matches_baseline_zz12_n8_free() {
        check_stream_build_matches_baseline::<ZZ12>(12, 8, true);
    }

    #[test]
    fn stream_build_matches_baseline_zz12_n8_rotation() {
        check_stream_build_matches_baseline::<ZZ12>(12, 8, false);
    }

    /// Force the volume-fan-out paths that a default-threshold run at
    /// small n never reaches, and which only grow at large n: a tiny
    /// per-writer buffer so each writer flushes into MANY run files
    /// (exercising cross-flush local dedup), plus every rat written
    /// through two different writers so the same record lands in
    /// different run files and the k-way merge must collapse it. The
    /// merged, rebuilt set must still equal the in-memory baseline
    /// exactly. This logic is ring-independent -- the rings only change
    /// the byte values, not the encode/sort/merge/dedup path -- so one
    /// ring suffices; what matters here is the run-file fan-out, not
    /// the choice of ring.
    #[test]
    fn stream_merge_dedups_across_many_runs_and_flushes() {
        use crate::rat_enum::stream::runs::{RunWriter, list_run_files};
        use crate::stringmatch::RatDafsa;

        let dir = tempdir();
        let runs_dir = dir.join(RUNS_SUBDIR);
        std::fs::create_dir_all(&runs_dir).unwrap();
        let prunes = Prunes::default();

        // Real canonical rats from the in-memory engine (ZZ12 n8 = 517).
        let ops = canonical_make_ops(true);
        let (baseline, _) = rat_enum_with::<ZZ12>(8, 1, ops, "baseline", "", false, &prunes);
        assert!(baseline.len() > 50, "need a non-trivial set to fan out");

        // threshold=7 -> dozens of flushes per writer; 3 writers; each
        // rat recorded into two distinct writers -> a duplicate sitting
        // in a different run file that only the global merge can dedup.
        {
            let mut writers: Vec<RunWriter> = (0..3)
                .map(|tid| RunWriter::with_threshold(&runs_dir, tid, 7))
                .collect();
            for (i, rat) in baseline.iter().enumerate() {
                writers[i % 3].record(rat);
                writers[(i + 1) % 3].record(rat);
            }
            // writers flush remaining buffers on Drop here.
        }

        let files = list_run_files(&runs_dir).unwrap();
        assert!(
            files.len() > 3,
            "tiny threshold should fan out into many run files, got {}",
            files.len()
        );

        let cert = merge_runs(&dir, 12, 8, 1, true).expect("merge_runs");
        assert_eq!(
            cert.unique_records as usize,
            baseline.len(),
            "k-way merge must collapse the duplicated records back to baseline cardinality"
        );

        let records = read_unique_records(&dir.join(UNIQUE_FILENAME))
            .unwrap()
            .map(|r| r.unwrap());
        let streamed = RatDafsa::from_sorted_unique_rats(records);
        let baseline_dafsa = RatDafsa::from_rats(baseline.iter().map(|v| v.as_slice()));

        assert_eq!(streamed.len(), baseline_dafsa.len(), "cardinality mismatch");
        let s: Vec<Vec<i8>> = streamed.iter().collect();
        let b: Vec<Vec<i8>> = baseline_dafsa.iter().collect();
        assert_eq!(s, b, "fan-out stream-merge set != in-memory baseline");
        for rat in &s {
            assert_eq!(
                streamed.index_of(rat.as_slice()),
                baseline_dafsa.index_of(rat.as_slice()),
                "index_of mismatch for {rat:?}"
            );
        }

        std::fs::remove_dir_all(&dir).ok();
    }

    /// Re-running the full pipeline (stream -> merge -> build) into
    /// the same output directory must produce byte-identical
    /// artifacts on the second invocation. Guards against any stage
    /// leaving stale state behind, file-creation races, or
    /// nondeterministic ordering inside the streaming dafsa builder.
    #[test]
    fn pipeline_idempotent_rerun_zz8_n8_free() {
        use crate::stringmatch::RatDafsa;

        let dir = tempdir();
        let prunes = Prunes::default();
        let ring = 8u8;
        let max_steps = 8;

        // First full run.
        stream_enum_parallel::<ZZ8>(max_steps, 1, 2, true, false, &prunes, &dir, None)
            .expect("stream pass 1");
        let cert1 = merge_runs(&dir, ring, max_steps, 1, true).expect("merge pass 1");
        let unique_bytes_1 = std::fs::read(dir.join(UNIQUE_FILENAME)).unwrap();
        let recs1: Vec<Vec<i8>> = read_unique_records(&dir.join(UNIQUE_FILENAME))
            .unwrap()
            .map(|r| r.unwrap())
            .collect();
        let dafsa1 = RatDafsa::from_sorted_unique_rats(recs1.iter().map(|v| v.as_slice()));
        let dafsa1_blocks_dir = dir.join("dafsa");
        std::fs::create_dir_all(&dafsa1_blocks_dir).unwrap();
        dafsa1
            .write_blocks(&dafsa1_blocks_dir, 8)
            .expect("build pass 1");
        let manifest_1 = std::fs::read(dafsa1_blocks_dir.join("block_index.json")).unwrap();

        // Second full run -- same output directory, same params. The
        // stream stage re-runs the DFS into the same runs/ dir
        // (which already has files from pass 1), so we wipe runs/
        // first to model the "user re-runs from scratch" workflow.
        std::fs::remove_dir_all(dir.join(super::RUNS_SUBDIR)).ok();
        stream_enum_parallel::<ZZ8>(max_steps, 1, 2, true, false, &prunes, &dir, None)
            .expect("stream pass 2");
        let cert2 = merge_runs(&dir, ring, max_steps, 1, true).expect("merge pass 2");
        let unique_bytes_2 = std::fs::read(dir.join(UNIQUE_FILENAME)).unwrap();
        let recs2: Vec<Vec<i8>> = read_unique_records(&dir.join(UNIQUE_FILENAME))
            .unwrap()
            .map(|r| r.unwrap())
            .collect();
        let dafsa2 = RatDafsa::from_sorted_unique_rats(recs2.iter().map(|v| v.as_slice()));
        dafsa2
            .write_blocks(&dafsa1_blocks_dir, 8)
            .expect("build pass 2");
        let manifest_2 = std::fs::read(dafsa1_blocks_dir.join("block_index.json")).unwrap();

        // unique.bin and certificate.blake3 must match across runs.
        assert_eq!(
            cert1.unique_blake3, cert2.unique_blake3,
            "certificate BLAKE3 differs across reruns"
        );
        assert_eq!(
            unique_bytes_1, unique_bytes_2,
            "unique.bin differs across reruns"
        );
        assert_eq!(cert1.unique_records, cert2.unique_records);

        // The block index manifest is structurally deterministic
        // (block IDs, states-per-block boundaries) and must round-trip
        // byte-for-byte on a same-params rerun.
        assert_eq!(
            manifest_1, manifest_2,
            "block_index.json differs across reruns"
        );

        // Spot-check the first block file too; if the manifest is
        // identical and the block writer is deterministic, every
        // block file (content-addressed by SHA-256) must exist on
        // disk under `blocks/`.
        let manifest: crate::stringmatch::dafsa::lazy::BlockManifest =
            serde_json::from_slice(&manifest_1).unwrap();
        assert!(!manifest.blocks.is_empty(), "no blocks emitted");
        let first = &manifest.blocks[0];
        let block_0_path = dafsa1_blocks_dir.join(manifest.block_filename(first));
        let block_0_bytes = std::fs::read(&block_0_path).unwrap();
        assert!(
            !block_0_bytes.is_empty(),
            "first block file missing: {block_0_path:?}"
        );

        let _ = std::fs::remove_dir_all(&dir);
    }

    /// Building without a prior `merge` step (so no `unique.bin`)
    /// must produce a clearly-typed I/O error rather than a panic or
    /// a silent wrong-shape DAFSA. The CLI uses this to print a
    /// helpful "run --mode merge first" message.
    #[test]
    fn read_unique_records_errors_when_missing() {
        let dir = tempdir();
        // No unique.bin at this path.
        let missing = dir.join(UNIQUE_FILENAME);
        let err = read_unique_records(&missing).unwrap_err();
        assert_eq!(err.kind(), std::io::ErrorKind::NotFound);
        let _ = std::fs::remove_dir_all(&dir);
    }
}
