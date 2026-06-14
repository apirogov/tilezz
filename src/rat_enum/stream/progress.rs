//! Live progress telemetry for Stage 1 (`--mode stream`).
//!
//! Picture: the stream stage is a pool of worker threads, each chewing
//! through a shared list of DFS seeds. Without instrumentation the run
//! is a black box for hours. This module surfaces two complementary
//! views:
//!
//!   * a periodic *heartbeat* -- a one-line snapshot of aggregate work
//!     rate, recursion depth, seeds dispatched/retired, and scratch
//!     growth, printed every `interval` while the run is live; and
//!   * a per-seed *cost distribution* -- how long each seed's subtree
//!     took and how many closures it yielded -- summarized at the end,
//!     so we can see whether the seeds are roughly equal work or
//!     heavy-tailed (the load-balancing / "are seeds equivalent?"
//!     question).
//!
//! Workers publish into a lock-free [`WorkerCell`] (one per worker,
//! cache-line padded to avoid false sharing); a monitor thread reads
//! them. Publishing is cadence-gated inside the DFS (every ~`2^20`
//! events) so the hot path stays plain integer adds -- see
//! `rat_enum_step`.

use std::path::Path;
use std::sync::atomic::{AtomicBool, AtomicU8, AtomicU32, AtomicU64, AtomicUsize, Ordering};
use std::time::{Duration, Instant};

use crate::rat_enum::seed::parallel::branch_factor;
use crate::rat_enum::stream::runs::list_run_files;

/// Worker lifecycle states, published as a `u8` so the monitor can
/// render "in-flight vs done" without a lock.
pub const STATE_IDLE: u8 = 0;
pub const STATE_RUNNING: u8 = 1;
pub const STATE_DONE: u8 = 2;

/// Flush cadence mask: a worker republishes its live snapshot whenever
/// its event counter crosses a multiple of `FLUSH_MASK + 1` (~1M). The
/// hot path only pays a single `& FLUSH_MASK == 0` test per branch.
pub const FLUSH_MASK: u64 = (1 << 20) - 1;

/// Lock-free per-worker published snapshot. Exactly one worker writes
/// its own cell; the monitor only reads. We use `Relaxed` everywhere --
/// this is a cheap live gauge, not a synchronization edge, and a
/// slightly stale read is fine for a heartbeat.
#[repr(align(64))]
#[derive(Default)]
pub struct WorkerCell {
    /// Total categorized DFS branch attempts so far (`DfsStats::total`).
    pub attempts: AtomicU64,
    /// Closures recorded so far by this worker.
    pub closed: AtomicU64,
    /// Current recursion depth (snake length) at last publish.
    pub depth: AtomicU32,
    /// Index of the seed this worker is currently draining.
    pub seed_idx: AtomicU32,
    /// Length of the current seed's prefix (so the depth above the seed
    /// can be subtracted off when reading the in-cylinder odometer).
    pub seed_len: AtomicU32,
    /// Run-relative wall-clock (ms since run start) at which the worker
    /// picked up its current seed. Powers the per-seed ETA.
    pub seed_start_ms: AtomicU64,
    /// Fraction of the current seed's cylinder swept so far, in parts
    /// per million (see [`odometer_fraction`]).
    pub progress_ppm: AtomicU32,
    /// One of [`STATE_IDLE`] / [`STATE_RUNNING`] / [`STATE_DONE`].
    pub state: AtomicU8,
}

impl WorkerCell {
    /// Publish a fresh snapshot (`Relaxed` stores). Callers gate this
    /// behind a [`FLUSH_MASK`] cadence test so the DFS hot path only
    /// pays the store cost ~once per million events.
    #[inline]
    pub fn publish(&self, attempts: u64, closed: u64, depth: u32, progress_ppm: u32) {
        self.attempts.store(attempts, Ordering::Relaxed);
        self.closed.store(closed, Ordering::Relaxed);
        self.depth.store(depth, Ordering::Relaxed);
        self.progress_ppm.store(progress_ppm, Ordering::Relaxed);
    }
}

/// Cost of draining one seed's DFS subtree, measured by the worker that
/// owned it: wall-clock and the number of closures it produced.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SeedCost {
    pub elapsed_ns: u64,
    pub closures: u64,
}

/// Distribution summary of per-seed costs. `skew = max / median`
/// elapsed is the headline number for "are the seeds roughly equal
/// work (skew ~1) or heavy-tailed (skew >> 1)?".
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SeedCostSummary {
    pub n: usize,
    pub min_ns: u64,
    pub median_ns: u64,
    pub p90_ns: u64,
    pub max_ns: u64,
    pub total_closures: u64,
    pub skew: f64,
}

impl SeedCostSummary {
    /// Summarize a batch of per-seed costs. Returns `None` for an empty
    /// batch (nothing to summarize).
    pub fn from_costs(costs: &[SeedCost]) -> Option<SeedCostSummary> {
        if costs.is_empty() {
            return None;
        }
        let mut elapsed: Vec<u64> = costs.iter().map(|c| c.elapsed_ns).collect();
        elapsed.sort_unstable();
        let min_ns = elapsed[0];
        let max_ns = *elapsed.last().unwrap();
        let median_ns = percentile(&elapsed, 0.5).unwrap_or(0);
        let p90_ns = percentile(&elapsed, 0.9).unwrap_or(0);
        let skew = if median_ns == 0 {
            0.0
        } else {
            max_ns as f64 / median_ns as f64
        };
        Some(SeedCostSummary {
            n: costs.len(),
            min_ns,
            median_ns,
            p90_ns,
            max_ns,
            total_closures: costs.iter().map(|c| c.closures).sum(),
            skew,
        })
    }
}

/// Nearest-rank percentile element of an ascending-sorted slice (`p` in
/// `[0, 1]`). No interpolation -- returns an actual element. Generic so
/// the per-seed cost summary (u64 ns), the live recursion-depth spread
/// (u32) and the cylinder-sweep spread (f64) all share one convention.
/// `None` for an empty slice.
fn percentile<T: Copy>(sorted: &[T], p: f64) -> Option<T> {
    if sorted.is_empty() {
        return None;
    }
    let n = sorted.len();
    let rank = (p * n as f64).ceil() as usize;
    Some(sorted[rank.saturating_sub(1).min(n - 1)])
}

/// Cumulative counters sampled at one instant; the monitor diffs two
/// of these to derive per-second rates.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Snapshot {
    pub attempts: u64,
    pub closed: u64,
    pub seeds_done: u64,
    pub runs_bytes: u64,
}

/// Per-second rates derived from two [`Snapshot`]s and the elapsed
/// interval between them.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Rates {
    pub attempts_per_s: f64,
    pub closed_per_s: f64,
    pub bytes_per_s: f64,
}

impl Rates {
    /// Per-second deltas from `prev` to `cur` over `dt`. A zero (or
    /// absurdly small) interval yields all-zero rates rather than NaN/inf.
    pub fn between(prev: &Snapshot, cur: &Snapshot, dt: Duration) -> Rates {
        let secs = dt.as_secs_f64();
        if secs <= 0.0 {
            return Rates {
                attempts_per_s: 0.0,
                closed_per_s: 0.0,
                bytes_per_s: 0.0,
            };
        }
        let rate = |now: u64, then: u64| now.saturating_sub(then) as f64 / secs;
        Rates {
            attempts_per_s: rate(cur.attempts, prev.attempts),
            closed_per_s: rate(cur.closed, prev.closed),
            bytes_per_s: rate(cur.runs_bytes, prev.runs_bytes),
        }
    }
}

/// Format a byte count with binary (IEC) units and one decimal place
/// above 1 KiB (e.g. `6.4 GiB`); plain `N B` below.
pub fn human_bytes(n: u64) -> String {
    const UNITS: [&str; 5] = ["B", "KiB", "MiB", "GiB", "TiB"];
    if n < 1024 {
        return format!("{n} B");
    }
    let mut val = n as f64;
    let mut unit = 0;
    while val >= 1024.0 && unit < UNITS.len() - 1 {
        val /= 1024.0;
        unit += 1;
    }
    format!("{val:.1} {}", UNITS[unit])
}

/// Format a duration as `HH:MM:SS` (hours may exceed two digits for
/// multi-day runs).
pub fn fmt_hms(d: Duration) -> String {
    let s = d.as_secs();
    format!("{:02}:{:02}:{:02}", s / 3600, (s % 3600) / 60, s % 60)
}

/// Adaptive short duration: `ns` / `us` / `ms` / `s`, one decimal above
/// the base unit. Used for the per-seed cost spread, whose values span
/// sub-millisecond (near-empty cylinders) to minutes (the giants) -- a
/// fixed `ms` unit would truncate the small end to `0`.
pub fn fmt_dur(d: Duration) -> String {
    let ns = d.as_nanos();
    if ns < 1_000 {
        format!("{ns} ns")
    } else if ns < 1_000_000 {
        format!("{:.1} us", ns as f64 / 1e3)
    } else if ns < 1_000_000_000 {
        format!("{:.1} ms", ns as f64 / 1e6)
    } else {
        format!("{:.1} s", ns as f64 / 1e9)
    }
}

/// Per-seed remaining-time estimate from the uniform-odometer model:
/// `on_seed * (1 - frac) / frac`. Returns `None` when `frac` is outside
/// `[0.005, 1.0)` -- too small to extrapolate, or (defensively) `>= 1`,
/// which would otherwise drive the divide negative and panic
/// `Duration::from_secs_f64`. Capped so a tiny `frac` can't yield an
/// absurd `Duration`. Address-space-uniform (see [`odometer_fraction`]):
/// rough early, trustworthy as `frac -> 1`.
pub fn seed_eta(on_seed: Duration, frac: f64) -> Option<Duration> {
    if !(0.005..1.0).contains(&frac) {
        return None;
    }
    let secs = (on_seed.as_secs_f64() * (1.0 - frac) / frac).min(1e9);
    Some(Duration::from_secs_f64(secs))
}

/// Fraction of one seed's cylinder swept so far, read as a base-`B`
/// odometer over the in-seed turn digits (`rel_path` = the path below
/// the seed prefix; `hturn`/`step` are the ring's). `B` and the digit
/// weights come from [`branch_factor`] -- the same fan-out the seed
/// split uses -- so the two can't drift. `idx(turn) = (hturn-1 - turn)
/// / step` is the turn's position in the fixed descending branch order,
/// so digit 0 is the lexicographic start. Returns a value in `[0, 1)`
/// -- the *address-space* (uniform-measure) fraction of the cylinder:
/// monotone in DFS progress but, because the tree is non-uniform, not
/// linear in work or leaf count.
pub fn odometer_fraction(rel_path: &[i8], hturn: i8, step: i8) -> f64 {
    let b = branch_factor(hturn, step) as f64;
    let hm1 = (hturn - 1) as i32;
    let step = step.max(1) as i32;
    let mut frac = 0.0;
    let mut scale = 1.0 / b;
    for &turn in rel_path {
        let idx = ((hm1 - turn as i32) / step) as f64;
        frac += idx * scale;
        scale /= b;
    }
    frac
}

/// Sum the on-disk size of the Stage-1 run files. Best-effort: any file
/// that vanished mid-stat (rotated/replaced) just contributes 0.
fn runs_dir_bytes(runs_dir: &Path) -> u64 {
    list_run_files(runs_dir)
        .map(|files| {
            files
                .iter()
                .filter_map(|p| std::fs::metadata(p).ok())
                .map(|m| m.len())
                .sum()
        })
        .unwrap_or(0)
}

/// Live heartbeat loop, run on its own thread. Wakes every `interval`,
/// aggregates the worker board into a [`Snapshot`], and prints one
/// stderr line with per-second rates, recursion-depth spread, seed
/// progress, and scratch growth. Exits one tick after `done` flips
/// (printing a final line so the closing state is visible).
///
/// The sleep is sliced so a `done` flip is noticed within ~200ms rather
/// than after a full `interval`.
#[allow(clippy::too_many_arguments)]
pub fn run_monitor(
    board: &[WorkerCell],
    next_idx: &AtomicUsize,
    completed: &AtomicUsize,
    seeds_total: usize,
    runs_dir: &Path,
    started: Instant,
    interval: Duration,
    done: &AtomicBool,
) {
    eprintln!(
        "[hb] monitoring {seeds_total} seeds across {} workers, interval {}s",
        board.len(),
        interval.as_secs_f64(),
    );
    let mut prev = Snapshot::default();
    let mut prev_t = started;
    loop {
        let mut slept = Duration::ZERO;
        while slept < interval && !done.load(Ordering::Relaxed) {
            let slice = Duration::from_millis(200).min(interval - slept);
            std::thread::sleep(slice);
            slept += slice;
        }
        let finished = done.load(Ordering::Relaxed);

        let now_ms = started.elapsed().as_millis() as u64;
        let mut attempts = 0u64;
        let mut closed = 0u64;
        let mut running = 0usize;
        let mut depths: Vec<u32> = Vec::with_capacity(board.len());
        let mut sweeps: Vec<f64> = Vec::with_capacity(board.len());
        // The "hardest" in-flight seed is the least-swept one -- the
        // giant a worker is currently grinding. Track it for the ETA.
        let mut hardest: Option<(u32, f64, u64)> = None; // (seed_idx, frac, start_ms)
        for c in board {
            attempts += c.attempts.load(Ordering::Relaxed);
            closed += c.closed.load(Ordering::Relaxed);
            if c.state.load(Ordering::Relaxed) == STATE_RUNNING {
                running += 1;
                depths.push(c.depth.load(Ordering::Relaxed));
                let f = c.progress_ppm.load(Ordering::Relaxed) as f64 / 1e6;
                sweeps.push(f);
                let idx = c.seed_idx.load(Ordering::Relaxed);
                let start = c.seed_start_ms.load(Ordering::Relaxed);
                if hardest.is_none_or(|(_, hf, _)| f < hf) {
                    hardest = Some((idx, f, start));
                }
            }
        }
        depths.sort_unstable();
        let d_p50 = percentile(&depths, 0.5).unwrap_or(0);
        let d_max = depths.last().copied().unwrap_or(0);
        sweeps.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let s_min = sweeps.first().copied().unwrap_or(0.0);
        let s_p50 = percentile(&sweeps, 0.5).unwrap_or(0.0);
        let s_max = sweeps.last().copied().unwrap_or(0.0);

        let now = Instant::now();
        let cur = Snapshot {
            attempts,
            closed,
            seeds_done: completed.load(Ordering::Relaxed) as u64,
            runs_bytes: runs_dir_bytes(runs_dir),
        };
        let rates = Rates::between(&prev, &cur, now.duration_since(prev_t));

        eprintln!(
            "[hb {}] seeds {}/{} disp, {} done ({} in-flight) | \
             closures {} (+{:.2}M/s) | attempts {} (+{:.0}M/s) | \
             depth p50={} max={} | runs {} (+{}/s)",
            fmt_hms(now.duration_since(started)),
            next_idx.load(Ordering::Relaxed).min(seeds_total),
            seeds_total,
            cur.seeds_done,
            running,
            cur.closed,
            rates.closed_per_s / 1e6,
            cur.attempts,
            rates.attempts_per_s / 1e6,
            d_p50,
            d_max,
            human_bytes(cur.runs_bytes),
            human_bytes(rates.bytes_per_s as u64),
        );
        // Per-seed cylinder sweep + ETA for the hardest in-flight seed.
        // The ETA is address-space-uniform (see `odometer_fraction`):
        // trustworthy as f -> 1, rough early.
        eprintln!(
            "         seed-sweep min={:.2}% p50={:.1}% max={:.1}% | {}",
            s_min * 100.0,
            s_p50 * 100.0,
            s_max * 100.0,
            match hardest {
                None => "no seeds in flight".to_string(),
                // `idx` is the dispatch index into the seed list (0-based).
                Some((idx, f, start)) => {
                    let on_seed = Duration::from_millis(now_ms.saturating_sub(start));
                    match seed_eta(on_seed, f) {
                        Some(eta) => format!(
                            "hardest seed {}/{}: {:.2}% swept, {} in, eta ~{}",
                            idx,
                            seeds_total,
                            f * 100.0,
                            fmt_hms(on_seed),
                            fmt_hms(eta),
                        ),
                        None => format!(
                            "hardest seed {}/{}: {:.3}% swept, {} in, eta --",
                            idx,
                            seeds_total,
                            f * 100.0,
                            fmt_hms(on_seed),
                        ),
                    }
                }
            },
        );

        prev = cur;
        prev_t = now;
        if finished {
            break;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn seed_cost_summary_reports_percentiles_and_skew() {
        // Elapsed 10,20,30,80; nearest-rank over n=4:
        //   median p=0.5 -> rank ceil(2.0)=2 -> idx 1 -> 20
        //   p90    p=0.9 -> rank ceil(3.6)=4 -> idx 3 -> 80
        let costs = [
            SeedCost {
                elapsed_ns: 30,
                closures: 3,
            },
            SeedCost {
                elapsed_ns: 10,
                closures: 1,
            },
            SeedCost {
                elapsed_ns: 80,
                closures: 8,
            },
            SeedCost {
                elapsed_ns: 20,
                closures: 2,
            },
        ];
        let s = SeedCostSummary::from_costs(&costs).expect("non-empty batch");
        assert_eq!(s.n, 4);
        assert_eq!(s.min_ns, 10);
        assert_eq!(s.median_ns, 20);
        assert_eq!(s.p90_ns, 80);
        assert_eq!(s.max_ns, 80);
        assert_eq!(s.total_closures, 14);
        assert!(
            (s.skew - 4.0).abs() < 1e-9,
            "skew = max/median = 80/20 = 4, got {}",
            s.skew
        );
    }

    #[test]
    fn seed_cost_summary_empty_is_none() {
        assert!(SeedCostSummary::from_costs(&[]).is_none());
    }

    #[test]
    fn rates_are_per_second_deltas() {
        let prev = Snapshot {
            attempts: 1_000,
            closed: 100,
            seeds_done: 1,
            runs_bytes: 1_000_000,
        };
        let cur = Snapshot {
            attempts: 5_000,
            closed: 300,
            seeds_done: 4,
            runs_bytes: 3_000_000,
        };
        let r = Rates::between(&prev, &cur, Duration::from_secs(2));
        assert!((r.attempts_per_s - 2000.0).abs() < 1e-9);
        assert!((r.closed_per_s - 100.0).abs() < 1e-9);
        assert!((r.bytes_per_s - 1_000_000.0).abs() < 1e-9);
    }

    #[test]
    fn rates_zero_interval_is_zero_not_nan() {
        let snap = Snapshot {
            attempts: 5_000,
            closed: 300,
            seeds_done: 4,
            runs_bytes: 3_000_000,
        };
        let r = Rates::between(&snap, &snap, Duration::ZERO);
        assert_eq!(r.attempts_per_s, 0.0);
        assert_eq!(r.closed_per_s, 0.0);
        assert_eq!(r.bytes_per_s, 0.0);
    }

    #[test]
    fn human_bytes_uses_binary_units() {
        assert_eq!(human_bytes(0), "0 B");
        assert_eq!(human_bytes(1023), "1023 B");
        assert_eq!(human_bytes(1024), "1.0 KiB");
        assert_eq!(human_bytes(1536), "1.5 KiB");
        assert_eq!(human_bytes(1024 * 1024), "1.0 MiB");
        assert_eq!(human_bytes(3 * 1024 * 1024 * 1024), "3.0 GiB");
    }

    #[test]
    fn fmt_hms_formats_hours_minutes_seconds() {
        assert_eq!(fmt_hms(Duration::from_secs(0)), "00:00:00");
        assert_eq!(fmt_hms(Duration::from_secs(61)), "00:01:01");
        assert_eq!(fmt_hms(Duration::from_secs(3723)), "01:02:03");
    }

    #[test]
    fn odometer_is_zero_at_the_lexicographic_start() {
        // ZZ7 via ring-14 step-2: directions [6,4,2,0,-2,-4,-6], B=7.
        // The first branch (turn = +6, idx 0) at every level is the
        // lexicographic start of the cylinder -> swept fraction 0.
        assert_eq!(odometer_fraction(&[], 7, 2), 0.0);
        assert_eq!(odometer_fraction(&[6], 7, 2), 0.0);
        assert_eq!(odometer_fraction(&[6, 6, 6], 7, 2), 0.0);
    }

    #[test]
    fn odometer_reads_digits_in_descending_branch_order() {
        // turn -6 is the LAST branch (idx (6-(-6))/2 = 6); at depth 0
        // that's 6/7 of the cylinder already behind us.
        assert!((odometer_fraction(&[-6], 7, 2) - 6.0 / 7.0).abs() < 1e-12);
        // turn +4 is the 2nd branch (idx 1) -> 1/7.
        assert!((odometer_fraction(&[4], 7, 2) - 1.0 / 7.0).abs() < 1e-12);
        // [6, 4]: 0 at depth 0, then idx1 at depth 1 -> 1/49.
        assert!((odometer_fraction(&[6, 4], 7, 2) - 1.0 / 49.0).abs() < 1e-12);
    }

    #[test]
    fn odometer_is_monotone_and_below_one() {
        // Deepening into later branches only increases the reading, and
        // it never reaches 1.0 (the final branch is still in progress).
        let a = odometer_fraction(&[6, 4], 7, 2);
        let b = odometer_fraction(&[6, -6], 7, 2);
        assert!(a < b, "later branch must read higher: {a} < {b}");
        let last = odometer_fraction(&[-6, -6], 7, 2);
        assert!(last < 1.0 && last > 0.97, "near the end but < 1: {last}");
    }

    #[test]
    fn odometer_handles_step_one() {
        // step=1, hturn=4 (hm1=3) -> directions [3,2,1,0,-1,-2,-3], B=7,
        // idx(v) = 3 - v.
        assert_eq!(odometer_fraction(&[3], 4, 1), 0.0);
        assert!((odometer_fraction(&[-3], 4, 1) - 6.0 / 7.0).abs() < 1e-12);
        assert!((odometer_fraction(&[2], 4, 1) - 1.0 / 7.0).abs() < 1e-12);
    }

    #[test]
    fn percentile_picks_nearest_rank_element() {
        let v = [10u64, 20, 30, 80]; // ascending
        assert_eq!(percentile(&v, 0.0), Some(10));
        assert_eq!(percentile(&v, 0.5), Some(20)); // ceil(2.0)=2 -> idx 1
        assert_eq!(percentile(&v, 0.9), Some(80)); // ceil(3.6)=4 -> idx 3
        assert_eq!(percentile(&v, 1.0), Some(80));
        assert_eq!(percentile::<u64>(&[], 0.5), None);
        assert_eq!(percentile(&[42u64], 0.5), Some(42));
    }

    #[test]
    fn seed_eta_extrapolates_and_guards_edges() {
        // Half-swept after 60s -> ~60s remaining.
        let e = seed_eta(Duration::from_secs(60), 0.5).unwrap();
        assert!((e.as_secs_f64() - 60.0).abs() < 1e-6);
        // Too-small fraction -> None (can't extrapolate).
        assert_eq!(seed_eta(Duration::from_secs(60), 0.0), None);
        assert_eq!(seed_eta(Duration::from_secs(60), 0.001), None);
        // frac >= 1 must NOT panic (the divide would go negative and
        // Duration::from_secs_f64 would panic).
        assert_eq!(seed_eta(Duration::from_secs(60), 1.0), None);
        assert_eq!(seed_eta(Duration::from_secs(60), 1.5), None);
        // Near-done -> small positive remaining.
        let near = seed_eta(Duration::from_secs(100), 0.99).unwrap();
        assert!(near.as_secs_f64() > 0.0 && near.as_secs_f64() < 2.0);
    }

    #[test]
    fn fmt_dur_is_adaptive() {
        assert_eq!(fmt_dur(Duration::from_nanos(500)), "500 ns");
        assert_eq!(fmt_dur(Duration::from_micros(400)), "400.0 us");
        assert_eq!(fmt_dur(Duration::from_micros(1500)), "1.5 ms");
        assert_eq!(fmt_dur(Duration::from_millis(2100)), "2.1 s");
        assert_eq!(fmt_dur(Duration::from_secs(56)), "56.0 s");
    }
}
