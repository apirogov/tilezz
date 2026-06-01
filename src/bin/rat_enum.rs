//! `rat_enum`: enumerate every simple polygon on a cyclotomic-ring
//! lattice with boundary length up to `n`. Output mode (Render/Bench/
//! ListSeeds) is controlled by `--mode`.
//!
//! # Memory profile and how to scale to large `n`
//!
//! The result set is a `HashSet<Vec<i8>>` of canonical (or
//! dihedral-canonical) angle sequences. The *raw* sequence data is
//! tiny (`n` bytes per polygon), but several effects multiply the
//! peak memory of a single process by ~50-100x:
//!
//! 1. **HashSet bucket overhead.** Hashbrown stores each entry with
//!    a control byte, a cached hash, and the `Vec<i8>` header (24
//!    bytes pointer/len/cap). The `Vec`'s heap allocation is
//!    malloc-rounded (~16-32 bytes for an n=13 payload). With load
//!    factor ~0.875, one fully-populated set is ~80 bytes per
//!    polygon -- already 6x the raw payload.
//!
//! 2. **Hashbrown grows by doubling.** Inside a resize from 2^k to
//!    2^(k+1) buckets, *both* arrays are alive while entries are
//!    copied. A set that just resized briefly holds ~3x its final
//!    structure size. When 16 worker threads are all near the same
//!    fill level (atomic-counter dispatch keeps work balanced),
//!    several can resize concurrently.
//!
//! 3. **Per-thread duplication during parallel DFS.** In
//!    `rat_enum_parallel`, each of `n_threads` workers maintains its
//!    own thread-local `HashSet`. These sets are live for the entire
//!    DFS phase and only consolidated at the end. Peak during DFS:
//!    `n_threads * per-thread-final-size * resize-factor`.
//!
//! 4. **Merge double-buffer.** Once workers complete, the main
//!    thread folds each local into a `merged` set via `extend`.
//!    The local being drained is freed only at end of the
//!    for-iteration; the others sit in their `JoinHandle`s until
//!    consumed. Right before the last drain: 1 unconsumed local,
//!    plus `merged` at ~final size, plus a transient resize buffer
//!    inside `merged` -- another ~2x at peak.
//!
//! 5. **glibc malloc arena retention.** Each worker thread's malloc
//!    arena keeps its high-water-mark pages mapped, even after the
//!    transient `Vec<i8>` allocations from inside the close-event
//!    closure are freed. This isn't huge per thread, but it does
//!    add up at scale.
//!
//! Empirically at ZZ12 dihedral n=13 (4.08M polygons, ~50 MB raw),
//! single-process 16-thread peak is **~5.5 GB** -- a ~100x blowup.
//!
//! # Recommendation: separate processes
//!
//! **Running N independent single-threaded processes is strictly
//! better than running one process with `--threads N`** for memory:
//!
//! * Each process holds only its share of the polygons (~1x final
//!   set size for that share, plus the resize transient: peak ~3x).
//! * No per-thread duplication: there's only one HashSet per
//!   process.
//! * No cross-thread merge double-buffer: each process's final set
//!   *is* the final set for that share.
//! * On process exit, the OS reclaims **every page** (arenas
//!   included). Sequential or wall-clock-parallel batches never
//!   share memory across batches.
//!
//! Mechanics: use `--mode list-seeds` to emit a list of work-unit
//! prefixes, then run one `--seed <prefix> --threads 1` process per
//! prefix. An external orchestrator (`xargs -P N`, GNU `parallel`,
//! Slurm, etc.) picks how many run concurrently -- and that becomes
//! the only knob that affects total wall time. Per-process peak
//! memory is *unchanged* by orchestrator parallelism, because each
//! process holds only its own share.
//!
//! Example for n=14 at ZZ12 with 16 cores:
//!
//! ```sh
//! ./target/release/rat_enum --ring 12 -n 14 --dihedral \
//!     --mode list-seeds > listing.txt
//! grep '^SEED ' listing.txt | sed 's/SEED //' | \
//!   xargs -I{} -P 16 \
//!     ./target/release/rat_enum --ring 12 -n 14 --dihedral \
//!     --seed "{}" --threads 1 \
//!   > all_rats.txt
//! { grep '^RAT ' listing.txt; cat all_rats.txt; } \
//!   | sort -u | grep -c '^RAT '
//! ```
//!
//! Peak memory: ~16 processes * ~50-100 MB each ≈ ~1 GB total,
//! versus ~30-40 GB for the equivalent single-process run.
//!
//! `--seed <prefix> --threads N>1` exists (it routes through
//! `enumerate_from_seed_parallel`) for environments where an
//! orchestrator gives one big process at a time, but it has the
//! same per-process memory profile as `--mode bench --threads N` --
//! it does NOT recover the single-threaded-per-process savings.

use std::collections::HashSet;
use std::fs::File;
use std::io::BufWriter;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::thread;
use std::time::Instant;

use clap::{Parser, ValueEnum};

use tilezz::cyclotomic::*;
use tilezz::geom::rat::{lex_min_rot, Rat};
use tilezz::geom::snake::{Snake, Turtle};
use tilezz::stringmatch::{repetition_factor, RatDafsa};
use tilezz::vis::animation::render_gif;
use tilezz::vis::draw::{MarkerStyle, TileStyle};
use tilezz::vis::plotutils::P64;
use tilezz::vis::scene::{Color, Fill, Scene, Stroke, TextStyle, Viewport};

static VERBOSE: Mutex<bool> = Mutex::new(false);

/// When `true`, the DFS streams each newly-discovered rat to stdout as
/// `RAT [...]` on the fly. This is the protocol used by `--seed` and
/// `--mode bench`; modes like `--mode dafsa` that write a single binary
/// artifact set this to `false` before enumerating to avoid millions of
/// useless stdout lines. Default `true` to preserve the existing
/// per-mode behaviour.
static STREAM_RAT_LINES: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(true);

// -------- DFS stats --------

#[derive(Default)]
struct DfsStats {
    closed: u64,
    intersected: u64,
    too_far: u64,
    recursed: u64,
    canonical_skip: u64,
    /// Branches eliminated by the modular reachability prune (set via
    /// `--mod-prune`). Always 0 when the prune is off, so the existing
    /// bench output stays comparable.
    mod_skip: u64,
    /// Branches eliminated by the coordinate-projection prune (set via
    /// passed first, so this reflects the ADDITIONAL pruning power on
    /// top of the modular prune.
    /// Branches eliminated by the closure-key prune (set via
    /// `--closure-key-prune`). Counted only when all prior prunes
    /// passed, so this reflects ADDITIONAL pruning power.
    closure_key_skip: u64,}

impl DfsStats {
    fn total(&self) -> u64 {
        self.closed
            + self.intersected
            + self.too_far
            + self.recursed
            + self.canonical_skip
            + self.mod_skip
            + self.closure_key_skip    }

    fn merge(&mut self, other: &DfsStats) {
        self.closed += other.closed;
        self.intersected += other.intersected;
        self.too_far += other.too_far;
        self.recursed += other.recursed;
        self.canonical_skip += other.canonical_skip;
        self.mod_skip += other.mod_skip;
        self.closure_key_skip += other.closure_key_skip;    }
}

impl std::fmt::Display for DfsStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let t = self.total();
        writeln!(f, "DFS stats ({} total direction attempts):", t)?;
        writeln!(
            f,
            "  canonical_skip:  {:>10} ({:>5.1}%)",
            self.canonical_skip,
            100.0 * self.canonical_skip as f64 / t as f64
        )?;
        writeln!(
            f,
            "  intersected:    {:>10} ({:>5.1}%)",
            self.intersected,
            100.0 * self.intersected as f64 / t as f64
        )?;
        writeln!(
            f,
            "  closed:         {:>10} ({:>5.1}%)",
            self.closed,
            100.0 * self.closed as f64 / t as f64
        )?;
        writeln!(
            f,
            "  recursed:       {:>10} ({:>5.1}%)",
            self.recursed,
            100.0 * self.recursed as f64 / t as f64
        )?;
        writeln!(
            f,
            "  too_far:        {:>10} ({:>5.1}%)",
            self.too_far,
            100.0 * self.too_far as f64 / t as f64
        )?;
        write!(
            f,
            "  mod_skip:       {:>10} ({:>5.1}%)",
            self.mod_skip,
            100.0 * self.mod_skip as f64 / t as f64
        )?;
        writeln!(
            f,
        )?;
        write!(
            f,
            "  closure_key_skip:{:>9} ({:>5.1}%)",
            self.closure_key_skip,
            100.0 * self.closure_key_skip as f64 / t as f64        )
    }
}

// -------- Modular reachability prune --------
//
// Optional, opt-in via `--mod-prune`. For each chosen modulus `m`, BFS
// the set R_r^(m) ⊆ (Z/m)^phi of mod-m displacements reachable by sums
// of EXACTLY r unit vectors, for r in 0..=max_steps. Then during DFS,
// after the canonical / too-far prunes but before `Snake::add`:
//
//   - compute candidate new displacement `new_pt` (already done for
//     the too-far check);
//   - for each modulus, look up `pack(new_pt mod m)` in R_{r-1}^(m);
//   - if any modulus misses, prune (no length-(r-1) suffix can sum
//     to `-new_pt` mod m, so no closure is possible).
//
// The check uses one u64 hash lookup per modulus per DFS node. The
// reachable sets are stored as `FxHashSet<u64>` where each `u64` packs
// the mod-m coordinates as a base-`m` integer (so a `[i64; PHI]`
// digest-vector becomes one `u64` -- saving an alloc per check and
// keeping the hash hot in cache).
//
// We use a global `OnceLock<ModularPrune>` so the existing DFS
// signatures don't grow an `Option<&ModularPrune>` argument; the DFS
// reads from the static and pays nothing when the prune isn't set.

/// Per-(modulus, remaining) reachable-displacement table for the
/// modular prune.
struct ModularPrune {
    /// PHI (= basis dimension) of the underlying ring. Used to pack
    /// candidate displacements into a `u64` key.
    phi: usize,
    /// Active moduli (those whose `m^PHI` fits the budget; smaller
    /// rings get more moduli).
    moduli: Vec<i64>,
    /// `tables[i][r]` = set of packed mod-`moduli[i]` displacements
    /// reachable by sums of exactly `r` unit vectors. A miss in
    /// any modulus is a prune.
    tables: Vec<Vec<rustc_hash::FxHashSet<u64>>>,
}

/// Bundle of all optional DFS prunes. Cloning is `Arc` clones, so
/// snapshotting at DFS entry and passing `&Prunes` through recursion
/// is cheap. Production: built once at `main`-time and stashed in
/// `PRUNES`; the DFS top-level reads it via `snapshot_prunes()`.
/// Tests: construct directly, pass to lower-level DFS functions
/// without going through the global -- avoids the OnceLock-set-once
/// limitation when running different combinations in the same
/// process.
#[derive(Clone, Default)]
struct Prunes {
    mod_prune: Option<std::sync::Arc<ModularPrune>>,
    closure_key_prune: Option<std::sync::Arc<ClosureKeyPrune>>,
}

/// Single global `Prunes` configured by `main` from CLI flags.
/// Tests bypass this; they build a local `Prunes` and pass it to
/// DFS directly via [`run_rat_enum_seqs_with_prunes`].
static PRUNES: Mutex<Option<Prunes>> = Mutex::new(None);

fn snapshot_prunes() -> Prunes {
    PRUNES.lock().unwrap().clone().unwrap_or_default()
}

/// Maximum `m^PHI` cell count we'll keep per modulus. ZZ16/20/24 at
/// m=4 sits right at this budget; ZZ32/60 at m=2 fits comfortably.
/// Larger moduli on higher-phi rings either saturate the prune (and
/// stop being useful) or blow memory; this cutoff trades both off.
const MOD_PRUNE_CELL_BUDGET: u64 = 1 << 16; // 65536

impl ModularPrune {
    /// Build the prune tables from a per-ring unit-vector list. `units`
    /// is `n_units` slices each of length `phi`, holding the integer-
    /// basis coefficients of `Units::unit(d)` for d in
    /// `(-hturn+1)..hturn`.
    ///
    /// `moduli_override`, when present, replaces the default candidate
    /// list (used by `--mod-prune-moduli` for A/B testing). Either way,
    /// candidates are filtered by [`MOD_PRUNE_CELL_BUDGET`].
    fn build(
        units: &[Vec<i64>],
        phi: usize,
        max_steps: usize,
        moduli_override: Option<&[i64]>,
    ) -> Self {
        // Default candidate list: small composites covering the
        // common low-frequency arithmetic constraints (parity, 3-fold,
        // 4-fold, 6-fold). A/B tested against the wider `{2..=16}`
        // set on ZZ8/12/16 at n up to 15: extending the set finds
        // marginally more skip opportunities but the per-node cost
        // of extra hash lookups overshoots the savings -- wall time
        // is equal or worse. The "right" knob is the number of
        // moduli (~4 lookups balance lookup cost vs prune-rate);
        // specific values within {2,3,4,5,6,...} barely matter.
        // Override with `--mod-prune-moduli` to experiment.
        let default_candidates: [i64; 4] = [2, 3, 4, 6];
        let candidates: &[i64] = moduli_override.unwrap_or(&default_candidates);
        let mut moduli: Vec<i64> = Vec::new();
        for &m in candidates {
            if m < 2 {
                continue;
            }
            let cells = (m as u64).checked_pow(phi as u32).unwrap_or(u64::MAX);
            if cells <= MOD_PRUNE_CELL_BUDGET {
                moduli.push(m);
            }
        }

        let mut tables: Vec<Vec<rustc_hash::FxHashSet<u64>>> = Vec::with_capacity(moduli.len());
        for &m in &moduli {
            tables.push(Self::build_one_modulus(units, phi, m, max_steps));
        }

        ModularPrune { phi, moduli, tables }
    }

    fn build_one_modulus(
        units: &[Vec<i64>],
        phi: usize,
        m: i64,
        max_steps: usize,
    ) -> Vec<rustc_hash::FxHashSet<u64>> {
        // We store CUMULATIVE reachability tables: `layers[r]` is the
        // set of mod-`m` displacements reachable by sums of AT MOST
        // `r` unit vectors. The prune check is "could we close in
        // any number of steps from 0..=remaining_after?" -- not
        // "exactly remaining_after", which would falsely reject any
        // rat closing strictly before `max_steps`.
        let total = (m as u64).pow(phi as u32);

        let units_mod: Vec<Vec<i64>> = units
            .iter()
            .map(|u| u.iter().map(|&c| c.rem_euclid(m)).collect())
            .collect();

        let mut layers: Vec<rustc_hash::FxHashSet<u64>> = Vec::with_capacity(max_steps + 1);
        // r=0: only the origin (sum of zero vectors).
        let mut cumulative: rustc_hash::FxHashSet<u64> = rustc_hash::FxHashSet::default();
        cumulative.insert(0u64);
        layers.push(cumulative.clone());

        // Wavefront state for the exact-r BFS step (kept as `Vec<i64>`
        // for simple add-and-reduce; packed only when folding into
        // `cumulative`).
        let mut current_vec: rustc_hash::FxHashSet<Vec<i64>> = rustc_hash::FxHashSet::default();
        current_vec.insert(vec![0i64; phi]);

        let mut saturated = cumulative.len() as u64 == total;
        for _r in 1..=max_steps {
            if saturated {
                layers.push(layers.last().unwrap().clone());
                continue;
            }
            let mut next: rustc_hash::FxHashSet<Vec<i64>> = rustc_hash::FxHashSet::default();
            next.reserve(current_vec.len() * units_mod.len());
            for v in &current_vec {
                for u in &units_mod {
                    let sum: Vec<i64> = v
                        .iter()
                        .zip(u.iter())
                        .map(|(a, b)| (a + b).rem_euclid(m))
                        .collect();
                    next.insert(sum);
                }
            }
            // Fold the new exact-r layer into the cumulative one.
            for v in &next {
                cumulative.insert(pack_coeffs(v, m));
            }
            layers.push(cumulative.clone());
            current_vec = next;
            if cumulative.len() as u64 == total {
                saturated = true;
            }
        }
        layers
    }

    /// Returns `true` if some sum of `remaining` unit vectors equals
    /// the negation of `disp` (mod m) for every active modulus. When
    /// `false`, no length-`remaining` suffix can close, so the caller
    /// prunes. Hot path: one `FxHashSet<u64>::contains` per modulus.
    #[inline]
    fn allows_closure(&self, disp: &[i64], remaining: usize) -> bool {
        // R_r is symmetric under negation (every unit vector u has -u
        // in the direction set, so sums and their negations form the
        // same set), so we can pack `disp` directly instead of `-disp`.
        for (i, &m) in self.moduli.iter().enumerate() {
            let key = pack_coeffs(disp, m);
            let table = match self.tables[i].get(remaining) {
                Some(t) => t,
                None => continue, // remaining beyond max_steps: skip this modulus
            };
            if !table.contains(&key) {
                return false;
            }
        }
        true
    }

    fn cell_counts(&self) -> Vec<u64> {
        self.moduli.iter().map(|&m| (m as u64).pow(self.phi as u32)).collect()
    }
}

/// Pack a coefficient vector mod `m` as a base-`m` integer. Safe up
/// to `m^PHI <= u64::MAX`. Caller guarantees by `MOD_PRUNE_CELL_BUDGET`.
#[inline]
fn pack_coeffs(coeffs: &[i64], m: i64) -> u64 {
    let mut key = 0u64;
    let mut mult = 1u64;
    let m_u = m as u64;
    for &c in coeffs {
        let v = c.rem_euclid(m) as u64;
        key += v * mult;
        mult = mult.wrapping_mul(m_u);
    }
    key
}

/// Build the modular prune for ring `ring` with target `max_steps`,
/// install into the global [`PRUNES`] bundle, and report the moduli
/// kept. Idempotent within a process; replaces any previously
/// installed modular prune. `moduli_override` lets the CLI A/B test
/// by forcing a specific modulus list.
fn install_mod_prune(ring: u8, max_steps: usize, moduli_override: Option<&[i64]>) {
    let (units, phi) = unit_vectors_for_ring(ring);
    let prune = ModularPrune::build(&units, phi, max_steps, moduli_override);
    let cells = prune.cell_counts();
    let summary: Vec<String> = prune
        .moduli
        .iter()
        .zip(cells.iter())
        .map(|(m, c)| format!("m={m} ({c} cells)"))
        .collect();
    eprintln!(
        "mod-prune: ring={ring} phi={phi} max_steps={max_steps}: {}",
        summary.join(", ")
    );
    let mut guard = PRUNES.lock().unwrap();
    let mut current = guard.take().unwrap_or_default();
    current.mod_prune = Some(std::sync::Arc::new(prune));
    *guard = Some(current);
}

/// Per-ring extraction of unit vectors as integer coefficient arrays.
/// Returns ALL `n` unit vectors `unit(0)..unit(n-1)` (not the DFS's
/// `n-1`-direction window): a partial rat's cumulative facing can
/// land on any absolute direction over enough turns, so closure
/// feasibility must be assessed against the full set of possible
/// unit vectors. Used by [`ModularPrune`] to feed its BFS over
/// modular displacements.
fn unit_vectors_for_ring(ring: u8) -> (Vec<Vec<i64>>, usize) {
    fn extract<ZZ: IsRing, const PHI: usize, F>(coeffs_fn: F) -> (Vec<Vec<i64>>, usize)
    where
        F: Fn(&ZZ) -> [i64; PHI],
    {
        let mut out = Vec::new();
        for d in 0..ZZ::turn() {
            let u: ZZ = <ZZ as Units>::unit(d);
            out.push(coeffs_fn(&u).to_vec());
        }
        (out, PHI)
    }
    match ring {
        4 => extract::<ZZ4, 2, _>(|x| x.int_coeffs()),
        8 => extract::<ZZ8, 4, _>(|x| x.int_coeffs()),
        10 => extract::<ZZ10, 4, _>(|x| x.int_coeffs()),
        12 => extract::<ZZ12, 4, _>(|x| x.int_coeffs()),
        16 => extract::<ZZ16, 8, _>(|x| x.int_coeffs()),
        20 => extract::<ZZ20, 8, _>(|x| x.int_coeffs()),
        24 => extract::<ZZ24, 8, _>(|x| x.int_coeffs()),
        32 => extract::<ZZ32, 16, _>(|x| x.int_coeffs()),
        60 => extract::<ZZ60, 16, _>(|x| x.int_coeffs()),
        _ => panic!("invalid ring selected"),
    }
}

// -------- Closure-key prune --------
//
// Pre-pass: DFS-enumerate every simple open snake up to length L,
// store the set of (endpoint, facing) pairs they reach (the "closure
// keys" K_{<=L}). During the main DFS, when remaining_after <= L, the
// candidate prefix can close iff its required suffix's (endpoint,
// facing) is in K_{<=L} -- otherwise no length-(<=remaining_after)
// continuation can complete it.
//
// This is strictly stronger than any modular projection (it sees
// exact lattice + facing info, not just a quotient). The cost is
// memory + pre-pass time + a ring multiplication per hot-path check
// (to compute the target suffix endpoint via
// target = -unit(-facing) * disp).

struct ClosureKeyPrune {
    /// Maximum tabulated suffix length. The prune fires only when
    /// `remaining_after <= max_l`.
    max_l: usize,
    /// `K_{<=L} = {(endpoint coords, facing) : reachable by some
    /// simple open snake of length 0..=L}`. Length 0 (empty snake,
    /// key (0, 0)) is included so already-closed candidates aren't
    /// incorrectly rejected.
    keys: rustc_hash::FxHashSet<(Vec<i64>, i8)>,
}

fn collect_closure_keys_dfs<ZZ: IsRing>(
    snake: &mut Snake<ZZ>,
    max_l: usize,
    keys: &mut rustc_hash::FxHashSet<(Vec<i64>, i8)>,
) {
    if snake.angles().len() >= max_l {
        return;
    }
    let turn = ZZ::turn();
    for direction in ((-ZZ::hturn() + 1)..ZZ::hturn()).rev() {
        if !snake.add(direction) {
            continue;
        }
        // Normalize facing to 0..n-1: `Snake::direction` returns
        // truncating-mod (can be negative for net-CW walks); the
        // hot-path lookup uses `rem_euclid` (always 0..n-1), so we
        // canonicalize here to make the hash keys compatible.
        let facing = snake.direction().rem_euclid(turn);
        keys.insert((snake.offset().int_coeffs_slice().to_vec(), facing));
        collect_closure_keys_dfs::<ZZ>(snake, max_l, keys);
        snake.pop();
    }
}

fn collect_closure_keys<ZZ: IsRing>(max_l: usize) -> rustc_hash::FxHashSet<(Vec<i64>, i8)> {
    let mut keys = rustc_hash::FxHashSet::default();
    // Empty snake's key: closure is already achieved.
    let phi = <ZZ as Units>::unit(0).int_coeffs_slice().len();
    keys.insert((vec![0i64; phi], 0));
    let mut snake: Snake<ZZ> = Snake::new();
    collect_closure_keys_dfs::<ZZ>(&mut snake, max_l, &mut keys);
    keys
}

fn install_closure_key_prune(ring: u8, max_l: usize) {
    let t0 = Instant::now();
    let keys = match ring {
        4 => collect_closure_keys::<ZZ4>(max_l),
        8 => collect_closure_keys::<ZZ8>(max_l),
        10 => collect_closure_keys::<ZZ10>(max_l),
        12 => collect_closure_keys::<ZZ12>(max_l),
        16 => collect_closure_keys::<ZZ16>(max_l),
        20 => collect_closure_keys::<ZZ20>(max_l),
        24 => collect_closure_keys::<ZZ24>(max_l),
        32 => collect_closure_keys::<ZZ32>(max_l),
        60 => collect_closure_keys::<ZZ60>(max_l),
        _ => panic!("invalid ring selected"),
    };
    eprintln!(
        "closure-key-prune: ring={ring} max_l={max_l}: {} distinct keys collected in {:?}",
        keys.len(),
        t0.elapsed(),
    );
    let mut guard = PRUNES.lock().unwrap();
    let mut current = guard.take().unwrap_or_default();
    current.closure_key_prune = Some(std::sync::Arc::new(ClosureKeyPrune { max_l, keys }));
    *guard = Some(current);
}

// --------

/// Pair of canonical-check and output-mapping functions that parameterise
/// the DFS.  Both rotation-canonical and dihedral-canonical enumeration
/// share the same core walk; they differ only in which pair of functions
/// they supply.
///
/// * `is_canonical` — prefix prune applied *before* `Snake::add`.
///   Returns false when the extended walk cannot be the lex-min rotation
///   (or dihedral image) of any closure it could grow into.
/// * `canonicalize` — applied to the chirality-normalised canonical
///   rotation at closure, producing the key inserted into the result
///   `HashSet`.  For rotation-canonical this is the identity; for
///   dihedral-canonical it picks the lex-min over rotations and
///   reversed-rotations.
#[derive(Clone, Copy)]
struct CanonicalOps {
    is_canonical: fn(&[i8], i8) -> bool,
    canonicalize: fn(&[i8]) -> Vec<i8>,
}

/// Bind `prefix ++ [new]` so that index `i` (for `i < prefix.len() + 1`)
/// resolves to the corresponding walk angle. Out-of-range indices are
/// not handled -- callers must keep accesses within `0..d` where
/// `d = prefix.len() + 1`.
fn walk_get(prefix: &[i8], new: i8, i: usize) -> i8 {
    if i < prefix.len() {
        prefix[i]
    } else {
        // Caller must guarantee i == prefix.len(). debug_assert documents
        // the boundary; release-mode behavior is still correct (returns
        // `new`) but a violation indicates a bug in the caller.
        debug_assert_eq!(i, prefix.len(), "walk_get: index past extended prefix");
        new
    }
}

/// For the extended walk `prefix ++ [new]` of length `d`, return
/// `false` if some rotation `k > 0` makes the rotation lex-smaller
/// than the identity within the fully-known comparable region. The
/// caller pairs this with a complementary loop for the dihedral case
/// (see [`is_dihedral_canonical_extended`]).
fn rotation_lex_min_violated(prefix: &[i8], new: i8) -> bool {
    let d = prefix.len() + 1;
    for k in 1..d {
        for i in 0..(d - k) {
            match walk_get(prefix, new, k + i).cmp(&walk_get(prefix, new, i)) {
                std::cmp::Ordering::Less => return true,
                std::cmp::Ordering::Greater => break,
                std::cmp::Ordering::Equal => continue,
            }
        }
    }
    false
}

/// Canonical-rotation prune for the open walk `prefix` virtually
/// extended by one more angle `new`. Lets the caller skip walks
/// that cannot be the lex-min cyclic rotation of any closure they
/// could grow into, *before* paying the cyclotomic intersect cost
/// of `Snake::add`.
///
/// Returns `false` when there's a rotation index `k > 0` such that
/// the rotation `[prefix++new][k..]` is strictly lex-less than
/// `[prefix++new][0..]` within the wrap-free comparable region; that
/// decision is permanent (no future angles can flip it) so the
/// eventual closed polygon's canonical rotation cannot start at
/// position 0 -- pruning is safe.
fn is_canonical_extended(prefix: &[i8], new: i8) -> bool {
    !rotation_lex_min_violated(prefix, new)
}

fn canonical_identity(seq: &[i8]) -> Vec<i8> {
    seq.to_vec()
}

/// Enumerate every simple polygon with boundary length up to
/// `max_steps` over the cyclotomic ring `ZZ`, in canonical-CCW form.
///
/// Depth-first walk: a single `Snake` is mutated in place via
/// [`Snake::add`] / [`Snake::pop`] across the recursion, so the inner
/// loop allocates only when a polygon actually closes (the canonical
/// [`Rat`] needs a fresh vec for hashing). Each polygon's `n` cyclic
/// walks collapse to a single canonical walk via the lex-min
/// rotation prune in [`is_canonical_extended`].
///
/// `step` restricts the DFS to directions that are multiples of `step`.
/// `step = 1` (default) walks every direction. `step = 2` on ZZ20 walks
/// only even-indexed directions, enumerating the ZZ10-equivalent
/// subset; `step = 2` on ZZ24 enumerates the ZZ12 subset; etc.
#[cfg(test)]
fn rat_enum<ZZ: IsRing>(max_steps: usize, step: i8) -> (Vec<Vec<i8>>, DfsStats) {
    rat_enum_with::<ZZ>(
        max_steps,
        step,
        CanonicalOps {
            is_canonical: is_canonical_extended,
            canonicalize: canonical_identity,
        },
        "enumeration",
        "",
        false,
        &Prunes::default(),
    )
}

/// Dihedral-canonical variant of [`rat_enum`].
///
/// Two-stage design (see module-level comment above for rationale):
///
/// 1. [`is_dihedral_canonical_extended`] prunes walk prefixes whose
///    complement rotation is lex-smaller, reducing the search tree
///    (~2.3x speedup at ZZ12 n=10).
///
/// 2. [`dihedral_canonical`] at closure maps the chirality-normalized
///    canonical rotation to the lex-min dihedral form, so both
///    surviving members of each chiral pair hash to the same key.
///
/// Returns the same set as running [`rat_enum`] and quotienting by
/// dihedral equivalence.
#[cfg(test)]
fn rat_enum_dihedral<ZZ: IsRing>(max_steps: usize, step: i8) -> (Vec<Vec<i8>>, DfsStats) {
    rat_enum_with::<ZZ>(
        max_steps,
        step,
        CanonicalOps {
            is_canonical: is_dihedral_canonical_extended,
            canonicalize: dihedral_canonical,
        },
        "dihedral enumeration",
        "dihedral ",
        false,
        &Prunes::default(),
    )
}

fn rat_enum_with<ZZ: IsRing>(
    max_steps: usize,
    step: i8,
    ops: CanonicalOps,
    label: &str,
    prefix: &str,
    paranoid: bool,
    prunes: &Prunes,
) -> (Vec<Vec<i8>>, DfsStats) {
    let mut result: HashSet<Vec<i8>> = HashSet::new();
    let mut snake: Snake<ZZ> = Snake::new();
    let mut stats = DfsStats::default();

    println!("-------- {label} started --------");
    if paranoid {
        println!("paranoid: per-step fresh-snake cross-check enabled");
    }
    rat_enum_step::<ZZ>(
        &mut snake,
        max_steps,
        step,
        &mut result,
        &mut stats,
        ops,
        paranoid,
        prunes,
    );
    println!(
        "-------- {label} completed --------\n{prefix}{} rats found",
        result.len()
    );

    let mut result: Vec<Vec<i8>> = result.into_iter().collect();
    result.sort_by_key(|x| x.len());
    (result, stats)
}

fn rat_enum_step<ZZ: IsRing>(
    snake: &mut Snake<ZZ>,
    max_steps: usize,
    step: i8,
    result: &mut HashSet<Vec<i8>>,
    stats: &mut DfsStats,
    ops: CanonicalOps,
    paranoid: bool,
    prunes: &Prunes,
) {
    let depth = snake.angles().len();
    if depth >= max_steps {
        return;
    }
    let remaining = (max_steps - depth) as i64;

    for direction in ((-ZZ::hturn() + 1)..ZZ::hturn()).rev() {
        if direction.rem_euclid(step) != 0 {
            continue;
        }
        // Canonical prune before the radius and geometry checks so
        // that rejected branches pay nothing.
        if !(ops.is_canonical)(snake.angles(), direction) {
            stats.canonical_skip += 1;
            continue;
        }

        // Early reachability prune: compute the next head position
        // before paying the cyclotomic intersect cost of Snake::add.
        // If the new point is too far from the origin to ever close,
        // skip this direction entirely.
        let new_pt = snake.offset()
            + <ZZ as Units>::unit(snake.direction()) * <ZZ as Units>::unit(direction);
        if !new_pt.is_zero() && !new_pt.within_radius(remaining) {
            stats.too_far += 1;
            continue;
        }

        // Optional modular reachability prune (set via `--mod-prune`).
        // After taking this direction, the snake will be at `new_pt`
        // with `remaining - 1` directions still to add for closure.
        // If no sum of (remaining-1) unit vectors equals `-new_pt`
        // (modulo any active modulus), closure is impossible.
        let remaining_after = (remaining as usize).saturating_sub(1);
        if let Some(mp) = prunes.mod_prune.as_deref() {
            if !mp.allows_closure(new_pt.int_coeffs_slice(), remaining_after) {
                stats.mod_skip += 1;
                continue;
            }
        }
        // Optional closure-key prune (set via `--closure-key-prune`).
        // Only fires when `remaining_after <= max_l` (otherwise the
        // tabulated suffix lengths aren't enough to cover the
        // closing range, and pruning would be unsound).
        if let Some(ck) = prunes.closure_key_prune.as_deref() {
            if remaining_after <= ck.max_l {
                let turn = ZZ::turn();
                let new_facing = (snake.direction() + direction).rem_euclid(turn);
                let neg_facing = (-new_facing).rem_euclid(turn);
                // target suffix endpoint = -unit(-new_facing) * new_pt.
                let target: ZZ = -(<ZZ as Units>::unit(neg_facing) * new_pt);
                let key = (target.int_coeffs_slice().to_vec(), neg_facing);
                if !ck.keys.contains(&key) {
                    stats.closure_key_skip += 1;
                    continue;
                }
            }
        }
        if !snake.add(direction) {
            stats.intersected += 1;
            continue;
        }
        if paranoid {
            // After each successful add(), replay the entire current
            // angle prefix in a fresh Snake. Any disagreement means
            // the stateful incremental check accepts a prefix the
            // from-scratch check rejects -- a Snake bug.
            let angles = snake.angles().to_vec();
            let mut fresh: Snake<ZZ> = Snake::new();
            for (i, &a) in angles.iter().enumerate() {
                assert!(
                    fresh.add(a),
                    "stateful snake accepted full prefix {:?} but fresh snake \
                     rejected angle {} at step {}",
                    &angles,
                    a,
                    i
                );
            }
            assert_eq!(
                fresh.is_closed(),
                snake.is_closed(),
                "fresh snake disagrees on is_closed for {:?}",
                &angles
            );
        }
        if snake.is_closed() {
            stats.closed += 1;
            let r = {
                let tmp = Rat::from_unchecked(snake);
                if tmp.chirality() > 0 {
                    tmp
                } else {
                    tmp.reversed()
                }
                .canonical()
            };
            let seq = (ops.canonicalize)(r.seq());
            if result.insert(seq.clone())
                && STREAM_RAT_LINES.load(std::sync::atomic::Ordering::Relaxed)
            {
                println!("RAT {seq:?}");
            }
        } else {
            stats.recursed += 1;
            rat_enum_step::<ZZ>(snake, max_steps, step, result, stats, ops, paranoid, prunes);
        }
        snake.pop();
    }
}

// -------- dihedral-canonical enumeration --------
//
// Two-stage design: prefix pruning for speed + output dedup for
// correctness.  The prefix prune alone is NOT sufficient to produce
// exactly one representative per dihedral class.  Here is why.
//
// A polygon's dihedral images include all n rotations and all n
// reflections.  For a chiral pair {P, mirror(P)}, both have distinct
// rotation-canonical forms.  To produce dihedral classes, we must
// keep exactly one.
//
// Stage 1 — is_dihedral_canonical_extended:
//
//   Compares the walk prefix against complement rotations (negated
//   angle sequences).  At depth 0 this prunes all positive first
//   angles (-a_0 < a_0), cutting the root branching factor roughly in
//   half.  At deeper depths it prunes additional branches where a
//   complement rotation is lex-smaller.
//
//   This is a HEURISTIC prefix prune, not a dihedral guarantee.  The
//   check guards the raw walk prefix, but the DFS output at closure
//   goes through chirality normalization (.reversed() if chirality <
//   0) then canonical rotation (.canonical()).  These transformations
//   produce a different sequence than what the prefix check was
//   comparing, so both members of some chiral pairs still reach
//   closure as distinct canonical-rotation outputs.
//
// Stage 2 — dihedral_canonical at closure:
//
//   Maps the chirality-normalized canonical rotation to the lex-min
//   over all rotations AND reversed-rotations.  Both canonical
//   rotations of a chiral pair map to the same dihedral-canonical
//   form, so the HashSet deduplicates them.  Without this, the output
//   would contain duplicate chiral pairs — valid polygons, but not
//   dihedral-canonical.
//
// The speedup (~2.3x at ZZ12 n=10) comes from stage 1 reducing the
// search tree.  Stage 2 is cheap (O(n^2) per closed polygon) and
// handles the correctness gap left by stage 1.

/// Heuristic dihedral prefix prune.  Like [`is_canonical_extended`]
/// but also compares the walk prefix against complement rotations
/// (negated angle sequences).  Returns `false` when any complement
/// rotation of the eventual closure has a lex-smaller prefix than the
/// identity walk, indicating the walk is suboptimal under the
/// dihedral group.
///
/// This is sound (never prunes a walk that could produce the
/// dihedral-min output) but incomplete (does not eliminate all
/// chiral duplicates -- see module-level comment above).  At depth 0
/// it prunes all positive first angles (`-a_0 < a_0`), roughly
/// halving the root branching factor.
///
/// # Index bounds
///
/// Complement rotation `k` compares `-a_{k - i}` with `a_i` at
/// positions `i = 0..=k`. Both sides are decided by the known walk
/// only when `0 <= k - i < d` and `0 <= i < d`; for the inner-loop
/// range `i = 0..=k` to stay within the prefix we need `k < d`.
/// Hence the outer loop is `for k in 0..d` (NOT `0..=d`) -- accessing
/// `walk_get(prefix, new, d)` would silently return `new` again for
/// a position whose value is genuinely undetermined, which would
/// over-prune.
fn is_dihedral_canonical_extended(prefix: &[i8], new: i8) -> bool {
    if rotation_lex_min_violated(prefix, new) {
        return false;
    }

    let d = prefix.len() + 1;
    for k in 0..d {
        for i in 0..=k {
            match (-walk_get(prefix, new, k - i)).cmp(&walk_get(prefix, new, i)) {
                std::cmp::Ordering::Less => return false,
                std::cmp::Ordering::Greater => break,
                std::cmp::Ordering::Equal => continue,
            }
        }
    }

    true
}

/// Correctness layer for the dihedral DFS.  Computes the lex-minimum
/// over all cyclic rotations and reversed rotations of `seq`.
///
/// The inputs are already chirality-normalized to CW by the DFS
/// before this function runs. Under that normalization, two mirror
/// images of the same polygon shape (an enantiomer pair) end up as
/// CW walks of the original and the mirror. Algebraically these CW
/// walks differ by SEQUENCE REVERSAL (one polygon walked CW from
/// vertex V0, vs its mirror walked CW from V0' = mirror of V0,
/// gives sequences related by reverse). Plain reversal -- not
/// revcomp -- is therefore the right relation here.
///
/// Both members of an enantiomer pair map to the same value, so the
/// HashSet deduplicates them. Without this step, the output would
/// contain both members of every non-achiral chiral pair.
fn dihedral_canonical(seq: &[i8]) -> Vec<i8> {
    let n = seq.len();
    let mut best: Vec<i8> = seq.to_vec();
    let mut rot: Vec<i8> = seq.to_vec();
    for _ in 0..n {
        if rot < best {
            best = rot.clone();
        }
        let mut rev = rot.clone();
        rev.reverse();
        if rev < best {
            best = rev.clone();
        }
        rot.rotate_left(1);
    }
    best
}

// -------- multi-threaded enumeration --------

/// Pick a DFS splitting depth such that the seed-walk produces
/// roughly `10 * n_threads` work units. With a branching factor of
/// `b` candidate directions per level, depth `d` enumerates at most
/// `b^d` seeds (the canonical-rotation + intersect + reachability
/// prunes knock that down further, but the raw count is the right
/// upper bound for sizing). We invert: `d = ceil(log_b(10 * threads))`.
///
/// `branching` is the per-level branching factor of the DFS:
/// `2 * hturn - 1` for a ZZ ring (the loop walks `(-hturn+1)..hturn`).
/// E.g. ZZ4 -> 3, ZZ12 -> 11, ZZ24 -> 23.
fn splitting_depth(n_threads: usize, branching: usize) -> usize {
    if n_threads <= 1 || branching <= 1 {
        return 0;
    }
    let target = (10 * n_threads) as f64;
    let depth = (target.ln() / (branching as f64).ln()).ceil() as usize;
    depth.max(1)
}

/// Output accumulators for [`collect_seeds`]: incomplete-walk prefixes
/// (seeds for worker threads), already-closed polygons, and DFS stats.
struct SeedGather<'a> {
    seeds: &'a mut Vec<Vec<i8>>,
    closed: &'a mut HashSet<Vec<i8>>,
    stats: &'a mut DfsStats,
}

/// Walk the existing DFS only down to `split_depth` and collect every
/// alive snake state (angle prefix) that survives the canonical-rotation
/// prune, the `Snake::add` self-intersect check, and the reachability
/// heuristic. These are the seeds the worker threads will pick up.
///
/// Polygons that already close at or above `split_depth` are recorded
/// directly into `gather.closed` -- they have no remaining work to delegate.
fn collect_seeds<ZZ: IsRing>(
    snake: &mut Snake<ZZ>,
    max_steps: usize,
    step: i8,
    split_depth: usize,
    gather: &mut SeedGather<'_>,
    ops: CanonicalOps,
    paranoid: bool,
    prunes: &Prunes,
) {
    let depth = snake.angles().len();
    if depth >= max_steps {
        return;
    }
    let remaining = (max_steps - depth) as i64;

    for direction in ((-ZZ::hturn() + 1)..ZZ::hturn()).rev() {
        if direction.rem_euclid(step) != 0 {
            continue;
        }
        if !(ops.is_canonical)(snake.angles(), direction) {
            gather.stats.canonical_skip += 1;
            continue;
        }

        // Early reachability prune (see `rat_enum_step`).
        let new_pt = snake.offset()
            + <ZZ as Units>::unit(snake.direction()) * <ZZ as Units>::unit(direction);
        if !new_pt.is_zero() && !new_pt.within_radius(remaining) {
            gather.stats.too_far += 1;
            continue;
        }

        // Modular reachability prune (see `rat_enum_step`).
        let remaining_after = (remaining as usize).saturating_sub(1);
        if let Some(mp) = prunes.mod_prune.as_deref() {
            if !mp.allows_closure(new_pt.int_coeffs_slice(), remaining_after) {
                gather.stats.mod_skip += 1;
                continue;
            }
        }
        if let Some(ck) = prunes.closure_key_prune.as_deref() {
            if remaining_after <= ck.max_l {
                let turn = ZZ::turn();
                let new_facing = (snake.direction() + direction).rem_euclid(turn);
                let neg_facing = (-new_facing).rem_euclid(turn);
                let target: ZZ = -(<ZZ as Units>::unit(neg_facing) * new_pt);
                let key = (target.int_coeffs_slice().to_vec(), neg_facing);
                if !ck.keys.contains(&key) {
                    gather.stats.closure_key_skip += 1;
                    continue;
                }
            }
        }
        if !snake.add(direction) {
            gather.stats.intersected += 1;
            continue;
        }
        if paranoid {
            let angles = snake.angles().to_vec();
            let mut fresh: Snake<ZZ> = Snake::new();
            for (i, &a) in angles.iter().enumerate() {
                assert!(
                    fresh.add(a),
                    "stateful snake accepted full prefix {:?} but fresh snake \
                     rejected angle {} at step {}",
                    &angles,
                    a,
                    i
                );
            }
            assert_eq!(
                fresh.is_closed(),
                snake.is_closed(),
                "fresh snake disagrees on is_closed for {:?}",
                &angles
            );
        }
        if snake.is_closed() {
            gather.stats.closed += 1;
            let r = {
                let tmp = Rat::from_unchecked(snake);
                if tmp.chirality() > 0 {
                    tmp
                } else {
                    tmp.reversed()
                }
                .canonical()
            };
            let seq = (ops.canonicalize)(r.seq());
            if gather.closed.insert(seq.clone())
                && STREAM_RAT_LINES.load(std::sync::atomic::Ordering::Relaxed)
            {
                println!("RAT {seq:?}");
            }
        } else {
            gather.stats.recursed += 1;
            if snake.angles().len() >= split_depth {
                gather.seeds.push(snake.angles().to_vec());
            } else {
                collect_seeds::<ZZ>(snake, max_steps, step, split_depth, gather, ops, paranoid, prunes);
            }
        }
        snake.pop();
    }
}

/// Parallel variant of [`rat_enum`]: splits the DFS at `split_depth`
/// (selected via [`splitting_depth`]), then hands the resulting alive
/// prefixes out to `n_threads` worker threads via a shared atomic
/// counter. Each worker keeps its own `HashSet` and the main thread
/// merges the per-worker sets at the end.
fn rat_enum_parallel<ZZ: IsRing>(
    max_steps: usize,
    step: i8,
    n_threads: usize,
    ops: CanonicalOps,
    label: &str,
    prefix: &str,
    paranoid: bool,
    prunes: &Prunes,
) -> (Vec<Vec<i8>>, DfsStats) {
    let hm1 = (ZZ::hturn() as usize).saturating_sub(1);
    let branching = 2 * (hm1 / step.max(1) as usize) + 1;
    let split_depth = splitting_depth(n_threads, branching);

    println!("-------- {label} started --------");
    if paranoid {
        println!("paranoid: per-step fresh-snake cross-check enabled");
    }
    println!("parallel: n_threads={n_threads} branching={branching} split_depth={split_depth}");

    let mut closed_main: HashSet<Vec<i8>> = HashSet::new();
    let mut seeds: Vec<Vec<i8>> = Vec::new();
    let mut seed_stats = DfsStats::default();
    {
        let mut snake: Snake<ZZ> = Snake::new();
        let mut gather = SeedGather {
            seeds: &mut seeds,
            closed: &mut closed_main,
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
    }
    println!("parallel: {} seed states collected", seeds.len());

    let (merged, worker_stats) = parallel_drain_seeds::<ZZ>(
        &seeds,
        closed_main,
        seed_stats,
        max_steps,
        step,
        n_threads,
        ops,
        paranoid,
        prunes,
    );

    println!(
        "-------- {label} completed --------\n{prefix}{} rats found",
        merged.len()
    );

    let mut result: Vec<Vec<i8>> = merged.into_iter().collect();
    result.sort_by_key(|x| x.len());
    (result, worker_stats)
}

/// Dispatch a collected list of seed prefixes across `n_threads`
/// worker threads via an atomic counter. Each worker takes seeds one
/// at a time, runs `rat_enum_step` from the seed's snake state, and
/// accumulates canonical sequences into a thread-local HashSet. At
/// the end the locals are folded into `closed_main` (which already
/// contains any polygons that closed during seed collection).
///
/// Used by both `rat_enum_parallel` (whole-tree enumeration, seeds
/// collected from the root) and `enumerate_from_seed_parallel`
/// (single-seed sub-tree enumeration, sub-seeds collected from a
/// given prefix).
#[allow(clippy::too_many_arguments)]
fn parallel_drain_seeds<ZZ: IsRing>(
    seeds: &[Vec<i8>],
    closed_main: HashSet<Vec<i8>>,
    seed_stats: DfsStats,
    max_steps: usize,
    step: i8,
    n_threads: usize,
    ops: CanonicalOps,
    paranoid: bool,
    prunes: &Prunes,
) -> (HashSet<Vec<i8>>, DfsStats) {
    let next_idx = AtomicUsize::new(0);
    let next_ref = &next_idx;

    thread::scope(|s| {
        let mut handles = Vec::with_capacity(n_threads);
        for _ in 0..n_threads {
            handles.push(s.spawn(move || -> (HashSet<Vec<i8>>, DfsStats) {
                let mut local: HashSet<Vec<i8>> = HashSet::new();
                let mut stats = DfsStats::default();
                loop {
                    let i = next_ref.fetch_add(1, Ordering::Relaxed);
                    if i >= seeds.len() {
                        break;
                    }
                    let mut snake: Snake<ZZ> = Snake::from_slice_unsafe(&seeds[i]);
                    rat_enum_step::<ZZ>(
                        &mut snake, max_steps, step, &mut local, &mut stats, ops, paranoid, prunes,
                    );
                }
                (local, stats)
            }));
        }
        let mut merged = closed_main;
        let mut total_stats = seed_stats;
        for h in handles {
            let (local, wstats) = h.join().expect("worker panic");
            merged.extend(local);
            total_stats.merge(&wstats);
        }
        (merged, total_stats)
    })
}

fn polygons<ZZ: IsRing>(rats: Vec<Vec<i8>>) -> Vec<Vec<P64>> {
    rats.into_iter()
        .map(|seq| Rat::<ZZ>::from_slice_unchecked(&seq).to_polyline_f64(Turtle::default()))
        .collect()
}

type EnumResult = (Vec<Vec<i8>>, DfsStats);

fn enumerate_dispatch<ZZ: IsRing>(
    max_steps: usize,
    step: i8,
    n_threads: usize,
    dihedral: bool,
    paranoid: bool,
) -> EnumResult {
    let ops = if dihedral {
        CanonicalOps {
            is_canonical: is_dihedral_canonical_extended,
            canonicalize: dihedral_canonical,
        }
    } else {
        CanonicalOps {
            is_canonical: is_canonical_extended,
            canonicalize: canonical_identity,
        }
    };

    let label = if dihedral {
        "dihedral enumeration"
    } else {
        "enumeration"
    };
    let prefix = if dihedral { "dihedral " } else { "" };

    let prunes = snapshot_prunes();
    if n_threads <= 1 {
        rat_enum_with::<ZZ>(max_steps, step, ops, label, prefix, paranoid, &prunes)
    } else {
        rat_enum_parallel::<ZZ>(max_steps, step, n_threads, ops, label, prefix, paranoid, &prunes)
    }
}

fn run_rat_enum_polylines(
    ring: u8,
    max_steps: usize,
    step: i8,
    n_threads: usize,
    dihedral: bool,
    paranoid: bool,
) -> Vec<Vec<P64>> {
    match ring {
        4 => polygons::<ZZ4>(
            enumerate_dispatch::<ZZ4>(max_steps, step, n_threads, dihedral, paranoid).0,
        ),
        8 => polygons::<ZZ8>(
            enumerate_dispatch::<ZZ8>(max_steps, step, n_threads, dihedral, paranoid).0,
        ),
        10 => polygons::<ZZ10>(
            enumerate_dispatch::<ZZ10>(max_steps, step, n_threads, dihedral, paranoid).0,
        ),
        12 => polygons::<ZZ12>(
            enumerate_dispatch::<ZZ12>(max_steps, step, n_threads, dihedral, paranoid).0,
        ),
        16 => polygons::<ZZ16>(
            enumerate_dispatch::<ZZ16>(max_steps, step, n_threads, dihedral, paranoid).0,
        ),
        20 => polygons::<ZZ20>(
            enumerate_dispatch::<ZZ20>(max_steps, step, n_threads, dihedral, paranoid).0,
        ),
        24 => polygons::<ZZ24>(
            enumerate_dispatch::<ZZ24>(max_steps, step, n_threads, dihedral, paranoid).0,
        ),
        32 => polygons::<ZZ32>(
            enumerate_dispatch::<ZZ32>(max_steps, step, n_threads, dihedral, paranoid).0,
        ),
        60 => polygons::<ZZ60>(
            enumerate_dispatch::<ZZ60>(max_steps, step, n_threads, dihedral, paranoid).0,
        ),
        _ => panic!("invalid ring selected"),
    }
}

/// True iff this canonical sequence equals the canonical form of its
/// mirror image. The mirror polygon's CCW angle sequence is just
/// `reverse(seq)` -- traversing the mirrored vertices in CCW order (in
/// original coordinates) hits the original vertices in reverse order,
/// and each turn keeps its sign because the mirror's orientation flip
/// and the CCW-vs-CW flip cancel. Achiral polygons survive a chirality
/// quotient as a single class.
fn is_achiral(canonical: &[i8]) -> bool {
    let mut reflected: Vec<i8> = canonical.iter().rev().copied().collect();
    let offset = lex_min_rot(&reflected);
    reflected.rotate_left(offset);
    reflected == canonical
}

/// Print per-boundary-length statistics over the enumerated canonical
/// sequences: total count, achiral count, and a histogram of cyclic
/// rotational symmetry orders.
fn print_stats(rats: &[Vec<i8>]) {
    use std::collections::BTreeMap;
    // length -> (total, achiral, rep_factor -> count)
    let mut per_len: BTreeMap<usize, (usize, usize, BTreeMap<usize, usize>)> = BTreeMap::new();
    for seq in rats {
        let n = seq.len();
        let entry = per_len.entry(n).or_default();
        entry.0 += 1;
        if is_achiral(seq) {
            entry.1 += 1;
        }
        let rf = repetition_factor(seq);
        *entry.2.entry(rf).or_insert(0) += 1;
    }
    println!("statistics by boundary length:");
    println!("  length |  total | achiral (%) | rotational symmetry histogram (rep_factor: count)");
    println!("  -------+--------+-------------+--------------------------------------------------");
    for (n, (total, achiral, hist)) in &per_len {
        let pct = if *total > 0 {
            100.0 * (*achiral as f64) / (*total as f64)
        } else {
            0.0
        };
        let mut hist_str = String::new();
        for (rf, cnt) in hist {
            if !hist_str.is_empty() {
                hist_str.push_str("  ");
            }
            hist_str.push_str(&format!("x{rf}={cnt}"));
        }
        println!("  {n:>6} | {total:>6} | {achiral:>5} ({pct:>5.1}%) | {hist_str}");
    }
}

fn run_rat_enum_seqs(
    ring: u8,
    max_steps: usize,
    step: i8,
    n_threads: usize,
    dihedral: bool,
    paranoid: bool,
) -> EnumResult {
    let f: fn(usize, i8, usize, bool, bool) -> EnumResult = match ring {
        4 => enumerate_dispatch::<ZZ4>,
        8 => enumerate_dispatch::<ZZ8>,
        10 => enumerate_dispatch::<ZZ10>,
        12 => enumerate_dispatch::<ZZ12>,
        16 => enumerate_dispatch::<ZZ16>,
        20 => enumerate_dispatch::<ZZ20>,
        24 => enumerate_dispatch::<ZZ24>,
        32 => enumerate_dispatch::<ZZ32>,
        60 => enumerate_dispatch::<ZZ60>,
        _ => panic!("invalid ring selected"),
    };
    f(max_steps, step, n_threads, dihedral, paranoid)
}

// ---- Seed partitioning ----
//
// External-orchestration entry points. `--list-seeds` walks the DFS
// down to `split_depth` and prints (1) "SEED <comma,sep,angles>"
// lines for alive prefixes and (2) "RAT [...]" lines for polygons
// that closed before reaching `split_depth`. `--seed <prefix>`
// resumes the DFS from the given prefix and prints RAT lines for
// what it finds. The two outputs together (one --list-seeds plus
// per-seed jobs) are mergeable via `sort -u` on the RAT lines to
// recover the one-shot result.

/// Build the `CanonicalOps` pair for either rotation-canonical or
/// dihedral-canonical enumeration.
fn make_ops(dihedral: bool) -> CanonicalOps {
    if dihedral {
        CanonicalOps {
            is_canonical: is_dihedral_canonical_extended,
            canonicalize: dihedral_canonical,
        }
    } else {
        CanonicalOps {
            is_canonical: is_canonical_extended,
            canonicalize: canonical_identity,
        }
    }
}

/// Walk the DFS down to `split_depth` and return
/// `(closed, alive_prefixes)`. `closed` contains polygons that
/// closed *during the seed walk* (perimeter <= split_depth). The
/// caller is responsible for handling both: closed ones go straight
/// into the final result; alive prefixes are the work units to be
/// dispatched to per-seed jobs.
fn collect_seed_prefixes<ZZ: IsRing>(
    max_steps: usize,
    step: i8,
    split_depth: usize,
    dihedral: bool,
) -> (HashSet<Vec<i8>>, Vec<Vec<i8>>) {
    let ops = make_ops(dihedral);
    let prunes = snapshot_prunes();
    let mut snake: Snake<ZZ> = Snake::new();
    let mut seeds: Vec<Vec<i8>> = Vec::new();
    let mut closed: HashSet<Vec<i8>> = HashSet::new();
    let mut stats = DfsStats::default();
    {
        let mut gather = SeedGather {
            seeds: &mut seeds,
            closed: &mut closed,
            stats: &mut stats,
        };
        collect_seeds::<ZZ>(
            &mut snake,
            max_steps,
            step,
            split_depth,
            &mut gather,
            ops,
            false,
            &prunes,
        );
    }
    (closed, seeds)
}

/// Resume the DFS from `seed` (a walked angle prefix) and return
/// the unique canonical sequences found below it.
///
/// * `n_threads == 1` (recommended): single-threaded fast path. Runs
///   `rat_enum_step` directly from the seed's snake state. Holds at
///   most one `HashSet` (this seed's share of the final set).
/// * `n_threads > 1`: sub-splits the seed's subtree at depth
///   `seed.len() + splitting_depth(n_threads, branching)` (clamped
///   to `max_steps`) to produce ~10*n_threads work units, then
///   dispatches them across worker threads via `parallel_drain_seeds`.
///
/// # Memory profile and why separate processes are strictly better
///
/// **Single-threaded:** ~1x the final set size at peak (one HashSet,
/// growing with resize transients to ~3x during a bucket-array
/// doubling).
///
/// **Multi-threaded:** ~5x the final set size at peak. Each of
/// `n_threads` workers holds a thread-local HashSet that lives for
/// the whole DFS phase, plus the merge double-buffers as the main
/// `closed_main` absorbs each local. This is the same memory profile
/// as the whole-tree parallel mode.
///
/// **If memory matters and you have multiple seeds to process,
/// running them as separate single-threaded processes (e.g.
/// `xargs -P 16 ... --seed {} --threads 1`) is strictly better:**
/// each process holds only its own share of the final set, and the
/// OS reclaims everything (allocator arenas included) when each
/// process exits. See the file-level docstring for the full
/// accounting. Use `n_threads > 1` only when an orchestrator
/// constrains you to one process at a time.
fn enumerate_from_seed<ZZ: IsRing>(
    max_steps: usize,
    step: i8,
    seed: &[i8],
    n_threads: usize,
    dihedral: bool,
    paranoid: bool,
) -> HashSet<Vec<i8>> {
    let ops = make_ops(dihedral);
    let prunes = snapshot_prunes();

    if n_threads <= 1 {
        // Single-threaded fast path. No sub-split, no worker spawn:
        // just walk the DFS from this seed's snake state directly.
        let mut snake: Snake<ZZ> = Snake::from_slice_unsafe(seed);
        let mut local: HashSet<Vec<i8>> = HashSet::new();
        let mut stats = DfsStats::default();
        rat_enum_step::<ZZ>(
            &mut snake, max_steps, step, &mut local, &mut stats, ops, paranoid, &prunes,
        );
        return local;
    }

    // Multi-threaded: sub-split the seed's subtree to ~10*n_threads
    // work units (capped at `max_steps`), then dispatch.
    let hm1 = (ZZ::hturn() as usize).saturating_sub(1);
    let branching = 2 * (hm1 / step.max(1) as usize) + 1;
    let sub_split = splitting_depth(n_threads, branching);
    let split_depth = (seed.len() + sub_split).min(max_steps);

    let mut closed_main: HashSet<Vec<i8>> = HashSet::new();
    let mut sub_seeds: Vec<Vec<i8>> = Vec::new();
    let mut seed_stats = DfsStats::default();
    {
        let mut snake: Snake<ZZ> = Snake::from_slice_unsafe(seed);
        let mut gather = SeedGather {
            seeds: &mut sub_seeds,
            closed: &mut closed_main,
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
            &prunes,
        );
    }

    let (merged, _) = parallel_drain_seeds::<ZZ>(
        &sub_seeds,
        closed_main,
        seed_stats,
        max_steps,
        step,
        n_threads,
        ops,
        paranoid,
        &prunes,
    );
    merged
}

type SeedListing = (HashSet<Vec<i8>>, Vec<Vec<i8>>);
type CollectSeedFn = fn(usize, i8, usize, bool) -> SeedListing;
type EnumFromSeedFn = fn(usize, i8, &[i8], usize, bool, bool) -> HashSet<Vec<i8>>;

fn dispatch_collect_seed_prefixes(
    ring: u8,
    max_steps: usize,
    step: i8,
    split_depth: usize,
    dihedral: bool,
) -> SeedListing {
    let f: CollectSeedFn = match ring {
        4 => collect_seed_prefixes::<ZZ4>,
        8 => collect_seed_prefixes::<ZZ8>,
        10 => collect_seed_prefixes::<ZZ10>,
        12 => collect_seed_prefixes::<ZZ12>,
        16 => collect_seed_prefixes::<ZZ16>,
        20 => collect_seed_prefixes::<ZZ20>,
        24 => collect_seed_prefixes::<ZZ24>,
        32 => collect_seed_prefixes::<ZZ32>,
        60 => collect_seed_prefixes::<ZZ60>,
        _ => panic!("invalid ring selected"),
    };
    f(max_steps, step, split_depth, dihedral)
}

fn dispatch_enumerate_from_seed(
    ring: u8,
    max_steps: usize,
    step: i8,
    seed: &[i8],
    n_threads: usize,
    dihedral: bool,
    paranoid: bool,
) -> HashSet<Vec<i8>> {
    let f: EnumFromSeedFn = match ring {
        4 => enumerate_from_seed::<ZZ4>,
        8 => enumerate_from_seed::<ZZ8>,
        10 => enumerate_from_seed::<ZZ10>,
        12 => enumerate_from_seed::<ZZ12>,
        16 => enumerate_from_seed::<ZZ16>,
        20 => enumerate_from_seed::<ZZ20>,
        24 => enumerate_from_seed::<ZZ24>,
        32 => enumerate_from_seed::<ZZ32>,
        60 => enumerate_from_seed::<ZZ60>,
        _ => panic!("invalid ring selected"),
    };
    f(max_steps, step, seed, n_threads, dihedral, paranoid)
}

// --------

#[derive(Copy, Clone, Debug, ValueEnum)]
enum Mode {
    /// Enumerate and render output (GIF)
    Render,
    /// Enumerate only and report elapsed time
    Bench,
    /// Walk the DFS down to `split_depth` (default 3) and print alive
    /// prefixes one per line as "SEED a,b,c", plus any polygon that
    /// closed before reaching that depth as "RAT [...]". Use with
    /// `--seed` to dispatch per-prefix jobs externally and merge.
    ListSeeds,
    /// Enumerate, then write a gzipped `RatDafsa` to the path given
    /// by `-o`. The on-disk format is `tilezz-rat-dafsa`; assigned
    /// external indices follow `(length asc, lex asc)` order. The
    /// length-prefix encoding used internally is described in
    /// `RatDafsa::JSON_SCHEMA_DOC`.
    Dafsa,
}

#[derive(Parser, Debug)]
#[command(version, about = "Compute all simple polygons over a ring with a maximal boundary length", long_about = None)]
struct Cli {
    #[arg(short = 'r', long)]
    ring: u8,

    #[arg(short = 'n', long)]
    max_steps: usize,

    #[arg(long, value_enum, default_value_t = Mode::Render)]
    mode: Mode,

    #[arg(
        short = 'o',
        long,
        help = "Output filename: GIF in --mode render, gzipped DAFSA JSON in --mode dafsa"
    )]
    filename: Option<String>,

    #[arg(short, long)]
    verbose: bool,

    /// In Bench mode, write a flamegraph SVG to this path (requires
    /// the `pprof` cargo feature).
    #[arg(long)]
    profile: Option<String>,

    /// In Bench mode, print a per-boundary-length breakdown of total /
    /// achiral / rotational-symmetry histogram after enumeration.
    #[arg(long)]
    stats: bool,

    /// Number of worker threads for the DFS.
    ///
    /// `1` (default) runs the original single-threaded path
    /// unchanged. `0` resolves to `num_cpus::get()`. Any other
    /// positive value is used as-is, capped at `num_cpus::get()`.
    #[arg(long, default_value_t = 1)]
    threads: usize,

    /// Restrict the DFS to directions that are multiples of `step`.
    ///
    /// `1` (default) walks every direction in `(-hturn+1)..hturn`.
    /// `step = 2` on ZZ20 walks only even-indexed directions, which
    /// enumerates exactly the ZZ10-equivalent polygons (b=9 vs 19).
    /// Similarly `step = 2` on ZZ24 -> ZZ12 subset, `step = 3` on
    /// ZZ24 -> ZZ8 subset, etc.
    #[arg(long, default_value_t = 1)]
    step: i8,

    /// Use dihedral-canonical enumeration: outputs one representative
    /// per chiral pair (lex-min over rotations and reflections).
    /// Faster than the default rotation-canonical DFS because
    /// complement-reflection pruning halves the search space at
    /// the root.
    #[arg(long)]
    dihedral: bool,

    /// After every successful `Snake::add`, replay the entire current
    /// angle prefix in a fresh `Snake::new()` and assert that fresh
    /// snake accepts the same prefix and reaches the same
    /// `is_closed` state. Panics on any disagreement -- traps cases
    /// where the stateful incremental check (post-pop grid state)
    /// admits prefixes that a from-scratch check would reject.
    /// Roughly O(n) extra work per step; expect ~2x slowdown.
    #[arg(long)]
    paranoid: bool,

    /// Enable the modular reachability prune. For each modulus
    /// `m in {2..=16}` whose state space `m^phi` fits the cell
    /// budget, precomputes the cumulative reachable mod-m
    /// displacement set and rejects any candidate direction whose
    /// post-extension displacement can't reach 0 mod m in any
    /// number of remaining steps. Cheap per-node (one packed hash
    /// lookup per modulus). See `--mode probe-modular` for the
    /// per-ring saturation curves and `--mod-prune-moduli` for
    /// A/B testing different modulus sets.
    #[arg(long)]
    mod_prune: bool,

    /// Comma-separated list of moduli to use with `--mod-prune`
    /// (default: 2..=16, filtered by ring's state-space budget).
    /// Useful for A/B testing: e.g. `--mod-prune-moduli 2,3,4,6`
    /// reproduces the original hardcoded set.
    #[arg(long, value_delimiter = ',', num_args = 1..)]
    mod_prune_moduli: Option<Vec<i64>>,

    /// Enable the closure-key prune: pre-enumerate every simple
    /// open snake up to length L (= `--closure-key-depth`) and
    /// store their (endpoint, facing) keys. During DFS, when the
    /// remaining-after-this-direction count is <= L, the candidate
    /// is pruned iff its required suffix's closure key isn't in
    /// the precomputed set. Strictly stronger than any modular
    /// projection (uses exact lattice + facing info). More memory
    /// and pre-pass cost; per-node cost is one ring mul + one
    /// hash lookup.
    #[arg(long)]
    closure_key_prune: bool,



    /// Resume the DFS from a specific seed prefix (comma-separated
    /// angles like `-5,3,-4`). When set, the DFS skips seed
    /// collection and resumes from this prefix; output is the same
    /// RAT lines as Bench mode. Honours `--threads`. Used together
    /// with `--mode list-seeds` to dispatch per-prefix jobs
    /// externally.
    #[arg(long, value_delimiter = ',', allow_hyphen_values = true, num_args = 1)]
    seed: Option<Vec<i8>>,

    /// DFS splitting depth for `--mode list-seeds` and for the
    /// parallel (`--threads > 1`) seed-collection step. Defaults to
    /// 3 in list-seeds mode (~50-300 seeds for typical n), or to
    /// the auto-picked value otherwise.
    #[arg(long)]
    split_depth: Option<usize>,
}

/// Resolve the `--threads` CLI value to an actual worker count.
/// `0` -> `num_cpus::get()`; otherwise cap at `num_cpus::get()`.
fn resolve_n_threads(requested: usize) -> usize {
    let max = num_cpus::get().max(1);
    if requested == 0 {
        max
    } else {
        requested.min(max)
    }
}

fn main() {
    let cli = Cli::parse();
    if cli.verbose {
        let mut verbose = VERBOSE.lock().unwrap();
        *verbose = true;
    }

    let n_threads = resolve_n_threads(cli.threads);

    if cli.mod_prune {
        install_mod_prune(
            cli.ring,
            cli.max_steps,
            cli.mod_prune_moduli.as_deref(),
        );
    }
    if cli.closure_key_prune {
        install_closure_key_prune(cli.ring, cli.closure_key_depth);
    }
    // `--mode list-seeds`: walk the DFS down to split_depth, print
    // alive prefixes as `SEED a,b,c` lines and any already-closed
    // polygons as `RAT [...]` lines, then exit. The split_depth
    // defaults to 3 (typically ~50-300 seeds at usable n).
    if matches!(cli.mode, Mode::ListSeeds) {
        let split_depth = cli.split_depth.unwrap_or(3);
        let (closed, seeds) = dispatch_collect_seed_prefixes(
            cli.ring,
            cli.max_steps,
            cli.step,
            split_depth,
            cli.dihedral,
        );
        for prefix in &seeds {
            let s = prefix
                .iter()
                .map(|a| a.to_string())
                .collect::<Vec<_>>()
                .join(",");
            println!("SEED {s}");
        }
        // Polygons that closed before reaching split_depth (typically
        // perimeter-3 only, e.g. the triangle at ZZ12). These would
        // be missed if the orchestrator only ran the SEED jobs.
        for seq in &closed {
            println!("RAT {seq:?}");
        }
        eprintln!(
            "list-seeds: {} alive prefixes at depth {}, {} polygons closed before split",
            seeds.len(),
            split_depth,
            closed.len(),
        );
        return;
    }

    // `--seed <prefix>`: skip seed collection, resume DFS from the
    // given prefix, print RAT lines. With `--threads 1` (default and
    // recommended), this is one work unit; the orchestrator runs
    // multiple `--seed` jobs as separate processes to use more cores.
    // With `--threads N>1`, the seed's subtree is itself sub-split
    // for in-process parallelism -- useful when an orchestrator
    // gives one big process at a time, but per-process memory is
    // ~N x higher. See the file-level docstring for why separate
    // processes are strictly better for memory.
    if let Some(seed) = cli.seed.as_deref() {
        let t0 = Instant::now();
        let rats = dispatch_enumerate_from_seed(
            cli.ring,
            cli.max_steps,
            cli.step,
            seed,
            n_threads,
            cli.dihedral,
            cli.paranoid,
        );
        let dt = t0.elapsed();
        for seq in &rats {
            println!("RAT {seq:?}");
        }
        eprintln!("seed {:?}: {} unique rats in {dt:?}", seed, rats.len());
        return;
    }

    match cli.mode {
        Mode::ListSeeds => unreachable!("handled above"),
        Mode::Dafsa => {
            let Some(filename) = cli.filename.as_deref() else {
                eprintln!("--mode dafsa requires -o <output path>");
                std::process::exit(2);
            };

            // The DFS streams each newly-found rat to stdout by
            // default (for `--seed` / `--mode bench`); in dafsa mode
            // the artifact is the binary file at `-o`, so silence the
            // streaming lines before enumerating.
            STREAM_RAT_LINES.store(false, std::sync::atomic::Ordering::Relaxed);

            let t0 = Instant::now();
            let (rats, _stats) = run_rat_enum_seqs(
                cli.ring,
                cli.max_steps,
                cli.step,
                n_threads,
                cli.dihedral,
                cli.paranoid,
            );
            eprintln!("enumerated {} rats in {:?}", rats.len(), t0.elapsed());

            // `RatDafsa::from_rats` handles the (length, lex) sort,
            // dedup, and length-prefix encoding internally; the
            // resulting external index of each rat is its position in
            // the (length, lex)-ordered sequence.
            let t1 = Instant::now();
            let dafsa = RatDafsa::from_rats(rats.iter().map(|r| r.as_slice()));
            eprintln!(
                "built RatDafsa ({} entries) in {:?}",
                dafsa.len(),
                t1.elapsed()
            );

            let t2 = Instant::now();
            let file = File::create(filename).expect("create output file");
            dafsa
                .write_json_gz(BufWriter::new(file))
                .expect("write gzipped RatDafsa");
            let bytes = std::fs::metadata(filename).map(|m| m.len()).unwrap_or(0);
            eprintln!("wrote {filename} ({bytes} bytes) in {:?}", t2.elapsed());

            if cli.stats {
                print_stats(&rats);
            }
        }
        Mode::Bench => {
            let profile = tilezz::util::profile::ProfileGuard::start(cli.profile.as_deref());

            let t0 = Instant::now();
            let (rats, stats) = run_rat_enum_seqs(
                cli.ring,
                cli.max_steps,
                cli.step,
                n_threads,
                cli.dihedral,
                cli.paranoid,
            );
            let dt = t0.elapsed();

            // Use the result so it can't be trivially optimized away.
            let total_boundary_len: usize = rats.iter().map(|s| s.len()).sum();

            println!(
                "benchmark: ring={} step={} max_steps={} -> {} unique rats (total boundary len={}) in {:?}",
                cli.ring,
                cli.step,
                cli.max_steps,
                rats.len(),
                total_boundary_len,
                dt
            );

            println!("{stats}");

            profile.finish();

            if cli.stats {
                print_stats(&rats);
            }
        }
        Mode::Render => {
            let rats: Vec<Vec<P64>> = run_rat_enum_polylines(
                cli.ring,
                cli.max_steps,
                cli.step,
                n_threads,
                cli.dihedral,
                cli.paranoid,
            );

            let Some(filename) = cli.filename else {
                return;
            };

            // Symmetric square viewport sized to fit every polygon
            // centered at the origin. Each polygon goes through
            // Rat::to_polyline_f64 with the default Turtle (origin,
            // direction 0), so we just need the max |x|/|y| across
            // all polylines.
            let mut max_abs: f64 = 1.0;
            for poly in &rats {
                for (x, y) in poly {
                    max_abs = max_abs.max(x.abs()).max(y.abs());
                }
            }
            let pad = 0.15 * max_abs;
            let r = max_abs + pad;
            let bounds = ((-r, -r), (r, r));

            // Style: filled yellow tile with thin black border,
            // red vertex markers labelled with their index.
            let style = TileStyle::filled(
                Fill::solid(Color::YELLOW.with_alpha(80)),
                Stroke::solid(Color::BLACK, 0.01 * r),
            )
            .with_vertex_marker(MarkerStyle::filled_circle(0.06 * r, Color::RED))
            .with_vertex_labels(TextStyle::new(0.04 * r, Color::WHITE).bold());

            let frames: Vec<Scene> = rats
                .iter()
                .enumerate()
                .map(|(ix, tile)| {
                    println!("rendering frame {}/{}", ix + 1, rats.len());
                    let mut scene = Scene::new().with_background(Color::WHITE);
                    scene.draw_tile(tile, &style);
                    scene
                })
                .collect();

            let w = 500u32;
            let vp = Viewport::square_for(w, bounds, 8);
            let gif_bytes = render_gif(&frames, &vp, 500).expect("render GIF");
            std::fs::write(&filename, gif_bytes).expect("write GIF");
            println!("wrote {filename}");
        }
    }
}

#[cfg(test)]
mod dihedral_tests {
    use super::*;
    use std::collections::HashMap;
    use tilezz::cyclotomic::geometry::intersect;
    use tilezz::cyclotomic::{IsRing, Units, ZZ12, ZZ4, ZZ8};

    /// Independent validator: reconstruct the polygon from `angles` in
    /// exact ring arithmetic and check three properties without any
    /// dependence on `Snake::can_add`:
    ///   (a) the walk closes at the origin;
    ///   (b) the n vertices visited (origin + n-1 intermediates) are
    ///       all pairwise distinct;
    ///   (c) no two non-adjacent boundary edges share any point
    ///       (proper crossing or shared interior).
    fn validate_simple_polygon<ZZ: IsRing>(angles: &[i8]) -> Result<(), String> {
        let n = angles.len();
        if n < 3 {
            return Err(format!("perimeter {n} < 3"));
        }
        // Reconstruct vertices.
        let mut pts: Vec<ZZ> = Vec::with_capacity(n + 1);
        pts.push(ZZ::zero());
        let mut dir: i64 = 0;
        for &a in angles {
            dir = (dir + a as i64).rem_euclid(ZZ::turn() as i64);
            let next = *pts.last().unwrap() + <ZZ as Units>::unit(dir as i8);
            pts.push(next);
        }
        // (a) closes.
        if !pts.last().unwrap().is_zero() {
            return Err(format!("does not close: last={:?}", pts.last().unwrap()));
        }
        pts.pop(); // drop the closing duplicate; now pts has the n cycle vertices

        // (b) distinct vertices. Use a HashMap to find duplicates and
        // report both indices.
        let mut seen: HashMap<ZZ, usize> = HashMap::new();
        for (i, &p) in pts.iter().enumerate() {
            if let Some(&j) = seen.get(&p) {
                return Err(format!("duplicate vertex: pts[{j}] == pts[{i}]"));
            }
            seen.insert(p, i);
        }

        // (c) no two non-adjacent edges share any point. Brute O(n^2)
        // check via `intersect`. Adjacent edges share an endpoint by
        // construction, so we skip cyclically-adjacent pairs.
        for i in 0..n {
            let a1 = pts[i];
            let a2 = pts[(i + 1) % n];
            for j in (i + 2)..n {
                // Skip the cyclic-adjacent pair (last edge vs first edge).
                if i == 0 && j == n - 1 {
                    continue;
                }
                let b1 = pts[j];
                let b2 = pts[(j + 1) % n];
                if intersect(&(a1, a2), &(b1, b2)) {
                    return Err(format!(
                        "edges {i} ({a1:?}->{a2:?}) and {j} ({b1:?}->{b2:?}) intersect/touch"
                    ));
                }
            }
        }

        Ok(())
    }

    /// Per-boundary-length counts the ZZ12 dihedral enumeration must
    /// match: OEIS A316192 "Number of self-avoiding polygons with
    /// perimeter n and sides = 1 ... not counting rotations and
    /// reflections as distinct" (the boundary path may nowhere touch
    /// or intersect itself).
    ///
    /// Reference values, indexed by perimeter `n`:
    ///   n=1..10 -> [0, 0, 1, 3, 4, 22, 69, 418, 2210, 14024]
    ///
    /// This test currently FAILS at n=9 (we report 2217 vs 2210)
    /// and n=10 (14124 vs 14024). The extras are polygons whose
    /// boundary has a T-touch -- a vertex of one edge sitting on the
    /// interior of another edge -- which `intersect_unit_segments`
    /// misses on ZZ12 because the strict CCW check treats
    /// wedge == 0 the same as wedge < 0. See the focused unit test
    /// `test_intersect_zz12_t_touch_endpoint_on_segment` in
    /// `cyclotomic::geometry::tests`.
    #[test]
    fn test_oeis_a316192_zz12() {
        const OEIS: &[(usize, usize)] = &[
            (3, 1),
            (4, 3),
            (5, 4),
            (6, 22),
            (7, 69),
            (8, 418),
            (9, 2210),
            (10, 14024),
        ];
        // Enumerate up to the largest n in the reference. The DFS
        // walks the dihedral-canonical variant; output is bucketed
        // per perimeter length.
        let max_n = OEIS.iter().map(|&(n, _)| n).max().unwrap();
        let (rats, _) = super::rat_enum_dihedral::<ZZ12>(max_n, 1);

        // Bucket per perimeter length.
        let mut by_len: std::collections::BTreeMap<usize, usize> =
            std::collections::BTreeMap::new();
        for seq in &rats {
            *by_len.entry(seq.len()).or_insert(0) += 1;
        }

        let mut mismatches: Vec<(usize, usize, usize)> = Vec::new();
        for &(n, expected) in OEIS {
            let got = by_len.get(&n).copied().unwrap_or(0);
            if got != expected {
                mismatches.push((n, got, expected));
            }
        }
        if !mismatches.is_empty() {
            for (n, got, expected) in &mismatches {
                eprintln!(
                    "n={n}: got {got}, expected (OEIS A316192) {expected}, diff {:+}",
                    *got as i64 - *expected as i64
                );
            }
            panic!(
                "dihedral ZZ12 enumeration differs from OEIS A316192 at {} perimeter length(s)",
                mismatches.len()
            );
        }
    }

    #[test]
    fn test_dihedral_output_polygons_are_simple_and_unique() {
        let (dihedral_rats, _) = super::rat_enum_dihedral::<ZZ12>(9, 1);

        // (i) Encodings are pairwise distinct (the HashSet via which
        // they were collected enforces this; assert defensively).
        let unique: std::collections::HashSet<Vec<i8>> = dihedral_rats.iter().cloned().collect();
        assert_eq!(
            unique.len(),
            dihedral_rats.len(),
            "duplicate sequences in dihedral output"
        );

        // (ii) Every polygon is geometrically simple per the
        // independent validator above. Collect failures so a single
        // run reports every case.
        let mut failures: Vec<(Vec<i8>, String)> = Vec::new();
        for seq in &dihedral_rats {
            if let Err(why) = validate_simple_polygon::<ZZ12>(seq) {
                failures.push((seq.clone(), why));
            }
        }
        if !failures.is_empty() {
            eprintln!("{} polygons failed independent validation:", failures.len());
            for (seq, why) in failures.iter().take(20) {
                eprintln!("  {seq:?} -- {why}");
            }
            panic!("non-simple polygons in dihedral enumeration output");
        }
    }

    /// For a single `(ring, max_steps)` pair, check that the
    /// dihedral DFS output equals the rotation DFS output quotiented
    /// by `dihedral_canonical`. Stage-1 over-pruning shows up as
    /// MISSING dihedral classes; stage-2 under-deduplication shows
    /// up as EXTRA. The body prints both lists before asserting so
    /// any future regression is diagnosable from the test log.
    fn check_quotient_match<ZZ: IsRing>(ring_label: &str, max_steps: usize) {
        let (all_rats, _) = super::rat_enum::<ZZ>(max_steps, 1);
        let mut expected: HashSet<Vec<i8>> = HashSet::new();
        for seq in &all_rats {
            expected.insert(super::dihedral_canonical(seq));
        }
        let (dihedral_rats, _) = super::rat_enum_dihedral::<ZZ>(max_steps, 1);
        let actual: HashSet<Vec<i8>> = dihedral_rats.into_iter().collect();

        if expected != actual {
            for s in expected.difference(&actual) {
                eprintln!("[{ring_label} n={max_steps}] MISSING: {s:?}");
            }
            for s in actual.difference(&expected) {
                eprintln!("[{ring_label} n={max_steps}] EXTRA:   {s:?}");
            }
        }
        assert_eq!(
            expected, actual,
            "{ring_label} n={max_steps}: dihedral DFS != rotation DFS / dihedral",
        );
        eprintln!(
            "[{ring_label} n={max_steps}] OK -- {} dihedral classes from {} rotation-canonical rats",
            actual.len(),
            all_rats.len(),
        );
    }

    /// Sweep multiple `(ring, max_steps)` pairs. ZZ4 and ZZ8 are
    /// included to catch ring-specific over-pruning: with only 3 or
    /// 7 candidate directions per step, the per-rotation comparisons
    /// are more constrained and an off-by-one in the prefix prune is
    /// easier to spot than at ZZ12.
    #[test]
    fn test_dihedral_enum_matches_dfs_quotient() {
        for n in [4, 5, 6, 7, 8, 9] {
            check_quotient_match::<ZZ12>("ZZ12", n);
        }
        for n in [4, 6, 8, 10, 12] {
            check_quotient_match::<ZZ8>("ZZ8", n);
        }
        for n in [4, 6, 8, 10, 12, 14] {
            check_quotient_match::<ZZ4>("ZZ4", n);
        }
    }

    /// Mechanics of the `--mode list-seeds` + per-seed dispatch:
    /// running the DFS in one shot must yield the same set of
    /// canonical sequences as collecting seeds at some `split_depth`,
    /// then enumerating from each seed and union-ing with the
    /// already-closed polygons. Tested across:
    ///   * small `max_steps` (5-7) for tractable runtime;
    ///   * dihedral on AND off (each uses a different `CanonicalOps`
    ///     pair and a different dedup hash);
    ///   * `split_depth` in 1-3 (1: every root-direction branch is a
    ///     seed; 3: small polygons close before split, exercising
    ///     the `closed` half of the merge).
    ///
    /// This is the regression target for the seed-partition API.
    /// Catches any future drift where collect_seeds and the per-seed
    /// runner stop agreeing on what "the same enumeration" means.
    #[test]
    fn test_seed_partitioning_matches_one_shot() {
        for &n in &[5usize, 6, 7] {
            for &dihedral in &[false, true] {
                for &split_depth in &[1usize, 2, 3] {
                    if split_depth >= n {
                        continue;
                    }
                    let (one_shot_seqs, _) = if dihedral {
                        super::rat_enum_dihedral::<ZZ12>(n, 1)
                    } else {
                        super::rat_enum::<ZZ12>(n, 1)
                    };
                    let one_shot: HashSet<Vec<i8>> = one_shot_seqs.into_iter().collect();

                    let (mut seeded, prefixes) =
                        super::collect_seed_prefixes::<ZZ12>(n, 1, split_depth, dihedral);
                    for prefix in &prefixes {
                        // Sweep both branches of the n_threads dispatch
                        // (single-threaded direct path AND sub-split +
                        // parallel_drain_seeds) within the same test --
                        // they must produce the same set.
                        for nthreads in &[1usize, 4] {
                            seeded.extend(super::enumerate_from_seed::<ZZ12>(
                                n, 1, prefix, *nthreads, dihedral, false,
                            ));
                        }
                    }

                    assert_eq!(
                        one_shot,
                        seeded,
                        "n={n} dihedral={dihedral} split_depth={split_depth}: \
                         seed-partitioned != one-shot \
                         ({} prefixes + {} pre-closed)",
                        prefixes.len(),
                        seeded.len() - prefixes.iter().map(|_| 0).sum::<usize>(),
                    );
                }
            }
        }
    }
}
