//! Optional DFS prunes layered on top of the baseline rat enumeration.
//!
//! Two prunes are currently supported, each opt-in:
//!
//! - [`modular`]: precomputed mod-`m` reachable-displacement tables
//!   (cheap, broad applicability, 10-376× depending on ring).
//! - [`closure_key`]: precomputed exact (endpoint, facing) closure
//!   keys for short suffix lengths (more expensive, strictly stronger
//!   than the modular prune, an additional 2-5× on top).
//!
//! Both are gathered into the [`Prunes`] bundle which the DFS reads
//! via `&Prunes` parameter. A single global `Mutex<Option<Prunes>>`
//! is also exposed so the binary can configure prunes once from CLI
//! flags and `snapshot_prunes()` makes a cheap Arc-clone copy at DFS
//! entry. Tests bypass the global and pass `Prunes` directly.

pub mod closure_key;
pub mod modular;
pub mod shadow;
pub mod units;

pub use closure_key::ClosureKeyPrune;
pub use modular::{MOD_PRUNE_CELL_BUDGET, ModularPrune};
pub use shadow::ShadowPrune;

use std::sync::Arc;
use std::sync::Mutex;

/// Bundle of all optional DFS prunes. Cloning is cheap (Arc fields),
/// so snapshotting at DFS entry and passing `&Prunes` through
/// recursion is essentially free. Production: built once at
/// `main`-time and stashed in `PRUNES`; the DFS top-level reads it
/// via `snapshot_prunes()`. Tests: construct directly, pass to
/// lower-level DFS functions without going through the global --
/// avoids the OnceLock-set-once limitation when running different
/// combinations in the same process.
#[derive(Clone, Default)]
pub struct Prunes {
    pub mod_prune: Option<Arc<ModularPrune>>,
    pub closure_key_prune: Option<Arc<ClosureKeyPrune>>,
    pub shadow_prune: Option<Arc<ShadowPrune>>,
}

/// Single global `Prunes` configured by `main` from CLI flags.
/// Tests bypass this; they build a local `Prunes` and pass it to
/// DFS directly.
static PRUNES: Mutex<Option<Prunes>> = Mutex::new(None);

/// Make a cheap (Arc-clone) snapshot of the currently installed prunes.
/// The DFS top-level calls this once and passes the snapshot through
/// recursion to avoid lock contention on the hot path.
pub fn snapshot_prunes() -> Prunes {
    PRUNES.lock().unwrap().clone().unwrap_or_default()
}

/// Build the modular prune for ring `ring` with target `max_steps`,
/// install into the global `PRUNES` bundle, and report the moduli
/// kept. Idempotent within a process; replaces any previously
/// installed modular prune. `moduli_override` lets the CLI A/B test
/// by forcing a specific modulus list.
pub fn install_mod_prune(ring: u8, max_steps: usize, moduli_override: Option<&[i64]>) {
    let (units, phi) = units::unit_vectors_for_ring(ring);
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
    current.mod_prune = Some(Arc::new(prune));
    *guard = Some(current);
}

/// Build the closure-key prune for ring `ring` up to suffix length
/// `max_l`, install into the global `PRUNES` bundle, report the
/// build time and key count. Idempotent.
pub fn install_closure_key_prune(ring: u8, max_l: usize) {
    let t0 = std::time::Instant::now();
    let keys = closure_key::collect_closure_keys_for_ring(ring, max_l);
    eprintln!(
        "closure-key-prune: ring={ring} max_l={max_l}: {} distinct keys collected in {:?}",
        keys.len(),
        t0.elapsed(),
    );
    let prune = ClosureKeyPrune { max_l, keys };
    let mut guard = PRUNES.lock().unwrap();
    let mut current = guard.take().unwrap_or_default();
    current.closure_key_prune = Some(Arc::new(prune));
    *guard = Some(current);
}

/// Build the shadow-radius prune for ring `ring`, install into the
/// global `PRUNES` bundle, and report the shadow places. Idempotent.
pub fn install_shadow_prune(ring: u8) {
    let prune = ShadowPrune::for_ring(ring);
    eprintln!(
        "shadow-prune: ring={ring}: {} shadow place(s), exponents {:?}",
        prune.exps.len(),
        prune.exps,
    );
    let mut guard = PRUNES.lock().unwrap();
    let mut current = guard.take().unwrap_or_default();
    current.shadow_prune = Some(Arc::new(prune));
    *guard = Some(current);
}
