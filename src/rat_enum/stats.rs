//! Per-DFS counters reported by `rat_enum_step` and its callers.
//!
//! Each surviving DFS branch falls into exactly one of:
//! - `closed`: snake reached the origin and was recorded as a rat.
//! - `intersected`: `Snake::add` rejected the candidate (self-intersection).
//! - `too_far`: Euclidean reachability prune fired (snake too far from origin).
//! - `shadow_skip`: shadow-radius reachability prune fired (too far in a
//!   non-physical archimedean place; archimedean half of `--reachability-prune`).
//! - `recursed`: snake extended; DFS descended further.
//! - `canonical_skip`: canonical-rotation prune fired (lex-min violation).
//! - `modular_skip`: modular reachability prune fired (set via `--reachability-prune`).
//! - `closure_table_skip`: closure-table prune fired (set via `--closure-table-prune`).
//!
//! These categories partition every direction attempted; their sum
//! is `total()`.

#[derive(Default)]
pub struct DfsStats {
    pub closed: u64,
    pub intersected: u64,
    pub too_far: u64,
    /// Branches eliminated by the shadow-radius prune (reachability in
    /// the non-physical archimedean places; archimedean half of
    /// `--reachability-prune`). Counted only when `within_radius` (physical
    /// place) already passed, so it reflects ADDITIONAL pruning beyond
    /// `too_far`.
    pub shadow_skip: u64,
    pub recursed: u64,
    pub canonical_skip: u64,
    /// Branches eliminated by the modular reachability prune (set
    /// via `--reachability-prune`). Always 0 when the prune is off, so the
    /// existing bench output stays comparable.
    pub modular_skip: u64,
    /// Branches eliminated by the closure-table prune (set via
    /// `--closure-table-prune`). Counted only when all prior prunes
    /// passed, so this reflects ADDITIONAL pruning power.
    pub closure_table_skip: u64,
}

impl DfsStats {
    pub fn total(&self) -> u64 {
        self.closed
            + self.intersected
            + self.too_far
            + self.shadow_skip
            + self.recursed
            + self.canonical_skip
            + self.modular_skip
            + self.closure_table_skip
    }

    pub fn merge(&mut self, other: &DfsStats) {
        self.closed += other.closed;
        self.intersected += other.intersected;
        self.too_far += other.too_far;
        self.shadow_skip += other.shadow_skip;
        self.recursed += other.recursed;
        self.canonical_skip += other.canonical_skip;
        self.modular_skip += other.modular_skip;
        self.closure_table_skip += other.closure_table_skip;
    }
}

impl std::fmt::Display for DfsStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let t = self.total();
        let pct = |n: u64| 100.0 * n as f64 / t as f64;
        writeln!(f, "DFS stats ({} total direction attempts):", t)?;
        writeln!(
            f,
            "  canonical_skip:   {:>10} ({:>5.1}%)",
            self.canonical_skip,
            pct(self.canonical_skip)
        )?;
        writeln!(
            f,
            "  intersected:      {:>10} ({:>5.1}%)",
            self.intersected,
            pct(self.intersected)
        )?;
        writeln!(
            f,
            "  closed:           {:>10} ({:>5.1}%)",
            self.closed,
            pct(self.closed)
        )?;
        writeln!(
            f,
            "  recursed:         {:>10} ({:>5.1}%)",
            self.recursed,
            pct(self.recursed)
        )?;
        writeln!(
            f,
            "  too_far:          {:>10} ({:>5.1}%)",
            self.too_far,
            pct(self.too_far)
        )?;
        writeln!(
            f,
            "  shadow_skip:      {:>10} ({:>5.1}%)",
            self.shadow_skip,
            pct(self.shadow_skip)
        )?;
        writeln!(
            f,
            "  modular_skip:         {:>10} ({:>5.1}%)",
            self.modular_skip,
            pct(self.modular_skip)
        )?;
        write!(
            f,
            "  closure_table_skip: {:>10} ({:>5.1}%)",
            self.closure_table_skip,
            pct(self.closure_table_skip)
        )
    }
}
