//! Modular reachability prune.
//!
//! Opt-in via `--mod-prune`. For each chosen modulus `m`, BFS the
//! set `R_{<=r}^(m) ⊆ (Z/m)^phi` of mod-`m` displacements reachable
//! by sums of AT MOST `r` unit vectors, for `r` in `0..=max_steps`.
//! Then during DFS, after the canonical / too-far prunes but before
//! `Snake::add`:
//!
//!   - compute the candidate new displacement `new_pt` (already done
//!     for the too-far check);
//!   - for each modulus, look up `pack(new_pt mod m)` in
//!     `R_{<=remaining_after}^(m)`;
//!   - if any modulus misses, prune (no continuation of length
//!     `<= remaining_after` can sum to `-new_pt` mod m, so closure
//!     is impossible).
//!
//! Hot path: one `FxHashSet<u64>::contains` per modulus. The
//! reachable sets are stored with each `u64` packing the mod-`m`
//! coordinates as a base-`m` integer (a `[i64; PHI]` digest-vector
//! becomes one `u64`), saving an alloc per check and keeping the
//! hash hot in cache.

/// Maximum `m^PHI` cell count we'll keep per modulus. ZZ16/20/24 at
/// m=4 sits right at this budget; ZZ32/60 at m=2 fits comfortably.
/// Larger moduli on higher-phi rings either saturate the prune (and
/// stop being useful) or blow memory; this cutoff trades both off.
pub const MOD_PRUNE_CELL_BUDGET: u64 = 1 << 16; // 65536

/// Per-(modulus, remaining) reachable-displacement table.
pub struct ModularPrune {
    /// PHI (= basis dimension) of the underlying ring. Used to pack
    /// candidate displacements into a `u64` key.
    pub phi: usize,
    /// Active moduli (those whose `m^PHI` fits the budget; smaller
    /// rings get more moduli).
    pub moduli: Vec<i64>,
    /// `tables[i][r]` = set of packed mod-`moduli[i]` displacements
    /// reachable by sums of AT MOST `r` unit vectors. A miss in
    /// any modulus is a prune.
    pub tables: Vec<Vec<rustc_hash::FxHashSet<u64>>>,
}

impl ModularPrune {
    /// Build the prune tables from a per-ring unit-vector list.
    /// `units` is `n_units` slices each of length `phi`, holding the
    /// integer-basis coefficients of `Units::unit(d)` for
    /// `d in 0..n`.
    ///
    /// `moduli_override`, when present, replaces the default
    /// candidate list (used by `--mod-prune-moduli` for A/B testing).
    /// Either way, candidates are filtered by [`MOD_PRUNE_CELL_BUDGET`].
    pub fn build(
        units: &[Vec<i64>],
        phi: usize,
        max_steps: usize,
        moduli_override: Option<&[i64]>,
    ) -> Self {
        // Default candidate list: small composites covering the
        // common low-frequency arithmetic constraints (parity,
        // 3-fold, 4-fold, 6-fold). A/B tested against the wider
        // `{2..=16}` set on ZZ8/12/16 at n up to 15: extending the
        // set finds marginally more skip opportunities but the
        // per-node cost of extra hash lookups overshoots the
        // savings -- wall time is equal or worse. The "right" knob
        // is the number of moduli (~4 lookups balance lookup cost
        // vs prune-rate); specific values within {2,3,4,5,6,...}
        // barely matter. Override with `--mod-prune-moduli` to
        // experiment.
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

        ModularPrune {
            phi,
            moduli,
            tables,
        }
    }

    fn build_one_modulus(
        units: &[Vec<i64>],
        phi: usize,
        m: i64,
        max_steps: usize,
    ) -> Vec<rustc_hash::FxHashSet<u64>> {
        // CUMULATIVE reachability tables: `layers[r]` is the set of
        // mod-`m` displacements reachable by sums of AT MOST `r`
        // unit vectors. The prune check is "could we close in any
        // number of steps from 0..=remaining_after?" -- NOT "exactly
        // remaining_after", which would falsely reject any rat
        // closing strictly before `max_steps`.
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

        // Wavefront state for the exact-r BFS step (kept as
        // `Vec<i64>` for simple add-and-reduce; packed only when
        // folding into `cumulative`).
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

    /// Returns `true` if some sum of at most `remaining` unit
    /// vectors equals the negation of `disp` (mod m) for every
    /// active modulus. When `false`, no length-`<=remaining` suffix
    /// can close, so the caller prunes. Hot path: one
    /// `FxHashSet<u64>::contains` per modulus.
    #[inline]
    pub fn allows_closure(&self, disp: &[i64], remaining: usize) -> bool {
        // R_r is symmetric under negation (every unit vector u has
        // -u in the direction set, so sums and their negations form
        // the same set), so we can pack `disp` directly instead of
        // `-disp`.
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

    pub fn cell_counts(&self) -> Vec<u64> {
        self.moduli
            .iter()
            .map(|&m| (m as u64).pow(self.phi as u32))
            .collect()
    }
}

/// Pack a coefficient vector mod `m` as a base-`m` integer. Safe up
/// to `m^PHI <= u64::MAX`. Caller guarantees this via
/// [`MOD_PRUNE_CELL_BUDGET`].
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
