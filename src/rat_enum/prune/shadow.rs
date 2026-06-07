//! Shadow-radius reachability prune: `within_radius` extended to the
//! non-physical archimedean places of `Z[zeta_n]`.
//!
//! The baseline `within_radius` prune enforces that the head can still
//! reach the origin in the *physical* embedding: `|sigma_1(head)| <=
//! remaining`, since each remaining unit step has modulus 1 there. But
//! a closing walk returns to 0 in **every** embedding simultaneously,
//! and every unit step is a root of unity (modulus 1) in every
//! embedding -- so the identical bound holds in each conjugate
//! ("shadow") place:
//!
//! ```text
//!     to close,  |sigma_g(head)| <= remaining   for every embedding g.
//! ```
//!
//! `within_radius` checks one of the `phi/2` archimedean places; this
//! prune checks the other `phi/2 - 1`. One Galois exponent represents
//! each non-physical conjugate pair: the units `2 <= g <= n/2` coprime
//! to `n` (ZZ12 -> {5}, ZZ24 -> {5,7,11}, ZZ32 -> {3,5,7,9,11,13,15}).
//!
//! This is the cyclotomic spacing bound used as a prune: `|sigma_1| *
//! |sigma_g| >= 1` (the field norm is a nonzero integer), so a head
//! physically near the origin necessarily has a large shadow -- which
//! `within_radius` cannot see but this prune kills.
//!
//! Relation to [`super::modular`]: this is the **same** prune as the
//! modular one -- the necessary condition "can the remaining unit-step
//! budget cancel the head" -- just checked at the other places of the
//! field: archimedean here, finite (mod m) there. They are two halves
//! of one reachability prune and share the `--reachability-prune` flag. They are
//! complementary, not redundant: the modular tables are `m^phi` cells
//! (budget-capped, so they shrink as `phi` grows), while this check is
//! `O(phi)` arithmetic with no precomputation -- so it carries the most
//! weight on exactly the high-`phi` rings where the modular half is
//! reduced to `m = 2`.

// `conj`, `re_sign`, and `zero` are all reachable through `IsRing`'s
// supertraits (Conj / ReImSign / IntRing), so no extra `use` is needed.
use crate::cyclotomic::{IsRing, Units};

/// The Galois exponents of the non-physical conjugate places of a ring,
/// reused per node by [`ShadowPrune::allows_closure`].
pub struct ShadowPrune {
    /// Ring order `n` (so that `unit(k)` indices reduce mod `n`).
    pub n: i64,
    /// One representative `g` per non-physical conjugate pair: the units
    /// `2 <= g <= n/2` coprime to `n`. Empty for the discrete rings
    /// (`ZZ4`, `ZZ6`) where there is no shadow place -- then the prune
    /// is a no-op.
    pub exps: Vec<i64>,
}

impl ShadowPrune {
    /// Build the shadow places for ring `n`.
    pub fn for_ring(ring: u8) -> Self {
        let n = ring as i64;
        let exps = (2..=n / 2).filter(|&g| gcd(g, n) == 1).collect();
        ShadowPrune { n, exps }
    }

    /// `true` if `point` is still within closing range `remaining` in
    /// every shadow place; `false` (prune) if some `|sigma_g(point)| >
    /// remaining`, since the remaining unit steps -- each of modulus 1
    /// in that place -- cannot cancel it back to 0.
    ///
    /// Uses the same radius convention as [`WithinRadius`] (the physical
    /// place), so this is its exact archimedean analogue. Exact, no
    /// f64: `|sigma_g(point)|^2 = sigma_g(point * conj(point))` (the
    /// physical norm-square `point*conj` lives in the totally-real
    /// subring, and `zeta -> zeta^g` carries it to the squared modulus
    /// in the g-th place), compared to `remaining^2` via `re_sign`.
    ///
    /// [`WithinRadius`]: crate::cyclotomic::WithinRadius
    #[inline]
    pub fn allows_closure<ZZ: IsRing>(&self, point: &ZZ, remaining: i64) -> bool {
        if self.exps.is_empty() {
            return true;
        }
        let nsq = *point * point.conj();
        let coeffs = nsq.int_coeffs_slice();
        let r_sq = <ZZ as From<i64>>::from(remaining * remaining);
        for &g in &self.exps {
            // sigma_g(nsq) = sum_j coeffs[j] * unit((g*j) mod n)
            let mut s = ZZ::zero();
            for (j, &c) in coeffs.iter().enumerate() {
                if c != 0 {
                    let k = (g * j as i64).rem_euclid(self.n) as i8;
                    s = s + <ZZ as Units>::unit(k).scale(c);
                }
            }
            // r_sq - |sigma_g(point)|^2 < 0  <=>  |sigma_g(point)| > remaining
            if (r_sq - s).re_sign() < 0 {
                return false;
            }
        }
        true
    }
}

fn gcd(mut a: i64, mut b: i64) -> i64 {
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a.abs()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic::{ZZ12, ZZ24, ZZ32};

    #[test]
    fn shadow_exps_are_the_nonphysical_conjugate_reps() {
        assert_eq!(ShadowPrune::for_ring(12).exps, vec![5]);
        assert_eq!(ShadowPrune::for_ring(24).exps, vec![5, 7, 11]);
        assert_eq!(ShadowPrune::for_ring(32).exps, vec![3, 5, 7, 9, 11, 13, 15]);
        // discrete rings have no shadow place -> no-op prune
        assert!(ShadowPrune::for_ring(4).exps.is_empty());
        assert!(ShadowPrune::for_ring(6).exps.is_empty());
    }

    /// For every reachable point and every radius, the prune must agree
    /// with the exact statement "|sigma_g(point)| <= radius for all g":
    /// it must never reject a point that genuinely closes within radius.
    /// We check it against a direct f64 evaluation of the g-th embedding.
    fn agrees_with_direct_embedding<ZZ: IsRing>(ring: u8) {
        let sp = ShadowPrune::for_ring(ring);
        let n = ring as i64;
        // a spread of reachable points: short sums of unit steps
        let dirs: Vec<ZZ> = (0..ring as i8).map(|d| <ZZ as Units>::unit(d)).collect();
        let mut pts: Vec<ZZ> = vec![ZZ::zero()];
        for &a in &dirs {
            for &b in &dirs {
                pts.push(a + b);
                pts.push(a + b + <ZZ as Units>::unit(0));
            }
        }
        for p in &pts {
            let coeffs = p.int_coeffs_slice();
            for radius in 0..6i64 {
                // direct: max over shadow places of |sigma_g(p)|, vs radius
                let mut ok = true;
                for &g in &sp.exps {
                    let (mut re, mut im) = (0.0f64, 0.0f64);
                    for (j, &c) in coeffs.iter().enumerate() {
                        let ang =
                            std::f64::consts::TAU * (g * j as i64).rem_euclid(n) as f64 / n as f64;
                        re += c as f64 * ang.cos();
                        im += c as f64 * ang.sin();
                    }
                    if (re * re + im * im).sqrt() > radius as f64 + 1e-9 {
                        ok = false;
                    }
                }
                assert_eq!(
                    sp.allows_closure(p, radius),
                    ok,
                    "ZZ{ring}: shadow prune disagrees with direct embedding for \
                     coeffs={coeffs:?} radius={radius}"
                );
            }
        }
    }

    #[test]
    fn allows_closure_matches_direct_embedding_zz12() {
        agrees_with_direct_embedding::<ZZ12>(12);
    }

    #[test]
    fn allows_closure_matches_direct_embedding_zz24() {
        agrees_with_direct_embedding::<ZZ24>(24);
    }

    #[test]
    fn allows_closure_matches_direct_embedding_zz32() {
        agrees_with_direct_embedding::<ZZ32>(32);
    }
}
