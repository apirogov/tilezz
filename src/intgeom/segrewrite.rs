//! View over the rewrite rules computed inside [`SegSeqBFS`].
//!
//! The BFS computes everything (LHS, MatchMeta, RHS, witness) per
//! rule; this module just provides:
//!   - lookups by LHS (`by_lhs`, `rules_by_lhs`),
//!   - per-rule access by id (passthrough to the BFS),
//!   - segment classification: terminals and the cursed attractor.

use rustc_hash::{FxHashMap, FxHashSet};

use crate::cyclotomic::{IsComplex, IsRingOrField, Units};
use crate::intgeom::patch::PatchMatch;
use crate::intgeom::segseqbfs::{lhs_seg_range, MatchMeta, RuleRecord, SegSeqBFS};

#[derive(Debug)]
pub struct SegRewriteTable {
    /// LHS sequence → MatchMeta → RHS (= the classic rewrite map).
    pub by_lhs: FxHashMap<Vec<usize>, FxHashMap<MatchMeta, [usize; 3]>>,
    /// LHS sequence → list of rule_ids into the BFS's `rules`.
    pub rules_by_lhs: FxHashMap<Vec<usize>, Vec<usize>>,
}

impl SegRewriteTable {
    pub fn build_from<T: IsComplex + IsRingOrField + Units>(
        seq_bfs: &SegSeqBFS<T>,
    ) -> Self {
        Self::build_with_stats(seq_bfs)
    }

    pub fn build_with_stats<T: IsComplex + IsRingOrField + Units>(
        seq_bfs: &SegSeqBFS<T>,
    ) -> Self {
        let mut by_lhs: FxHashMap<Vec<usize>, FxHashMap<MatchMeta, [usize; 3]>> =
            FxHashMap::default();
        let mut rules_by_lhs: FxHashMap<Vec<usize>, Vec<usize>> =
            FxHashMap::default();
        for (rule_id, rule) in seq_bfs.rules().iter().enumerate() {
            by_lhs
                .entry(rule.lhs.clone())
                .or_default()
                .insert(rule.meta, rule.rhs);
            rules_by_lhs
                .entry(rule.lhs.clone())
                .or_default()
                .push(rule_id);
        }
        SegRewriteTable { by_lhs, rules_by_lhs }
    }

    /// Distinct LHS sequences.
    pub fn num_lhs(&self) -> usize {
        self.by_lhs.len()
    }

    /// Total distinct (LHS, MatchMeta, RHS) rules.
    pub fn num_rules(&self) -> usize {
        self.by_lhs.values().map(|m| m.len()).sum()
    }

    /// All rewrites firing on a given LHS, if any.
    pub fn rules_for_lhs(&self, lhs: &[usize]) -> Option<&FxHashMap<MatchMeta, [usize; 3]>> {
        self.by_lhs.get(lhs)
    }

    /// Rule ids (indices into `seq_bfs.rules()`) for a given LHS.
    pub fn rule_ids_for_lhs(&self, lhs: &[usize]) -> Option<&[usize]> {
        self.rules_by_lhs.get(lhs).map(|v| v.as_slice())
    }
}

/// Classification analysis over the rule set.
///
/// The "cursed" concept is defined as a *per-rule* least fixed
/// point. Intuitively: applying a non-cursed rule must land you in a
/// state where you can keep rewriting non-cursedly.
///
/// 1. **Terminal segment**: a seg id that appears in some rule's RHS
///    but never in any rule's LHS. Once introduced, it cannot be
///    rewritten away.
///
/// 2. **Cursed rule** (least fixed point):
///    - Base case: any rule whose RHS contains a terminal segment.
///    - Inductive case: a rule R is cursed iff its own trial T_R
///      (= the un-normalized patch you get by applying R to its
///      canonical witness) has NO non-cursed LHS covering R.rhs at
///      the RHS position. Concretely: every LHS substring of T_R's
///      boundary that satisfies `q ≤ rhs_idx` AND `q + L > rhs_idx`
///      (= would rewrite at least r's first segment) is itself
///      cursed.
///    - This is checked iteratively until the cursed set stabilizes.
///
/// 3. **Cursed LHS**: every rule with this LHS is cursed.
///
/// 4. **Dead RHS**: every rule producing this RHS is cursed (= no
///    non-cursed rule can produce it). Derived as a final pass over
///    `rules.rhs` partition by cursed-status.
///
/// 5. **Nonlocal rule** (independent of cursedness): rule succeeds on
///    its canonical witness but `add_tile` fails on at least one
///    other patch containing the same LHS. Context beyond L matters
///    for applicability. Computed only over non-cursed rules.
///
/// `context_sensitive_rhs` is left empty under the per-rule definition
/// (it's not a natural derived concept) and may be reintroduced later
/// if the use case clarifies.
#[derive(Debug)]
pub struct RuleAnalysis {
    pub terminals: FxHashSet<usize>,
    pub dead_rhs: FxHashSet<[usize; 3]>,
    pub context_sensitive_rhs: FxHashSet<[usize; 3]>,
    pub cursed_rules: FxHashSet<usize>,
    pub cursed_lhs: FxHashSet<Vec<usize>>,
    pub nonlocal_rules: FxHashSet<usize>,
}

impl RuleAnalysis {
    pub fn build<T: IsComplex + IsRingOrField + Units>(
        seq_bfs: &SegSeqBFS<T>,
        table: &SegRewriteTable,
    ) -> Self {
        let rules: &[RuleRecord] = seq_bfs.rules();
        // Collect seg ids that appear in any LHS and any RHS.
        let mut lhs_segs: FxHashSet<usize> = FxHashSet::default();
        let mut rhs_segs: FxHashSet<usize> = FxHashSet::default();
        for r in rules {
            for &s in &r.lhs {
                lhs_segs.insert(s);
            }
            for &s in &r.rhs {
                rhs_segs.insert(s);
            }
        }
        let terminals: FxHashSet<usize> =
            rhs_segs.difference(&lhs_segs).copied().collect();

        // For each rule, collect the set of LHS-ids that cover its
        // RHS from the left on the rule's own trial. This is the
        // per-rule "alive evidence": a rule is alive iff at least one
        // of its covering LHSes is *non-cursed*. We compute the
        // evidence once, then iterate the cursed set to a fixed point.
        let lhs_index: FxHashMap<Vec<usize>, usize> = table
            .rules_by_lhs
            .keys()
            .enumerate()
            .map(|(i, k)| (k.clone(), i))
            .collect();
        let n_lhses = lhs_index.len();
        let max_lhs_len = lhs_index.keys().map(|l| l.len()).max().unwrap_or(0);

        // Closure: gather all LHS-ids in our set that cover `center`
        // (= satisfy q ≤ center AND q+L > center) within the cyclic
        // probe window around `center` of length up to 2*(l-1)+3.
        let r_len = 3usize;
        let probe_pad = max_lhs_len.saturating_sub(1);
        let covering_lhses = |seg_ids: &[usize], center: usize, out: &mut FxHashSet<usize>| {
            let k = seg_ids.len();
            if k < r_len {
                return;
            }
            let max_total_pad = k.saturating_sub(r_len);
            let pad_left = probe_pad.min(max_total_pad);
            let pad_right =
                probe_pad.min(max_total_pad.saturating_sub(pad_left));
            let probe_len = pad_left + r_len + pad_right;
            let probe_start = (center + k - pad_left) % k;
            let probe: Vec<usize> = (0..probe_len)
                .map(|i| seg_ids[(probe_start + i) % k])
                .collect();
            let p_in_probe = pad_left;
            for q in 0..=p_in_probe {
                let max_len = (probe_len - q).min(max_lhs_len);
                for len in 1..=max_len {
                    if q + len <= p_in_probe {
                        continue;
                    }
                    let sub: Vec<usize> = probe[q..q + len].to_vec();
                    if let Some(&id) = lhs_index.get(&sub) {
                        out.insert(id);
                    }
                }
            }
        };

        // Per-rule alive evidence: for each rule, the set of LHS-ids
        // covering its RHS on its own trial.
        let mut alive_by_rule: Vec<FxHashSet<usize>> = (0..rules.len())
            .map(|_| FxHashSet::default())
            .collect();

        // For each rule, build the trial and collect covering LHSes.
        for (rule_id, rule) in rules.iter().enumerate() {
            let src_patch = &seq_bfs.patches()[rule.witness_patch_id];
            let src_n = src_patch.boundary_len();
            let src_juncs: Vec<usize> =
                (0..src_n).filter(|&i| src_patch.is_junction(i)).collect();
            let (start_seg, _) =
                lhs_seg_range(&src_juncs, src_n, &rule.witness_pm);
            let j_outer_cw_src = src_juncs[start_seg];
            let ccw_anchor_src =
                (rule.witness_pm.start_a + rule.witness_pm.len) % src_n;
            let j_outer_cw_trial =
                (j_outer_cw_src + src_n - ccw_anchor_src) % src_n;

            let mut trial = src_patch.clone();
            if !trial.add_tile(&rule.witness_pm) {
                continue;
            }
            let trial_n = trial.boundary_len();
            let trial_juncs: Vec<usize> =
                (0..trial_n).filter(|&i| trial.is_junction(i)).collect();
            let trial_k = trial_juncs.len();
            if trial_k < 3 {
                continue;
            }
            let mut trial_sids: Vec<usize> = Vec::with_capacity(trial_k);
            for j in 0..trial_k {
                let cw_cj = trial
                    .coarse_junction_at(trial_juncs[j])
                    .expect("known junction");
                let ccw_cj = trial
                    .coarse_junction_at(trial_juncs[(j + 1) % trial_k])
                    .expect("known junction");
                let id = seq_bfs
                    .seg_bfs()
                    .lookup_seg(&cw_cj, &ccw_cj)
                    .expect("seg in DB");
                trial_sids.push(id);
            }
            let rhs_idx = match trial_juncs
                .iter()
                .position(|&p| p == j_outer_cw_trial)
            {
                Some(i) => i,
                None => continue,
            };
            debug_assert_eq!(trial_sids[rhs_idx], rule.rhs[0]);
            covering_lhses(&trial_sids, rhs_idx, &mut alive_by_rule[rule_id]);
        }

        // Step A: compute nonlocal_rules FIRST (over all rules), so
        // the fixed-point can exclude nonlocal rules from "alive
        // evidence" — a nonlocal rule is unreliable, applying it can
        // fail on the very patch we'd want it on, so we shouldn't
        // count it as a continuation guarantee.
        let known_lhses_ref: FxHashSet<&Vec<usize>> =
            table.rules_by_lhs.keys().collect();
        let mut lhs_occurrences: FxHashMap<Vec<usize>, Vec<(usize, usize)>> =
            FxHashMap::default();
        for (patch_id, patch) in seq_bfs.patches().iter().enumerate() {
            let n = patch.boundary_len();
            if n == 0 {
                continue;
            }
            let juncs: Vec<usize> =
                (0..n).filter(|&i| patch.is_junction(i)).collect();
            let k = juncs.len();
            if k < 2 {
                continue;
            }
            let mut seg_ids: Vec<usize> = Vec::with_capacity(k);
            for j in 0..k {
                let cw_cj = patch
                    .coarse_junction_at(juncs[j])
                    .expect("known junction");
                let ccw_cj = patch
                    .coarse_junction_at(juncs[(j + 1) % k])
                    .expect("known junction");
                let id = seq_bfs
                    .seg_bfs()
                    .lookup_seg(&cw_cj, &ccw_cj)
                    .expect("seg in DB");
                seg_ids.push(id);
            }
            let max_len = max_lhs_len.min(k);
            for start in 0..k {
                let mut subseq: Vec<usize> = Vec::with_capacity(max_len);
                for off in 0..max_len {
                    subseq.push(seg_ids[(start + off) % k]);
                    if known_lhses_ref.contains(&subseq) {
                        lhs_occurrences
                            .entry(subseq.clone())
                            .or_default()
                            .push((patch_id, juncs[start]));
                    }
                }
            }
        }

        let mut nonlocal_rules: FxHashSet<usize> = FxHashSet::default();
        for (rule_id, r) in rules.iter().enumerate() {
            let r_patch = &seq_bfs.patches()[r.witness_patch_id];
            let r_n = r_patch.boundary_len();
            let r_juncs: Vec<usize> =
                (0..r_n).filter(|&i| r_patch.is_junction(i)).collect();
            let (r_start_seg, _) = lhs_seg_range(&r_juncs, r_n, &r.witness_pm);
            let r_canonical = (r.witness_patch_id, r_juncs[r_start_seg]);

            let occurrences = match lhs_occurrences.get(&r.lhs) {
                Some(v) => v,
                None => continue,
            };
            for &(patch_id, j_cw) in occurrences {
                if (patch_id, j_cw) == r_canonical {
                    continue;
                }
                let w_patch = &seq_bfs.patches()[patch_id];
                let n_w = w_patch.boundary_len();
                let start_a = (j_cw + r.meta.start_edge_in_lhs) % n_w;
                let pm = PatchMatch {
                    start_a,
                    len: r.meta.match_len,
                    start_b: r.meta.tile_offset,
                    tile_id: r.meta.tile_id,
                };
                let mut trial = w_patch.clone();
                if !trial.add_tile(&pm) {
                    nonlocal_rules.insert(rule_id);
                    break;
                }
            }
        }

        // Step B: per-rule fixed-point iteration.
        //   - A rule R is "unusable" iff cursed OR nonlocal — we
        //     can't rely on it for continuation.
        //   - An LHS is "usable" iff some rule for it is non-unusable.
        //   - A rule R is cursed iff R.rhs contains a terminal OR
        //     every covering LHS in alive_by_rule[R] is non-usable.
        //   - Iterate until stable. (`nonlocal_rules` is fixed.)
        let mut cursed_rules: FxHashSet<usize> = FxHashSet::default();
        for (rule_id, r) in rules.iter().enumerate() {
            if r.rhs.iter().any(|s| terminals.contains(s)) {
                cursed_rules.insert(rule_id);
            }
        }
        let mut usable_lhs_ids: FxHashSet<usize> = FxHashSet::default();
        let mut iterations = 0usize;
        loop {
            iterations += 1;
            // Recompute usable_lhs_ids: LHS usable iff some rule for
            // it is neither cursed nor nonlocal.
            usable_lhs_ids.clear();
            for (lhs, ids) in &table.rules_by_lhs {
                let any_usable = ids.iter().any(|id| {
                    !cursed_rules.contains(id) && !nonlocal_rules.contains(id)
                });
                if any_usable {
                    if let Some(&lid) = lhs_index.get(lhs) {
                        usable_lhs_ids.insert(lid);
                    }
                }
            }
            // Recompute cursed_rules.
            let mut new_cursed_rules: FxHashSet<usize> = FxHashSet::default();
            for (rule_id, r) in rules.iter().enumerate() {
                if r.rhs.iter().any(|s| terminals.contains(s)) {
                    new_cursed_rules.insert(rule_id);
                    continue;
                }
                let any_alive = alive_by_rule[rule_id]
                    .iter()
                    .any(|id| usable_lhs_ids.contains(id));
                if !any_alive {
                    new_cursed_rules.insert(rule_id);
                }
            }
            if new_cursed_rules == cursed_rules {
                break;
            }
            cursed_rules = new_cursed_rules;
        }

        // For output `cursed_lhs`, use the stricter "all rules cursed"
        // definition (not "all rules unusable") — nonlocal is reported
        // separately.
        let mut cursed_lhs_ids: FxHashSet<usize> = FxHashSet::default();
        for (lhs, ids) in &table.rules_by_lhs {
            if ids.iter().all(|id| cursed_rules.contains(id)) {
                if let Some(&lid) = lhs_index.get(lhs) {
                    cursed_lhs_ids.insert(lid);
                }
            }
        }

        // Final dead_rhs: an RHS is dead iff every rule with that
        // RHS is cursed (= no non-cursed rule produces it).
        let mut rhs_has_noncursed_rule: FxHashMap<[usize; 3], bool> =
            FxHashMap::default();
        for (rule_id, r) in rules.iter().enumerate() {
            let entry = rhs_has_noncursed_rule.entry(r.rhs).or_insert(false);
            if !cursed_rules.contains(&rule_id) {
                *entry = true;
            }
        }
        let mut dead_rhs: FxHashSet<[usize; 3]> = FxHashSet::default();
        for (rhs, has_noncursed) in &rhs_has_noncursed_rule {
            if !*has_noncursed {
                dead_rhs.insert(*rhs);
            }
        }

        // Reify final cursed_lhs as Vec<usize> keys.
        let id_to_lhs: Vec<Vec<usize>> = {
            let mut v = vec![Vec::new(); n_lhses];
            for (lhs, &id) in &lhs_index {
                v[id] = lhs.clone();
            }
            v
        };
        let cursed_lhs: FxHashSet<Vec<usize>> = cursed_lhs_ids
            .iter()
            .map(|&id| id_to_lhs[id].clone())
            .collect();

        // context_sensitive_rhs is harder to define cleanly under
        // the fixed-point definition (since "alive" now means
        // "some non-cursed LHS covers it"); leave as an empty set
        // for now.
        let context_sensitive_rhs: FxHashSet<[usize; 3]> = FxHashSet::default();

        eprintln!("cursed-set fixed-point converged in {} iterations", iterations);

        // Sanity: by construction, every non-cursed rule's trial has
        // at least one USABLE covering LHS.
        for (rule_id, _) in rules.iter().enumerate() {
            if cursed_rules.contains(&rule_id) {
                continue;
            }
            let any_alive = alive_by_rule[rule_id]
                .iter()
                .any(|id| usable_lhs_ids.contains(id));
            assert!(
                any_alive,
                "non-cursed rule {} has no usable covering LHS on its trial",
                rule_id,
            );
        }

        RuleAnalysis {
            terminals,
            dead_rhs,
            context_sensitive_rhs,
            cursed_rules,
            cursed_lhs,
            nonlocal_rules,
        }
    }

    pub fn is_terminal(&self, seg: usize) -> bool {
        self.terminals.contains(&seg)
    }

    pub fn is_cursed_rule(&self, rule_id: usize) -> bool {
        self.cursed_rules.contains(&rule_id)
    }

    pub fn is_cursed_lhs(&self, lhs: &[usize]) -> bool {
        self.cursed_lhs.contains(lhs)
    }

    pub fn is_nonlocal_rule(&self, rule_id: usize) -> bool {
        self.nonlocal_rules.contains(&rule_id)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic::ZZ12;
    use crate::intgeom::rat::Rat;
    use crate::intgeom::tiles;
    use crate::intgeom::tileset::TileSet;
    use std::sync::{Arc, OnceLock};

    fn hex_tileset() -> Arc<TileSet<ZZ12>> {
        Arc::new(TileSet::new(vec![
            Rat::try_from(&tiles::hexagon::<ZZ12>()).unwrap(),
        ]))
    }

    fn spectre_tileset() -> Arc<TileSet<ZZ12>> {
        Arc::new(TileSet::new(vec![
            Rat::try_from(&tiles::spectre::<ZZ12>()).unwrap(),
        ]))
    }

    fn mixed_tileset() -> Arc<TileSet<ZZ12>> {
        Arc::new(TileSet::new(vec![
            Rat::try_from(&tiles::square::<ZZ12>()).unwrap(),
            Rat::try_from(&tiles::hexagon::<ZZ12>()).unwrap(),
        ]))
    }

    /// Cached SegSeqBFS for hex — running it once across all tests
    /// avoids the ~25s rebuild.
    fn hex_seq_bfs() -> &'static SegSeqBFS<ZZ12> {
        static BFS: OnceLock<SegSeqBFS<ZZ12>> = OnceLock::new();
        BFS.get_or_init(|| {
            SegSeqBFS::run(hex_tileset(), 1_000_000, 1_000_000).unwrap()
        })
    }

    #[test]
    fn hex_rewrite_rules_extracted() {
        let seq_bfs = hex_seq_bfs();
        let table = SegRewriteTable::build_with_stats(seq_bfs);
        eprintln!(
            "hex rewrites: distinct_lhs={}, total_rules={}, bfs_rules={}",
            table.num_lhs(),
            table.num_rules(),
            seq_bfs.num_rules()
        );
        assert!(table.num_rules() > 0);
        assert_eq!(
            seq_bfs.num_rules(),
            table.num_rules(),
            "table rules must equal BFS rules — one canonical witness per rule"
        );
        // Sanity: RHS is always 3.
        for metas in table.by_lhs.values() {
            for rhs in metas.values() {
                assert_eq!(rhs.len(), 3);
            }
        }
    }

    #[test]
    fn hex_save_load_roundtrip() {
        let bfs = hex_seq_bfs();
        let tmp = std::env::temp_dir().join(format!(
            "tilezz_hex_seq_bfs_{}.bin",
            std::process::id()
        ));
        bfs.save_to(&tmp).expect("save");
        let bytes = std::fs::metadata(&tmp).expect("metadata").len();
        eprintln!("hex snapshot size: {} bytes ({:.1} MB)", bytes, bytes as f64 / 1_048_576.0);
        let loaded = SegSeqBFS::<ZZ12>::load_from(hex_tileset(), &tmp).expect("load");
        let _ = std::fs::remove_file(&tmp);
        assert_eq!(loaded.num_sequences(), bfs.num_sequences());
        assert_eq!(loaded.num_rules(), bfs.num_rules());
        assert_eq!(loaded.num_witness_patches(), bfs.num_witness_patches());
        // Spot-check: same rules in same order, same canonical witness.
        for (a, b) in loaded.rules().iter().zip(bfs.rules().iter()) {
            assert_eq!(a.lhs, b.lhs);
            assert_eq!(a.meta, b.meta);
            assert_eq!(a.rhs, b.rhs);
            assert_eq!(a.witness_patch_id, b.witness_patch_id);
            assert_eq!(a.witness_pm, b.witness_pm);
        }
        // Analysis on loaded should match analysis on original.
        let t1 = SegRewriteTable::build_with_stats(bfs);
        let a1 = RuleAnalysis::build(bfs, &t1);
        let t2 = SegRewriteTable::build_with_stats(&loaded);
        let a2 = RuleAnalysis::build(&loaded, &t2);
        assert_eq!(t1.num_lhs(), t2.num_lhs());
        assert_eq!(t1.num_rules(), t2.num_rules());
        assert_eq!(a1.terminals, a2.terminals);
        assert_eq!(a1.cursed_rules, a2.cursed_rules);
        assert_eq!(a1.cursed_lhs, a2.cursed_lhs);
    }

    #[test]
    fn hex_rule_analysis() {
        let seq_bfs = hex_seq_bfs();
        let table = SegRewriteTable::build_with_stats(seq_bfs);
        let analysis = RuleAnalysis::build(seq_bfs, &table);
        eprintln!(
            "hex: terminals={}, dead_rhs={}, ctx_sens_rhs={}, cursed_rules={}/{}, cursed_lhs={}/{}, nonlocal_rules={}/{}",
            analysis.terminals.len(),
            analysis.dead_rhs.len(),
            analysis.context_sensitive_rhs.len(),
            analysis.cursed_rules.len(),
            table.num_rules(),
            analysis.cursed_lhs.len(),
            table.num_lhs(),
            analysis.nonlocal_rules.len(),
            table.num_rules(),
        );
        // Every cursed LHS has only cursed rules.
        for lhs in &analysis.cursed_lhs {
            let ids = table.rule_ids_for_lhs(lhs).expect("lhs registered");
            for id in ids {
                assert!(
                    analysis.cursed_rules.contains(id),
                    "cursed LHS has a non-cursed rule"
                );
            }
        }
        // Every cursed rule's RHS contains a terminal OR is dead.
        for &id in &analysis.cursed_rules {
            let rule = &seq_bfs.rules()[id];
            assert!(
                rule.rhs.iter().any(|s| analysis.terminals.contains(s))
                    || analysis.dead_rhs.contains(&rule.rhs),
                "cursed rule has non-terminal non-dead RHS"
            );
        }
        // No terminal appears in any LHS.
        for rule in seq_bfs.rules() {
            for s in &rule.lhs {
                assert!(
                    !analysis.terminals.contains(s),
                    "terminal seg {} appears in an LHS",
                    s
                );
            }
        }
    }

    #[test]
    #[ignore]
    fn spectre_rewrite_rules_extracted() {
        let seq_bfs = SegSeqBFS::run(spectre_tileset(), 10_000_000, 10_000_000).unwrap();
        let table = SegRewriteTable::build_with_stats(&seq_bfs);
        let analysis = RuleAnalysis::build(&seq_bfs, &table);
        eprintln!(
            "spectre: sequences(bfs)={}, distinct_lhs(rules)={}, total_rules={}, bfs_rules={}",
            seq_bfs.num_sequences(),
            table.num_lhs(),
            table.num_rules(),
            seq_bfs.num_rules()
        );
        eprintln!(
            "spectre: terminals={}, dead_rhs={}, ctx_sens_rhs={}, cursed_rules={}/{}, cursed_lhs={}/{}, nonlocal_rules={}/{}",
            analysis.terminals.len(),
            analysis.dead_rhs.len(),
            analysis.context_sensitive_rhs.len(),
            analysis.cursed_rules.len(),
            table.num_rules(),
            analysis.cursed_lhs.len(),
            table.num_lhs(),
            analysis.nonlocal_rules.len(),
            table.num_rules(),
        );
        assert_eq!(seq_bfs.num_rules(), table.num_rules(),
            "table rules must equal BFS rules");
        assert_eq!(seq_bfs.num_sequences(), table.num_lhs(),
            "BFS sequences must match table LHS count");
        for metas in table.by_lhs.values() {
            for rhs in metas.values() {
                assert_eq!(rhs.len(), 3);
            }
        }
    }

    #[test]
    #[ignore]
    fn mixed_rewrite_rules_extracted() {
        let seq_bfs = SegSeqBFS::run(mixed_tileset(), 20_000_000, 20_000_000).unwrap();
        let table = SegRewriteTable::build_with_stats(&seq_bfs);
        let analysis = RuleAnalysis::build(&seq_bfs, &table);
        eprintln!(
            "mixed: sequences(bfs)={}, distinct_lhs(rules)={}, total_rules={}, bfs_rules={}",
            seq_bfs.num_sequences(),
            table.num_lhs(),
            table.num_rules(),
            seq_bfs.num_rules()
        );
        eprintln!(
            "mixed: terminals={}, dead_rhs={}, ctx_sens_rhs={}, cursed_rules={}/{}, cursed_lhs={}/{}, nonlocal_rules={}/{}",
            analysis.terminals.len(),
            analysis.dead_rhs.len(),
            analysis.context_sensitive_rhs.len(),
            analysis.cursed_rules.len(),
            table.num_rules(),
            analysis.cursed_lhs.len(),
            table.num_lhs(),
            analysis.nonlocal_rules.len(),
            table.num_rules(),
        );
        assert_eq!(seq_bfs.num_rules(), table.num_rules(),
            "table rules must equal BFS rules");
        assert_eq!(seq_bfs.num_sequences(), table.num_lhs(),
            "BFS sequences must match table LHS count");
        for metas in table.by_lhs.values() {
            for rhs in metas.values() {
                assert_eq!(rhs.len(), 3);
            }
        }
    }
}
