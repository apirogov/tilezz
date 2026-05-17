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
/// Three escalating "blocked" categories at the segment / RHS / rule
/// / LHS levels:
///
/// 1. **Terminal segment**: a seg id that appears in some rule's RHS
///    but never in any rule's LHS. Once introduced, it cannot be
///    rewritten away — there is no rule that pattern-matches it.
///
/// 2. **Dead RHS**: the 3-seg RHS, considered as a substring of the
///    boundary of any patch where it appears (= the rule's own trial
///    *and* any registered witness whose boundary happens to contain
///    it), cannot be acted on by any rule. Concretely: on every such
///    patch, no LHS substring of the boundary at some position
///    `[q, q+L)` satisfies `q ≤ p` AND `q + L > p`, where `p` is the
///    RHS's start position. (= no LHS "begins at or before the RHS
///    and rewrites at least the RHS's first segment".)
///
///    Containing a terminal seg trivially implies dead — but the
///    converse is not always true. In principle a tileset could have
///    non-terminal-containing RHSes that are still dead: every
///    individual seg is in some LHS, but no length-`L'` LHS covers
///    this specific 3-seg combination from the left. For both
///    spectre and hex, no such cases occur (dead_rhs ⊆ rhs-with-
///    terminal), but the check is kept for tilesets where it might.
///
/// 3. **Context-sensitive RHS**: alive on some patches that contain
///    it but dead on others. Tells us aliveness depends on the
///    surrounding context, not just on the RHS itself.
///
/// 4. **Cursed rule**: RHS contains a terminal OR RHS is dead. Rule
///    is a one-way trip into something we can't undo.
///
/// 5. **Cursed LHS**: every rule with this LHS is cursed. From this
///    LHS there is no escape that avoids inserting irreducibility.
///
/// Independent of cursedness:
///
/// 6. **Nonlocal rule**: rule succeeds on its canonical witness but
///    `add_tile` fails on at least one other patch that contains the
///    same LHS. The LHS alone doesn't determine the rule's
///    applicability — context beyond L matters. Computed only over
///    non-cursed rules.
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

        // Dead RHS detection (probe-based, following the user's
        // description literally).
        //
        // For each RHS r, for each witness W with r appearing in W's
        // boundary at position p:
        //   - Build a "probe" = the boundary cyclic substring of
        //     length up to (2*(l-1) + 3), centered around r with up
        //     to l-1 segs of left and right context (clipped at the
        //     point where the extension would wrap past r from the
        //     other side, on small patches).
        //   - Search the probe for an LHS substring at position
        //     [q, q+L) with q ≤ p_in_probe AND q + L > p_in_probe —
        //     i.e. the LHS starts at or before r and would rewrite
        //     at least r's first segment.
        //   - If any such LHS exists in any probe → r is alive.
        let known_lhses: FxHashSet<&Vec<usize>> =
            table.rules_by_lhs.keys().collect();
        let max_lhs_len = table
            .rules_by_lhs
            .keys()
            .map(|l| l.len())
            .max()
            .unwrap_or(0);

        // For each RHS, track two flags: has_alive (seen on some
        // patch at a left-coverable position) and has_dead (seen on
        // some patch at a position not left-coverable).
        let r_len = 3usize;
        let probe_pad = max_lhs_len.saturating_sub(1);
        let known_rhses: FxHashSet<[usize; 3]> =
            rules.iter().map(|r| r.rhs).collect();

        // Closure: given a patch's seg_ids and a probe center index,
        // build the probe and check whether the center is alive.
        let probe_is_alive = |seg_ids: &[usize], center: usize| -> bool {
            let k = seg_ids.len();
            if k < r_len {
                return false;
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
                    if known_lhses.contains(&sub) {
                        return true;
                    }
                }
            }
            false
        };

        let mut rhs_flags: FxHashMap<[usize; 3], (bool, bool)> =
            FxHashMap::default(); // (has_alive, has_dead)

        // Pass 1: each rule's trial — RHS at known position on the
        // un-normalized trial.
        for rule in rules {
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

            let alive = probe_is_alive(&trial_sids, rhs_idx);
            let entry = rhs_flags.entry(rule.rhs).or_insert((false, false));
            if alive {
                entry.0 = true;
            } else {
                entry.1 = true;
            }
        }

        // Pass 2: every registered witness — scan its boundary for
        // length-3 substrings that match a known RHS, check each.
        for patch in seq_bfs.patches() {
            let n = patch.boundary_len();
            if n == 0 {
                continue;
            }
            let juncs: Vec<usize> =
                (0..n).filter(|&i| patch.is_junction(i)).collect();
            let k = juncs.len();
            if k < 3 {
                continue;
            }
            let mut sids: Vec<usize> = Vec::with_capacity(k);
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
                sids.push(id);
            }
            for s in 0..k {
                let triple = [sids[s], sids[(s + 1) % k], sids[(s + 2) % k]];
                if !known_rhses.contains(&triple) {
                    continue;
                }
                let alive = probe_is_alive(&sids, s);
                let entry = rhs_flags.entry(triple).or_insert((false, false));
                if alive {
                    entry.0 = true;
                } else {
                    entry.1 = true;
                }
            }
        }

        // Classify each known RHS.
        let mut dead_rhs: FxHashSet<[usize; 3]> = FxHashSet::default();
        let mut context_sensitive_rhs: FxHashSet<[usize; 3]> = FxHashSet::default();
        for rhs in &known_rhses {
            match rhs_flags.get(rhs).copied().unwrap_or((false, false)) {
                (false, _) => {
                    dead_rhs.insert(*rhs);
                }
                (true, true) => {
                    context_sensitive_rhs.insert(*rhs);
                }
                (true, false) => {}
            }
        }

        // Cursed rule = RHS contains terminal OR RHS is dead.
        let mut cursed_rules: FxHashSet<usize> = FxHashSet::default();
        for (rule_id, r) in rules.iter().enumerate() {
            if r.rhs.iter().any(|s| terminals.contains(s))
                || dead_rhs.contains(&r.rhs)
            {
                cursed_rules.insert(rule_id);
            }
        }

        // Cursed LHS: every rule with that LHS is cursed.
        let mut cursed_lhs: FxHashSet<Vec<usize>> = FxHashSet::default();
        for (lhs, ids) in &table.rules_by_lhs {
            if ids.iter().all(|id| cursed_rules.contains(id)) {
                cursed_lhs.insert(lhs.clone());
            }
        }

        // Nonlocal rules: a rule R with LHS L is nonlocal iff there
        // exists ANY patch in the BFS whose boundary contains L at
        // some position where R's match fails `add_tile` — not just
        // canonical witnesses for rules with L. We build an inverted
        // index `lhs_occurrences: L -> [(patch_id, j_outer_cw_pos)]`
        // by scanning every witness patch's cyclic seg sequence for
        // every contiguous subsequence that matches a known LHS.
        let known_lhses: FxHashSet<&Vec<usize>> = table.rules_by_lhs.keys().collect();
        let max_lhs_len = table.rules_by_lhs.keys().map(|l| l.len()).max().unwrap_or(0);
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
                let cw_cj = patch.coarse_junction_at(juncs[j]).expect("known junction");
                let ccw_cj = patch.coarse_junction_at(juncs[(j + 1) % k]).expect("known junction");
                let id = seq_bfs
                    .seg_bfs()
                    .lookup_seg(&cw_cj, &ccw_cj)
                    .expect("seg in DB");
                seg_ids.push(id);
            }
            // Enumerate cyclic subsequences of length 1..=min(max_lhs_len, k)
            // starting at each segment index.
            let max_len = max_lhs_len.min(k);
            for start in 0..k {
                let mut subseq: Vec<usize> = Vec::with_capacity(max_len);
                for off in 0..max_len {
                    subseq.push(seg_ids[(start + off) % k]);
                    if known_lhses.contains(&subseq) {
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
            // Cursed rules are already known unusable (= they inject a
            // terminal). Don't bother testing locality for them.
            if cursed_rules.contains(&rule_id) {
                continue;
            }
            // Compute R's own canonical L-position to skip it.
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
