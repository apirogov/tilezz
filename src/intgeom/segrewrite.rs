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
/// - **Terminals** = segment ids that appear in some rule's RHS but
///   never in any rule's LHS. By definition there's no rewrite that
///   replaces them; once present they stay.
/// - **Cursed rule** = a rule whose RHS contains a terminal segment.
///   Applying it injects a terminal that cannot be eliminated later.
/// - **Cursed LHS** = an LHS for which *every* applicable rule is
///   cursed. From such an LHS there's no escape that avoids inserting
///   a terminal.
/// - **Nonlocal rule** = a rule that succeeds on its own canonical
///   witness but fails (`add_tile` returns false) on at least one
///   *other* witness with the same LHS. The LHS alone is not
///   sufficient to know whether the rule is applicable — context
///   beyond L matters. Strictly weaker than "cursed" (= just
///   unreliable).
#[derive(Debug)]
pub struct RuleAnalysis {
    pub terminals: FxHashSet<usize>,
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

        // Cursed rules: any RHS seg is a terminal.
        let mut cursed_rules: FxHashSet<usize> = FxHashSet::default();
        for (rule_id, r) in rules.iter().enumerate() {
            if r.rhs.iter().any(|s| terminals.contains(s)) {
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
            "hex: terminals={}, cursed_rules={}/{}, cursed_lhs={}/{}, nonlocal_rules={}/{}",
            analysis.terminals.len(),
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
        // Every cursed rule's RHS contains some terminal.
        for &id in &analysis.cursed_rules {
            let rule = &seq_bfs.rules()[id];
            assert!(
                rule.rhs.iter().any(|s| analysis.terminals.contains(s)),
                "cursed rule has no terminal in RHS"
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
            "spectre: terminals={}, cursed_rules={}/{}, cursed_lhs={}/{}, nonlocal_rules={}/{}",
            analysis.terminals.len(),
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
            "mixed: terminals={}, cursed_rules={}/{}, cursed_lhs={}/{}, nonlocal_rules={}/{}",
            analysis.terminals.len(),
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
