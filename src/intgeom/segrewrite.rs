//! View over the rewrite rules computed inside [`SegSeqBFS`].
//!
//! The BFS computes everything (LHS, MatchMeta, RHS, witness) per
//! rule; this module just provides:
//!   - lookups by LHS (`by_lhs`, `rules_by_lhs`),
//!   - per-rule access by id (passthrough to the BFS),
//!   - segment classification: terminals and the cursed attractor.

use rustc_hash::{FxHashMap, FxHashSet};

use crate::cyclotomic::{IsComplex, IsRingOrField, Units};
use crate::intgeom::segseqbfs::{MatchMeta, RuleRecord, SegSeqBFS};

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
#[derive(Debug)]
pub struct RuleAnalysis {
    pub terminals: FxHashSet<usize>,
    pub cursed_rules: FxHashSet<usize>,
    pub cursed_lhs: FxHashSet<Vec<usize>>,
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

        RuleAnalysis { terminals, cursed_rules, cursed_lhs }
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
    fn hex_rule_analysis() {
        let seq_bfs = hex_seq_bfs();
        let table = SegRewriteTable::build_with_stats(seq_bfs);
        let analysis = RuleAnalysis::build(seq_bfs, &table);
        eprintln!(
            "hex: terminals={}, cursed_rules={}/{}, cursed_lhs={}/{}",
            analysis.terminals.len(),
            analysis.cursed_rules.len(),
            table.num_rules(),
            analysis.cursed_lhs.len(),
            table.num_lhs(),
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
            "spectre: terminals={}, cursed_rules={}/{}, cursed_lhs={}/{}",
            analysis.terminals.len(),
            analysis.cursed_rules.len(),
            table.num_rules(),
            analysis.cursed_lhs.len(),
            table.num_lhs(),
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
            "mixed: terminals={}, cursed_rules={}/{}, cursed_lhs={}/{}",
            analysis.terminals.len(),
            analysis.cursed_rules.len(),
            table.num_rules(),
            analysis.cursed_lhs.len(),
            table.num_lhs(),
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
