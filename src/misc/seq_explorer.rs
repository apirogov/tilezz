use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::sync::Arc;
use std::time::Instant;

use rustc_hash::{FxHashMap, FxHashSet};

use crate::cyclotomic::{IsComplex, IsRingOrField, Units};
use crate::intgeom::matchtypes::{BpSeed, MatchFinder};
use crate::intgeom::rat::Rat;
use crate::intgeom::snake::Snake;
use crate::intgeom::tileset::TileSet;

pub fn replay_glue<T: IsComplex + IsRingOrField + Units>(
    source: &Rat<T>,
    tile: &Rat<T>,
    start_a: usize,
    start_b: usize,
    len: usize,
) -> Result<Rat<T>, String> {
    source
        .try_glue_precomputed((start_a as i64, len, start_b as i64), tile, true)
        .map_err(|e| e.to_string())
}

struct TrieNode {
    children: HashMap<i8, usize>,
    witness: Option<usize>,
}

impl TrieNode {
    fn new() -> Self {
        TrieNode {
            children: HashMap::new(),
            witness: None,
        }
    }
}

struct SeqTrie {
    nodes: Vec<TrieNode>,
}

impl SeqTrie {
    fn new() -> Self {
        SeqTrie {
            nodes: vec![TrieNode::new()],
        }
    }

    fn insert_cyclic_subseqs(
        &mut self,
        seq: &[i8],
        k: usize,
        witness_idx: usize,
    ) -> (usize, Vec<bool>) {
        let n = seq.len();
        if n == 0 {
            return (0, Vec::new());
        }
        let max_len = k.min(n);
        let mut new_count = 0;
        let mut active = vec![false; n];

        for start in 0..n {
            let mut node = 0;
            let mut any_new = false;
            for l in 0..max_len {
                let angle = seq[(start + l) % n];
                let next = self.nodes[node].children.get(&angle).copied();
                node = match next {
                    Some(child) => child,
                    None => {
                        let new_node = self.nodes.len();
                        self.nodes.push(TrieNode::new());
                        self.nodes[node].children.insert(angle, new_node);
                        self.nodes[new_node].witness = Some(witness_idx);
                        new_count += 1;
                        any_new = true;
                        new_node
                    }
                };
            }
            if any_new {
                active[start] = true;
            }
        }

        (new_count, active)
    }

    fn len(&self) -> usize {
        self.nodes.len() - 1
    }

    fn sequences_by_witness(&self) -> BTreeMap<usize, Vec<Vec<i8>>> {
        let mut result: BTreeMap<usize, Vec<Vec<i8>>> = BTreeMap::new();
        self.dfs_collect(0, &mut Vec::new(), &mut result);
        result
    }

    fn dfs_collect(
        &self,
        node: usize,
        path: &mut Vec<i8>,
        result: &mut BTreeMap<usize, Vec<Vec<i8>>>,
    ) {
        if let Some(w) = self.nodes[node].witness {
            result.entry(w).or_default().push(path.clone());
        }
        let mut sorted_children: Vec<(i8, usize)> = self.nodes[node]
            .children
            .iter()
            .map(|(&a, &c)| (a, c))
            .collect();
        sorted_children.sort_by_key(|(a, _)| *a);
        for (angle, child) in sorted_children {
            path.push(angle);
            self.dfs_collect(child, path, result);
            path.pop();
        }
    }
}

fn junction_active_mask(
    raw_active: &[bool],
    k: usize,
    provenance: &Provenance,
    source_len: usize,
) -> Vec<bool> {
    let n = raw_active.len();
    match provenance {
        Provenance::Seed { .. } => vec![true; n],
        Provenance::Glue { len, .. } => {
            let cw_j = source_len - len;
            let ccw_j = 0usize;
            let mut active = vec![false; n];
            let check_radius = k.min(n);
            let dilate = k.saturating_sub(1).min(n);

            for &j in &[ccw_j, cw_j] {
                let mut induced = false;
                for d in 0..check_radius {
                    if raw_active[(j + n - d) % n] {
                        induced = true;
                        break;
                    }
                }
                if induced {
                    for d in 0..=dilate {
                        active[(j + n - d) % n] = true;
                        active[(j + d) % n] = true;
                    }
                }
            }
            active
        }
    }
}

#[derive(Debug, Clone)]
pub enum Provenance {
    Seed {
        tile_idx: usize,
    },
    Glue {
        source_rat_id: usize,
        tile_idx: usize,
        start_a: usize,
        start_b: usize,
        len: usize,
    },
}

pub struct SeqExplorer<T: IsComplex> {
    tileset: Arc<TileSet<T>>,
    trie: SeqTrie,
    rats: Vec<Rat<T>>,
    provenances: Vec<Provenance>,
    rat_to_id: FxHashMap<Rat<T>, usize>,
    witnesses: Vec<usize>,
    max_subseq_len: usize,
    total_matches_skipped: usize,
    total_matches_considered: usize,
}

impl<T: IsComplex + IsRingOrField + Units> SeqExplorer<T> {
    pub fn new(tileset: Arc<TileSet<T>>) -> Self {
        let t0 = Instant::now();
        let k = tileset.rats().iter().map(|r| r.len()).max().unwrap_or(0);

        let mut trie = SeqTrie::new();
        let mut rats: Vec<Rat<T>> = Vec::new();
        let mut provenances: Vec<Provenance> = Vec::new();
        let mut rat_to_id: FxHashMap<Rat<T>, usize> = FxHashMap::default();
        let mut witnesses: Vec<usize> = Vec::new();
        let mut pending: BTreeSet<usize> = BTreeSet::new();
        let mut active_map: HashMap<usize, Vec<bool>> = HashMap::new();

        for tile_idx in 0..tileset.num_tiles() {
            let rat = tileset.rat(tile_idx).clone();
            if let std::collections::hash_map::Entry::Vacant(e) = rat_to_id.entry(rat.clone()) {
                let rat_id = rats.len();
                e.insert(rat_id);
                rats.push(rat.clone());
                provenances.push(Provenance::Seed { tile_idx });

                let witness_idx = witnesses.len();
                let (new_count, raw_active) = trie.insert_cyclic_subseqs(rat.seq(), k, witness_idx);
                if new_count > 0 {
                    let n = rat.seq().len();
                    let mask =
                        junction_active_mask(&raw_active, k, &Provenance::Seed { tile_idx }, n);
                    active_map.insert(rat_id, mask);
                    witnesses.push(rat_id);
                    pending.insert(rat_id);
                }
            }
        }

        eprintln!(
            "  Seeds: {} rats, {} subseqs, {} pending",
            rats.len(),
            trie.len(),
            pending.len(),
        );

        let mut layer = 0;
        let mut total_skipped = 0usize;
        let mut total_considered = 0usize;

        // The B-side (`tileset`) is fixed across all layers, so its
        // bit-parallel masks can be precomputed once and reused for
        // every batch. The A-side (`batch_ts`) changes per layer.
        let bp_seed = BpSeed::new(Arc::clone(&tileset));

        while !pending.is_empty() {
            layer += 1;
            let batch_ids: Vec<usize> = pending.into_iter().collect();
            pending = BTreeSet::new();

            let batch_rats: Vec<Rat<T>> = batch_ids.iter().map(|&id| rats[id].clone()).collect();
            let batch_ts = Arc::new(TileSet::new(batch_rats));
            assert_eq!(batch_ts.num_tiles(), batch_ids.len());

            let global_from_batch: Vec<usize> = {
                let ts = &batch_ts;
                (0..ts.num_tiles())
                    .map(|batch_idx| {
                        let batch_rat = ts.rat(batch_idx);
                        *rat_to_id.get(batch_rat).expect("batch rat in registry")
                    })
                    .collect()
            };

            let batch_active: Vec<Vec<bool>> = (0..batch_ts.num_tiles())
                .map(|batch_idx| {
                    let global_id = global_from_batch[batch_idx];
                    active_map.get(&global_id).cloned().unwrap_or_default()
                })
                .collect();

            let mf = MatchFinder::crossing_with_seed(Arc::clone(&batch_ts), bp_seed.clone());

            let mut new_rat_entries: Vec<(Rat<T>, Provenance)> = Vec::new();
            let mut seen_new: FxHashSet<Rat<T>> = FxHashSet::default();
            let mut seen_invalid: FxHashSet<Rat<T>> = FxHashSet::default();
            let layer_skipped = 0usize;
            let mut layer_considered = 0usize;

            for batch_idx in 0..batch_ts.num_tiles() {
                let _n = batch_ts.rat(batch_idx).len();
                let active = &batch_active[batch_idx];
                let global_source = global_from_batch[batch_idx];

                for tile_idx in 0..tileset.num_tiles() {
                    for (glued, mt) in mf.valid_matches_filtered(batch_idx, tile_idx, active) {
                        layer_considered += 1;

                        if rat_to_id.contains_key(&glued)
                            || seen_new.contains(&glued)
                            || seen_invalid.contains(&glued)
                        {
                            continue;
                        }

                        if Snake::<T>::try_from(glued.seq()).is_err() {
                            seen_invalid.insert(glued);
                            continue;
                        }

                        seen_new.insert(glued.clone());
                        new_rat_entries.push((
                            glued,
                            Provenance::Glue {
                                source_rat_id: global_source,
                                tile_idx,
                                start_a: mt.start_a,
                                start_b: mt.start_b,
                                len: mt.len,
                            },
                        ));
                    }
                }
            }

            total_skipped += layer_skipped;
            total_considered += layer_considered;

            let new_count_total = new_rat_entries.len();
            let mut len_min = usize::MAX;
            let mut len_max = 0usize;

            for (rat, _) in &new_rat_entries {
                let len = rat.len();
                len_min = len_min.min(len);
                len_max = len_max.max(len);
            }

            let mut new_from_layer = 0usize;
            for (rat, prov) in new_rat_entries {
                let rat_id = rats.len();
                rats.push(rat.clone());
                provenances.push(prov);
                rat_to_id.insert(rat.clone(), rat_id);

                let witness_idx = witnesses.len();
                let (new_count, raw_active) = trie.insert_cyclic_subseqs(rat.seq(), k, witness_idx);
                if new_count > 0 {
                    let source_len = match &provenances[rat_id] {
                        Provenance::Seed { .. } => rat.len(),
                        Provenance::Glue { source_rat_id, .. } => rats[*source_rat_id].len(),
                    };
                    let mask =
                        junction_active_mask(&raw_active, k, &provenances[rat_id], source_len);
                    active_map.insert(rat_id, mask);
                    witnesses.push(rat_id);
                    pending.insert(rat_id);
                    new_from_layer += 1;
                }
            }

            eprintln!(
                "  Layer {}: batch={} new_rats={} new_wit={} trie={} pending={} lens=[{},{}] skip={}/{}",
                layer,
                batch_ts.num_tiles(),
                new_count_total,
                new_from_layer,
                trie.len(),
                pending.len(),
                len_min,
                len_max,
                layer_skipped,
                layer_skipped + layer_considered,
            );
        }

        eprintln!(
            "  Done: {} layers, {} rats, {} subseqs, {} witnesses, skip={}/{} time={:.2?}",
            layer,
            rats.len(),
            trie.len(),
            witnesses.len(),
            total_skipped,
            total_skipped + total_considered,
            t0.elapsed(),
        );

        SeqExplorer {
            tileset,
            trie,
            rats,
            provenances,
            rat_to_id,
            witnesses,
            max_subseq_len: k,
            total_matches_skipped: total_skipped,
            total_matches_considered: total_considered,
        }
    }

    pub fn num_subseqs(&self) -> usize {
        self.trie.len()
    }

    pub fn num_witnesses(&self) -> usize {
        self.witnesses.len()
    }

    pub fn num_known_rats(&self) -> usize {
        self.rats.len()
    }

    pub fn max_subseq_len(&self) -> usize {
        self.max_subseq_len
    }

    pub fn rats(&self) -> &[Rat<T>] {
        &self.rats
    }

    pub fn provenances(&self) -> &[Provenance] {
        &self.provenances
    }

    pub fn witnesses(&self) -> &[usize] {
        &self.witnesses
    }

    pub fn witness_rat_ids(&self) -> &[usize] {
        &self.witnesses
    }

    pub fn provenance(&self, rat_id: usize) -> &Provenance {
        &self.provenances[rat_id]
    }

    pub fn rat(&self, id: usize) -> &Rat<T> {
        &self.rats[id]
    }

    pub fn rat_id_of(&self, rat: &Rat<T>) -> Option<usize> {
        self.rat_to_id.get(rat).copied()
    }

    pub fn sequences_by_witness(&self) -> BTreeMap<usize, Vec<Vec<i8>>> {
        self.trie.sequences_by_witness()
    }

    pub fn tileset(&self) -> &Arc<TileSet<T>> {
        &self.tileset
    }

    pub fn stats(&self) -> (usize, usize) {
        (self.total_matches_skipped, self.total_matches_considered)
    }

    pub fn export_collected(&self) -> BTreeMap<usize, Vec<Vec<i8>>> {
        let by_idx = self.trie.sequences_by_witness();
        let mut result = BTreeMap::new();
        for (witness_local_idx, seqs) in by_idx {
            let rat_id = self.witnesses[witness_local_idx];
            result.insert(rat_id, seqs);
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic::ZZ12;
    use crate::intgeom::snake::Snake;
    use crate::intgeom::tiles;

    fn hex_ts() -> Arc<TileSet<ZZ12>> {
        let hex: Snake<ZZ12> = tiles::hexagon();
        let rat = Rat::try_from(&hex).unwrap();
        Arc::new(TileSet::new(vec![rat]))
    }

    fn sq_ts() -> Arc<TileSet<ZZ12>> {
        let sq: Snake<ZZ12> = tiles::square();
        let rat = Rat::try_from(&sq).unwrap();
        Arc::new(TileSet::new(vec![rat]))
    }

    #[test]
    fn hex_seq_explorer() {
        let explorer = SeqExplorer::new(hex_ts());
        eprintln!(
            "[hex] subseqs={} witnesses={} rats={} k={}",
            explorer.num_subseqs(),
            explorer.num_witnesses(),
            explorer.num_known_rats(),
            explorer.max_subseq_len(),
        );
        assert!(explorer.num_subseqs() > 0);
        assert!(explorer.num_known_rats() > 0);
        assert!(explorer.num_witnesses() <= explorer.num_subseqs());
    }

    #[test]
    fn square_seq_explorer() {
        let explorer = SeqExplorer::new(sq_ts());
        eprintln!(
            "[square] subseqs={} witnesses={} rats={} k={}",
            explorer.num_subseqs(),
            explorer.num_witnesses(),
            explorer.num_known_rats(),
            explorer.max_subseq_len(),
        );
        assert!(explorer.num_subseqs() > 0);
        assert!(explorer.num_known_rats() > 0);
        assert!(explorer.num_witnesses() <= explorer.num_subseqs());
    }

    fn mixed_ts() -> Arc<TileSet<ZZ12>> {
        let sq: Snake<ZZ12> = tiles::square();
        let hex: Snake<ZZ12> = tiles::hexagon();
        let sq_rat = Rat::try_from(&sq).unwrap();
        let hex_rat = Rat::try_from(&hex).unwrap();
        Arc::new(TileSet::new(vec![sq_rat, hex_rat]))
    }

    #[test]
    #[ignore]
    fn mixed_seq_explorer() {
        let explorer = SeqExplorer::new(mixed_ts());
        eprintln!(
            "[mixed] subseqs={} witnesses={} rats={} k={}",
            explorer.num_subseqs(),
            explorer.num_witnesses(),
            explorer.num_known_rats(),
            explorer.max_subseq_len(),
        );
    }
}
