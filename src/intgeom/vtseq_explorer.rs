use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::sync::Arc;

use crate::cyclotomic::{IsComplex, IsRingOrField, Units};
use crate::intgeom::patch::{GrowingPatch, PatchMatch};
use crate::intgeom::rat::Rat;
use crate::intgeom::tileset::TileSet;

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
        Provenance::Extend { patch_match, .. } => {
            let len = patch_match.len;
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

fn match_touches_active(start_a: usize, mlen: usize, n: usize, active: &[bool]) -> bool {
    if mlen == 0 || n == 0 {
        return false;
    }
    let before = (start_a + n - 1) % n;
    if active[before] {
        return true;
    }
    let after = (start_a + mlen) % n;
    if active[after] {
        return true;
    }
    for i in 0..mlen {
        if active[(start_a + i) % n] {
            return true;
        }
    }
    false
}

#[derive(Debug, Clone)]
pub enum Provenance {
    Seed {
        tile_idx: usize,
    },
    Extend {
        source_witness_id: usize,
        patch_match: PatchMatch,
    },
}

pub struct PatchSeqExplorer<T: IsComplex> {
    tileset: Arc<TileSet<T>>,
    trie: SeqTrie,
    witnesses: Vec<GrowingPatch<T>>,
    provenances: Vec<Provenance>,
    rat_to_id: BTreeMap<Rat<T>, usize>,
    witness_indices: Vec<usize>,
    max_subseq_len: usize,
    total_matches_skipped: usize,
    total_matches_considered: usize,
}

impl<T: IsComplex + IsRingOrField + Units> PatchSeqExplorer<T> {
    pub fn new(tileset: Arc<TileSet<T>>) -> Self {
        let k = tileset.rats().iter().map(|r| r.len()).max().unwrap_or(0);

        let mut trie = SeqTrie::new();
        let mut witnesses: Vec<GrowingPatch<T>> = Vec::new();
        let mut provenances: Vec<Provenance> = Vec::new();
        let mut rat_to_id: BTreeMap<Rat<T>, usize> = BTreeMap::new();
        let mut witness_indices: Vec<usize> = Vec::new();
        let mut pending: BTreeSet<usize> = BTreeSet::new();
        let mut active_map: HashMap<usize, Vec<bool>> = HashMap::new();

        for tile_idx in 0..tileset.num_tiles() {
            let tile_rat = tileset.rat(tile_idx).clone();
            if rat_to_id.contains_key(&tile_rat) {
                continue;
            }

            let seed = GrowingPatch::new(Arc::clone(&tileset), tile_idx);
            let rat_id = witnesses.len();
            rat_to_id.insert(tile_rat.clone(), rat_id);
            provenances.push(Provenance::Seed { tile_idx });
            witnesses.push(seed);

            let witness_idx = witness_indices.len();
            let seq = tile_rat.seq();
            let (new_count, raw_active) = trie.insert_cyclic_subseqs(seq, k, witness_idx);
            if new_count > 0 {
                let n = seq.len();
                let mask = junction_active_mask(&raw_active, k, &Provenance::Seed { tile_idx }, n);
                active_map.insert(rat_id, mask);
                witness_indices.push(rat_id);
                pending.insert(rat_id);
            }
        }

        for tile_idx in 0..tileset.num_tiles() {
            let seed = GrowingPatch::new(Arc::clone(&tileset), tile_idx);
            for pm in seed.get_all_matches() {
                let mut gp = seed.clone_for_mutation();
                if gp.add_tile(&pm).is_none() || !gp.is_growing() {
                    continue;
                }
                let rat = gp.to_rat();
                if rat_to_id.contains_key(&rat) {
                    continue;
                }

                let rat_id = witnesses.len();
                rat_to_id.insert(rat.clone(), rat_id);
                provenances.push(Provenance::Seed { tile_idx });

                let seq = rat.seq();
                let witness_idx = witness_indices.len();
                let (new_count, raw_active) = trie.insert_cyclic_subseqs(seq, k, witness_idx);
                if new_count > 0 {
                    let n = seq.len();
                    let mask =
                        junction_active_mask(&raw_active, k, &Provenance::Seed { tile_idx }, n);
                    active_map.insert(rat_id, mask);
                    witness_indices.push(rat_id);
                    witnesses.push(gp);
                    pending.insert(rat_id);
                }
            }
        }

        eprintln!(
            "  Seeds: {} patches, {} subseqs, {} pending",
            witnesses.len(),
            trie.len(),
            pending.len(),
        );

        let mut layer = 0u32;
        let mut total_skipped = 0usize;
        let mut total_considered = 0usize;

        while !pending.is_empty() {
            layer += 1;
            let batch_ids: Vec<usize> = pending.into_iter().collect();
            pending = BTreeSet::new();

            let batch_active: Vec<Vec<bool>> = batch_ids
                .iter()
                .map(|&id| active_map.get(&id).cloned().unwrap_or_default())
                .collect();

            let mut new_entries: Vec<(GrowingPatch<T>, Provenance, usize)> = Vec::new();
            let mut seen_new: BTreeSet<Rat<T>> = BTreeSet::new();
            let mut layer_skipped = 0usize;
            let mut layer_considered = 0usize;

            for (batch_pos, &source_id) in batch_ids.iter().enumerate() {
                let patch = &witnesses[source_id];
                let n = patch.boundary_len();

                let matches = patch.get_all_matches();
                for pm in &matches {
                    if n > 0 {
                        let active = &batch_active[batch_pos];
                        if !match_touches_active(pm.start_a, pm.len, n, active) {
                            layer_skipped += 1;
                            continue;
                        }
                    }
                    layer_considered += 1;

                    let mut gp = patch.clone_for_mutation();
                    if gp.add_tile(pm).is_none() || !gp.is_growing() {
                        continue;
                    }

                    let rat = gp.to_rat();
                    if !rat_to_id.contains_key(&rat) && seen_new.insert(rat.clone()) {
                        let source_len = if n > 0 {
                            n
                        } else {
                            match &provenances[source_id] {
                                Provenance::Seed { tile_idx } => tileset.rat(*tile_idx).seq().len(),
                                Provenance::Extend { .. } => n,
                            }
                        };
                        new_entries.push((
                            gp,
                            Provenance::Extend {
                                source_witness_id: source_id,
                                patch_match: pm.clone(),
                            },
                            source_len,
                        ));
                    }
                }
            }

            total_skipped += layer_skipped;
            total_considered += layer_considered;

            let new_count_total = new_entries.len();
            let mut new_from_layer = 0usize;

            for (gp, prov, source_len) in new_entries {
                let rat = gp.to_rat();
                let rat_id = witnesses.len();
                rat_to_id.insert(rat.clone(), rat_id);
                provenances.push(prov);
                witnesses.push(gp);

                let witness_idx = witness_indices.len();
                let seq = rat.seq();
                let (new_count, raw_active) = trie.insert_cyclic_subseqs(seq, k, witness_idx);
                if new_count > 0 {
                    let mask =
                        junction_active_mask(&raw_active, k, &provenances[rat_id], source_len);
                    active_map.insert(rat_id, mask);
                    witness_indices.push(rat_id);
                    pending.insert(rat_id);
                    new_from_layer += 1;
                }
            }

            eprintln!(
                "  Layer {}: batch={} new_patches={} new_wit={} trie={} pending={} skip={}/{}",
                layer,
                batch_ids.len(),
                new_count_total,
                new_from_layer,
                trie.len(),
                pending.len(),
                layer_skipped,
                layer_skipped + layer_considered,
            );
        }

        eprintln!(
            "  Done: {} layers, {} patches, {} subseqs, {} witnesses, skip={}/{}",
            layer,
            witnesses.len(),
            trie.len(),
            witness_indices.len(),
            total_skipped,
            total_skipped + total_considered,
        );

        PatchSeqExplorer {
            tileset,
            trie,
            witnesses,
            provenances,
            rat_to_id,
            witness_indices,
            max_subseq_len: k,
            total_matches_skipped: total_skipped,
            total_matches_considered: total_considered,
        }
    }

    pub fn num_subseqs(&self) -> usize {
        self.trie.len()
    }

    pub fn num_witnesses(&self) -> usize {
        self.witness_indices.len()
    }

    pub fn num_patches(&self) -> usize {
        self.witnesses.len()
    }

    pub fn max_subseq_len(&self) -> usize {
        self.max_subseq_len
    }

    pub fn witnesses(&self) -> &[GrowingPatch<T>] {
        &self.witnesses
    }

    pub fn provenances(&self) -> &[Provenance] {
        &self.provenances
    }

    pub fn witness_indices(&self) -> &[usize] {
        &self.witness_indices
    }

    pub fn sequences_by_witness(&self) -> BTreeMap<usize, Vec<Vec<i8>>> {
        self.trie.sequences_by_witness()
    }

    pub fn stats(&self) -> (usize, usize) {
        (self.total_matches_skipped, self.total_matches_considered)
    }

    pub fn tileset(&self) -> &Arc<TileSet<T>> {
        &self.tileset
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic::ZZ12;
    use crate::intgeom::rat::Rat;
    use crate::intgeom::tiles;
    use std::sync::Arc;

    #[test]
    fn hex_patch_seq_explorer() {
        let hex: crate::intgeom::snake::Snake<ZZ12> = tiles::hexagon();
        let rat = Rat::try_from(&hex).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let explorer = PatchSeqExplorer::new(ts);
        eprintln!(
            "[hex] subseqs={} witnesses={} patches={} k={}",
            explorer.num_subseqs(),
            explorer.num_witnesses(),
            explorer.num_patches(),
            explorer.max_subseq_len(),
        );
        assert!(explorer.num_subseqs() > 0);
        assert!(explorer.num_patches() > 0);
        assert!(explorer.num_witnesses() <= explorer.num_subseqs());
    }

    #[test]
    fn square_patch_seq_explorer() {
        let sq: crate::intgeom::snake::Snake<ZZ12> = tiles::square();
        let rat = Rat::try_from(&sq).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let explorer = PatchSeqExplorer::new(ts);
        eprintln!(
            "[square] subseqs={} witnesses={} patches={} k={}",
            explorer.num_subseqs(),
            explorer.num_witnesses(),
            explorer.num_patches(),
            explorer.max_subseq_len(),
        );
        assert!(explorer.num_subseqs() > 0);
        assert!(explorer.num_patches() > 0);
        assert!(explorer.num_witnesses() <= explorer.num_subseqs());
    }
}
