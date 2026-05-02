use std::hash::{Hash, Hasher};
use std::sync::Arc;

use rustc_hash::{FxHashMap, FxHasher};

use crate::cyclotomic::{IsComplex, IsRingOrField, Units};
use crate::intgeom::matchtypes::MatchTypeIndex;
use crate::intgeom::patch::{compute_glue_angles, EdgeInfo, GrowingPatch};
use crate::intgeom::snake::Snake;
use crate::intgeom::tileset::TileSet;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OpenNeighborhoodType {
    pub id: usize,
    pub central_tile_id: usize,
    pub gap_start: u8,
    pub gap_len: u8,
    pub angles: Vec<i8>,
    pub edges: Vec<EdgeInfo>,
    pub inner_chains: Vec<Vec<EdgeInfo>>,
}

impl OpenNeighborhoodType {
    pub fn gap_len(&self) -> usize {
        self.gap_len as usize
    }
}

fn dedup_hash(
    central_tile_id: usize,
    angles: &[i8],
    edges: &[EdgeInfo],
    gap_start: u8,
    gap_len: u8,
) -> u64 {
    let mut h = FxHasher::default();
    central_tile_id.hash(&mut h);
    angles.hash(&mut h);
    edges.hash(&mut h);
    gap_start.hash(&mut h);
    gap_len.hash(&mut h);
    h.finish()
}

struct SeenIndex {
    by_hash: FxHashMap<u64, Vec<usize>>,
}

impl SeenIndex {
    fn new() -> Self {
        SeenIndex {
            by_hash: FxHashMap::default(),
        }
    }

    fn insert(
        &mut self,
        all_types: &[OpenNeighborhoodType],
        central_tile_id: usize,
        ca: &[i8],
        ce: &[EdgeInfo],
        gap_start: u8,
        gap_len: u8,
    ) -> bool {
        let hash = dedup_hash(central_tile_id, ca, ce, gap_start, gap_len);
        let indices = self.by_hash.entry(hash).or_default();
        for &idx in indices.iter() {
            let t = &all_types[idx];
            if t.central_tile_id == central_tile_id
                && t.gap_start == gap_start
                && t.gap_len == gap_len
                && t.angles == ca
                && t.edges == ce
            {
                return false;
            }
        }
        indices.push(all_types.len());
        true
    }
}

fn find_best_rotation(angles: &[i8], edges: &[EdgeInfo]) -> usize {
    let n = angles.len();
    if n == 0 {
        return 0;
    }
    let mut best = 0;
    for r in 1..n {
        for i in 0..n {
            let idx_r = (r + i) % n;
            let idx_b = (best + i) % n;
            match angles[idx_r].cmp(&angles[idx_b]) {
                std::cmp::Ordering::Less => {
                    best = r;
                    break;
                }
                std::cmp::Ordering::Greater => break,
                std::cmp::Ordering::Equal => match edges[idx_r].cmp(&edges[idx_b]) {
                    std::cmp::Ordering::Less => {
                        best = r;
                        break;
                    }
                    std::cmp::Ordering::Greater => break,
                    std::cmp::Ordering::Equal => continue,
                },
            }
        }
    }
    best
}

fn canonicalize(angles: &[i8], edges: &[EdgeInfo]) -> (Vec<i8>, Vec<EdgeInfo>, usize) {
    let rot = find_best_rotation(angles, edges);
    let mut a = angles.to_vec();
    let mut e = edges.to_vec();
    a.rotate_left(rot);
    e.rotate_left(rot);
    (a, e, rot)
}

fn in_consumed_range(pos: usize, start_a: usize, match_len: usize, n: usize) -> bool {
    if match_len == 0 || n == 0 {
        return false;
    }
    if match_len >= n {
        return true;
    }
    let end = (start_a + match_len) % n;
    if end == 0 {
        pos >= start_a
    } else if start_a < end {
        pos >= start_a && pos < end
    } else {
        pos >= start_a || pos < end
    }
}

fn compute_new_edges_seed(
    tile_id: usize,
    n: usize,
    pm: &crate::intgeom::patch::PatchMatch,
    m: usize,
) -> Vec<EdgeInfo> {
    let ccw_pos = (pm.start_a + pm.len) % n;
    let seg_len_old = n - pm.len;
    let seg_len_new = m - pm.len;
    let mut edges = Vec::with_capacity(seg_len_old + seg_len_new);
    for i in 0..seg_len_old {
        edges.push(EdgeInfo {
            tile_id,
            tile_offset: (ccw_pos + i) % n,
        });
    }
    edges.push(EdgeInfo {
        tile_id: pm.tile_id,
        tile_offset: pm.start_b,
    });
    for k in 1..seg_len_new {
        edges.push(EdgeInfo {
            tile_id: pm.tile_id,
            tile_offset: (pm.start_b + pm.len + k) % m,
        });
    }
    edges
}

fn compute_new_edges_growing(
    edges: &[EdgeInfo],
    n: usize,
    pm: &crate::intgeom::patch::PatchMatch,
    m: usize,
) -> Vec<EdgeInfo> {
    let ccw_pos = (pm.start_a + pm.len) % n;
    let seg_len_old = n - pm.len;
    let seg_len_new = m - pm.len;
    let mut new_edges = Vec::with_capacity(seg_len_old + seg_len_new);
    for i in 0..seg_len_old {
        new_edges.push(edges[(ccw_pos + i) % n]);
    }
    new_edges.push(EdgeInfo {
        tile_id: pm.tile_id,
        tile_offset: pm.start_b,
    });
    for k in 1..seg_len_new {
        new_edges.push(EdgeInfo {
            tile_id: pm.tile_id,
            tile_offset: (pm.start_b + pm.len + k) % m,
        });
    }
    new_edges
}

struct BfsState {
    angles: Vec<i8>,
    edges: Vec<EdgeInfo>,
    gap_start: u8,
    gap_len: u8,
    central_tile_id: usize,
    central_n: usize,
}

pub struct OpenNeighborhoodIndex<T: IsComplex> {
    tileset: Arc<TileSet<T>>,
    entries: Vec<OpenNeighborhoodType>,
}

fn count_covered_gap_edges(
    gap_len: u8,
    start_a: usize,
    match_len: usize,
    boundary_n: usize,
) -> usize {
    (0..gap_len as usize)
        .filter(|&p| in_consumed_range(p, start_a, match_len, boundary_n))
        .count()
}

impl<T: IsComplex + IsRingOrField + Units> OpenNeighborhoodIndex<T> {
    pub fn new(tileset: Arc<TileSet<T>>) -> Self {
        let match_index = Arc::new(MatchTypeIndex::new(Arc::clone(&tileset)));
        let mut all_types: Vec<OpenNeighborhoodType> = Vec::new();
        let mut seen = SeenIndex::new();
        let mut current_level: Vec<BfsState> = Vec::new();

        for tile_id in 0..tileset.num_tiles() {
            let n = tileset.rat(tile_id).seq().len();
            assert!(n <= 32, "tile edge count exceeds 32");
            let seed_angles: Vec<i8> = tileset.rat(tile_id).seq().to_vec();
            let seed_edges: Vec<EdgeInfo> = (0..n)
                .map(|i| EdgeInfo {
                    tile_id,
                    tile_offset: i,
                })
                .collect();

            let seed_candidates =
                GrowingPatch::compute_all_candidates(&match_index, &seed_angles, &seed_edges);
            eprintln!(
                "  tile {}: {} edges, {} seed candidates",
                tile_id,
                n,
                seed_candidates.iter().map(|v| v.len()).sum::<usize>()
            );

            let mut seed_types = 0;
            for matches in &seed_candidates {
                for pm in matches {
                    let m = tileset.rat(pm.tile_id).seq().len();
                    let new_angles = match compute_glue_angles::<T>(&seed_angles, pm, &tileset) {
                        Some(a) => a,
                        None => continue,
                    };
                    let snake = match Snake::<T>::try_from(new_angles.as_slice()) {
                        Ok(s) => s,
                        Err(_) => continue,
                    };
                    if !snake.is_closed() {
                        continue;
                    }

                    let new_edges = compute_new_edges_seed(tile_id, n, pm, m);
                    let covered_count = pm.len;
                    let new_gap_start = ((pm.start_a + covered_count) % n) as u8;
                    let new_gap_len = (n - covered_count) as u8;
                    if new_gap_len == 0 {
                        continue;
                    }

                    let (ca, ce, _rot) = canonicalize(&new_angles, &new_edges);
                    if !seen.insert(&all_types, tile_id, &ca, &ce, new_gap_start, new_gap_len) {
                        continue;
                    }

                    all_types.push(OpenNeighborhoodType {
                        id: 0,
                        central_tile_id: tile_id,
                        gap_start: new_gap_start,
                        gap_len: new_gap_len,
                        angles: ca,
                        edges: ce,
                        inner_chains: vec![],
                    });
                    seed_types += 1;

                    current_level.push(BfsState {
                        angles: new_angles,
                        edges: new_edges,
                        gap_start: new_gap_start,
                        gap_len: new_gap_len,
                        central_tile_id: tile_id,
                        central_n: n,
                    });
                }
            }
            eprintln!(
                "  tile {} seeds done: {} types, {} queue, {} seen",
                tile_id,
                seed_types,
                current_level.len(),
                seen.by_hash.len()
            );
        }

        let mut explored: usize = 0;
        while !current_level.is_empty() {
            let mut next_level: Vec<BfsState> = Vec::new();
            for state in current_level {
                if state.gap_len <= 1 {
                    continue;
                }
                let boundary_n = state.angles.len();
                explored += 1;
                if explored.is_multiple_of(5000) {
                    eprintln!(
                        "  explored {}, next {}, seen {}, types {}",
                        explored,
                        next_level.len(),
                        seen.by_hash.len(),
                        all_types.len(),
                    );
                }

                let candidates =
                    GrowingPatch::compute_all_candidates(&match_index, &state.angles, &state.edges);

                for matches in &candidates {
                    for pm in matches {
                        if !in_consumed_range(0, pm.start_a, pm.len, boundary_n) {
                            continue;
                        }

                        let covered_count =
                            count_covered_gap_edges(state.gap_len, pm.start_a, pm.len, boundary_n);
                        if covered_count == 0 {
                            continue;
                        }

                        let m = tileset.rat(pm.tile_id).seq().len();
                        let new_angles = match compute_glue_angles::<T>(&state.angles, pm, &tileset)
                        {
                            Some(a) => a,
                            None => continue,
                        };
                        let snake = match Snake::<T>::try_from(new_angles.as_slice()) {
                            Ok(s) => s,
                            Err(_) => continue,
                        };
                        if !snake.is_closed() {
                            continue;
                        }

                        let new_edges = compute_new_edges_growing(&state.edges, boundary_n, pm, m);
                        let new_gap_start =
                            (state.gap_start as usize + covered_count) % state.central_n;
                        let new_gap_len = state.gap_len as usize - covered_count;
                        if new_gap_len == 0 {
                            continue;
                        }

                        let (ca, ce, _rot) = canonicalize(&new_angles, &new_edges);
                        if !seen.insert(
                            &all_types,
                            state.central_tile_id,
                            &ca,
                            &ce,
                            new_gap_start as u8,
                            new_gap_len as u8,
                        ) {
                            continue;
                        }

                        all_types.push(OpenNeighborhoodType {
                            id: 0,
                            central_tile_id: state.central_tile_id,
                            gap_start: new_gap_start as u8,
                            gap_len: new_gap_len as u8,
                            angles: ca,
                            edges: ce,
                            inner_chains: vec![],
                        });

                        next_level.push(BfsState {
                            angles: new_angles,
                            edges: new_edges,
                            gap_start: new_gap_start as u8,
                            gap_len: new_gap_len as u8,
                            central_tile_id: state.central_tile_id,
                            central_n: state.central_n,
                        });
                    }
                }
            }
            eprintln!(
                "  level done: explored {}, next_level {}, total types {}",
                explored,
                next_level.len(),
                all_types.len(),
            );
            current_level = next_level;
        }

        for nt in &mut all_types {
            if nt.inner_chains.is_empty() {
                nt.inner_chains = vec![vec![]; nt.angles.len()];
            }
        }

        for (i, nt) in all_types.iter_mut().enumerate() {
            nt.id = i + 1;
        }

        eprintln!("  total: {} types", all_types.len());

        OpenNeighborhoodIndex {
            tileset,
            entries: all_types,
        }
    }

    pub fn entries(&self) -> &[OpenNeighborhoodType] {
        &self.entries
    }

    pub fn num_types(&self) -> usize {
        self.entries.len()
    }

    pub fn tileset(&self) -> &Arc<TileSet<T>> {
        &self.tileset
    }

    pub fn validate(&self) -> Vec<usize> {
        let match_index = Arc::new(MatchTypeIndex::new(Arc::clone(&self.tileset)));
        let mut invalid = Vec::new();
        for nhood in &self.entries {
            if !self.validate_one(nhood, &match_index) {
                invalid.push(nhood.id);
            }
        }
        invalid
    }

    fn validate_one(
        &self,
        nhood: &OpenNeighborhoodType,
        match_index: &Arc<MatchTypeIndex<T>>,
    ) -> bool {
        if nhood.gap_len == 0 {
            return false;
        }
        if nhood.angles.len() != nhood.edges.len() || nhood.angles.is_empty() {
            return false;
        }

        let candidates =
            GrowingPatch::compute_all_candidates(match_index, &nhood.angles, &nhood.edges);
        let patch = match GrowingPatch::from_parts(
            Arc::clone(match_index),
            nhood.angles.clone(),
            nhood.edges.clone(),
            candidates,
            nhood.inner_chains.clone(),
        ) {
            Some(p) => p,
            None => return false,
        };

        if patch.boundary_len() != nhood.angles.len() {
            return false;
        }
        if patch.edges().len() != nhood.edges.len() {
            return false;
        }
        for (i, (a, b)) in patch.angles().iter().zip(nhood.angles.iter()).enumerate() {
            if a != b {
                eprintln!("  angle mismatch at {}: {} vs {}", i, a, b);
                return false;
            }
        }
        for (i, (a, b)) in patch.edges().iter().zip(nhood.edges.iter()).enumerate() {
            if a != b {
                eprintln!(
                    "  edge mismatch at {}: {}.{} vs {}.{}",
                    i, a.tile_id, a.tile_offset, b.tile_id, b.tile_offset
                );
                return false;
            }
        }

        let central_count = patch
            .edges()
            .iter()
            .filter(|ei| ei.tile_id == nhood.central_tile_id)
            .count();
        if central_count < nhood.gap_len as usize {
            eprintln!(
                "  central tile edge count {} < gap_len {}",
                central_count, nhood.gap_len
            );
            return false;
        }
        true
    }

    pub fn write_collection(&self, out: &mut impl std::io::Write) -> std::io::Result<()> {
        for nhood in &self.entries {
            let n = nhood.angles.len();
            let angles_str = nhood
                .angles
                .iter()
                .map(|a| a.to_string())
                .collect::<Vec<_>>()
                .join(" ");
            let edges_str = nhood
                .edges
                .iter()
                .map(|e| format!("{}.{}", e.tile_id, e.tile_offset))
                .collect::<Vec<_>>()
                .join(" ");
            let inner_chains_str = nhood
                .inner_chains
                .iter()
                .map(|chain| {
                    if chain.is_empty() {
                        "-".to_string()
                    } else {
                        chain
                            .iter()
                            .map(|e| format!("{}.{}", e.tile_id, e.tile_offset))
                            .collect::<Vec<_>>()
                            .join(",")
                    }
                })
                .collect::<Vec<_>>()
                .join(" ");

            writeln!(
                out,
                "NTYPE {} {} {} {} {} {} {} {}",
                nhood.id,
                nhood.central_tile_id,
                nhood.gap_start,
                nhood.gap_len,
                n,
                angles_str,
                edges_str,
                inner_chains_str,
            )?;
        }
        Ok(())
    }

    pub fn parse_file(tileset: Arc<TileSet<T>>, input: &str) -> Result<Self, String> {
        let mut entries = Vec::new();
        for line in input.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.is_empty() || parts[0] != "NTYPE" {
                continue;
            }
            if parts.len() < 6 {
                return Err(format!("NTYPE line too short: {}", line));
            }

            let id: usize = parts[1]
                .parse()
                .map_err(|_| format!("bad id: {}", parts[1]))?;
            let central_tile_id: usize = parts[2]
                .parse()
                .map_err(|_| format!("bad central_tile_id: {}", parts[2]))?;
            let gap_start: u8 = parts[3]
                .parse()
                .map_err(|_| format!("bad gap_start: {}", parts[3]))?;
            let gap_len: u8 = parts[4]
                .parse()
                .map_err(|_| format!("bad gap_len: {}", parts[4]))?;
            let n: usize = parts[5]
                .parse()
                .map_err(|_| format!("bad n: {}", parts[5]))?;

            if parts.len() < 6 + n {
                return Err(format!(
                    "NTYPE {} expected {} angles, got {} fields",
                    id,
                    n,
                    parts.len() - 6
                ));
            }
            let angles: Vec<i8> = parts[6..6 + n]
                .iter()
                .map(|s| s.parse::<i8>())
                .collect::<Result<Vec<_>, _>>()
                .map_err(|e| format!("bad angle: {}", e))?;

            let edges_start = 6 + n;
            if parts.len() < edges_start + n {
                return Err(format!(
                    "NTYPE {} expected {} edges, got {} fields",
                    id,
                    n,
                    parts.len() - edges_start
                ));
            }
            let edges: Vec<EdgeInfo> = parts[edges_start..edges_start + n]
                .iter()
                .map(|s| parse_edge_info(s))
                .collect::<Result<Vec<_>, _>>()?;

            let inner_start = edges_start + n;
            let inner_chains: Vec<Vec<EdgeInfo>> = if inner_start < parts.len() {
                parts[inner_start..]
                    .iter()
                    .map(|s| {
                        if *s == "-" {
                            Ok(vec![])
                        } else {
                            s.split(',')
                                .map(parse_edge_info)
                                .collect::<Result<Vec<_>, _>>()
                        }
                    })
                    .collect::<Result<Vec<_>, _>>()?
            } else {
                vec![vec![]; n]
            };

            if inner_chains.len() != n {
                return Err(format!(
                    "NTYPE {} expected {} inner chains, got {}",
                    id,
                    n,
                    inner_chains.len()
                ));
            }

            entries.push(OpenNeighborhoodType {
                id,
                central_tile_id,
                gap_start,
                gap_len,
                angles,
                edges,
                inner_chains,
            });
        }
        Ok(OpenNeighborhoodIndex { tileset, entries })
    }
}

fn parse_edge_info(s: &str) -> Result<EdgeInfo, String> {
    let parts: Vec<&str> = s.split('.').collect();
    if parts.len() != 2 {
        return Err(format!("bad edge info '{}', expected tile_id.offset", s));
    }
    let tile_id: usize = parts[0]
        .parse()
        .map_err(|_| format!("bad tile_id: {}", parts[0]))?;
    let tile_offset: usize = parts[1]
        .parse()
        .map_err(|_| format!("bad tile_offset: {}", parts[1]))?;
    Ok(EdgeInfo {
        tile_id,
        tile_offset,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic::ZZ12;
    use crate::intgeom::rat::Rat;
    use crate::intgeom::tiles;
    use std::collections::BTreeMap;

    #[test]
    fn test_square_neighborhood_types() {
        let sq: crate::intgeom::snake::Snake<ZZ12> = tiles::square();
        let rat = Rat::try_from(&sq).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = OpenNeighborhoodIndex::new(ts);
        eprintln!("square: {} neighborhood types", idx.num_types());

        let mut by_gap_len: BTreeMap<usize, Vec<&OpenNeighborhoodType>> = BTreeMap::new();
        for nhood in idx.entries() {
            by_gap_len.entry(nhood.gap_len()).or_default().push(nhood);
        }
        for (gl, types) in &by_gap_len {
            eprintln!("  gap_len={}: {} types", gl, types.len());
        }

        let invalid = idx.validate();
        assert!(
            invalid.is_empty(),
            "found {} invalid: {:?}",
            invalid.len(),
            &invalid[..invalid.len().min(10)],
        );
    }

    #[test]
    fn test_hex_neighborhood_summary() {
        let hex: crate::intgeom::snake::Snake<ZZ12> = tiles::hexagon();
        let rat = Rat::try_from(&hex).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = OpenNeighborhoodIndex::new(ts);

        let mut by_gap: BTreeMap<usize, BTreeMap<usize, usize>> = BTreeMap::new();
        for nhood in idx.entries() {
            *by_gap
                .entry(nhood.gap_len())
                .or_default()
                .entry(nhood.angles.len())
                .or_insert(0) += 1;
        }
        eprintln!("hex: {} total types", idx.num_types());
        for (gl, lens) in &by_gap {
            for (bl, cnt) in lens {
                eprintln!("  gap_len={} boundary_len={}: {} types", gl, bl, cnt);
            }
        }

        let invalid = idx.validate();
        assert!(
            invalid.is_empty(),
            "found {} invalid: {:?}",
            invalid.len(),
            &invalid[..invalid.len().min(10)],
        );
    }

    #[test]
    fn test_square_neighborhood_roundtrip() {
        let sq: crate::intgeom::snake::Snake<ZZ12> = tiles::square();
        let rat = Rat::try_from(&sq).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = OpenNeighborhoodIndex::new(Arc::clone(&ts));

        let mut buf = Vec::new();
        idx.write_collection(&mut buf).unwrap();
        let output = String::from_utf8(buf).unwrap();

        let idx2 = OpenNeighborhoodIndex::parse_file(Arc::clone(&ts), &output).unwrap();
        assert_eq!(idx.entries().len(), idx2.entries().len());

        for (a, b) in idx.entries().iter().zip(idx2.entries().iter()) {
            assert_eq!(a.id, b.id);
            assert_eq!(a.central_tile_id, b.central_tile_id);
            assert_eq!(a.gap_start, b.gap_start);
            assert_eq!(a.gap_len, b.gap_len);
            assert_eq!(a.angles, b.angles);
            assert_eq!(a.edges, b.edges);
            assert_eq!(a.inner_chains, b.inner_chains);
        }
    }

    #[test]
    fn test_spectre_neighborhood_summary() {
        let sp: crate::intgeom::snake::Snake<ZZ12> = tiles::spectre();
        let rat = Rat::try_from(&sp).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = OpenNeighborhoodIndex::new(ts);

        let mut by_gap: BTreeMap<usize, BTreeMap<usize, usize>> = BTreeMap::new();
        for nhood in idx.entries() {
            *by_gap
                .entry(nhood.gap_len())
                .or_default()
                .entry(nhood.angles.len())
                .or_insert(0) += 1;
        }
        eprintln!("spectre: {} total types", idx.num_types());
        let mut by_gap_total: BTreeMap<usize, usize> = BTreeMap::new();
        for nhood in idx.entries() {
            *by_gap_total.entry(nhood.gap_len()).or_insert(0) += 1;
        }
        for (gl, cnt) in &by_gap_total {
            eprintln!("  gap_len={}: {} types", gl, cnt);
        }

        let invalid = idx.validate();
        assert!(
            invalid.is_empty(),
            "found {} invalid: {:?}",
            invalid.len(),
            &invalid[..invalid.len().min(10)],
        );
    }
}
