use std::collections::{HashSet, VecDeque};
use std::sync::Arc;

use crate::cyclotomic::{IsComplex, IsRingOrField, Units};
use crate::intgeom::matchtypes::MatchTypeIndex;
use crate::intgeom::patch::{EdgeInfo, GrowingPatch};
use crate::intgeom::tileset::TileSet;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OpenNeighborhoodType {
    pub id: usize,
    pub central_tile_id: usize,
    pub gap_bitmask: u32,
    pub gap_positions: Vec<usize>,
    pub angles: Vec<i8>,
    pub edges: Vec<EdgeInfo>,
    pub inner_chains: Vec<Vec<EdgeInfo>>,
}

impl OpenNeighborhoodType {
    pub fn gap_len(&self) -> usize {
        self.gap_bitmask.count_ones() as usize
    }
}

type DedupKey = (
    usize,
    Vec<i8>,
    Vec<EdgeInfo>,
    Vec<Vec<EdgeInfo>>,
    Vec<usize>,
);

fn find_best_rotation(angles: &[i8], edges: &[EdgeInfo], inner_chains: &[Vec<EdgeInfo>]) -> usize {
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
                    std::cmp::Ordering::Equal => {
                        match inner_chains[idx_r].cmp(&inner_chains[idx_b]) {
                            std::cmp::Ordering::Less => {
                                best = r;
                                break;
                            }
                            std::cmp::Ordering::Greater => break,
                            std::cmp::Ordering::Equal => continue,
                        }
                    }
                },
            }
        }
    }
    best
}

fn canonicalize_boundary(
    angles: &[i8],
    edges: &[EdgeInfo],
    inner_chains: &[Vec<EdgeInfo>],
) -> (Vec<i8>, Vec<EdgeInfo>, Vec<Vec<EdgeInfo>>, usize) {
    let n = angles.len();
    if n == 0 {
        return (vec![], vec![], vec![], 0);
    }
    let rot = find_best_rotation(angles, edges, inner_chains);
    let mut a = angles.to_vec();
    let mut e = edges.to_vec();
    let mut ic = inner_chains.to_vec();
    a.rotate_left(rot);
    e.rotate_left(rot);
    ic.rotate_left(rot);
    (a, e, ic, rot)
}

fn rotate_positions(positions: &[usize], rot: usize, n: usize) -> Vec<usize> {
    let mut result: Vec<usize> = positions.iter().map(|&p| (p + n - rot) % n).collect();
    result.sort();
    result
}

fn has_one_contiguous_run(bitmask: u32, n: usize) -> bool {
    if bitmask == 0 {
        return false;
    }
    let mut transitions = 0;
    for i in 0..n {
        let curr = (bitmask >> i) & 1;
        let next = (bitmask >> ((i + 1) % n)) & 1;
        if curr != next {
            transitions += 1;
        }
    }
    transitions == 2
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

fn update_positions(
    positions: &[usize],
    start_a: usize,
    match_len: usize,
    old_n: usize,
) -> Vec<usize> {
    let ccw_pos = (start_a + match_len) % old_n;
    let mut new_positions = Vec::new();
    for &p in positions {
        if in_consumed_range(p, start_a, match_len, old_n) {
            continue;
        }
        let new_p = if p >= ccw_pos {
            p - ccw_pos
        } else {
            (old_n - ccw_pos) + p
        };
        new_positions.push(new_p);
    }
    new_positions.sort();
    new_positions
}

fn compute_gap_bitmask(edges: &[EdgeInfo], positions: &[usize], central_tile_id: usize) -> u32 {
    let mut bitmask = 0u32;
    for &p in positions {
        let ei = edges[p];
        if ei.tile_id == central_tile_id {
            bitmask |= 1u32 << ei.tile_offset;
        }
    }
    bitmask
}

pub struct OpenNeighborhoodIndex<T: IsComplex> {
    tileset: Arc<TileSet<T>>,
    entries: Vec<OpenNeighborhoodType>,
}

impl<T: IsComplex + IsRingOrField + Units> OpenNeighborhoodIndex<T> {
    pub fn new(tileset: Arc<TileSet<T>>) -> Self {
        let mut all_types: Vec<OpenNeighborhoodType> = Vec::new();
        let mut seen: HashSet<DedupKey> = HashSet::new();
        let mut queue: VecDeque<(GrowingPatch<T>, Vec<usize>, usize)> = VecDeque::new();

        for tile_id in 0..tileset.num_tiles() {
            let n = tileset.rat(tile_id).seq().len();
            assert!(n <= 32, "tile edge count exceeds 32");
            let initial_positions: Vec<usize> = (0..n).collect();

            let seed = GrowingPatch::new(Arc::clone(&tileset), tile_id);
            let seed_matches = seed.get_all_matches();
            eprintln!(
                "  tile {}: {} edges, {} seed matches",
                tile_id,
                n,
                seed_matches.len()
            );

            let mut seed_types = 0;
            for pm in seed_matches {
                let mut gp = seed.clone_for_mutation();
                if gp.add_tile(&pm).is_none() || !gp.is_growing() {
                    continue;
                }

                let new_positions = update_positions(&initial_positions, pm.start_a, pm.len, n);

                if new_positions.is_empty() {
                    continue;
                }

                let bitmask = compute_gap_bitmask(gp.edges(), &new_positions, tile_id);

                let (ca, ce, ci, rot) =
                    canonicalize_boundary(gp.angles(), gp.edges(), gp.inner_chains());
                let cp = rotate_positions(&new_positions, rot, ca.len());

                let key = (tile_id, ca.clone(), ce.clone(), ci.clone(), cp.clone());
                if !seen.insert(key) {
                    continue;
                }

                if has_one_contiguous_run(bitmask, n) {
                    all_types.push(OpenNeighborhoodType {
                        id: 0,
                        central_tile_id: tile_id,
                        gap_bitmask: bitmask,
                        gap_positions: cp,
                        angles: ca,
                        edges: ce,
                        inner_chains: ci,
                    });
                    seed_types += 1;
                    queue.push_back((gp, new_positions, tile_id));
                }
            }
            eprintln!(
                "  tile {} seeds done: {} queue, {} types, {} seen",
                tile_id,
                queue.len(),
                seed_types,
                seen.len()
            );
        }

        let mut explored: usize = 0;
        let max_states = 100_000;
        while let Some((gp, positions, tile_id)) = queue.pop_front() {
            if seen.len() > max_states {
                eprintln!(
                    "  WARNING: exceeded {} states, stopping BFS for tile {}",
                    max_states, tile_id
                );
                break;
            }
            let n = tileset.rat(tile_id).seq().len();
            let boundary_n = gp.boundary_len();
            let matches = gp.get_all_matches();
            explored += 1;
            if explored <= 10 || explored.is_multiple_of(1000) {
                eprintln!(
                    "  explored {}, queue {}, seen {}, types {}, matches {}",
                    explored,
                    queue.len(),
                    seen.len(),
                    all_types.len(),
                    matches.len()
                );
            }

            for pm in matches {
                let touches_central = positions
                    .iter()
                    .any(|&p| in_consumed_range(p, pm.start_a, pm.len, boundary_n));
                if !touches_central {
                    continue;
                }

                let new_positions = update_positions(&positions, pm.start_a, pm.len, boundary_n);

                if new_positions.is_empty() {
                    continue;
                }

                let mut gp2 = gp.clone_for_mutation();
                if gp2.add_tile(&pm).is_none() || !gp2.is_growing() {
                    continue;
                }

                let bitmask = compute_gap_bitmask(gp2.edges(), &new_positions, tile_id);

                let (ca, ce, ci, rot) =
                    canonicalize_boundary(gp2.angles(), gp2.edges(), gp2.inner_chains());
                let cp = rotate_positions(&new_positions, rot, ca.len());

                let key = (tile_id, ca.clone(), ce.clone(), ci.clone(), cp.clone());
                if !seen.insert(key) {
                    continue;
                }

                if has_one_contiguous_run(bitmask, n) {
                    all_types.push(OpenNeighborhoodType {
                        id: 0,
                        central_tile_id: tile_id,
                        gap_bitmask: bitmask,
                        gap_positions: cp,
                        angles: ca,
                        edges: ce,
                        inner_chains: ci,
                    });
                    queue.push_back((gp2, new_positions, tile_id));
                }
            }
        }

        for (i, nt) in all_types.iter_mut().enumerate() {
            nt.id = i + 1;
        }

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
        let n = self.tileset.rat(nhood.central_tile_id).seq().len();

        if nhood.gap_bitmask == 0 {
            return false;
        }

        if !has_one_contiguous_run(nhood.gap_bitmask, n) {
            return false;
        }

        if nhood.angles.len() != nhood.edges.len()
            || nhood.inner_chains.len() != nhood.angles.len()
            || nhood.angles.is_empty()
        {
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

        let boundary_n = patch.boundary_len();
        if boundary_n != nhood.angles.len() {
            return false;
        }

        if nhood.gap_positions.is_empty() {
            return false;
        }

        for &p in &nhood.gap_positions {
            if p >= boundary_n {
                return false;
            }
            let ei = patch.edges()[p];
            if ei.tile_id != nhood.central_tile_id {
                return false;
            }
            let bit = 1u32 << ei.tile_offset;
            if nhood.gap_bitmask & bit == 0 {
                return false;
            }
        }

        let observed_bitmask =
            compute_gap_bitmask(patch.edges(), &nhood.gap_positions, nhood.central_tile_id);
        if observed_bitmask != nhood.gap_bitmask {
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
            let gap_pos_str = nhood
                .gap_positions
                .iter()
                .map(|p| p.to_string())
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
                nhood.gap_bitmask,
                n,
                angles_str,
                edges_str,
                gap_pos_str,
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
            let gap_bitmask: u32 = parts[3]
                .parse()
                .map_err(|_| format!("bad gap_bitmask: {}", parts[3]))?;
            let n: usize = parts[4]
                .parse()
                .map_err(|_| format!("bad n: {}", parts[4]))?;

            if parts.len() < 5 + n {
                return Err(format!(
                    "NTYPE {} expected {} angles, got {} fields",
                    id,
                    n,
                    parts.len() - 5
                ));
            }
            let angles: Vec<i8> = parts[5..5 + n]
                .iter()
                .map(|s| s.parse::<i8>())
                .collect::<Result<Vec<_>, _>>()
                .map_err(|e| format!("bad angle: {}", e))?;

            let edges_start = 5 + n;
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

            let gap_len = gap_bitmask.count_ones() as usize;
            let gap_pos_start = edges_start + n;
            if parts.len() < gap_pos_start + gap_len {
                return Err(format!(
                    "NTYPE {} expected {} gap positions, got {} fields",
                    id,
                    gap_len,
                    parts.len() - gap_pos_start
                ));
            }
            let gap_positions: Vec<usize> = parts[gap_pos_start..gap_pos_start + gap_len]
                .iter()
                .map(|s| s.parse::<usize>())
                .collect::<Result<Vec<_>, _>>()
                .map_err(|e| format!("bad gap_position: {}", e))?;

            let inner_start = gap_pos_start + gap_len;
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
                gap_bitmask,
                gap_positions,
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

    #[test]
    fn test_has_one_contiguous_run() {
        assert!(has_one_contiguous_run(0b000111, 6));
        assert!(has_one_contiguous_run(0b110000, 6));
        assert!(has_one_contiguous_run(0b100001, 6));
        assert!(has_one_contiguous_run(0b000001, 6));
        assert!(has_one_contiguous_run(0b011111, 6));
        assert!(!has_one_contiguous_run(0b000000, 6));
        assert!(!has_one_contiguous_run(0b010010, 6));
        assert!(!has_one_contiguous_run(0b101010, 6));
    }

    #[test]
    fn test_hexagon_neighborhood_types() {
        let hex: crate::intgeom::snake::Snake<ZZ12> = tiles::hexagon();
        let rat = Rat::try_from(&hex).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = OpenNeighborhoodIndex::new(ts);
        assert!(idx.num_types() > 0, "should find some neighborhood types");
        eprintln!("hex: {} neighborhood types", idx.num_types());

        let invalid = idx.validate();
        assert!(
            invalid.is_empty(),
            "found {} invalid neighborhood types: {:?}",
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
            assert_eq!(a.gap_bitmask, b.gap_bitmask);
            assert_eq!(a.gap_positions, b.gap_positions);
            assert_eq!(a.angles, b.angles);
            assert_eq!(a.edges, b.edges);
            assert_eq!(a.inner_chains, b.inner_chains);
        }
    }
}
