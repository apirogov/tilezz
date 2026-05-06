use std::hash::{Hash, Hasher};
use std::sync::Arc;

use rustc_hash::{FxHashMap, FxHasher};

use crate::cyclotomic::{IsComplex, IsRingOrField, Units};
use crate::intgeom::matchtypes::MatchTypeIndex;
use crate::intgeom::patch::{
    compute_glue_angles, cyclic_range_contains, full_vertex_type_from, update_inner_chains,
    EdgeInfo, GrowingPatch, PatchMatch, VertexType,
};
use crate::intgeom::snake::Snake;
use crate::intgeom::tileset::TileSet;

const MAX_BOUNDARY_SIZE: usize = 80;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NeighborhoodType {
    pub id: usize,
    pub central_tile_id: usize,
    pub central_anchor_edge: usize,
    pub gap_len: u8,
    pub context_vertex_types: Vec<VertexType>,
    pub angles: Vec<i8>,
    pub edges: Vec<EdgeInfo>,
    pub inner_chains: Vec<Vec<EdgeInfo>>,
    pub gap_start: usize,
}

impl NeighborhoodType {
    pub fn gap_len(&self) -> usize {
        self.gap_len as usize
    }
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

fn count_covered_gap_edges(gap_len: usize, start_a: usize, match_len: usize, n: usize) -> usize {
    (0..gap_len)
        .filter(|&p| in_consumed_range(p, start_a, match_len, n))
        .count()
}

fn extract_covered_vertex_types(
    edges: &[EdgeInfo],
    inner_chains: &[Vec<EdgeInfo>],
    gap_start: usize,
    gap_len: usize,
) -> Vec<VertexType> {
    let n = edges.len();
    let covered_len = n - gap_len;
    if covered_len < 2 {
        return Vec::new();
    }
    let num_inner = covered_len - 1;
    let mut result = Vec::with_capacity(num_inner);
    for k in 1..=num_inner {
        let pos = (gap_start + n - k) % n;
        result.push(full_vertex_type_from(edges, inner_chains, pos));
    }
    result
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

struct CanonicalState {
    angles: Vec<i8>,
    edges: Vec<EdgeInfo>,
    inner_chains: Vec<Vec<EdgeInfo>>,
    gap_start: usize,
}

fn canonicalize_state(
    angles: &[i8],
    edges: &[EdgeInfo],
    inner_chains: &[Vec<EdgeInfo>],
) -> CanonicalState {
    let rot = find_best_rotation(angles, edges);
    let n = angles.len();
    let mut a = angles.to_vec();
    let mut e = edges.to_vec();
    let mut ic = inner_chains.to_vec();
    a.rotate_left(rot);
    e.rotate_left(rot);
    ic.rotate_left(rot);
    CanonicalState {
        angles: a,
        edges: e,
        inner_chains: ic,
        gap_start: rot % n,
    }
}

fn dedup_hash(
    central_tile_id: usize,
    ca: &[i8],
    ce: &[EdgeInfo],
    ci: &[Vec<EdgeInfo>],
    gap_len: u8,
) -> u64 {
    let mut h = FxHasher::default();
    central_tile_id.hash(&mut h);
    ca.hash(&mut h);
    ce.hash(&mut h);
    ci.hash(&mut h);
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
        all_types: &[NeighborhoodType],
        central_tile_id: usize,
        canon: &CanonicalState,
        gap_len: u8,
    ) -> bool {
        let hash = dedup_hash(
            central_tile_id,
            &canon.angles,
            &canon.edges,
            &canon.inner_chains,
            gap_len,
        );
        let indices = self.by_hash.entry(hash).or_default();
        for &idx in indices.iter() {
            let t = &all_types[idx];
            if t.central_tile_id == central_tile_id
                && t.gap_len == gap_len
                && t.angles == canon.angles
                && t.edges == canon.edges
                && t.inner_chains == canon.inner_chains
            {
                return false;
            }
        }
        indices.push(all_types.len());
        true
    }
}

fn compute_new_edges_seed(tile_id: usize, n: usize, pm: &PatchMatch, m: usize) -> Vec<EdgeInfo> {
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
    pm: &PatchMatch,
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

struct NtBfsState {
    angles: Vec<i8>,
    edges: Vec<EdgeInfo>,
    inner_chains: Vec<Vec<EdgeInfo>>,
    gap_len: u8,
    central_tile_id: usize,
}

pub struct NeighborhoodIndex<T: IsComplex> {
    tileset: Arc<TileSet<T>>,
    entries: Vec<NeighborhoodType>,
}

impl<T: IsComplex + IsRingOrField + Units> NeighborhoodIndex<T> {
    pub fn new(tileset: Arc<TileSet<T>>) -> Self {
        let match_index = Arc::new(MatchTypeIndex::new(Arc::clone(&tileset)));
        let mut all_types: Vec<NeighborhoodType> = Vec::new();
        let mut seen = SeenIndex::new();
        let mut current_level: Vec<NtBfsState> = Vec::new();

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
            let seed_inner: Vec<Vec<EdgeInfo>> = vec![vec![]; n];

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

                    let new_n = new_angles.len();
                    let new_edges = compute_new_edges_seed(tile_id, n, pm, m);
                    let new_inner = update_inner_chains(&seed_inner, &seed_edges, pm, new_n);

                    let gap_len = n - pm.len;
                    if gap_len == 0 || gap_len > u8::MAX as usize {
                        continue;
                    }

                    let canon = canonicalize_state(&new_angles, &new_edges, &new_inner);
                    if !seen.insert(&all_types, tile_id, &canon, gap_len as u8) {
                        continue;
                    }

                    let context_vts = extract_covered_vertex_types(
                        &canon.edges,
                        &canon.inner_chains,
                        canon.gap_start,
                        gap_len,
                    );
                    let central_anchor = canon.edges[canon.gap_start].tile_offset;

                    all_types.push(NeighborhoodType {
                        id: 0,
                        central_tile_id: tile_id,
                        central_anchor_edge: central_anchor,
                        gap_len: gap_len as u8,
                        context_vertex_types: context_vts,
                        angles: canon.angles.clone(),
                        edges: canon.edges.clone(),
                        inner_chains: canon.inner_chains.clone(),
                        gap_start: canon.gap_start,
                    });
                    seed_types += 1;

                    current_level.push(NtBfsState {
                        angles: new_angles,
                        edges: new_edges,
                        inner_chains: new_inner,
                        gap_len: gap_len as u8,
                        central_tile_id: tile_id,
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
        let mut level_num: usize = 0;
        while !current_level.is_empty() {
            let mut next_level: Vec<NtBfsState> = Vec::new();
            let mut stat_candidates: usize = 0;
            let mut stat_touching: usize = 0;
            let mut stat_glue: usize = 0;
            let mut stat_snake: usize = 0;
            let mut stat_dedup: usize = 0;
            let mut stat_update: usize = 0;
            let mut stat_max_bsize: usize = 0;

            for state in &current_level {
                let boundary_n = state.angles.len();
                if boundary_n > MAX_BOUNDARY_SIZE {
                    continue;
                }
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

                let all_candidates =
                    GrowingPatch::compute_all_candidates(&match_index, &state.angles, &state.edges);

                stat_candidates += all_candidates.iter().map(|v| v.len()).sum::<usize>();

                let cw_prev = boundary_n - 1;

                struct DeferredUpdate {
                    pm: PatchMatch,
                    new_angles: Vec<i8>,
                    new_n: usize,
                }
                let mut deferred_updates: Vec<DeferredUpdate> = Vec::new();
                let mut had_close_success = false;

                for pm in all_candidates.iter().flatten() {
                    let covers_ccw = cyclic_range_contains(pm.start_a, pm.len, 0, boundary_n);
                    let covers_cw = cyclic_range_contains(pm.start_a, pm.len, cw_prev, boundary_n);
                    if !covers_ccw && !covers_cw {
                        continue;
                    }
                    stat_touching += 1;

                    let covered_count = count_covered_gap_edges(
                        state.gap_len as usize,
                        pm.start_a,
                        pm.len,
                        boundary_n,
                    );
                    let new_gap_len = state.gap_len as usize - covered_count;
                    if new_gap_len == 0 || new_gap_len > u8::MAX as usize {
                        continue;
                    }
                    let is_update = covered_count == 0;

                    let m = tileset.rat(pm.tile_id).seq().len();
                    let new_angles = match compute_glue_angles::<T>(&state.angles, pm, &tileset) {
                        Some(a) => a,
                        None => continue,
                    };
                    stat_glue += 1;

                    let snake = match Snake::<T>::try_from(new_angles.as_slice()) {
                        Ok(s) => s,
                        Err(_) => continue,
                    };
                    if !snake.is_closed() {
                        continue;
                    }
                    stat_snake += 1;

                    let new_n = new_angles.len();
                    if new_n > MAX_BOUNDARY_SIZE {
                        continue;
                    }

                    if is_update {
                        deferred_updates.push(DeferredUpdate {
                            pm: pm.clone(),
                            new_angles,
                            new_n,
                        });
                        continue;
                    }

                    let new_edges = compute_new_edges_growing(&state.edges, boundary_n, pm, m);
                    let new_inner =
                        update_inner_chains(&state.inner_chains, &state.edges, pm, new_n);

                    let canon = canonicalize_state(&new_angles, &new_edges, &new_inner);
                    if !seen.insert(&all_types, state.central_tile_id, &canon, new_gap_len as u8) {
                        continue;
                    }
                    stat_dedup += 1;
                    had_close_success = true;

                    let context_vts = extract_covered_vertex_types(
                        &canon.edges,
                        &canon.inner_chains,
                        canon.gap_start,
                        new_gap_len,
                    );
                    let central_anchor = canon.edges[canon.gap_start].tile_offset;

                    all_types.push(NeighborhoodType {
                        id: 0,
                        central_tile_id: state.central_tile_id,
                        central_anchor_edge: central_anchor,
                        gap_len: new_gap_len as u8,
                        context_vertex_types: context_vts,
                        angles: canon.angles.clone(),
                        edges: canon.edges.clone(),
                        inner_chains: canon.inner_chains.clone(),
                        gap_start: canon.gap_start,
                    });

                    next_level.push(NtBfsState {
                        angles: new_angles,
                        edges: new_edges,
                        inner_chains: new_inner,
                        gap_len: new_gap_len as u8,
                        central_tile_id: state.central_tile_id,
                    });
                }

                if !had_close_success && state.gap_len >= 2 {
                    for du in deferred_updates {
                        stat_update += 1;
                        let m = tileset.rat(du.pm.tile_id).seq().len();
                        let new_edges =
                            compute_new_edges_growing(&state.edges, boundary_n, &du.pm, m);
                        let new_inner = update_inner_chains(
                            &state.inner_chains,
                            &state.edges,
                            &du.pm,
                            du.new_n,
                        );

                        let canon = canonicalize_state(&du.new_angles, &new_edges, &new_inner);
                        if !seen.insert(&all_types, state.central_tile_id, &canon, state.gap_len) {
                            continue;
                        }
                        stat_dedup += 1;

                        let context_vts = extract_covered_vertex_types(
                            &canon.edges,
                            &canon.inner_chains,
                            canon.gap_start,
                            state.gap_len as usize,
                        );
                        let central_anchor = canon.edges[canon.gap_start].tile_offset;

                        all_types.push(NeighborhoodType {
                            id: 0,
                            central_tile_id: state.central_tile_id,
                            central_anchor_edge: central_anchor,
                            gap_len: state.gap_len,
                            context_vertex_types: context_vts,
                            angles: canon.angles.clone(),
                            edges: canon.edges.clone(),
                            inner_chains: canon.inner_chains.clone(),
                            gap_start: canon.gap_start,
                        });

                        next_level.push(NtBfsState {
                            angles: du.new_angles,
                            edges: new_edges,
                            inner_chains: new_inner,
                            gap_len: state.gap_len,
                            central_tile_id: state.central_tile_id,
                        });
                    }
                }

                stat_max_bsize = stat_max_bsize.max(boundary_n);
            }
            level_num += 1;
            eprintln!(
                "  level {}: explored {}, next {}, types {} | cands={} touch={} glue={} snake={} update={} new={} maxbs={}",
                level_num,
                explored,
                next_level.len(),
                all_types.len(),
                stat_candidates,
                stat_touching,
                stat_glue,
                stat_snake,
                stat_update,
                stat_dedup,
                stat_max_bsize,
            );
            current_level = next_level;
        }

        for (i, nt) in all_types.iter_mut().enumerate() {
            nt.id = i + 1;
        }

        eprintln!("  total: {} types", all_types.len());

        NeighborhoodIndex {
            tileset,
            entries: all_types,
        }
    }

    pub fn entries(&self) -> &[NeighborhoodType] {
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

    fn validate_one(&self, nhood: &NeighborhoodType, match_index: &Arc<MatchTypeIndex<T>>) -> bool {
        if nhood.gap_len == 0 {
            return false;
        }
        if nhood.angles.len() != nhood.edges.len() || nhood.angles.is_empty() {
            return false;
        }
        if nhood.inner_chains.len() != nhood.angles.len() {
            return false;
        }

        let patch = match GrowingPatch::from_parts(
            Arc::clone(match_index),
            nhood.angles.clone(),
            nhood.edges.clone(),
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

        let n = nhood.edges.len();
        let gs = nhood.gap_start;
        let gl = nhood.gap_len as usize;
        let covered_len = n - gl;
        if covered_len == 0 {
            eprintln!("  no covered edges");
            return false;
        }

        let expected_vts = extract_covered_vertex_types(&nhood.edges, &nhood.inner_chains, gs, gl);
        if expected_vts != nhood.context_vertex_types {
            eprintln!("  context vertex types mismatch");
            return false;
        }

        if nhood.edges[gs].tile_offset != nhood.central_anchor_edge {
            eprintln!(
                "  central_anchor_edge mismatch: found {} vs claimed {}",
                nhood.edges[gs].tile_offset, nhood.central_anchor_edge
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
            let vts_str = nhood
                .context_vertex_types
                .iter()
                .map(format_vtype)
                .collect::<Vec<_>>()
                .join(" ");

            writeln!(
                out,
                "NTYPE {} {} {} {} {} {} {} {} {} {}",
                nhood.id,
                nhood.central_tile_id,
                nhood.gap_start,
                nhood.gap_len,
                nhood.central_anchor_edge,
                n,
                angles_str,
                edges_str,
                inner_chains_str,
                vts_str,
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
            if parts.len() < 7 {
                return Err(format!("NTYPE line too short: {}", line));
            }

            let id: usize = parts[1]
                .parse()
                .map_err(|_| format!("bad id: {}", parts[1]))?;
            let central_tile_id: usize = parts[2]
                .parse()
                .map_err(|_| format!("bad central_tile_id: {}", parts[2]))?;
            let gap_start: usize = parts[3]
                .parse()
                .map_err(|_| format!("bad gap_start: {}", parts[3]))?;
            let gap_len: u8 = parts[4]
                .parse()
                .map_err(|_| format!("bad gap_len: {}", parts[4]))?;
            let central_anchor_edge: usize = parts[5]
                .parse()
                .map_err(|_| format!("bad central_anchor_edge: {}", parts[5]))?;
            let n: usize = parts[6]
                .parse()
                .map_err(|_| format!("bad n: {}", parts[6]))?;

            if parts.len() < 7 + n {
                return Err(format!(
                    "NTYPE {} expected {} angles, got {} fields",
                    id,
                    n,
                    parts.len() - 7
                ));
            }
            let angles: Vec<i8> = parts[7..7 + n]
                .iter()
                .map(|s| s.parse::<i8>())
                .collect::<Result<Vec<_>, _>>()
                .map_err(|e| format!("bad angle: {}", e))?;

            let edges_start = 7 + n;
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
            let inner_chains: Vec<Vec<EdgeInfo>> = if parts.len() >= inner_start + n {
                parts[inner_start..inner_start + n]
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

            let vts_start = inner_start + n;
            let context_vertex_types = if parts.len() > vts_start {
                parts[vts_start..]
                    .iter()
                    .map(|s| parse_vtype(s))
                    .collect::<Result<Vec<_>, _>>()?
            } else {
                extract_covered_vertex_types(&edges, &inner_chains, gap_start, gap_len as usize)
            };

            entries.push(NeighborhoodType {
                id,
                central_tile_id,
                central_anchor_edge,
                gap_len,
                context_vertex_types,
                angles,
                edges,
                inner_chains,
                gap_start,
            });
        }
        Ok(NeighborhoodIndex { tileset, entries })
    }
}

fn format_vtype(vt: &VertexType) -> String {
    let inner_str = if vt.inner.is_empty() {
        "-".to_string()
    } else {
        vt.inner
            .iter()
            .map(|e| format!("{}.{}", e.tile_id, e.tile_offset))
            .collect::<Vec<_>>()
            .join(",")
    };
    format!(
        "{}.{};{};{}.{}",
        vt.cw.tile_id, vt.cw.tile_offset, inner_str, vt.ccw.tile_id, vt.ccw.tile_offset
    )
}

fn parse_vtype(s: &str) -> Result<VertexType, String> {
    let parts: Vec<&str> = s.split(';').collect();
    if parts.len() != 3 {
        return Err(format!("bad vertex type '{}', expected cw;inner;ccw", s));
    }
    let cw = parse_edge_info(parts[0])?;
    let inner = if parts[1] == "-" {
        vec![]
    } else {
        parts[1]
            .split(',')
            .map(parse_edge_info)
            .collect::<Result<Vec<_>, _>>()?
    };
    let ccw = parse_edge_info(parts[2])?;
    Ok(VertexType { cw, inner, ccw })
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
        let idx = NeighborhoodIndex::new(ts);
        eprintln!("square: {} neighborhood types", idx.num_types());

        let mut by_gap_len: BTreeMap<usize, Vec<&NeighborhoodType>> = BTreeMap::new();
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
        let idx = NeighborhoodIndex::new(ts);

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
        let idx = NeighborhoodIndex::new(Arc::clone(&ts));

        let mut buf = Vec::new();
        idx.write_collection(&mut buf).unwrap();
        let output = String::from_utf8(buf).unwrap();

        let idx2 = NeighborhoodIndex::parse_file(Arc::clone(&ts), &output).unwrap();
        assert_eq!(idx.entries().len(), idx2.entries().len());

        for (a, b) in idx.entries().iter().zip(idx2.entries().iter()) {
            assert_eq!(a.id, b.id);
            assert_eq!(a.central_tile_id, b.central_tile_id);
            assert_eq!(a.gap_start, b.gap_start);
            assert_eq!(a.gap_len, b.gap_len);
            assert_eq!(a.central_anchor_edge, b.central_anchor_edge);
            assert_eq!(a.angles, b.angles);
            assert_eq!(a.edges, b.edges);
            assert_eq!(a.inner_chains, b.inner_chains);
            assert_eq!(a.context_vertex_types, b.context_vertex_types);
        }
    }

    #[test]
    fn test_spectre_neighborhood_summary() {
        let sp: crate::intgeom::snake::Snake<ZZ12> = tiles::spectre();
        let rat = Rat::try_from(&sp).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = NeighborhoodIndex::new(ts);

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
