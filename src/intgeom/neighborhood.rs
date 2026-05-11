use std::collections::VecDeque;
use std::sync::Arc;

use rustc_hash::FxHashMap;

use crate::cyclotomic::{IsComplex, IsRingOrField, Units};
use crate::intgeom::matchtypes::MatchTypeIndex;
use crate::intgeom::patch::{EdgeInfo, GrowingPatch, PatchMatch, VertexType};
use crate::intgeom::tileset::TileSet;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct NeighborhoodType {
    pub central_tile_id: usize,
    pub cw_anchor_on_central: usize,
    pub cw_anchor_on_context: usize,
    pub vt_seq: Vec<VertexType>,
}

pub const NT_CLOSED_ID: usize = 0;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum NtKind {
    Dead,
    Undead,
    Blessed,
    Free,
}

#[derive(Clone, Debug)]
pub struct NtTransition {
    pub src_id: usize,
    pub dst_id: usize,
    pub tile_id: usize,
    pub tile_offset: usize,
}

pub struct NeighborhoodIndex<T: IsComplex> {
    tileset: Arc<TileSet<T>>,
    entries: Vec<NeighborhoodType>,
    transitions: Vec<NtTransition>,
}

impl<T: IsComplex + IsRingOrField + Units> NeighborhoodIndex<T> {
    pub fn new(tileset: Arc<TileSet<T>>) -> Self {
        let match_index = Arc::new(MatchTypeIndex::new(Arc::clone(&tileset)));
        let mut entries: Vec<NeighborhoodType> = Vec::new();
        let mut transitions: Vec<NtTransition> = Vec::new();
        // IDs are 1-based so NT_CLOSED_ID = 0 is unambiguous: entries[i] has
        // id i + 1. `seen` and transitions store ids, not indices.
        let mut seen: FxHashMap<NeighborhoodType, usize> = FxHashMap::default();
        let mut queue: VecDeque<NeighborhoodType> = VecDeque::new();

        for id in 1..=match_index.num_types() {
            let mt = match_index.get(id);
            let mut patch = GrowingPatch::new(Arc::clone(match_index.tileset()), mt.tile_a);
            let pm = PatchMatch {
                start_a: mt.start_a,
                len: mt.len,
                start_b: mt.start_b,
                tile_id: mt.tile_b,
            };
            if patch.add_tile(&pm).is_none() {
                continue;
            }
            patch.normalize();
            let patch_n = patch.boundary_len();

            // Collect all ctx-ctx junctions with their positions, sorted by
            // boundary position. We then filter per third_pm to only those
            // incident with the match (= within the covered segment).
            let all_juncs: Vec<(usize, VertexType)> = (0..patch_n)
                .filter_map(|i| patch.junction_vertex_type_at(i).map(|vt| (i, vt)))
                .collect();

            for third_pm in patch.get_all_matches() {
                // get_all_matches returns edge-compatible glues but for
                // non-convex tiles (e.g. spectre) the central might still
                // collide with the existing two-tile patch. Try the glue on
                // a clone first; skip if it fails.
                let mut trial = patch.clone_for_mutation();
                if trial.add_tile(&third_pm).is_none() {
                    continue;
                }

                let central_tile_id = third_pm.tile_id;
                let cw_anchor_on_central = third_pm.start_b;

                // vt_seq is the junctions incident with the central match:
                // those at positions in [anchor, anchor + match_len], i.e.
                // CCW distance from the anchor is at most match_len.
                let filtered: Vec<&(usize, VertexType)> = all_juncs
                    .iter()
                    .filter(|(pos, _)| (pos + patch_n - third_pm.start_a) % patch_n <= third_pm.len)
                    .collect();
                if filtered.is_empty() {
                    continue;
                }
                let first_junc = filtered[0].0;
                let vt_seq: Vec<VertexType> = filtered.iter().map(|(_, vt)| vt.clone()).collect();
                // cw_anchor_on_context = CW distance from first_junc back to
                // the anchor (equivalently, CCW distance from anchor to
                // first_junc). The first VT sits at or CCW of the anchor, so
                // this distance is small and invariant under BFS growth
                // (which only adds edges CCW of the last VT, never between
                // the anchor and first_junc).
                let cw_anchor_on_context = (first_junc + patch_n - third_pm.start_a) % patch_n;
                let nt = NeighborhoodType {
                    central_tile_id,
                    cw_anchor_on_central,
                    cw_anchor_on_context,
                    vt_seq,
                };
                if seen.contains_key(&nt) {
                    continue;
                }
                // The trial.add_tile above verified the central glues into
                // the normalized two-tile patch. Also verify the abstract NT
                // reconstructs from its vt_seq (needed for non-convex tiles
                // where seed-gen and reconstruct can diverge).
                if !nt_is_valid(&nt, &match_index) {
                    continue;
                }
                let id = entries.len() + 1;
                seen.insert(nt.clone(), id);
                entries.push(nt.clone());
                queue.push_back(nt);
            }
        }

        let mut explored = 0usize;
        while let Some(state) = queue.pop_front() {
            explored += 1;
            let src_id = *seen.get(&state).expect("dequeued state must be in seen");
            for outcome in explore_one(&state, &match_index) {
                let dst_id = match outcome.kind {
                    OutcomeKind::Closed => NT_CLOSED_ID,
                    OutcomeKind::Open { nt: new_nt, .. } => {
                        if let Some(&id) = seen.get(&new_nt) {
                            id
                        } else {
                            let id = entries.len() + 1;
                            seen.insert(new_nt.clone(), id);
                            entries.push(new_nt.clone());
                            queue.push_back(new_nt);
                            id
                        }
                    }
                };
                transitions.push(NtTransition {
                    src_id,
                    dst_id,
                    tile_id: outcome.petal_pm.tile_id,
                    tile_offset: outcome.petal_pm.start_b,
                });
            }
        }

        eprintln!(
            "  total: {} types, {} transitions, {} explored",
            entries.len(),
            transitions.len(),
            explored,
        );

        NeighborhoodIndex {
            tileset,
            entries,
            transitions,
        }
    }

    pub fn entries(&self) -> &[NeighborhoodType] {
        &self.entries
    }

    pub fn num_types(&self) -> usize {
        self.entries.len()
    }

    pub fn transitions(&self) -> &[NtTransition] {
        &self.transitions
    }

    pub fn tileset(&self) -> &Arc<TileSet<T>> {
        &self.tileset
    }

    pub fn classify_all(&self) -> Vec<NtKind> {
        // Transition ids are 1-based; convert with `- 1` when indexing into
        // these per-entry vectors. NT_CLOSED_ID (= 0) is the closed sentinel
        // and never appears as a src_id.
        let n = self.entries.len();
        let mut succ_sets: Vec<Vec<usize>> = vec![vec![]; n];
        let mut has_outgoing = vec![false; n];

        for t in &self.transitions {
            debug_assert!(t.src_id != NT_CLOSED_ID, "closed cannot be a src");
            let src_idx = t.src_id - 1;
            has_outgoing[src_idx] = true;
            if t.dst_id != NT_CLOSED_ID {
                succ_sets[src_idx].push(t.dst_id - 1);
            }
        }

        let mut cursed = vec![false; n];
        for i in 0..n {
            cursed[i] = !has_outgoing[i];
        }

        let mut changed = true;
        while changed {
            changed = false;
            for i in 0..n {
                if cursed[i] {
                    continue;
                }
                let succs = &succ_sets[i];
                if succs.is_empty() {
                    continue;
                }
                if succs.iter().all(|&s| cursed[s]) {
                    cursed[i] = true;
                    changed = true;
                }
            }
        }

        let mut blessed = vec![false; n];
        let mut has_closed_or_blessed_succ = vec![false; n];

        for t in &self.transitions {
            let src_idx = t.src_id - 1;
            if t.dst_id == NT_CLOSED_ID {
                has_closed_or_blessed_succ[src_idx] = true;
            }
        }

        changed = true;
        while changed {
            changed = false;
            for i in 0..n {
                if blessed[i] || cursed[i] {
                    continue;
                }
                let mut any = false;
                let src_id_i = i + 1;
                let all_good = self
                    .transitions
                    .iter()
                    .filter(|t| t.src_id == src_id_i)
                    .all(|t| {
                        any = true;
                        t.dst_id == NT_CLOSED_ID || blessed[t.dst_id - 1]
                    });
                if any && all_good {
                    blessed[i] = true;
                    has_closed_or_blessed_succ[i] = true;
                    changed = true;
                }
            }
        }

        self.entries
            .iter()
            .enumerate()
            .map(|(i, _)| {
                if cursed[i] && !has_outgoing[i] {
                    NtKind::Dead
                } else if cursed[i] {
                    NtKind::Undead
                } else if blessed[i] {
                    NtKind::Blessed
                } else {
                    NtKind::Free
                }
            })
            .collect()
    }

    /// Return the 1-based ids of entries that fail to reconstruct or that
    /// cannot accept the central tile at the stored anchor. Empty result
    /// means every entry is well-formed.
    pub fn validate(&self) -> Vec<usize> {
        let match_index = Arc::new(MatchTypeIndex::new(Arc::clone(&self.tileset)));
        self.entries
            .iter()
            .enumerate()
            .filter_map(|(idx, nt)| {
                if nt_is_valid(nt, &match_index) {
                    None
                } else {
                    Some(idx + 1)
                }
            })
            .collect()
    }

    /// Serialize the entries + transitions to a text format.
    ///
    /// Format (one record per line, whitespace-separated tokens):
    ///
    /// ```text
    /// NTYPE <id> <kind> <central_id> <ac> <ax> <n_vts> <vt0> <vt1> ...
    /// TRANS <src_id> <dst_id> <tile_id> <tile_offset>
    /// ```
    ///
    /// Each `<vti>` encodes a single VertexType as
    /// `<cw_id>.<cw_off>|<inner_list>|<ccw_id>.<ccw_off>` where
    /// `<inner_list>` is either `-` (empty) or a comma-separated list of
    /// `<tile_id>.<tile_offset>` edges. `<kind>` is one of `dead`,
    /// `undead`, `blessed`, `free`. The tileset is not serialized: callers
    /// must pass a matching tileset to `parse_file`.
    pub fn write_collection(&self, out: &mut impl std::io::Write) -> std::io::Result<()> {
        let kinds = self.classify_all();
        for (idx, nt) in self.entries.iter().enumerate() {
            let id = idx + 1;
            let kind = match kinds[idx] {
                NtKind::Dead => "dead",
                NtKind::Undead => "undead",
                NtKind::Blessed => "blessed",
                NtKind::Free => "free",
            };
            write!(
                out,
                "NTYPE {} {} {} {} {} {}",
                id,
                kind,
                nt.central_tile_id,
                nt.cw_anchor_on_central,
                nt.cw_anchor_on_context,
                nt.vt_seq.len(),
            )?;
            for vt in &nt.vt_seq {
                let inner = if vt.inner.is_empty() {
                    "-".to_string()
                } else {
                    vt.inner
                        .iter()
                        .map(|e| format!("{}.{}", e.tile_id, e.tile_offset))
                        .collect::<Vec<_>>()
                        .join(",")
                };
                write!(
                    out,
                    " {}.{}|{}|{}.{}",
                    vt.cw.tile_id, vt.cw.tile_offset, inner, vt.ccw.tile_id, vt.ccw.tile_offset,
                )?;
            }
            writeln!(out)?;
        }
        for t in &self.transitions {
            writeln!(
                out,
                "TRANS {} {} {} {}",
                t.src_id, t.dst_id, t.tile_id, t.tile_offset,
            )?;
        }
        Ok(())
    }

    /// Parse a collection previously written by [`write_collection`]. The
    /// tileset must match the one used to produce the file (we do not store
    /// or verify it here).
    pub fn parse_file(tileset: Arc<TileSet<T>>, input: &str) -> Result<Self, String> {
        let mut entries: Vec<NeighborhoodType> = Vec::new();
        let mut transitions: Vec<NtTransition> = Vec::new();
        for (lineno, line) in input.lines().enumerate() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let mut tok = line.split_whitespace();
            let kw = tok
                .next()
                .ok_or_else(|| format!("line {}: empty record", lineno + 1))?;
            match kw {
                "NTYPE" => {
                    let nt = parse_ntype_line(&mut tok, lineno + 1)?;
                    entries.push(nt);
                }
                "TRANS" => {
                    let parse_usize = |s: Option<&str>, what: &str| -> Result<usize, String> {
                        s.ok_or_else(|| format!("line {}: missing {}", lineno + 1, what))?
                            .parse::<usize>()
                            .map_err(|e| format!("line {}: bad {}: {}", lineno + 1, what, e))
                    };
                    let src_id = parse_usize(tok.next(), "src_id")?;
                    let dst_id = parse_usize(tok.next(), "dst_id")?;
                    let tile_id = parse_usize(tok.next(), "tile_id")?;
                    let tile_offset = parse_usize(tok.next(), "tile_offset")?;
                    transitions.push(NtTransition {
                        src_id,
                        dst_id,
                        tile_id,
                        tile_offset,
                    });
                }
                other => {
                    return Err(format!(
                        "line {}: unknown record kind `{}`",
                        lineno + 1,
                        other
                    ));
                }
            }
        }
        Ok(NeighborhoodIndex {
            tileset,
            entries,
            transitions,
        })
    }
}

fn parse_edge_info(s: &str) -> Result<EdgeInfo, String> {
    let (id, off) = s
        .split_once('.')
        .ok_or_else(|| format!("bad edge `{}`: expected `id.offset`", s))?;
    let tile_id: usize = id
        .parse()
        .map_err(|e| format!("bad edge tile_id `{}`: {}", id, e))?;
    let tile_offset: usize = off
        .parse()
        .map_err(|e| format!("bad edge tile_offset `{}`: {}", off, e))?;
    Ok(EdgeInfo {
        tile_id,
        tile_offset,
    })
}

fn parse_vt_token(tok: &str) -> Result<VertexType, String> {
    let parts: Vec<&str> = tok.split('|').collect();
    if parts.len() != 3 {
        return Err(format!(
            "bad VT `{}`: expected `cw|inner|ccw`, got {} parts",
            tok,
            parts.len()
        ));
    }
    let cw = parse_edge_info(parts[0])?;
    let inner = if parts[1] == "-" {
        Vec::new()
    } else {
        parts[1]
            .split(',')
            .map(parse_edge_info)
            .collect::<Result<Vec<_>, _>>()?
    };
    let ccw = parse_edge_info(parts[2])?;
    Ok(VertexType { cw, inner, ccw })
}

fn parse_ntype_line(
    tok: &mut std::str::SplitWhitespace<'_>,
    lineno: usize,
) -> Result<NeighborhoodType, String> {
    let mut next = |what: &str| -> Result<String, String> {
        tok.next()
            .map(|s| s.to_string())
            .ok_or_else(|| format!("line {}: missing {}", lineno, what))
    };
    let _id = next("id")?; // id is implicit in order; we ignore it
    let _kind = next("kind")?; // kind is informational; recomputed via classify_all
    let central_tile_id: usize = next("central_id")?
        .parse()
        .map_err(|e| format!("line {}: bad central_id: {}", lineno, e))?;
    let cw_anchor_on_central: usize = next("ac")?
        .parse()
        .map_err(|e| format!("line {}: bad ac: {}", lineno, e))?;
    let cw_anchor_on_context: usize = next("ax")?
        .parse()
        .map_err(|e| format!("line {}: bad ax: {}", lineno, e))?;
    let n_vts: usize = next("n_vts")?
        .parse()
        .map_err(|e| format!("line {}: bad n_vts: {}", lineno, e))?;
    let mut vt_seq = Vec::with_capacity(n_vts);
    for i in 0..n_vts {
        let s = next(&format!("vt[{}]", i))?;
        vt_seq.push(parse_vt_token(&s).map_err(|e| format!("line {}: {}", lineno, e))?);
    }
    Ok(NeighborhoodType {
        central_tile_id,
        cw_anchor_on_central,
        cw_anchor_on_context,
        vt_seq,
    })
}

#[allow(dead_code)]
struct AugmentedContext<T: IsComplex> {
    augmented: GrowingPatch<T>,
    gap_start: usize,
    gap_len: usize,
    central_ptid: usize,
}

#[allow(dead_code)]
struct FrontierInfo {
    dist_to_frontier: usize,
}

#[allow(dead_code)]
struct ExploreOutcome<T: IsComplex> {
    petal_pm: PatchMatch,
    kind: OutcomeKind<T>,
}

#[allow(dead_code, clippy::large_enum_variant)]
enum OutcomeKind<T: IsComplex> {
    Closed,
    Open {
        nt: NeighborhoodType,
        trial: GrowingPatch<T>,
        central_ptid: usize,
    },
}

#[allow(dead_code)]
fn attach_central<T: IsComplex + IsRingOrField + Units>(
    nt: &NeighborhoodType,
    match_index: &Arc<MatchTypeIndex<T>>,
) -> Option<AugmentedContext<T>> {
    let (mut context, junc_positions) =
        GrowingPatch::construct_witness_from_vt_sequence(&nt.vt_seq, Arc::clone(match_index))?;
    let first_junc = junc_positions[0];
    let ctx_n = context.boundary_len();
    let anchor_pos = (first_junc + ctx_n - nt.cw_anchor_on_context) % ctx_n;

    let pm = PatchMatch {
        start_a: anchor_pos,
        len: {
            let tileset = match_index.tileset();
            let central_seq = tileset.rat(nt.central_tile_id).seq();
            crate::intgeom::patch::forward_match_length(
                context.angles(),
                anchor_pos,
                central_seq,
                nt.cw_anchor_on_central,
            )
        },
        start_b: nt.cw_anchor_on_central,
        tile_id: nt.central_tile_id,
    };

    let tileset = match_index.tileset();
    let central_len = tileset.rat(nt.central_tile_id).seq().len();
    let ctx_n = context.boundary_len();

    context.add_tile(&pm)?;

    let aug_n = context.boundary_len();
    let gap_len = central_len - pm.len;
    let seg_len_old = ctx_n - pm.len;
    let gap_start = seg_len_old % aug_n;
    let central_ptid = context.patch_tile_ids()[gap_start];

    Some(AugmentedContext {
        augmented: context,
        gap_start,
        gap_len,
        central_ptid,
    })
}

#[allow(dead_code)]
fn find_gap_frontier<T: IsComplex + IsRingOrField + Units>(
    aug: &AugmentedContext<T>,
) -> Option<FrontierInfo> {
    let patch = &aug.augmented;
    let n = patch.boundary_len();

    let gap_end = (aug.gap_start + aug.gap_len) % n;

    let mut frontier_pos = gap_end;
    for _ in 0..n {
        if patch.is_junction(frontier_pos) {
            break;
        }
        frontier_pos = (frontier_pos + 1) % n;
    }
    if !patch.is_junction(frontier_pos) {
        return None;
    }

    let mut dist_to_frontier = aug.gap_start;
    for i in (0..aug.gap_start).rev() {
        if patch.is_junction(i) {
            dist_to_frontier = aug.gap_start - i - 1;
            break;
        }
    }

    Some(FrontierInfo { dist_to_frontier })
}

#[allow(dead_code)]
fn explore_step<T: IsComplex + IsRingOrField + Units>(
    nt: &NeighborhoodType,
    match_index: &Arc<MatchTypeIndex<T>>,
) -> Option<(AugmentedContext<T>, FrontierInfo, Vec<ExploreOutcome<T>>)> {
    let aug = attach_central(nt, match_index)?;
    let frontier = find_gap_frontier(&aug)?;
    let patch = &aug.augmented;
    let n = patch.boundary_len();
    let central_ptid = aug.central_ptid;

    let gap_end = (aug.gap_start + aug.gap_len) % n;
    let mut frontier_pos = gap_end;
    for _ in 0..n {
        if patch.is_junction(frontier_pos) {
            break;
        }
        frontier_pos = (frontier_pos + 1) % n;
    }

    let candidates = GrowingPatch::compute_candidates_covering_position(
        patch.match_index(),
        patch.angles(),
        patch.edges(),
        frontier_pos,
    );

    let mut outcomes = Vec::new();
    for petal_pm in &candidates {
        let mut trial = aug.augmented.clone_for_mutation();
        if trial.add_tile(petal_pm).is_none() {
            continue;
        }

        let kind = if find_remaining_gap(&trial, central_ptid).is_none() {
            OutcomeKind::Closed
        } else {
            OutcomeKind::Open {
                nt: nt.clone(),
                trial,
                central_ptid,
            }
        };
        outcomes.push(ExploreOutcome {
            petal_pm: petal_pm.clone(),
            kind,
        });
    }

    Some((aug, frontier, outcomes))
}

#[allow(dead_code)]
fn find_remaining_gap<T: IsComplex + IsRingOrField + Units>(
    trial: &GrowingPatch<T>,
    central_ptid: usize,
) -> Option<(usize, usize)> {
    let n = trial.boundary_len();
    let ptids = trial.patch_tile_ids();
    let first = (0..n).find(|&i| ptids[i] == central_ptid)?;
    let mut len = 1;
    while len < n && ptids[(first + len) % n] == central_ptid {
        len += 1;
    }
    Some((first, len))
}

#[allow(dead_code)]
struct AttachedContext<T: IsComplex> {
    aug: GrowingPatch<T>,
    central_ptid: usize,
    frontier_pos_on_aug: usize,
    frontier_is_junction_in_ctx: bool,
    last_covered_ctx_edge: EdgeInfo,
    /// Number of context boundary edges on `aug` (= `ctx_n - match_len`).
    /// On `aug`, positions `[0, gap_start)` are surviving context edges
    /// and `[gap_start, aug_n)` are the central tile gap edges.
    gap_start: usize,
    aug_n: usize,
}

#[allow(dead_code)]
fn build_attached_context<T: IsComplex + IsRingOrField + Units>(
    nt: &NeighborhoodType,
    match_index: &Arc<MatchTypeIndex<T>>,
) -> Option<AttachedContext<T>> {
    let (ctx, junc_positions) =
        GrowingPatch::construct_witness_from_vt_sequence(&nt.vt_seq, Arc::clone(match_index))?;
    let ctx_n = ctx.boundary_len();
    if ctx_n == 0 {
        return None;
    }
    let first_junc = junc_positions[0];
    let anchor_pos = (first_junc + ctx_n - nt.cw_anchor_on_context) % ctx_n;

    let tileset = match_index.tileset();
    let central_seq = tileset.rat(nt.central_tile_id).seq();
    let central_n = central_seq.len();
    let match_len = crate::intgeom::patch::forward_match_length(
        ctx.angles(),
        anchor_pos,
        central_seq,
        nt.cw_anchor_on_central,
    );
    if match_len == 0 || match_len >= central_n {
        return None;
    }

    let frontier_on_ctx = (anchor_pos + match_len) % ctx_n;
    let frontier_is_junction_in_ctx = ctx.is_junction(frontier_on_ctx);
    let last_covered_ctx_edge = ctx.edges()[(anchor_pos + match_len - 1) % ctx_n];

    let central_pm = PatchMatch {
        start_a: anchor_pos,
        len: match_len,
        start_b: nt.cw_anchor_on_central,
        tile_id: nt.central_tile_id,
    };
    let mut aug = ctx.clone_for_mutation();
    aug.add_tile(&central_pm)?;

    let aug_n = aug.boundary_len();
    let gap_start = ctx_n - match_len;
    let central_ptid = aug.patch_tile_ids()[gap_start];

    // gap_end == 0 by construction; advance CCW until a junction is found.
    let mut frontier_pos_on_aug = 0usize;
    let mut found = false;
    for _ in 0..aug_n {
        if aug.is_junction(frontier_pos_on_aug) {
            found = true;
            break;
        }
        frontier_pos_on_aug = (frontier_pos_on_aug + 1) % aug_n;
    }
    if !found {
        return None;
    }

    Some(AttachedContext {
        aug,
        central_ptid,
        frontier_pos_on_aug,
        frontier_is_junction_in_ctx,
        last_covered_ctx_edge,
        gap_start,
        aug_n,
    })
}

#[allow(dead_code)]
fn explore_one<T: IsComplex + IsRingOrField + Units>(
    nt: &NeighborhoodType,
    match_index: &Arc<MatchTypeIndex<T>>,
) -> Vec<ExploreOutcome<T>> {
    let ac = match build_attached_context(nt, match_index) {
        Some(a) => a,
        None => return Vec::new(),
    };

    let candidates = GrowingPatch::compute_candidates_covering_position(
        ac.aug.match_index(),
        ac.aug.angles(),
        ac.aug.edges(),
        ac.frontier_pos_on_aug,
    );

    let mut results = Vec::new();
    for petal_pm in candidates {
        let mut trial = ac.aug.clone_for_mutation();
        if trial.add_tile(&petal_pm).is_none() {
            continue;
        }

        // Closed: petal absorbed all central edges.
        if !trial.patch_tile_ids().contains(&ac.central_ptid) {
            results.push(ExploreOutcome {
                petal_pm,
                kind: OutcomeKind::Closed,
            });
            continue;
        }

        // Count how many of the petal's matched edges (on aug) are in the
        // gap region [gap_start, aug_n). Those edges are NOT matched in the
        // ctx-only view (the central tile is absent there) and become
        // surviving boundary edges of the petal in ctx-only.
        //
        // In aug, the petal's matched edges are at petal tile offsets
        // (start_b - 1), (start_b - 2), ..., (start_b - mlen) (mod petal_n);
        // the first surviving is at offset start_b. In ctx-only the matched
        // offsets are only the ones whose corresponding aug positions are
        // context positions (i.e. in [0, gap_start)). The remaining matched
        // offsets (gap-consumed) become the first surviving offsets going
        // CCW from the old frontier in ctx-only, so the ctx-only first
        // surviving offset is `start_b - petal_gap_consumed`.
        let petal_n = match_index.tileset().rat(petal_pm.tile_id).seq().len();
        let mut petal_gap_consumed = 0usize;
        for i in 0..petal_pm.len {
            let p = (petal_pm.start_a + i) % ac.aug_n;
            if p >= ac.gap_start {
                petal_gap_consumed += 1;
            }
        }
        let petal_edge = EdgeInfo {
            tile_id: petal_pm.tile_id,
            tile_offset: (petal_pm.start_b + petal_n - petal_gap_consumed) % petal_n,
        };

        let mut new_vt_seq = nt.vt_seq.clone();
        if ac.frontier_is_junction_in_ctx {
            // Frontier vertex was already a ctx-ctx junction (last VT in vt_seq).
            // The petal becomes the new ccw-most tile at this junction; the old
            // ccw tile is pushed into `inner`.
            let last = new_vt_seq.last_mut().expect("vt_seq is non-empty");
            let old_ccw = last.ccw;
            last.inner.push(old_ccw);
            last.ccw = petal_edge;
        } else {
            // Frontier vertex was not a ctx-ctx junction; the petal creates a
            // new junction there. cw side is the last covered ctx edge.
            new_vt_seq.push(VertexType {
                cw: ac.last_covered_ctx_edge,
                inner: vec![],
                ccw: petal_edge,
            });
        }

        let new_nt = NeighborhoodType {
            central_tile_id: nt.central_tile_id,
            cw_anchor_on_central: nt.cw_anchor_on_central,
            cw_anchor_on_context: nt.cw_anchor_on_context,
            vt_seq: new_vt_seq,
        };

        // Verify the abstract vt_seq update produced a valid NT: it must
        // reconstruct + glue the central. The petal glue itself was already
        // validated on the augmented patch (`trial.add_tile` succeeded), but
        // the abstract Case A / Case B update may not always capture the
        // geometry (e.g. for non-convex tiles like spectre).
        if !nt_is_valid(&new_nt, match_index) {
            continue;
        }

        results.push(ExploreOutcome {
            petal_pm,
            kind: OutcomeKind::Open {
                nt: new_nt,
                trial,
                central_ptid: ac.central_ptid,
            },
        });
    }
    results
}

fn nt_is_valid<T: IsComplex + IsRingOrField + Units>(
    nt: &NeighborhoodType,
    match_index: &Arc<MatchTypeIndex<T>>,
) -> bool {
    let Some((mut ctx, jp)) =
        GrowingPatch::construct_witness_from_vt_sequence(&nt.vt_seq, Arc::clone(match_index))
    else {
        return false;
    };
    let ctx_n = ctx.boundary_len();
    if ctx_n == 0 || jp.is_empty() {
        return false;
    }
    let first_junc = jp[0];
    let anchor = (first_junc + ctx_n - nt.cw_anchor_on_context) % ctx_n;
    let tileset = match_index.tileset();
    let central_seq = tileset.rat(nt.central_tile_id).seq();
    let central_n = central_seq.len();
    let mlen = crate::intgeom::patch::forward_match_length(
        ctx.angles(),
        anchor,
        central_seq,
        nt.cw_anchor_on_central,
    );
    if mlen == 0 || mlen >= central_n {
        return false;
    }
    let pm = PatchMatch {
        start_a: anchor,
        len: mlen,
        start_b: nt.cw_anchor_on_central,
        tile_id: nt.central_tile_id,
    };
    ctx.add_tile(&pm).is_some()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic::ZZ12;
    use crate::intgeom::matchtypes::MatchTypeIndex;
    use crate::intgeom::patch::GrowingPatch;
    use crate::intgeom::rat::Rat;
    use crate::intgeom::tiles;

    fn square_tileset() -> Arc<TileSet<ZZ12>> {
        let sq = tiles::square::<ZZ12>();
        let rat = Rat::try_from(&sq).unwrap();
        Arc::new(TileSet::new(vec![rat]))
    }

    fn hex_tileset() -> Arc<TileSet<ZZ12>> {
        let hex = tiles::hexagon::<ZZ12>();
        let rat = Rat::try_from(&hex).unwrap();
        Arc::new(TileSet::new(vec![rat]))
    }

    #[test]
    fn square_seed_count() {
        let idx = NeighborhoodIndex::new(square_tileset());
        assert!(idx.num_types() > 0, "expected non-empty seed collection");
        for nt in idx.entries() {
            assert_eq!(nt.central_tile_id, 0, "single-tile tileset");
            assert!(!nt.vt_seq.is_empty(), "seeds must have non-empty vt_seq");
        }
    }

    #[test]
    fn hex_seed_count() {
        let idx = NeighborhoodIndex::new(hex_tileset());
        assert!(idx.num_types() > 0, "expected non-empty seed collection");
        for nt in idx.entries() {
            assert_eq!(nt.central_tile_id, 0, "single-tile tileset");
            assert!(!nt.vt_seq.is_empty(), "seeds must have non-empty vt_seq");
        }
    }

    #[test]
    fn seeds_have_no_duplicate_keys() {
        let idx = NeighborhoodIndex::new(hex_tileset());
        let mut keys = std::collections::HashSet::new();
        for (i, nt) in idx.entries().iter().enumerate() {
            assert!(
                keys.insert((
                    nt.central_tile_id,
                    nt.cw_anchor_on_central,
                    nt.cw_anchor_on_context,
                    nt.vt_seq.clone(),
                )),
                "[{}] duplicate seed key: central={}, anchor={}, ctx_anchor={}, vt_seq={:?}",
                i,
                nt.central_tile_id,
                nt.cw_anchor_on_central,
                nt.cw_anchor_on_context,
                nt.vt_seq,
            );
        }
    }

    fn assert_roundtrip<T: IsComplex + IsRingOrField + Units>(idx: &NeighborhoodIndex<T>) {
        let mut buf: Vec<u8> = Vec::new();
        idx.write_collection(&mut buf).expect("write");
        let text = std::str::from_utf8(&buf).expect("utf8");
        let parsed = NeighborhoodIndex::parse_file(Arc::clone(idx.tileset()), text).expect("parse");
        assert_eq!(
            parsed.num_types(),
            idx.num_types(),
            "entries count mismatch"
        );
        for (i, (a, b)) in idx
            .entries()
            .iter()
            .zip(parsed.entries().iter())
            .enumerate()
        {
            assert_eq!(a, b, "entry {} mismatch", i);
        }
        assert_eq!(
            parsed.transitions().len(),
            idx.transitions().len(),
            "transitions count mismatch"
        );
        for (i, (a, b)) in idx
            .transitions()
            .iter()
            .zip(parsed.transitions().iter())
            .enumerate()
        {
            assert_eq!(a.src_id, b.src_id, "transition {} src", i);
            assert_eq!(a.dst_id, b.dst_id, "transition {} dst", i);
            assert_eq!(a.tile_id, b.tile_id, "transition {} tile_id", i);
            assert_eq!(a.tile_offset, b.tile_offset, "transition {} tile_offset", i);
        }
    }

    #[test]
    fn square_roundtrip_collection() {
        let idx = NeighborhoodIndex::new(square_tileset());
        assert_roundtrip(&idx);
    }

    #[test]
    fn hex_roundtrip_collection() {
        let idx = NeighborhoodIndex::new(hex_tileset());
        assert_roundtrip(&idx);
    }

    fn spectre_tileset() -> Arc<TileSet<ZZ12>> {
        let sp = tiles::spectre::<ZZ12>();
        let rat = Rat::try_from(&sp).unwrap();
        Arc::new(TileSet::new(vec![rat]))
    }

    #[test]
    fn validate_returns_empty_for_square_and_hex() {
        let sq = NeighborhoodIndex::new(square_tileset());
        let bad = sq.validate();
        assert!(
            bad.is_empty(),
            "square has {} invalid entries (first few: {:?})",
            bad.len(),
            &bad[..bad.len().min(5)]
        );
        let hex = NeighborhoodIndex::new(hex_tileset());
        let bad = hex.validate();
        assert!(
            bad.is_empty(),
            "hex has {} invalid entries (first few: {:?})",
            bad.len(),
            &bad[..bad.len().min(5)]
        );
    }

    #[test]
    fn spectre_validates() {
        let idx = NeighborhoodIndex::new(spectre_tileset());
        let bad = idx.validate();
        assert!(
            bad.is_empty(),
            "spectre has {}/{} invalid entries (first few: {:?})",
            bad.len(),
            idx.num_types(),
            &bad[..bad.len().min(5)]
        );
        // Sanity: spectre tiles the plane, so the BFS must reach closed
        // configurations from at least some open NTs.
        let n_closed = idx
            .transitions()
            .iter()
            .filter(|t| t.dst_id == NT_CLOSED_ID)
            .count();
        assert!(
            n_closed > 0,
            "spectre BFS should produce closed transitions"
        );
        let kinds = idx.classify_all();
        let n_blessed = kinds.iter().filter(|k| **k == NtKind::Blessed).count();
        assert!(n_blessed > 0, "spectre should have Blessed entries");
    }

    #[test]
    fn classify_all_returns_one_per_entry() {
        for idx in [
            NeighborhoodIndex::new(square_tileset()),
            NeighborhoodIndex::new(hex_tileset()),
        ] {
            let kinds = idx.classify_all();
            assert_eq!(kinds.len(), idx.num_types());
            // BFS reaches a frontier where every petal closes, so there must
            // be at least one Blessed entry (= all outgoing transitions go to
            // closed).
            let n_blessed = kinds.iter().filter(|k| **k == NtKind::Blessed).count();
            assert!(n_blessed > 0, "expected at least one Blessed entry");
        }
    }

    #[test]
    fn bfs_produces_transitions() {
        let sq_idx = NeighborhoodIndex::new(square_tileset());
        assert!(sq_idx.num_types() > 0);
        assert!(
            !sq_idx.transitions().is_empty(),
            "BFS should produce transitions"
        );
        // src_id is a 1-based entry id; dst_id is either 1-based or
        // NT_CLOSED_ID (= 0).
        for t in sq_idx.transitions() {
            assert!(t.src_id >= 1 && t.src_id <= sq_idx.num_types());
            assert!(t.dst_id == NT_CLOSED_ID || t.dst_id <= sq_idx.num_types());
        }

        let hex_idx = NeighborhoodIndex::new(hex_tileset());
        assert!(hex_idx.num_types() > 0);
        assert!(!hex_idx.transitions().is_empty());
        for t in hex_idx.transitions() {
            assert!(t.src_id >= 1 && t.src_id <= hex_idx.num_types());
            assert!(t.dst_id == NT_CLOSED_ID || t.dst_id <= hex_idx.num_types());
        }
    }

    fn validate_seeds<T: IsComplex + IsRingOrField + Units>(
        idx: &NeighborhoodIndex<T>,
    ) -> Vec<String> {
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        let mut errors = Vec::new();
        for nt in idx.entries() {
            let Some((mut ctx, junc_positions)) =
                GrowingPatch::construct_witness_from_vt_sequence(&nt.vt_seq, Arc::clone(&mi))
            else {
                errors.push(format!(
                    "nt c={} ac={} ax={}: reconstruct failed",
                    nt.central_tile_id, nt.cw_anchor_on_central, nt.cw_anchor_on_context
                ));
                continue;
            };
            let ctx_n = ctx.boundary_len();
            let first_junc = junc_positions[0];
            let anchor_pos = (first_junc + ctx_n - nt.cw_anchor_on_context) % ctx_n;
            let found = ctx.get_all_matches().into_iter().find(|pm| {
                pm.tile_id == nt.central_tile_id
                    && pm.start_a == anchor_pos
                    && pm.start_b == nt.cw_anchor_on_central
            });
            let Some(pm) = found else {
                errors.push(format!(
                    "nt c={} ac={} ax={}: no match at anchor positions",
                    nt.central_tile_id, nt.cw_anchor_on_central, nt.cw_anchor_on_context
                ));
                continue;
            };
            if ctx.add_tile(&pm).is_none() {
                errors.push(format!(
                    "nt c={} ac={} ax={}: add_tile failed",
                    nt.central_tile_id, nt.cw_anchor_on_central, nt.cw_anchor_on_context
                ));
            }
        }
        errors
    }

    #[test]
    fn square_seeds_validate() {
        let idx = NeighborhoodIndex::new(square_tileset());
        let errors = validate_seeds(&idx);
        assert!(
            errors.is_empty(),
            "validation errors:\n{}",
            errors.join("\n")
        );
    }

    #[test]
    fn hex_seeds_validate() {
        let idx = NeighborhoodIndex::new(hex_tileset());
        let errors = validate_seeds(&idx);
        assert!(
            errors.is_empty(),
            "validation errors:\n{}",
            errors.join("\n")
        );
    }

    #[test]
    fn attach_central_square_seeds() {
        let idx = NeighborhoodIndex::new(square_tileset());
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        for (i, nt) in idx.entries().iter().enumerate() {
            let aug = attach_central(nt, &mi)
                .unwrap_or_else(|| panic!("seed {}: attach_central failed for nt {:?}", i, nt));
            assert!(aug.gap_len > 0, "seed {}: gap_len should be > 0", i);
            assert!(
                aug.gap_start < aug.augmented.boundary_len(),
                "seed {}: gap_start out of range",
                i
            );
            let tileset = mi.tileset();
            let central_len = tileset.rat(nt.central_tile_id).seq().len();
            assert!(
                aug.gap_len < central_len,
                "seed {}: gap_len {} should be < central_len {}",
                i,
                aug.gap_len,
                central_len
            );
        }
    }

    #[test]
    fn attach_central_hex_seeds() {
        let idx = NeighborhoodIndex::new(hex_tileset());
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        for (i, nt) in idx.entries().iter().enumerate() {
            let aug = attach_central(nt, &mi)
                .unwrap_or_else(|| panic!("seed {}: attach_central failed for nt {:?}", i, nt));
            assert!(aug.gap_len > 0, "seed {}: gap_len should be > 0", i);
        }
    }

    #[test]
    fn find_gap_frontier_square_seeds() {
        let idx = NeighborhoodIndex::new(square_tileset());
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        for (i, nt) in idx.entries().iter().enumerate() {
            let (aug, frontier, petal_outcomes) =
                explore_step(nt, &mi).unwrap_or_else(|| panic!("seed {}: explore_step failed", i));
            assert!(
                frontier.dist_to_frontier < aug.augmented.boundary_len(),
                "seed {}: dist_to_frontier out of range",
                i
            );
            assert!(
                !petal_outcomes.is_empty(),
                "seed {}: should have at least one successful petal",
                i
            );
        }
    }

    #[test]
    fn find_gap_frontier_hex_seeds() {
        let idx = NeighborhoodIndex::new(hex_tileset());
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        for (i, nt) in idx.entries().iter().enumerate() {
            let (aug, frontier, petal_outcomes) =
                explore_step(nt, &mi).unwrap_or_else(|| panic!("seed {}: explore_step failed", i));
            assert!(
                frontier.dist_to_frontier < aug.augmented.boundary_len(),
                "seed {}: dist_to_frontier out of range",
                i
            );
            assert!(
                !petal_outcomes.is_empty(),
                "seed {}: should have at least one successful petal",
                i
            );
        }
    }

    #[test]
    fn gap_edges_are_central_tile() {
        let idx = NeighborhoodIndex::new(square_tileset());
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        for (i, nt) in idx.entries().iter().enumerate() {
            let aug = attach_central(nt, &mi)
                .unwrap_or_else(|| panic!("seed {}: attach_central failed", i));
            let edges = aug.augmented.edges();
            let n = aug.augmented.boundary_len();
            let ptids = aug.augmented.patch_tile_ids();
            let central_ptid = ptids[aug.gap_start];
            for k in 0..aug.gap_len {
                let pos = (aug.gap_start + k) % n;
                assert_eq!(
                    edges[pos].tile_id, nt.central_tile_id,
                    "seed {}: gap edge at pos {} should be central tile",
                    i, pos
                );
                assert_eq!(
                    ptids[pos], central_ptid,
                    "seed {}: gap edge at pos {} should have same ptid",
                    i, pos
                );
            }
        }
    }

    #[test]
    fn frontier_adjacent_to_gap() {
        let idx = NeighborhoodIndex::new(hex_tileset());
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        for (i, nt) in idx.entries().iter().enumerate() {
            let (aug, _frontier, _petal_outcomes) =
                explore_step(nt, &mi).unwrap_or_else(|| panic!("seed {}: explore_step failed", i));
            let n = aug.augmented.boundary_len();
            let gap_end = (aug.gap_start + aug.gap_len) % n;
            assert!(
                aug.augmented.is_junction(gap_end),
                "seed {}: gap_end {} should be a junction",
                i,
                gap_end
            );
        }
    }

    #[test]
    fn explore_step_square_has_petals() {
        let idx = NeighborhoodIndex::new(square_tileset());
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        let nt = &idx.entries()[0];
        let (_aug, _frontier, petal_outcomes) =
            explore_step(nt, &mi).expect("explore_step on first seed");
        assert!(
            !petal_outcomes.is_empty(),
            "should have at least one petal outcome"
        );
    }

    #[test]
    fn find_remaining_gap_after_petal() {
        let idx = NeighborhoodIndex::new(square_tileset());
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        for (i, nt) in idx.entries().iter().enumerate() {
            let (aug, _frontier, petal_outcomes) =
                explore_step(nt, &mi).unwrap_or_else(|| panic!("seed {}: explore_step failed", i));
            let central_ptid = aug.central_ptid;
            for outcome in &petal_outcomes {
                match &outcome.kind {
                    OutcomeKind::Closed => {}
                    OutcomeKind::Open {
                        trial,
                        central_ptid: cptid,
                        ..
                    } => {
                        let gap = find_remaining_gap(trial, *cptid);
                        let n = trial.boundary_len();
                        let ptids = trial.patch_tile_ids();
                        let central_count = ptids.iter().filter(|&&id| id == central_ptid).count();
                        if central_count == 0 {
                            assert!(
                                gap.is_none(),
                                "seed {}: gap should be None when fully consumed",
                                i
                            );
                        } else {
                            let (gs, gl) = gap.unwrap_or_else(|| {
                                panic!(
                                    "seed {}: gap should exist ({} central edges)",
                                    i, central_count
                                )
                            });
                            assert_eq!(gl, central_count, "seed {}: gap length mismatch", i);
                            for k in 0..gl {
                                assert_eq!(
                                    ptids[(gs + k) % n],
                                    central_ptid,
                                    "seed {}: gap edge at offset {} should be central_ptid",
                                    i,
                                    k
                                );
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn explore_one_hex_seeds_reconstruct() {
        let idx = NeighborhoodIndex::new(hex_tileset());
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        let mut checked = 0usize;
        let mut errors = Vec::new();
        for (i, nt) in idx.entries().iter().enumerate() {
            for outcome in explore_one(nt, &mi) {
                let new_nt = match &outcome.kind {
                    OutcomeKind::Closed => continue,
                    OutcomeKind::Open { nt, .. } => nt.clone(),
                };
                let (mut ctx, junc_pos) = GrowingPatch::construct_witness_from_vt_sequence(
                    &new_nt.vt_seq,
                    Arc::clone(&mi),
                )
                .unwrap_or_else(|| {
                    panic!(
                        "seed {} result: reconstruct failed for new_nt {:?}",
                        i, new_nt
                    )
                });
                let first_junc = junc_pos[0];
                let ctx_n = ctx.boundary_len();
                let anchor_pos = (first_junc + ctx_n - new_nt.cw_anchor_on_context) % ctx_n;
                let found = ctx.get_all_matches().into_iter().find(|pm| {
                    pm.tile_id == new_nt.central_tile_id
                        && pm.start_a == anchor_pos
                        && pm.start_b == new_nt.cw_anchor_on_central
                });
                let Some(pm) = found else {
                    errors.push(format!(
                        "seed {}: no match c={} ac={} ax={} (ctx_n={})",
                        i,
                        new_nt.central_tile_id,
                        new_nt.cw_anchor_on_central,
                        new_nt.cw_anchor_on_context,
                        ctx.boundary_len()
                    ));
                    continue;
                };
                if ctx.add_tile(&pm).is_none() {
                    errors.push(format!(
                        "seed {}: add_tile failed for new_nt {:?}",
                        i, new_nt
                    ));
                    continue;
                }
                checked += 1;
            }
        }
        eprintln!("validated {} new NTs from hex seeds", checked);
        assert!(checked > 0);
        if !errors.is_empty() {
            for e in &errors[..errors.len().min(10)] {
                eprintln!("  {}", e);
            }
            panic!("{} reattach errors", errors.len());
        }
    }

    #[test]
    fn dist_to_frontier_computation() {
        let idx = NeighborhoodIndex::new(square_tileset());
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        for (i, nt) in idx.entries().iter().enumerate() {
            let (aug, frontier, _petal_outcomes) =
                explore_step(nt, &mi).unwrap_or_else(|| panic!("seed {}: explore_step failed", i));
            assert!(
                frontier.dist_to_frontier < aug.augmented.boundary_len(),
                "seed {}: dist_to_frontier out of range",
                i
            );
        }
    }

    #[test]
    fn explore_one_square_seeds_reconstruct() {
        let idx = NeighborhoodIndex::new(square_tileset());
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        let mut checked = 0usize;
        let mut errors = Vec::new();
        for (i, nt) in idx.entries().iter().enumerate() {
            for outcome in explore_one(nt, &mi) {
                let new_nt = match &outcome.kind {
                    OutcomeKind::Closed => continue,
                    OutcomeKind::Open { nt, .. } => nt.clone(),
                };
                let (mut ctx, junc_pos) = GrowingPatch::construct_witness_from_vt_sequence(
                    &new_nt.vt_seq,
                    Arc::clone(&mi),
                )
                .unwrap_or_else(|| {
                    panic!(
                        "seed {} result: reconstruct failed for new_nt {:?}",
                        i, new_nt
                    )
                });
                let first_junc = junc_pos[0];
                let ctx_n = ctx.boundary_len();
                let anchor_pos = (first_junc + ctx_n - new_nt.cw_anchor_on_context) % ctx_n;
                let found = ctx.get_all_matches().into_iter().find(|pm| {
                    pm.tile_id == new_nt.central_tile_id
                        && pm.start_a == anchor_pos
                        && pm.start_b == new_nt.cw_anchor_on_central
                });
                let Some(pm) = found else {
                    errors.push(format!(
                        "seed {}: no match c={} ac={} ax={} (ctx_n={})",
                        i,
                        new_nt.central_tile_id,
                        new_nt.cw_anchor_on_central,
                        new_nt.cw_anchor_on_context,
                        ctx.boundary_len()
                    ));
                    continue;
                };
                if ctx.add_tile(&pm).is_none() {
                    errors.push(format!(
                        "seed {}: add_tile failed for new_nt {:?}",
                        i, new_nt
                    ));
                    continue;
                }
                checked += 1;
            }
        }
        eprintln!("validated {} new NTs from square seeds", checked);
        assert!(checked > 0);
        if !errors.is_empty() {
            for e in &errors[..errors.len().min(10)] {
                eprintln!("  {}", e);
            }
            panic!("{} reattach errors", errors.len());
        }
    }
}
