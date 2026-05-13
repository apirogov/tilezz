use std::collections::VecDeque;
use std::sync::Arc;

use rustc_hash::FxHashMap;

use crate::cyclotomic::{IsComplex, IsRingOrField, Units};
use crate::intgeom::matchtypes::MatchTypeIndex;
use crate::intgeom::patch::{EdgeInfo, GrowingPatch, OpenVertexType, PatchMatch, TransitionSide};
use crate::intgeom::tileset::TileSet;

/// A neighborhood type: a central tile, a context patch matched against a
/// contiguous segment of the central's perimeter, and the open gap that
/// remains on the central.
///
/// The matched segment is the CCW range `[cw_anchor_on_central,
/// ccw_anchor_on_central)` on the central tile (mod `central_n`); its
/// length is `(ccw_anchor_on_central - cw_anchor_on_central) mod
/// central_n`. The corresponding covered segment on the context boundary
/// is described by `vt_seq`: the sequence of context-context junctions
/// incident with the match, in CCW order from `vt_seq[0]` (the CCW-most
/// junction reachable CW from the CW frontier) to `vt_seq.last()` (the
/// CW-most junction reachable CCW from the CCW frontier).
///
/// Frontier positions on context are parameterized by short distances
/// from `vt_seq`'s endpoints:
/// * `cw_anchor_on_context` = CW distance from `vt_seq[0]`'s junction to
///   the CW frontier vertex on context.
/// * `ccw_anchor_on_context` = CCW distance from `vt_seq.last()`'s
///   junction to the CCW frontier vertex on context.
///
/// `num_ctx_tiles` records the number of context tiles in the matched
/// patch; it grows by exactly 1 per BFS step regardless of direction, so
/// it gives the natural DAG layering for both CCW and CW transitions.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct NeighborhoodType {
    pub central_tile_id: usize,
    pub cw_anchor_on_central: usize,
    pub ccw_anchor_on_central: usize,
    pub cw_anchor_on_context: usize,
    pub ccw_anchor_on_context: usize,
    pub num_ctx_tiles: usize,
    pub vt_seq: Vec<OpenVertexType>,
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
    pub side: TransitionSide,
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
            if !patch.add_tile(&pm) {
                continue;
            }
            patch.normalize();
            let patch_n = patch.boundary_len();

            // Collect all ctx-ctx junctions with their positions, sorted by
            // boundary position. We then filter per third_pm to only those
            // incident with the match (= within the covered segment).
            let all_juncs: Vec<(usize, OpenVertexType)> = (0..patch_n)
                .filter_map(|i| patch.junction_vertex_type_at(i).map(|vt| (i, vt)))
                .collect();

            for third_pm in patch.get_all_matches() {
                // get_all_matches returns edge-compatible glues but for
                // non-convex tiles (e.g. spectre) the central might still
                // collide with the existing two-tile patch. Try the glue on
                // a clone first; skip if it fails.
                let mut trial = patch.clone();
                if !trial.add_tile(&third_pm) {
                    continue;
                }

                let central_tile_id = third_pm.tile_id;
                let central_n = match_index.tileset().rat(central_tile_id).seq().len();
                let cw_anchor_on_central = third_pm.start_b;
                // ccw_anchor_on_central is the central tile offset at the
                // CCW anchor vertex. Per the shared-edge correspondence
                // (CCW on context = CW on central), the CCW anchor vertex
                // on context maps to tile offset `start_b - match_len` mod
                // central_n.
                let ccw_anchor_on_central =
                    (cw_anchor_on_central + central_n - third_pm.len) % central_n;

                // vt_seq is the junctions incident with the central match:
                // those at positions in [anchor, anchor + match_len], i.e.
                // CCW distance from the anchor is at most match_len.
                let filtered: Vec<&(usize, OpenVertexType)> = all_juncs
                    .iter()
                    .filter(|(pos, _)| (pos + patch_n - third_pm.start_a) % patch_n <= third_pm.len)
                    .collect();
                if filtered.is_empty() {
                    continue;
                }
                let first_junc = filtered[0].0;
                let last_junc = filtered[filtered.len() - 1].0;
                let vt_seq: Vec<OpenVertexType> =
                    filtered.iter().map(|(_, vt)| vt.clone()).collect();
                // cw_anchor_on_context = CW distance from first_junc to the
                // CW anchor on context. ccw_anchor_on_context = CCW distance
                // from last_junc to the CCW anchor on context. Both stay
                // small as the BFS only adds edges OUTSIDE [cw_anchor,
                // ccw_anchor].
                let cw_anchor_on_context = (first_junc + patch_n - third_pm.start_a) % patch_n;
                let ccw_anchor_pos = (third_pm.start_a + third_pm.len) % patch_n;
                let ccw_anchor_on_context = (ccw_anchor_pos + patch_n - last_junc) % patch_n;
                let nt = NeighborhoodType {
                    central_tile_id,
                    cw_anchor_on_central,
                    ccw_anchor_on_central,
                    cw_anchor_on_context,
                    ccw_anchor_on_context,
                    num_ctx_tiles: 2,
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
                    side: outcome.side,
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
        // Worklist propagation over the reverse adjacency. O(n + m) total.
        // IDs are 1-based; NT_CLOSED_ID (= 0) is the closed sentinel and
        // never appears as a src_id.
        let n = self.entries.len();
        let mut succ_count = vec![0usize; n]; // non-closed successors of i
        let mut preds: Vec<Vec<usize>> = vec![vec![]; n]; // reverse edges (non-closed dst -> src)
        let mut has_outgoing = vec![false; n];
        for t in &self.transitions {
            debug_assert!(t.src_id != NT_CLOSED_ID, "closed cannot be a src");
            let src_idx = t.src_id - 1;
            has_outgoing[src_idx] = true;
            if t.dst_id != NT_CLOSED_ID {
                succ_count[src_idx] += 1;
                preds[t.dst_id - 1].push(src_idx);
            }
        }

        // Cursed: every non-closed successor is cursed, AND there is at
        // least one non-closed successor or no outgoing at all.
        let mut cursed = vec![false; n];
        let mut remaining = succ_count.clone();
        let mut queue: VecDeque<usize> = VecDeque::new();
        for i in 0..n {
            if !has_outgoing[i] {
                cursed[i] = true;
                queue.push_back(i);
            }
        }
        while let Some(v) = queue.pop_front() {
            for &p in &preds[v] {
                if cursed[p] {
                    continue;
                }
                remaining[p] -= 1;
                if remaining[p] == 0 && succ_count[p] > 0 {
                    cursed[p] = true;
                    queue.push_back(p);
                }
            }
        }

        // Blessed: every outgoing transition leads to Closed or Blessed.
        // Seed with entries whose only successors are closed.
        let mut blessed = vec![false; n];
        let mut remaining = succ_count.clone();
        let mut queue: VecDeque<usize> = VecDeque::new();
        for i in 0..n {
            if has_outgoing[i] && succ_count[i] == 0 {
                blessed[i] = true;
                queue.push_back(i);
            }
        }
        while let Some(v) = queue.pop_front() {
            for &p in &preds[v] {
                if blessed[p] {
                    continue;
                }
                remaining[p] -= 1;
                if remaining[p] == 0 {
                    blessed[p] = true;
                    queue.push_back(p);
                }
            }
        }

        (0..n)
            .map(|i| {
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
    /// Each `<vti>` encodes a single OpenVertexType as
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
                "NTYPE {} {} {} {} {} {} {} {} {}",
                id,
                kind,
                nt.central_tile_id,
                nt.cw_anchor_on_central,
                nt.ccw_anchor_on_central,
                nt.cw_anchor_on_context,
                nt.ccw_anchor_on_context,
                nt.num_ctx_tiles,
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
            let side = match t.side {
                TransitionSide::Cw => "cw",
                TransitionSide::Ccw => "ccw",
            };
            writeln!(
                out,
                "TRANS {} {} {} {} {}",
                t.src_id, t.dst_id, side, t.tile_id, t.tile_offset,
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
                    let side_tok = tok
                        .next()
                        .ok_or_else(|| format!("line {}: missing side", lineno + 1))?;
                    let side = match side_tok {
                        "cw" => TransitionSide::Cw,
                        "ccw" => TransitionSide::Ccw,
                        other => {
                            return Err(format!(
                                "line {}: unknown transition side `{}`",
                                lineno + 1,
                                other
                            ));
                        }
                    };
                    let tile_id = parse_usize(tok.next(), "tile_id")?;
                    let tile_offset = parse_usize(tok.next(), "tile_offset")?;
                    transitions.push(NtTransition {
                        src_id,
                        dst_id,
                        side,
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

fn parse_vt_token(tok: &str) -> Result<OpenVertexType, String> {
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
    Ok(OpenVertexType { cw, inner, ccw })
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
    let cw_anchor_on_central: usize = next("cw_ac")?
        .parse()
        .map_err(|e| format!("line {}: bad cw_ac: {}", lineno, e))?;
    let ccw_anchor_on_central: usize = next("ccw_ac")?
        .parse()
        .map_err(|e| format!("line {}: bad ccw_ac: {}", lineno, e))?;
    let cw_anchor_on_context: usize = next("cw_ax")?
        .parse()
        .map_err(|e| format!("line {}: bad cw_ax: {}", lineno, e))?;
    let ccw_anchor_on_context: usize = next("ccw_ax")?
        .parse()
        .map_err(|e| format!("line {}: bad ccw_ax: {}", lineno, e))?;
    let num_ctx_tiles: usize = next("num_ctx_tiles")?
        .parse()
        .map_err(|e| format!("line {}: bad num_ctx_tiles: {}", lineno, e))?;
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
        ccw_anchor_on_central,
        cw_anchor_on_context,
        ccw_anchor_on_context,
        num_ctx_tiles,
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
    side: TransitionSide,
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

    if !context.add_tile(&pm) {
        return None;
    }

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
        let mut trial = aug.augmented.clone();
        if !trial.add_tile(petal_pm) {
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
            side: TransitionSide::Ccw,
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
    /// CCW frontier position on aug (gap_end vertex, always a junction).
    frontier_pos_on_aug: usize,
    /// CW frontier (= anchor) position on aug = gap_start.
    cw_frontier_pos_on_aug: usize,
    frontier_is_junction_in_ctx: bool,
    cw_anchor_is_junction_in_ctx: bool,
    last_covered_ctx_edge: EdgeInfo,
    /// The first covered ctx edge (ctx.edges()[cw_anchor_pos]).
    first_covered_ctx_edge: EdgeInfo,
    /// Number of context boundary edges on `aug` (= `ctx_n - match_len`).
    /// On `aug`, positions `[0, gap_start)` are surviving context edges
    /// and `[gap_start, aug_n)` are the central tile gap edges.
    gap_start: usize,
    aug_n: usize,
    /// Old `match_len = (ccw_anchor_on_central - cw_anchor_on_central) mod
    /// central_n`. Cached here so the BFS step can compute the new
    /// match_len = old + petal_gap_consumed without recomputing.
    match_len: usize,
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
    let cw_anchor_is_junction_in_ctx = ctx.is_junction(anchor_pos);
    let last_covered_ctx_edge = ctx.edges()[(anchor_pos + match_len - 1) % ctx_n];
    let first_covered_ctx_edge = ctx.edges()[anchor_pos];

    let central_pm = PatchMatch {
        start_a: anchor_pos,
        len: match_len,
        start_b: nt.cw_anchor_on_central,
        tile_id: nt.central_tile_id,
    };
    let mut aug = ctx.clone();
    if !aug.add_tile(&central_pm) {
        return None;
    }

    let aug_n = aug.boundary_len();
    let gap_start = ctx_n - match_len;
    let central_ptid = aug.patch_tile_ids()[gap_start];

    // CCW frontier on aug: gap_end == 0 by construction.
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
    // CW frontier on aug: gap_start vertex (always a junction since the
    // central tile begins there).
    let cw_frontier_pos_on_aug = gap_start;

    Some(AttachedContext {
        aug,
        central_ptid,
        frontier_pos_on_aug,
        cw_frontier_pos_on_aug,
        frontier_is_junction_in_ctx,
        cw_anchor_is_junction_in_ctx,
        last_covered_ctx_edge,
        first_covered_ctx_edge,
        gap_start,
        aug_n,
        match_len,
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

    let mut results = Vec::new();
    let target_edge = ac.frontier_pos_on_aug;
    let ccw_candidates = GrowingPatch::compute_candidates_covering_position(
        ac.aug.match_index(),
        ac.aug.angles(),
        ac.aug.edges(),
        ac.frontier_pos_on_aug,
    );
    for petal_pm in ccw_candidates {
        if !match_absorbs_edge(&petal_pm, target_edge, ac.aug_n) {
            continue;
        }
        if let Some(outcome) = try_step_ccw(nt, &ac, petal_pm, match_index) {
            results.push(outcome);
        }
    }
    results
}

/// Whether the cyclic range of `pm.len` edges starting at `pm.start_a` in
/// a boundary of length `n` includes the edge at position `target_edge`.
fn match_absorbs_edge(pm: &PatchMatch, target_edge: usize, n: usize) -> bool {
    if pm.len == 0 {
        return false;
    }
    let end_inclusive = (pm.start_a + pm.len - 1) % n;
    if pm.start_a <= end_inclusive {
        target_edge >= pm.start_a && target_edge <= end_inclusive
    } else {
        target_edge >= pm.start_a || target_edge <= end_inclusive
    }
}

/// Try to apply a CCW petal to `nt` via the candidate `petal_pm`. Returns
/// `None` if the glue collides or the resulting vt_seq doesn't reconstruct.
fn try_step_ccw<T: IsComplex + IsRingOrField + Units>(
    nt: &NeighborhoodType,
    ac: &AttachedContext<T>,
    petal_pm: PatchMatch,
    match_index: &Arc<MatchTypeIndex<T>>,
) -> Option<ExploreOutcome<T>> {
    let mut trial = ac.aug.clone();
    if !trial.add_tile(&petal_pm) {
        return None;
    }

    // Closed: petal absorbed all central edges.
    if !trial.patch_tile_ids().contains(&ac.central_ptid) {
        return Some(ExploreOutcome {
            side: TransitionSide::Ccw,
            petal_pm,
            kind: OutcomeKind::Closed,
        });
    }

    let petal_n = match_index.tileset().rat(petal_pm.tile_id).seq().len();

    // The OLD CCW anchor vertex on aug is at position frontier_pos_on_aug
    // = 0. In the petal-vertex coordinate system, aug vertex (start_a + i)
    // ↔ petal vertex (start_b - i). So OLD CCW anchor ↔ petal vertex
    // (start_b - i_anchor) where i_anchor = (0 - start_a) mod aug_n. The
    // new ccw edge of vt_seq[-1] starts at this petal vertex going CCW =
    // petal edge[start_b - i_anchor].
    let i_anchor = (ac.aug_n - petal_pm.start_a) % ac.aug_n;
    let petal_edge = EdgeInfo {
        tile_id: petal_pm.tile_id,
        tile_offset: (petal_pm.start_b + petal_n - i_anchor) % petal_n,
    };

    let mut new_vt_seq = nt.vt_seq.clone();
    if ac.frontier_is_junction_in_ctx {
        // CCW Case A: extend last VT — petal becomes new ccw, old ccw goes
        // to the END of inner.
        let last = new_vt_seq.last_mut().expect("vt_seq is non-empty");
        let old_ccw = last.ccw;
        last.inner.push(old_ccw);
        last.ccw = petal_edge;
    } else {
        // CCW Case B: append new VT. cw = last covered ctx edge.
        new_vt_seq.push(OpenVertexType {
            cw: ac.last_covered_ctx_edge,
            inner: vec![],
            ccw: petal_edge,
        });
    }

    let new_nt = try_construct_nt_from_cw(
        nt.central_tile_id,
        nt.cw_anchor_on_central,
        nt.cw_anchor_on_context,
        nt.num_ctx_tiles + 1,
        new_vt_seq,
        match_index,
    )?;

    Some(ExploreOutcome {
        side: TransitionSide::Ccw,
        petal_pm,
        kind: OutcomeKind::Open {
            nt: new_nt,
            trial,
            central_ptid: ac.central_ptid,
        },
    })
}

/// Re-extract the canonical (vt_seq, cw_anchor_on_context, ccw_anchor_on_context)
/// from a reconstructed ctx + known anchor positions. Walks the covered
/// segment on ctx and collects junctions at each position via
/// `junction_vertex_type_at`. The resulting vt_seq is canonical — it
/// matches what `seed_gen` would have produced for the same geometric
/// state, regardless of which path the BFS took to reach it.
fn canonicalize_vt_seq_on_ctx<T: IsComplex + IsRingOrField + Units>(
    ctx: &GrowingPatch<T>,
    cw_anchor_pos: usize,
    ccw_anchor_pos: usize,
    ctx_n: usize,
) -> Option<(Vec<OpenVertexType>, usize, usize)> {
    let match_len = (ccw_anchor_pos + ctx_n - cw_anchor_pos) % ctx_n;
    let mut juncs: Vec<(usize, OpenVertexType)> = Vec::new();
    for off in 0..=match_len {
        let pos = (cw_anchor_pos + off) % ctx_n;
        if let Some(vt) = ctx.junction_vertex_type_at(pos) {
            juncs.push((pos, vt));
        }
    }
    if juncs.is_empty() {
        return None;
    }
    let first_junc = juncs[0].0;
    let last_junc = juncs[juncs.len() - 1].0;
    let cw_on_ctx = (first_junc + ctx_n - cw_anchor_pos) % ctx_n;
    let ccw_on_ctx = (ccw_anchor_pos + ctx_n - last_junc) % ctx_n;
    let vt_seq = juncs.into_iter().map(|(_, vt)| vt).collect();
    Some((vt_seq, cw_on_ctx, ccw_on_ctx))
}

/// Reconstruct an NT from its vt_seq + CW-side anchor fields, deriving the
/// CCW-side anchors from the geometry, and verifying the central can be
/// glued. The stored `vt_seq` is re-extracted from the reconstructed ctx
/// so that BFS dedup uses a canonical encoding regardless of how the NT
/// was assembled. Returns None if reconstruction fails, the match length
/// is 0, the match closes the central, or the central glue collides.
fn try_construct_nt_from_cw<T: IsComplex + IsRingOrField + Units>(
    central_tile_id: usize,
    cw_anchor_on_central: usize,
    cw_anchor_on_context: usize,
    num_ctx_tiles: usize,
    vt_seq: Vec<OpenVertexType>,
    match_index: &Arc<MatchTypeIndex<T>>,
) -> Option<NeighborhoodType> {
    let (mut ctx, jp) =
        GrowingPatch::construct_witness_from_vt_sequence(&vt_seq, Arc::clone(match_index))?;
    let ctx_n = ctx.boundary_len();
    if ctx_n == 0 || jp.is_empty() {
        return None;
    }
    let first_junc = jp[0];
    let cw_anchor_pos = (first_junc + ctx_n - cw_anchor_on_context) % ctx_n;
    let tileset = match_index.tileset();
    let central_seq = tileset.rat(central_tile_id).seq();
    let central_n = central_seq.len();
    let match_len = crate::intgeom::patch::forward_match_length(
        ctx.angles(),
        cw_anchor_pos,
        central_seq,
        cw_anchor_on_central,
    );
    if match_len == 0 || match_len >= central_n {
        return None;
    }
    let ccw_anchor_on_central = (cw_anchor_on_central + central_n - match_len) % central_n;
    let ccw_anchor_pos = (cw_anchor_pos + match_len) % ctx_n;
    let (canon_vt_seq, canon_cw_on_ctx, canon_ccw_on_ctx) =
        canonicalize_vt_seq_on_ctx(&ctx, cw_anchor_pos, ccw_anchor_pos, ctx_n)?;
    let pm = PatchMatch {
        start_a: cw_anchor_pos,
        len: match_len,
        start_b: cw_anchor_on_central,
        tile_id: central_tile_id,
    };
    if !ctx.add_tile(&pm) {
        return None;
    }
    Some(NeighborhoodType {
        central_tile_id,
        cw_anchor_on_central,
        ccw_anchor_on_central,
        cw_anchor_on_context: canon_cw_on_ctx,
        ccw_anchor_on_context: canon_ccw_on_ctx,
        num_ctx_tiles,
        vt_seq: canon_vt_seq,
    })
}

fn nt_is_valid<T: IsComplex + IsRingOrField + Units>(
    nt: &NeighborhoodType,
    match_index: &Arc<MatchTypeIndex<T>>,
) -> bool {
    let Some(reconstructed) = try_construct_nt_from_cw(
        nt.central_tile_id,
        nt.cw_anchor_on_central,
        nt.cw_anchor_on_context,
        nt.num_ctx_tiles,
        nt.vt_seq.clone(),
        match_index,
    ) else {
        return false;
    };
    reconstructed.ccw_anchor_on_central == nt.ccw_anchor_on_central
        && reconstructed.ccw_anchor_on_context == nt.ccw_anchor_on_context
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic::ZZ12;
    use crate::intgeom::matchtypes::MatchTypeIndex;
    use crate::intgeom::patch::GrowingPatch;
    use crate::intgeom::rat::Rat;
    use crate::intgeom::tiles;
    use std::sync::OnceLock;

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

    /// Shared, lazy-initialized neighborhood indices. Tests should call
    /// these instead of `NeighborhoodIndex::new(...)` to avoid rebuilding
    /// expensive indices for every test. CCW-only direction is the
    /// established working baseline.
    fn square_idx() -> &'static NeighborhoodIndex<ZZ12> {
        static SQUARE_IDX: OnceLock<NeighborhoodIndex<ZZ12>> = OnceLock::new();
        SQUARE_IDX.get_or_init(|| NeighborhoodIndex::new(square_tileset()))
    }

    fn hex_idx() -> &'static NeighborhoodIndex<ZZ12> {
        static HEX_IDX: OnceLock<NeighborhoodIndex<ZZ12>> = OnceLock::new();
        HEX_IDX.get_or_init(|| NeighborhoodIndex::new(hex_tileset()))
    }

    fn spectre_idx() -> &'static NeighborhoodIndex<ZZ12> {
        static SPECTRE_IDX: OnceLock<NeighborhoodIndex<ZZ12>> = OnceLock::new();
        SPECTRE_IDX.get_or_init(|| NeighborhoodIndex::new(spectre_tileset()))
    }

    #[test]
    fn square_seed_count() {
        let idx = square_idx();
        assert!(idx.num_types() > 0, "expected non-empty seed collection");
        for nt in idx.entries() {
            assert_eq!(nt.central_tile_id, 0, "single-tile tileset");
            assert!(!nt.vt_seq.is_empty(), "seeds must have non-empty vt_seq");
        }
    }

    #[test]
    fn hex_seed_count() {
        let idx = hex_idx();
        assert!(idx.num_types() > 0, "expected non-empty seed collection");
        for nt in idx.entries() {
            assert_eq!(nt.central_tile_id, 0, "single-tile tileset");
            assert!(!nt.vt_seq.is_empty(), "seeds must have non-empty vt_seq");
        }
    }

    #[test]
    fn seeds_have_no_duplicate_keys() {
        let idx = hex_idx();
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
        let idx = square_idx();
        assert_roundtrip(idx);
    }

    #[test]
    fn hex_roundtrip_collection() {
        let idx = hex_idx();
        assert_roundtrip(idx);
    }

    fn spectre_tileset() -> Arc<TileSet<ZZ12>> {
        let sp = tiles::spectre::<ZZ12>();
        let rat = Rat::try_from(&sp).unwrap();
        Arc::new(TileSet::new(vec![rat]))
    }

    #[test]
    fn validate_returns_empty_for_square_and_hex() {
        let sq = square_idx();
        let bad = sq.validate();
        assert!(
            bad.is_empty(),
            "square has {} invalid entries (first few: {:?})",
            bad.len(),
            &bad[..bad.len().min(5)]
        );
        let hex = hex_idx();
        let bad = hex.validate();
        assert!(
            bad.is_empty(),
            "hex has {} invalid entries (first few: {:?})",
            bad.len(),
            &bad[..bad.len().min(5)]
        );
    }

    /// Extract a canonical fingerprint of an NT's geometric shape: the
    /// lex_min_rot-normalized aug boundary angle sequence + central tile
    /// id + cw_anchor_on_central. Two NTs with the same fingerprint are
    /// the same geometric state regardless of vt_seq encoding.
    fn nt_boundary_fingerprint<T: IsComplex + IsRingOrField + Units>(
        nt: &NeighborhoodType,
        mi: &Arc<MatchTypeIndex<T>>,
    ) -> Option<(usize, usize, Vec<i8>)> {
        let (mut ctx, jp) =
            GrowingPatch::construct_witness_from_vt_sequence(&nt.vt_seq, Arc::clone(mi))?;
        let ctx_n = ctx.boundary_len();
        let first_junc = *jp.first()?;
        let cw_anchor_pos = (first_junc + ctx_n - nt.cw_anchor_on_context) % ctx_n;
        let pm = PatchMatch {
            start_a: cw_anchor_pos,
            len: (nt.cw_anchor_on_central + mi.tileset().rat(nt.central_tile_id).seq().len()
                - nt.ccw_anchor_on_central)
                % mi.tileset().rat(nt.central_tile_id).seq().len(),
            start_b: nt.cw_anchor_on_central,
            tile_id: nt.central_tile_id,
        };
        if !ctx.add_tile(&pm) {
            return None;
        }
        let angles = ctx.angles().to_vec();
        let rot = crate::intgeom::rat::lex_min_rot(&angles);
        let mut canon = angles;
        canon.rotate_left(rot);
        Some((nt.central_tile_id, nt.cw_anchor_on_central, canon))
    }

    /// Completeness check for CCW-only BFS: for every NT in the BFS set,
    /// brute-force enumerate all valid CCW petal candidates via direct
    /// `get_match` calls on the aug boundary. For each candidate, the
    /// new aug boundary (= trial.angles() after add_tile) is normalized
    /// via lex_min_rot and compared against the fingerprint set of all
    /// BFS NTs. Each candidate must yield either:
    ///   - a closed configuration (= petal absorbs all of central), OR
    ///   - an NT whose canonical boundary matches one in our set.
    ///
    /// This uses NO logic from the BFS step (try_step_*) — it's a pure
    /// brute-force independent check.
    #[test]
    fn spectre_ccw_only_is_complete() {
        let ccw = spectre_idx();
        let tileset = ccw.tileset().clone();
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(&tileset)));
        // The fingerprint of an NT's geometric shape = (central_tile_id,
        // cw_anchor_on_central, lex_min_rot of aug boundary angles).
        let bfs_fps: FxHashMap<(usize, usize, Vec<i8>), ()> = ccw
            .entries()
            .iter()
            .filter_map(|nt| nt_boundary_fingerprint(nt, &mi).map(|fp| (fp, ())))
            .collect();
        let mut missing_count = 0usize;
        let mut closed_count = 0usize;
        let mut already_in_set = 0usize;
        let mut total_petals = 0usize;
        for nt in ccw.entries() {
            let Some(ac) = build_attached_context(nt, &mi) else {
                continue;
            };
            let target_edge = ac.frontier_pos_on_aug;
            let n_aug = ac.aug_n;
            let aug_rat = crate::intgeom::rat::Rat::from_slice_unchecked(ac.aug.angles());
            let mut seen: FxHashMap<(usize, usize, usize, usize), ()> = FxHashMap::default();
            for tile_b in 0..mi.tileset().num_tiles() {
                let other = mi.tileset().rat(tile_b);
                let n_b = other.seq().len();
                for pos in 0..n_aug {
                    for b_off in 0..n_b {
                        let (ns, len, ne) = aug_rat.get_match((pos as i64, b_off as i64), other);
                        if len < 1 {
                            continue;
                        }
                        let ns_u = ns.rem_euclid(n_aug as i64) as usize;
                        let ne_u = ne.rem_euclid(n_b as i64) as usize;
                        let pm = PatchMatch {
                            start_a: ns_u,
                            len,
                            start_b: ne_u,
                            tile_id: tile_b,
                        };
                        if !match_absorbs_edge(&pm, target_edge, n_aug) {
                            continue;
                        }
                        if seen.contains_key(&(ns_u, len, ne_u, tile_b)) {
                            continue;
                        }
                        let mut trial = ac.aug.clone();
                        if !trial.add_tile(&pm) {
                            continue;
                        }
                        seen.insert((ns_u, len, ne_u, tile_b), ());
                        total_petals += 1;
                        if !trial.patch_tile_ids().contains(&ac.central_ptid) {
                            closed_count += 1;
                            continue;
                        }
                        // Brute-force compute the new fingerprint directly
                        // from the trial's boundary angles. NO try_step_*
                        // logic — purely the geometric result.
                        let trial_angles = trial.angles().to_vec();
                        let rot = crate::intgeom::rat::lex_min_rot(&trial_angles);
                        let mut canon = trial_angles;
                        canon.rotate_left(rot);
                        let new_fp = (nt.central_tile_id, nt.cw_anchor_on_central, canon);
                        if bfs_fps.contains_key(&new_fp) {
                            already_in_set += 1;
                        } else {
                            missing_count += 1;
                            if missing_count <= 3 {
                                eprintln!(
                                    "MISSING: parent num_ctx={} pm(start_a={} len={} start_b={} tile_id={})",
                                    nt.num_ctx_tiles,
                                    pm.start_a, pm.len, pm.start_b, pm.tile_id
                                );
                            }
                        }
                    }
                }
            }
        }
        eprintln!(
            "CCW completeness (brute-force, NO BFS step logic): total_petals={} closed={} already_in_set={} MISSING={}",
            total_petals, closed_count, already_in_set, missing_count
        );
        assert_eq!(
            missing_count, 0,
            "CCW BFS is incomplete: {} states missing",
            missing_count
        );
    }

    #[test]
    fn spectre_validates() {
        let idx = spectre_idx();
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
        for idx in [square_idx(), hex_idx()] {
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
        let sq_idx = square_idx();
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

        let hex_idx = hex_idx();
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
            if !ctx.add_tile(&pm) {
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
        let idx = square_idx();
        let errors = validate_seeds(idx);
        assert!(
            errors.is_empty(),
            "validation errors:\n{}",
            errors.join("\n")
        );
    }

    #[test]
    fn hex_seeds_validate() {
        let idx = hex_idx();
        let errors = validate_seeds(idx);
        assert!(
            errors.is_empty(),
            "validation errors:\n{}",
            errors.join("\n")
        );
    }

    #[test]
    fn attach_central_square_seeds() {
        let idx = square_idx();
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
        let idx = hex_idx();
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        for (i, nt) in idx.entries().iter().enumerate() {
            let aug = attach_central(nt, &mi)
                .unwrap_or_else(|| panic!("seed {}: attach_central failed for nt {:?}", i, nt));
            assert!(aug.gap_len > 0, "seed {}: gap_len should be > 0", i);
        }
    }

    #[test]
    fn find_gap_frontier_square_seeds() {
        let idx = square_idx();
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
        let idx = hex_idx();
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
        let idx = square_idx();
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
        let idx = hex_idx();
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
        let idx = square_idx();
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
        let idx = square_idx();
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
        let idx = hex_idx();
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
                if !ctx.add_tile(&pm) {
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
        let idx = square_idx();
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
        let idx = square_idx();
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
                if !ctx.add_tile(&pm) {
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
