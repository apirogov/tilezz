//! Neighborhood types: catalog of local *interfaces* between a tile
//! and any legal context patch that could meet it along one
//! contiguous edge run.
//!
//! # What we are enumerating
//!
//! For a given tileset, this module produces every distinct way a
//! tile (the **central**) can meet a surrounding patch of other tiles
//! (the **context**) along a single contiguous matched segment.
//! What matters for the interface is only the central + the shape of
//! the matched edge run; the *outer* perimeter of the context (the
//! part not touching the central) is irrelevant for "what tile can
//! grow next here". The output is therefore the set of all possible
//! such interfaces, classified by whether they extend to a fully
//! surrounded central or trap somewhere along the way.
//!
//! # Why the abstraction, not pure boundary enumeration
//!
//! A naive "enumerate all length-≤k boundary strings that can occur
//! as a match interface" approach has to separately filter realizable
//! from junk strings - a costly secondary step. Growing the context
//! one tile at a time via BFS gives realizability for free: every NT
//! is built from a known-realizable predecessor, so no candidate has
//! to be checked against geometric legality post-hoc.
//!
//! The trade-off: the BFS catalog is **finer** than the interface
//! equivalence quotient. The same interface can be reached by
//! multiple distinct growth histories (different sequences of petal
//! attachments), and each history becomes its own catalog entry. The
//! true interface quotient is a coarser, derived view.
//!
//! # Classification → forbidden patterns
//!
//! Every NT gets a `Dead`/`Undead`/`Blessed`/`Free` reachability
//! kind via [`NeighborhoodIndex::classify_all`]. The Dead and Undead
//! interfaces are the actionable output: they are the local boundary
//! patterns a tiling can never recover from. Downstream consumers
//! treat the set of Dead/Undead match-side boundary strings as
//! **forbidden patterns** that must not appear in any extending
//! glue. Together with the dead/undead VTs (a finer per-vertex
//! constraint, see `vertextypes.rs`), these constraints prune the
//! tile-growth search.
//!
//! See [`NeighborhoodType`] for the data layout and
//! [`NeighborhoodIndex::new`] for the BFS algorithm.

use std::collections::VecDeque;
use std::sync::Arc;

use rustc_hash::FxHashMap;

use crate::cyclotomic::{IsComplex, IsRingOrField, Units};
use crate::intgeom::matchtypes::MatchTypeIndex;
use crate::intgeom::patch::{EdgeInfo, GrowingPatch, OpenVertexType, PatchMatch, TransitionSide};
use crate::intgeom::tileset::TileSet;

/// One catalog entry of a local **interface** between a tile and any
/// legal context patch sharing one contiguous matched edge run with
/// it. The geometrically-significant content is the central + the
/// shape of the matched segment; everything else about the context
/// (its outer perimeter, internal tile arrangement, total tile
/// count) is incidental to the interface but is retained here to
/// support the BFS that generates the catalog.
///
/// # Interface vs BFS entry
///
/// **A `NeighborhoodType` is a BFS state, not an interface class.**
/// The catalog can contain multiple BFS entries that describe the
/// same interface, reached via different growth histories. The
/// interface-equivalence quotient (= same `central_tile_id`, same
/// `cw_anchor_on_central`, same canonical matched-segment boundary)
/// is a coarser view, computable from this catalog but not
/// represented as a separate type. Downstream consumers exporting
/// forbidden-pattern sets work at the interface level; the BFS works
/// at the entry level for combinatorial cleanliness (every step
/// adds exactly one tile, no canonicalization-collapse cases to
/// reason about).
///
/// # Field semantics
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
///
/// # Derived (redundant) fields
///
/// The CCW-side anchors are cached derivations of the CW-side anchors
/// plus the geometry; they are stored to avoid recomputation per
/// BFS step but are not independent state:
///
/// * `ccw_anchor_on_central = (cw_anchor_on_central - match_len) mod
///   central_n`, where
///   `match_len = (cw_anchor_on_central - ccw_anchor_on_central) mod
///   central_n` (the same equation in the other direction).
/// * `ccw_anchor_on_context` is determined by walking the
///   reconstructed context to find the CCW-most junction within the
///   covered segment - see [`try_construct_nt_from_cw`] for the
///   canonical re-derivation.
///
/// Constructing a `NeighborhoodType` directly (as the parser and
/// `seed_phase` do) can produce inconsistent CCW fields. The
/// invariant - "stored CCW fields equal what the geometry implies" -
/// is checked by [`nt_is_valid`] and is part of `validate`'s
/// contract.
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
    /// Build the full open-NT catalog for `tileset` by BFS over the
    /// neighborhood-transition graph.
    ///
    /// # Catalog vs interface quotient
    ///
    /// The catalog is **finer than the interface-equivalence
    /// quotient**: the same interface (= same central +
    /// same canonical matched-segment boundary) can be reached by
    /// multiple growth histories, and the BFS records each as a
    /// distinct entry. The interface quotient is therefore a derived
    /// view, not directly stored. This is intentional - the BFS gets
    /// "every step adds exactly one tile" as a clean invariant in
    /// exchange for that over-counting. Consumers wanting the
    /// quotient compute it post-hoc by grouping entries by their
    /// canonical interface fingerprint.
    ///
    /// # Algorithm
    ///
    /// 1. **Seed generation.** For every match-type `(tile_a,
    ///    start_a, tile_b, start_b, len)` in [`MatchTypeIndex`]:
    ///    build a two-tile patch by gluing `tile_b` to `tile_a` and
    ///    normalize. For every additional candidate from
    ///    [`GrowingPatch::get_all_matches`], try gluing a third
    ///    ("central") tile to a clone (non-convex tiles can collide;
    ///    skip on failure). Filter the patch's junctions to those
    ///    incident with the match (CCW distance from the match
    ///    anchor is at most the match length) - that's `vt_seq`.
    ///    Validate via [`nt_is_valid`] and dedup against the BFS's
    ///    `PartialEq` dedup key (which is the full
    ///    `NeighborhoodType`, not the interface).
    /// 2. **BFS growth.** Dequeue each NT and run [`explore_one`] to
    ///    enumerate one-petal CCW successors. New NTs are enqueued;
    ///    transitions are recorded with the right `TransitionSide`.
    ///    Every step grows `num_ctx_tiles` by exactly 1, so the
    ///    transition graph is a DAG.
    ///
    /// # CW direction
    ///
    /// CW-direction exploration is **not implemented**. The CCW-only
    /// BFS is verified complete by the brute-force regression test
    /// `spectre_ccw_only_is_complete`: starting from every NT in the
    /// catalog, every direct `get_match` candidate at the CCW
    /// frontier either closes the configuration or produces an NT
    /// already in the catalog. Earlier attempts at `try_step_cw`
    /// produced spurious NTs because the CCW Case-A/Case-B
    /// `vt_seq` mutation breaks down when the OLD CW anchor's
    /// junction status flips between old ctx and new ctx (= when
    /// the petal contributes zero angle at that vertex); a clean
    /// re-derivation is needed before reintroducing CW exploration.
    pub fn new(tileset: Arc<TileSet<T>>) -> Self {
        let match_index = Arc::new(MatchTypeIndex::new(Arc::clone(&tileset)));
        // IDs are 1-based so NT_CLOSED_ID = 0 is unambiguous: entries[i]
        // has id i + 1. `seen` and transitions store ids, not indices.
        let mut state = BfsState::default();
        seed_phase(&mut state, &match_index);
        bfs_phase(&mut state, &match_index);

        let BfsState {
            entries,
            transitions,
            ..
        } = state;

        NeighborhoodIndex {
            tileset,
            entries,
            transitions,
        }
    }

    /// All [`NeighborhoodType`] entries in BFS-discovery order.
    /// `entries()[id - 1]` is the entry for 1-based id `id`.
    pub fn entries(&self) -> &[NeighborhoodType] {
        &self.entries
    }

    /// Number of NTs in the catalog. Ids run `1..=num_types()`;
    /// [`NT_CLOSED_ID`] (= 0) is reserved as the closed-configuration
    /// sentinel.
    pub fn num_types(&self) -> usize {
        self.entries.len()
    }

    /// All recorded BFS transitions, in BFS-emission order.
    /// `src_id` is always a real 1-based entry id; `dst_id` is either
    /// a 1-based id or [`NT_CLOSED_ID`].
    pub fn transitions(&self) -> &[NtTransition] {
        &self.transitions
    }

    /// The tileset this catalog was built from.
    pub fn tileset(&self) -> &Arc<TileSet<T>> {
        &self.tileset
    }

    /// Classify every entry as Dead / Undead / Blessed / Free by
    /// walking the transition graph.
    ///
    /// * `Dead` - no outgoing transitions.
    /// * `Undead` - has transitions, every reachable path terminates
    ///   at a Dead/Undead (never at the Closed sink).
    /// * `Blessed` - every outgoing transition reaches Closed or
    ///   another Blessed.
    /// * `Free` - has at least one Closed-reaching path and at least
    ///   one Dead-terminating path.
    ///
    /// # Downstream use: forbidden patterns
    ///
    /// The Dead and Undead kinds are the actionable output of the
    /// classification. Their match-side boundary strings are
    /// **forbidden local boundary patterns** for the tileset: any
    /// tile growth that would create one of these interfaces is
    /// guaranteed to never complete to a full corona, so consumers
    /// growing tilings prune candidate glues against this set. The
    /// vertex-level analog lives in `vertextypes.rs` (Dead/Undead
    /// VTs are even smaller local forbidden patterns; every
    /// Dead/Undead NT contains at least one Dead/Undead VT, so the
    /// VT antichain is the minimal-pattern compression of the
    /// constraint set).
    ///
    /// # Algorithm
    ///
    /// O(n + m) via worklist propagation on the reverse adjacency.
    pub fn classify_all(&self) -> Vec<NtKind> {
        // Transition ids are 1-based; convert with `- 1` when indexing into
        // these per-entry vectors. NT_CLOSED_ID (= 0) is the closed sentinel
        // and never appears as a src_id.
        let n = self.entries.len();
        let mut succ_count = vec![0usize; n]; // non-closed successors of i
        let mut preds: Vec<Vec<usize>> = vec![vec![]; n]; // reverse edges (non-closed dst -> src)
        let mut has_outgoing = vec![false; n];
        let mut has_closing = vec![false; n];
        for t in &self.transitions {
            debug_assert!(t.src_id != NT_CLOSED_ID, "closed cannot be a src");
            let src_idx = t.src_id - 1;
            has_outgoing[src_idx] = true;
            if t.dst_id == NT_CLOSED_ID {
                has_closing[src_idx] = true;
            } else {
                succ_count[src_idx] += 1;
                preds[t.dst_id - 1].push(src_idx);
            }
        }

        // Cursed: every non-closed successor is cursed, AND there is at
        // least one non-closed successor or no outgoing at all, AND no
        // closing transition exists (a closing transition is an edge to
        // the Closed sink, which is not cursed - an NT with a closing
        // transition has an escape path and is therefore not cursed).
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
                if remaining[p] == 0 && succ_count[p] > 0 && !has_closing[p] {
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
                TransitionSide::Both => "both",
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
    ///
    /// Validation: NTYPE ids are 1-based and dense (`1, 2, ...,
    /// num_entries`), TRANS `src_id`/`dst_id` are in range, and the
    /// recorded NTYPE `kind` agrees with what [`Self::classify_all`]
    /// re-derives from the parsed transitions. Mismatched kinds catch
    /// stale serialized catalogs after a classification change.
    pub fn parse_file(tileset: Arc<TileSet<T>>, input: &str) -> Result<Self, String> {
        let mut entries: Vec<NeighborhoodType> = Vec::new();
        let mut transitions: Vec<NtTransition> = Vec::new();
        let mut recorded_kinds: Vec<NtKind> = Vec::new();
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
                    let (parsed_id, parsed_kind, nt) = parse_ntype_line(&mut tok, lineno + 1)?;
                    let expected_id = entries.len() + 1;
                    if parsed_id != expected_id {
                        return Err(format!(
                            "line {}: NTYPE id sequence broken (expected {}, got {})",
                            lineno + 1,
                            expected_id,
                            parsed_id
                        ));
                    }
                    entries.push(nt);
                    recorded_kinds.push(parsed_kind);
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
                        "both" => TransitionSide::Both,
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
                    if src_id < 1 || src_id > entries.len() {
                        return Err(format!(
                            "line {}: TRANS src_id {} out of range (have {} NTYPE entries)",
                            lineno + 1,
                            src_id,
                            entries.len()
                        ));
                    }
                    if dst_id != NT_CLOSED_ID && dst_id > entries.len() {
                        return Err(format!(
                            "line {}: TRANS dst_id {} out of range (have {} NTYPE entries, NT_CLOSED_ID = {})",
                            lineno + 1,
                            dst_id,
                            entries.len(),
                            NT_CLOSED_ID
                        ));
                    }
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
        let idx = NeighborhoodIndex {
            tileset,
            entries,
            transitions,
        };
        // Cross-check the recorded kinds against what classify_all
        // re-derives from the transitions. Catches stale serialized
        // catalogs (e.g. written before a classifier fix).
        let computed_kinds = idx.classify_all();
        for (i, (recorded, computed)) in recorded_kinds
            .iter()
            .zip(computed_kinds.iter())
            .enumerate()
        {
            if recorded != computed {
                return Err(format!(
                    "NTYPE id {} kind mismatch: file says {:?}, classify_all derives {:?}",
                    i + 1,
                    recorded,
                    computed
                ));
            }
        }
        Ok(idx)
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

/// Returns `(parsed_id, parsed_kind, parsed_nt)` so the caller can
/// validate the id sequence is dense + 1-based and cross-check the
/// recorded kind against `classify_all`.
fn parse_ntype_line(
    tok: &mut std::str::SplitWhitespace<'_>,
    lineno: usize,
) -> Result<(usize, NtKind, NeighborhoodType), String> {
    let mut next = |what: &str| -> Result<String, String> {
        tok.next()
            .map(|s| s.to_string())
            .ok_or_else(|| format!("line {}: missing {}", lineno, what))
    };
    let id_str = next("id")?;
    let parsed_id: usize = id_str
        .parse()
        .map_err(|e| format!("line {}: bad id: {}", lineno, e))?;
    let kind_str = next("kind")?;
    let parsed_kind = match kind_str.as_str() {
        "dead" => NtKind::Dead,
        "undead" => NtKind::Undead,
        "blessed" => NtKind::Blessed,
        "free" => NtKind::Free,
        other => {
            return Err(format!("line {}: unknown kind `{}`", lineno, other));
        }
    };
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
    Ok((
        parsed_id,
        parsed_kind,
        NeighborhoodType {
            central_tile_id,
            cw_anchor_on_central,
            ccw_anchor_on_central,
            cw_anchor_on_context,
            ccw_anchor_on_context,
            num_ctx_tiles,
            vt_seq,
        },
    ))
}

/// Mutable state shared across the seed and BFS phases of
/// [`NeighborhoodIndex::new`].
struct BfsState {
    /// All discovered NTs in BFS order. `entries[i]` has 1-based id
    /// `i + 1`.
    entries: Vec<NeighborhoodType>,
    /// All recorded transitions between NTs (and to the Closed sink).
    transitions: Vec<NtTransition>,
    /// `seen[nt]` = 1-based id of the NT in `entries`.
    seen: FxHashMap<NeighborhoodType, usize>,
    /// BFS frontier.
    queue: VecDeque<NeighborhoodType>,
}

impl Default for BfsState {
    fn default() -> Self {
        BfsState {
            entries: Vec::new(),
            transitions: Vec::new(),
            seen: FxHashMap::default(),
            queue: VecDeque::new(),
        }
    }
}

/// Phase 1: enumerate every two-tile patch from `MatchTypeIndex`, try
/// every third-tile glue against it, and accept each successful glue
/// as a seed NT (after a self-consistency check via [`nt_is_valid`]).
fn seed_phase<T: IsComplex + IsRingOrField + Units>(
    state: &mut BfsState,
    match_index: &Arc<MatchTypeIndex<T>>,
) {
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

        // Collect all ctx-ctx junctions with their positions. Filtered
        // per third_pm to only those incident with the match.
        let all_juncs: Vec<(usize, OpenVertexType)> = (0..patch_n)
            .filter_map(|i| patch.junction_vertex_type_at(i).map(|vt| (i, vt)))
            .collect();

        for third_pm in patch.get_all_matches() {
            // get_all_matches returns edge-compatible glues but for
            // non-convex tiles (e.g. spectre) the central might still
            // collide with the existing two-tile patch. Try the glue
            // on a clone first; skip if it fails.
            let mut trial = patch.clone();
            if !trial.add_tile(&third_pm) {
                continue;
            }
            if let Some(nt) =
                build_seed_nt(&third_pm, patch_n, &all_juncs, match_index)
            {
                if state.seen.contains_key(&nt) {
                    continue;
                }
                // trial.add_tile above verified the central glues into
                // the normalized two-tile patch. Also verify the
                // abstract NT reconstructs from its vt_seq (needed for
                // non-convex tiles where seed-gen and reconstruct can
                // diverge).
                if !nt_is_valid(&nt, match_index) {
                    continue;
                }
                let id = state.entries.len() + 1;
                state.seen.insert(nt.clone(), id);
                state.entries.push(nt.clone());
                state.queue.push_back(nt);
            }
        }
    }
}

/// Construct a candidate seed NT from `third_pm`'s glue position on a
/// two-tile patch. Returns `None` if no junction is incident with the
/// match (= the central tile shares no junction with the patch).
fn build_seed_nt<T: IsComplex + IsRingOrField + Units>(
    third_pm: &PatchMatch,
    patch_n: usize,
    all_juncs: &[(usize, OpenVertexType)],
    match_index: &Arc<MatchTypeIndex<T>>,
) -> Option<NeighborhoodType> {
    let central_tile_id = third_pm.tile_id;
    let central_n = match_index.tileset().rat(central_tile_id).seq().len();
    let cw_anchor_on_central = third_pm.start_b;
    // Per the shared-edge correspondence (CCW on context = CW on
    // central), the CCW anchor vertex on context maps to tile offset
    // `start_b - match_len` mod central_n.
    let ccw_anchor_on_central =
        (cw_anchor_on_central + central_n - third_pm.len) % central_n;

    // vt_seq is the junctions incident with the central match: those
    // at positions in [anchor, anchor + match_len], i.e. CCW distance
    // from the anchor is at most match_len.
    let filtered: Vec<&(usize, OpenVertexType)> = all_juncs
        .iter()
        .filter(|(pos, _)| (pos + patch_n - third_pm.start_a) % patch_n <= third_pm.len)
        .collect();
    if filtered.is_empty() {
        return None;
    }
    let first_junc = filtered[0].0;
    let last_junc = filtered[filtered.len() - 1].0;
    let vt_seq: Vec<OpenVertexType> = filtered.iter().map(|(_, vt)| vt.clone()).collect();
    // cw_anchor_on_context = CW distance from first_junc to the CW
    // anchor on context. ccw_anchor_on_context = CCW distance from
    // last_junc to the CCW anchor. Both stay small as the BFS only
    // adds edges OUTSIDE [cw_anchor, ccw_anchor].
    let cw_anchor_on_context = (first_junc + patch_n - third_pm.start_a) % patch_n;
    let ccw_anchor_pos = (third_pm.start_a + third_pm.len) % patch_n;
    let ccw_anchor_on_context = (ccw_anchor_pos + patch_n - last_junc) % patch_n;
    Some(NeighborhoodType {
        central_tile_id,
        cw_anchor_on_central,
        ccw_anchor_on_central,
        cw_anchor_on_context,
        ccw_anchor_on_context,
        num_ctx_tiles: 2,
        vt_seq,
    })
}

/// Phase 2: dequeue each NT, enumerate CCW one-petal successors via
/// [`explore_one`], dedup against `seen`, and record transitions.
/// Every step grows `num_ctx_tiles` by 1, so the transition graph is
/// a DAG.
fn bfs_phase<T: IsComplex + IsRingOrField + Units>(
    state: &mut BfsState,
    match_index: &Arc<MatchTypeIndex<T>>,
) {
    while let Some(nt) = state.queue.pop_front() {
        let src_id = *state.seen.get(&nt).expect("dequeued state must be in seen");
        for outcome in explore_one(&nt, match_index) {
            let dst_id = match outcome.kind {
                OutcomeKind::Closed => NT_CLOSED_ID,
                OutcomeKind::Open(new_nt) => {
                    if let Some(&id) = state.seen.get(&new_nt) {
                        id
                    } else {
                        let id = state.entries.len() + 1;
                        state.seen.insert(new_nt.clone(), id);
                        state.entries.push(new_nt.clone());
                        state.queue.push_back(new_nt);
                        id
                    }
                }
            };
            state.transitions.push(NtTransition {
                src_id,
                dst_id,
                side: outcome.side,
                tile_id: outcome.petal_pm.tile_id,
                tile_offset: outcome.petal_pm.start_b,
            });
        }
    }
}

struct ExploreOutcome {
    side: TransitionSide,
    petal_pm: PatchMatch,
    kind: OutcomeKind,
}

enum OutcomeKind {
    Closed,
    Open(NeighborhoodType),
}

struct AttachedContext<T: IsComplex> {
    aug: GrowingPatch<T>,
    central_ptid: usize,
    /// CCW frontier position on aug (gap_end vertex, always a junction).
    frontier_pos_on_aug: usize,
    /// Whether the CCW frontier vertex was already a junction on the
    /// pre-glue context boundary. Drives CCW Case A vs Case B in
    /// [`try_step_ccw`].
    frontier_is_junction_in_ctx: bool,
    /// Last context boundary edge covered by the central tile.
    /// Used as the `cw` edge of a freshly-appended vt in Case B.
    last_covered_ctx_edge: EdgeInfo,
    /// Number of context boundary edges on `aug` (= `ctx_n - match_len`).
    /// On `aug`, positions `[0, gap_start)` are surviving context edges
    /// and `[gap_start, aug_n)` are the central tile gap edges. Read
    /// only by tests; the production CCW step uses `frontier_pos_on_aug`
    /// instead.
    #[cfg_attr(not(test), allow(dead_code))]
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

fn explore_one<T: IsComplex + IsRingOrField + Units>(
    nt: &NeighborhoodType,
    match_index: &Arc<MatchTypeIndex<T>>,
) -> Vec<ExploreOutcome> {
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
///
/// **Edge-inclusive** check: `target_edge in [pm.start_a, pm.start_a +
/// pm.len - 1]` (mod n). Compare to [`crate::intgeom::patch::cyclic_range_contains`]
/// which is **vertex-inclusive** (covers `len + 1` vertices). Don't
/// substitute one for the other: e.g. for `start=25, len=1, n=26`,
/// `match_absorbs_edge(target=0) = false` (edge 0 is not in
/// `[25, 25]`) but `cyclic_range_contains(0) = true` (vertex 0 is
/// the CCW endpoint of edge 25).
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
) -> Option<ExploreOutcome> {
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
        kind: OutcomeKind::Open(new_nt),
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

/// Check that an NT is self-consistent: reconstruct from the CW-side
/// fields and verify the stored CCW-side fields agree with what the
/// geometry implies. Used as a post-hoc validator (see
/// [`NeighborhoodIndex::validate`]) and as a seed-acceptance filter
/// during BFS construction.
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

    /// Helper for the baked-in count tests: assert
    /// `(num_types, transitions, dead, undead, blessed, free,
    /// closed_transitions)` exactly. Catches silent drift in
    /// seed generation, BFS exploration, or classification.
    #[allow(clippy::too_many_arguments)]
    fn assert_counts<T: IsComplex + IsRingOrField + Units>(
        idx: &NeighborhoodIndex<T>,
        expected_types: usize,
        expected_transitions: usize,
        expected_closed: usize,
        expected_dead: usize,
        expected_undead: usize,
        expected_blessed: usize,
        expected_free: usize,
    ) {
        let kinds = idx.classify_all();
        let dead = kinds.iter().filter(|k| **k == NtKind::Dead).count();
        let undead = kinds.iter().filter(|k| **k == NtKind::Undead).count();
        let blessed = kinds.iter().filter(|k| **k == NtKind::Blessed).count();
        let free = kinds.iter().filter(|k| **k == NtKind::Free).count();
        let closed_trans = idx
            .transitions()
            .iter()
            .filter(|t| t.dst_id == NT_CLOSED_ID)
            .count();

        assert_eq!(idx.num_types(), expected_types, "num_types");
        assert_eq!(idx.transitions().len(), expected_transitions, "transitions");
        assert_eq!(closed_trans, expected_closed, "closed_trans");
        assert_eq!(dead, expected_dead, "dead");
        assert_eq!(undead, expected_undead, "undead");
        assert_eq!(blessed, expected_blessed, "blessed");
        assert_eq!(free, expected_free, "free");
        assert_eq!(dead + undead + blessed + free, expected_types, "kinds sum");
    }

    #[test]
    fn square_counts() {
        // Square tiles the plane convexly; every NT is Blessed (every
        // continuation closes).
        assert_counts(square_idx(), 109184, 436736, 327680, 0, 0, 109184, 0);
    }

    #[test]
    fn hex_counts() {
        assert_counts(hex_idx(), 55944, 335664, 279936, 0, 0, 55944, 0);
    }

    #[test]
    fn spectre_counts() {
        // Spectre is non-convex; exercises every kind.
        assert_counts(spectre_idx(), 10114, 11297, 1453, 6001, 1615, 584, 1914);
    }

    #[test]
    fn seeds_are_well_formed() {
        for idx in [
            square_idx() as &NeighborhoodIndex<ZZ12>,
            hex_idx(),
            spectre_idx(),
        ] {
            assert!(idx.num_types() > 0, "expected non-empty seed collection");
            let num_tiles = idx.tileset().num_tiles();
            for nt in idx.entries() {
                assert!(
                    nt.central_tile_id < num_tiles,
                    "central_tile_id {} out of range (have {} tiles)",
                    nt.central_tile_id,
                    num_tiles
                );
                assert!(!nt.vt_seq.is_empty(), "seeds must have non-empty vt_seq");
            }
        }
    }

    /// No two catalog entries are equal by `PartialEq`. This is the
    /// dedup contract of the BFS — not interface-level
    /// canonicalization. The catalog **does** contain multiple
    /// entries per interface class by design (different growth
    /// histories produce distinct `NeighborhoodType`s with the same
    /// match-side interface); see the module doc on the
    /// interface-quotient view.
    #[test]
    fn entries_have_no_partial_eq_duplicates() {
        for (name, idx) in [
            ("square", square_idx() as &NeighborhoodIndex<ZZ12>),
            ("hex", hex_idx()),
            ("spectre", spectre_idx()),
        ] {
            let mut keys: std::collections::HashSet<&NeighborhoodType> =
                std::collections::HashSet::new();
            for (i, nt) in idx.entries().iter().enumerate() {
                assert!(
                    keys.insert(nt),
                    "{}: duplicate entry at id {} (vt_seq={:?})",
                    name,
                    i + 1,
                    nt.vt_seq
                );
            }
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

    /// parse_file must reject a serialized catalog whose NTYPE kind
    /// disagrees with what classify_all re-derives from the parsed
    /// transitions. Catches stale catalogs from before a classifier
    /// fix.
    /// Pure unit test for [`match_absorbs_edge`]. The function appears
    /// in both production (`explore_one` filter) and in the brute
    /// side of `spectre_ccw_only_is_complete`, so a bug here would
    /// not be caught by the existing brute test (circular).
    /// Exhaustive cross-check against a brute reference for moderate
    /// `n`. The reference walks the edge set explicitly.
    #[test]
    fn match_absorbs_edge_unit() {
        fn brute(start: usize, len: usize, target: usize, n: usize) -> bool {
            if len == 0 || n == 0 {
                return false;
            }
            (0..len).any(|i| (start + i) % n == target)
        }

        for n in [1usize, 2, 5, 13, 26] {
            for start in 0..n {
                // len > n is a precondition violation (a match
                // cannot cover more edges than the boundary has).
                // Don't test those inputs.
                for len in 0..=n {
                    for target in 0..n {
                        let pm = PatchMatch {
                            start_a: start,
                            len,
                            start_b: 0,
                            tile_id: 0,
                        };
                        let got = match_absorbs_edge(&pm, target, n);
                        let want = brute(start, len, target, n);
                        assert_eq!(
                            got, want,
                            "n={n} start={start} len={len} target={target}"
                        );
                    }
                }
            }
        }

        // Pin the wrap regression case directly: edge n-1, len 1.
        let pm = PatchMatch {
            start_a: 25,
            len: 1,
            start_b: 0,
            tile_id: 0,
        };
        assert!(match_absorbs_edge(&pm, 25, 26));
        assert!(!match_absorbs_edge(&pm, 0, 26));

        // Pin the wrap-spanning case: starts at n-1, len 2, covers
        // edges n-1 and 0.
        let pm = PatchMatch {
            start_a: 25,
            len: 2,
            start_b: 0,
            tile_id: 0,
        };
        assert!(match_absorbs_edge(&pm, 25, 26));
        assert!(match_absorbs_edge(&pm, 0, 26));
        assert!(!match_absorbs_edge(&pm, 1, 26));
    }

    #[test]
    fn parse_file_rejects_kind_mismatch() {
        let idx = square_idx();
        let mut buf: Vec<u8> = Vec::new();
        idx.write_collection(&mut buf).expect("write");
        let text = std::str::from_utf8(&buf).expect("utf8");
        // Square's NTs are all Blessed; swap the first one's kind
        // to "dead" - that must fail to parse.
        let corrupted = text.replacen(" blessed ", " dead ", 1);
        assert_ne!(corrupted, text, "test fixture did not introduce a swap");
        let err = match NeighborhoodIndex::parse_file(Arc::clone(idx.tileset()), &corrupted) {
            Ok(_) => panic!("parser should reject kind mismatch"),
            Err(e) => e,
        };
        assert!(
            err.contains("kind mismatch"),
            "expected kind-mismatch error, got: {err}"
        );
    }

    fn spectre_tileset() -> Arc<TileSet<ZZ12>> {
        let sp = tiles::spectre::<ZZ12>();
        let rat = Rat::try_from(&sp).unwrap();
        Arc::new(TileSet::new(vec![rat]))
    }

    #[test]
    fn validate_returns_empty_for_all_tilesets() {
        for (name, idx) in [
            ("square", square_idx() as &NeighborhoodIndex<ZZ12>),
            ("hex", hex_idx()),
            ("spectre", spectre_idx()),
        ] {
            let bad = idx.validate();
            assert!(
                bad.is_empty(),
                "{} has {}/{} invalid entries (first few: {:?})",
                name,
                bad.len(),
                idx.num_types(),
                &bad[..bad.len().min(5)]
            );
        }
    }

    /// Canonical fingerprint of the **augmented patch boundary** of
    /// an NT (= context + central glued together, lex_min_rot of the
    /// resulting cyclic angle sequence, paired with the central tile
    /// id and cw_anchor_on_central).
    ///
    /// **This is not the interface fingerprint.** The interface
    /// fingerprint would only canonicalize the matched edge run on
    /// the context boundary (the linear segment shared with the
    /// central), discarding everything else. The fingerprint here is
    /// finer than the interface (it distinguishes patches whose
    /// non-matched perimeters differ) and coarser than `PartialEq`
    /// on `NeighborhoodType` (it ignores growth-history-dependent
    /// fields like `num_ctx_tiles`).
    ///
    /// Used by the brute-force completeness + correctness tests as
    /// an independent canonical form for comparing trial patches to
    /// catalog entries.
    fn aug_boundary_fingerprint<T: IsComplex + IsRingOrField + Units>(
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

    /// **Canonical-class fingerprint of an NT**: the equivalence
    /// class that exactly preserves the four-kind classification.
    /// Two NTs with the same value here are physically identical
    /// configurations (modulo translation and reflection-free
    /// rotation in the plane), so their futures of growth and hence
    /// their kinds are identical.
    ///
    /// Components:
    ///   - `central_tile_id` - which tile is the central.
    ///   - `cw_anchor_on_central` - which edge of the central is the
    ///     CW endpoint of the matched segment (= identifies the
    ///     "first central-side gap edge" at the anchor).
    ///   - The **canonical (lex_min_rot) aug-boundary angle sequence** -
    ///     captures the shape of the entire post-glue (ctx + central)
    ///     patch.
    ///   - The **CW anchor's position on that canonical boundary** -
    ///     pins where the central is attached on the canonical
    ///     shape. Without this, two NTs differing only by where the
    ///     central is glued on the same shape would collapse to the
    ///     same class but have different futures (different
    ///     neighboring un-matched-ctx geometry near the anchor).
    ///
    /// Empirically (`diagnose_canonical_class_compression`):
    /// kind-constancy across this fingerprint holds on hex, square,
    /// and spectre (0 violations on all three). Compression vs the
    /// BFS catalog ranges from ~1x (spectre, near-canonical) to
    /// ~2300x (hex, massively redundant).
    fn canonical_class_fingerprint<T: IsComplex + IsRingOrField + Units>(
        nt: &NeighborhoodType,
        mi: &Arc<MatchTypeIndex<T>>,
    ) -> Option<(usize, usize, Vec<i8>, usize)> {
        let ac = build_attached_context(nt, mi)?;
        let aug_n = ac.aug_n;
        let angles = ac.aug.angles().to_vec();
        let rot = crate::intgeom::rat::lex_min_rot(&angles);
        let mut canon = angles;
        canon.rotate_left(rot);
        // canon[i] = aug.angles()[(i + rot) % aug_n]. The CW anchor
        // vertex on aug lives at index ac.gap_start (= first gap
        // edge's starting vertex). Its position on the canonical
        // boundary is (gap_start - rot) mod aug_n.
        let cw_anchor_in_canon = (ac.gap_start + aug_n - rot) % aug_n;
        Some((
            nt.central_tile_id,
            nt.cw_anchor_on_central,
            canon,
            cw_anchor_in_canon,
        ))
    }

    /// Count catalog entries that share a canonical-class fingerprint
    /// but disagree on kind. Returns `(violations, first_violation,
    /// total_classes)`.
    fn count_canonical_class_kind_violations<T: IsComplex + IsRingOrField + Units>(
        idx: &NeighborhoodIndex<T>,
    ) -> (usize, Option<String>, usize) {
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        let kinds = idx.classify_all();
        let mut fp_to_kind: FxHashMap<
            (usize, usize, Vec<i8>, usize),
            (NtKind, usize),
        > = FxHashMap::default();
        let mut violations = 0usize;
        let mut first_violation: Option<String> = None;
        for (i, nt) in idx.entries().iter().enumerate() {
            let Some(fp) = canonical_class_fingerprint(nt, &mi) else {
                continue;
            };
            let kind = kinds[i];
            if let Some(&(prev_kind, prev_id)) = fp_to_kind.get(&fp) {
                if kind != prev_kind {
                    violations += 1;
                    if first_violation.is_none() {
                        first_violation = Some(format!(
                            "entry id {} kind={:?} shares canonical class with id {} kind={:?}",
                            i + 1,
                            kind,
                            prev_id,
                            prev_kind
                        ));
                    }
                }
            } else {
                fp_to_kind.insert(fp, (kind, i + 1));
            }
        }
        (violations, first_violation, fp_to_kind.len())
    }

    /// Pin the classification-preserving invariant: every catalog
    /// entry's kind is determined by its
    /// [`canonical_class_fingerprint`]. Two entries with the same
    /// fingerprint must classify the same. Holds on convex tilesets
    /// trivially (lots of redundancy in the BFS) and on non-convex
    /// tilesets non-trivially (the BFS is near-canonical but the
    /// frontier-position component matters).
    #[test]
    fn hex_kind_constant_per_canonical_class() {
        let (violations, msg, _) = count_canonical_class_kind_violations(hex_idx());
        assert_eq!(violations, 0, "{:?}", msg);
    }

    #[test]
    fn spectre_kind_constant_per_canonical_class() {
        let (violations, msg, _) = count_canonical_class_kind_violations(spectre_idx());
        assert_eq!(violations, 0, "{:?}", msg);
    }

    /// Square is the largest tileset; gated for cost.
    #[test]
    #[ignore]
    fn square_kind_constant_per_canonical_class() {
        let (violations, msg, _) = count_canonical_class_kind_violations(square_idx());
        assert_eq!(violations, 0, "{:?}", msg);
    }

    /// Data-only diagnostic: prints (entries, canonical-class count,
    /// compression ratio) per tileset. Useful for tracking the BFS's
    /// canonicalization quality over time without pinning fragile
    /// exact-count regressions. Run via
    /// `cargo test --release intgeom::neighborhood::tests::diagnose_canonical_class_compression -- --ignored --nocapture`.
    #[test]
    #[ignore]
    fn diagnose_canonical_class_compression() {
        for (name, idx) in [
            ("square", square_idx() as &NeighborhoodIndex<ZZ12>),
            ("hex", hex_idx()),
            ("spectre", spectre_idx()),
        ] {
            let (violations, _, n_classes) = count_canonical_class_kind_violations(idx);
            let n = idx.num_types();
            eprintln!(
                "{name}: entries={} canonical_classes={} ratio={:.2}x violations={}",
                n,
                n_classes,
                n as f64 / n_classes as f64,
                violations
            );
        }
    }

    /// Combined completeness + transition correctness check, using
    /// `rat.get_match` directly on each NT's augmented boundary as
    /// the independent reference.
    ///
    /// # Fingerprint level
    ///
    /// Comparison is at the **augmented-boundary fingerprint** level
    /// (see [`aug_boundary_fingerprint`]). That fingerprint is
    /// stricter than the interface fingerprint and looser than
    /// `PartialEq` on `NeighborhoodType`: it distinguishes patches
    /// whose full perimeters differ but treats two BFS entries with
    /// the same augmented patch as equivalent (= the BFS may produce
    /// multiple entries per geometric state via different
    /// `num_ctx_tiles` / `vt_seq` encodings, and we want to ignore
    /// that here).
    ///
    /// # Properties
    ///
    /// For each src, for each canonical petal placement
    /// (brute-enumerated via `rat.get_match`) absorbing the frontier
    /// edge:
    ///
    ///   1. **Completeness**: the trial's augmented-boundary
    ///      fingerprint is in the catalog, or the petal closed the
    ///      configuration. Catches "BFS missed a reachable state".
    ///   2. **Transition fingerprint coverage**: the set of distinct
    ///      destination fingerprints recorded for `(src, tile_id,
    ///      tile_offset)` includes the trial's destination
    ///      fingerprint. Catches "BFS recorded a transition pointing
    ///      to a geometrically-wrong destination".
    ///
    /// Failure modes NOT caught by this check:
    ///   - A bug at *interface* granularity where two distinct
    ///     interfaces collapse to the same augmented-boundary
    ///     fingerprint by accident (unlikely; geometrically rare).
    ///   - A bug at *id* granularity where two entries with the same
    ///     augmented fingerprint differ in some other respect that
    ///     matters to downstream consumers. (Currently no consumer
    ///     cares about anything other than the interface, so this
    ///     wouldn't matter.)
    ///
    /// Uses NO logic from `explore_one` / `try_step_*` /
    /// `compute_candidates_covering_position` - pure brute force.
    fn assert_brute_completeness_and_transition_correctness<T>(idx: &NeighborhoodIndex<T>)
    where
        T: IsComplex + IsRingOrField + Units,
    {
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        // fp_set: set of canonical fingerprints in the catalog.
        // id_to_fp: per-entry fingerprint for transition cross-check.
        let mut fp_set: FxHashMap<(usize, usize, Vec<i8>), ()> = FxHashMap::default();
        let mut id_to_fp: Vec<Option<(usize, usize, Vec<i8>)>> =
            Vec::with_capacity(idx.num_types());
        for nt in idx.entries() {
            let fp = aug_boundary_fingerprint(nt, &mi);
            if let Some(ref fp) = fp {
                fp_set.insert(fp.clone(), ());
            }
            id_to_fp.push(fp);
        }
        // recorded_by_key: (src_id, tile_id, tile_offset) -> set of
        // destination fingerprints (or "CLOSED" sentinel = None).
        type DstFp = Option<(usize, usize, Vec<i8>)>;
        let mut recorded_by_key: FxHashMap<
            (usize, usize, usize),
            std::collections::HashSet<DstFp>,
        > = FxHashMap::default();
        for t in idx.transitions() {
            let dst_fp: DstFp = if t.dst_id == NT_CLOSED_ID {
                None
            } else {
                id_to_fp[t.dst_id - 1].clone()
            };
            recorded_by_key
                .entry((t.src_id, t.tile_id, t.tile_offset))
                .or_default()
                .insert(dst_fp);
        }

        for (i, nt) in idx.entries().iter().enumerate() {
            let src_id = i + 1;
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
                        let (ns, len, ne) =
                            aug_rat.get_match((pos as i64, b_off as i64), other);
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

                        let trial_dst_fp: DstFp =
                            if !trial.patch_tile_ids().contains(&ac.central_ptid) {
                                None
                            } else {
                                let trial_angles = trial.angles().to_vec();
                                let rot = crate::intgeom::rat::lex_min_rot(&trial_angles);
                                let mut canon = trial_angles;
                                canon.rotate_left(rot);
                                let fp =
                                    (nt.central_tile_id, nt.cw_anchor_on_central, canon);
                                assert!(
                                    fp_set.contains_key(&fp),
                                    "completeness violation: src={} brute glue \
                                     (tile={}, start_a={}, len={}, start_b={}) \
                                     produces state not in catalog",
                                    src_id, tile_b, ns_u, len, ne_u
                                );
                                Some(fp)
                            };

                        let key = (src_id, tile_b, ne_u);
                        let recorded_set = recorded_by_key.get(&key);
                        assert!(
                            recorded_set.is_some_and(|s| s.contains(&trial_dst_fp)),
                            "transition fingerprint mismatch: src={} brute glue \
                             (tile={}, start_a={}, len={}, start_b={}) -> {} but \
                             catalog has no transition with that destination \
                             fingerprint for (src, tile, tile_offset)",
                            src_id, tile_b, ns_u, len, ne_u,
                            if trial_dst_fp.is_none() { "CLOSED".to_string() }
                            else { "Open(fp)".to_string() }
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn spectre_brute_complete_and_correct() {
        assert_brute_completeness_and_transition_correctness(spectre_idx());
    }

    #[test]
    fn hex_brute_complete_and_correct() {
        assert_brute_completeness_and_transition_correctness(hex_idx());
    }

    /// Square has 109k entries and 437k transitions; gated behind
    /// `--ignored` because runtime is order minutes.
    #[test]
    #[ignore]
    fn square_brute_complete_and_correct() {
        assert_brute_completeness_and_transition_correctness(square_idx());
    }



    /// Re-derive the four-kind classification predicates from
    /// `idx.transitions()` and assert they agree with `classify_all`.
    ///
    /// Cross-checks each kind via global reachability computed
    /// independently of `classify_all`'s algorithm:
    ///   - `reaches_closed[i]` = exists a forward path from i to
    ///     `NT_CLOSED_ID`
    ///   - `reaches_dead[i]` = exists a forward path from i to a
    ///     no-outgoing entry
    ///
    /// The four kinds are pinned by these reachability properties:
    ///   - Dead: no outgoing transitions
    ///   - Undead: has outgoing AND NOT reaches_closed
    ///   - Blessed: has outgoing AND NOT reaches_dead
    ///   - Free: has outgoing AND reaches_closed AND reaches_dead
    ///
    /// Catches the V2b shape (closing-as-escape blindspot) AND the
    /// "Free is the trash bucket" shape where Free becomes vacuous.
    fn assert_classify_invariants<T: IsComplex + IsRingOrField + Units>(
        idx: &NeighborhoodIndex<T>,
    ) {
        let kinds = idx.classify_all();
        let n = idx.num_types();

        let mut succs: Vec<Vec<usize>> = vec![Vec::new(); n];
        let mut preds: Vec<Vec<usize>> = vec![Vec::new(); n];
        let mut has_outgoing = vec![false; n];
        let mut closing_seeds: Vec<usize> = Vec::new();
        for t in idx.transitions() {
            let src = t.src_id - 1;
            has_outgoing[src] = true;
            if t.dst_id == NT_CLOSED_ID {
                closing_seeds.push(src);
            } else {
                succs[src].push(t.dst_id - 1);
                preds[t.dst_id - 1].push(src);
            }
        }

        let propagate = |seeds: Vec<usize>| -> Vec<bool> {
            let mut reached = vec![false; n];
            let mut queue: VecDeque<usize> = VecDeque::new();
            for s in seeds {
                if !reached[s] {
                    reached[s] = true;
                    queue.push_back(s);
                }
            }
            while let Some(v) = queue.pop_front() {
                for &p in &preds[v] {
                    if !reached[p] {
                        reached[p] = true;
                        queue.push_back(p);
                    }
                }
            }
            reached
        };

        let reaches_closed = propagate(closing_seeds);
        // Dead seeds for the reaches_dead reachability: NTs with no
        // outgoing transitions are themselves Dead. A path to Dead
        // means a path that traps. (Note: a Dead NT trivially reaches
        // itself, so reaches_dead[Dead] = true for all Dead.)
        let dead_seeds: Vec<usize> = (0..n).filter(|&i| !has_outgoing[i]).collect();
        let reaches_dead = propagate(dead_seeds);

        for (i, &kind) in kinds.iter().enumerate() {
            let id = i + 1;
            let expected = if !has_outgoing[i] {
                NtKind::Dead
            } else if !reaches_closed[i] {
                NtKind::Undead
            } else if !reaches_dead[i] {
                NtKind::Blessed
            } else {
                NtKind::Free
            };
            assert_eq!(
                kind, expected,
                "NT id {id} classify_all says {:?} but reachability says {:?} \
                 (has_outgoing={} reaches_closed={} reaches_dead={})",
                kind, expected, has_outgoing[i], reaches_closed[i], reaches_dead[i]
            );
        }
    }

    #[test]
    fn square_classify_invariants() {
        assert_classify_invariants(square_idx());
    }

    #[test]
    fn hex_classify_invariants() {
        assert_classify_invariants(hex_idx());
    }

    #[test]
    fn spectre_classify_invariants() {
        assert_classify_invariants(spectre_idx());
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

    /// Generic over the tileset: build_attached_context succeeds for
    /// every seed, the gap is non-empty, and gap_start is in range.
    fn assert_seed_geometry<T: IsComplex + IsRingOrField + Units>(
        idx: &NeighborhoodIndex<T>,
    ) {
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        for (i, nt) in idx.entries().iter().enumerate() {
            let ac = build_attached_context(nt, &mi)
                .unwrap_or_else(|| panic!("seed {}: build_attached_context failed for nt {:?}", i, nt));
            let gap_len = ac.aug_n - ac.gap_start;
            assert!(gap_len > 0, "seed {}: gap should be non-empty", i);
            assert!(ac.gap_start < ac.aug_n, "seed {}: gap_start out of range", i);
            let central_len = idx.tileset().rat(nt.central_tile_id).seq().len();
            assert!(
                gap_len < central_len,
                "seed {}: gap_len {} should be < central_len {}",
                i,
                gap_len,
                central_len
            );
            assert!(
                !explore_one(nt, &mi).is_empty(),
                "seed {}: should have at least one successful petal",
                i
            );
        }
    }

    #[test]
    fn square_seed_geometry() {
        assert_seed_geometry(square_idx());
    }

    #[test]
    fn hex_seed_geometry() {
        assert_seed_geometry(hex_idx());
    }

    /// Brute-enumerate every legal 3-tile seed configuration via
    /// `rat.get_match` directly (independent of `MatchTypeIndex` and
    /// `patch.get_all_matches`, which `seed_phase` uses internally),
    /// canonicalize each, and assert the resulting set equals the
    /// catalog's seed entries (= entries with `num_ctx_tiles == 2`).
    ///
    /// Catches seed-completeness regressions of the V2c shape: if
    /// the production match enumerator misses canonical matches, the
    /// seed phase would silently under-produce, and downstream
    /// closure tests (which iterate over existing entries) couldn't
    /// recover the missing seeds.
    fn assert_seeds_match_brute<T: IsComplex + IsRingOrField + Units>(
        tileset: Arc<TileSet<T>>,
        idx: &NeighborhoodIndex<T>,
    ) {
        type Class = (usize, usize, Vec<i8>, usize);
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(&tileset)));

        // Brute side: enumerate canonical 3-tile patches via direct
        // get_match, applying the same filters seed_phase does.
        let mut brute: FxHashMap<Class, ()> = FxHashMap::default();
        let num_tiles = tileset.num_tiles();
        for a in 0..num_tiles {
            let tile_a = tileset.rat(a);
            let n_a = tile_a.seq().len();
            for b in 0..num_tiles {
                let tile_b = tileset.rat(b);
                let n_b = tile_b.seq().len();
                let mut seen_ab: FxHashMap<(usize, usize, usize), ()> = FxHashMap::default();
                for ia in 0..n_a {
                    for ib in 0..n_b {
                        let (ns, len, ne) = tile_a.get_match((ia as i64, ib as i64), tile_b);
                        if len == 0 {
                            continue;
                        }
                        let ns_u = ns.rem_euclid(n_a as i64) as usize;
                        let ne_u = ne.rem_euclid(n_b as i64) as usize;
                        if seen_ab.contains_key(&(ns_u, len, ne_u)) {
                            continue;
                        }
                        seen_ab.insert((ns_u, len, ne_u), ());
                        if !crate::intgeom::matchtypes::junction_gap_nonnegative(
                            tile_a.seq(),
                            ns_u,
                            len,
                            tile_b.seq(),
                            ne_u,
                        ) {
                            continue;
                        }
                        if tile_a
                            .try_glue_precomputed((ns, len, ne), tile_b, true)
                            .is_err()
                        {
                            continue;
                        }
                        // Build the 2-tile patch.
                        let mut patch = GrowingPatch::new(Arc::clone(&tileset), a);
                        let pm1 = PatchMatch {
                            start_a: ns_u,
                            len,
                            start_b: ne_u,
                            tile_id: b,
                        };
                        if !patch.add_tile(&pm1) {
                            continue;
                        }
                        patch.normalize();
                        let patch_n = patch.boundary_len();
                        let patch_rat = crate::intgeom::rat::Rat::from_slice_unchecked(
                            patch.angles(),
                        );

                        // Brute third-tile candidates over the 2-tile boundary.
                        let mut seen_third: FxHashMap<(usize, usize, usize, usize), ()> =
                            FxHashMap::default();
                        for c in 0..num_tiles {
                            let tile_c = tileset.rat(c);
                            let n_c = tile_c.seq().len();
                            for ic in 0..patch_n {
                                for jc in 0..n_c {
                                    let (ns3, len3, ne3) = patch_rat.get_match(
                                        (ic as i64, jc as i64),
                                        tile_c,
                                    );
                                    if len3 == 0 {
                                        continue;
                                    }
                                    let ns3_u = ns3.rem_euclid(patch_n as i64) as usize;
                                    let ne3_u = ne3.rem_euclid(n_c as i64) as usize;
                                    if seen_third.contains_key(&(ns3_u, len3, ne3_u, c)) {
                                        continue;
                                    }
                                    seen_third.insert((ns3_u, len3, ne3_u, c), ());
                                    let pm3 = PatchMatch {
                                        start_a: ns3_u,
                                        len: len3,
                                        start_b: ne3_u,
                                        tile_id: c,
                                    };
                                    let mut trial = patch.clone();
                                    if !trial.add_tile(&pm3) {
                                        continue;
                                    }
                                    // Filter: match must touch at least one
                                    // junction of the 2-tile patch (else
                                    // build_seed_nt returns None).
                                    let touches_junction = (0..patch_n).any(|jp| {
                                        if !patch.is_junction(jp) {
                                            return false;
                                        }
                                        (jp + patch_n - ns3_u) % patch_n <= len3
                                    });
                                    if !touches_junction {
                                        continue;
                                    }
                                    // Compute canonical class fingerprint.
                                    let aug_angles = trial.angles().to_vec();
                                    let aug_n = aug_angles.len();
                                    let rot = crate::intgeom::rat::lex_min_rot(&aug_angles);
                                    let mut canon = aug_angles;
                                    canon.rotate_left(rot);
                                    let gap_start = patch_n - len3;
                                    let cw_anchor_in_canon =
                                        (gap_start + aug_n - rot) % aug_n;
                                    brute.insert(
                                        (c, ne3_u, canon, cw_anchor_in_canon),
                                        (),
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }

        // Catalog side: canonical class of every entry with
        // num_ctx_tiles == 2 (= seed entries).
        let mut catalog: FxHashMap<Class, ()> = FxHashMap::default();
        for nt in idx.entries() {
            if nt.num_ctx_tiles != 2 {
                continue;
            }
            if let Some(fp) = canonical_class_fingerprint(nt, &mi) {
                catalog.insert(fp, ());
            }
        }

        let brute_only: Vec<&Class> = brute.keys().filter(|k| !catalog.contains_key(*k)).collect();
        let catalog_only: Vec<&Class> = catalog.keys().filter(|k| !brute.contains_key(*k)).collect();
        assert!(
            brute_only.is_empty(),
            "seed completeness: {} canonical seed classes found by brute but missing from catalog",
            brute_only.len()
        );
        assert!(
            catalog_only.is_empty(),
            "seed soundness: {} catalog seed classes not produced by brute",
            catalog_only.len()
        );
    }

    #[test]
    fn hex_seeds_match_brute() {
        assert_seeds_match_brute(Arc::clone(hex_idx().tileset()), hex_idx());
    }

    #[test]
    fn square_seeds_match_brute() {
        assert_seeds_match_brute(Arc::clone(square_idx().tileset()), square_idx());
    }

    #[test]
    fn spectre_seeds_match_brute() {
        assert_seeds_match_brute(Arc::clone(spectre_idx().tileset()), spectre_idx());
    }

    #[test]
    fn gap_edges_are_central_tile() {
        let idx = square_idx();
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        for (i, nt) in idx.entries().iter().enumerate() {
            let ac = build_attached_context(nt, &mi)
                .unwrap_or_else(|| panic!("seed {}: build_attached_context failed", i));
            let edges = ac.aug.edges();
            let ptids = ac.aug.patch_tile_ids();
            let central_ptid = ptids[ac.gap_start];
            for k in ac.gap_start..ac.aug_n {
                assert_eq!(
                    edges[k].tile_id, nt.central_tile_id,
                    "seed {}: gap edge at pos {} should be central tile",
                    i, k
                );
                assert_eq!(
                    ptids[k], central_ptid,
                    "seed {}: gap edge at pos {} should have same ptid",
                    i, k
                );
            }
        }
    }

    #[test]
    fn ccw_frontier_at_gap_end() {
        // After `build_attached_context`, the boundary is rotated so
        // gap_end == 0. The CCW frontier vertex must therefore be at
        // position 0 - confirming gap_end is itself a junction.
        let idx = hex_idx();
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        for (i, nt) in idx.entries().iter().enumerate() {
            let ac = build_attached_context(nt, &mi)
                .unwrap_or_else(|| panic!("seed {}: build_attached_context failed", i));
            assert_eq!(
                ac.frontier_pos_on_aug, 0,
                "seed {}: gap_end (pos 0 on aug) should be the CCW frontier junction",
                i
            );
            assert!(ac.aug.is_junction(0), "seed {}: position 0 must be junction", i);
        }
    }

    #[test]
    fn explore_one_square_has_petals() {
        let idx = square_idx();
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        let nt = &idx.entries()[0];
        assert!(
            !explore_one(nt, &mi).is_empty(),
            "should have at least one petal outcome"
        );
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
                    OutcomeKind::Open(nt) => nt.clone(),
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
        assert!(checked > 0, "should have validated at least one new NT");
        assert!(
            errors.is_empty(),
            "{} reattach errors (first few): {}",
            errors.len(),
            errors[..errors.len().min(10)].join("; ")
        );
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
                    OutcomeKind::Open(nt) => nt.clone(),
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
        assert!(checked > 0, "should have validated at least one new NT");
        assert!(
            errors.is_empty(),
            "{} reattach errors (first few): {}",
            errors.len(),
            errors[..errors.len().min(10)].join("; ")
        );
    }
}
