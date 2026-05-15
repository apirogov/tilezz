//! Neighborhood types: catalog of local *interfaces* between a tile
//! and any legal context patch that could meet it along one
//! contiguous edge run, plus the post-absorption "surrounded tile"
//! states that arise once the central tile's edges are all
//! consumed.
//!
//! # Two phases
//!
//! The catalog has two kinds of entries (unified as [`NtEntry`]):
//!
//! - **Phase 1 — [`NeighborhoodType`]**: the central tile is still
//!   on the patch boundary along a matched segment, with CW/CCW
//!   frontier vertices where context tiles can still be added.
//!   This is the original "open NT" model.
//!
//! - **Phase 2 — [`SurroundedTile`]**: a phase-1 step absorbed the
//!   last central edge, so the central is now interior; only its
//!   vertices may remain incident with the boundary (= an "open
//!   vertex" where the angle gap hasn't closed). The full
//!   `GrowingPatch` boundary state is stored, in normalized form;
//!   see [`SurroundedTile::to_patch`] for reconstruction.
//!
//! Phase-2 entries either:
//!   - **closed** (no open vertex remaining; central fully
//!     surrounded — terminal Blessed state), or
//!   - **open** (one open vertex remains; phase-2 BFS attempts to
//!     close it by gluing further petals).
//!
//! # Why grow incrementally
//!
//! A naive "enumerate all length-≤k boundary strings" approach has
//! to separately filter realizable from junk strings — a costly
//! secondary step. Growing the patch one tile at a time via BFS
//! gives realizability for free: every entry is built from a
//! known-realizable predecessor, so no candidate has to be checked
//! against geometric legality post-hoc.
//!
//! The trade-off: the BFS catalog is **finer** than the interface
//! equivalence quotient. The same configuration can be reached by
//! multiple distinct growth histories, and each history becomes its
//! own catalog entry. The true interface quotient is a coarser,
//! derived view; the canonicalisation pass that runs after each
//! glue (via [`GrowingPatch::normalize`]) deduplicates aggressively
//! but not exhaustively.
//!
//! # Classification → forbidden patterns
//!
//! Every entry gets a `Dead`/`Undead`/`Blessed`/`Free` reachability
//! kind via [`NeighborhoodIndex::classify_all`]. The Dead and
//! Undead interfaces are the actionable output: they are the local
//! patterns a tiling can never recover from. Downstream consumers
//! treat them as **forbidden patterns** that must not appear in any
//! extending glue. Together with the dead/undead VTs (a finer
//! per-vertex constraint, see `vertextypes.rs`), these constraints
//! prune the tile-growth search.
//!
//! See [`NeighborhoodType`] / [`SurroundedTile`] for entry layout
//! and [`NeighborhoodIndex::new`] for the BFS algorithm.

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

/// Phase-2 entry: the central tile has had all of its edges absorbed
/// by surrounding tiles ("surrounded"). Stores the full `GrowingPatch`
/// boundary state, in **normalized** form (= rotated to lex-min
/// boundary angles, `patch_tile_ids` densified to `0..k` by
/// CCW-first-occurrence). Two SurroundedTiles describing geometrically
/// equivalent patches compare equal via derived `Eq`/`Hash`.
///
/// `is_closed = true` when the corona fully encloses the central
/// (no boundary vertex incident with central). The entry is a
/// terminal Blessed state — it has no outgoing transitions, and
/// [`NeighborhoodIndex::classify_all`] recognises it as an
/// escape destination for the entries that transition into it.
///
/// `is_closed = false` when an open vertex remains on the boundary.
/// `open_vertex_pos` is the boundary position of that open vertex
/// in the normalized boundary (= where phase-2 BFS attaches petals).
///
/// To reconstruct as a `GrowingPatch`, use [`SurroundedTile::to_patch`].
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct SurroundedTile {
    pub central_tile_id: usize,
    pub is_closed: bool,
    /// Open vertex's boundary position in the normalized boundary
    /// (only meaningful when `!is_closed`; set to 0 by convention
    /// for closed entries).
    pub open_vertex_pos: usize,
    /// Number of distinct boundary tile instances (= `next_tile_id`
    /// after [`GrowingPatch::normalize`]). Fully-interior tiles are
    /// not counted here — see normalize's note.
    pub num_tile_instances: usize,
    pub angles: Vec<i8>,
    pub edges: Vec<EdgeInfo>,
    pub inner_chains: Vec<Vec<EdgeInfo>>,
    pub patch_tile_ids: Vec<usize>,
}

impl SurroundedTile {
    /// Reconstruct the surrounded-tile patch as a `GrowingPatch`. Returns
    /// `None` only on `from_parts` invariant violations (= shouldn't
    /// happen for SurroundedTiles produced by the BFS).
    ///
    /// In debug builds, asserts the internal-consistency invariants
    /// callers must uphold:
    /// - Per-position arrays (`angles`, `edges`, `inner_chains`,
    ///   `patch_tile_ids`) are the same length.
    /// - `num_tile_instances` equals `max(patch_tile_ids) + 1` (= the
    ///   dense-renumbering output of [`GrowingPatch::normalize`]).
    /// - `open_vertex_pos < angles.len()` when `!is_closed`.
    pub fn to_patch<T: IsComplex + IsRingOrField + Units>(
        &self,
        match_index: Arc<MatchTypeIndex<T>>,
    ) -> Option<GrowingPatch<T>> {
        let n = self.angles.len();
        debug_assert_eq!(self.edges.len(), n, "edges length mismatch");
        debug_assert_eq!(self.inner_chains.len(), n, "inner_chains length mismatch");
        debug_assert_eq!(
            self.patch_tile_ids.len(),
            n,
            "patch_tile_ids length mismatch"
        );
        debug_assert_eq!(
            self.num_tile_instances,
            self.patch_tile_ids.iter().max().map_or(0, |m| m + 1),
            "num_tile_instances doesn't match patch_tile_ids"
        );
        debug_assert!(
            self.is_closed || self.open_vertex_pos < n,
            "open_vertex_pos {} out of range for boundary len {}",
            self.open_vertex_pos,
            n
        );
        GrowingPatch::from_parts(
            match_index,
            self.angles.clone(),
            self.edges.clone(),
            self.inner_chains.clone(),
            self.patch_tile_ids.clone(),
            self.num_tile_instances,
        )
    }
}

/// A catalog entry: either a phase-1 NT (central still on the
/// patch boundary along a matched segment) or a phase-2 SurroundedTile
/// (central fully edge-absorbed, possibly with an open vertex).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum NtEntry {
    Phase1(NeighborhoodType),
    Phase2(SurroundedTile),
}

impl NtEntry {
    pub fn central_tile_id(&self) -> usize {
        match self {
            NtEntry::Phase1(nt) => nt.central_tile_id,
            NtEntry::Phase2(st) => st.central_tile_id,
        }
    }

    /// For phase-1 entries: the matched-segment `vt_seq`. For phase-2
    /// entries: an empty slice (the surrounded-tile representation is
    /// now boundary state, not a junction sequence; query
    /// [`SurroundedTile::to_patch`] and `junction_vertex_type_at` for
    /// per-junction info).
    pub fn vt_seq(&self) -> &[OpenVertexType] {
        match self {
            NtEntry::Phase1(nt) => &nt.vt_seq,
            NtEntry::Phase2(_) => &[],
        }
    }

    pub fn as_phase1(&self) -> Option<&NeighborhoodType> {
        match self {
            NtEntry::Phase1(nt) => Some(nt),
            _ => None,
        }
    }

    pub fn as_phase2(&self) -> Option<&SurroundedTile> {
        match self {
            NtEntry::Phase2(st) => Some(st),
            _ => None,
        }
    }
}

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
    entries: Vec<NtEntry>,
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
    /// 2. **BFS growth.** Dequeue each entry id and dispatch:
    ///    - [`NtEntry::Phase1`]: enumerate one-petal CCW successors
    ///      via [`explore_phase1`]. Each successor is either another
    ///      phase-1 NT (central edges still remain) or a phase-2
    ///      [`SurroundedTile`] (central edges all absorbed).
    ///    - [`NtEntry::Phase2`] (open only): explore one petal at
    ///      the open vertex via [`explore_phase2`]; closed phase-2
    ///      entries are terminal Blessed states with no outgoing
    ///      transitions ([`Self::classify_all`] recognises them
    ///      directly as escape destinations).
    ///    New entries are enqueued; transitions are recorded.
    ///    Every step adds exactly one tile to the underlying patch,
    ///    so the transition graph is a DAG.
    ///
    /// # CCW-only exploration
    ///
    /// Phase-1 BFS only explores in the CCW direction; the
    /// `spectre_ccw_only_is_complete` regression test confirms this
    /// is complete (= every brute-force petal candidate at the CCW
    /// frontier either closes the central or produces an NT already
    /// in the catalog). Phase-2 BFS likewise only attaches at the
    /// open vertex's CCW edge.
    pub fn new(tileset: Arc<TileSet<T>>) -> Self {
        let match_index = Arc::new(MatchTypeIndex::new(Arc::clone(&tileset)));
        // IDs are 1-based: entries[i] has id i + 1. `seen` and
        // transitions store ids, not indices.
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

    /// All catalog entries in BFS-discovery order.
    /// `entries()[id - 1]` is the entry for 1-based id `id`.
    pub fn entries(&self) -> &[NtEntry] {
        &self.entries
    }

    /// Number of catalog entries (phase 1 + phase 2). Ids run
    /// `1..=num_types()`.
    pub fn num_types(&self) -> usize {
        self.entries.len()
    }

    /// Number of phase-1 (open-NT) entries.
    pub fn num_phase1(&self) -> usize {
        self.entries
            .iter()
            .filter(|e| matches!(e, NtEntry::Phase1(_)))
            .count()
    }

    /// Number of phase-2 (SurroundedTile) entries.
    pub fn num_phase2(&self) -> usize {
        self.entries
            .iter()
            .filter(|e| matches!(e, NtEntry::Phase2(_)))
            .count()
    }

    /// All recorded BFS transitions, in BFS-emission order. Both
    /// `src_id` and `dst_id` are 1-based entry ids.
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
    /// * `Dead` — no outgoing transitions, not a closed phase-2 terminal.
    /// * `Undead` — has outgoing transitions, every reachable path
    ///   eventually traps at a Dead/Undead entry (never reaches a
    ///   closed phase-2 terminal).
    /// * `Blessed` — every outgoing transition reaches a closed
    ///   phase-2 terminal or another Blessed entry. Closed phase-2
    ///   terminals are themselves Blessed.
    /// * `Free` — has at least one path that reaches a closed
    ///   terminal AND at least one that traps at Dead/Undead.
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
        // these per-entry vectors.
        //
        // Closed phase-2 SurroundedTile entries are "terminal-good" —
        // they have no outgoing transitions but represent a valid
        // closure (= central fully surrounded). They are the *only*
        // Blessed seeds. All Blessed and Undead status reaches other
        // entries via standard reverse-edge propagation.
        let n = self.entries.len();
        let mut is_closed_terminal = vec![false; n];
        for (i, entry) in self.entries.iter().enumerate() {
            if let NtEntry::Phase2(st) = entry {
                if st.is_closed {
                    is_closed_terminal[i] = true;
                }
            }
        }

        let mut succ_count = vec![0usize; n];
        let mut preds: Vec<Vec<usize>> = vec![vec![]; n];
        let mut has_outgoing = vec![false; n];
        for t in &self.transitions {
            let src_idx = t.src_id - 1;
            let dst_idx = t.dst_id - 1;
            has_outgoing[src_idx] = true;
            succ_count[src_idx] += 1;
            preds[dst_idx].push(src_idx);
        }

        // Cursed propagation: an entry becomes cursed when all of its
        // successors are cursed. Closed terminals are never cursed
        // (they're seeded as Blessed instead), so any entry with a
        // closed-terminal successor has remaining > 0 permanently —
        // the "has-closing-escape" exemption is implicit.
        let mut cursed = vec![false; n];
        let mut remaining = succ_count.clone();
        let mut queue: VecDeque<usize> = VecDeque::new();
        for i in 0..n {
            if !has_outgoing[i] && !is_closed_terminal[i] {
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

        // Blessed propagation: seed with closed terminals only.
        // Standard reverse-edge propagation marks each predecessor
        // Blessed once all of its successors are Blessed.
        let mut blessed = vec![false; n];
        let mut remaining = succ_count.clone();
        let mut queue: VecDeque<usize> = VecDeque::new();
        for i in 0..n {
            if is_closed_terminal[i] {
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

    /// Return the 1-based ids of entries that fail self-consistency:
    /// phase-1 entries that don't reconstruct from their CW-side
    /// fields. Phase-2 entries are validated by attempting
    /// [`SurroundedTile::to_patch`] — `from_parts` enforces the
    /// shape's basic invariants.
    pub fn validate(&self) -> Vec<usize> {
        let match_index = Arc::new(MatchTypeIndex::new(Arc::clone(&self.tileset)));
        let num_tiles = self.tileset.num_tiles();
        self.entries
            .iter()
            .enumerate()
            .filter_map(|(idx, entry)| {
                let ok = match entry {
                    NtEntry::Phase1(nt) => nt_is_valid(nt, &match_index),
                    NtEntry::Phase2(st) => {
                        st.central_tile_id < num_tiles
                            && !st.angles.is_empty()
                            && st.to_patch(Arc::clone(&match_index)).is_some()
                    }
                };
                if ok {
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
    /// NTYPE <id> <kind> <central_id> <cwc> <ccwc> <cwx> <ccwx> <nctx> <n_vts> <vt0> ...
    /// STYPE <id> <kind> <central_id> <is_closed:0|1> <open_vertex_pos> <num_tile_instances> <n>
    ///   <angle0> ... <angle_{n-1}>
    ///   <edge0> ... <edge_{n-1}>
    ///   <ptid0> ... <ptid_{n-1}>
    ///   <ichain0> ... <ichain_{n-1}>
    /// TRANS <src_id> <dst_id> <side> <tile_id> <tile_offset>
    /// ```
    ///
    /// STYPE records use a multi-line format because the patch state
    /// is dense; the inner-chain entries are encoded
    /// `<count>:<edge0>,<edge1>,...` or just `-` for empty. Each edge is
    /// `<tile_id>.<tile_offset>`.
    ///
    /// NTYPE: each `<vti>` encodes a single OpenVertexType as
    /// `<cw_id>.<cw_off>|<inner_list>|<ccw_id>.<ccw_off>` where
    /// `<inner_list>` is either `-` (empty) or a comma-separated list of
    /// `<tile_id>.<tile_offset>` edges. `<kind>` is one of `dead`,
    /// `undead`, `blessed`, `free`. The tileset is not serialized: callers
    /// must pass a matching tileset to `parse_file`.
    pub fn write_collection(&self, out: &mut impl std::io::Write) -> std::io::Result<()> {
        let kinds = self.classify_all();
        for (idx, entry) in self.entries.iter().enumerate() {
            let id = idx + 1;
            let kind = match kinds[idx] {
                NtKind::Dead => "dead",
                NtKind::Undead => "undead",
                NtKind::Blessed => "blessed",
                NtKind::Free => "free",
            };
            match entry {
                NtEntry::Phase1(nt) => {
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
                    write_vt_seq(out, &nt.vt_seq)?;
                    writeln!(out)?;
                }
                NtEntry::Phase2(st) => {
                    let n = st.angles.len();
                    write!(
                        out,
                        "STYPE {} {} {} {} {} {} {}",
                        id,
                        kind,
                        st.central_tile_id,
                        if st.is_closed { 1 } else { 0 },
                        st.open_vertex_pos,
                        st.num_tile_instances,
                        n,
                    )?;
                    for &a in &st.angles {
                        write!(out, " {}", a)?;
                    }
                    for e in &st.edges {
                        write!(out, " {}.{}", e.tile_id, e.tile_offset)?;
                    }
                    for &p in &st.patch_tile_ids {
                        write!(out, " {}", p)?;
                    }
                    for ic in &st.inner_chains {
                        if ic.is_empty() {
                            write!(out, " -")?;
                        } else {
                            write!(out, " {}", ic.len())?;
                            for e in ic {
                                write!(out, ",{}.{}", e.tile_id, e.tile_offset)?;
                            }
                        }
                    }
                    writeln!(out)?;
                }
            }
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
    /// Validation: NTYPE/STYPE ids are 1-based and dense (`1, 2, ...,
    /// num_entries`), TRANS `src_id`/`dst_id` are in range, and the
    /// recorded entry `kind` agrees with what [`Self::classify_all`]
    /// re-derives from the parsed transitions. Mismatched kinds catch
    /// stale serialized catalogs after a classification change.
    pub fn parse_file(tileset: Arc<TileSet<T>>, input: &str) -> Result<Self, String> {
        let mut entries: Vec<NtEntry> = Vec::new();
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
                    entries.push(NtEntry::Phase1(nt));
                    recorded_kinds.push(parsed_kind);
                }
                "STYPE" => {
                    let (parsed_id, parsed_kind, st) = parse_stype_line(&mut tok, lineno + 1)?;
                    let expected_id = entries.len() + 1;
                    if parsed_id != expected_id {
                        return Err(format!(
                            "line {}: STYPE id sequence broken (expected {}, got {})",
                            lineno + 1,
                            expected_id,
                            parsed_id
                        ));
                    }
                    entries.push(NtEntry::Phase2(st));
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
                    if dst_id < 1 || dst_id > entries.len() {
                        return Err(format!(
                            "line {}: TRANS dst_id {} out of range (have {} entries; \
                             dst_id must be 1-based)",
                            lineno + 1,
                            dst_id,
                            entries.len()
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
        for (i, (recorded, computed)) in
            recorded_kinds.iter().zip(computed_kinds.iter()).enumerate()
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

fn parse_stype_line(
    tok: &mut std::str::SplitWhitespace<'_>,
    lineno: usize,
) -> Result<(usize, NtKind, SurroundedTile), String> {
    let mut next = |what: &str| -> Result<String, String> {
        tok.next()
            .map(|s| s.to_string())
            .ok_or_else(|| format!("line {}: missing {}", lineno, what))
    };
    let parsed_id: usize = next("id")?
        .parse()
        .map_err(|e| format!("line {}: bad id: {}", lineno, e))?;
    let parsed_kind = match next("kind")?.as_str() {
        "dead" => NtKind::Dead,
        "undead" => NtKind::Undead,
        "blessed" => NtKind::Blessed,
        "free" => NtKind::Free,
        other => return Err(format!("line {}: unknown kind `{}`", lineno, other)),
    };
    let central_tile_id: usize = next("central_id")?
        .parse()
        .map_err(|e| format!("line {}: bad central_id: {}", lineno, e))?;
    let is_closed_raw: u32 = next("is_closed")?
        .parse()
        .map_err(|e| format!("line {}: bad is_closed: {}", lineno, e))?;
    let is_closed = is_closed_raw != 0;
    let open_vertex_pos: usize = next("open_vertex_pos")?
        .parse()
        .map_err(|e| format!("line {}: bad open_vertex_pos: {}", lineno, e))?;
    let num_tile_instances: usize = next("num_tile_instances")?
        .parse()
        .map_err(|e| format!("line {}: bad num_tile_instances: {}", lineno, e))?;
    let n: usize = next("n")?
        .parse()
        .map_err(|e| format!("line {}: bad n: {}", lineno, e))?;
    let mut angles = Vec::with_capacity(n);
    for i in 0..n {
        let s = next(&format!("angle[{}]", i))?;
        angles.push(
            s.parse::<i8>()
                .map_err(|e| format!("line {}: bad angle[{}]: {}", lineno, i, e))?,
        );
    }
    let mut edges = Vec::with_capacity(n);
    for i in 0..n {
        let s = next(&format!("edge[{}]", i))?;
        edges.push(parse_edge_info(&s).map_err(|e| format!("line {}: {}", lineno, e))?);
    }
    let mut patch_tile_ids = Vec::with_capacity(n);
    for i in 0..n {
        let s = next(&format!("ptid[{}]", i))?;
        patch_tile_ids.push(
            s.parse::<usize>()
                .map_err(|e| format!("line {}: bad ptid[{}]: {}", lineno, i, e))?,
        );
    }
    let mut inner_chains = Vec::with_capacity(n);
    for i in 0..n {
        let s = next(&format!("ichain[{}]", i))?;
        if s == "-" {
            inner_chains.push(Vec::new());
        } else {
            let mut parts = s.splitn(2, ',');
            let count_str = parts
                .next()
                .ok_or_else(|| format!("line {}: bad ichain[{}]: empty", lineno, i))?;
            let _count: usize = count_str
                .parse()
                .map_err(|e| format!("line {}: bad ichain[{}] count: {}", lineno, i, e))?;
            let rest = parts.next().unwrap_or("");
            let chain: Vec<EdgeInfo> = if rest.is_empty() {
                Vec::new()
            } else {
                rest.split(',')
                    .map(parse_edge_info)
                    .collect::<Result<Vec<_>, _>>()
                    .map_err(|e| format!("line {}: {}", lineno, e))?
            };
            inner_chains.push(chain);
        }
    }
    Ok((
        parsed_id,
        parsed_kind,
        SurroundedTile {
            central_tile_id,
            is_closed,
            open_vertex_pos,
            num_tile_instances,
            angles,
            edges,
            inner_chains,
            patch_tile_ids,
        },
    ))
}

fn write_vt_seq(out: &mut impl std::io::Write, vt_seq: &[OpenVertexType]) -> std::io::Result<()> {
    for vt in vt_seq {
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
    Ok(())
}

/// Mutable state shared across the seed and BFS phases of
/// [`NeighborhoodIndex::new`].
struct BfsState {
    /// All discovered entries in BFS order. `entries[i]` has 1-based id
    /// `i + 1`. Phase-1 and phase-2 entries share one id space.
    entries: Vec<NtEntry>,
    /// All recorded transitions between entries.
    transitions: Vec<NtTransition>,
    /// `seen[entry]` = 1-based id of the entry in `entries`.
    seen: FxHashMap<NtEntry, usize>,
    /// BFS frontier of entry ids waiting to be explored.
    queue: VecDeque<usize>,
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
            if let Some(nt) = build_seed_nt(&third_pm, patch_n, &all_juncs, match_index) {
                let entry = NtEntry::Phase1(nt);
                if state.seen.contains_key(&entry) {
                    continue;
                }
                let NtEntry::Phase1(ref nt_ref) = entry else {
                    unreachable!()
                };
                // trial.add_tile above verified the central glues into
                // the normalized two-tile patch. Also verify the
                // abstract NT reconstructs from its vt_seq (needed for
                // non-convex tiles where seed-gen and reconstruct can
                // diverge).
                if !nt_is_valid(nt_ref, match_index) {
                    continue;
                }
                let id = state.entries.len() + 1;
                state.seen.insert(entry.clone(), id);
                state.entries.push(entry);
                state.queue.push_back(id);
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
    let ccw_anchor_on_central = (cw_anchor_on_central + central_n - third_pm.len) % central_n;

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

/// BFS growth: dequeue each entry id, enumerate one-petal successors
/// (phase-1 step for `NtEntry::Phase1`, phase-2 step for open
/// `NtEntry::Phase2`), dedup against `seen`, and record transitions.
/// Phase-1 entries grow `num_ctx_tiles` by 1 per step; phase-2 entries
/// grow corona by 1 per step. The transition graph is a DAG.
///
/// Closed `SurroundedTile` entries (`is_closed = true`) are pushed
/// into `entries` but not into the BFS queue — they are terminal.
/// [`NeighborhoodIndex::classify_all`] recognises them directly as
/// Blessed escape destinations.
fn bfs_phase<T: IsComplex + IsRingOrField + Units>(
    state: &mut BfsState,
    match_index: &Arc<MatchTypeIndex<T>>,
) {
    while let Some(src_id) = state.queue.pop_front() {
        let entry = state.entries[src_id - 1].clone();
        let outcomes = match &entry {
            NtEntry::Phase1(nt) => explore_phase1(nt, match_index),
            NtEntry::Phase2(st) => explore_phase2(st, match_index),
        };
        for outcome in outcomes {
            let new_entry = match outcome.kind {
                OutcomeKind::Phase1(nt) => NtEntry::Phase1(nt),
                OutcomeKind::Phase2(st) => NtEntry::Phase2(st),
            };
            let needs_bfs = matches!(
                &new_entry,
                NtEntry::Phase1(_)
                    | NtEntry::Phase2(SurroundedTile {
                        is_closed: false,
                        ..
                    })
            );
            let dst_id = if let Some(&id) = state.seen.get(&new_entry) {
                id
            } else {
                let id = state.entries.len() + 1;
                state.seen.insert(new_entry.clone(), id);
                state.entries.push(new_entry);
                // Phase-1 entries and open phase-2 entries get further
                // BFS exploration. Closed phase-2 entries are terminal
                // (= they carry the closed-corona data; no outgoing
                // transitions, recognised as escape by classify_all).
                if needs_bfs {
                    state.queue.push_back(id);
                }
                id
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

/// The result of one BFS step: a new (or existing) catalog entry to
/// transition into. Phase-1 entries describe states where the central
/// is still on the boundary; phase-2 entries (`SurroundedTile`) describe
/// states where the central has been edge-absorbed, with a flag for
/// whether the surrounding corona has closed the final open vertex.
enum OutcomeKind {
    Phase1(NeighborhoodType),
    Phase2(SurroundedTile),
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
    /// and `[gap_start, aug_n)` are the central tile gap edges. Used
    /// by [`try_step_ccw`] to identify the CW frontier ctx edge for
    /// the phase-2 transition's `is_closed` check.
    gap_start: usize,
    aug_n: usize,
}

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

fn explore_phase1<T: IsComplex + IsRingOrField + Units>(
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

/// Phase-2 BFS step from a `SurroundedTile`. Reconstructs the corona
/// patch via [`SurroundedTile::to_patch`] (which calls
/// [`GrowingPatch::from_parts`] directly — no `vt_seq` reconstruction,
/// so the closed-corona / fully-interior-tile pitfalls of
/// `construct_witness_from_vt_sequence` are avoided). Attempts to
/// attach a petal at the open vertex (= `st.open_vertex_pos`).
///
/// For each successful glue, builds the next `SurroundedTile`.
/// `is_closed` on the result is true iff the petal's match absorbed
/// both edges incident to the old open vertex (= the petal seals the
/// angle around it).
fn explore_phase2<T: IsComplex + IsRingOrField + Units>(
    st: &SurroundedTile,
    match_index: &Arc<MatchTypeIndex<T>>,
) -> Vec<ExploreOutcome> {
    debug_assert!(
        !st.is_closed,
        "explore_phase2 should not run on closed SurroundedTiles"
    );
    let corona = match st.to_patch(Arc::clone(match_index)) {
        Some(p) => p,
        None => return Vec::new(),
    };
    let corona_n = corona.boundary_len();
    if corona_n == 0 {
        return Vec::new();
    }
    let open_vertex_pos = st.open_vertex_pos;
    let ccw_edge = open_vertex_pos;
    let cw_edge = (open_vertex_pos + corona_n - 1) % corona_n;

    let candidates = GrowingPatch::compute_candidates_covering_position(
        corona.match_index(),
        corona.angles(),
        corona.edges(),
        ccw_edge,
    );
    let mut results = Vec::new();
    for petal_pm in candidates {
        if !match_absorbs_edge(&petal_pm, ccw_edge, corona_n) {
            continue;
        }
        let mut trial = corona.clone();
        if !trial.add_tile(&petal_pm) {
            continue;
        }
        let is_closed = match_absorbs_edge(&petal_pm, cw_edge, corona_n);
        let open_vertex_trial_pos = if is_closed {
            None
        } else {
            let ccw_on_corona = (petal_pm.start_a + petal_pm.len) % corona_n;
            Some((open_vertex_pos + corona_n - ccw_on_corona) % corona_n)
        };
        let new_st = build_surrounded_tile_from_trial(
            trial,
            st.central_tile_id,
            is_closed,
            open_vertex_trial_pos,
        );
        results.push(ExploreOutcome {
            side: TransitionSide::Ccw,
            petal_pm,
            kind: OutcomeKind::Phase2(new_st),
        });
    }
    results
}

/// Build a `SurroundedTile` from a trial patch. Normalizes the patch
/// (= lex-min rotation + dense `0..k` `patch_tile_ids`) and stores
/// its boundary state directly. For open entries, maps the original
/// open-vertex position into the normalized boundary via the rotation
/// offset that `normalize` reports.
fn build_surrounded_tile_from_trial<T: IsComplex + IsRingOrField + Units>(
    mut trial: GrowingPatch<T>,
    central_tile_id: usize,
    is_closed: bool,
    open_vertex_trial_pos: Option<usize>,
) -> SurroundedTile {
    let n = trial.boundary_len();
    let rot = trial.normalize();
    let open_vertex_pos = match open_vertex_trial_pos {
        Some(p) if n > 0 => (p + n - rot) % n,
        _ => 0,
    };
    let patch_tile_ids = trial.patch_tile_ids().to_vec();
    let num_tile_instances = patch_tile_ids.iter().max().map_or(0, |m| m + 1);
    SurroundedTile {
        central_tile_id,
        is_closed,
        open_vertex_pos,
        num_tile_instances,
        angles: trial.angles().to_vec(),
        edges: trial.edges().to_vec(),
        inner_chains: trial.inner_chains().to_vec(),
        patch_tile_ids,
    }
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

    // Central edges fully absorbed → emit a phase-2 SurroundedTile.
    //
    // The CCW frontier ctx edge (aug position 0) is always in the
    // petal's match by the candidate filter, so the CCW anchor
    // (aug vertex 0) is always absorbed when the petal absorbs all
    // central. The CW anchor (aug vertex `gap_start`) survives iff
    // its outer ctx edge (= aug edge `gap_start - 1`) is not also
    // absorbed.
    //
    // - **Closed phase 2** (petal absorbs CW outer ctx too): both
    //   anchors gone; central is fully surrounded; terminal.
    // - **Open phase 2** (CW outer ctx not absorbed): CW anchor
    //   survives as the sole frontier vertex touching the central;
    //   phase-2 BFS continues there.
    if !trial.patch_tile_ids().contains(&ac.central_ptid) {
        let cw_outer_edge = (ac.gap_start + ac.aug_n - 1) % ac.aug_n;
        let is_closed = match_absorbs_edge(&petal_pm, cw_outer_edge, ac.aug_n);
        let canonical_open_trial_pos = if is_closed {
            None
        } else {
            let ccw_pos_on_aug = (petal_pm.start_a + petal_pm.len) % ac.aug_n;
            Some((ac.gap_start + ac.aug_n - ccw_pos_on_aug) % ac.aug_n)
        };
        let st = build_surrounded_tile_from_trial(
            trial,
            nt.central_tile_id,
            is_closed,
            canonical_open_trial_pos,
        );
        return Some(ExploreOutcome {
            side: TransitionSide::Ccw,
            petal_pm,
            kind: OutcomeKind::Phase2(st),
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
        kind: OutcomeKind::Phase1(new_nt),
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
        expected_closed_phase2: usize,
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
        // Closed phase-2 SurroundedTile entries (= the terminal-good
        // states that classify_all treats as escape destinations).
        let closed_phase2 = idx
            .entries()
            .iter()
            .filter(|e| matches!(e, NtEntry::Phase2(st) if st.is_closed))
            .count();

        assert_eq!(idx.num_types(), expected_types, "num_types");
        assert_eq!(idx.transitions().len(), expected_transitions, "transitions");
        assert_eq!(closed_phase2, expected_closed_phase2, "closed_phase2");
        assert_eq!(dead, expected_dead, "dead");
        assert_eq!(undead, expected_undead, "undead");
        assert_eq!(blessed, expected_blessed, "blessed");
        assert_eq!(free, expected_free, "free");
        assert_eq!(dead + undead + blessed + free, expected_types, "kinds sum");
    }

    #[test]
    #[cfg_attr(debug_assertions, ignore = "release-only: square idx build ~21s in release, much slower in debug")]
    fn square_counts() {
        // Square tiles convexly; every continuation closes. All
        // 240 256 entries are Blessed: 109 184 phase-1 + 131 072
        // phase-2 (split 65 536 closed terminals + 65 536 open
        // intermediates that phase-2 BFS eventually closes).
        // Unlike hex, square produces both open and closed phase-2
        // entries from phase-1 → phase-2 transitions.
        assert_counts(square_idx(), 240256, 698880, 65536, 0, 0, 240256, 0);
    }

    #[test]
    fn hex_counts() {
        // Hex tiles perfectly: every central-edge-absorption is a
        // closed corona (no open vertex). All 46656 phase-2 entries
        // are closed terminals, classified directly as Blessed.
        assert_counts(hex_idx(), 102600, 335664, 46656, 0, 0, 102600, 0);
    }

    #[test]
    fn spectre_counts() {
        // Spectre is non-convex; exercises every kind. Phase-1
        // entries classify the same as before phase-2 was introduced
        // (dead=6001, undead=1615, free=1914), confirming the
        // phase-1 BFS graph is structurally unchanged. Phase-2 adds
        // 445 entries (251 closed terminals + 194 open
        // intermediates), all Blessed (the closed terminals as
        // direct Blessed seeds; the open ones via propagation
        // through their successors).
        assert_counts(spectre_idx(), 10559, 11742, 251, 6001, 1615, 1029, 1914);
    }

    #[test]
    #[cfg_attr(debug_assertions, ignore = "release-only: iterates square_idx (~21s build in release)")]
    fn seeds_are_well_formed() {
        for idx in [
            square_idx() as &NeighborhoodIndex<ZZ12>,
            hex_idx(),
            spectre_idx(),
        ] {
            assert!(idx.num_types() > 0, "expected non-empty seed collection");
            let num_tiles = idx.tileset().num_tiles();
            for entry in idx.entries() {
                assert!(
                    entry.central_tile_id() < num_tiles,
                    "central_tile_id {} out of range (have {} tiles)",
                    entry.central_tile_id(),
                    num_tiles
                );
                match entry {
                    NtEntry::Phase1(nt) => assert!(
                        !nt.vt_seq.is_empty(),
                        "phase-1 seeds must have non-empty vt_seq"
                    ),
                    NtEntry::Phase2(st) => assert!(
                        !st.angles.is_empty(),
                        "phase-2 entries must have non-empty boundary"
                    ),
                }
            }
        }
    }

    /// No two catalog entries are equal by `PartialEq`. This is the
    /// dedup contract of the BFS — not interface-level
    /// canonicalization. The catalog **does** contain multiple
    /// entries per interface class by design (different growth
    /// histories produce distinct entries with the same match-side
    /// interface); see the module doc on the interface-quotient view.
    #[test]
    #[cfg_attr(debug_assertions, ignore = "release-only: iterates square_idx (~21s build in release)")]
    fn entries_have_no_partial_eq_duplicates() {
        for (name, idx) in [
            ("square", square_idx() as &NeighborhoodIndex<ZZ12>),
            ("hex", hex_idx()),
            ("spectre", spectre_idx()),
        ] {
            let mut keys: std::collections::HashSet<&NtEntry> = std::collections::HashSet::new();
            for (i, entry) in idx.entries().iter().enumerate() {
                assert!(
                    keys.insert(entry),
                    "{}: duplicate entry at id {} (vt_seq={:?})",
                    name,
                    i + 1,
                    entry.vt_seq()
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
    #[cfg_attr(debug_assertions, ignore = "release-only: square idx build ~21s in release, much slower in debug")]
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
    /// in both production (`explore_phase1` / `explore_phase2` filter)
    /// and in the brute side of `spectre_ccw_only_is_complete`, so a
    /// bug here would not be caught by the existing brute test
    /// (circular).
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
                        assert_eq!(got, want, "n={n} start={start} len={len} target={target}");
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
    #[cfg_attr(debug_assertions, ignore = "release-only: validate() reconstructs all square phase-2 patches (~28s in release)")]
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
        let mut fp_to_kind: FxHashMap<(usize, usize, Vec<i8>, usize), (NtKind, usize)> =
            FxHashMap::default();
        let mut violations = 0usize;
        let mut first_violation: Option<String> = None;
        for (i, entry) in idx.entries().iter().enumerate() {
            let Some(nt) = entry.as_phase1() else {
                continue;
            };
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
    /// Uses NO logic from `explore_phase1` / `explore_phase2` /
    /// `try_step_*` / `compute_candidates_covering_position` —
    /// pure brute force.
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
        for entry in idx.entries() {
            let fp = match entry {
                NtEntry::Phase1(nt) => aug_boundary_fingerprint(nt, &mi),
                // Phase-2 entries have no aug boundary to fingerprint
                // (central is interior). Skip; their classification
                // correctness is checked separately.
                NtEntry::Phase2(_) => None,
            };
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
            // id_to_fp[i] is None for phase-2 destinations (= closed
            // and open phase-2 entries don't have an aug-boundary
            // fingerprint), so transitions to phase-2 surface as
            // dst_fp = None.
            let dst_fp: DstFp = id_to_fp[t.dst_id - 1].clone();
            recorded_by_key
                .entry((t.src_id, t.tile_id, t.tile_offset))
                .or_default()
                .insert(dst_fp);
        }

        for (i, entry) in idx.entries().iter().enumerate() {
            let src_id = i + 1;
            let Some(nt) = entry.as_phase1() else {
                continue;
            };
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

                        let trial_dst_fp: DstFp =
                            if !trial.patch_tile_ids().contains(&ac.central_ptid) {
                                None
                            } else {
                                let trial_angles = trial.angles().to_vec();
                                let rot = crate::intgeom::rat::lex_min_rot(&trial_angles);
                                let mut canon = trial_angles;
                                canon.rotate_left(rot);
                                let fp = (nt.central_tile_id, nt.cw_anchor_on_central, canon);
                                assert!(
                                    fp_set.contains_key(&fp),
                                    "completeness violation: src={} brute glue \
                                     (tile={}, start_a={}, len={}, start_b={}) \
                                     produces state not in catalog",
                                    src_id,
                                    tile_b,
                                    ns_u,
                                    len,
                                    ne_u
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
                            src_id,
                            tile_b,
                            ns_u,
                            len,
                            ne_u,
                            if trial_dst_fp.is_none() {
                                "CLOSED".to_string()
                            } else {
                                "Open(fp)".to_string()
                            }
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
    ///   - `reaches_closed[i]` = exists a forward path from i to a
    ///     closed phase-2 terminal
    ///   - `reaches_dead[i]` = exists a forward path from i to a
    ///     stuck (no-outgoing, not closed-terminal) entry
    ///
    /// The four kinds are pinned by these reachability properties:
    ///   - Dead: no outgoing, not closed-terminal
    ///   - Undead: has outgoing AND NOT reaches_closed
    ///   - Blessed: has outgoing AND NOT reaches_dead (or IS
    ///     closed-terminal)
    ///   - Free: has outgoing AND reaches_closed AND reaches_dead
    ///
    /// Catches the V2b shape (closing-as-escape blindspot) AND the
    /// "Free is the trash bucket" shape where Free becomes vacuous.
    fn assert_classify_invariants<T: IsComplex + IsRingOrField + Units>(
        idx: &NeighborhoodIndex<T>,
    ) {
        let kinds = idx.classify_all();
        let n = idx.num_types();

        // Identify closed phase-2 entries (= terminal Blessed states).
        let is_closed_terminal: Vec<bool> = idx
            .entries()
            .iter()
            .map(|e| matches!(e, NtEntry::Phase2(st) if st.is_closed))
            .collect();

        let mut preds: Vec<Vec<usize>> = vec![Vec::new(); n];
        let mut has_outgoing = vec![false; n];
        for t in idx.transitions() {
            let src = t.src_id - 1;
            has_outgoing[src] = true;
            preds[t.dst_id - 1].push(src);
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

        // "reaches_closed": there's a forward path from i to a
        // closed-terminal entry (Blessed escape).
        let closing_seeds: Vec<usize> = (0..n).filter(|&i| is_closed_terminal[i]).collect();
        let reaches_closed = propagate(closing_seeds);

        // "reaches_dead": there's a forward path from i to an entry
        // that's stuck (= has no outgoing AND is not closed-terminal).
        let dead_seeds: Vec<usize> = (0..n)
            .filter(|&i| !has_outgoing[i] && !is_closed_terminal[i])
            .collect();
        let reaches_dead = propagate(dead_seeds);

        for (i, &kind) in kinds.iter().enumerate() {
            let id = i + 1;
            // Closed-terminal entries are always Blessed (no
            // outgoing but reachable-as-closed trivially).
            let expected = if is_closed_terminal[i] {
                NtKind::Blessed
            } else if !has_outgoing[i] {
                NtKind::Dead
            } else if !reaches_closed[i] {
                NtKind::Undead
            } else if !reaches_dead[i] {
                NtKind::Blessed
            } else {
                NtKind::Free
            };
            assert_eq!(
                kind,
                expected,
                "NT id {id} classify_all says {:?} but reachability says {:?} \
                 (has_outgoing={} reaches_closed={} reaches_dead={} closed_terminal={})",
                kind,
                expected,
                has_outgoing[i],
                reaches_closed[i],
                reaches_dead[i],
                is_closed_terminal[i]
            );
        }
    }

    #[test]
    #[cfg_attr(debug_assertions, ignore = "release-only: square idx build ~21s in release, much slower in debug")]
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
        // Both src_id and dst_id are 1-based entry ids.
        for t in sq_idx.transitions() {
            assert!(t.src_id >= 1 && t.src_id <= sq_idx.num_types());
            assert!(t.dst_id >= 1 && t.dst_id <= sq_idx.num_types());
        }

        let hex_idx = hex_idx();
        assert!(hex_idx.num_types() > 0);
        assert!(!hex_idx.transitions().is_empty());
        for t in hex_idx.transitions() {
            assert!(t.src_id >= 1 && t.src_id <= hex_idx.num_types());
            assert!(t.dst_id >= 1 && t.dst_id <= hex_idx.num_types());
        }
    }

    fn validate_seeds<T: IsComplex + IsRingOrField + Units>(
        idx: &NeighborhoodIndex<T>,
    ) -> Vec<String> {
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        let mut errors = Vec::new();
        for entry in idx.entries() {
            let Some(nt) = entry.as_phase1() else {
                continue;
            };
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
    #[cfg_attr(debug_assertions, ignore = "release-only: square idx build ~21s in release, much slower in debug")]
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
    fn assert_seed_geometry<T: IsComplex + IsRingOrField + Units>(idx: &NeighborhoodIndex<T>) {
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        for (i, entry) in idx.entries().iter().enumerate() {
            let Some(nt) = entry.as_phase1() else {
                continue;
            };
            let ac = build_attached_context(nt, &mi).unwrap_or_else(|| {
                panic!("seed {}: build_attached_context failed for nt {:?}", i, nt)
            });
            let gap_len = ac.aug_n - ac.gap_start;
            assert!(gap_len > 0, "seed {}: gap should be non-empty", i);
            assert!(
                ac.gap_start < ac.aug_n,
                "seed {}: gap_start out of range",
                i
            );
            let central_len = idx.tileset().rat(nt.central_tile_id).seq().len();
            assert!(
                gap_len < central_len,
                "seed {}: gap_len {} should be < central_len {}",
                i,
                gap_len,
                central_len
            );
            assert!(
                !explore_phase1(nt, &mi).is_empty(),
                "seed {}: should have at least one successful petal",
                i
            );
        }
    }

    #[test]
    #[cfg_attr(debug_assertions, ignore = "release-only: square idx build ~21s in release, much slower in debug")]
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
                        let patch_rat =
                            crate::intgeom::rat::Rat::from_slice_unchecked(patch.angles());

                        // Brute third-tile candidates over the 2-tile boundary.
                        let mut seen_third: FxHashMap<(usize, usize, usize, usize), ()> =
                            FxHashMap::default();
                        for c in 0..num_tiles {
                            let tile_c = tileset.rat(c);
                            let n_c = tile_c.seq().len();
                            for ic in 0..patch_n {
                                for jc in 0..n_c {
                                    let (ns3, len3, ne3) =
                                        patch_rat.get_match((ic as i64, jc as i64), tile_c);
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
                                    let cw_anchor_in_canon = (gap_start + aug_n - rot) % aug_n;
                                    brute.insert((c, ne3_u, canon, cw_anchor_in_canon), ());
                                }
                            }
                        }
                    }
                }
            }
        }

        // Catalog side: canonical class of every entry with
        // num_ctx_tiles == 2 (= seed entries). Phase-2 entries have no
        // num_ctx_tiles; skip them.
        let mut catalog: FxHashMap<Class, ()> = FxHashMap::default();
        for entry in idx.entries() {
            let Some(nt) = entry.as_phase1() else {
                continue;
            };
            if nt.num_ctx_tiles != 2 {
                continue;
            }
            if let Some(fp) = canonical_class_fingerprint(nt, &mi) {
                catalog.insert(fp, ());
            }
        }

        let brute_only: Vec<&Class> = brute.keys().filter(|k| !catalog.contains_key(*k)).collect();
        let catalog_only: Vec<&Class> =
            catalog.keys().filter(|k| !brute.contains_key(*k)).collect();
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
    #[cfg_attr(debug_assertions, ignore = "release-only: square idx build ~21s in release, much slower in debug")]
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
        for (i, entry) in idx.entries().iter().enumerate() {
            let Some(nt) = entry.as_phase1() else {
                continue;
            };
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
        for (i, entry) in idx.entries().iter().enumerate() {
            let Some(nt) = entry.as_phase1() else {
                continue;
            };
            let ac = build_attached_context(nt, &mi)
                .unwrap_or_else(|| panic!("seed {}: build_attached_context failed", i));
            assert_eq!(
                ac.frontier_pos_on_aug, 0,
                "seed {}: gap_end (pos 0 on aug) should be the CCW frontier junction",
                i
            );
            assert!(
                ac.aug.is_junction(0),
                "seed {}: position 0 must be junction",
                i
            );
        }
    }

    #[test]
    fn explore_one_square_has_petals() {
        let idx = square_idx();
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        let nt = idx.entries()[0]
            .as_phase1()
            .expect("first entry should be phase-1");
        assert!(
            !explore_phase1(nt, &mi).is_empty(),
            "should have at least one petal outcome"
        );
    }

    #[test]
    fn explore_one_hex_seeds_reconstruct() {
        let idx = hex_idx();
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        let mut checked = 0usize;
        let mut errors = Vec::new();
        for (i, entry) in idx.entries().iter().enumerate() {
            let Some(nt) = entry.as_phase1() else {
                continue;
            };
            for outcome in explore_phase1(nt, &mi) {
                let new_nt = match &outcome.kind {
                    OutcomeKind::Phase2(_) => continue,
                    OutcomeKind::Phase1(nt) => nt.clone(),
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
        for (i, entry) in idx.entries().iter().enumerate() {
            let Some(nt) = entry.as_phase1() else {
                continue;
            };
            for outcome in explore_phase1(nt, &mi) {
                let new_nt = match &outcome.kind {
                    OutcomeKind::Phase2(_) => continue,
                    OutcomeKind::Phase1(nt) => nt.clone(),
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
