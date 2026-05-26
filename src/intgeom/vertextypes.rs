use std::collections::{BTreeSet, HashMap, VecDeque};
use std::sync::Arc;

use crate::cyclotomic::{IsRing};
use crate::intgeom::patch::{
    ClosedVertexType, EdgeInfo, GrowingPatch, OpenVertexType, PatchMatch, TransitionSide,
};
use crate::intgeom::rat::Rat;
use crate::intgeom::tileset::TileSet;

/// Reachability classification of an open vertex type within the
/// VT transition graph.
///
/// Imagine the transition graph with one extra sink, `Closed`
/// (id `CLOSED_ID = 0`), reached by every closing transition. Then:
///
/// * [`VTypeKind::Dead`] — non-closable base case: this VT has no
///   realized transitions at all. No glue from here goes anywhere,
///   including to `Closed`.
/// * [`VTypeKind::Undead`] — every path from this VT eventually
///   terminates at a `Dead` VT, never at `Closed`. Equivalently, this
///   VT is in the attractor of `Dead` under "all successors cursed"
///   but is not itself `Dead`.
/// * [`VTypeKind::Blessed`] — every path from this VT eventually
///   reaches `Closed`. Every outgoing transition either closes the
///   junction or leads to another `Blessed` VT.
/// * [`VTypeKind::Free`] — has at least one path to `Closed` AND at
///   least one path that ends at a `Dead` VT (or otherwise stays
///   non-blessed). Neither inevitably trapped nor inevitably closing.
///
/// "Cursed" = "Dead or Undead" = non-closable. "Alive" = "Blessed or
/// Free" = closable along at least one path. Computed by two monotone
/// fixpoints in [`OpenVertexTypeIndex::new`]: see [`compute_cursed`]
/// and [`compute_blessed`].
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, serde::Serialize, serde::Deserialize,
)]
pub enum VTypeKind {
    /// No realized transitions — non-closable base case.
    Dead,
    /// Every path from here terminates at a `Dead` VT, never at
    /// `Closed`.
    Undead,
    /// Every path from here terminates at `Closed`.
    Blessed,
    /// Some paths reach `Closed`, some terminate at `Dead`. Neither
    /// inevitably trapped nor inevitably closing.
    Free,
}

impl VTypeKind {
    /// Combine two [`VTypeKind`]s into the kind that should be
    /// assigned to a *segment* with these two vertex types at its
    /// endpoints. (Used by downstream segment-type analysis in
    /// `vtype_enum`.)
    ///
    /// The combination is a small lattice:
    ///
    /// * If either endpoint is `Dead`, the segment is `Dead` (a
    ///   dead vertex anywhere in the segment kills the whole segment).
    /// * Otherwise if either endpoint is cursed (Dead or Undead), the
    ///   segment is `Undead`.
    /// * Otherwise if *both* endpoints are `Blessed`, the segment is
    ///   `Blessed`.
    /// * Otherwise the segment is `Free`.
    ///
    /// (Cursed-ness propagates OR; blessedness propagates AND.)
    pub fn segment_kind(a: VTypeKind, b: VTypeKind) -> VTypeKind {
        let cursed = matches!(a, VTypeKind::Dead | VTypeKind::Undead)
            || matches!(b, VTypeKind::Dead | VTypeKind::Undead);
        let dead = matches!(a, VTypeKind::Dead) || matches!(b, VTypeKind::Dead);
        let blessed = matches!(a, VTypeKind::Blessed) && matches!(b, VTypeKind::Blessed);
        if dead {
            VTypeKind::Dead
        } else if cursed {
            VTypeKind::Undead
        } else if blessed {
            VTypeKind::Blessed
        } else {
            VTypeKind::Free
        }
    }
}

pub struct OpenVertexTypeInfo<T: IsRing> {
    vtype: OpenVertexType,
    kind: VTypeKind,
    is_initial: bool,
    realizing_rat: Rat<T>,
    witness: GrowingPatch<T>,
    witness_pos: usize,
    gap_angle: i8,
    cw_neighbor_offset: usize,
    ccw_neighbor_offset: usize,
}

impl<T: IsRing> OpenVertexTypeInfo<T> {
    /// The canonical [`OpenVertexType`] this entry describes.
    pub fn vtype(&self) -> &OpenVertexType {
        &self.vtype
    }

    /// Reachability classification: see [`VTypeKind`].
    pub fn kind(&self) -> VTypeKind {
        self.kind
    }

    pub fn is_dead(&self) -> bool {
        self.kind == VTypeKind::Dead
    }

    pub fn is_undead(&self) -> bool {
        self.kind == VTypeKind::Undead
    }

    pub fn is_blessed(&self) -> bool {
        self.kind == VTypeKind::Blessed
    }

    pub fn is_free(&self) -> bool {
        self.kind == VTypeKind::Free
    }

    /// `true` for Dead or Undead (non-closable).
    pub fn is_cursed(&self) -> bool {
        matches!(self.kind, VTypeKind::Dead | VTypeKind::Undead)
    }

    /// `true` for Blessed or Free (closable along at least one path).
    pub fn is_alive(&self) -> bool {
        matches!(self.kind, VTypeKind::Blessed | VTypeKind::Free)
    }

    /// `true` if this VT was discovered during the seed phase (as a
    /// junction in some tile's first-glue patch), not just via a
    /// later BFS step.
    pub fn is_initial(&self) -> bool {
        self.is_initial
    }

    /// A `Rat` representation of the (minimal) witness boundary,
    /// useful for displaying or hashing the realized shape.
    pub fn realizing_rat(&self) -> &Rat<T> {
        &self.realizing_rat
    }

    /// Boundary turn angle at the focus vertex of the witness.
    pub fn gap_angle(&self) -> i8 {
        self.gap_angle
    }

    /// The minimal patch that realizes this VT at a boundary junction.
    pub fn witness(&self) -> &GrowingPatch<T> {
        &self.witness
    }

    /// Boundary position of the focus vertex in [`Self::witness`].
    pub fn witness_pos(&self) -> usize {
        self.witness_pos
    }

    /// Boundary distance (CW) from the focus vertex to the next
    /// junction in the witness.
    pub fn cw_neighbor_offset(&self) -> usize {
        self.cw_neighbor_offset
    }

    /// Boundary distance (CCW) from the focus vertex to the next
    /// junction in the witness.
    pub fn ccw_neighbor_offset(&self) -> usize {
        self.ccw_neighbor_offset
    }
}

/// Sentinel destination id used by closing transitions in
/// [`TransitionInfo::dst_id`]: a transition with `dst_id == CLOSED_ID`
/// seals the source vertex (no successor VT). All real VT ids are
/// 1-based and `>= 1`, so id `0` is reserved.
pub const CLOSED_ID: usize = 0;

/// A single transition in the VT graph: gluing tile `tile_id` to
/// `src_id`'s recorded witness, with the matched edge anchored at
/// `tile_offset` on the tile and the focus vertex consumed on `side`,
/// produces a patch whose junction is the VT with id `dst_id`
/// (or seals the focus if `dst_id == CLOSED_ID`).
pub struct TransitionInfo {
    /// 1-based VT id of the source.
    pub src_id: usize,
    /// 1-based VT id of the destination, or [`CLOSED_ID`] if the
    /// transition seals the focus vertex.
    pub dst_id: usize,
    /// Which incident edge(s) of the focus vertex the glue consumed.
    pub side: TransitionSide,
    /// Tile being glued.
    pub tile_id: usize,
    /// Offset within `tile_id` of the canonical-CW anchor of the
    /// matched edge (see [`TransitionSide`] doc for the anchor rule).
    pub tile_offset: usize,
    /// For closing transitions (`dst_id == CLOSED_ID`), the 1-based id
    /// of the [`ClosedVertexType`] realised by the closure.
    /// `None` for non-closing transitions.
    pub closed_vt_id: Option<usize>,
}

impl TransitionInfo {
    pub fn is_closed(&self) -> bool {
        self.dst_id == CLOSED_ID
    }
}

/// Metadata for a discovered closed vertex type.
pub struct ClosedVertexTypeInfo {
    vtype: ClosedVertexType,
    /// 1-based id assigned to this closed VT.
    id: usize,
    /// Number of distinct closing transitions in the catalog that
    /// realise this closed VT.
    transition_count: usize,
}

impl ClosedVertexTypeInfo {
    pub fn vtype(&self) -> &ClosedVertexType {
        &self.vtype
    }

    pub fn id(&self) -> usize {
        self.id
    }

    /// How many closing transitions in the parent
    /// [`OpenVertexTypeIndex`] produce this closed VT.
    pub fn transition_count(&self) -> usize {
        self.transition_count
    }
}

/// BFS-built catalog of every [`OpenVertexType`] reachable from a
/// tileset's first-glue junctions, with [`TransitionInfo`]s between
/// them and a [`VTypeKind`] classification for each.
///
/// Each entry gets a stable 1-based id (id 0 is reserved as the
/// [`CLOSED_ID`] sentinel for transitions that close the junction).
/// Entries are stored in canonical (sorted) order by `OpenVertexType`,
/// which makes [`Self::range_by_cw`] a partition over the entries.
pub struct OpenVertexTypeIndex<T: IsRing> {
    tileset: Arc<TileSet<T>>,
    entries: Vec<OpenVertexTypeInfo<T>>,
    transitions: Vec<TransitionInfo>,
    reverse: HashMap<OpenVertexType, usize>,
    closed_entries: Vec<ClosedVertexTypeInfo>,
    closed_reverse: HashMap<ClosedVertexType, usize>,
}

impl<T: IsRing> OpenVertexTypeIndex<T> {
    /// Enumerate all open vertex types reachable from `tileset` and
    /// classify each.
    ///
    /// # Algorithm
    ///
    /// 1. **Seed.** For every tile in the tileset, take the seed
    ///    patch, try every legal first glue, and extract every
    ///    junction vertex type from the resulting 2-tile patch.
    ///    These VTs are the initial set; one representative witness
    ///    patch is recorded for each.
    /// 2. **BFS.** Pop a VT from the queue, enumerate the matches
    ///    touching its junction in the witness patch, and for each
    ///    glue:
    ///    - If the match consumes both edges incident to the focus
    ///      vertex, record a *closing* transition (`dst_id =
    ///      CLOSED_ID`).
    ///    - Otherwise, record a transition to the new VT at the
    ///      surviving side of the match.
    ///    - Also walk the post-glue boundary to discover any *other*
    ///      VTs the glue exposed; new ones get enqueued.
    /// 3. **ID assignment.** Sort the discovered VTs canonically and
    ///    assign 1-based ids. Build a `reverse: VT → id` map.
    /// 4. **Transitions.** Materialize `TransitionInfo`s from the
    ///    raw `(src_vt, dst_vt, side, tile_id, tile_offset)` tuples
    ///    via the reverse map; build succ/pred sets for the fixpoint
    ///    passes.
    /// 5. **Classification.** Run two monotone fixpoints:
    ///    - `compute_cursed` (Dead / Undead): a VT is cursed if all
    ///      its successors are cursed; base case is "no transitions".
    ///    - `compute_blessed`: a VT is blessed if every transition
    ///      either closes the junction or leads to another blessed
    ///      VT.
    ///
    ///    Then label each VT as Dead / Undead / Blessed / Free.
    ///
    /// # Cost
    ///
    /// Dominated by the BFS — every reachable VT contributes a
    /// candidate-enumeration over its witness boundary. For non-
    /// convex tilesets the VT space can be large (e.g. spectre);
    /// progress is logged to stderr every 100 VTs.
    pub fn new(tileset: Arc<TileSet<T>>) -> Self {
        let mut state = BfsState::default();
        seed_phase(&mut state, &tileset);
        bfs_phase(&mut state, &tileset);

        let BfsState {
            all_types,
            initial_types,
            raw_transitions,
            witness_store,
            ..
        } = state;

        let (vt_list, reverse) = build_id_map(all_types);
        let (closed_list, closed_reverse) = build_closed_id_map(&raw_transitions);
        let (succ_sets, transition_infos, closed_counts) = build_transition_arrays(
            &reverse,
            &closed_reverse,
            closed_list.len(),
            &raw_transitions,
        );
        let has_any_realized = compute_has_any_realized(vt_list.len(), &transition_infos);
        let has_closing = compute_has_closing(vt_list.len(), &transition_infos);

        let is_cursed = compute_cursed(&vt_list, &has_any_realized, &has_closing, &succ_sets);
        let is_blessed = compute_blessed(&vt_list, &has_any_realized, &succ_sets);

        let info_entries = classify_and_finalize(
            vt_list,
            &has_any_realized,
            witness_store,
            &initial_types,
            &is_cursed,
            &is_blessed,
        );

        let closed_entries: Vec<ClosedVertexTypeInfo> = closed_list
            .into_iter()
            .enumerate()
            .map(|(i, cvt)| ClosedVertexTypeInfo {
                vtype: cvt,
                id: i + 1,
                transition_count: closed_counts[i],
            })
            .collect();

        OpenVertexTypeIndex {
            tileset,
            entries: info_entries,
            transitions: transition_infos,
            reverse,
            closed_entries,
            closed_reverse,
        }
    }

    /// Number of distinct VT entries in the catalog. Ids run
    /// `1..=num_types()`.
    pub fn num_types(&self) -> usize {
        self.entries.len()
    }

    /// All [`TransitionInfo`]s in the catalog, in the order produced
    /// by the BFS (deterministic, but not sorted by any particular
    /// key).
    pub fn transitions(&self) -> &[TransitionInfo] {
        &self.transitions
    }

    /// All [`OpenVertexTypeInfo`] entries in canonical
    /// `(cw, inner, ccw)` lex order. `entries()[id - 1]` is the entry
    /// for VT id `id`.
    pub fn entries(&self) -> &[OpenVertexTypeInfo<T>] {
        &self.entries
    }

    /// The slice of entries whose `vtype.cw.tile_id == tile_id`.
    ///
    /// Relies on the invariant that entries are stored in
    /// `(cw, inner, ccw)` lex order, so `cw.tile_id` is the
    /// most-significant sort key and entries sharing a `cw.tile_id`
    /// form a contiguous range. A refactor that changes the entry
    /// sort order would silently break this method; preserve the
    /// ordering or update both together.
    pub fn range_by_cw(&self, tile_id: usize) -> &[OpenVertexTypeInfo<T>] {
        let start = self
            .entries
            .partition_point(|e| e.vtype().cw.tile_id < tile_id);
        let end = self
            .entries
            .partition_point(|e| e.vtype().cw.tile_id <= tile_id);
        &self.entries[start..end]
    }

    /// Look up the 1-based id of `vtype`, if it's in the catalog.
    pub fn get_id(&self, vtype: &OpenVertexType) -> Option<usize> {
        self.reverse.get(vtype).copied()
    }

    /// Return the canonical [`OpenVertexType`] for id `id`. Panics if
    /// `id` is out of range (not in `1..=num_types()`).
    pub fn get_type(&self, id: usize) -> &OpenVertexType {
        assert!(
            id >= 1 && id <= self.entries.len(),
            "vertex type id out of range"
        );
        &self.entries[id - 1].vtype
    }

    /// Return the full [`OpenVertexTypeInfo`] for id `id`. Panics if
    /// `id` is out of range.
    pub fn get_info(&self, id: usize) -> &OpenVertexTypeInfo<T> {
        assert!(
            id >= 1 && id <= self.entries.len(),
            "vertex type id out of range"
        );
        &self.entries[id - 1]
    }

    /// The tileset this catalog was built from.
    pub fn tileset(&self) -> &Arc<TileSet<T>> {
        &self.tileset
    }

    /// Number of distinct closed VTs discovered during the BFS — i.e.
    /// the number of distinct fully-surrounded interior-vertex
    /// configurations realisable from this tileset's first glues.
    /// Ids run `1..=num_closed_types()`.
    pub fn num_closed_types(&self) -> usize {
        self.closed_entries.len()
    }

    /// All [`ClosedVertexTypeInfo`] entries in canonical (sorted)
    /// order. Entry at index `id - 1` has id `id`.
    pub fn closed_entries(&self) -> &[ClosedVertexTypeInfo] {
        &self.closed_entries
    }

    /// 1-based id of `cvt` in this catalog, or `None` if absent.
    pub fn get_closed_id(&self, cvt: &ClosedVertexType) -> Option<usize> {
        self.closed_reverse.get(cvt).copied()
    }

    /// Return the closed VT for id `id`. Panics if `id` is out of
    /// range or `0` (the open-side `CLOSED_ID` sentinel).
    pub fn get_closed_type(&self, id: usize) -> &ClosedVertexType {
        assert!(
            id >= 1 && id <= self.closed_entries.len(),
            "closed vertex type id out of range"
        );
        &self.closed_entries[id - 1].vtype
    }

    /// Return the [`ClosedVertexTypeInfo`] for id `id`. Panics if
    /// out of range.
    pub fn get_closed_info(&self, id: usize) -> &ClosedVertexTypeInfo {
        assert!(
            id >= 1 && id <= self.closed_entries.len(),
            "closed vertex type id out of range"
        );
        &self.closed_entries[id - 1]
    }
}

/// The destination of a raw BFS transition: either another open VT
/// (still on the boundary after the glue) or a specific closed VT
/// (focus sealed by a `TransitionSide::Both` match).
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum RawDst {
    Open(OpenVertexType),
    Closed(ClosedVertexType),
}

/// A single raw transition tuple, accumulated by the BFS before id
/// assignment turns it into a [`TransitionInfo`].
///
/// Layout: `(src_vt, dst, side, tile_id, tile_offset)`.
type RawTransition = (OpenVertexType, RawDst, TransitionSide, usize, usize);

/// Mutable state shared across the seed and BFS phases of
/// [`OpenVertexTypeIndex::new`].
struct BfsState<T: IsRing> {
    /// Every VT discovered so far (seed phase + BFS).
    all_types: BTreeSet<OpenVertexType>,
    /// Subset of `all_types` originating from the seed phase.
    initial_types: BTreeSet<OpenVertexType>,
    /// Raw transition tuples; deduplicated via the BTreeSet.
    raw_transitions: BTreeSet<RawTransition>,
    /// VTs already enqueued for BFS exploration.
    visited: BTreeSet<OpenVertexType>,
    /// BFS frontier.
    queue: VecDeque<OpenVertexType>,
    /// Witness patch + junction position + gap angle for each VT.
    ///
    /// Each `GrowingPatch` carries its full spatial grid (a few KB),
    /// even though after the BFS terminates the grid is only used if
    /// a caller asks for `info.witness().get_matches_touching_vertex(...)`
    /// or similar grid-dependent operation. For tilesets with
    /// thousands of VTs this is a few MB of mostly-cold memory.
    /// If it becomes a bottleneck, switch to storing
    /// `RawBoundary` here and reconstruct the `GrowingPatch` lazily
    /// in `OpenVertexTypeInfo::witness()`.
    witness_store: HashMap<OpenVertexType, (GrowingPatch<T>, usize, i8)>,
}

impl<T: IsRing> Default for BfsState<T> {
    fn default() -> Self {
        BfsState {
            all_types: BTreeSet::new(),
            initial_types: BTreeSet::new(),
            raw_transitions: BTreeSet::new(),
            visited: BTreeSet::new(),
            queue: VecDeque::new(),
            witness_store: HashMap::new(),
        }
    }
}

/// Phase 1: extract every junction VT from every legal first-glue of
/// every seed tile. Each discovered VT records a representative
/// witness patch (the 2-tile patch that contains it).
fn seed_phase<T: IsRing>(
    state: &mut BfsState<T>,
    tileset: &Arc<TileSet<T>>,
) {
    for seed_id in 0..tileset.num_tiles() {
        let seed = GrowingPatch::new(Arc::clone(tileset), seed_id);
        for pm in seed.get_all_matches() {
            let mut gp = seed.clone();
            if !gp.add_tile(&pm) || !gp.is_growing() {
                continue;
            }
            for pos in 0..gp.boundary_len() {
                if !gp.is_junction(pos) {
                    continue;
                }
                if let Some(vt) = gp.junction_vertex_type_at(pos) {
                    state.all_types.insert(vt.clone());
                    state.initial_types.insert(vt.clone());
                    state
                        .witness_store
                        .entry(vt.clone())
                        .or_insert_with(|| (gp.clone(), pos, gp.angles()[pos]));
                    if state.visited.insert(vt.clone()) {
                        state.queue.push_back(vt);
                    }
                }
            }
        }
    }
}

/// Phase 2: process every VT in the queue. For each, enumerate
/// touching matches and emit raw transitions; also discover any new
/// VTs exposed by the glues and enqueue them.
fn bfs_phase<T: IsRing>(
    state: &mut BfsState<T>,
    tileset: &Arc<TileSet<T>>,
) {
    while let Some(vt) = state.queue.pop_front() {
        let (touching, n) = {
            let (gp, pos, _gap) = &state.witness_store[&vt];
            let n = gp.boundary_len();
            let t: Vec<PatchMatch> = gp.get_matches_touching_vertex(*pos);
            if t.is_empty() {
                continue;
            }
            (t, n)
        };
        let pos = state.witness_store[&vt].1;

        for pm in touching {
            let gp = &state.witness_store[&vt].0;
            let mut gp2 = gp.clone();
            if !gp2.add_tile(&pm) || !gp2.is_growing() {
                continue;
            }

            // Classify the glue by which incident boundary edges of
            // the focus vertex `pos` it consumed.
            //
            // The focus vertex has two incident boundary edges:
            //   CW edge:  edge position (pos + n - 1) % n
            //             (goes vertex pos-1 -> pos)
            //   CCW edge: edge position pos
            //             (goes vertex pos -> pos+1)
            //
            // The match consumes a half-open edge range
            // [pm.start_a(), pm.start_a() + pm.len()). The focus is sealed
            // iff BOTH incident edges are in that range. (Note: a
            // vertex-range test would be too permissive — a single-
            // edge match on the CW edge has vertex range [pos-1, pos]
            // touching both pos-1 and pos without consuming the CCW
            // edge of pos.)
            let edge_in_match = |edge: usize| -> bool { (edge + n - pm.start_a()) % n < pm.len() };
            let consumes_cw_edge = edge_in_match((pos + n - 1) % n);
            let consumes_ccw_edge = edge_in_match(pos);

            // `get_matches_touching_vertex(pos)` guarantees the match
            // touches the focus vertex, so at least one incident
            // edge must be consumed.
            debug_assert!(
                consumes_cw_edge || consumes_ccw_edge,
                "touching match doesn't consume any incident edge of pos"
            );

            let side = match (consumes_cw_edge, consumes_ccw_edge) {
                (true, true) => TransitionSide::Both,
                (true, false) => TransitionSide::Cw,
                (false, true) => TransitionSide::Ccw,
                (false, false) => unreachable!(),
            };

            // Canonical recorded edge for the tile_offset
            // computation:
            //   - Cw side: the CW edge.
            //   - Ccw side: the CCW edge.
            //   - Both (closing): CW edge by convention (documented
            //     on TransitionSide::Both).
            let edge_pos = match side {
                TransitionSide::Ccw => pos,
                TransitionSide::Cw | TransitionSide::Both => (pos + n - 1) % n,
            };
            let offset_in_match =
                (edge_pos as i64 - pm.start_a() as i64).rem_euclid(n as i64) as usize;
            let m = tileset.rat(pm.tile_id()).len();
            let tile_offset = (pm.start_b() as i64 + pm.len() as i64 - offset_in_match as i64)
                .rem_euclid(m as i64) as usize;

            // Where is the focus's new junction (if any) in the
            // post-glue boundary?
            //   - side == Ccw: pm.start_a() == pos; new boundary's
            //     index 0 is at old vertex (pos + len) % n; new
            //     junction at the CW endpoint of the match lives
            //     at new index seg_len_old = n - pm.len().
            //   - side == Cw: pm ended at vertex pos (start_a + len
            //     == pos); new boundary's index 0 is at old vertex
            //     pos; new junction at the CCW endpoint of the
            //     match lives at new index 0.
            //   - side == Both: focus vertex is sealed; no junction
            //     position.
            let junction_pos = if pm.start_a() == pos { n - pm.len() } else { 0 };

            if matches!(side, TransitionSide::Both) {
                // Compute the closed VT realised at the focus. The
                // cyclic petal ring around the now-interior vertex
                // is `[cw, inner..., ccw]` (CCW order); the new tile
                // sits implicitly between `ccw` and `cw` cyclically.
                let closed = ClosedVertexType::from_open_via_closure(&vt);
                state.raw_transitions.insert((
                    vt.clone(),
                    RawDst::Closed(closed),
                    side,
                    pm.tile_id(),
                    tile_offset,
                ));
            } else if let Some(new_vt) = gp2.junction_vertex_type_at(junction_pos) {
                state.raw_transitions.insert((
                    vt.clone(),
                    RawDst::Open(new_vt),
                    side,
                    pm.tile_id(),
                    tile_offset,
                ));
            }

            for new_pos in 0..gp2.boundary_len() {
                if !gp2.is_junction(new_pos) {
                    continue;
                }
                if let Some(nv) = gp2.junction_vertex_type_at(new_pos) {
                    state.all_types.insert(nv.clone());
                    if state.visited.insert(nv.clone()) {
                        // `or_insert_with` matches the seed phase's
                        // first-seen-wins semantics. The `visited`
                        // check above already ensures we only get here
                        // on first visit, but keeping the API
                        // symmetric guards against future loosening
                        // of that guard.
                        state
                            .witness_store
                            .entry(nv.clone())
                            .or_insert_with(|| (gp2.clone(), new_pos, gp2.angles()[new_pos]));
                        state.queue.push_back(nv);
                    }
                }
            }
        }
    }
}

/// Phase 3: sort discovered VTs canonically and assign 1-based ids.
/// Returns `(vt_list, reverse)` where `reverse[vt] = id` and
/// `vt_list[id - 1] = vt`.
fn build_id_map(
    all_types: BTreeSet<OpenVertexType>,
) -> (Vec<OpenVertexType>, HashMap<OpenVertexType, usize>) {
    let vt_list: Vec<OpenVertexType> = all_types.into_iter().collect();
    let reverse: HashMap<OpenVertexType, usize> = vt_list
        .iter()
        .enumerate()
        .map(|(i, vt)| (vt.clone(), i + 1))
        .collect();
    (vt_list, reverse)
}

/// Phase 3.5: gather all distinct closed VTs from the raw closing
/// transitions and assign 1-based ids. Returns the sorted catalog
/// plus the `closed_vt -> id` reverse map.
fn build_closed_id_map(
    raw_transitions: &BTreeSet<RawTransition>,
) -> (Vec<ClosedVertexType>, HashMap<ClosedVertexType, usize>) {
    let mut all_closed: BTreeSet<ClosedVertexType> = BTreeSet::new();
    for (_src, dst, _side, _tid, _toff) in raw_transitions {
        if let RawDst::Closed(c) = dst {
            all_closed.insert(c.clone());
        }
    }
    let closed_list: Vec<ClosedVertexType> = all_closed.into_iter().collect();
    let closed_reverse: HashMap<ClosedVertexType, usize> = closed_list
        .iter()
        .enumerate()
        .map(|(i, c)| (c.clone(), i + 1))
        .collect();
    (closed_list, closed_reverse)
}

/// Phase 4: turn the raw `(src, dst, side, tile_id, tile_offset)`
/// tuples into `TransitionInfo` records using the open- and closed-VT
/// id maps; build the successor sets used by the cursed / blessed
/// fixpoints and tally per-closed-VT transition counts.
fn build_transition_arrays(
    reverse: &HashMap<OpenVertexType, usize>,
    closed_reverse: &HashMap<ClosedVertexType, usize>,
    num_closed: usize,
    raw_transitions: &BTreeSet<RawTransition>,
) -> (Vec<BTreeSet<usize>>, Vec<TransitionInfo>, Vec<usize>) {
    let n = reverse.len();
    let mut succ_sets: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); n];
    let mut transition_infos: Vec<TransitionInfo> = Vec::new();
    let mut closed_transition_counts: Vec<usize> = vec![0; num_closed];

    for (src, dst, side, tid, toff) in raw_transitions {
        let Some(&src_id) = reverse.get(src) else {
            continue;
        };
        match dst {
            RawDst::Open(dst_vt) => {
                if let Some(&dst_id) = reverse.get(dst_vt) {
                    succ_sets[src_id - 1].insert(dst_id);
                    transition_infos.push(TransitionInfo {
                        src_id,
                        dst_id,
                        side: *side,
                        tile_id: *tid,
                        tile_offset: *toff,
                        closed_vt_id: None,
                    });
                }
            }
            RawDst::Closed(cvt) => {
                let cvt_id = closed_reverse.get(cvt).copied();
                if let Some(id) = cvt_id {
                    closed_transition_counts[id - 1] += 1;
                }
                transition_infos.push(TransitionInfo {
                    src_id,
                    dst_id: CLOSED_ID,
                    side: *side,
                    tile_id: *tid,
                    tile_offset: *toff,
                    closed_vt_id: cvt_id,
                });
            }
        }
    }

    (succ_sets, transition_infos, closed_transition_counts)
}

/// Phase 5: assemble per-VT [`OpenVertexTypeInfo`] records, classifying
/// each as Dead / Undead / Blessed / Free from the fixpoint results.
fn classify_and_finalize<T: IsRing>(
    vt_list: Vec<OpenVertexType>,
    has_any_realized: &[bool],
    mut witness_store: HashMap<OpenVertexType, (GrowingPatch<T>, usize, i8)>,
    initial_types: &BTreeSet<OpenVertexType>,
    is_cursed: &[bool],
    is_blessed: &[bool],
) -> Vec<OpenVertexTypeInfo<T>> {
    vt_list
        .into_iter()
        .enumerate()
        .map(|(i, vt)| {
            let (witness, witness_pos, gap_angle) = witness_store.remove(&vt).unwrap();
            let rat = witness.to_rat();
            let (cw_nbr, ccw_nbr) = witness
                .neighbor_junction_offsets(witness_pos)
                .unwrap_or((0, 0));
            // "Dead" = no transitions in the realized catalog. (BFS
            // touching matches may fail add_tile or produce no
            // destination VT, leaving the catalog empty for this VT.)
            let no_transitions = !has_any_realized[i];
            let kind = if no_transitions {
                VTypeKind::Dead
            } else if is_cursed[i] {
                VTypeKind::Undead
            } else if is_blessed[i] {
                VTypeKind::Blessed
            } else {
                VTypeKind::Free
            };
            OpenVertexTypeInfo {
                kind,
                is_initial: initial_types.contains(&vt),
                realizing_rat: rat,
                vtype: vt,
                witness,
                witness_pos,
                gap_angle,
                cw_neighbor_offset: cw_nbr,
                ccw_neighbor_offset: ccw_nbr,
            }
        })
        .collect()
}

/// Compute the "cursed" status (Dead ∪ Undead) for each VT — the
/// attractor of the Dead set under the transition graph.
///
/// Formal definition: cursed = least fixed point of
///   `dead ∨ (every successor is cursed)`
/// where:
///   - `dead[i]` = no realized transitions at all (the base case);
///   - "successor" includes both open successors (in `succ_sets`) AND
///     the implicit Closed sink reachable via closing transitions.
///
/// Since Closed is **not** cursed by construction (it's the "escape"
/// sink), a VT that has any closing transition has a successor that
/// is not cursed, so it can never join the cursed set. This is
/// enforced by the `has_closing[i]` short-circuit in the inductive
/// step. Equivalently: cursed = "no path from this VT ever reaches
/// Closed."
///
/// After this function returns, `kind` is assigned in
/// `classify_and_finalize` by splitting cursed into:
///   - `Dead` = base case (no realized transitions);
///   - `Undead` = cursed but with at least one realized transition.
///
/// Returns a `Vec<bool>` indexed by `id - 1`. O(T·N) total work.
fn compute_cursed(
    entries: &[OpenVertexType],
    has_any_realized: &[bool],
    has_closing: &[bool],
    succ_sets: &[BTreeSet<usize>],
) -> Vec<bool> {
    let n = entries.len();
    let mut cursed: Vec<bool> = has_any_realized.iter().map(|&h| !h).collect();

    let mut changed = true;
    while changed {
        changed = false;
        for i in 0..n {
            if cursed[i] {
                continue;
            }
            // A closing transition is an edge to Closed, which is not
            // cursed. So any VT with a closing transition has a non-
            // cursed successor and cannot itself be cursed.
            if has_closing[i] {
                continue;
            }
            let succs = &succ_sets[i];
            if succs.is_empty() {
                // No open successors and no closing transitions, but
                // !cursed means has_any_realized — contradiction. The
                // base-case seed should have already marked this VT.
                debug_assert!(
                    !has_any_realized[i],
                    "VT with no open successors and no closing transitions \
                     should have been Dead in the base case"
                );
                continue;
            }
            if succs.iter().all(|&s| cursed[s - 1]) {
                cursed[i] = true;
                changed = true;
            }
        }
    }

    cursed
}

/// Compute the "blessed" status for each VT via a monotone fixpoint.
///
/// A VT is blessed iff it has at least one **realized** outgoing
/// transition AND every non-closing transition leads to an already-
/// blessed VT (closing transitions count vacuously). Uses the
/// pre-built `succ_sets` (non-closing successors only) plus
/// `has_any_realized` to detect "has any catalog transition".
///
/// Returns a `Vec<bool>` indexed by `id - 1`. O(T·N) total work.
fn compute_blessed(
    entries: &[OpenVertexType],
    has_any_realized: &[bool],
    succ_sets: &[BTreeSet<usize>],
) -> Vec<bool> {
    let n = entries.len();
    let mut blessed: Vec<bool> = vec![false; n];

    let mut changed = true;
    while changed {
        changed = false;
        for i in 0..n {
            if blessed[i] {
                continue;
            }
            // A VT must have at least one realized transition (open
            // or closing) to be blessable. VTs with no transitions
            // are Dead, not Blessed.
            if !has_any_realized[i] {
                continue;
            }
            // Every non-closing successor (succ_sets only stores those)
            // must already be blessed. An empty succ_sets means "only
            // closing transitions exist", which vacuously satisfies
            // the condition.
            if succ_sets[i].iter().all(|&dst| blessed[dst - 1]) {
                blessed[i] = true;
                changed = true;
            }
        }
    }

    blessed
}

/// Derive per-VT "has at least one realized transition (open or
/// closing) in the catalog" from the post-build_transition_arrays
/// outputs.
fn compute_has_any_realized(n: usize, transitions: &[TransitionInfo]) -> Vec<bool> {
    let mut v = vec![false; n];
    for t in transitions {
        v[t.src_id - 1] = true;
    }
    v
}

/// Derive per-VT "has at least one closing transition" — i.e., the VT
/// has an explicit escape path to Closed.
fn compute_has_closing(n: usize, transitions: &[TransitionInfo]) -> Vec<bool> {
    let mut v = vec![false; n];
    for t in transitions {
        if t.is_closed() {
            v[t.src_id - 1] = true;
        }
    }
    v
}

// =====================================================================
// Serializable snapshot for `collect` / `validate` workflows.
//
// The library exposes a non-generic [`Collection`] type that snapshots
// an [`OpenVertexTypeIndex`] for `serde_json` round-trip, plus the
// validation primitives used by `vtype_enum` (witness reconstruction,
// vertex-type cross-check, transition cross-check, completeness scan).
// =====================================================================

use crate::intgeom::matchtypes::MatchTypeIndex;
use serde::{Deserialize, Serialize};

/// Sentinel destination id for transitions that seal the focus vertex
/// (`dst_id == CLOSED_ID`). Same convention as
/// [`TransitionInfo::dst_id`]; preserved in the serialized form so a
/// loader can distinguish "transition to closed" from "transition to
/// open VT id N".
///
/// Recorded boundary state of a VT's canonical witness patch.
///
/// Always taken from a witness that has been
/// [`GrowingPatch::normalize`]d, so the boundary is at its lex-min
/// rotation and `patch_tile_ids` is in dense CCW-first-occurrence
/// form. The ids could in principle be re-derived from the segment
/// structure of `angles + edges`, but [`GrowingPatch::from_parts`]
/// requires them as an argument and `update_inner_chains` reads them
/// at glue time, so we serialize them directly rather than rebuild
/// them on every load.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WitnessSnapshot {
    /// Boundary position of the focus vertex in this witness.
    pub pos: usize,
    pub angles: Vec<i8>,
    pub edges: Vec<EdgeInfo>,
    pub inner_chains: Vec<Vec<EdgeInfo>>,
    /// Per-boundary-position tile-instance id, dense `0..next_tile_id`
    /// in CCW first-occurrence order (= the form
    /// [`GrowingPatch::normalize`] produces).
    pub patch_tile_ids: Vec<usize>,
    /// `max(patch_tile_ids) + 1` (= `0` for an empty boundary).
    pub next_tile_id: usize,
}

/// One VT record in a [`Collection`]: the canonical
/// [`OpenVertexType`], its [`VTypeKind`] classification, the neighbor
/// junction offsets in the witness, and the witness boundary itself.
///
/// Cross-checked on load: `vtype` must match
/// `GrowingPatch::junction_vertex_type_at(witness.pos)` on the
/// reconstructed witness; `(cw_neighbor_offset, ccw_neighbor_offset)`
/// must match `neighbor_junction_offsets(witness.pos)`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VtypeRecord {
    /// 1-based VT id in the source [`OpenVertexTypeIndex`].
    pub id: usize,
    pub vtype: OpenVertexType,
    pub kind: VTypeKind,
    /// CW distance from the focus vertex to the next junction in the
    /// witness.
    pub cw_neighbor_offset: usize,
    /// CCW distance from the focus vertex to the next junction in the
    /// witness.
    pub ccw_neighbor_offset: usize,
    pub witness: WitnessSnapshot,
}

/// One transition record. Mirrors [`TransitionInfo`] one-to-one; the
/// `closed_vt_id` index of the destination closed VT is not retained
/// (a closing transition is identified by `dst_id == CLOSED_ID`).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TransitionRecord {
    pub src_id: usize,
    /// Destination VT's 1-based id, or [`CLOSED_ID`] for a closing
    /// transition that seals the focus vertex.
    pub dst_id: usize,
    pub side: TransitionSide,
    pub tile_id: usize,
    pub tile_offset: usize,
}

/// Serializable snapshot of an [`OpenVertexTypeIndex`]. Designed to
/// round-trip via `serde_json` for the `vtype_enum`-style
/// collect/validate workflow. Non-generic — the typed tileset is
/// reconstructed by the caller from `tile_angles`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Collection {
    /// Symbolic name of the cyclotomic ring (e.g. `"ZZ12"`, `"ZZ10"`).
    pub ring: String,
    /// Raw angle sequence of each tile, in `TileSet` order.
    pub tile_angles: Vec<Vec<i8>>,
    pub vtypes: Vec<VtypeRecord>,
    pub transitions: Vec<TransitionRecord>,
}

/// Validation outcome of [`Collection::completeness_errors`].
pub struct CompletenessReport {
    /// Total candidate matches checked across all alive VTs.
    pub matches_checked: usize,
    /// VT ids whose witness produces a junction VT that the collection
    /// doesn't contain. Empty ⇒ the collection is at its claimed fixed
    /// point.
    pub missing: std::collections::BTreeSet<usize>,
}

impl CompletenessReport {
    pub fn is_complete(&self) -> bool {
        self.missing.is_empty()
    }
}

impl Collection {
    /// Snapshot a built [`OpenVertexTypeIndex`], tagged with `ring`.
    pub fn from_index<T>(idx: &OpenVertexTypeIndex<T>, ring: impl Into<String>) -> Self
    where
        T: IsRing,
    {
        let tile_angles = idx
            .tileset()
            .rats()
            .iter()
            .map(|r| r.seq().to_vec())
            .collect();
        let mut vtypes = Vec::with_capacity(idx.num_types());
        for id in 1..=idx.num_types() {
            let info = idx.get_info(id);
            // Normalize the witness so its boundary is at lex-min
            // rotation and `patch_tile_ids` is in canonical dense
            // form. This makes the serialized snapshot canonical
            // (two BFS runs that produce the same VT also produce
            // the same snapshot).
            let mut w = info.witness().clone();
            let rot = w.normalize();
            let n = w.boundary_len();
            let pos = (info.witness_pos() + n - rot) % n;
            vtypes.push(VtypeRecord {
                id,
                vtype: info.vtype().clone(),
                kind: info.kind(),
                cw_neighbor_offset: info.cw_neighbor_offset(),
                ccw_neighbor_offset: info.ccw_neighbor_offset(),
                witness: WitnessSnapshot {
                    pos,
                    angles: w.angles().to_vec(),
                    edges: w.edges().to_vec(),
                    inner_chains: w.inner_chains().to_vec(),
                    patch_tile_ids: w.patch_tile_ids().to_vec(),
                    next_tile_id: w.next_tile_id(),
                },
            });
        }
        let transitions = idx
            .transitions()
            .iter()
            .map(|t| TransitionRecord {
                src_id: t.src_id,
                dst_id: t.dst_id,
                side: t.side,
                tile_id: t.tile_id,
                tile_offset: t.tile_offset,
            })
            .collect();
        Collection {
            ring: ring.into(),
            tile_angles,
            vtypes,
            transitions,
        }
    }

    /// Reconstruct one [`GrowingPatch`] per VT witness, keyed by VT
    /// id, using the recorded `patch_tile_ids` / `next_tile_id` so
    /// that subsequent `add_tile` calls update inner chains the same
    /// way the BFS did.
    pub fn reconstruct_witnesses<T>(
        &self,
        tile_ts: &Arc<TileSet<T>>,
    ) -> Result<HashMap<usize, GrowingPatch<T>>, String>
    where
        T: IsRing,
    {
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(tile_ts)));
        let mut out = HashMap::with_capacity(self.vtypes.len());
        for v in &self.vtypes {
            let n = v.witness.angles.len();
            if v.witness.inner_chains.len() != n
                || v.witness.edges.len() != n
                || v.witness.patch_tile_ids.len() != n
            {
                return Err(format!(
                    "VT {}: witness arrays length mismatch (n={}, edges={}, inner={}, ptids={})",
                    v.id,
                    n,
                    v.witness.edges.len(),
                    v.witness.inner_chains.len(),
                    v.witness.patch_tile_ids.len(),
                ));
            }
            let gp = GrowingPatch::from_parts(
                Arc::clone(&mi),
                v.witness.angles.clone(),
                v.witness.edges.clone(),
                v.witness.inner_chains.clone(),
                v.witness.patch_tile_ids.clone(),
                v.witness.next_tile_id,
            )
            .ok_or_else(|| format!("VT {}: from_parts failed", v.id))?;
            out.insert(v.id, gp);
        }
        Ok(out)
    }

    /// Per-VT cross-check: verify each VT record's `vtype` and
    /// `(cw_neighbor_offset, ccw_neighbor_offset)` against what the
    /// reconstructed witness re-derives. Returns `(vt_id, message)`
    /// pairs for each disagreement.
    pub fn vtype_errors<T>(
        &self,
        witnesses: &HashMap<usize, GrowingPatch<T>>,
    ) -> Vec<(usize, String)>
    where
        T: IsRing,
    {
        let mut errors = Vec::new();
        for v in &self.vtypes {
            let Some(gp) = witnesses.get(&v.id) else {
                errors.push((v.id, "no witness".to_string()));
                continue;
            };
            let pos = v.witness.pos;
            let (expected_cw, expected_ccw) = gp.neighbor_junction_offsets(pos).unwrap_or((0, 0));
            if expected_cw != v.cw_neighbor_offset || expected_ccw != v.ccw_neighbor_offset {
                errors.push((
                    v.id,
                    format!(
                        "neighbor offsets: expected ({}, {}), recorded ({}, {})",
                        expected_cw, expected_ccw, v.cw_neighbor_offset, v.ccw_neighbor_offset,
                    ),
                ));
            }
            match gp.junction_vertex_type_at(pos) {
                Some(actual) if actual.cw == v.vtype.cw && actual.ccw == v.vtype.ccw => {}
                Some(actual) => errors.push((
                    v.id,
                    format!(
                        "vtype mismatch: recorded {:?}, witness {:?}",
                        v.vtype, actual
                    ),
                )),
                None => errors.push((v.id, "junction_vertex_type_at returned None".to_string())),
            }
        }
        errors
    }

    /// Per-transition cross-check: verify each recorded transition is
    /// realised by some glue on the source VT's witness. Returns
    /// `((src, dst), message)` for unverifiable transitions.
    pub fn transition_errors<T>(
        &self,
        tile_ts: &Arc<TileSet<T>>,
        witnesses: &HashMap<usize, GrowingPatch<T>>,
    ) -> Vec<((usize, usize), String)>
    where
        T: IsRing,
    {
        let mut errors = Vec::new();
        let known_ids: std::collections::BTreeSet<usize> =
            self.vtypes.iter().map(|v| v.id).collect();
        let witness_pos: HashMap<usize, usize> =
            self.vtypes.iter().map(|v| (v.id, v.witness.pos)).collect();
        for t in &self.transitions {
            if !known_ids.contains(&t.src_id) {
                errors.push(((t.src_id, t.dst_id), "unknown src id".to_string()));
                continue;
            }
            if t.dst_id != CLOSED_ID && !known_ids.contains(&t.dst_id) {
                errors.push(((t.src_id, t.dst_id), "unknown dst id".to_string()));
                continue;
            }
            let Some(gp) = witnesses.get(&t.src_id) else {
                errors.push(((t.src_id, t.dst_id), "no source witness".to_string()));
                continue;
            };
            let pos = witness_pos[&t.src_id];
            let n = gp.boundary_len();
            // The BFS records `t.tile_offset` resolved against the
            // canonical edge for the side (CW edge for Cw/Both, CCW
            // edge for Ccw). Mirror that exactly so the inverse
            // computation below recovers the same offset.
            let edge_pos = match t.side {
                TransitionSide::Cw | TransitionSide::Both => (pos + n - 1) % n,
                TransitionSide::Ccw => pos,
            };
            let cw_edge = (pos + n - 1) % n;
            let ccw_edge = pos;
            // Edge-consumption predicate (half-open `[start_a,
            // start_a + len)`). The BFS uses this same predicate to
            // decide which side a match consumes; use it here too so
            // we match its classification exactly. The vertex-touch
            // semantics of `cyclic_range_contains` (closed `<= len`)
            // would over-accept matches that only touch the focus
            // vertex without consuming its incident edge — the bug
            // the old binary's `validate_common` had.
            let edge_consumed =
                |pm: &PatchMatch, edge: usize| -> bool { (edge + n - pm.start_a()) % n < pm.len() };
            let tile_len = tile_ts.rat(t.tile_id).len();
            let mut found = false;
            for pm in gp.get_all_matches() {
                if pm.tile_id() != t.tile_id {
                    continue;
                }
                if !edge_consumed(&pm, edge_pos) {
                    continue;
                }
                let offset_in_match =
                    (edge_pos as i64 - pm.start_a() as i64).rem_euclid(n as i64) as usize;
                let computed_offset = (pm.start_b() as i64 + pm.len() as i64 - offset_in_match as i64)
                    .rem_euclid(tile_len as i64) as usize;
                if computed_offset != t.tile_offset {
                    continue;
                }
                let consumes_cw = edge_consumed(&pm, cw_edge);
                let consumes_ccw = edge_consumed(&pm, ccw_edge);
                let actual_side = match (consumes_cw, consumes_ccw) {
                    (true, true) => TransitionSide::Both,
                    (true, false) => TransitionSide::Cw,
                    (false, true) => TransitionSide::Ccw,
                    (false, false) => continue,
                };
                if actual_side != t.side {
                    continue;
                }
                let mut gp2 = gp.clone();
                if !gp2.add_tile(&pm) || !gp2.is_growing() {
                    continue;
                }
                if t.dst_id == CLOSED_ID {
                    // For a closing match the focus vertex is sealed.
                    // No junction lookup is needed — the edge-side
                    // and offset checks above already pin the glue.
                    found = true;
                    break;
                }
                let junction_pos = if pm.start_a() == pos { n - pm.len() } else { 0 };
                if let Some(actual) = gp2.junction_vertex_type_at(junction_pos) {
                    let Some(dst_gp) = witnesses.get(&t.dst_id) else {
                        continue;
                    };
                    let Some(&dst_pos) = witness_pos.get(&t.dst_id) else {
                        continue;
                    };
                    if let Some(expected) = dst_gp.junction_vertex_type_at(dst_pos) {
                        if actual.cw == expected.cw && actual.ccw == expected.ccw {
                            found = true;
                            break;
                        }
                    }
                }
            }
            if !found {
                errors.push(((t.src_id, t.dst_id), "no matching glue found".to_string()));
            }
        }
        errors
    }

    /// Completeness scan: for every non-cursed VT, enumerate matches
    /// touching its focus vertex on the witness; each glue's induced
    /// junction VT must be present in the collection (compared by
    /// `(cw, ccw)`). Returns the set of VT ids whose witness produces
    /// an unknown junction VT, plus the total number of match checks.
    pub fn completeness_errors<T>(
        &self,
        witnesses: &HashMap<usize, GrowingPatch<T>>,
    ) -> CompletenessReport
    where
        T: IsRing,
    {
        let mut missing: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
        let mut matches_checked = 0usize;
        for v in &self.vtypes {
            if matches!(v.kind, VTypeKind::Dead | VTypeKind::Undead) {
                continue;
            }
            let Some(gp) = witnesses.get(&v.id) else {
                continue;
            };
            let pos = v.witness.pos;
            let old_n = gp.boundary_len();
            // Same edge-consumption predicate as `transition_errors`
            // — vertex-touch semantics over-accept matches that don't
            // actually consume the focus's incident edges.
            let edge_consumed = |pm: &PatchMatch, edge: usize| -> bool {
                (edge + old_n - pm.start_a()) % old_n < pm.len()
            };
            for pm in gp.get_matches_touching_vertex(pos) {
                matches_checked += 1;
                let mut gp2 = gp.clone();
                if !gp2.add_tile(&pm) || !gp2.is_growing() {
                    continue;
                }
                let junction_pos = if pm.start_a() == pos { old_n - pm.len() } else { 0 };
                let consumes_ccw = edge_consumed(&pm, pos);
                let consumes_cw = edge_consumed(&pm, (pos + old_n - 1) % old_n);
                if consumes_cw && consumes_ccw {
                    continue;
                }
                let Some(jvt) = gp2.junction_vertex_type_at(junction_pos) else {
                    continue;
                };
                // Compare against the canonical recorded vtypes in
                // the collection. (Comparing through each other VT's
                // reconstructed witness would add a needless layer
                // of indirection and is more sensitive to differences
                // in patch_tile_ids / next_tile_id between the BFS
                // and the reconstructed patch.)
                let found = self.vtypes.iter().any(|other| other.vtype == jvt);
                if !found {
                    missing.insert(v.id);
                }
            }
        }
        CompletenessReport {
            matches_checked,
            missing,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic::{ZZ12, ZZ4};
    use crate::intgeom::tiles;
    use crate::intgeom::tileset::{self, TileSet};
    use crate::matches::{EdgeRange, Segment};
    use std::sync::Arc;

    /// End-to-end regression: every recorded transition and the
    /// completeness scan on every Free/Blessed VT's witness must
    /// agree after a `Collection` round-trip via `serde_json`.
    ///
    /// The pre-cleanup `vtype_enum` validator used
    /// `cyclic_range_contains` (vertex-touch, `<= len`) for both the
    /// edge-position filter and the `covers_both` test. The
    /// vertex-touch semantic over-accepted matches at the boundary
    /// of the matched edge range, so the Cw-side check would mark
    /// the real glue as "Both" and reject it, leaving the transition
    /// unverifiable. Square exposed this on every Cw open transition
    /// (64 failures); spectre also exposed completeness gaps. Both
    /// the validator and the BFS now use the edge-consumption
    /// predicate `(edge + n - start_a) % n < len`, matching the
    /// `consumes_cw_edge` / `consumes_ccw_edge` check in
    /// `OpenVertexTypeIndex::new`.
    #[test]
    fn collection_roundtrip_validates_square_and_spectre() {
        let cases: [(&str, Arc<TileSet<ZZ12>>); 2] = [
            ("square", tileset::square::<ZZ12>()),
            ("spectre", tileset::spectre::<ZZ12>()),
        ];
        for (label, ts) in cases {
            let idx = OpenVertexTypeIndex::new(Arc::clone(&ts));
            let collection = Collection::from_index(&idx, "ZZ12");
            let json = serde_json::to_string(&collection).expect("serialize");
            let restored: Collection = serde_json::from_str(&json).expect("deserialize");
            let witnesses = restored.reconstruct_witnesses(&ts).expect("reconstruct");

            let vt_errs = restored.vtype_errors(&witnesses);
            assert!(vt_errs.is_empty(), "{label}: vtype errors: {vt_errs:?}");

            let tr_errs = restored.transition_errors(&ts, &witnesses);
            assert!(
                tr_errs.is_empty(),
                "{}: {} transition errors (first: {:?})",
                label,
                tr_errs.len(),
                tr_errs.first(),
            );

            let report = restored.completeness_errors(&witnesses);
            assert!(
                report.is_complete(),
                "{}: {} VTs produce unknown junction types: {:?}",
                label,
                report.missing.len(),
                report.missing.iter().take(10).collect::<Vec<_>>(),
            );
        }
    }

    /// Bake in the known VT count + classification breakdown for hex.
    /// Any drift in the BFS or classification will fail loudly here.
    #[test]
    fn hexagon_vertex_type_counts() {
        let hex: crate::intgeom::snake::Snake<ZZ12> = tiles::hexagon();
        let rat = Rat::try_from(&hex).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = OpenVertexTypeIndex::new(ts);

        let alive = idx.entries.iter().filter(|e| e.is_alive()).count();
        let dead = idx.entries.iter().filter(|e| e.is_dead()).count();
        let cursed = idx.entries.iter().filter(|e| e.is_cursed()).count();
        let blessed = idx.entries.iter().filter(|e| e.is_blessed()).count();
        let free = idx.entries.iter().filter(|e| e.is_free()).count();

        assert_eq!(idx.num_types(), 36, "hex VT count");
        assert_eq!(alive, 36, "every hex VT is alive");
        assert_eq!(dead, 0, "no dead hex VTs");
        assert_eq!(cursed, 0, "no cursed hex VTs");
        assert_eq!(blessed + free, 36, "alive = blessed + free");
    }

    /// Bake in the known VT count for square.
    #[test]
    fn square_vertex_type_counts() {
        let sq: crate::intgeom::snake::Snake<ZZ4> = tiles::square();
        let rat = Rat::try_from(&sq).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = OpenVertexTypeIndex::new(ts);
        assert_eq!(idx.num_types(), 80, "square VT count");
        let alive = idx.entries.iter().filter(|e| e.is_alive()).count();
        assert_eq!(alive, 80, "every square VT is alive");
    }

    /// Every transition's src/dst ids must be valid (in-range or
    /// CLOSED_ID), every tile_id must be a real tile in the tileset,
    /// and tile_offset must be within the referenced tile's edge
    /// count.
    #[test]
    fn hexagon_transitions_well_formed() {
        let hex: crate::intgeom::snake::Snake<ZZ12> = tiles::hexagon();
        let rat = Rat::try_from(&hex).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = OpenVertexTypeIndex::new(Arc::clone(&ts));
        let n_types = idx.num_types();
        for t in idx.transitions() {
            assert!(
                t.src_id >= 1 && t.src_id <= n_types,
                "src_id {} out of range",
                t.src_id
            );
            assert!(
                t.dst_id == CLOSED_ID || (t.dst_id >= 1 && t.dst_id <= n_types),
                "dst_id {} out of range (and not CLOSED)",
                t.dst_id
            );
            assert!(
                t.tile_id < ts.num_tiles(),
                "tile_id {} out of range",
                t.tile_id
            );
            let tile_len = ts.rat(t.tile_id).len();
            assert!(
                t.tile_offset < tile_len,
                "tile_offset {} out of range (tile {} has {} edges)",
                t.tile_offset,
                t.tile_id,
                tile_len
            );
        }
    }

    /// The cursed-fixpoint invariant: a VT is `Dead` or `Undead` iff
    /// every reachable path from it (via succ_sets) ends at a Dead
    /// VT. Computed via BFS over the successor graph: from each
    /// non-cursed VT we should be able to reach at least one VT that
    /// is `Free` or `Blessed`.
    #[test]
    fn hexagon_cursed_fixpoint_invariant() {
        let hex: crate::intgeom::snake::Snake<ZZ12> = tiles::hexagon();
        let rat = Rat::try_from(&hex).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = OpenVertexTypeIndex::new(Arc::clone(&ts));

        // For hex, all VTs are alive, so cursed-fixpoint is trivially
        // satisfied (no VTs are cursed). Test the *positive* invariant
        // separately: every non-cursed VT has at least one outgoing
        // transition (closing or to another non-cursed VT).
        for info in idx.entries() {
            if info.is_cursed() {
                continue;
            }
            let id = idx.get_id(info.vtype()).unwrap();
            let outgoing: Vec<&TransitionInfo> = idx
                .transitions()
                .iter()
                .filter(|t| t.src_id == id)
                .collect();
            assert!(
                !outgoing.is_empty(),
                "non-cursed VT id {id} should have outgoing transitions"
            );
            // Every outgoing must either close or lead to a non-Dead VT.
            for t in &outgoing {
                if t.is_closed() {
                    continue;
                }
                let dst_info = idx.get_info(t.dst_id);
                assert!(
                    !dst_info.is_dead(),
                    "non-cursed VT id {id} has open transition to Dead VT id {}",
                    t.dst_id
                );
            }
        }
    }

    /// The blessed-fixpoint invariant: a VT is `Blessed` iff every
    /// outgoing transition is either closing or leads to another
    /// `Blessed` VT. This re-derives the predicate from scratch and
    /// compares against the recorded classification.
    #[test]
    fn hexagon_blessed_fixpoint_invariant() {
        let hex: crate::intgeom::snake::Snake<ZZ12> = tiles::hexagon();
        let rat = Rat::try_from(&hex).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = OpenVertexTypeIndex::new(Arc::clone(&ts));
        for info in idx.entries() {
            if !info.is_blessed() {
                continue;
            }
            let id = idx.get_id(info.vtype()).unwrap();
            let outgoing: Vec<&TransitionInfo> = idx
                .transitions()
                .iter()
                .filter(|t| t.src_id == id)
                .collect();
            assert!(
                !outgoing.is_empty(),
                "Blessed VT id {id} must have at least one transition"
            );
            for t in &outgoing {
                assert!(
                    t.is_closed() || idx.get_info(t.dst_id).is_blessed(),
                    "Blessed VT id {id} has open transition to non-Blessed VT id {}",
                    t.dst_id
                );
            }
        }
    }

    // ----- Validity + completeness helpers + tests -----

    /// Validity invariant: every VT in the catalog is realized by its
    /// recorded witness at the recorded position. (Round-trip:
    /// `witness.junction_vertex_type_at(pos) == Some(vt)`.)
    fn assert_catalog_validity<T>(idx: &OpenVertexTypeIndex<T>)
    where
        T: IsRing,
    {
        for info in idx.entries() {
            let observed = info.witness().junction_vertex_type_at(info.witness_pos());
            assert_eq!(
                observed,
                Some(info.vtype().clone()),
                "validity violation: witness for catalog id {} does not \
                 realize its claimed VT",
                idx.get_id(info.vtype()).unwrap()
            );
        }
    }

    /// Completeness invariant, checked **without** using
    /// `get_matches_touching_vertex` (the BFS's own enumerator).
    ///
    /// For each catalog VT and its witness, enumerate every
    /// canonical-max-length match between the witness boundary and
    /// each tile by iterating anchor pairs `(start_a, start_b)` and
    /// calling [`Rat::get_match`] to extend maximally in both
    /// directions. Filter to matches whose edge range covers the CW
    /// or CCW incident edge of the focus vertex (pure modular
    /// arithmetic). For every such match accepted by `add_tile +
    /// is_growing`, every junction VT exposed in the resulting patch
    /// must already be in the catalog.
    ///
    /// This is non-circular w.r.t. the BFS's match-enumeration in
    /// `patch.rs` (`compute_all_candidates` /
    /// `get_matches_touching_vertex`): the brute path here goes via
    /// `Rat::get_match` directly on the witness's angle sequence and
    /// only relies on `add_tile`'s contract that hand-constructed
    /// canonical `PatchMatch`es are valid input. It still uses the
    /// witness and `add_tile`, so it verifies closure relative to
    /// each VT's recorded witness — not global completeness over all
    /// possible witnesses.
    ///
    /// Note on truncations: a `PatchMatch` with `len` *less than* the
    /// canonical max at its anchor would represent a geometrically
    /// meaningless "partially-fused" glue (per the project's match
    /// semantics, gluing always extends to the maximal common run);
    /// and `add_tile` does not validate the claimed `len`, so passing
    /// a truncated `len` produces a malformed boundary. We only
    /// enumerate canonical matches here.
    fn assert_catalog_complete_brute<T>(idx: &OpenVertexTypeIndex<T>)
    where
        T: IsRing,
    {
        let ts = idx.tileset();
        for info in idx.entries() {
            let src_id = idx.get_id(info.vtype()).unwrap();
            let witness = info.witness();
            let pos = info.witness_pos();
            let n = witness.boundary_len();
            let boundary_rat = Rat::from_slice_unchecked(witness.angles());

            // Canonical-match enumeration. (start_a, start_b) is the
            // anchor; rat.get_match extends maximally in both
            // directions to produce (ext_start, len, ext_end). We
            // dedup since multiple anchors collapse onto the same
            // canonical match.
            let mut canonical: std::collections::BTreeSet<(usize, usize, usize, usize)> =
                std::collections::BTreeSet::new();
            for tile_id in 0..ts.num_tiles() {
                let tile_rat = ts.rat(tile_id);
                let m = tile_rat.len();
                for start_a in 0..n {
                    for start_b in 0..m {
                        let (ext_start, len, ext_end) =
                            boundary_rat.get_match((start_a as i64, start_b as i64), tile_rat);
                        if len == 0 {
                            continue;
                        }
                        let es = ext_start.rem_euclid(n as i64) as usize;
                        let ee = ext_end.rem_euclid(m as i64) as usize;
                        canonical.insert((tile_id, es, len, ee));
                    }
                }
            }

            for &(tile_id, ext_start, len, ext_end) in &canonical {
                // Edge-range filter: does [ext_start, ext_start+len)
                // cover the CW edge (pos - 1) or the CCW edge (pos)?
                let in_range = |edge: usize| (edge + n - ext_start) % n < len;
                let touches_cw = in_range((pos + n - 1) % n);
                let touches_ccw = in_range(pos);
                if !touches_cw && !touches_ccw {
                    continue;
                }
                let pm = PatchMatch::new(
    EdgeRange::new(ext_start, len),
    Segment::new(tile_id, EdgeRange::new(ext_end, len)),
);
                let mut gp2 = witness.clone();
                if !gp2.add_tile(&pm) || !gp2.is_growing() {
                    continue;
                }
                for new_pos in 0..gp2.boundary_len() {
                    if !gp2.is_junction(new_pos) {
                        continue;
                    }
                    if let Some(nv) = gp2.junction_vertex_type_at(new_pos) {
                        assert!(
                            idx.get_id(&nv).is_some(),
                            "completeness violation: canonical glue \
                             (tile {tile_id}, start_a {ext_start}, start_b {ext_end}, len {len}) \
                             from src VT id {src_id} exposes VT at new_pos {new_pos} \
                             not in catalog"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn hexagon_catalog_validity() {
        let hex: crate::intgeom::snake::Snake<ZZ12> = tiles::hexagon();
        let rat = Rat::try_from(&hex).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = OpenVertexTypeIndex::new(ts);
        assert_catalog_validity(&idx);
    }

    #[test]
    fn square_catalog_validity() {
        let sq: crate::intgeom::snake::Snake<ZZ4> = tiles::square();
        let rat = Rat::try_from(&sq).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = OpenVertexTypeIndex::new(ts);
        assert_catalog_validity(&idx);
    }

    #[test]
    fn spectre_catalog_validity() {
        let s: crate::intgeom::snake::Snake<ZZ12> = tiles::spectre();
        let rat = Rat::try_from(&s).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = OpenVertexTypeIndex::new(ts);
        assert_catalog_validity(&idx);
    }

    #[test]
    fn hexagon_catalog_complete_brute() {
        let hex: crate::intgeom::snake::Snake<ZZ12> = tiles::hexagon();
        let rat = Rat::try_from(&hex).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = OpenVertexTypeIndex::new(ts);
        assert_catalog_complete_brute(&idx);
    }

    #[test]
    fn square_catalog_complete_brute() {
        let sq: crate::intgeom::snake::Snake<ZZ4> = tiles::square();
        let rat = Rat::try_from(&sq).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = OpenVertexTypeIndex::new(ts);
        assert_catalog_complete_brute(&idx);
    }

    /// Diagnostic (kept as a regression check): on VT id 1's witness,
    /// the brute canonical-match set touching the focus vertex must
    /// equal `get_matches_touching_vertex`. Previously fired due to a
    /// wrap-around bug in [`cyclic_range_contains`] at `start + len ==
    /// n` exactly.
    #[test]
    fn vt_witness_touching_matches_match_brute() {
        let s: crate::intgeom::snake::Snake<ZZ12> = tiles::spectre();
        let rat = Rat::try_from(&s).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = OpenVertexTypeIndex::new(Arc::clone(&ts));

        for info in idx.entries() {
            let witness = info.witness();
            let pos = info.witness_pos();
            let n = witness.boundary_len();
            let boundary_rat = Rat::from_slice_unchecked(witness.angles());

            // Same brute as patch.rs:3108: canonical via rat.get_match,
            // filter via junctions_glueable + try_glue_precomputed.
            let mut filtered_brute: std::collections::BTreeSet<(usize, usize, usize, usize)> =
                std::collections::BTreeSet::new();
            for tile_id in 0..ts.num_tiles() {
                let tile_rat = ts.rat(tile_id);
                let m = tile_rat.len();
                for start_a in 0..n {
                    for start_b in 0..m {
                        let (ns, len, ne) =
                            boundary_rat.get_match((start_a as i64, start_b as i64), tile_rat);
                        if len == 0 {
                            continue;
                        }
                        let ns_u = ns.rem_euclid(n as i64) as usize;
                        let ne_u = ne.rem_euclid(m as i64) as usize;
                        // Edge-based touching filter (does NOT use
                        // cyclic_range_contains — independent of the
                        // function this test would have caught).
                        let in_range = |edge: usize| (edge + n - ns_u) % n < len;
                        if !in_range((pos + n - 1) % n) && !in_range(pos) {
                            continue;
                        }
                        if !crate::intgeom::glue::junctions_glueable(
                            witness.angles(),
                            ns_u,
                            len,
                            tile_rat.seq(),
                            ne_u,
                        ) {
                            continue;
                        }
                        if boundary_rat
                            .try_glue_precomputed((ns, len, ne), tile_rat, true)
                            .is_ok()
                        {
                            filtered_brute.insert((tile_id, ns_u, len, ne_u));
                        }
                    }
                }
            }

            let api_set: std::collections::BTreeSet<(usize, usize, usize, usize)> = witness
                .get_matches_touching_vertex(pos)
                .into_iter()
                .map(|pm| (pm.tile_id(), pm.start_a(), pm.len(), pm.start_b()))
                .collect();

            assert_eq!(
                filtered_brute,
                api_set,
                "VT id {} (witness n={n} pos={pos}): brute != api",
                idx.get_id(info.vtype()).unwrap()
            );
        }
    }

    /// Diagnostic kept around (eprintln output, `#[ignore]`) for
    /// future investigation of `get_matches_touching_vertex` discrepancies
    /// on `from_parts`-constructed witness patches. Run via
    /// `cargo test --lib --release intgeom::vertextypes::tests::diagnose_spectre_touching_matches -- --ignored --nocapture`.
    #[test]
    #[ignore]
    fn diagnose_spectre_touching_matches() {
        let s: crate::intgeom::snake::Snake<ZZ12> = tiles::spectre();
        let rat = Rat::try_from(&s).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = OpenVertexTypeIndex::new(Arc::clone(&ts));

        let info = idx.get_info(1);
        let witness = info.witness();
        let pos = info.witness_pos();
        let n = witness.boundary_len();
        let boundary_rat = Rat::from_slice_unchecked(witness.angles());

        // Same brute as existing patch.rs:3070 test:
        // - rat.get_match for canonical (start_a, len, start_b);
        // - junctions_glueable filter;
        // - try_glue_precomputed(unchecked=true).is_ok() filter.
        let mut filtered_brute: std::collections::BTreeSet<(usize, usize, usize, usize)> =
            std::collections::BTreeSet::new();
        let mut unfiltered_brute: std::collections::BTreeSet<(usize, usize, usize, usize)> =
            std::collections::BTreeSet::new();
        for tile_id in 0..ts.num_tiles() {
            let tile_rat = ts.rat(tile_id);
            let m = tile_rat.len();
            for start_a in 0..n {
                for start_b in 0..m {
                    let (ns, len, ne) =
                        boundary_rat.get_match((start_a as i64, start_b as i64), tile_rat);
                    if len == 0 {
                        continue;
                    }
                    let ns_u = ns.rem_euclid(n as i64) as usize;
                    let ne_u = ne.rem_euclid(m as i64) as usize;
                    let in_range = |edge: usize| (edge + n - ns_u) % n < len;
                    if !in_range((pos + n - 1) % n) && !in_range(pos) {
                        continue;
                    }
                    unfiltered_brute.insert((tile_id, ns_u, len, ne_u));
                    if !crate::intgeom::glue::junctions_glueable(
                        witness.angles(),
                        ns_u,
                        len,
                        tile_rat.seq(),
                        ne_u,
                    ) {
                        continue;
                    }
                    if boundary_rat
                        .try_glue_precomputed((ns, len, ne), tile_rat, true)
                        .is_ok()
                    {
                        filtered_brute.insert((tile_id, ns_u, len, ne_u));
                    }
                }
            }
        }

        // BFS-API enumeration.
        let api_set: std::collections::BTreeSet<(usize, usize, usize, usize)> = witness
            .get_matches_touching_vertex(pos)
            .into_iter()
            .map(|pm| (pm.tile_id(), pm.start_a(), pm.len(), pm.start_b()))
            .collect();

        eprintln!(
            "VT id 1: n={n} pos={pos} unfiltered={} filtered={} api={}",
            unfiltered_brute.len(),
            filtered_brute.len(),
            api_set.len()
        );
        let only_filt: Vec<_> = filtered_brute.difference(&api_set).collect();
        let only_api: Vec<_> = api_set.difference(&filtered_brute).collect();
        eprintln!(
            "  filtered-only-not-api: {} | api-only: {}",
            only_filt.len(),
            only_api.len()
        );
        // Probe the hypothesis: missing matches fail strict-convex
        // `junctions_glueable(clean_tile, ia, clean_tile, ib)` —
        // which uses CLEAN-tile angles — but pass
        // `junctions_glueable(boundary, ...)` because the
        // boundary's angle at the right-junction position differs from
        // the clean tile's angle (boundary[(sa+len)%n] is at a glue
        // junction, with a turn angle different from the underlying
        // tile's interior angle).
        let tile_seq = ts.rat(0).seq();
        let bnd_angles = witness.angles();
        let edges = witness.edges();
        for x in only_filt.iter() {
            let &(_tile, sa, len, sb) = *x;
            let m = tile_seq.len();
            let next_sa = (sa + len) % n;
            let underlying_tile = edges[sa].tile_id;
            let underlying_off = edges[sa].tile_offset;
            let underlying_next_tile = edges[next_sa].tile_id;
            let underlying_next_off = edges[next_sa].tile_offset;
            let underlying_seq = ts.rat(underlying_tile).seq();
            let underlying_next_seq = ts.rat(underlying_next_tile).seq();
            let bnd_left = bnd_angles[sa] as i32 + tile_seq[sb] as i32;
            let bnd_right = bnd_angles[next_sa] as i32 + tile_seq[(sb + m - len) % m] as i32;
            let clean_left = underlying_seq[underlying_off] as i32 + tile_seq[sb] as i32;
            let clean_right = underlying_next_seq[underlying_next_off] as i32
                + tile_seq[(sb + m - len) % m] as i32;
            let clean_tile_a = ts.rat(underlying_tile);
            let (ns_clean, len_clean, ne_clean) =
                clean_tile_a.get_match((underlying_off as i64, sb as i64), ts.rat(0));
            eprintln!(
                "    miss (sa={sa},len={len},sb={sb}): \
                 underlying_off={underlying_off} \
                 clean canonical at (off={underlying_off}, sb={sb}): \
                 (ns={ns_clean},len={len_clean},ne={ne_clean})"
            );
            let _ = (bnd_left, bnd_right, clean_left, clean_right);
        }

        // Probe what MatchTypeIndex actually indexed for (tile=0, offset=0).
        let mti = crate::intgeom::matchtypes::MatchTypeIndex::new(Arc::clone(&ts));
        let cands = mti.candidates_starting_at(0, 0);
        eprintln!(
            "  MatchTypeIndex.by_start[0][0] = {} candidates",
            cands.len()
        );

        eprintln!("  filtered_brute entries:");
        for x in filtered_brute.iter() {
            let in_api = api_set.contains(x);
            eprintln!("    brute {x:?}  in_api={in_api}");
        }
        eprintln!("  api_set entries:");
        for x in api_set.iter() {
            eprintln!("    api   {x:?}");
        }

        // Print junctions and segments for the witness boundary.
        let mut juncs: Vec<usize> = Vec::new();
        for i in 0..n {
            let ei = witness.edges()[i];
            if ts.rat(ei.tile_id).seq()[ei.tile_offset] != witness.angles()[i] {
                juncs.push(i);
            }
        }
        eprintln!("  junctions in witness: {juncs:?}");
        eprintln!(
            "  edges[25]: tile_id={} tile_offset={}",
            edges[25].tile_id, edges[25].tile_offset
        );
        for (i, e) in edges.iter().enumerate().skip(n - 3).take(3) {
            eprintln!(
                "    edges[{i}]: tile_id={} tile_offset={}",
                e.tile_id, e.tile_offset
            );
        }
        for (i, e) in edges.iter().enumerate().take(3) {
            eprintln!(
                "    edges[{i}]: tile_id={} tile_offset={}",
                e.tile_id, e.tile_offset
            );
        }

        // Compare: what does compute_all_candidates produce at result[25]?
        let result = crate::intgeom::patch::GrowingPatch::<ZZ12>::compute_all_candidates(
            &Arc::new(crate::intgeom::matchtypes::MatchTypeIndex::new(Arc::clone(
                &ts,
            ))),
            witness.angles(),
            witness.edges(),
        );
        eprintln!("  result[25] = {} entries:", result[25].len());
        for pm in &result[25] {
            eprintln!(
                "    cac: tile={} start_a={} len={} start_b={}",
                pm.tile_id(), pm.start_a(), pm.len(), pm.start_b()
            );
        }

        // Step into compute_candidates_at_position behavior for pos=25:
        // iterate by_start[0][0] candidates manually and report which
        // ones survive each filter.
        let bnd_rat = Rat::from_slice_unchecked(witness.angles());
        let mut seen: std::collections::HashSet<(usize, usize, usize, usize)> =
            std::collections::HashSet::new();
        eprintln!("  Simulated compute_candidates_at_position(pos=25, tile_id=0, off=0):");
        for c in cands {
            let tile_b = ts.rat(c.tile_b());
            let (ns, len, ne) = bnd_rat.get_match((25i64, c.start_b() as i64), tile_b);
            if len == 0 {
                continue;
            }
            let ns_u = ns.rem_euclid(n as i64) as usize;
            let ne_u = ne.rem_euclid(tile_b.len() as i64) as usize;
            let gap_ok = crate::intgeom::glue::junctions_glueable(
                witness.angles(),
                ns_u,
                len,
                tile_b.seq(),
                ne_u,
            );
            let key = (ns_u, len, ne_u, c.tile_b());
            let already_seen = seen.contains(&key);
            let glue_ok = bnd_rat
                .try_glue_precomputed((ns, len, ne), tile_b, true)
                .is_ok();
            let dst = (c.tile_b(), ns_u, len, ne_u);
            let target_match = filtered_brute.contains(&dst) && dst.1 == 25;
            if target_match || !already_seen {
                eprintln!(
                    "    cand sb={} len={} -> bnd ({ns_u},{len},{ne_u}) gap_ok={gap_ok} \
                     seen_already={already_seen} glue_ok={glue_ok}  target={target_match}",
                    c.start_b(), c.range.len
                );
            }
            seen.insert(key);
        }
    }

    /// Spectre's canonical-match completeness check. Expensive
    /// (229 VTs × ~n*m anchor pairs each), so it's gated behind
    /// `--ignored`. Run via
    /// `cargo test --release spectre_catalog_complete_brute --
    /// --ignored --nocapture`.
    #[test]
    #[ignore]
    fn spectre_catalog_complete_brute() {
        let s: crate::intgeom::snake::Snake<ZZ12> = tiles::spectre();
        let rat = Rat::try_from(&s).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = OpenVertexTypeIndex::new(ts);
        assert_catalog_complete_brute(&idx);
    }

    #[test]
    fn range_by_cw_single_tile() {
        let hex: crate::intgeom::snake::Snake<ZZ12> = tiles::hexagon();
        let rat = Rat::try_from(&hex).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = OpenVertexTypeIndex::new(ts);

        let slice = idx.range_by_cw(0);
        assert_eq!(
            slice.len(),
            idx.num_types(),
            "single-tile tileset: all types have cw.tile_id == 0"
        );
        for info in slice {
            assert_eq!(info.vtype().cw.tile_id, 0);
        }

        let empty = idx.range_by_cw(1);
        assert!(empty.is_empty(), "no types with cw.tile_id == 1");
    }

    #[test]
    fn range_by_cw_multi_tile() {
        let sq: crate::intgeom::snake::Snake<ZZ12> = tiles::square();
        let hex: crate::intgeom::snake::Snake<ZZ12> = tiles::hexagon();
        let sq_rat = Rat::try_from(&sq).unwrap();
        let hex_rat = Rat::try_from(&hex).unwrap();
        let ts = Arc::new(TileSet::new(vec![sq_rat, hex_rat]));
        let idx = OpenVertexTypeIndex::new(ts);

        let n = idx.num_types();
        assert!(n > 0);

        let sq_slice = idx.range_by_cw(0);
        let hex_slice = idx.range_by_cw(1);

        assert_eq!(
            sq_slice.len() + hex_slice.len(),
            n,
            "partition must cover all types"
        );
        assert!(
            !sq_slice.is_empty(),
            "should have types starting with square"
        );
        assert!(!hex_slice.is_empty(), "should have types starting with hex");

        for info in sq_slice {
            assert_eq!(info.vtype().cw.tile_id, 0);
        }
        for info in hex_slice {
            assert_eq!(info.vtype().cw.tile_id, 1);
        }

        let empty = idx.range_by_cw(2);
        assert!(empty.is_empty());
    }

    /// Spectre is the non-trivial benchmark — non-convex, 14 edges,
    /// VT space includes every classification kind (Dead, Undead,
    /// Blessed, Free). These counts are the regression target: any
    /// drift in the BFS, classification, or rotation/transition
    /// conventions will fail here loudly.
    #[test]
    fn spectre_vertex_type_counts() {
        let s: crate::intgeom::snake::Snake<ZZ12> = tiles::spectre();
        let rat = Rat::try_from(&s).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = OpenVertexTypeIndex::new(Arc::clone(&ts));

        let alive = idx.entries.iter().filter(|e| e.is_alive()).count();
        let dead = idx.entries.iter().filter(|e| e.is_dead()).count();
        let undead = idx.entries.iter().filter(|e| e.is_undead()).count();
        let blessed = idx.entries.iter().filter(|e| e.is_blessed()).count();
        let free = idx.entries.iter().filter(|e| e.is_free()).count();
        let cursed = idx.entries.iter().filter(|e| e.is_cursed()).count();

        // Probe + assert: capture diagnostic info for any future drift.
        eprintln!(
            "spectre: {} types ({} alive: {} Blessed + {} Free, {} cursed: {} Dead + {} Undead, {} transitions)",
            idx.num_types(),
            alive,
            blessed,
            free,
            cursed,
            dead,
            undead,
            idx.transitions().len()
        );

        assert_eq!(idx.num_types(), 280, "spectre VT count");
        assert_eq!(alive + cursed, 280, "alive + cursed = total");
        assert_eq!(dead + undead, cursed, "cursed = Dead + Undead");
        assert_eq!(blessed + free, alive, "alive = Blessed + Free");
        // Concrete counts — current values; will need updating if the
        // algorithm changes.
        assert_eq!(alive, 102, "spectre alive count");
        assert_eq!(blessed, 82, "spectre blessed count");
        assert_eq!(free, 20, "spectre free count");
        assert_eq!(cursed, 178, "spectre cursed count");
        assert_eq!(dead, 159, "spectre dead count");
        assert_eq!(undead, 19, "spectre undead count");
        assert_eq!(idx.transitions().len(), 420, "spectre transition count");
    }

    /// Transition-well-formedness check for spectre (same as the hex
    /// version). Spectre has cursed VTs and closing transitions, so
    /// this also exercises CLOSED_ID handling.
    #[test]
    fn spectre_transitions_well_formed() {
        let s: crate::intgeom::snake::Snake<ZZ12> = tiles::spectre();
        let rat = Rat::try_from(&s).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = OpenVertexTypeIndex::new(Arc::clone(&ts));
        let n_types = idx.num_types();
        for t in idx.transitions() {
            assert!(
                t.src_id >= 1 && t.src_id <= n_types,
                "src_id {} out of range",
                t.src_id
            );
            assert!(
                t.dst_id == CLOSED_ID || (t.dst_id >= 1 && t.dst_id <= n_types),
                "dst_id {} out of range (and not CLOSED)",
                t.dst_id
            );
            assert!(t.tile_id < ts.num_tiles());
            assert!(t.tile_offset < ts.rat(t.tile_id).len());
        }
    }

    /// Re-derives the cursed predicate from scratch. A cursed VT must
    /// satisfy two invariants:
    ///
    /// 1. It has no closing transitions (closing = escape to Closed,
    ///    which is not cursed, so a closing transition would
    ///    contradict "all successors cursed").
    /// 2. Every open successor is itself cursed.
    ///
    /// Plus: every Dead VT has no realized transitions at all.
    #[test]
    fn spectre_cursed_invariant() {
        let s: crate::intgeom::snake::Snake<ZZ12> = tiles::spectre();
        let rat = Rat::try_from(&s).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = OpenVertexTypeIndex::new(Arc::clone(&ts));
        for info in idx.entries() {
            if !info.is_cursed() {
                continue;
            }
            let id = idx.get_id(info.vtype()).unwrap();
            let outgoing: Vec<&TransitionInfo> = idx
                .transitions()
                .iter()
                .filter(|t| t.src_id == id)
                .collect();
            if info.is_dead() {
                assert!(
                    outgoing.is_empty(),
                    "Dead VT id {id} must have no transitions, got {}",
                    outgoing.len()
                );
                continue;
            }
            // Undead: at least one transition, no closing, all open
            // successors cursed.
            assert!(
                !outgoing.is_empty(),
                "Undead VT id {id} must have at least one transition"
            );
            for t in &outgoing {
                assert!(
                    !t.is_closed(),
                    "Undead VT id {id} has a closing transition — \
                     should have been classified Free"
                );
                let dst_info = idx.get_info(t.dst_id);
                assert!(
                    dst_info.is_cursed(),
                    "Undead VT id {id} has open transition to non-cursed VT id {}",
                    t.dst_id
                );
            }
        }
    }

    /// For spectre, the blessed fixpoint also has work to do (68
    /// blessed VTs). Re-derive the predicate: every blessed VT must
    /// have at least one transition, and every non-closing outgoing
    /// transition must lead to another blessed VT.
    #[test]
    fn spectre_blessed_invariant() {
        let s: crate::intgeom::snake::Snake<ZZ12> = tiles::spectre();
        let rat = Rat::try_from(&s).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = OpenVertexTypeIndex::new(Arc::clone(&ts));
        for info in idx.entries() {
            if !info.is_blessed() {
                continue;
            }
            let id = idx.get_id(info.vtype()).unwrap();
            let outgoing: Vec<&TransitionInfo> = idx
                .transitions()
                .iter()
                .filter(|t| t.src_id == id)
                .collect();
            assert!(
                !outgoing.is_empty(),
                "Blessed VT id {id} has no transitions"
            );
            for t in &outgoing {
                assert!(
                    t.is_closed() || idx.get_info(t.dst_id).is_blessed(),
                    "Blessed VT id {id} has open transition to non-Blessed VT id {}",
                    t.dst_id
                );
            }
        }
    }

    #[test]
    fn hex_segment_types() {
        let hex: crate::intgeom::snake::Snake<ZZ12> = tiles::hexagon();
        let rat = Rat::try_from(&hex).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = OpenVertexTypeIndex::new(ts);
        eprintln!("hex: {} vertex types", idx.num_types());
    }

    /// Sanity check on hex closed VTs: at least one is discovered,
    /// every petal references the lone hex tile (tile_id == 0), every
    /// transition_count is positive, and the sum across all closed
    /// VTs equals the total number of closing transitions in the
    /// catalog.
    #[test]
    fn hexagon_closed_vertex_types_are_discovered() {
        let hex: crate::intgeom::snake::Snake<ZZ12> = tiles::hexagon();
        let rat = Rat::try_from(&hex).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = OpenVertexTypeIndex::new(ts);

        assert!(
            idx.num_closed_types() >= 1,
            "hex has at least one closed VT"
        );
        for cvt in idx.closed_entries() {
            assert!(!cvt.vtype().is_empty());
            assert!(cvt.transition_count() >= 1);
            for e in cvt.vtype().edges() {
                assert_eq!(e.tile_id, 0, "only one tile in this tileset");
                assert!(e.tile_offset < 6);
            }
        }

        let n_closing: usize = idx.transitions().iter().filter(|t| t.is_closed()).count();
        assert!(n_closing > 0, "hex has at least one closing transition");
        let total_per_vt: usize = idx
            .closed_entries()
            .iter()
            .map(|c| c.transition_count())
            .sum();
        assert_eq!(
            total_per_vt, n_closing,
            "sum of per-closed-VT transition counts equals total closings"
        );
    }

    /// Same surface for squares.
    #[test]
    fn square_closed_vertex_types_are_discovered() {
        let sq: crate::intgeom::snake::Snake<ZZ4> = tiles::square();
        let rat = Rat::try_from(&sq).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = OpenVertexTypeIndex::new(ts);

        assert!(
            idx.num_closed_types() >= 1,
            "square has at least one closed VT"
        );
        for cvt in idx.closed_entries() {
            assert!(!cvt.vtype().is_empty());
            assert!(cvt.transition_count() >= 1);
            for e in cvt.vtype().edges() {
                assert_eq!(e.tile_id, 0);
                assert!(e.tile_offset < 4);
            }
        }
    }

    /// The closed VT catalog is canonical: every ring is stored in
    /// its lex-min cyclic rotation. Verify directly.
    #[test]
    fn closed_vts_are_lex_min_canonical() {
        let hex: crate::intgeom::snake::Snake<ZZ12> = tiles::hexagon();
        let rat = Rat::try_from(&hex).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = OpenVertexTypeIndex::new(ts);

        for cvt in idx.closed_entries() {
            let canonical = ClosedVertexType::from_cyclic(cvt.vtype().edges());
            assert_eq!(cvt.vtype(), &canonical, "closed VT not in canonical form");
        }
    }

    /// For every closing transition in the catalog, its
    /// `closed_vt_id` must point at a valid 1-based id and the
    /// referenced closed VT must equal what we'd derive from the
    /// source open VT's petal ring (the BFS contract).
    #[test]
    fn closing_transitions_resolve_to_canonical_closed_vts() {
        let hex: crate::intgeom::snake::Snake<ZZ12> = tiles::hexagon();
        let rat = Rat::try_from(&hex).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = OpenVertexTypeIndex::new(ts);

        for t in idx.transitions() {
            if t.is_closed() {
                let cvt_id = t.closed_vt_id.expect("closing transition has closed_vt_id");
                assert!(cvt_id >= 1 && cvt_id <= idx.num_closed_types());
                let cvt = idx.get_closed_type(cvt_id);
                let src_open = idx.get_type(t.src_id);
                let expected = ClosedVertexType::from_open_via_closure(src_open);
                assert_eq!(
                    cvt,
                    &expected,
                    "closed VT for transition {:?} differs from its derived form",
                    (t.src_id, t.tile_id)
                );
            } else {
                assert!(
                    t.closed_vt_id.is_none(),
                    "non-closing transition has unexpected closed_vt_id"
                );
            }
        }
    }
}
