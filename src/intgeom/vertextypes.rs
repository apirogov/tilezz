use std::collections::{BTreeSet, HashMap, VecDeque};
use std::sync::Arc;

use crate::cyclotomic::{IsComplex, IsRingOrField, Units};
use crate::intgeom::patch::{GrowingPatch, OpenVertexType, PatchMatch, TransitionSide};
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
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
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

pub struct OpenVertexTypeInfo<T: IsComplex> {
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

impl<T: IsComplex> OpenVertexTypeInfo<T> {
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
}

impl TransitionInfo {
    pub fn is_closed(&self) -> bool {
        self.dst_id == CLOSED_ID
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
pub struct OpenVertexTypeIndex<T: IsComplex> {
    tileset: Arc<TileSet<T>>,
    entries: Vec<OpenVertexTypeInfo<T>>,
    transitions: Vec<TransitionInfo>,
    reverse: HashMap<OpenVertexType, usize>,
}

impl<T: IsComplex + IsRingOrField + Units> OpenVertexTypeIndex<T> {
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
        let (succ_sets, transition_infos) =
            build_transition_arrays(&reverse, &raw_transitions);
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

        OpenVertexTypeIndex {
            tileset,
            entries: info_entries,
            transitions: transition_infos,
            reverse,
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
}

/// A single raw transition tuple, accumulated by the BFS before id
/// assignment turns it into a [`TransitionInfo`].
///
/// Layout: `(src_vt, dst_vt_or_None_for_closing, side, tile_id, tile_offset)`.
type RawTransition = (
    OpenVertexType,
    Option<OpenVertexType>,
    TransitionSide,
    usize,
    usize,
);

/// Mutable state shared across the seed and BFS phases of
/// [`OpenVertexTypeIndex::new`].
struct BfsState<T: IsComplex> {
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
    witness_store: HashMap<OpenVertexType, (GrowingPatch<T>, usize, i8)>,
}

impl<T: IsComplex> Default for BfsState<T> {
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
fn seed_phase<T: IsComplex + IsRingOrField + Units>(
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
fn bfs_phase<T: IsComplex + IsRingOrField + Units>(
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
            // [pm.start_a, pm.start_a + pm.len). The focus is sealed
            // iff BOTH incident edges are in that range. (Note: a
            // vertex-range test would be too permissive — a single-
            // edge match on the CW edge has vertex range [pos-1, pos]
            // touching both pos-1 and pos without consuming the CCW
            // edge of pos.)
            let edge_in_match = |edge: usize| -> bool { (edge + n - pm.start_a) % n < pm.len };
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
                (edge_pos as i64 - pm.start_a as i64).rem_euclid(n as i64) as usize;
            let m = tileset.rat(pm.tile_id).len();
            let tile_offset = (pm.start_b as i64 + pm.len as i64 - offset_in_match as i64)
                .rem_euclid(m as i64) as usize;

            // Where is the focus's new junction (if any) in the
            // post-glue boundary?
            //   - side == Ccw: pm.start_a == pos; new boundary's
            //     index 0 is at old vertex (pos + len) % n; new
            //     junction at the CW endpoint of the match lives
            //     at new index seg_len_old = n - pm.len.
            //   - side == Cw: pm ended at vertex pos (start_a + len
            //     == pos); new boundary's index 0 is at old vertex
            //     pos; new junction at the CCW endpoint of the
            //     match lives at new index 0.
            //   - side == Both: focus vertex is sealed; no junction
            //     position.
            let junction_pos = if pm.start_a == pos { n - pm.len } else { 0 };

            if matches!(side, TransitionSide::Both) {
                state
                    .raw_transitions
                    .insert((vt.clone(), None, side, pm.tile_id, tile_offset));
            } else if let Some(new_vt) = gp2.junction_vertex_type_at(junction_pos) {
                state.raw_transitions.insert((
                    vt.clone(),
                    Some(new_vt),
                    side,
                    pm.tile_id,
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
                        state.witness_store.entry(nv.clone()).or_insert_with(|| {
                            (gp2.clone(), new_pos, gp2.angles()[new_pos])
                        });
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

/// Phase 4: turn the raw `(src, dst, side, tile_id, tile_offset)`
/// tuples into `TransitionInfo` records using the id map; build the
/// successor sets used by the cursed / blessed fixpoints.
fn build_transition_arrays(
    reverse: &HashMap<OpenVertexType, usize>,
    raw_transitions: &BTreeSet<RawTransition>,
) -> (Vec<BTreeSet<usize>>, Vec<TransitionInfo>) {
    let n = reverse.len();
    let mut succ_sets: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); n];
    let mut transition_infos: Vec<TransitionInfo> = Vec::new();

    for (src, dst, side, tid, toff) in raw_transitions {
        let Some(&src_id) = reverse.get(src) else {
            continue;
        };
        match dst {
            Some(dst_vt) => {
                if let Some(&dst_id) = reverse.get(dst_vt) {
                    succ_sets[src_id - 1].insert(dst_id);
                    transition_infos.push(TransitionInfo {
                        src_id,
                        dst_id,
                        side: *side,
                        tile_id: *tid,
                        tile_offset: *toff,
                    });
                }
            }
            None => {
                transition_infos.push(TransitionInfo {
                    src_id,
                    dst_id: CLOSED_ID,
                    side: *side,
                    tile_id: *tid,
                    tile_offset: *toff,
                });
            }
        }
    }

    (succ_sets, transition_infos)
}

/// Phase 5: assemble per-VT [`OpenVertexTypeInfo`] records, classifying
/// each as Dead / Undead / Blessed / Free from the fixpoint results.
fn classify_and_finalize<T: IsComplex + IsRingOrField + Units>(
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic::{ZZ12, ZZ4};
    use crate::intgeom::tiles;
    use crate::intgeom::tileset::TileSet;
    use std::sync::Arc;

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
            let outgoing: Vec<&TransitionInfo> =
                idx.transitions().iter().filter(|t| t.src_id == id).collect();
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
            let outgoing: Vec<&TransitionInfo> =
                idx.transitions().iter().filter(|t| t.src_id == id).collect();
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
        T: IsComplex + IsRingOrField + Units,
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
        T: IsComplex + IsRingOrField + Units,
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
                let pm = PatchMatch {
                    tile_id,
                    start_a: ext_start,
                    start_b: ext_end,
                    len,
                };
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
            // filter via junction_gap_nonnegative + try_glue_precomputed.
            let mut filtered_brute: std::collections::BTreeSet<
                (usize, usize, usize, usize),
            > = std::collections::BTreeSet::new();
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
                        let in_range =
                            |edge: usize| (edge + n - ns_u) % n < len;
                        if !in_range((pos + n - 1) % n) && !in_range(pos) {
                            continue;
                        }
                        if !crate::intgeom::matchtypes::junction_gap_nonnegative(
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

            let api_set: std::collections::BTreeSet<(usize, usize, usize, usize)> =
                witness
                    .get_matches_touching_vertex(pos)
                    .into_iter()
                    .map(|pm| (pm.tile_id, pm.start_a, pm.len, pm.start_b))
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
        // - junction_gap_nonnegative filter;
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
                    if !crate::intgeom::matchtypes::junction_gap_nonnegative(
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
            .map(|pm| (pm.tile_id, pm.start_a, pm.len, pm.start_b))
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
        // `is_single_edge_candidate(clean_tile, ia, clean_tile, ib)` —
        // which uses CLEAN-tile angles — but pass
        // `junction_gap_nonnegative(boundary, ...)` because the
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
            let bnd_right = bnd_angles[next_sa] as i32
                + tile_seq[(sb + m - len) % m] as i32;
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
        eprintln!("  edges[25]: tile_id={} tile_offset={}", edges[25].tile_id, edges[25].tile_offset);
        for i in (n - 3)..=(n - 1) {
            eprintln!(
                "    edges[{i}]: tile_id={} tile_offset={}",
                edges[i].tile_id, edges[i].tile_offset
            );
        }
        for i in 0..3 {
            eprintln!(
                "    edges[{i}]: tile_id={} tile_offset={}",
                edges[i].tile_id, edges[i].tile_offset
            );
        }

        // Compare: what does compute_all_candidates produce at result[25]?
        let result = crate::intgeom::patch::GrowingPatch::<ZZ12>::compute_all_candidates(
            &Arc::new(crate::intgeom::matchtypes::MatchTypeIndex::new(Arc::clone(&ts))),
            witness.angles(),
            witness.edges(),
        );
        eprintln!("  result[25] = {} entries:", result[25].len());
        for pm in &result[25] {
            eprintln!(
                "    cac: tile={} start_a={} len={} start_b={}",
                pm.tile_id, pm.start_a, pm.len, pm.start_b
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
            let tile_b = ts.rat(c.tile_b);
            let (ns, len, ne) = bnd_rat.get_match((25i64, c.start_b as i64), tile_b);
            if len == 0 {
                continue;
            }
            let ns_u = ns.rem_euclid(n as i64) as usize;
            let ne_u = ne.rem_euclid(tile_b.len() as i64) as usize;
            let gap_ok = crate::intgeom::matchtypes::junction_gap_nonnegative(
                witness.angles(),
                ns_u,
                len,
                tile_b.seq(),
                ne_u,
            );
            let key = (ns_u, len, ne_u, c.tile_b);
            let already_seen = seen.contains(&key);
            let glue_ok = bnd_rat
                .try_glue_precomputed((ns, len, ne), tile_b, true)
                .is_ok();
            let dst = (c.tile_b, ns_u, len, ne_u);
            let target_match = filtered_brute.contains(&dst) && dst.1 == 25;
            if target_match || !already_seen {
                eprintln!(
                    "    cand sb={} len={} -> bnd ({ns_u},{len},{ne_u}) gap_ok={gap_ok} \
                     seen_already={already_seen} glue_ok={glue_ok}  target={target_match}",
                    c.start_b, c.len
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
    ///   1. It has no closing transitions (closing = escape to Closed,
    ///      which is not cursed, so a closing transition would
    ///      contradict "all successors cursed").
    ///   2. Every open successor is itself cursed.
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
            let outgoing: Vec<&TransitionInfo> =
                idx.transitions().iter().filter(|t| t.src_id == id).collect();
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
            let outgoing: Vec<&TransitionInfo> =
                idx.transitions().iter().filter(|t| t.src_id == id).collect();
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
}
