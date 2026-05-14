use std::collections::{BTreeSet, HashMap, VecDeque};
use std::sync::Arc;

use crate::cyclotomic::{IsComplex, IsRingOrField, Units};
use crate::intgeom::patch::{GrowingPatch, OpenVertexType, PatchMatch, TransitionSide};
use crate::intgeom::rat::Rat;
use crate::intgeom::tileset::TileSet;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
enum HasTransitions {
    Yes,
    No,
}

/// Reachability classification of an open vertex type within the
/// VT transition graph.
///
/// Computed by [`OpenVertexTypeIndex::new`] via two monotone fixpoint
/// passes (one for dead/undead, one for blessed) over the transition
/// graph.
///
/// The four kinds partition the VT space:
///
/// * [`VTypeKind::Dead`] — no outgoing transition exists. Every glue
///   that touches this vertex either fails (geometric collision, bad
///   angle) or no such glue exists at all. This is the base case for
///   the cursed fixpoint.
/// * [`VTypeKind::Undead`] — reachable but trapped. Has at least one
///   outgoing transition, but every successor — transitively — is
///   cursed (Dead or Undead). No glue path from here leads to a Free
///   or Blessed state.
/// * [`VTypeKind::Blessed`] — every outgoing transition either closes
///   the junction (`dst_id == CLOSED_ID`) or leads to another Blessed
///   VT (fixpoint). Geometrically: every continuation from this VT
///   eventually seals the vertex into the interior of a patch.
/// * [`VTypeKind::Free`] — alive (not cursed) and not blessed. Has at
///   least one continuation that's still in play — neither known to
///   dead-end nor known to inevitably close.
///
/// "Cursed" is shorthand for "Dead or Undead"; "alive" for "Free or
/// Blessed". See [`OpenVertexTypeInfo::is_cursed`] / `is_alive`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum VTypeKind {
    /// No outgoing transitions — every glue touching this vertex
    /// fails or no candidate exists.
    Dead,
    /// Reachable but trapped: every successor is itself cursed.
    Undead,
    /// Every continuation either closes the vertex or stays Blessed.
    Blessed,
    /// Alive and not Blessed — at least one continuation is still
    /// in play.
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
    #[allow(dead_code)]
    has_transitions: HasTransitions,
    #[allow(dead_code)]
    successors: Vec<usize>,
    #[allow(dead_code)]
    predecessors: Vec<usize>,
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
    pub fn vtype(&self) -> &OpenVertexType {
        &self.vtype
    }

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

    pub fn is_cursed(&self) -> bool {
        matches!(self.kind, VTypeKind::Dead | VTypeKind::Undead)
    }

    pub fn is_alive(&self) -> bool {
        matches!(self.kind, VTypeKind::Blessed | VTypeKind::Free)
    }

    pub fn is_initial(&self) -> bool {
        self.is_initial
    }

    pub fn realizing_rat(&self) -> &Rat<T> {
        &self.realizing_rat
    }

    pub fn gap_angle(&self) -> i8 {
        self.gap_angle
    }

    pub fn witness(&self) -> &GrowingPatch<T> {
        &self.witness
    }

    pub fn witness_pos(&self) -> usize {
        self.witness_pos
    }

    pub fn cw_neighbor_offset(&self) -> usize {
        self.cw_neighbor_offset
    }

    pub fn ccw_neighbor_offset(&self) -> usize {
        self.ccw_neighbor_offset
    }
}

pub const CLOSED_ID: usize = 0;

pub struct TransitionInfo {
    pub src_id: usize,
    pub dst_id: usize,
    pub side: TransitionSide,
    pub tile_id: usize,
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
        let mut all_types: BTreeSet<OpenVertexType> = BTreeSet::new();
        let mut initial_types: BTreeSet<OpenVertexType> = BTreeSet::new();
        let mut transition_map: HashMap<OpenVertexType, HasTransitions> = HashMap::new();
        let mut raw_transitions: BTreeSet<(
            OpenVertexType,
            Option<OpenVertexType>,
            TransitionSide,
            usize,
            usize,
        )> = BTreeSet::new();

        let mut visited: BTreeSet<OpenVertexType> = BTreeSet::new();
        let mut queue: VecDeque<OpenVertexType> = VecDeque::new();
        let mut witness_store: HashMap<OpenVertexType, (GrowingPatch<T>, usize, i8)> =
            HashMap::new();

        for seed_id in 0..tileset.num_tiles() {
            eprintln!(
                "[VT BFS] seeding tile {}/{}...",
                seed_id,
                tileset.num_tiles()
            );
            let seed = GrowingPatch::new(Arc::clone(&tileset), seed_id);
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
                        all_types.insert(vt.clone());
                        initial_types.insert(vt.clone());
                        witness_store
                            .entry(vt.clone())
                            .or_insert_with(|| (gp.clone(), pos, gp.angles()[pos]));
                        if visited.insert(vt.clone()) {
                            queue.push_back(vt);
                        }
                    }
                }
            }
        }

        let mut bfs_processed: usize = 0;
        let bfs_start = std::time::Instant::now();

        while let Some(vt) = queue.pop_front() {
            if transition_map.get(&vt) == Some(&HasTransitions::No) {
                continue;
            }

            bfs_processed += 1;
            if bfs_processed.is_multiple_of(100) || queue.is_empty() {
                eprintln!(
                    "[VT BFS] processed {} | queue {} | types {} | {:.1?}",
                    bfs_processed,
                    queue.len(),
                    visited.len(),
                    bfs_start.elapsed(),
                );
            }

            let touching = {
                let (gp, pos, _gap) = &witness_store[&vt];
                let n = gp.boundary_len();
                let t: Vec<PatchMatch> = gp.get_matches_touching_vertex(*pos);
                if t.is_empty() {
                    transition_map.insert(vt, HasTransitions::No);
                    continue;
                }
                transition_map
                    .entry(vt.clone())
                    .or_insert(HasTransitions::Yes);
                (t, n)
            };
            let (touching, n) = touching;
            let pos = witness_store[&vt].1;

            for pm in touching {
                let gp = &witness_store[&vt].0;
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
                    raw_transitions.insert((vt.clone(), None, side, pm.tile_id, tile_offset));
                } else if let Some(new_vt) = gp2.junction_vertex_type_at(junction_pos) {
                    raw_transitions.insert((
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
                        all_types.insert(nv.clone());
                        if visited.insert(nv.clone()) {
                            witness_store
                                .insert(nv.clone(), (gp2.clone(), new_pos, gp2.angles()[new_pos]));
                            queue.push_back(nv);
                        }
                    }
                }
            }
        }

        let entries: Vec<OpenVertexType> = all_types.into_iter().collect();
        let reverse: HashMap<OpenVertexType, usize> = entries
            .iter()
            .enumerate()
            .map(|(i, vt)| (vt.clone(), i + 1))
            .collect();

        let mut succ_sets: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); entries.len()];
        let mut pred_sets: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); entries.len()];

        let mut transition_infos: Vec<TransitionInfo> = Vec::new();
        for (src, dst, side, tid, toff) in &raw_transitions {
            if let Some(&src_id) = reverse.get(src) {
                match dst {
                    Some(dst_vt) => {
                        if let Some(&dst_id) = reverse.get(dst_vt) {
                            succ_sets[src_id - 1].insert(dst_id);
                            pred_sets[dst_id - 1].insert(src_id);
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
        }

        let is_cursed = compute_cursed(&entries, &transition_map, &succ_sets);
        let is_blessed = compute_blessed(&entries, &transition_infos);

        let info_entries: Vec<OpenVertexTypeInfo<T>> = entries
            .into_iter()
            .enumerate()
            .map(|(i, vt)| {
                let id = i + 1;
                let (witness, witness_pos, gap_angle) = witness_store.remove(&vt).unwrap();
                let rat = witness.to_rat();
                let (cw_nbr, ccw_nbr) = witness
                    .neighbor_junction_offsets(witness_pos)
                    .unwrap_or((0, 0));
                let no_transitions = transition_map.get(&vt) == Some(&HasTransitions::No);
                let kind = if no_transitions {
                    VTypeKind::Dead
                } else if is_cursed[&id] {
                    VTypeKind::Undead
                } else if is_blessed[&id] {
                    VTypeKind::Blessed
                } else {
                    VTypeKind::Free
                };
                OpenVertexTypeInfo {
                    has_transitions: if no_transitions {
                        HasTransitions::No
                    } else {
                        HasTransitions::Yes
                    },
                    successors: succ_sets[i].iter().copied().collect(),
                    predecessors: pred_sets[i].iter().copied().collect(),
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
            .collect();

        OpenVertexTypeIndex {
            tileset,
            entries: info_entries,
            transitions: transition_infos,
            reverse,
        }
    }

    pub fn num_types(&self) -> usize {
        self.entries.len()
    }

    pub fn transitions(&self) -> &[TransitionInfo] {
        &self.transitions
    }

    pub fn entries(&self) -> &[OpenVertexTypeInfo<T>] {
        &self.entries
    }

    pub fn range_by_cw(&self, tile_id: usize) -> &[OpenVertexTypeInfo<T>] {
        let start = self
            .entries
            .partition_point(|e| e.vtype().cw.tile_id < tile_id);
        let end = self
            .entries
            .partition_point(|e| e.vtype().cw.tile_id <= tile_id);
        &self.entries[start..end]
    }

    pub fn get_id(&self, vtype: &OpenVertexType) -> Option<usize> {
        self.reverse.get(vtype).copied()
    }

    pub fn get_type(&self, id: usize) -> &OpenVertexType {
        assert!(
            id >= 1 && id <= self.entries.len(),
            "vertex type id out of range"
        );
        &self.entries[id - 1].vtype
    }

    pub fn get_info(&self, id: usize) -> &OpenVertexTypeInfo<T> {
        assert!(
            id >= 1 && id <= self.entries.len(),
            "vertex type id out of range"
        );
        &self.entries[id - 1]
    }

    pub fn tileset(&self) -> &Arc<TileSet<T>> {
        &self.tileset
    }
}

fn compute_cursed(
    entries: &[OpenVertexType],
    transition_map: &HashMap<OpenVertexType, HasTransitions>,
    succ_sets: &[BTreeSet<usize>],
) -> HashMap<usize, bool> {
    let n = entries.len();
    let mut cursed: HashMap<usize, bool> = HashMap::with_capacity(n);

    for (i, vt) in entries.iter().enumerate() {
        let id = i + 1;
        cursed.insert(id, transition_map.get(vt) == Some(&HasTransitions::No));
    }

    let mut changed = true;
    while changed {
        changed = false;
        for (i, _vt) in entries.iter().enumerate() {
            let id = i + 1;
            if cursed[&id] {
                continue;
            }
            let succs = &succ_sets[i];
            if succs.is_empty() {
                continue;
            }
            if succs.iter().all(|s| cursed[s]) {
                cursed.insert(id, true);
                changed = true;
            }
        }
    }

    cursed
}

fn compute_blessed(
    entries: &[OpenVertexType],
    transitions: &[TransitionInfo],
) -> HashMap<usize, bool> {
    let n = entries.len();
    let mut blessed: HashMap<usize, bool> = HashMap::with_capacity(n);

    for (i, _vt) in entries.iter().enumerate() {
        let id = i + 1;
        blessed.insert(id, false);
    }

    let mut changed = true;
    while changed {
        changed = false;
        for (i, _vt) in entries.iter().enumerate() {
            let id = i + 1;
            if blessed[&id] {
                continue;
            }
            let has_any = transitions.iter().any(|t| t.src_id == id);
            if !has_any {
                continue;
            }
            let all_closed_or_blessed = transitions
                .iter()
                .filter(|t| t.src_id == id)
                .all(|t| t.is_closed() || blessed.get(&t.dst_id).copied().unwrap_or(false));
            if all_closed_or_blessed {
                blessed.insert(id, true);
                changed = true;
            }
        }
    }

    blessed
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic::{ZZ12, ZZ4};
    use crate::intgeom::tiles;
    use crate::intgeom::tileset::TileSet;
    use std::sync::Arc;

    #[test]
    fn hexagon_vertex_type_counts() {
        let hex: crate::intgeom::snake::Snake<ZZ12> = tiles::hexagon();
        let rat = Rat::try_from(&hex).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = OpenVertexTypeIndex::new(ts);
        assert!(idx.num_types() > 0, "should discover some vertex types");
        let alive = idx.entries.iter().filter(|e| e.is_alive()).count();
        let dead = idx.entries.iter().filter(|e| e.is_dead()).count();
        let cursed = idx.entries.iter().filter(|e| e.is_cursed()).count();
        eprintln!(
            "hex: {} types ({} alive, {} dead, {} cursed)",
            idx.num_types(),
            alive,
            dead,
            cursed
        );
    }

    #[test]
    fn square_vertex_type_counts() {
        let sq: crate::intgeom::snake::Snake<ZZ4> = tiles::square();
        let rat = Rat::try_from(&sq).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = OpenVertexTypeIndex::new(ts);
        assert!(idx.num_types() > 0);
        let alive = idx.entries.iter().filter(|e| e.is_alive()).count();
        eprintln!("square: {} types ({} alive)", idx.num_types(), alive);
    }

    #[test]
    fn hexagon_witness_consistency() {
        let hex: crate::intgeom::snake::Snake<ZZ12> = tiles::hexagon();
        let rat = Rat::try_from(&hex).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = OpenVertexTypeIndex::new(ts);
        for info in &idx.entries {
            let vt = info.witness().junction_vertex_type_at(info.witness_pos());
            assert_eq!(
                vt,
                Some(info.vtype.clone()),
                "witness should realize its claimed type"
            );
        }
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

    #[test]
    fn hex_segment_types() {
        let hex: crate::intgeom::snake::Snake<ZZ12> = tiles::hexagon();
        let rat = Rat::try_from(&hex).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let idx = OpenVertexTypeIndex::new(ts);
        eprintln!("hex: {} vertex types", idx.num_types());
    }
}
