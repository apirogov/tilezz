use rustc_hash::FxHashSet;
use std::sync::Arc;

use crate::cyclotomic::geometry::intersect;
use crate::cyclotomic::{IsComplex, IsRingOrField, Units};
use crate::intgeom::angles;
use crate::intgeom::grid::UnitSquareGrid;
use crate::intgeom::matchtypes::{
    is_single_edge_candidate, junction_gap_nonnegative, CandidateMatch, MatchTypeIndex,
};
use crate::intgeom::rat::Rat;
use crate::intgeom::snake::Snake;
use crate::intgeom::tileset::TileSet;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct EdgeInfo {
    pub tile_id: usize,
    pub tile_offset: usize,
}

/// Auxiliary per-vertex info on a patch boundary.
///
/// Each boundary position has an angle and two adjacent edges (cw, ccw).
/// This is instrumental infrastructure — the interesting structure lives
/// in [`OpenVertexType`] at actual junction vertices.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct PatchVertexInfo {
    pub angle: i8,
    pub cw: EdgeInfo,
    pub ccw: EdgeInfo,
}

/// A maximal contiguous range `[patch_start, patch_end)` of boundary
/// positions whose edges all come from a single tile instance.
///
/// `tile_id` is the tile that segment belongs to; `offset_start` is the
/// tile-offset of `edges[patch_start]`. Within the segment,
/// `edges[patch_start + k].tile_offset == (offset_start + k) mod
/// tile_len`.
///
/// **Cyclic-vs-linear caveat.** The boundary is cyclic but
/// `TileSegment`s are produced from the linear array `[0, n)`. If
/// position 0 is not a junction, a single cyclic tile-instance run that
/// straddles position 0 is split into two linear `TileSegment`s — one
/// at the start `[0, first_junction)` and one at the end
/// `[last_junction, n)`. Both have the same `tile_id`, and their
/// offsets are cyclically continuous (the last segment's offsets
/// continue, modulo tile length, into the first segment's). Callers
/// that need cyclic tile instances rather than linear segments must
/// stitch these two halves themselves.
pub struct TileSegment {
    pub patch_start: usize,
    pub patch_end: usize,
    pub tile_id: usize,
    pub offset_start: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PatchMatch {
    pub start_a: usize,
    pub len: usize,
    pub start_b: usize,
    pub tile_id: usize,
}

/// An **open** junction vertex: the arrangement of tiles meeting at a
/// boundary vertex that is *not* fully surrounded.
///
/// At a junction (where the boundary angle differs from the tile's internal
/// angle), multiple tile instances converge. `cw` is the edge arriving from
/// the CW direction, `ccw` is the edge departing in the CCW direction, and
/// `inner` lists any enclosed tile edges between them.
///
/// "Open" means the vertex is on the patch boundary — there is still room
/// for additional tiles to be glued in (i.e. the cumulative angle around
/// the vertex is strictly less than a full turn). Equivalently, the boundary
/// angle at the vertex is **not** ±half-turn — a half-turn boundary angle
/// would mean the boundary doubles back on itself (a degenerate pinched
/// vertex), which is excluded by construction.
///
/// In contrast, a *closed* vertex type (not represented by this struct)
/// would describe a fully-surrounded interior vertex, where the petals
/// cover the entire turn and the vertex no longer lies on any boundary.
#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct OpenVertexType {
    pub cw: EdgeInfo,
    pub inner: Vec<EdgeInfo>,
    pub ccw: EdgeInfo,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum TransitionSide {
    Cw,
    Ccw,
}

/// An incrementally grown, edge-to-edge tiling of a connected,
/// **hole-free** region of the plane.
///
/// # Invariants
///
/// * **Hole-free** (simply connected). The patch's interior is
///   topologically a disk: every interior point is reachable from every
///   other along a path that doesn't cross the boundary, and the boundary
///   itself is a single simple closed polygon. There are no enclosed
///   "empty" regions between tiles.
///
/// * **Edge-to-edge**. Each tile edge either lies on the patch boundary
///   or coincides exactly with another tile's edge. Partial-edge overlaps
///   ("T-junctions" along an edge) are not allowed.
///
/// These invariants are upheld by [`GrowingPatch::add_tile`] (geometric
/// collision check via the spatial grid, plus the angle-matching
/// constraints) and by [`GrowingPatch::construct_witness_from_vt_sequence`]
/// (which additionally rejects ±hturn boundary junctions). Several
/// pieces of internal code rely on hole-freeness for correctness; the
/// relevant comments cite this invariant where it matters.
#[derive(Clone)]
pub struct GrowingPatch<T: IsComplex> {
    match_index: Arc<MatchTypeIndex<T>>,
    state: PatchState<T>,
}

#[derive(Clone)]
#[allow(clippy::large_enum_variant)]
enum PatchState<T: IsComplex> {
    Seed {
        tile_id: usize,
        cached_matches: Vec<PatchMatch>,
    },
    /// Active patch with a boundary.
    ///
    /// # Grid bookkeeping
    ///
    /// The spatial grid uses **stable edge IDs** so incremental glue
    /// updates only touch the changed segments:
    ///
    /// * `seg_data[id] -> (p1, p2)` for each edge ever created.
    /// * `boundary_edge_ids[pos] -> id` for the current boundary.
    /// * `next_edge_id` is the next free ID; bumped on each new edge.
    /// * `grid` indexes `seg_data` by ID.
    ///
    /// Removed (matched-away) edges leave their entry in `seg_data`
    /// rather than being shifted out — keeping live IDs stable across
    /// glues is what lets `unregister_segment` work in O(1). Consequently
    /// `seg_data.len()` grows with **total glue history**, not with
    /// current boundary length. For typical workflows (short-lived
    /// patches in BFS enumeration; small handwritten fixtures) the
    /// overhead is small. For long-lived patches grown over many glues,
    /// pass through [`GrowingPatch::from_parts`] to rebuild the grid
    /// from `angles`/`edges`/`inner_chains`/`patch_tile_ids` —
    /// `from_parts` re-traces positions, builds a fresh grid, and
    /// resets `boundary_edge_ids` to `(0..n)`, dropping all tombstones.
    Growing {
        angles: Vec<i8>,
        edges: Vec<EdgeInfo>,
        candidates_by_start: Option<Vec<Vec<PatchMatch>>>,
        inner_chains: Vec<Vec<EdgeInfo>>,
        positions: Vec<T>,
        grid: UnitSquareGrid,
        seg_data: Vec<(T, T)>,
        boundary_edge_ids: Vec<usize>,
        next_edge_id: usize,
        patch_tile_ids: Vec<usize>,
        next_tile_id: usize,
    },
}

/// Returns true if position `pos` on the boundary is a junction vertex,
/// i.e. the boundary angle differs from the tile's internal angle at the
/// outgoing edge.
pub(crate) fn is_junction_at<T: IsComplex + IsRingOrField + Units>(
    angles: &[i8],
    edges: &[EdgeInfo],
    tileset: &TileSet<T>,
    pos: usize,
) -> bool {
    let ei = edges[pos];
    tileset.rat(ei.tile_id).seq()[ei.tile_offset] != angles[pos]
}

/// Raw vertex type extraction at a boundary position (no junction check).
///
/// Constructs the (cw, inner, ccw) triple from boundary edge/inner-chain data.
/// Used only by tests; prefer [`GrowingPatch::junction_vertex_type_at`] which
/// returns `None` at non-junction positions.
#[cfg(test)]
pub(crate) fn vertex_type_raw_from(
    edges: &[EdgeInfo],
    inner_chains: &[Vec<EdgeInfo>],
    pos: usize,
) -> OpenVertexType {
    let n = edges.len();
    OpenVertexType {
        cw: edges[(pos + n - 1) % n],
        inner: inner_chains[pos].clone(),
        ccw: edges[pos],
    }
}

/// Compute the maximum match length when gluing `other` onto `self_angles`
/// starting at `self_start`, with the match endpoint at `other_junction`
/// on the other side.
///
/// The two sides are walked outward from a shared anchor edge: `self`
/// advances CCW from `self_start`, `other` advances CW from
/// `other_junction`. Edges are compatible when `self_angles[k] ==
/// -other_angles[mirror(k)]` (the angles meet head-to-tail with opposite
/// signs, since one side is traversed CCW and the other CW).
///
/// Length is at least 1 (the anchor edge itself always matches by
/// construction) and stops as soon as the angle pair disagrees, capped
/// at `min(self_len, other_len)`.
pub(crate) fn forward_match_length(
    self_angles: &[i8],
    self_start: usize,
    other_angles: &[i8],
    other_junction: usize,
) -> usize {
    let n = self_angles.len();
    let m = other_angles.len();
    let max_len = n.min(m);
    let mut len = 1;
    for i in 1..max_len {
        let xi = self_angles[(self_start + i) % n];
        let yi = -other_angles[(other_junction + m - i) % m];
        if xi != yi {
            break;
        }
        len += 1;
    }
    len
}

pub fn update_inner_chains(
    old_inner: &[Vec<EdgeInfo>],
    old_edges: &[EdgeInfo],
    pm: &PatchMatch,
    new_n: usize,
    old_ptids: &[usize],
) -> Vec<Vec<EdgeInfo>> {
    let n = old_edges.len();
    let seg_len_old = n - pm.len;
    let ccw_pos = (pm.start_a + pm.len) % n;
    let cw_end_matched = (pm.start_a + pm.len - 1) % n;

    let mut new_inner = vec![Vec::new(); new_n];

    for i in 1..seg_len_old {
        new_inner[i] = old_inner[(ccw_pos + i) % n].clone();
    }

    let consumed_cw_ptid = old_ptids[cw_end_matched];
    let ccw_survivor_ptid = old_ptids[ccw_pos];
    let mut chain_ccw: Vec<EdgeInfo> = old_inner[ccw_pos].clone();
    if consumed_cw_ptid != ccw_survivor_ptid {
        chain_ccw.push(old_edges[cw_end_matched]);
    }
    new_inner[0] = chain_ccw;

    let consumed_ccw_ptid = old_ptids[pm.start_a];
    let cw_survivor_ptid = old_ptids[(pm.start_a + n - 1) % n];
    let mut chain_cw: Vec<EdgeInfo> = old_inner[pm.start_a].clone();
    if consumed_ccw_ptid != cw_survivor_ptid {
        chain_cw.push(old_edges[pm.start_a]);
    }
    new_inner[seg_len_old] = chain_cw;

    new_inner
}

/// Raw boundary representation used during witness construction.
///
/// Like [`GrowingPatch`], a `RawBoundary` describes the outline of a
/// hole-free, edge-to-edge tiling — it carries the same per-position
/// angle / edge / inner-chain / patch-tile-id data, but without the
/// spatial-grid bookkeeping needed for incremental glue checks. The
/// hole-free invariant applies here too: code that consumes a
/// `RawBoundary` may assume the underlying region is simply connected.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct RawBoundary {
    pub angles: Vec<i8>,
    pub edges: Vec<EdgeInfo>,
    pub inner_chains: Vec<Vec<EdgeInfo>>,
    pub patch_tile_ids: Vec<usize>,
}

#[cfg(test)]
pub(crate) fn raw_is_junction<T: IsComplex + IsRingOrField + Units>(
    boundary: &RawBoundary,
    tileset: &TileSet<T>,
    pos: usize,
) -> bool {
    is_junction_at(&boundary.angles, &boundary.edges, tileset, pos)
}

#[cfg(test)]
pub(crate) fn next_junction_on_raw_boundary<T: IsComplex + IsRingOrField + Units>(
    boundary: &RawBoundary,
    tileset: &TileSet<T>,
    from_pos: usize,
) -> Option<usize> {
    let n = boundary.angles.len();
    for step in 1..=n {
        let pos = (from_pos + step) % n;
        if raw_is_junction(boundary, tileset, pos) {
            return Some(pos);
        }
    }
    None
}

/// Result of [`glue_match_to_raw_boundary`].
///
/// `new_junc_pos` is the index, in the resulting `boundary`, of the
/// junction at the CCW end of the match — i.e. the first new-tile edge.
/// Equivalently, it equals `old_survivor_len` (see the **rotation
/// convention** documented on [`glue_match_to_raw_boundary`]).
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct RawGlueResult {
    pub boundary: RawBoundary,
    pub new_junc_pos: usize,
    pub match_len: usize,
    pub old_survivor_len: usize,
    pub new_tile_survivor_len: usize,
}

/// Glue one tile to a raw boundary at a boundary junction.
///
/// The `target.tile_offset` is the tile-side junction offset used as the
/// `start_b`/`other_end` argument to `glue_raw_angles`. Surviving target edges
/// are therefore `tile_offset, tile_offset + 1, ...`, modulo the tile edge
/// count. This helper is intentionally shared by VT witness construction and NT
/// enumeration so that edge-offset conventions cannot diverge again.
pub(crate) fn glue_tile_to_raw_boundary<T: IsComplex + IsRingOrField + Units>(
    boundary: &RawBoundary,
    junc_pos: usize,
    target: EdgeInfo,
    tileset: &TileSet<T>,
    first_step: bool,
    new_tile_id: usize,
) -> Option<RawGlueResult> {
    let tile_seq = tileset.rat(target.tile_id).seq();
    let tile_junc = target.tile_offset;
    let mlen = forward_match_length(&boundary.angles, junc_pos, tile_seq, tile_junc);
    glue_match_to_raw_boundary::<T>(
        boundary,
        &PatchMatch {
            start_a: junc_pos,
            len: mlen,
            start_b: tile_junc,
            tile_id: target.tile_id,
        },
        tileset,
        first_step,
        new_tile_id,
    )
}

/// Glue one match onto a raw boundary.
///
/// # Rotation convention
///
/// The returned boundary is **rotated** so that index 0 is the first
/// surviving old edge immediately CCW of the match (i.e. position
/// `(pm.start_a + pm.len) % old_n` in the old boundary). The layout is:
///
/// ```text
/// new_edges[0 .. old_survivor_len]           = surviving old edges
/// new_edges[old_survivor_len .. new_n]       = new tile's surviving edges
/// ```
///
/// `RawGlueResult::new_junc_pos = old_survivor_len` is the index of the
/// CCW junction (first new-tile edge) in the new boundary. The CW
/// junction is at index 0 of the new boundary.
///
/// **Caller consequence.** A position `j` tracked in the *old* boundary
/// that survives the glue (i.e. lies outside the matched range) maps to
/// `(j + old_n - ccw_pos) % old_n` in the new boundary, where
/// `ccw_pos = (pm.start_a + pm.len) % old_n`. [`flower_petal_glue`]
/// applies this remap automatically for the `prior_juncs` it tracks
/// across multiple glues; bespoke callers in `neighborhood.rs` use the
/// returned `new_junc_pos` as their CCW-side anchor and derive other
/// positions relative to it, so they don't need an explicit remap.
///
/// This rotation is mandatory in some form because the boundary length
/// changes (`old_n → old_n - pm.len + new_tile_seg_len`); positions in
/// the new array necessarily differ from positions in the old one. The
/// "anchor at CCW survivor" choice keeps the formula uniform — a single
/// modular subtraction — and matches what's most natural for callers
/// that just glued at a junction and want the next junction's index.
pub(crate) fn glue_match_to_raw_boundary<T: IsComplex + IsRingOrField + Units>(
    boundary: &RawBoundary,
    pm: &PatchMatch,
    tileset: &TileSet<T>,
    first_step: bool,
    new_tile_id: usize,
) -> Option<RawGlueResult> {
    let tile_seq = tileset.rat(pm.tile_id).seq();
    let m = tile_seq.len();
    let n = boundary.angles.len();
    let mlen = pm.len;

    if mlen == 0 || mlen > n || mlen > m {
        return None;
    }

    let seg_len_old = n - mlen;
    let seg_len_new = m - mlen;
    if seg_len_new == 0 {
        return None;
    }
    let new_n = seg_len_old + seg_len_new;

    let (new_edges, new_patch_tile_ids) = build_glued_edges(
        &boundary.edges,
        &boundary.patch_tile_ids,
        pm,
        m,
        new_tile_id,
    );

    let new_inner = if first_step {
        vec![vec![]; new_n]
    } else {
        update_inner_chains(
            &boundary.inner_chains,
            &boundary.edges,
            pm,
            new_n,
            &boundary.patch_tile_ids,
        )
    };

    let gr =
        angles::glue_raw_angles::<T>(&boundary.angles, tile_seq, pm.start_a, pm.len, pm.start_b)?;

    Some(RawGlueResult {
        boundary: RawBoundary {
            angles: gr.angles,
            edges: new_edges,
            inner_chains: new_inner,
            patch_tile_ids: new_patch_tile_ids,
        },
        new_junc_pos: seg_len_old,
        match_len: mlen,
        old_survivor_len: seg_len_old,
        new_tile_survivor_len: seg_len_new,
    })
}

/// Compute the sequence of boundary angles at a junction vertex as petals
/// are stacked CCW from the CW side.
///
/// Returns a vector of length `1 + |inner| + 1` where:
///  * `result[0]` is the boundary angle with only the CW tile present —
///    that is, the internal angle of `vtype.cw`'s tile at the junction
///    vertex (the corner just CCW of `cw.tile_offset`).
///  * Each subsequent entry is the new boundary angle after adding one
///    more petal, computed as `next = petal_angle + prev - hturn`.
///    Geometrically: each petal eats away `hturn - petal_angle` from the
///    remaining boundary angle.
///  * `result.last()` is the boundary angle once all petals (inner +
///    ccw) are in place — this is what the witness boundary should show
///    at the junction position.
///
/// For convex tiles (positive internal angles), the sequence is
/// monotonically non-increasing. This is asserted by
/// `assert_junction_angle_sequence_valid` in tests.
pub fn junction_angle_sequence<T: IsComplex + IsRingOrField + Units>(
    vtype: &OpenVertexType,
    tileset: &TileSet<T>,
) -> Vec<i8> {
    let seed_seq = tileset.rat(vtype.cw.tile_id).seq();
    let n_seed = seed_seq.len();
    let junc_vertex = (vtype.cw.tile_offset + 1) % n_seed;
    let mut result = vec![seed_seq[junc_vertex]];

    let mut petals = vtype.inner.clone();
    petals.push(vtype.ccw);

    for petal in &petals {
        let petal_seq = tileset.rat(petal.tile_id).seq();
        let petal_angle = petal_seq[petal.tile_offset];
        let prev = *result.last().unwrap();
        let next = angles::normalize_angle::<T>(petal_angle + prev - T::hturn());
        result.push(next);
    }

    result
}

/// Glue a sequence of "petal" tiles onto a raw boundary at a single junction.
///
/// Each call to [`glue_tile_to_raw_boundary`] rotates the boundary per the
/// convention documented on [`glue_match_to_raw_boundary`]. Any junction
/// positions in `prior_juncs` that the caller is tracking across this
/// call get remapped via `(j + n_old - ccw_pos) % n_old` after each glue,
/// so the caller sees them at their new indices in the final boundary.
///
/// Returns `(new_boundary, junc_pos_after_last_glue)` on success.
fn flower_petal_glue<T: IsComplex + IsRingOrField + Units>(
    mut boundary: RawBoundary,
    mut junc_pos: usize,
    targets: &[EdgeInfo],
    tileset: &TileSet<T>,
    mut first_step: bool,
    next_tile_id: &mut usize,
    prior_juncs: &mut [usize],
) -> Option<(RawBoundary, usize)> {
    for target in targets {
        let tile_id = *next_tile_id;
        *next_tile_id += 1;
        let n_old = boundary.angles.len();
        let result = glue_tile_to_raw_boundary::<T>(
            &boundary, junc_pos, *target, tileset, first_step, tile_id,
        )?;

        // Remap previously tracked junction positions per the rotation
        // convention documented on `glue_match_to_raw_boundary`.
        let ccw_pos = (junc_pos + result.match_len) % n_old;
        for j in prior_juncs.iter_mut() {
            *j = (*j + n_old - ccw_pos) % n_old;
        }

        first_step = false;
        boundary = result.boundary;
        junc_pos = result.new_junc_pos;
    }

    Some((boundary, junc_pos))
}

impl<T: IsComplex + IsRingOrField + Units> GrowingPatch<T> {
    pub fn new(tileset: Arc<TileSet<T>>, seed_tile_id: usize) -> Self {
        let match_index = Arc::new(MatchTypeIndex::new(Arc::clone(&tileset)));
        let seed_matches =
            Self::compute_seed_matches(&match_index, match_index.tileset(), seed_tile_id);
        GrowingPatch {
            match_index,
            state: PatchState::Seed {
                tile_id: seed_tile_id,
                cached_matches: seed_matches,
            },
        }
    }

    /// Reconstruct a `Growing` patch from its core boundary data.
    ///
    /// Returns `None` if the input slices are inconsistent (empty, or
    /// mismatched lengths between `angles`, `edges`, `inner_chains`,
    /// `patch_tile_ids`).
    ///
    /// **Caveats — derived state is rebuilt, not restored.** This
    /// constructor recomputes the spatial grid (`grid`, `seg_data`,
    /// `boundary_edge_ids`, `next_edge_id`) from `angles` via
    /// `trace_boundary_positions` + `build_boundary_grid`. Specifically:
    /// `boundary_edge_ids` is set to `(0..n).collect()` regardless of
    /// any pre-existing identification scheme, and `candidates_by_start`
    /// is left as `None` (will be filled on demand). Use this when you
    /// have a serialized boundary and want a fresh `GrowingPatch`, not
    /// to faithfully restore a snapshot of every internal field.
    pub fn from_parts(
        match_index: Arc<MatchTypeIndex<T>>,
        angles: Vec<i8>,
        edges: Vec<EdgeInfo>,
        inner_chains: Vec<Vec<EdgeInfo>>,
        patch_tile_ids: Vec<usize>,
        next_tile_id: usize,
    ) -> Option<Self> {
        if angles.is_empty()
            || angles.len() != edges.len()
            || inner_chains.len() != angles.len()
            || patch_tile_ids.len() != angles.len()
        {
            return None;
        }
        let n = angles.len();
        let positions = trace_boundary_positions::<T>(&angles);
        let (grid, seg_data, boundary_edge_ids, next_edge_id) =
            build_boundary_grid::<T>(&positions, n);
        Some(GrowingPatch {
            match_index,
            state: PatchState::Growing {
                angles,
                edges,
                candidates_by_start: None,
                inner_chains,
                positions,
                grid,
                seg_data,
                boundary_edge_ids,
                next_edge_id,
                patch_tile_ids,
                next_tile_id,
            },
        })
    }

    /// Construct the smallest patch that realises a single
    /// [`OpenVertexType`] at one of its boundary junction positions.
    ///
    /// Returns the resulting patch and the junction position `wpos` within
    /// it. Returns `None` if the witness cannot be realised — including
    /// the case where the requested vertex would close into a degenerate
    /// (±hturn / pinched) boundary, which would violate the open-VT
    /// contract.
    pub fn construct_minimal_witness(
        vtype: &OpenVertexType,
        match_index: Arc<MatchTypeIndex<T>>,
    ) -> Option<(Self, usize)> {
        let (patch, juncs) =
            Self::construct_witness_from_vt_sequence(std::slice::from_ref(vtype), match_index)?;
        Some((patch, juncs[0]))
    }

    /// Construct a patch whose boundary realises the given sequence of
    /// [`OpenVertexType`]s as consecutive junctions, in CCW order.
    ///
    /// Returns the patch plus the boundary positions of each junction in
    /// `vtypes` order.
    ///
    /// Returns `None` if any of the requested junctions cannot be realised
    /// as an open boundary vertex — specifically, if a construction step
    /// would produce a ±hturn (pinched / degenerate) boundary angle at a
    /// tracked junction position. Since the function's contract is to
    /// deliver a chain of *open* boundary junctions, such a configuration
    /// is rejected rather than returned silently.
    pub fn construct_witness_from_vt_sequence(
        vtypes: &[OpenVertexType],
        match_index: Arc<MatchTypeIndex<T>>,
    ) -> Option<(Self, Vec<usize>)> {
        if vtypes.is_empty() {
            return None;
        }
        Self::construct_witness_from_vt_sequence_inner(vtypes, 0, &match_index)
    }

    fn construct_witness_from_vt_sequence_inner(
        vtypes: &[OpenVertexType],
        start: usize,
        match_index: &Arc<MatchTypeIndex<T>>,
    ) -> Option<(Self, Vec<usize>)> {
        let tileset = match_index.tileset();
        let n_vts = vtypes.len();
        let seed_id = vtypes[start].cw.tile_id;
        let seed_seq = tileset.rat(seed_id).seq();
        let n_seed = seed_seq.len();

        let mut raw = RawBoundary {
            angles: seed_seq.to_vec(),
            edges: (0..n_seed)
                .map(|i| EdgeInfo {
                    tile_id: seed_id,
                    tile_offset: i,
                })
                .collect(),
            inner_chains: vec![vec![]; n_seed],
            patch_tile_ids: vec![0; n_seed],
        };
        let mut next_tile_id = 1usize;
        let mut junc_positions = Vec::new();

        let vt0 = &vtypes[start];
        let junc_pos = (vt0.cw.tile_offset + 1) % n_seed;
        let mut all_targets = vt0.inner.clone();
        all_targets.push(vt0.ccw);

        let (new_raw, new_junc) = flower_petal_glue::<T>(
            raw,
            junc_pos,
            &all_targets,
            tileset.as_ref(),
            true,
            &mut next_tile_id,
            &mut junc_positions,
        )?;
        // Open-VT invariant: the reconstructed junction must not be ±hturn.
        // A ±hturn boundary angle means the boundary doubles back at this
        // vertex (the vertex is closed/degenerate), which contradicts the
        // input being a sequence of open boundary junctions.
        if new_raw.angles[new_junc].abs() == T::hturn() {
            return None;
        }
        raw = new_raw;
        junc_positions.push(new_junc);

        for step in 1..n_vts {
            let k_prev = (start + step - 1) % n_vts;
            let k = (start + step) % n_vts;
            let vt = &vtypes[k];
            let vt_prev = &vtypes[k_prev];

            let tile_size = tileset.rat(vt_prev.ccw.tile_id).seq().len();
            let offset = ((vt.cw.tile_offset as isize - vt_prev.ccw.tile_offset as isize
                + tile_size as isize) as usize)
                % tile_size
                + 1;
            let prev_junc = junc_positions[step - 1];
            let n = raw.edges.len();
            let junc_k = (prev_junc + offset) % n;

            let angle_seq = junction_angle_sequence(vt, tileset.as_ref());
            let mut all_targets = vt.inner.clone();
            all_targets.push(vt.ccw);

            let boundary_angle = raw.angles[junc_k];
            let skip = angle_seq
                .iter()
                .position(|&a| a == boundary_angle)
                .unwrap_or(0);
            let targets = &all_targets[skip..];

            if targets.is_empty() {
                junc_positions.push(junc_k);
                continue;
            }

            let (new_raw, new_junc) = flower_petal_glue::<T>(
                raw,
                junc_k,
                targets,
                tileset.as_ref(),
                false,
                &mut next_tile_id,
                &mut junc_positions,
            )?;
            // Open-VT invariant: see comment on the initial glue above.
            if new_raw.angles[new_junc].abs() == T::hturn() {
                return None;
            }
            raw = new_raw;
            junc_positions.push(new_junc);
        }

        let patch = GrowingPatch::from_parts(
            Arc::clone(match_index),
            raw.angles,
            raw.edges,
            raw.inner_chains,
            raw.patch_tile_ids,
            next_tile_id,
        )?;

        Some((patch, junc_positions))
    }

    pub fn compute_candidates_covering_position(
        match_index: &Arc<MatchTypeIndex<T>>,
        angles: &[i8],
        edges: &[EdgeInfo],
        target: usize,
    ) -> Vec<PatchMatch> {
        let n = angles.len();
        let rat = Rat::from_slice_unchecked(angles);
        let tileset = match_index.tileset();
        let max_tile_len = tileset.rats().iter().map(|r| r.len()).max().unwrap_or(0);
        let mut global_seen: FxHashSet<(usize, usize, usize, usize)> = FxHashSet::default();
        let mut result: Vec<PatchMatch> = Vec::new();

        let segments = compute_segments(angles, edges, tileset);
        let junctions_set: FxHashSet<usize> = compute_junctions(angles, edges, tileset)
            .into_iter()
            .collect();

        for offset in 0..max_tile_len.min(n) {
            let pos = (target + n - offset) % n;

            if let Some(segment) = segments
                .iter()
                .find(|s| pos >= s.patch_start && pos < s.patch_end)
            {
                let tile_id = segment.tile_id;
                let tile_offset = (segment.offset_start + (pos - segment.patch_start))
                    % tileset.rat(tile_id).len();
                let mut tmp_by_start: Vec<Vec<PatchMatch>> = vec![Vec::new(); n];
                compute_candidates_at_position(
                    pos,
                    angles,
                    &rat,
                    n,
                    tile_id,
                    tile_offset,
                    match_index,
                    tileset,
                    &mut global_seen,
                    &mut tmp_by_start,
                );
                for pm in tmp_by_start.into_iter().flatten() {
                    if cyclic_range_contains(pm.start_a, pm.len, target, n) {
                        result.push(pm);
                    }
                }
            }

            if junctions_set.contains(&pos) {
                enumerate_junction_candidates_at(
                    pos,
                    angles,
                    &rat,
                    n,
                    tileset,
                    &mut global_seen,
                    |start_a, len| cyclic_range_contains(start_a, len, target, n),
                    |pm| result.push(pm),
                );
            }
        }

        result
    }

    pub fn compute_all_candidates(
        match_index: &Arc<MatchTypeIndex<T>>,
        angles: &[i8],
        edges: &[EdgeInfo],
    ) -> Vec<Vec<PatchMatch>> {
        let n = angles.len();
        let rat = Rat::from_slice_unchecked(angles);
        let tileset = match_index.tileset();
        let mut result = vec![Vec::new(); n];

        let segments = compute_segments(angles, edges, tileset);
        let junctions = compute_junctions(angles, edges, tileset);

        let mut global_seen: FxHashSet<(usize, usize, usize, usize)> = FxHashSet::default();

        for segment in &segments {
            let tile_id = segment.tile_id;
            let seg_len = segment.patch_end - segment.patch_start;
            for local_k in 0..seg_len {
                let tile_offset = (segment.offset_start + local_k) % tileset.rat(tile_id).len();
                let patch_pos = segment.patch_start + local_k;
                compute_candidates_at_position(
                    patch_pos,
                    angles,
                    &rat,
                    n,
                    tile_id,
                    tile_offset,
                    match_index,
                    tileset,
                    &mut global_seen,
                    &mut result,
                );
            }
        }

        for &junc_idx in &junctions {
            enumerate_junction_candidates_at(
                junc_idx,
                angles,
                &rat,
                n,
                tileset,
                &mut global_seen,
                |_, _| true,
                |pm| result[pm.start_a].push(pm),
            );
        }

        result
    }

    pub fn is_growing(&self) -> bool {
        matches!(self.state, PatchState::Growing { .. })
    }

    pub fn get_all_matches(&self) -> Vec<PatchMatch> {
        match &self.state {
            PatchState::Seed { cached_matches, .. } => cached_matches.clone(),
            PatchState::Growing {
                angles,
                edges,
                candidates_by_start,
                ..
            } => match candidates_by_start {
                Some(cbs) => cbs.iter().flatten().cloned().collect(),
                None => Self::compute_all_candidates(&self.match_index, angles, edges)
                    .into_iter()
                    .flatten()
                    .collect(),
            },
        }
    }

    pub fn get_matches_for_tile(&self, tile_id: usize) -> impl Iterator<Item = PatchMatch> + '_ {
        self.get_all_matches()
            .into_iter()
            .filter(move |m| m.tile_id == tile_id)
    }

    pub fn get_matches_touching_vertex(&self, vertex_index: usize) -> Vec<PatchMatch> {
        match &self.state {
            PatchState::Seed { cached_matches, .. } => {
                let n = self.match_index.tileset().rat(0).len();
                cached_matches
                    .iter()
                    .filter(|pm| cyclic_range_contains(pm.start_a, pm.len, vertex_index, n))
                    .cloned()
                    .collect()
            }
            PatchState::Growing {
                candidates_by_start,
                angles,
                edges,
                ..
            } => {
                let n = angles.len();
                let computed: Vec<Vec<PatchMatch>>;
                let source: &[Vec<PatchMatch>] = match candidates_by_start {
                    Some(cbs) => cbs,
                    None => {
                        computed = Self::compute_all_candidates(&self.match_index, angles, edges);
                        &computed
                    }
                };
                let k = self
                    .match_index
                    .tileset()
                    .rats()
                    .iter()
                    .map(|r| r.len())
                    .max()
                    .unwrap_or(0)
                    .min(n);
                let mut result = Vec::new();
                for offset in 0..=k {
                    let start = (vertex_index + n - offset) % n;
                    for pm in &source[start] {
                        if cyclic_range_contains(pm.start_a, pm.len, vertex_index, n) {
                            result.push(pm.clone());
                        }
                    }
                }
                result
            }
        }
    }

    pub fn ensure_candidates_materialized(&mut self) {
        if let PatchState::Growing {
            angles,
            edges,
            candidates_by_start,
            ..
        } = &mut self.state
        {
            if candidates_by_start.is_none() {
                *candidates_by_start = Some(Self::compute_all_candidates(
                    &self.match_index,
                    angles,
                    edges,
                ));
            }
        }
    }

    pub fn num_tiles(&self) -> usize {
        self.match_index.tileset().num_tiles()
    }

    pub fn boundary_len(&self) -> usize {
        match &self.state {
            PatchState::Growing { angles, .. } => angles.len(),
            _ => 0,
        }
    }

    pub fn angles(&self) -> &[i8] {
        match &self.state {
            PatchState::Growing { angles, .. } => angles,
            _ => &[],
        }
    }

    pub fn to_rat(&self) -> Rat<T> {
        match &self.state {
            PatchState::Seed { tile_id, .. } => self.match_index.tileset().rat(*tile_id).clone(),
            PatchState::Growing { angles, .. } => Rat::from_slice_unchecked(angles),
        }
    }

    pub fn tileset(&self) -> &Arc<TileSet<T>> {
        self.match_index.tileset()
    }

    pub fn match_index(&self) -> &Arc<MatchTypeIndex<T>> {
        &self.match_index
    }

    pub fn edges(&self) -> &[EdgeInfo] {
        match &self.state {
            PatchState::Growing { edges, .. } => edges,
            _ => &[],
        }
    }

    pub fn inner_chains(&self) -> &[Vec<EdgeInfo>] {
        match &self.state {
            PatchState::Growing { inner_chains, .. } => inner_chains,
            _ => &[],
        }
    }

    pub fn patch_tile_ids(&self) -> &[usize] {
        match &self.state {
            PatchState::Growing { patch_tile_ids, .. } => patch_tile_ids,
            _ => &[],
        }
    }

    pub fn next_tile_id(&self) -> usize {
        match &self.state {
            PatchState::Growing { next_tile_id, .. } => *next_tile_id,
            _ => 0,
        }
    }

    pub fn normalize(&mut self) {
        let (
            angles,
            edges,
            inner_chains,
            positions,
            grid,
            seg_data,
            boundary_edge_ids,
            next_edge_id,
            mut patch_tile_ids,
            _next_tile_id,
        ) = match &mut self.state {
            PatchState::Growing {
                angles,
                edges,
                inner_chains,
                positions,
                grid,
                seg_data,
                boundary_edge_ids,
                next_edge_id,
                patch_tile_ids,
                next_tile_id,
                candidates_by_start: _,
            } => (
                std::mem::take(angles),
                std::mem::take(edges),
                std::mem::take(inner_chains),
                std::mem::take(positions),
                std::mem::take(grid),
                std::mem::take(seg_data),
                std::mem::take(boundary_edge_ids),
                *next_edge_id,
                std::mem::take(patch_tile_ids),
                *next_tile_id,
            ),
            _ => return,
        };

        let n = angles.len();
        if n == 0 {
            return;
        }

        let rot = crate::intgeom::rat::lex_min_rot(&angles);

        let mut angles = angles;
        let mut edges = edges;
        let mut inner_chains = inner_chains;
        let mut positions = positions;
        let mut boundary_edge_ids = boundary_edge_ids;

        if rot != 0 {
            angles.rotate_left(rot);
            edges.rotate_left(rot);
            inner_chains.rotate_left(rot);
            patch_tile_ids.rotate_left(rot);
            boundary_edge_ids.rotate_left(rot);
            // `positions` has n + 1 entries (closing vertex repeated). Rotate
            // only the first n; the closing entry is always positions[0].
            let last_idx = positions.len() - 1;
            positions.truncate(last_idx);
            positions.rotate_left(rot);
            let first = positions[0];
            positions.push(first);
        }

        let mut remap: rustc_hash::FxHashMap<usize, usize> = rustc_hash::FxHashMap::default();
        let mut next = 0usize;
        for id in &mut patch_tile_ids {
            let new_id = *remap.entry(*id).or_insert_with(|| {
                let v = next;
                next += 1;
                v
            });
            *id = new_id;
        }

        self.state = PatchState::Growing {
            angles,
            edges,
            candidates_by_start: None,
            inner_chains,
            positions,
            grid,
            seg_data,
            boundary_edge_ids,
            next_edge_id,
            patch_tile_ids,
            next_tile_id: next,
        };
    }

    pub fn candidates_by_start(&self) -> &[Vec<PatchMatch>] {
        match &self.state {
            PatchState::Growing {
                candidates_by_start,
                ..
            } => candidates_by_start.as_deref().unwrap_or(&[]),
            _ => &[],
        }
    }

    pub fn vertex_type_at(&self, i: usize) -> Option<PatchVertexInfo> {
        let n = self.boundary_len();
        if n == 0 || i >= n {
            return None;
        }
        let edges = match &self.state {
            PatchState::Growing { edges, .. } => edges,
            _ => return None,
        };
        let angles = self.angles();
        Some(PatchVertexInfo {
            angle: angles[i],
            cw: edges[(i + n - 1) % n],
            ccw: edges[i],
        })
    }

    /// Return the junction vertex type at position `i`, or `None` if not a junction.
    pub fn junction_vertex_type_at(&self, i: usize) -> Option<OpenVertexType> {
        let n = self.boundary_len();
        if n == 0 || i >= n {
            return None;
        }
        match &self.state {
            PatchState::Growing {
                edges,
                inner_chains,
                ..
            } => {
                if !self.is_junction(i) {
                    return None;
                }
                let n = edges.len();
                Some(OpenVertexType {
                    cw: edges[(i + n - 1) % n],
                    inner: inner_chains[i].clone(),
                    ccw: edges[i],
                })
            }
            _ => None,
        }
    }

    pub fn is_junction(&self, i: usize) -> bool {
        let edges = match &self.state {
            PatchState::Growing { edges, .. } => edges,
            _ => return false,
        };
        if edges.is_empty() {
            return false;
        }
        is_junction_at(self.angles(), edges, self.match_index.tileset(), i)
    }

    pub fn neighbor_junction_offsets(&self, pos: usize) -> Option<(usize, usize)> {
        let (edges, n) = match &self.state {
            PatchState::Growing { edges, .. } => (edges, edges.len()),
            _ => return None,
        };
        if n == 0 || pos >= n {
            return None;
        }

        let mut j_cw = (pos + n - 1) % n;
        while j_cw != pos && !self.is_junction(j_cw) {
            j_cw = (j_cw + n - 1) % n;
        }

        let mut j_ccw = (pos + 1) % n;
        while j_ccw != pos && !self.is_junction(j_ccw) {
            j_ccw = (j_ccw + 1) % n;
        }

        let cw_offset = edges[j_cw].tile_offset;
        let ccw_prev = (j_ccw + n - 1) % n;
        let ccw_edge = edges[ccw_prev];
        let tile_len = self.match_index.tileset().rat(ccw_edge.tile_id).len();
        let ccw_offset = (ccw_edge.tile_offset + 1) % tile_len;

        Some((cw_offset, ccw_offset))
    }

    pub fn junction_vertices(&self) -> Vec<(usize, PatchVertexInfo)> {
        let n = self.boundary_len();
        let mut result = Vec::new();
        for i in 0..n {
            if self.is_junction(i) {
                result.push((i, self.vertex_type_at(i).unwrap()));
            }
        }
        result
    }

    /// Partition the boundary into [`TileSegment`]s — maximal contiguous
    /// runs of edges from the same tile instance.
    ///
    /// See [`TileSegment`] for the cyclic-vs-linear caveat: when position
    /// 0 is not at a junction, one cyclic tile-instance run is split into
    /// two linear segments at the array seam.
    pub fn tile_segments(&self) -> Vec<TileSegment> {
        let edges = match &self.state {
            PatchState::Growing { edges, .. } => edges,
            _ => return vec![],
        };
        let n = edges.len();
        if n == 0 {
            return vec![];
        }
        compute_segments(self.angles(), edges, self.match_index.tileset())
    }

    pub fn add_tile(&mut self, pm: &PatchMatch) -> bool {
        match &self.state {
            PatchState::Seed { tile_id, .. } => {
                let seed_id = *tile_id;
                self.init_from_first_add(seed_id, pm)
            }
            PatchState::Growing { .. } => self.add_tile_growing(pm),
        }
    }

    fn init_from_first_add(&mut self, seed_id: usize, pm: &PatchMatch) -> bool {
        let tileset = self.match_index.tileset();
        let seed_rat = tileset.rat(seed_id);
        let n = seed_rat.seq().len();
        let m = tileset.rat(pm.tile_id).seq().len();
        let mlen = pm.len;

        if mlen == 0 || mlen > n || mlen > m {
            return false;
        }

        let (_ns, seed_len, _ne) = seed_rat.get_match(
            (pm.start_a as i64, pm.start_b as i64),
            tileset.rat(pm.tile_id),
        );
        if seed_len == 0 {
            return false;
        }

        let seed_angles = seed_rat.seq().to_vec();
        let new_angles = match compute_glue_angles::<T>(&seed_angles, pm, tileset) {
            Some(a) => a,
            None => return false,
        };

        let seg_len_old = n - mlen;
        let seg_len_new = m - mlen;
        let new_len = seg_len_old + seg_len_new;
        if new_len == 0 {
            return false;
        }

        // Geometric self-intersection check.
        //
        // On the first add we don't yet have a `UnitSquareGrid` to query
        // incrementally, so we build the candidate boundary from scratch
        // and run Snake's batch validator. `add_tile_growing` uses the
        // incremental `check_segment_clear` path instead (it maintains the
        // grid across glues). Both ultimately use the same `intersect` +
        // `UnitSquareGrid` primitive, and their agreement is cross-checked
        // by `add_tile_decision_agrees_with_snake_on_spectre`.
        if Snake::<T>::try_from(new_angles.as_slice()).is_err() {
            return false;
        }

        let positions = trace_boundary_positions::<T>(&new_angles);
        let (grid, seg_data, boundary_edge_ids, next_edge_id) =
            build_boundary_grid::<T>(&positions, new_len);

        // Synthesize the seed's "old edges" so we can share build_glued_edges
        // with the growing path. The seed occupies patch_tile_id 0; the
        // glued tile becomes patch_tile_id 1.
        let seed_old_edges: Vec<EdgeInfo> = (0..n)
            .map(|i| EdgeInfo {
                tile_id: seed_id,
                tile_offset: i,
            })
            .collect();
        let seed_old_ptids = vec![0usize; n];
        let (edges, patch_tile_ids) = build_glued_edges(&seed_old_edges, &seed_old_ptids, pm, m, 1);

        debug_assert_eq!(edges.len(), new_len);

        let inner_chains = vec![vec![]; new_len];

        self.state = PatchState::Growing {
            angles: new_angles,
            edges,
            candidates_by_start: None,
            inner_chains,
            positions,
            grid,
            seg_data,
            boundary_edge_ids,
            next_edge_id,
            patch_tile_ids,
            next_tile_id: 2,
        };

        true
    }

    fn add_tile_growing(&mut self, pm: &PatchMatch) -> bool {
        // === Pre-checks against immutable state ===
        //
        // Paths 1, 2, 4 are invariant violations that legitimate callers
        // (i.e. `get_all_matches`) never produce: bounds (1) are enforced by
        // `forward_match_length`, the ±hturn rejection (2) is already done
        // by `try_glue_precomputed` inside `get_all_matches`, and full
        // closure (4) is geometrically impossible since the patch interior
        // is already filled. The `debug_assert!`s catch any bugs in tests;
        // the `return false`s are release-mode safety nets that leave
        // `self.state` untouched (no mem::take has happened yet).
        let n = self.boundary_len();
        let mlen = pm.len;
        let m = self.match_index.tileset().rat(pm.tile_id).seq().len();

        debug_assert!(
            mlen > 0 && mlen <= n && mlen <= m,
            "add_tile_growing: invalid mlen={mlen} (n={n}, m={m}); \
             get_all_matches() should never produce this"
        );
        if mlen == 0 || mlen > n || mlen > m {
            return false;
        }

        let new_angles =
            match compute_glue_angles::<T>(self.angles(), pm, self.match_index.tileset()) {
                Some(a) => a,
                None => {
                    debug_assert!(
                        false,
                        "compute_glue_angles returned None; get_all_matches() \
                         should already filter ±hturn glues via try_glue_precomputed"
                    );
                    return false;
                }
            };

        let seg_len_old = n - mlen;
        let seg_len_new = m - mlen;
        let new_len = seg_len_old + seg_len_new;
        debug_assert!(
            new_len > 0,
            "add_tile_growing: new_len == 0 is geometrically impossible — \
             would require placing the new tile entirely inside the existing patch"
        );
        if new_len == 0 {
            return false;
        }

        // === Take state. Path 3 (geometric collision) is the only remaining
        // rollback path; it mutates `grid`, so we must restore on failure. ===
        let (
            old_angles,
            edges,
            old_inner,
            positions,
            mut grid,
            seg_data,
            boundary_edge_ids,
            next_edge_id,
            patch_tile_ids,
            next_tile_id,
        ) = match &mut self.state {
            PatchState::Growing {
                angles,
                edges,
                candidates_by_start: _,
                inner_chains,
                positions,
                grid,
                seg_data,
                boundary_edge_ids,
                next_edge_id,
                patch_tile_ids,
                next_tile_id,
            } => (
                std::mem::take(angles),
                std::mem::take(edges),
                std::mem::take(inner_chains),
                std::mem::take(positions),
                std::mem::take(grid),
                std::mem::take(seg_data),
                std::mem::take(boundary_edge_ids),
                *next_edge_id,
                std::mem::take(patch_tile_ids),
                *next_tile_id,
            ),
            _ => return false,
        };

        let ccw_pos = (pm.start_a + mlen) % n;

        let removed_ids: Vec<usize> = (0..mlen)
            .map(|i| boundary_edge_ids[(pm.start_a + i) % n])
            .collect();

        for &id in &removed_ids {
            unregister_segment::<T>(&mut grid, &seg_data, id);
        }

        // Trace the new tile's boundary positions, starting from the CW
        // junction and walking the new-tile half of `new_angles`.
        let cw_junction = positions[pm.start_a];
        let ccw_junction = positions[ccw_pos];
        let initial_dir = {
            let prev = (pm.start_a + n - 1) % n;
            dir_of_edge::<T>(positions[prev], positions[pm.start_a])
        };
        let new_tile_positions =
            trace_polyline_from::<T>(cw_junction, initial_dir, &new_angles[seg_len_old..]);
        debug_assert_eq!(
            *new_tile_positions.last().unwrap(),
            ccw_junction,
            "new tile trace should end at CCW junction"
        );

        // Geometric collision check (path 3): verify the new tile's
        // segments don't cross any surviving boundary segment. The new
        // tile's final endpoint is allowed to coincide with `ccw_junction`
        // (that's where it rejoins the existing boundary).
        let all_clear =
            segments_all_clear::<T>(&grid, &seg_data, &new_tile_positions, ccw_junction);

        if !all_clear {
            for &id in &removed_ids {
                let (p1, p2) = seg_data[id];
                register_segment::<T>(&mut grid, p1, p2, id);
            }
            self.restore_growing(
                old_angles,
                edges,
                old_inner,
                positions,
                grid,
                seg_data,
                boundary_edge_ids,
                next_edge_id,
                patch_tile_ids,
                next_tile_id,
            );
            return false;
        }

        let mut new_positions = Vec::with_capacity(new_len + 1);
        for i in 0..=seg_len_old {
            new_positions.push(positions[(ccw_pos + i) % n]);
        }
        for p in new_tile_positions.iter().take(seg_len_new + 1).skip(1) {
            new_positions.push(*p);
        }

        let mut new_boundary_edge_ids = Vec::with_capacity(new_len);
        let mut new_seg_data = seg_data;
        let mut next_id = next_edge_id;

        for i in 0..seg_len_old {
            new_boundary_edge_ids.push(boundary_edge_ids[(ccw_pos + i) % n]);
        }

        for i in 0..seg_len_new {
            let id = next_id;
            next_id += 1;
            new_boundary_edge_ids.push(id);
            let p1 = new_tile_positions[i];
            let p2 = new_tile_positions[i + 1];
            new_seg_data.push((p1, p2));
            register_segment::<T>(&mut grid, p1, p2, id);
        }

        let (new_edges, new_patch_tile_ids) =
            build_glued_edges(&edges, &patch_tile_ids, pm, m, next_tile_id);
        debug_assert_eq!(new_edges.len(), new_len);

        let new_inner = update_inner_chains(&old_inner, &edges, pm, new_len, &patch_tile_ids);

        self.state = PatchState::Growing {
            angles: new_angles,
            edges: new_edges,
            candidates_by_start: None,
            inner_chains: new_inner,
            positions: new_positions,
            grid,
            seg_data: new_seg_data,
            boundary_edge_ids: new_boundary_edge_ids,
            next_edge_id: next_id,
            patch_tile_ids: new_patch_tile_ids,
            next_tile_id: next_tile_id + 1,
        };

        true
    }

    #[allow(clippy::too_many_arguments)]
    fn restore_growing(
        &mut self,
        angles: Vec<i8>,
        edges: Vec<EdgeInfo>,
        inner_chains: Vec<Vec<EdgeInfo>>,
        positions: Vec<T>,
        grid: UnitSquareGrid,
        seg_data: Vec<(T, T)>,
        boundary_edge_ids: Vec<usize>,
        next_edge_id: usize,
        patch_tile_ids: Vec<usize>,
        next_tile_id: usize,
    ) {
        self.state = PatchState::Growing {
            angles,
            edges,
            candidates_by_start: None,
            inner_chains,
            positions,
            grid,
            seg_data,
            boundary_edge_ids,
            next_edge_id,
            patch_tile_ids,
            next_tile_id,
        };
    }

    fn compute_seed_matches(
        match_index: &MatchTypeIndex<T>,
        tileset: &TileSet<T>,
        seed_tile_id: usize,
    ) -> Vec<PatchMatch> {
        let seed = tileset.rat(seed_tile_id);
        let seed_seq = seed.seq();
        let n = seed_seq.len();
        let mut seen: FxHashSet<(usize, usize, usize, usize)> = FxHashSet::default();
        let mut matches = Vec::new();

        for offset in 0..n {
            for cand in match_index.candidates_starting_at(seed_tile_id, offset) {
                let tile_b = tileset.rat(cand.tile_b);
                let (ns, len, ne) = seed.get_match((offset as i64, cand.start_b as i64), tile_b);
                if len == 0 {
                    continue;
                }
                let ns_u = ns.rem_euclid(n as i64) as usize;
                let ne_u = ne.rem_euclid(tile_b.len() as i64) as usize;
                if !junction_gap_nonnegative(seed_seq, ns_u, len, tile_b.seq(), ne_u) {
                    continue;
                }
                if let Ok(_glued) = seed.try_glue_precomputed((ns, len, ne), tile_b, true) {
                    let key = (ns_u, len, ne_u, cand.tile_b);
                    if seen.insert(key) {
                        matches.push(PatchMatch {
                            start_a: ns_u,
                            len,
                            start_b: ne_u,
                            tile_id: cand.tile_b,
                        });
                    }
                }
            }
        }

        matches
    }
}

/// Build the post-glue `(edges, patch_tile_ids)` vectors.
///
/// Surviving old boundary edges come first (`seg_len_old` entries starting
/// at `(pm.start_a + pm.len) % n`), then the new tile's surviving edges
/// (`seg_len_new` entries starting at `pm.start_b`). The new tile's entries
/// are tagged with `new_tile_id` in the returned ptid vector.
fn build_glued_edges(
    old_edges: &[EdgeInfo],
    old_patch_tile_ids: &[usize],
    pm: &PatchMatch,
    m_tile: usize,
    new_tile_id: usize,
) -> (Vec<EdgeInfo>, Vec<usize>) {
    let n = old_edges.len();
    let mlen = pm.len;
    let seg_len_old = n - mlen;
    let seg_len_new = m_tile - mlen;
    let ccw_pos = (pm.start_a + mlen) % n;
    let new_len = seg_len_old + seg_len_new;

    let mut new_edges = Vec::with_capacity(new_len);
    let mut new_ptids = Vec::with_capacity(new_len);

    for i in 0..seg_len_old {
        new_edges.push(old_edges[(ccw_pos + i) % n]);
        new_ptids.push(old_patch_tile_ids[(ccw_pos + i) % n]);
    }
    new_edges.push(EdgeInfo {
        tile_id: pm.tile_id,
        tile_offset: pm.start_b,
    });
    new_ptids.push(new_tile_id);
    for k in 1..seg_len_new {
        new_edges.push(EdgeInfo {
            tile_id: pm.tile_id,
            tile_offset: (pm.start_b + k) % m_tile,
        });
        new_ptids.push(new_tile_id);
    }
    (new_edges, new_ptids)
}

/// Partition the boundary into [`TileSegment`]s — maximal contiguous runs
/// of edges from the same tile instance (no junction between them).
///
/// The segment-break condition is the **canonical** junction check
/// [`is_junction_at`] — i.e. `angles[i] != tile.seq()[edges[i].tile_offset]`.
/// Every glue updates the boundary angle at the new junction
/// (`compute_glue_angles` produces an angle that for non-degenerate tiles
/// is provably different from either tile's natural internal angle at the
/// junction position), so this check reliably detects every tile-instance
/// boundary.
///
/// Note: a simpler offset-arithmetic heuristic (same `tile_id` + contiguous
/// `tile_offset` from the segment start) is **not** used here. It would
/// over-segment on intra-tile wrap-arounds (a single tile instance whose
/// surviving boundary edges have offsets like `2,3,4,5,0`) and could
/// under-segment in pathological same-tile-id glues whose offsets happen
/// to line up. The angle-based check is canonical and avoids both.
///
/// Cyclic-vs-linear caveat: when position 0 is not at a junction, a single
/// cyclic tile-instance run is split into two linear segments at the array
/// seam (`[0, first_junction)` and `[last_junction, n)`). See
/// [`TileSegment`] for details.
fn compute_segments<T: IsComplex + IsRingOrField + Units>(
    angles: &[i8],
    edges: &[EdgeInfo],
    tileset: &TileSet<T>,
) -> Vec<TileSegment> {
    let n = edges.len();
    if n == 0 {
        return vec![];
    }
    let mut segments: Vec<TileSegment> = Vec::new();
    let mut seg_start = 0;
    for i in 1..=n {
        let break_here = i == n || is_junction_at(angles, edges, tileset, i);
        if break_here {
            segments.push(TileSegment {
                patch_start: seg_start,
                patch_end: i,
                tile_id: edges[seg_start].tile_id,
                offset_start: edges[seg_start].tile_offset,
            });
            seg_start = i;
        }
    }
    segments
}

fn compute_junctions<T: IsComplex + IsRingOrField + Units>(
    angles: &[i8],
    edges: &[EdgeInfo],
    tileset: &TileSet<T>,
) -> Vec<usize> {
    (0..edges.len())
        .filter(|&i| is_junction_at(angles, edges, tileset, i))
        .collect()
}

#[allow(clippy::too_many_arguments)]
fn compute_candidates_at_position<T: IsComplex + IsRingOrField + Units>(
    pos: usize,
    angles: &[i8],
    rat: &Rat<T>,
    n: usize,
    tile_id: usize,
    tile_offset: usize,
    match_index: &MatchTypeIndex<T>,
    tileset: &TileSet<T>,
    global_seen: &mut FxHashSet<(usize, usize, usize, usize)>,
    result: &mut [Vec<PatchMatch>],
) {
    for cand in match_index.candidates_starting_at(tile_id, tile_offset) {
        append_match_candidate(pos, angles, rat, n, cand, tileset, global_seen, result);
    }
}

#[allow(clippy::too_many_arguments)]
fn append_match_candidate<T: IsComplex + IsRingOrField + Units>(
    pos: usize,
    angles: &[i8],
    rat: &Rat<T>,
    n: usize,
    cand: &CandidateMatch,
    tileset: &TileSet<T>,
    seen: &mut FxHashSet<(usize, usize, usize, usize)>,
    result: &mut [Vec<PatchMatch>],
) {
    let tile_b = tileset.rat(cand.tile_b);
    let (ns, len, ne) = rat.get_match((pos as i64, cand.start_b as i64), tile_b);
    if len == 0 {
        return;
    }
    let ns_u = ns.rem_euclid(n as i64) as usize;
    let ne_u = ne.rem_euclid(tile_b.len() as i64) as usize;
    if !junction_gap_nonnegative(angles, ns_u, len, tile_b.seq(), ne_u) {
        return;
    }
    let key = (ns_u, len, ne_u, cand.tile_b);
    if !seen.insert(key) {
        return;
    }
    if rat
        .try_glue_precomputed((ns, len, ne), tile_b, true)
        .is_ok()
    {
        result[ns_u].push(PatchMatch {
            start_a: ns_u,
            len,
            start_b: ne_u,
            tile_id: cand.tile_b,
        });
    }
}

/// Enumerate single-edge match candidates that anchor at junction position `pos`.
///
/// For each `(tile_id_b, ib)` pair, checks single-edge compatibility, dedups
/// via `seen`, applies the `keep(start_a, len)` filter before the relatively
/// expensive `try_glue_precomputed`, and hands surviving matches to `emit`.
#[allow(clippy::too_many_arguments)]
fn enumerate_junction_candidates_at<T: IsComplex + IsRingOrField + Units>(
    pos: usize,
    angles: &[i8],
    rat: &Rat<T>,
    n: usize,
    tileset: &TileSet<T>,
    seen: &mut FxHashSet<(usize, usize, usize, usize)>,
    keep: impl Fn(usize, usize) -> bool,
    mut emit: impl FnMut(PatchMatch),
) {
    for tile_id_b in 0..tileset.num_tiles() {
        let tile_b = tileset.rat(tile_id_b);
        let b_seq = tile_b.seq();
        let m = b_seq.len();
        for ib in 0..m {
            if !is_single_edge_candidate(angles, pos, b_seq, ib) {
                continue;
            }
            let (ns, len, ne) = rat.get_match((pos as i64, ib as i64), tile_b);
            if len != 1 {
                continue;
            }
            let ns_u = ns.rem_euclid(n as i64) as usize;
            let ne_u = ne.rem_euclid(m as i64) as usize;
            let key = (ns_u, len, ne_u, tile_id_b);
            if !seen.insert(key) {
                continue;
            }
            if !keep(ns_u, len) {
                continue;
            }
            if rat
                .try_glue_precomputed((ns, len, ne), tile_b, true)
                .is_ok()
            {
                emit(PatchMatch {
                    start_a: ns_u,
                    len,
                    start_b: ne_u,
                    tile_id: tile_id_b,
                });
            }
        }
    }
}

/// Returns true if vertex `index` is touched by a match of `len` edges
/// starting at edge position `start` on a cyclic boundary of length `n`.
///
/// A match of length `len` starting at edge `start` covers edges
/// `[start, start+1, ..., start+len-1]` and therefore *touches* the
/// `len + 1` boundary **vertices** `[start, start+1, ..., start+len]`
/// (the CW endpoint of the first edge through the CCW endpoint of the
/// last edge). This function is inclusive on both ends — both the CW
/// vertex (`start`) and the CCW vertex (`start + len`) are considered
/// to be touched.
///
/// Used by `compute_candidates_covering_position` and
/// `get_matches_touching_vertex` to find matches incident with a given
/// boundary vertex.
pub fn cyclic_range_contains(start: usize, len: usize, index: usize, n: usize) -> bool {
    if len == 0 || n == 0 {
        return false;
    }
    let end = start + len;
    if end <= n {
        index >= start && index <= end
    } else {
        index >= start || index <= end % n
    }
}

/// Compute the new boundary angles after gluing `pm.tile_id`'s tile onto
/// the boundary `angles` along the match described by `pm`.
///
/// Wraps [`angles::glue_raw_angles`] and additionally rejects glues that
/// would produce a ±half-turn at either of the two new junction angles
/// (a half-turn boundary angle means the boundary doubles back on itself
/// — a degenerate pinched vertex, which is not a valid patch boundary).
///
/// Returns `None` if either the raw glue would be empty or a ±hturn
/// junction would result. This is the live-patch glue path; the parallel
/// raw-boundary path used by witness construction (`glue_match_to_raw_boundary`)
/// is intentionally permissive at this level — see
/// [`construct_witness_from_vt_sequence`] for the open-VT enforcement.
pub fn compute_glue_angles<T: IsComplex + IsRingOrField + Units>(
    angles: &[i8],
    pm: &PatchMatch,
    tileset: &Arc<TileSet<T>>,
) -> Option<Vec<i8>> {
    let other_seq = tileset.rat(pm.tile_id).seq();
    let gr = angles::glue_raw_angles::<T>(angles, other_seq, pm.start_a, pm.len, pm.start_b)?;
    if let (Some(a_yx), Some(a_xy)) = (gr.a_yx, gr.a_xy) {
        if a_yx.abs() == T::hturn() || a_xy.abs() == T::hturn() {
            return None;
        }
    }
    Some(gr.angles)
}

/// Trace a polyline of unit edges starting at `start` facing
/// `initial_dir`, applying each angle as a turn.
///
/// Returns `angles.len() + 1` positions (the starting vertex plus one
/// after each edge).
fn trace_polyline_from<T: IsComplex + IsRingOrField + Units>(
    start: T,
    initial_dir: i8,
    angles: &[i8],
) -> Vec<T> {
    let mut positions = Vec::with_capacity(angles.len() + 1);
    positions.push(start);
    let mut dir = initial_dir;
    for &a in angles {
        dir = (dir as i64 + a as i64).rem_euclid(T::turn() as i64) as i8;
        let last = *positions.last().unwrap();
        positions.push(last + T::unit(dir));
    }
    positions
}

fn trace_boundary_positions<T: IsComplex + IsRingOrField + Units>(angles: &[i8]) -> Vec<T> {
    trace_polyline_from(T::zero(), 0, angles)
}

fn build_boundary_grid<T: IsComplex + IsRingOrField + Units>(
    positions: &[T],
    n: usize,
) -> (UnitSquareGrid, Vec<(T, T)>, Vec<usize>, usize) {
    let mut grid = UnitSquareGrid::new();
    let mut seg_data = Vec::with_capacity(n);
    let boundary_edge_ids: Vec<usize> = (0..n).collect();
    for i in 0..n {
        let p1 = positions[i];
        let p2 = positions[i + 1];
        seg_data.push((p1, p2));
        let cells = UnitSquareGrid::seg_neighborhood_of(p1, p2);
        for &cell in &cells {
            grid.add(cell, i);
        }
    }
    (grid, seg_data, boundary_edge_ids, n)
}

fn register_segment<T: IsComplex + IsRingOrField + Units>(
    grid: &mut UnitSquareGrid,
    p1: T,
    p2: T,
    id: usize,
) {
    let cells = UnitSquareGrid::seg_neighborhood_of(p1, p2);
    for &cell in &cells {
        grid.add(cell, id);
    }
}

fn unregister_segment<T: IsComplex + IsRingOrField + Units>(
    grid: &mut UnitSquareGrid,
    seg_data: &[(T, T)],
    id: usize,
) {
    let (p1, p2) = seg_data[id];
    let cells = UnitSquareGrid::seg_neighborhood_of(p1, p2);
    for &cell in &cells {
        grid.remove(cell, id);
    }
}

fn check_segment_clear<T: IsComplex + IsRingOrField + Units>(
    grid: &UnitSquareGrid,
    seg_data: &[(T, T)],
    p1: T,
    p2: T,
    allowed_endpoint: Option<T>,
) -> bool {
    let cells = UnitSquareGrid::seg_neighborhood_of(p1, p2);
    let near_ids = grid.get_cells(&cells);
    for &id in &near_ids {
        let (x, y) = seg_data[id];
        let is_allowed = allowed_endpoint == Some(p2);
        if !is_allowed && (p2 == x || p2 == y) {
            return false;
        }
        if intersect(&(p1, p2), &(x, y)) {
            return false;
        }
    }
    true
}

/// Walk consecutive segments of a polyline, checking each against the
/// boundary grid via [`check_segment_clear`]. The last segment's CCW
/// endpoint is permitted to coincide with `allowed_last_endpoint` (this
/// is where the new tile rejoins the existing boundary).
fn segments_all_clear<T: IsComplex + IsRingOrField + Units>(
    grid: &UnitSquareGrid,
    seg_data: &[(T, T)],
    positions: &[T],
    allowed_last_endpoint: T,
) -> bool {
    let n_edges = positions.len().saturating_sub(1);
    for i in 0..n_edges {
        let allowed = (i == n_edges - 1).then_some(allowed_last_endpoint);
        if !check_segment_clear::<T>(grid, seg_data, positions[i], positions[i + 1], allowed) {
            return false;
        }
    }
    true
}

/// Find the direction `dir` such that `from + T::unit(dir) == to`.
///
/// Caller invariant: `(from, to)` must be a unit-length boundary edge, i.e.
/// `to - from` is a unit vector in the cyclotomic ring. This is enforced
/// upstream by [`trace_boundary_positions`] and the live patch growth
/// path (every boundary edge is constructed as a single unit step). A
/// failure here means upstream broke the unit-edge invariant — an
/// internal bug, not bad input.
fn dir_of_edge<T: IsComplex + IsRingOrField + Units>(from: T, to: T) -> i8 {
    let d = to - from;
    for dir in 0..T::turn() {
        if T::unit(dir) == d {
            return dir;
        }
    }
    unreachable!("dir_of_edge: ({from:?} -> {to:?}) is not a unit vector");
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic::{ZZ12, ZZ4};
    use crate::intgeom::matchtypes::MatchTypeIndex;
    use crate::intgeom::snake::Snake;
    use crate::intgeom::tiles;
    use std::collections::BTreeMap;

    fn square_patch() -> GrowingPatch<ZZ4> {
        let sq: Snake<ZZ4> = tiles::square();
        let rat = Rat::try_from(&sq).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        GrowingPatch::new(ts, 0)
    }

    fn hex_patch() -> GrowingPatch<ZZ12> {
        let hex: Snake<ZZ12> = tiles::hexagon();
        let rat = Rat::try_from(&hex).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        GrowingPatch::new(ts, 0)
    }

    #[test]
    fn seed_patch_has_no_boundary() {
        let gp = square_patch();
        assert!(!gp.is_growing());
        assert_eq!(gp.boundary_len(), 0);
        assert_eq!(gp.edges().len(), 0);
    }

    /// Cloning a `Seed` must preserve the seed-match set (the seed
    /// boundary never changes, so the matches computed at construction
    /// are still authoritative for the clone). Regression test for an
    /// earlier latent bug in the (now-removed) `clone_for_mutation`
    /// method, which used to reset `cached_matches` to empty.
    #[test]
    fn clone_preserves_seed_matches() {
        let gp = hex_patch();
        let original = gp.get_all_matches();
        assert!(
            !original.is_empty(),
            "fixture: seed should have non-empty matches"
        );
        let clone = gp.clone();
        let cloned_matches = clone.get_all_matches();
        assert_eq!(
            cloned_matches, original,
            "cloned Seed must report the same matches as the original"
        );
    }

    #[test]
    fn first_add_produces_growing() {
        let mut gp = hex_patch();
        let pm = gp.get_all_matches()[0].clone();
        assert!(gp.add_tile(&pm), "first add");

        assert!(gp.is_growing());
        assert_eq!(gp.boundary_len(), 12 - 2 * pm.len);
        assert_eq!(gp.edges().len(), gp.boundary_len());
        assert_eq!(gp.angles().len(), gp.boundary_len());
    }

    #[test]
    fn junction_vertex_ids_nonempty_after_each_add() {
        let mut gp = hex_patch();
        let mut step = 0;
        // Recompute candidates after each add — pms from a stale patch
        // state are not valid input to `add_tile` once the boundary changes.
        while step < 3 {
            let candidates = gp.get_all_matches();
            let pm = match candidates.first() {
                Some(pm) => pm.clone(),
                None => break,
            };
            if !gp.add_tile(&pm) {
                break;
            }
            assert!(
                !gp.edges().is_empty(),
                "step {step}: edges should not be empty"
            );
            assert!(
                !gp.junction_vertices().is_empty(),
                "step {step}: should have junction vertices"
            );
            step += 1;
        }
        assert!(step > 0, "expected at least one successful add");
    }

    #[test]
    fn hexagon_all_36_matches_produce_valid_bi_hexes() {
        let gp = hex_patch();
        let matches = gp.get_all_matches();
        assert_eq!(matches.len(), 36, "hex self-matches = 36");

        for pm in &matches {
            let mut gp2 = hex_patch();
            assert!(gp2.add_tile(pm), "first add should succeed for pm {:?}", pm);
            assert!(gp2.is_growing());
            assert_eq!(gp2.boundary_len(), 12 - 2 * pm.len);
            assert_eq!(gp2.edges().len(), gp2.boundary_len());

            let rat = gp2.to_rat();
            assert!(
                Snake::<ZZ12>::try_from(rat.seq()).is_ok(),
                "valid snake for pm {:?}",
                pm
            );
        }
    }

    #[test]
    fn square_all_16_matches_produce_valid_bi_squares() {
        let gp = square_patch();
        let matches = gp.get_all_matches();
        assert_eq!(matches.len(), 16, "square self-matches = 16");

        for pm in &matches {
            let mut gp2 = square_patch();
            assert!(gp2.add_tile(pm), "first add should succeed for pm {:?}", pm);
            assert!(gp2.is_growing());
            assert_eq!(gp2.boundary_len(), 8 - 2 * pm.len);

            let rat = gp2.to_rat();
            assert!(
                Snake::<ZZ4>::try_from(rat.seq()).is_ok(),
                "valid snake for pm {:?}",
                pm
            );
        }
    }

    #[test]
    fn to_rat_matches_direct_glue_for_all_matches() {
        let gp = hex_patch();
        let matches = gp.get_all_matches();
        let ts = gp.tileset().clone();

        for pm in &matches {
            let mut gp2 = GrowingPatch::<ZZ12>::new(Arc::clone(&ts), 0);
            assert!(gp2.add_tile(pm), "first add");
            let rat = gp2.to_rat();

            let seed_rat = ts.rat(0);
            let new_rat = ts.rat(pm.tile_id);
            let glued = seed_rat.try_glue((pm.start_a as i64, pm.start_b as i64), new_rat);
            match glued {
                Ok(g) => assert_eq!(rat.seq(), g.seq(), "mismatch for pm {:?}", pm),
                Err(e) => panic!("glue failed for pm {:?}: {}", pm, e),
            }
        }
    }

    #[test]
    fn edges_self_consistent() {
        let gp_sq: GrowingPatch<ZZ4> = square_patch();
        for pm in gp_sq.get_all_matches() {
            let mut gp2 = gp_sq.clone();
            if !gp2.add_tile(&pm) || !gp2.is_growing() {
                continue;
            }
            verify_edges_consistency(&gp2, gp2.tileset(), &format!("bi-sq pm {:?}", pm));
        }
        let gp_hex: GrowingPatch<ZZ12> = hex_patch();
        for pm in gp_hex.get_all_matches() {
            let mut gp2 = gp_hex.clone();
            if !gp2.add_tile(&pm) || !gp2.is_growing() {
                continue;
            }
            verify_edges_consistency(&gp2, gp2.tileset(), &format!("bi-hex pm {:?}", pm));
            for pm2 in gp2.get_all_matches() {
                let mut gp3 = gp2.clone();
                if gp3.add_tile(&pm2) && gp3.is_growing() {
                    verify_edges_consistency(&gp3, gp3.tileset(), "3-hex");
                }
            }
        }
    }

    /// Assert two angle sequences are equal up to cyclic rotation (i.e.
    /// describe the same closed shape). Stronger than the
    /// sort-the-multisets comparison: two unrelated sequences with the
    /// same angle counts would pass a sorted equality but fail this.
    fn assert_same_cyclic_shape(a: &[i8], b: &[i8], label: &str) {
        assert_eq!(
            a.len(),
            b.len(),
            "{label}: angle sequences have different lengths ({} vs {})",
            a.len(),
            b.len(),
        );
        if a.is_empty() {
            return;
        }
        let mut a_canon = a.to_vec();
        let a_rot = crate::intgeom::rat::lex_min_rot(&a_canon);
        a_canon.rotate_left(a_rot);
        let mut b_canon = b.to_vec();
        let b_rot = crate::intgeom::rat::lex_min_rot(&b_canon);
        b_canon.rotate_left(b_rot);
        assert_eq!(
            a_canon, b_canon,
            "{label}: angle sequences are not cyclic rotations of each other"
        );
    }

    fn verify_edges_consistency<T: IsComplex + IsRingOrField + Units>(
        gp: &GrowingPatch<T>,
        ts: &Arc<TileSet<T>>,
        label: &str,
    ) {
        let n = gp.boundary_len();
        assert!(n > 0, "[{}] patch should be growing", label);
        let edges = gp.edges();
        assert_eq!(edges.len(), n, "[{}] edges length", label);

        for (i, edge) in edges.iter().enumerate().take(n) {
            assert!(
                edge.tile_id < ts.num_tiles(),
                "[{}] pos {}: invalid tile_id {}",
                label,
                i,
                edge.tile_id
            );
            let tile_len = ts.rat(edge.tile_id).len();
            assert!(
                edge.tile_offset < tile_len,
                "[{}] pos {}: invalid offset {} for tile {} (len {})",
                label,
                i,
                edge.tile_offset,
                edge.tile_id,
                tile_len
            );
        }

        for i in 0..n {
            let j = (i + 1) % n;
            if edges[i].tile_id == edges[j].tile_id && !gp.is_junction(i) && !gp.is_junction(j) {
                let tile_len = ts.rat(edges[i].tile_id).len();
                let expected_next = (edges[i].tile_offset + 1) % tile_len;
                assert_eq!(
                    edges[j].tile_offset, expected_next,
                    "[{}] pos {}->{}: same-tile continuation expected offset {} got {}",
                    label, i, j, expected_next, edges[j].tile_offset
                );
            }
        }

        let angles = gp.angles();
        assert_eq!(angles.len(), n, "[{}] angles length", label);
    }

    /// For each junction position on `glued`, build the minimal VT witness
    /// and assert that extracting the VT from the witness yields the original.
    /// Returns the number of junctions checked (asserts at least one).
    fn assert_minimal_witness_roundtrips_for<T: IsComplex + IsRingOrField + Units>(
        glued: &GrowingPatch<T>,
        mi: &Arc<MatchTypeIndex<T>>,
        label: &str,
    ) {
        let mut checked = 0;
        for pos in 0..glued.boundary_len() {
            let vt = match glued.junction_vertex_type_at(pos) {
                Some(vt) => vt,
                None => continue,
            };
            let (witness, wpos) = GrowingPatch::construct_minimal_witness(&vt, Arc::clone(mi))
                .unwrap_or_else(|| {
                    panic!("{label}: construct_minimal_witness failed at pos={pos} vt={vt:?}")
                });
            let reconstructed = vertex_type_raw_from(witness.edges(), witness.inner_chains(), wpos);
            assert_eq!(
                reconstructed, vt,
                "{label}: roundtrip failed at pos={pos} for vt={vt:?}",
            );
            checked += 1;
        }
        assert!(checked > 0, "{label}: expected at least one junction");
    }

    /// For each junction position on `brute`, construct the minimal witness,
    /// locate the matching position in the witness boundary, and assert the
    /// witness/brute VTs agree. Also asserts that the witness angle multiset
    /// equals the brute-force angle multiset.
    fn assert_witness_matches_brute_force<T: IsComplex + IsRingOrField + Units>(
        brute: &GrowingPatch<T>,
        mi: &Arc<MatchTypeIndex<T>>,
        label: &str,
    ) {
        let brute_angles = brute.angles().to_vec();
        let brute_edges = brute.edges().to_vec();
        let brute_inner = brute.inner_chains().to_vec();

        for pos in 0..brute.boundary_len() {
            let vt = match brute.junction_vertex_type_at(pos) {
                Some(vt) => vt,
                None => continue,
            };
            let (witness, _wpos) = GrowingPatch::construct_minimal_witness(&vt, Arc::clone(mi))
                .unwrap_or_else(|| panic!("{label}: witness construction failed at pos={pos}"));
            let w_edges = witness.edges();
            let w_inner = witness.inner_chains();

            let mut found = false;
            for wpos in 0..witness.boundary_len() {
                let wvt = vertex_type_raw_from(w_edges, w_inner, wpos);
                if wvt == vt {
                    let brute_vt = vertex_type_raw_from(&brute_edges, &brute_inner, pos);
                    assert_eq!(
                        wvt, brute_vt,
                        "{label}: witness VT != brute-force VT at pos={pos}"
                    );
                    found = true;
                    break;
                }
            }
            assert!(
                found,
                "{label}: no matching position in witness for vt={vt:?} at pos={pos}"
            );

            // For these fixtures (bi-hex and bi-square), the brute patch
            // *is* the minimum-witness shape, so witness and brute should
            // describe the same closed shape — same boundary up to
            // cyclic rotation.
            assert_same_cyclic_shape(
                witness.angles(),
                &brute_angles,
                &format!("{label}: witness vs brute"),
            );
        }
    }

    /// For each junction on `glued`, assert that `junction_angle_sequence`
    /// (a) ends at the witness junction angle, and (b) is monotonically
    /// non-increasing. Returns the number of junctions checked (asserts at
    /// least one).
    fn assert_junction_angle_sequence_valid<T: IsComplex + IsRingOrField + Units>(
        glued: &GrowingPatch<T>,
        mi: &Arc<MatchTypeIndex<T>>,
        label: &str,
    ) {
        let tileset = mi.tileset();
        let mut checked = 0;
        for pos in 0..glued.boundary_len() {
            let vt = match glued.junction_vertex_type_at(pos) {
                Some(vt) => vt,
                None => continue,
            };
            let angles = junction_angle_sequence::<T>(&vt, tileset.as_ref());
            let (witness, wpos) =
                GrowingPatch::construct_minimal_witness(&vt, Arc::clone(mi)).expect("witness");
            assert_eq!(
                *angles.last().unwrap(),
                witness.angles()[wpos],
                "{label}: last angle should match witness junction angle for vt={vt:?}",
            );
            assert!(
                angles[0] > 0,
                "{label}: seed angle should be positive for vt={vt:?} (convex-tile invariant)",
            );
            for i in 1..angles.len() {
                assert!(
                    angles[i] <= angles[i - 1],
                    "{label}: angles should be monotone decreasing at i={i} for vt={vt:?}: {angles:?}",
                );
            }
            checked += 1;
        }
        assert!(checked > 0, "{label}: expected at least one junction");
    }

    #[test]
    fn edges_mixed_consistency() {
        let hex_snake: Snake<ZZ12> = tiles::hexagon();
        let sq_snake: Snake<ZZ12> = tiles::square();
        let hex_rat = Rat::try_from(&hex_snake).unwrap();
        let sq_rat = Rat::try_from(&sq_snake).unwrap();
        let ts = Arc::new(TileSet::new(vec![hex_rat, sq_rat]));

        for seed_id in 0..ts.num_tiles() {
            let gp = GrowingPatch::<ZZ12>::new(Arc::clone(&ts), seed_id);
            for pm in gp.get_all_matches() {
                let mut gp2 = GrowingPatch::new(Arc::clone(&ts), seed_id);
                if gp2.add_tile(&pm) && gp2.is_growing() {
                    verify_edges_consistency(
                        &gp2,
                        &ts,
                        &format!("mixed seed={} pm {:?}", seed_id, pm),
                    );
                }
            }
        }
    }

    #[test]
    fn brute_force_squares_up_to_4_tiles() {
        let sq: Snake<ZZ4> = tiles::square();
        let rat = Rat::try_from(&sq).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let patches = brute_force_patches(&ts, 4);

        let mut by_tiles: BTreeMap<usize, (usize, usize)> = BTreeMap::new();
        for ways in patches.values() {
            let n = ways[0].len() + 1;
            let e = by_tiles.entry(n).or_insert((0, 0));
            e.0 += 1;
            e.1 += ways.len();
        }

        assert_eq!(
            by_tiles.get(&1).map(|(s, _)| *s).unwrap_or(0),
            1,
            "1 mono-square"
        );
        assert_eq!(by_tiles.get(&2), Some(&(1, 16)), "1 bi-square, 16 ways");
        assert!(
            by_tiles.get(&3).map(|(s, _)| *s).unwrap_or(0) >= 2,
            "at least 2 tri-squares"
        );
    }

    #[test]
    fn brute_force_hexagons_up_to_3_tiles() {
        let hex: Snake<ZZ12> = tiles::hexagon();
        let rat = Rat::try_from(&hex).unwrap();
        let ts = Arc::new(TileSet::new(vec![rat]));
        let patches = brute_force_patches(&ts, 3);

        let mut by_tiles: BTreeMap<usize, usize> = BTreeMap::new();
        for ways in patches.values() {
            let n = ways[0].len() + 1;
            by_tiles.entry(n).and_modify(|c| *c += 1).or_insert(1);
        }

        assert_eq!(by_tiles.get(&1).copied().unwrap_or(0), 1, "1 mono-hex");
        assert_eq!(by_tiles.get(&2).copied().unwrap_or(0), 1, "1 bi-hex");
        assert!(
            by_tiles.get(&3).copied().unwrap_or(0) >= 1,
            "at least 1 tri-hex"
        );
    }

    fn brute_force_recurse<T: IsComplex + IsRingOrField + Units>(
        gp: &mut GrowingPatch<T>,
        history: &mut Vec<PatchMatch>,
        max_tiles: usize,
        results: &mut BTreeMap<Rat<T>, Vec<Vec<PatchMatch>>>,
    ) {
        let num_tiles = history.len() + 1;
        let rat = gp.to_rat();
        results.entry(rat).or_default().push(history.clone());

        if num_tiles >= max_tiles || !gp.is_growing() {
            return;
        }

        for pm in &gp.get_all_matches() {
            let mut gp2 = gp.clone();
            if gp2.add_tile(pm) {
                history.push(pm.clone());
                brute_force_recurse(&mut gp2, history, max_tiles, results);
                history.pop();
            }
        }
    }

    fn brute_force_patches<T: IsComplex + IsRingOrField + Units>(
        ts: &Arc<TileSet<T>>,
        max_tiles: usize,
    ) -> BTreeMap<Rat<T>, Vec<Vec<PatchMatch>>> {
        let mut results: BTreeMap<Rat<T>, Vec<Vec<PatchMatch>>> = BTreeMap::new();
        results
            .entry(ts.rat(0).clone())
            .or_default()
            .push(Vec::new());

        let seed_matches = GrowingPatch::new(Arc::clone(ts), 0).get_all_matches();
        for pm in &seed_matches {
            let mut gp = GrowingPatch::new(Arc::clone(ts), 0);
            assert!(gp.add_tile(pm), "first add");
            let mut history = vec![pm.clone()];
            brute_force_recurse(&mut gp, &mut history, max_tiles, &mut results);
        }

        results
    }

    #[test]
    fn inner_chains_empty_after_first_glue() {
        let gp = hex_patch();
        let pm = gp.get_all_matches()[0].clone();
        let mut gp2 = gp.clone();
        assert!(gp2.add_tile(&pm), "first add");
        let inner = match &gp2.state {
            PatchState::Growing { inner_chains, .. } => inner_chains.clone(),
            _ => panic!("expected Growing"),
        };
        for (i, chain) in inner.iter().enumerate() {
            assert!(
                chain.is_empty(),
                "inner chain at position {i} should be empty after first glue, got {chain:?}"
            );
        }
    }

    #[test]
    fn inner_chains_grow_on_second_glue() {
        let gp = hex_patch();
        let first_match = gp.get_all_matches()[0].clone();
        let mut gp2 = gp.clone();
        assert!(gp2.add_tile(&first_match), "first add");

        let candidates = gp2.get_all_matches();
        let second = candidates
            .iter()
            .find(|pm| pm.len == 1)
            .expect("need len-1 match");
        let mut gp3 = gp2.clone();
        assert!(gp3.add_tile(second), "second add");

        let n = gp3.boundary_len();
        let edges = gp3.edges();
        let ptids = gp3.patch_tile_ids();
        let inners = gp3.inner_chains();

        for pos in 0..n {
            let prev = (pos + n - 1) % n;
            let cw_ptid = ptids[prev];
            let ccw_ptid = ptids[pos];
            for entry in &inners[pos] {
                assert_ne!(
                    entry.tile_id, edges[prev].tile_id,
                    "inner at {pos} should not be from CW tile"
                );
                assert_ne!(
                    entry.tile_id, edges[pos].tile_id,
                    "inner at {pos} should not be from CCW tile"
                );
                assert!(
                    cw_ptid != ccw_ptid || inners[pos].is_empty(),
                    "when CW and CCW have same ptid, inner should be empty at {pos}"
                );
            }
        }
    }

    #[test]
    fn junction_vertex_type_roundtrip_after_first_glue() {
        let gp = hex_patch();
        let pm = gp.get_all_matches()[0].clone();
        let mut gp2 = gp.clone();
        assert!(gp2.add_tile(&pm), "first add");
        let n = gp2.boundary_len();
        let mut junction_count = 0;
        for i in 0..n {
            if let Some(vt) = gp2.junction_vertex_type_at(i) {
                assert!(vt.inner.is_empty(), "inner should be empty at pos {i}");
                junction_count += 1;
            }
        }
        assert!(junction_count > 0, "should have at least one junction");
    }

    #[test]
    fn construct_minimal_witness_hex_roundtrip() {
        let gp = hex_patch();
        let mi = gp.match_index().clone();
        for pm in &gp.get_all_matches() {
            let mut glued = gp.clone();
            assert!(glued.add_tile(pm), "glue should succeed");
            assert_minimal_witness_roundtrips_for(&glued, &mi, &format!("hex pm {:?}", pm));
        }
    }

    #[test]
    fn construct_minimal_witness_square_roundtrip() {
        let gp = square_patch();
        let mi = gp.match_index().clone();
        for pm in &gp.get_all_matches() {
            let mut glued = gp.clone();
            assert!(glued.add_tile(pm), "glue should succeed");
            assert_minimal_witness_roundtrips_for(&glued, &mi, &format!("square pm {:?}", pm));
        }
    }

    #[test]
    fn construct_minimal_witness_hex_with_inner() {
        let gp = hex_patch();
        let mi = gp.match_index().clone();
        let first = gp.get_all_matches()[0].clone();
        let mut gp2 = gp.clone();
        assert!(gp2.add_tile(&first), "first add");

        let len1_match = gp2
            .get_all_matches()
            .into_iter()
            .find(|pm| pm.len == 1)
            .expect("need len-1 match");
        let mut gp3 = gp2.clone();
        assert!(gp3.add_tile(&len1_match), "second add");

        assert_minimal_witness_roundtrips_for(&gp3, &mi, "hex two-glue with inner");
    }

    #[test]
    fn get_matches_touching_vertex_lazy_matches_eager() {
        let ts: Arc<TileSet<ZZ12>> =
            Arc::new(TileSet::new(
                vec![Rat::try_from(&tiles::spectre()).unwrap()],
            ));
        let mut gp = GrowingPatch::new(Arc::clone(&ts), 0);
        let pm = gp.get_all_matches().into_iter().next().unwrap();
        assert!(gp.add_tile(&pm));
        gp.ensure_candidates_materialized();

        let n = gp.boundary_len();
        for target in 0..n {
            let lazy = {
                let gp2 = gp.clone();
                gp2.get_matches_touching_vertex(target)
            };
            let eager = gp.get_matches_touching_vertex(target);
            let mut lazy_sorted = lazy;
            lazy_sorted.sort_by_key(|pm| (pm.start_a, pm.len, pm.start_b, pm.tile_id));
            let mut eager_sorted = eager;
            eager_sorted.sort_by_key(|pm| (pm.start_a, pm.len, pm.start_b, pm.tile_id));
            assert_eq!(
                lazy_sorted,
                eager_sorted,
                "Mismatch at target={target}: lazy={}, eager={}",
                lazy_sorted.len(),
                eager_sorted.len()
            );
        }
    }

    #[test]
    fn compute_candidates_covering_position_matches_full_enumeration() {
        let ts: Arc<TileSet<ZZ12>> =
            Arc::new(TileSet::new(
                vec![Rat::try_from(&tiles::spectre()).unwrap()],
            ));
        let mi: Arc<MatchTypeIndex<ZZ12>> = Arc::new(MatchTypeIndex::new(Arc::clone(&ts)));
        let mut gp = GrowingPatch::new(Arc::clone(&ts), 0);
        let pm = gp.get_all_matches().into_iter().next().unwrap();
        assert!(gp.add_tile(&pm));

        let all_cands = GrowingPatch::compute_all_candidates(&mi, gp.angles(), gp.edges());
        let n = gp.angles().len();
        let sort_key = |pm: &PatchMatch| (pm.start_a, pm.len, pm.start_b, pm.tile_id);

        for target in 0..n {
            let mut covering = GrowingPatch::compute_candidates_covering_position(
                &mi,
                gp.angles(),
                gp.edges(),
                target,
            );

            // Ground truth: every match in the full enumeration that touches
            // `target`. Compared as multisets via sorting.
            let mut touching_truth: Vec<PatchMatch> = all_cands
                .iter()
                .flatten()
                .filter(|pm| cyclic_range_contains(pm.start_a, pm.len, target, n))
                .cloned()
                .collect();

            covering.sort_by_key(sort_key);
            touching_truth.sort_by_key(sort_key);
            assert_eq!(
                covering, touching_truth,
                "covering vs touching-from-all mismatch at target={target}",
            );
        }
    }

    /// Snapshot of externally observable patch state plus a probe of the
    /// internal spatial grid via candidate accept/reject classification.
    /// Two patches with equal snapshots behave identically against further
    /// `add_tile` attempts — grid corruption would show up as a different
    /// reject set even when angles/edges/etc. are still equal.
    fn classify_candidates<T: IsComplex + IsRingOrField + Units>(
        gp: &GrowingPatch<T>,
    ) -> Vec<(PatchMatch, bool)> {
        let mut results: Vec<(PatchMatch, bool)> = gp
            .get_all_matches()
            .into_iter()
            .map(|pm| {
                let mut trial = gp.clone();
                let ok = trial.add_tile(&pm);
                (pm, ok)
            })
            .collect();
        results.sort_by_key(|(pm, _)| (pm.start_a, pm.len, pm.start_b, pm.tile_id));
        results
    }

    /// Snapshot every publicly observable component of a growing patch,
    /// plus the candidate classification (which doubles as a grid probe).
    #[allow(clippy::type_complexity)]
    fn snapshot_growing<T: IsComplex + IsRingOrField + Units>(
        gp: &GrowingPatch<T>,
    ) -> (
        Vec<i8>,
        Vec<EdgeInfo>,
        Vec<Vec<EdgeInfo>>,
        Vec<usize>,
        usize,
        usize,
        Vec<(PatchMatch, bool)>,
    ) {
        (
            gp.angles().to_vec(),
            gp.edges().to_vec(),
            gp.inner_chains().to_vec(),
            gp.patch_tile_ids().to_vec(),
            gp.next_tile_id(),
            gp.boundary_len(),
            classify_candidates(gp),
        )
    }

    /// The only legitimate rejection path in `add_tile_growing` is the
    /// geometric collision check (`check_segment_clear`) — paths 1, 2, 4
    /// are invariants that legitimate callers (`get_all_matches`) never
    /// violate. This test exercises path 3 against a 2-spectre patch and
    /// asserts that the full patch state (plus a grid probe via candidate
    /// classification) is byte-identical after the failed `add_tile`.
    #[test]
    fn add_tile_failure_leaves_state_unchanged() {
        let ts: Arc<TileSet<ZZ12>> =
            Arc::new(TileSet::new(
                vec![Rat::try_from(&tiles::spectre()).unwrap()],
            ));
        let mut gp = GrowingPatch::new(Arc::clone(&ts), 0);
        let first = gp.get_all_matches().into_iter().next().unwrap();
        assert!(gp.add_tile(&first), "fixture setup");
        let before = snapshot_growing(&gp);
        let failing_pm = before
            .6
            .iter()
            .find(|(_, ok)| !*ok)
            .map(|(pm, _)| pm.clone())
            .expect("expected at least one colliding candidate");
        assert!(
            !gp.add_tile(&failing_pm),
            "must reject a colliding candidate",
        );
        assert_eq!(
            snapshot_growing(&gp),
            before,
            "state changed after a geometrically-rejected pm",
        );
    }

    /// `get_all_matches()` returns edge-compatible candidates without
    /// checking spatial overlap (it only filters via single-edge
    /// compatibility and angle math). For non-convex tiles like spectre,
    /// some of those candidates would self-intersect with existing tiles,
    /// and `add_tile`'s `check_segment_clear` path is the safety net that
    /// catches them. This test pins that behavior: at least one returned
    /// candidate must be rejected, and at least one must be accepted (so
    /// we know the candidate list is non-trivial).
    #[test]
    fn add_tile_rejects_geometrically_invalid_candidate() {
        let ts: Arc<TileSet<ZZ12>> =
            Arc::new(TileSet::new(
                vec![Rat::try_from(&tiles::spectre()).unwrap()],
            ));
        let mut gp = GrowingPatch::new(Arc::clone(&ts), 0);
        let first = gp.get_all_matches().into_iter().next().unwrap();
        assert!(gp.add_tile(&first));

        let candidates = gp.get_all_matches();
        let (mut accepted, mut rejected) = (0usize, 0usize);
        for pm in &candidates {
            let mut trial = gp.clone();
            if trial.add_tile(pm) {
                accepted += 1;
            } else {
                rejected += 1;
            }
        }
        assert!(
            rejected > 0,
            "expected at least one geometrically-invalid candidate to be rejected; \
             all {} candidates were accepted",
            candidates.len()
        );
        assert!(
            accepted > 0,
            "expected at least one valid candidate to be accepted; all {} rejected",
            candidates.len()
        );
    }

    /// Cross-check between the two independent geometric implementations:
    /// GrowingPatch's incremental check (which maintains a UnitSquareGrid
    /// across glues that remove multiple segments and add new ones with an
    /// allowed-endpoint exception) versus Snake's batch validator (which
    /// walks the resulting boundary segment-by-segment from origin and
    /// checks each new segment against the previously visited ones).
    ///
    /// Both ultimately use the same `intersect` + `UnitSquareGrid` primitive
    /// but compose it differently. They must agree on accept/reject for
    /// every candidate.
    ///
    /// Skips candidates that would produce ±hturn (Snake panics on hturn,
    /// and `compute_glue_angles` would have already rejected them at the
    /// add_tile level — the two paths trivially agree there).
    #[test]
    fn add_tile_decision_agrees_with_snake_on_spectre() {
        let ts: Arc<TileSet<ZZ12>> =
            Arc::new(TileSet::new(
                vec![Rat::try_from(&tiles::spectre()).unwrap()],
            ));
        let mut gp = GrowingPatch::new(Arc::clone(&ts), 0);
        let first = gp.get_all_matches().into_iter().next().unwrap();
        assert!(gp.add_tile(&first), "fixture setup");

        let candidates = gp.get_all_matches();
        let tileset = gp.tileset().clone();
        let mut compared = 0usize;
        let mut discrepancies: Vec<(PatchMatch, bool, bool)> = Vec::new();

        for pm in &candidates {
            let new_angles = match compute_glue_angles::<ZZ12>(gp.angles(), pm, &tileset) {
                Some(a) => a,
                None => continue,
            };
            let snake_ok = Snake::<ZZ12>::try_from(new_angles.as_slice()).is_ok();
            let mut trial = gp.clone();
            let gp_ok = trial.add_tile(pm);
            if snake_ok != gp_ok {
                discrepancies.push((pm.clone(), snake_ok, gp_ok));
            }
            compared += 1;
        }

        assert!(compared > 0, "expected non-zero candidates to compare");
        assert!(
            discrepancies.is_empty(),
            "Snake and add_tile disagreed on {} of {} candidates: {:?}",
            discrepancies.len(),
            compared,
            discrepancies
        );
    }

    /// After every successful `add_tile`, the resulting boundary should
    /// be a valid (non-self-intersecting) closed Snake polygon. Spectre
    /// is the right fixture because it has a non-convex shape — most of
    /// the candidate boundaries are non-trivial.
    #[test]
    fn growing_patch_boundary_validates_as_snake_through_growth() {
        let ts: Arc<TileSet<ZZ12>> =
            Arc::new(TileSet::new(
                vec![Rat::try_from(&tiles::spectre()).unwrap()],
            ));
        let mut gp = GrowingPatch::new(Arc::clone(&ts), 0);
        let mut step = 0usize;
        while step < 4 {
            let pm = match gp.get_all_matches().first() {
                Some(pm) => pm.clone(),
                None => break,
            };
            if !gp.add_tile(&pm) {
                break;
            }
            let angles = gp.angles().to_vec();
            let snake = Snake::<ZZ12>::try_from(angles.as_slice());
            assert!(
                snake.is_ok(),
                "step {step}: GrowingPatch's boundary failed Snake validation: angles={angles:?}"
            );
            assert!(
                snake.unwrap().is_closed(),
                "step {step}: GrowingPatch's boundary should close as a polygon"
            );
            step += 1;
        }
        assert!(step > 0, "expected at least one successful add");
    }

    /// Brute-force candidate enumeration independent of `MatchTypeIndex`.
    ///
    /// `compute_all_candidates` (and therefore `get_all_matches`) relies on
    /// the pre-computed `MatchTypeIndex::candidates_starting_at` index for
    /// the segment path, and direct iteration for the junction path. This
    /// test brute-forces every `(tile_id_b, ib, start_a)` triple and
    /// applies the same downstream filters
    /// (`junction_gap_nonnegative`, `try_glue_precomputed`), so a mismatch
    /// against `get_all_matches()` would indicate a bug in either the
    /// index or the segment/junction routing.
    #[test]
    fn get_all_matches_matches_brute_force_on_spectre() {
        let ts: Arc<TileSet<ZZ12>> =
            Arc::new(TileSet::new(
                vec![Rat::try_from(&tiles::spectre()).unwrap()],
            ));
        let mut gp = GrowingPatch::new(Arc::clone(&ts), 0);
        let first = gp.get_all_matches().into_iter().next().unwrap();
        assert!(gp.add_tile(&first), "fixture setup");

        let n = gp.boundary_len();
        let rat = Rat::from_slice_unchecked(gp.angles());

        let mut brute: std::collections::BTreeSet<(usize, usize, usize, usize)> =
            std::collections::BTreeSet::new();
        for tile_id_b in 0..ts.num_tiles() {
            let tile_b = ts.rat(tile_id_b);
            let b_seq = tile_b.seq();
            let m_tile = b_seq.len();
            for ib in 0..m_tile {
                for start_a in 0..n {
                    let (ns, len, ne) = rat.get_match((start_a as i64, ib as i64), tile_b);
                    if len == 0 {
                        continue;
                    }
                    let ns_u = ns.rem_euclid(n as i64) as usize;
                    let ne_u = ne.rem_euclid(m_tile as i64) as usize;
                    if !crate::intgeom::matchtypes::junction_gap_nonnegative(
                        gp.angles(),
                        ns_u,
                        len,
                        b_seq,
                        ne_u,
                    ) {
                        continue;
                    }
                    if rat
                        .try_glue_precomputed((ns, len, ne), tile_b, true)
                        .is_ok()
                    {
                        brute.insert((ns_u, len, ne_u, tile_id_b));
                    }
                }
            }
        }

        let from_api: std::collections::BTreeSet<(usize, usize, usize, usize)> = gp
            .get_all_matches()
            .into_iter()
            .map(|pm| (pm.start_a, pm.len, pm.start_b, pm.tile_id))
            .collect();

        assert_eq!(
            brute, from_api,
            "brute-force candidate set differs from get_all_matches()"
        );
    }

    /// Like `get_all_matches_matches_brute_force_on_spectre` but for
    /// `get_matches_touching_vertex`: brute-force enumerate all matches,
    /// filter by `cyclic_range_contains(start_a, len, v, n)` for each
    /// vertex `v`, and compare against the per-vertex fast path.
    #[test]
    fn get_matches_touching_vertex_matches_brute_force_on_spectre() {
        let ts: Arc<TileSet<ZZ12>> =
            Arc::new(TileSet::new(
                vec![Rat::try_from(&tiles::spectre()).unwrap()],
            ));
        let mut gp = GrowingPatch::new(Arc::clone(&ts), 0);
        let first = gp.get_all_matches().into_iter().next().unwrap();
        assert!(gp.add_tile(&first), "fixture setup");
        gp.ensure_candidates_materialized();

        let n = gp.boundary_len();
        let rat = Rat::from_slice_unchecked(gp.angles());

        let mut brute_matches: Vec<PatchMatch> = Vec::new();
        for tile_id_b in 0..ts.num_tiles() {
            let tile_b = ts.rat(tile_id_b);
            let b_seq = tile_b.seq();
            let m_tile = b_seq.len();
            for ib in 0..m_tile {
                for start_a in 0..n {
                    let (ns, len, ne) = rat.get_match((start_a as i64, ib as i64), tile_b);
                    if len == 0 {
                        continue;
                    }
                    let ns_u = ns.rem_euclid(n as i64) as usize;
                    let ne_u = ne.rem_euclid(m_tile as i64) as usize;
                    if !crate::intgeom::matchtypes::junction_gap_nonnegative(
                        gp.angles(),
                        ns_u,
                        len,
                        b_seq,
                        ne_u,
                    ) {
                        continue;
                    }
                    if rat
                        .try_glue_precomputed((ns, len, ne), tile_b, true)
                        .is_ok()
                    {
                        brute_matches.push(PatchMatch {
                            start_a: ns_u,
                            len,
                            start_b: ne_u,
                            tile_id: tile_id_b,
                        });
                    }
                }
            }
        }
        // Dedup the brute set (the (start_a, ib) double-counts hit the same
        // canonical match).
        let brute_set: std::collections::BTreeSet<(usize, usize, usize, usize)> = brute_matches
            .iter()
            .map(|pm| (pm.start_a, pm.len, pm.start_b, pm.tile_id))
            .collect();

        for target in 0..n {
            let touching_brute: std::collections::BTreeSet<(usize, usize, usize, usize)> =
                brute_set
                    .iter()
                    .copied()
                    .filter(|(start_a, len, _, _)| cyclic_range_contains(*start_a, *len, target, n))
                    .collect();
            let touching_api: std::collections::BTreeSet<(usize, usize, usize, usize)> = gp
                .get_matches_touching_vertex(target)
                .into_iter()
                .map(|pm| (pm.start_a, pm.len, pm.start_b, pm.tile_id))
                .collect();
            assert_eq!(
                touching_brute, touching_api,
                "mismatch at target={target}: brute={touching_brute:?} api={touching_api:?}"
            );
        }
    }

    /// `vertex_type_at(i)` for a `Growing` patch should return
    /// `Some(PatchVertexInfo { angle: angles[i], cw: edges[i-1], ccw: edges[i] })`
    /// for every in-range `i`, and `None` for out-of-range / Seed state.
    #[test]
    fn vertex_type_at_returns_consistent_info() {
        let mut gp = hex_patch();
        let pm = gp
            .get_all_matches()
            .into_iter()
            .find(|p| p.len == 1)
            .expect("len-1 hex match");
        assert!(gp.add_tile(&pm), "fixture");

        let n = gp.boundary_len();
        let angles = gp.angles().to_vec();
        let edges = gp.edges().to_vec();
        for i in 0..n {
            let info = gp.vertex_type_at(i).expect("Some for in-range i");
            assert_eq!(info.angle, angles[i], "angle mismatch at {i}");
            assert_eq!(info.cw, edges[(i + n - 1) % n], "cw edge at {i}");
            assert_eq!(info.ccw, edges[i], "ccw edge at {i}");
        }
        assert!(gp.vertex_type_at(n).is_none(), "out-of-range returns None");
        assert!(
            gp.vertex_type_at(n + 7).is_none(),
            "far out-of-range returns None"
        );

        // Seed state has no boundary; vertex_type_at always returns None.
        let seed = hex_patch();
        assert!(!seed.is_growing());
        assert!(seed.vertex_type_at(0).is_none(), "seed has no boundary");
    }

    /// `neighbor_junction_offsets(pos)` returns offsets into the CW and CCW
    /// neighbouring junctions' tile sequences. The returned values must
    /// (a) be within the relevant tile's length and (b) correctly identify
    /// the CW junction's edge and the (ccw_prev + 1) offset of the CCW
    /// junction's preceding edge.
    #[test]
    fn neighbor_junction_offsets_returns_valid_offsets() {
        let mut gp = hex_patch();
        let pm = gp
            .get_all_matches()
            .into_iter()
            .find(|p| p.len == 1)
            .expect("len-1 hex match");
        assert!(gp.add_tile(&pm), "fixture");
        let n = gp.boundary_len();
        let edges = gp.edges().to_vec();
        let ts = gp.tileset().clone();

        for pos in 0..n {
            let (cw_off, ccw_off) = gp
                .neighbor_junction_offsets(pos)
                .expect("Some for valid pos");

            // Walk CW to the nearest junction (possibly == pos itself).
            let mut j_cw = (pos + n - 1) % n;
            while j_cw != pos && !gp.is_junction(j_cw) {
                j_cw = (j_cw + n - 1) % n;
            }
            let cw_tile_len = ts.rat(edges[j_cw].tile_id).len();
            assert!(cw_off < cw_tile_len, "cw_off out of range at pos {pos}");
            assert_eq!(
                cw_off, edges[j_cw].tile_offset,
                "cw_off should be the CW junction's tile_offset at pos {pos}",
            );

            // Walk CCW to the nearest junction.
            let mut j_ccw = (pos + 1) % n;
            while j_ccw != pos && !gp.is_junction(j_ccw) {
                j_ccw = (j_ccw + 1) % n;
            }
            let ccw_prev_edge = edges[(j_ccw + n - 1) % n];
            let ccw_tile_len = ts.rat(ccw_prev_edge.tile_id).len();
            assert!(ccw_off < ccw_tile_len, "ccw_off out of range at pos {pos}");
            assert_eq!(
                ccw_off,
                (ccw_prev_edge.tile_offset + 1) % ccw_tile_len,
                "ccw_off should be (ccw_prev edge's offset + 1) at pos {pos}",
            );
        }

        // Out-of-range and Seed state both return None.
        assert!(gp.neighbor_junction_offsets(n).is_none());
        let seed = hex_patch();
        assert!(seed.neighbor_junction_offsets(0).is_none());
    }

    /// `tile_segments()` should:
    /// (a) cover the boundary contiguously (segments concatenated from
    /// 0 to `n` with no gaps),
    /// (b) have segment boundaries at exactly the junction positions,
    /// (c) within each segment, `tile_id` is constant and `tile_offset`
    /// advances by 1 modulo the tile's edge count.
    #[test]
    fn tile_segments_partitions_boundary() {
        let mut gp = hex_patch();
        let pm = gp
            .get_all_matches()
            .into_iter()
            .find(|p| p.len == 1)
            .expect("len-1 hex match");
        assert!(gp.add_tile(&pm), "fixture");
        let n = gp.boundary_len();
        let edges = gp.edges().to_vec();
        let segs = gp.tile_segments();

        // Contiguous partition.
        assert_eq!(
            segs.first().map(|s| s.patch_start),
            Some(0),
            "first segment starts at 0"
        );
        assert_eq!(
            segs.last().map(|s| s.patch_end),
            Some(n),
            "last segment ends at n"
        );
        for w in segs.windows(2) {
            assert_eq!(
                w[0].patch_end, w[1].patch_start,
                "segments must be contiguous"
            );
        }

        // Consistent tile_id and contiguous offsets within each segment.
        for seg in &segs {
            let tile_id = seg.tile_id;
            let tile_len = gp.tileset().rat(tile_id).len();
            for k in 0..(seg.patch_end - seg.patch_start) {
                let pos = seg.patch_start + k;
                assert_eq!(edges[pos].tile_id, tile_id, "tile_id at pos {pos}");
                assert_eq!(
                    edges[pos].tile_offset,
                    (seg.offset_start + k) % tile_len,
                    "tile_offset at pos {pos}",
                );
            }
        }

        // A position is a segment start iff it is position 0 (the
        // linear-partition seam, always present) or a junction.
        let expected_starts: std::collections::BTreeSet<usize> = std::iter::once(0)
            .chain((0..n).filter(|&i| gp.is_junction(i)))
            .collect();
        let actual_starts: std::collections::BTreeSet<usize> =
            segs.iter().map(|s| s.patch_start).collect();
        assert_eq!(
            actual_starts, expected_starts,
            "segment starts must equal {{0}} ∪ junctions"
        );
    }

    #[test]
    fn construct_minimal_witness_hex_boundary_matches_brute_force() {
        let gp = hex_patch();
        let mi = gp.match_index().clone();
        for pm in &gp.get_all_matches() {
            let mut brute = gp.clone();
            assert!(brute.add_tile(pm), "brute glue");
            assert_witness_matches_brute_force(&brute, &mi, &format!("hex pm {:?}", pm));
        }
    }

    #[test]
    fn construct_minimal_witness_square_boundary_matches_brute_force() {
        let gp = square_patch();
        let mi = gp.match_index().clone();
        for pm in &gp.get_all_matches() {
            let mut brute = gp.clone();
            assert!(brute.add_tile(pm), "brute glue");
            assert_witness_matches_brute_force(&brute, &mi, &format!("square pm {:?}", pm));
        }
    }

    #[test]
    fn construct_minimal_witness_spectre_roundtrip() {
        let ts: Arc<TileSet<ZZ12>> =
            Arc::new(TileSet::new(
                vec![Rat::try_from(&tiles::spectre()).unwrap()],
            ));
        let mut gp = GrowingPatch::new(Arc::clone(&ts), 0);
        let pm = gp.get_all_matches().into_iter().next().unwrap();
        assert!(gp.add_tile(&pm), "first spectre glue");
        let mi = gp.match_index().clone();
        assert_minimal_witness_roundtrips_for(&gp, &mi, "spectre first-glue");
    }

    #[test]
    fn forward_match_length_hex_basic() {
        let hex: Snake<ZZ12> = tiles::hexagon();
        let rat = Rat::try_from(&hex).unwrap();
        let seq = rat.seq();

        assert_eq!(forward_match_length(seq, 0, seq, 0), 1);
        assert_eq!(forward_match_length(seq, 3, seq, 3), 1);
        assert_eq!(forward_match_length(seq, 0, seq, 1), 1);

        let boundary: Vec<i8> = vec![-2, 2, 2, 2, 2, -2, 2, 2, 2, 2];
        assert_eq!(forward_match_length(&boundary, 5, seq, 0), 1);
        assert_eq!(forward_match_length(&boundary, 0, seq, 0), 1);
    }

    #[test]
    fn forward_match_length_square_basic() {
        let sq: Snake<ZZ4> = tiles::square();
        let rat = Rat::try_from(&sq).unwrap();
        let seq = rat.seq();

        assert_eq!(forward_match_length(seq, 0, seq, 0), 1);
        assert_eq!(forward_match_length(seq, 2, seq, 2), 1);
    }

    #[test]
    fn glue_raw_angles_hex_self_glue() {
        let hex: Snake<ZZ12> = tiles::hexagon();
        let rat = Rat::try_from(&hex).unwrap();
        let seq = rat.seq().to_vec();

        let result = angles::glue_raw_angles::<ZZ12>(&seq, &seq, 0, 1, 0);
        assert!(result.is_some());
        let gr = result.unwrap();
        assert_eq!(gr.angles.len(), 10);
        assert_eq!(gr.a_yx, Some(-2));
        assert_eq!(gr.a_xy, Some(-2));
    }

    #[test]
    fn glue_raw_angles_matches_rat_glue() {
        let hex: Snake<ZZ12> = tiles::hexagon();
        let rat = Rat::try_from(&hex).unwrap();
        let seq = rat.seq();

        let rat_result = rat.try_glue((0, 0), &rat).expect("rat glue");
        let raw_result = angles::glue_raw_angles::<ZZ12>(seq, seq, 0, 1, 0).expect("raw glue");

        // Both glue paths must produce the same boundary up to cyclic
        // rotation (they may pick different starting positions).
        assert_same_cyclic_shape(rat_result.seq(), &raw_result.angles, "rat vs raw glue");
    }

    #[test]
    fn test_junction_angle_sequence_hex() {
        let gp = hex_patch();
        let mi = gp.match_index().clone();
        for pm in &gp.get_all_matches() {
            let mut glued = gp.clone();
            assert!(glued.add_tile(pm), "glue");
            assert_junction_angle_sequence_valid(&glued, &mi, &format!("hex pm {:?}", pm));
        }
    }

    #[test]
    fn construct_witness_from_vt_sequence_single_vt_roundtrip() {
        let gp = hex_patch();
        let mi = gp.match_index().clone();

        let mut gp = gp;
        let pm = gp
            .get_all_matches()
            .into_iter()
            .find(|pm| pm.len == 1)
            .expect("len-1 match");
        assert!(gp.add_tile(&pm), "first glue");

        let vt = gp.junction_vertex_type_at(0).expect("junction at 0");

        let (minimal, _wpos) =
            GrowingPatch::construct_minimal_witness(&vt, mi.clone()).expect("minimal witness");

        let (reconstructed, _junc_positions) =
            GrowingPatch::construct_witness_from_vt_sequence(std::slice::from_ref(&vt), mi)
                .expect("reconstruction");

        // construct_minimal_witness delegates to
        // construct_witness_from_vt_sequence for single-element input,
        // so the two outputs must be byte-identical, not just congruent.
        assert_eq!(minimal.angles(), reconstructed.angles());
        assert_eq!(minimal.edges(), reconstructed.edges());
        assert_eq!(minimal.inner_chains(), reconstructed.inner_chains());
    }

    /// Build a 5-hexagon plus-shaped cross: a central hex with four
    /// hexagons attached on alternating sides. The resulting patch has
    /// 18 boundary edges and 6 junctions arranged symmetrically.
    ///
    /// Construction picks each glue by *resulting boundary length*:
    /// start with a bi-hex (10 edges), then grow to 14 → 16 → 18.
    /// `five_hex_cross_structure` verifies the boundary symmetry,
    /// junction count, and tile-id pattern.
    ///
    /// FIXME: replace with explicit tile-placement construction once
    /// that API exists. The current recipe depends on which pm is
    /// returned first by `get_all_matches()` for each target
    /// boundary_len; if two pms produced the same length and the
    /// iteration order changed, we'd silently build a different
    /// (equally-shaped) cross — caught by the structure test, but
    /// avoidable with a direct geometry-based fixture.
    fn five_hex_cross() -> GrowingPatch<ZZ12> {
        let mut gp = hex_patch();
        let first = gp
            .get_all_matches()
            .into_iter()
            .find(|pm| pm.len == 1)
            .expect("seed has len-1 match");
        assert!(gp.add_tile(&first), "glue tile 2");
        assert_eq!(gp.boundary_len(), 10, "bi-hex should have 10 edges");

        let target_sequence = [14usize, 16, 18];
        for (step, &target_n) in target_sequence.iter().enumerate() {
            let matches = gp.get_all_matches();
            let mut found = false;
            for pm in &matches {
                let mut trial = gp.clone();
                if !trial.add_tile(pm) {
                    continue;
                }
                if trial.boundary_len() == target_n {
                    assert!(gp.add_tile(pm), "glue failed at step {}", step + 2);
                    assert_eq!(
                        gp.boundary_len(),
                        target_n,
                        "step {}: boundary_len should be {target_n} after glue",
                        step + 2
                    );
                    found = true;
                    break;
                }
            }
            assert!(
                found,
                "no match giving boundary_len={target_n} at step {}",
                step + 2
            );
        }
        gp
    }

    #[test]
    fn five_hex_cross_structure() {
        let gp = five_hex_cross();
        let n = gp.boundary_len();
        assert_eq!(n, 18);

        let angles = gp.angles();
        assert_eq!(&angles[..9], &angles[9..], "boundary should be symmetric");

        let junctions: Vec<usize> = (0..n).filter(|&i| gp.is_junction(i)).collect();
        assert_eq!(junctions.len(), 6);

        let mut segs: Vec<usize> = Vec::new();
        for w in junctions.windows(2) {
            segs.push(w[1] - w[0]);
        }
        segs.push(n - junctions[5] + junctions[0]);
        assert_eq!(
            segs,
            vec![1, 4, 4, 1, 4, 4],
            "junction offsets should be 1,4,4,1,4,4"
        );

        for i in 0..n {
            let prev = (i + n - 1) % n;
            let id = gp.patch_tile_ids()[i];
            let prev_id = gp.patch_tile_ids()[prev];
            if gp.is_junction(i) {
                assert_ne!(
                    id, prev_id,
                    "junction at {i} should have distinct patch_tile_ids"
                );
            }
        }

        let mut run_start = 0;
        let mut runs: Vec<(usize, usize)> = Vec::new();
        for i in 1..=n {
            if i == n || gp.patch_tile_ids()[i] != gp.patch_tile_ids()[run_start] {
                runs.push((gp.patch_tile_ids()[run_start], i - run_start));
                run_start = i;
            }
        }
        assert_eq!(runs.len(), 6, "should have 6 runs of patch_tile_ids");
        let center_runs: Vec<&(usize, usize)> = runs.iter().filter(|(id, _)| *id == 0).collect();
        assert_eq!(
            center_runs.len(),
            2,
            "center tile should appear in exactly 2 runs"
        );
        assert_eq!(center_runs[0].1, 1, "each center run should be 1 edge");
        assert_eq!(center_runs[1].1, 1, "each center run should be 1 edge");
    }

    #[test]
    fn reconstruct_five_hex_cross() {
        let gp = five_hex_cross();
        let mi = gp.match_index().clone();

        let n = gp.boundary_len();
        let mut vt_seq: Vec<OpenVertexType> = Vec::new();
        for i in 0..n {
            if gp.is_junction(i) {
                let vt = gp.junction_vertex_type_at(i).unwrap();
                assert!(
                    vt.inner.is_empty(),
                    "hex boundary junctions should have empty inner"
                );
                vt_seq.push(vt);
            }
        }
        assert_eq!(vt_seq.len(), 6);

        let result = GrowingPatch::construct_witness_from_vt_sequence(&vt_seq, mi);
        let (reconstructed, _junc_positions) = result.expect("reconstruction should succeed");

        assert_same_cyclic_shape(
            gp.angles(),
            reconstructed.angles(),
            "5-hex-cross: reconstructed vs original",
        );
        assert_eq!(
            reconstructed.boundary_len(),
            gp.boundary_len(),
            "boundary length should match"
        );

        let recon_juncs: Vec<usize> = (0..reconstructed.boundary_len())
            .filter(|&i| reconstructed.is_junction(i))
            .collect();
        assert_eq!(recon_juncs.len(), 6, "should have 6 junctions");
    }

    #[test]
    fn next_junction_on_raw_boundary_finds_all_junctions() {
        let gp = hex_patch();
        let ts = gp.tileset().clone();

        let mut gp = gp;
        let pm = gp
            .get_all_matches()
            .into_iter()
            .find(|pm| pm.len == 1)
            .expect("len-1 match");
        assert!(gp.add_tile(&pm), "first glue");

        let raw = RawBoundary {
            angles: gp.angles().to_vec(),
            edges: gp.edges().to_vec(),
            inner_chains: gp.inner_chains().to_vec(),
            patch_tile_ids: gp.patch_tile_ids().to_vec(),
        };

        let n = raw.angles.len();
        let mut junctions: Vec<usize> = Vec::new();
        for i in 0..n {
            if raw_is_junction(&raw, ts.as_ref(), i) {
                junctions.push(i);
            }
        }
        assert_eq!(junctions.len(), 2, "two-hex should have 2 junctions");

        let j1 = next_junction_on_raw_boundary(&raw, ts.as_ref(), junctions[0])
            .expect("should find next junction");
        assert_eq!(j1, junctions[1], "should find the other junction");

        let j0 = next_junction_on_raw_boundary(&raw, ts.as_ref(), junctions[1])
            .expect("should wrap around");
        assert_eq!(j0, junctions[0], "should wrap to first junction");
    }

    #[test]
    fn test_junction_angle_sequence_square() {
        let gp = square_patch();
        let mi = gp.match_index().clone();
        for pm in &gp.get_all_matches() {
            let mut glued = gp.clone();
            assert!(glued.add_tile(pm), "glue");
            assert_junction_angle_sequence_valid(&glued, &mi, &format!("square pm {:?}", pm));
        }
    }

    #[test]
    fn normalize_five_hex_cross() {
        let gp = five_hex_cross();
        let mut gp2 = gp.clone();
        gp2.normalize();

        assert_eq!(gp2.boundary_len(), 18);

        let ptids = gp2.patch_tile_ids();
        let mut seen = std::collections::HashSet::new();
        for &id in ptids {
            seen.insert(id);
        }
        let max_id = *seen.iter().max().unwrap();
        assert_eq!(
            seen.len(),
            max_id + 1,
            "ptids should be 0..=max with no gaps"
        );
        assert_eq!(gp2.next_tile_id(), seen.len());

        let angles = gp2.angles();
        let min_angle = *angles.iter().min().unwrap();
        assert_eq!(
            angles[0], min_angle,
            "normalized boundary should start at lex-min angle"
        );
    }

    #[test]
    fn normalize_idempotent() {
        let gp = five_hex_cross();
        let mut gp1 = gp.clone();
        gp1.normalize();
        let snap1 = (
            gp1.angles().to_vec(),
            gp1.edges().to_vec(),
            gp1.patch_tile_ids().to_vec(),
        );
        gp1.normalize();
        let snap2 = (
            gp1.angles().to_vec(),
            gp1.edges().to_vec(),
            gp1.patch_tile_ids().to_vec(),
        );
        assert_eq!(snap1, snap2, "normalize should be idempotent");
    }

    fn t_tetromino_angles() -> Vec<i8> {
        let snake: Snake<ZZ4> = tiles::tetromino_T();
        let rat = Rat::try_from(&snake).unwrap();
        rat.seq().to_vec()
    }

    /// Build the T-tetromino — 4 unit squares in a T shape:
    ///
    /// ```text
    ///     +---+
    ///     |   |
    /// +---+   +---+
    /// |             |
    /// +---+---+---+
    /// ```
    ///
    /// Glues three squares onto the seed by picking each candidate
    /// match by its `(start_a, len, start_b)` triple. Per-step boundary
    /// length assertions (square seed: 4 edges; bi-square: 6; tri-square:
    /// 8; T: 10) catch the case where the match-finder semantics drift
    /// such that the selected triple produces a different shape.
    ///
    /// FIXME: replace with explicit tile-placement construction once
    /// that API exists. The current recipe is opaque: a reader has to
    /// reverse-engineer the intended shape from three triples.
    fn t_tetromino() -> GrowingPatch<ZZ4> {
        let mut gp = square_patch();
        // Glue 1: square + square along seed's edge 0 → bi-square (6 edges).
        let ms: Vec<_> = gp.get_all_matches();
        let pm1 = ms
            .iter()
            .find(|pm| pm.start_a == 0 && pm.len == 1 && pm.start_b == 0)
            .expect("first match");
        assert!(gp.add_tile(pm1), "glue 1");
        assert_eq!(gp.boundary_len(), 6, "bi-square should have 6 edges");

        // Glue 2: third square along the bi-square's edge 0 → tri-square (8 edges).
        let ms: Vec<_> = gp.get_all_matches();
        let pm2 = ms
            .iter()
            .find(|pm| pm.start_a == 0 && pm.len == 1 && pm.start_b == 1)
            .expect("second match");
        assert!(gp.add_tile(pm2), "glue 2");
        assert_eq!(gp.boundary_len(), 8, "tri-square should have 8 edges");

        // Glue 3: fourth square to form the T-stem → T-tetromino (10 edges).
        let ms: Vec<_> = gp.get_all_matches();
        let pm3 = ms
            .iter()
            .find(|pm| pm.start_a == 0 && pm.len == 1 && pm.start_b == 1)
            .expect("third match");
        assert!(gp.add_tile(pm3), "glue 3");
        assert_eq!(gp.boundary_len(), 10, "T-tetromino should have 10 edges");
        gp
    }

    #[test]
    fn reconstruct_t_tetromino() {
        let gp = t_tetromino();
        let mi = gp.match_index().clone();
        let n = gp.boundary_len();
        assert_eq!(n, 10);

        let ref_angles = t_tetromino_angles();
        assert_same_cyclic_shape(
            gp.angles(),
            &ref_angles,
            "built patch should be the T tetromino shape",
        );

        let mut vt_seq: Vec<OpenVertexType> = Vec::new();
        for i in 0..n {
            if gp.is_junction(i) {
                let vt = gp.junction_vertex_type_at(i).unwrap();
                vt_seq.push(vt);
            }
        }
        assert!(!vt_seq.is_empty(), "T should have junctions");

        let has_inner = vt_seq.iter().any(|vt| !vt.inner.is_empty());
        assert!(
            has_inner,
            "T tetromino should have junctions with non-empty inner"
        );

        let result = GrowingPatch::construct_witness_from_vt_sequence(&vt_seq, mi);
        let (reconstructed, _junc_positions) = result.expect("reconstruction should succeed");

        assert_eq!(
            reconstructed.boundary_len(),
            n,
            "boundary length should match"
        );

        assert_same_cyclic_shape(
            reconstructed.angles(),
            &ref_angles,
            "reconstructed angles should match T",
        );

        let recon_juncs: Vec<usize> = (0..reconstructed.boundary_len())
            .filter(|&i| reconstructed.is_junction(i))
            .collect();
        assert_eq!(
            recon_juncs.len(),
            vt_seq.len(),
            "junction count should match"
        );
    }
}
