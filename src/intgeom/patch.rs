use rustc_hash::FxHashSet;
use std::sync::Arc;

use crate::cyclotomic::geometry::intersect_unit_segments;
use crate::cyclotomic::{IsRingOrField, Units};
use crate::intgeom::angles;
use crate::intgeom::grid::UnitSquareGrid;
use crate::intgeom::matchtypes::{
    is_single_edge_candidate, junction_gap_nonnegative, CandidateMatch, MatchTypeIndex,
};
use crate::intgeom::rat::Rat;
use crate::intgeom::snake::Snake;
use crate::intgeom::tileset::TileSet;

/// Identifies a single tile edge by `(tile_id, tile_offset)`.
///
/// `tile_id` indexes into the patch's [`TileSet`]; `tile_offset` is the
/// edge index within that tile's boundary (`0..tile.len()`).
///
/// On a [`GrowingPatch`] boundary, an `EdgeInfo` at position `i` says
/// "the boundary edge at position `i` is edge `tile_offset` of an
/// instance of tile `tile_id`". Different boundary positions can share
/// the same `EdgeInfo` if the same tile shape appears multiple times
/// in the patch (each tile *instance* gets a separate `patch_tile_id`,
/// tracked elsewhere in `GrowingPatch`).
#[derive(
    Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, serde::Serialize, serde::Deserialize,
)]
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
#[allow(dead_code)] // only used in patch.rs tests today; kept as crate-internal API
pub(crate) struct PatchVertexInfo {
    pub(crate) angle: i8,
    pub(crate) cw: EdgeInfo,
    pub(crate) ccw: EdgeInfo,
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
pub(crate) struct TileSegment {
    pub(crate) patch_start: usize,
    pub(crate) patch_end: usize,
    pub(crate) tile_id: usize,
    pub(crate) offset_start: usize,
}

/// A description of a candidate glue between the current patch
/// boundary and a new tile.
///
/// Returned by [`GrowingPatch::get_all_matches`] and consumed by
/// [`GrowingPatch::add_tile`].
///
/// # Edge-offset convention
///
/// * `start_a` is the first **matched** boundary position on the
///   patch side: the match consumes edges `start_a, start_a+1, ...,
///   start_a+len-1` (modulo current boundary length).
/// * `start_b` is the first **surviving** edge on the new tile's
///   side, just past the match. The new tile's surviving edges run
///   `start_b, start_b+1, ..., start_b + (tile_len - len) - 1`
///   (modulo `tile_len`).
/// * `tile_id` indexes into the patch's [`TileSet`] — the shape of
///   the new tile.
///
/// (Same convention as the lower-level [`MatchType`](crate::intgeom::matchtypes::MatchType).)
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
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
/// In contrast, a [`ClosedVertexType`] describes a fully-surrounded
/// interior vertex, where the petals cover the entire turn and the
/// vertex no longer lies on any boundary.
#[derive(
    Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, serde::Serialize, serde::Deserialize,
)]
pub struct OpenVertexType {
    pub cw: EdgeInfo,
    pub inner: Vec<EdgeInfo>,
    pub ccw: EdgeInfo,
}

/// A **closed** junction vertex: a fully-surrounded interior vertex
/// whose petals cover the entire turn (no boundary gap).
///
/// Stored as the cyclic CCW sequence of [`EdgeInfo`]s around the
/// vertex, canonicalised to its lex-min cyclic rotation so two
/// configurations that differ only by which petal is listed first
/// compare equal. Two consecutive entries (cyclically) flank a single
/// incident tile-instance at the vertex.
///
/// # Construction
///
/// Built from the petal sequence of an [`OpenVertexType`] at the
/// moment a closing glue (`TransitionSide::Both` — both incident
/// boundary edges consumed by the same match) seals the focus
/// vertex: the closed petal ring is `[cw, inner[0], ..., inner[m-1],
/// ccw]`, with the new tile contributing implicitly between `ccw`
/// and `cw` cyclically. See [`Self::from_open_via_closure`].
///
/// # Identity caveat
///
/// Two distinct *new tiles* that close the same `OpenVertexType` and
/// produce the same cyclic edge ring compare equal under this type —
/// the new tile's identity is implicit in the `(ccw, cw)` cyclic
/// adjacency. If the call site cares about which closing tile was
/// used, it must track that separately (the transition that produced
/// the closure already records `tile_id` and `tile_offset`).
#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct ClosedVertexType {
    edges: Vec<EdgeInfo>,
}

impl ClosedVertexType {
    /// Build from a raw cyclic edge sequence; the result is stored in
    /// lex-min cyclic rotation.
    pub fn from_cyclic(petals: &[EdgeInfo]) -> Self {
        ClosedVertexType {
            edges: canonical_cyclic_rotation(petals),
        }
    }

    /// Build the closed VT realised when the given open VT is sealed
    /// by a closing glue. The petal ring is `[cw, inner..., ccw]`
    /// (CCW around the focus); the new tile sits implicitly between
    /// `ccw` and `cw` cyclically.
    pub fn from_open_via_closure(open: &OpenVertexType) -> Self {
        let mut petals = Vec::with_capacity(open.inner.len() + 2);
        petals.push(open.cw);
        petals.extend(open.inner.iter().copied());
        petals.push(open.ccw);
        Self::from_cyclic(&petals)
    }

    /// The petals, in CCW order, starting at the lex-min rotation.
    pub fn edges(&self) -> &[EdgeInfo] {
        &self.edges
    }

    /// Number of petals (= number of incident tile-instances at the
    /// vertex, since each pair of consecutive petals cyclically
    /// flanks one tile).
    pub fn len(&self) -> usize {
        self.edges.len()
    }

    pub fn is_empty(&self) -> bool {
        self.edges.is_empty()
    }
}

/// Return the lex-min cyclic rotation of `seq`. Compares rotations
/// lexicographically by their full `n`-element unrolling.
pub(crate) fn canonical_cyclic_rotation<T: Ord + Clone>(seq: &[T]) -> Vec<T> {
    let n = seq.len();
    if n <= 1 {
        return seq.to_vec();
    }
    let mut best = 0usize;
    for start in 1..n {
        // Compare rotation `start` vs current best.
        let mut cmp = std::cmp::Ordering::Equal;
        for i in 0..n {
            let a = &seq[(start + i) % n];
            let b = &seq[(best + i) % n];
            cmp = a.cmp(b);
            if cmp != std::cmp::Ordering::Equal {
                break;
            }
        }
        if cmp == std::cmp::Ordering::Less {
            best = start;
        }
    }
    (0..n).map(|i| seq[(best + i) % n].clone()).collect()
}

/// A coarser variant of [`OpenVertexType`] that ignores the interior
/// tile arrangement at the junction and instead records the
/// **boundary angle** at the vertex.
///
/// At a junction vertex, the boundary angle (= what the patch stores
/// in its `angles` vec) is `full_turn − sum(incident interior tile
/// angles)`. It captures the boundary's local "shape" at the junction
/// — what matters for future tile attachments — without depending on
/// which specific tiles are buried in the interior.
///
/// Two junctions with the same `(cw_edge, ccw_edge, angle)` admit the
/// same outgoing tile attachments, even if their interior populations
/// differ. This makes `CoarseJunction` the right equivalence for seg
/// enumeration: it conflates configurations whose only difference is
/// invisible from outside the boundary.
#[derive(
    Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, serde::Serialize, serde::Deserialize,
)]
pub struct CoarseJunction {
    pub cw_edge: EdgeInfo,
    pub ccw_edge: EdgeInfo,
    pub angle: i8,
}

/// Which incident edge(s) of a focus junction vertex a transition
/// consumed.
///
/// Used by both VT and NT enumeration. At an open vertex with two
/// incident boundary edges (CW and CCW), a glue can:
///
/// - consume only the CW edge → [`TransitionSide::Cw`],
/// - consume only the CCW edge → [`TransitionSide::Ccw`],
/// - consume **both** incident edges, sealing the vertex →
///   [`TransitionSide::Both`].
///
/// When the side is `Both`, the destination of the transition is
/// always the closed-vertex sentinel (`vertextypes::CLOSED_ID`); the
/// transition's `tile_offset` is normalized to the **CW edge's**
/// matching tile-offset by convention.
///
/// For NT transitions (in `neighborhood`), only `Cw` and `Ccw` arise
/// — NT enumeration steps the central tile along one direction at a
/// time and never produces a `Both` step.
#[derive(
    Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, serde::Serialize, serde::Deserialize,
)]
pub enum TransitionSide {
    /// The glue consumed only the CW edge of the focus vertex.
    Cw,
    /// The glue consumed only the CCW edge of the focus vertex.
    Ccw,
    /// The glue consumed both incident edges, sealing the focus
    /// vertex. The transition's `tile_offset` is canonicalized to
    /// the CW edge.
    Both,
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
pub struct GrowingPatch<T: IsRingOrField> {
    match_index: Arc<MatchTypeIndex<T>>,
    state: PatchState<T>,
}

#[derive(Clone)]
#[allow(clippy::large_enum_variant)]
enum PatchState<T: IsRingOrField> {
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
fn is_junction_at<T: IsRingOrField + Units>(
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
fn vertex_type_raw_from(
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

fn update_inner_chains(
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

    let consumed_ccw_ptid = old_ptids[pm.start_a];
    let cw_survivor_ptid = old_ptids[(pm.start_a + n - 1) % n];
    let mut chain_cw: Vec<EdgeInfo> = old_inner[pm.start_a].clone();
    if consumed_ccw_ptid != cw_survivor_ptid {
        chain_cw.push(old_edges[pm.start_a]);
    }

    // Keystone case: when the petal contributes zero surviving edges
    // (`seg_len_new == 0`, hence `new_n == seg_len_old`), its CW and
    // CCW match endpoints collapse to a single new boundary vertex at
    // index 0. Both chains apply there; concat with chain_cw first
    // (going CW boundary edge → interior → CCW boundary edge passes
    // through the pm.start_a-side first, then the ccw_pos-side).
    // The petal's own perimeter contributions are not represented in
    // inner_chains (same convention as the normal-glue path which
    // also only records consumed-adjacency tile edges).
    if new_n == seg_len_old {
        debug_assert_eq!(seg_len_old, new_n, "keystone glue invariant");
        let mut merged = chain_cw;
        merged.extend(chain_ccw);
        new_inner[0] = merged;
    } else {
        new_inner[0] = chain_ccw;
        new_inner[seg_len_old] = chain_cw;
    }

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
struct RawBoundary {
    angles: Vec<i8>,
    edges: Vec<EdgeInfo>,
    inner_chains: Vec<Vec<EdgeInfo>>,
    patch_tile_ids: Vec<usize>,
}

#[cfg(test)]
fn raw_is_junction<T: IsRingOrField + Units>(
    boundary: &RawBoundary,
    tileset: &TileSet<T>,
    pos: usize,
) -> bool {
    is_junction_at(&boundary.angles, &boundary.edges, tileset, pos)
}

#[cfg(test)]
fn next_junction_on_raw_boundary<T: IsRingOrField + Units>(
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
struct RawGlueResult {
    boundary: RawBoundary,
    new_junc_pos: usize,
    match_len: usize,
    old_survivor_len: usize,
    new_tile_survivor_len: usize,
}

/// Glue one tile to a raw boundary at a boundary junction.
///
/// The `target.tile_offset` is the tile-side junction offset used as the
/// `start_b`/`other_end` argument to `glue_raw_angles`. Surviving target edges
/// are therefore `tile_offset, tile_offset + 1, ...`, modulo the tile edge
/// count. This helper is intentionally shared by VT witness construction and NT
/// enumeration so that edge-offset conventions cannot diverge again.
fn glue_tile_to_raw_boundary<T: IsRingOrField + Units>(
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
/// # Caller contract
///
/// Inherits the precondition from [`angles::glue_raw_angles`]: `pm`
/// must describe a **real** match. The sibling helper
/// [`glue_tile_to_raw_boundary`] satisfies this by computing
/// `pm.len = forward_match_length(...)` (canonical max) before
/// constructing the `PatchMatch`. Callers passing `pm` directly are
/// responsible for ensuring the same invariant — there is no
/// post-glue Snake/grid check on this path (the result is a
/// `RawBoundary` consumed by [`GrowingPatch::from_parts`], which
/// trusts the caller).
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
fn glue_match_to_raw_boundary<T: IsRingOrField + Units>(
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
fn junction_angle_sequence<T: IsRingOrField + Units>(
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
fn flower_petal_glue<T: IsRingOrField + Units>(
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

impl<T: IsRingOrField + Units> GrowingPatch<T> {
    /// Construct a fresh `GrowingPatch` seeded with one tile from
    /// `tileset` (the tile at `seed_tile_id`).
    ///
    /// The patch starts in the `Seed` state — no boundary yet, only
    /// the seed tile. The first call to [`Self::add_tile`] transitions
    /// to `Growing` and produces the initial boundary.
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
    ///
    /// # Precondition 1: vtypes is a LINEAR sequence, not a cyclic description
    ///
    /// `vtypes` must describe a **partial / open** boundary segment, not a
    /// cyclic description of an entire closed boundary. The construction
    /// glues tiles one at a time at each successive junction without
    /// recognising any closure between `vtypes[n-1]` and `vtypes[0]`. If
    /// the supplied sequence happens to encode a closed loop, the chain
    /// will keep adding tiles past where the loop would close and produce
    /// a geometrically invalid (self-intersecting) patch.
    ///
    /// For example: the 6 outer-corner junctions of a 7-hex full corona
    /// form a cyclic boundary trace (= they describe the full closed
    /// outer perimeter). Passing those 6 vtypes to this function does
    /// **not** rebuild the corona's 18-edge boundary — it produces a
    /// 30-edge self-intersecting spiral of 7 tiles, with the 7th tile
    /// overlapping the seed.
    ///
    /// Valid use cases pass partial vt_seqs: phase-1 NT matched-segment
    /// junctions, single open-VT witnesses, and similar.
    ///
    /// # Precondition 2: vtypes must capture every tile of the realising patch
    ///
    /// Each [`OpenVertexType`] carries the tile edges meeting at one
    /// boundary junction via its `cw`, `inner`, and `ccw` fields. This
    /// covers every tile *incident with the boundary* at that vertex.
    /// **Tiles that are fully interior to the patch — i.e., have no
    /// vertex on the boundary, so they appear in no junction's `inner`
    /// field — are not captured by any vtype and will not be placed by
    /// this function.**
    ///
    /// For the 7-hex full corona, this also bites: the central hex is
    /// fully interior with respect to the outer boundary and is in no
    /// outer junction's `inner`. Without the central's geometric
    /// constraint, `flower_petal_glue` lays the outer hexes out as a
    /// curving chain rather than a corona — the chain spirals inward
    /// instead of maintaining the spacing the central would force.
    ///
    /// # Summary of caller obligations
    ///
    /// - `vtypes` describes a partial / open boundary segment.
    /// - Every tile of the realising patch has at least one vertex on
    ///   the boundary segment described by `vtypes`, captured via
    ///   `cw` / `ccw` / `inner` of some entry.
    ///
    /// Both conditions hold for minimal-witness uses (one open-VT, a
    /// chain of phase-1 NT matched-segment junctions with the central
    /// appearing in `inner` where it meets the boundary). Neither
    /// holds for phase-2 *closed* SurroundedTile coronas, where the
    /// outer boundary is cyclic and the central is fully interior.
    /// The function does **not** validate either precondition and
    /// will silently return a bogus patch on violation; callers must
    /// guarantee.
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

    /// `true` if the patch has a boundary (one or more `add_tile`
    /// calls have succeeded); `false` while it's still just a seed.
    pub fn is_growing(&self) -> bool {
        matches!(self.state, PatchState::Growing { .. })
    }

    /// All legal `add_tile` candidates for the current boundary, in
    /// arbitrary order. Each candidate is edge-compatible and passes
    /// the angle-math check, but geometric self-intersection is *not*
    /// pre-filtered — `add_tile` does that final check.
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

    #[allow(dead_code)]
    pub(crate) fn get_matches_for_tile(
        &self,
        tile_id: usize,
    ) -> impl Iterator<Item = PatchMatch> + '_ {
        self.get_all_matches()
            .into_iter()
            .filter(move |m| m.tile_id == tile_id)
    }

    /// All legal `add_tile` candidates whose match touches the boundary
    /// vertex at `vertex_index`. A match "touches" a vertex if the
    /// vertex lies in the closed range `[start_a, start_a + len]` of
    /// the matched edges (see [`cyclic_range_contains`]).
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

    /// All `add_tile` candidates whose matched edge run absorbs at
    /// least one edge in the cyclic closed-closed edge range
    /// `[start_edge, ..., end_edge]`.
    ///
    /// The range is interpreted CCW: if `start_edge <= end_edge` it's
    /// the linear segment `[start_edge, end_edge]`, otherwise it
    /// wraps through `n-1`/`0` and covers `[start_edge, ..., n-1]`
    /// and `[0, ..., end_edge]`. Length is always `(end_edge -
    /// start_edge + n) % n + 1` edges. To query a single edge, pass
    /// `start_edge == end_edge`. To query the full boundary, pass
    /// any `start_edge` with `end_edge = (start_edge + n - 1) % n`.
    ///
    /// A match `(start_a, len)` "absorbs" edges `[start_a, start_a +
    /// len - 1]` (mod `n`). The result is every match whose absorbed
    /// run shares at least one edge with the queried range. Returns
    /// an empty vec for empty boundaries.
    pub fn get_matches_in_edge_range(&self, start_edge: usize, end_edge: usize) -> Vec<PatchMatch> {
        let n = self.boundary_len();
        if n == 0 {
            return Vec::new();
        }
        let range_len = (end_edge + n - start_edge) % n + 1;
        self.get_all_matches()
            .into_iter()
            .filter(|pm| cyclic_arcs_overlap(start_edge, range_len, pm.start_a, pm.len, n))
            .collect()
    }

    #[allow(dead_code)]
    pub(crate) fn ensure_candidates_materialized(&mut self) {
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

    /// Number of distinct tile shapes in the underlying tileset (not
    /// the number of tile instances in the patch).
    pub fn num_tiles(&self) -> usize {
        self.match_index.tileset().num_tiles()
    }

    /// Length of the current boundary in edges. Returns 0 for a `Seed`
    /// state (use [`Self::is_growing`] to check).
    pub fn boundary_len(&self) -> usize {
        match &self.state {
            PatchState::Growing { angles, .. } => angles.len(),
            _ => 0,
        }
    }

    /// The boundary's cyclic angle sequence (one entry per boundary
    /// vertex). Length matches [`Self::boundary_len`]; empty for a
    /// `Seed` state.
    pub fn angles(&self) -> &[i8] {
        match &self.state {
            PatchState::Growing { angles, .. } => angles,
            _ => &[],
        }
    }

    /// Materialise the current patch as a `Rat` describing its
    /// boundary. For a `Seed` state, returns the seed tile.
    pub fn to_rat(&self) -> Rat<T> {
        match &self.state {
            PatchState::Seed { tile_id, .. } => self.match_index.tileset().rat(*tile_id).clone(),
            PatchState::Growing { angles, .. } => Rat::from_slice_unchecked(angles),
        }
    }

    /// Reference to the underlying tileset.
    pub fn tileset(&self) -> &Arc<TileSet<T>> {
        self.match_index.tileset()
    }

    /// Reference to the underlying `MatchTypeIndex` used for candidate
    /// enumeration. Shared via `Arc` so callers can pass it to
    /// associated functions like
    /// [`Self::construct_witness_from_vt_sequence`].
    pub fn match_index(&self) -> &Arc<MatchTypeIndex<T>> {
        &self.match_index
    }

    /// Per-boundary-position [`EdgeInfo`]: which tile shape and which
    /// of its edges occupies each boundary position. Length matches
    /// [`Self::boundary_len`].
    pub fn edges(&self) -> &[EdgeInfo] {
        match &self.state {
            PatchState::Growing { edges, .. } => edges,
            _ => &[],
        }
    }

    /// Per-boundary-position lists of interior tile edges meeting at
    /// that vertex. Non-empty only at junction vertices that have
    /// multiple tiles converging (i.e. an [`OpenVertexType`] with
    /// non-empty `inner`).
    pub fn inner_chains(&self) -> &[Vec<EdgeInfo>] {
        match &self.state {
            PatchState::Growing { inner_chains, .. } => inner_chains,
            _ => &[],
        }
    }

    /// Per-boundary-position **patch tile id**: a fresh monotonic id
    /// assigned to each tile *instance* placed in the patch.
    /// Distinct from `EdgeInfo::tile_id` (which identifies the tile
    /// *shape*): two boundary positions can share a tile_id but have
    /// different patch_tile_ids if they belong to different instances
    /// of the same shape.
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

    /// Rotate the boundary to its canonical (lex-min) starting
    /// position and remap `patch_tile_ids` to dense `0..k` order so
    /// two patches with the same shape and same tile arrangement
    /// compare equal. Useful before hashing or storing as a key.
    ///
    /// Returns the rotation offset applied (= `lex_min_rot(angles)`
    /// computed on the pre-normalize state). A position `p` in the
    /// pre-normalize boundary maps to `(p + n - rot) % n` in the
    /// post-normalize boundary. Returns `0` when no rotation is
    /// applied — including: `Seed` state, empty boundary, and an
    /// already-normalized patch (where `lex_min_rot` returns 0).
    pub fn normalize(&mut self) -> usize {
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
            _ => return 0,
        };

        let n = angles.len();
        if n == 0 {
            return 0;
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
        rot
    }

    #[allow(dead_code)]
    pub(crate) fn candidates_by_start(&self) -> &[Vec<PatchMatch>] {
        match &self.state {
            PatchState::Growing {
                candidates_by_start,
                ..
            } => candidates_by_start.as_deref().unwrap_or(&[]),
            _ => &[],
        }
    }

    #[allow(dead_code)]
    pub(crate) fn vertex_type_at(&self, i: usize) -> Option<PatchVertexInfo> {
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

    /// Coarser variant of [`Self::junction_vertex_type_at`]: returns
    /// `Some(CoarseJunction)` if `i` is a junction, with `cw_edge` /
    /// `ccw_edge` populated and `angle` equal to the boundary angle
    /// at `i`. Drops the interior-tile list — see [`CoarseJunction`]
    /// for the equivalence semantics.
    pub fn coarse_junction_at(&self, i: usize) -> Option<CoarseJunction> {
        let n = self.boundary_len();
        if n == 0 || i >= n {
            return None;
        }
        match &self.state {
            PatchState::Growing { edges, .. } => {
                if !self.is_junction(i) {
                    return None;
                }
                Some(CoarseJunction {
                    cw_edge: edges[(i + n - 1) % n],
                    ccw_edge: edges[i],
                    angle: self.angles()[i],
                })
            }
            _ => None,
        }
    }

    /// `true` if boundary position `i` is a junction vertex — i.e.
    /// the boundary angle there differs from the local tile's natural
    /// internal angle, meaning multiple tiles meet at this vertex.
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

    /// At boundary position `pos`, find the nearest junctions in the
    /// CW and CCW directions and return their relevant tile-offsets
    /// (`(cw_offset, ccw_offset)`). Returns `None` for a `Seed`-state
    /// patch or for `pos >= boundary_len`.
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

    #[allow(dead_code)]
    pub(crate) fn junction_vertices(&self) -> Vec<(usize, PatchVertexInfo)> {
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
    #[allow(dead_code)]
    pub(crate) fn tile_segments(&self) -> Vec<TileSegment> {
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

    /// Attempt to glue the tile described by `pm` onto the current
    /// boundary. Returns `true` on success; `false` if the glue is
    /// rejected (geometric collision, ±hturn junction, or invalid
    /// match dimensions). On failure the patch state is unchanged.
    ///
    /// If the patch is in `Seed` state, the first successful call
    /// transitions it to `Growing` and constructs the initial
    /// boundary.
    ///
    /// # Caller contract
    ///
    /// `pm` must be a **canonical** match obtained from one of:
    /// [`Self::get_all_matches`], [`Self::get_matches_touching_vertex`],
    /// or — for hand construction — built via [`Rat::get_match`] /
    /// [`forward_match_length`] and the boundary's own angle sequence.
    /// `add_tile` does **not** validate that `(pm.start_a, pm.len,
    /// pm.start_b)` describes a real match: it only checks ±hturn
    /// junction degeneracy and post-glue self-intersection. A bogus
    /// interval that happens to produce a non-self-intersecting
    /// polyline will be silently accepted with a geometrically
    /// nonsensical boundary. A `debug_assert!` in
    /// [`angles::glue_raw_angles`] catches this in test builds.
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
    // Keystone glue (seg_len_new == 0): no surviving petal edges, so
    // do not push any new tile edge. For seg_len_new >= 1 the first
    // pushed edge has tile_offset = pm.start_b and subsequent edges
    // follow in petal-CCW order.
    if seg_len_new > 0 {
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
fn compute_segments<T: IsRingOrField + Units>(
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

fn compute_junctions<T: IsRingOrField + Units>(
    angles: &[i8],
    edges: &[EdgeInfo],
    tileset: &TileSet<T>,
) -> Vec<usize> {
    (0..edges.len())
        .filter(|&i| is_junction_at(angles, edges, tileset, i))
        .collect()
}

#[allow(clippy::too_many_arguments)]
fn compute_candidates_at_position<T: IsRingOrField + Units>(
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
fn append_match_candidate<T: IsRingOrField + Units>(
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
fn enumerate_junction_candidates_at<T: IsRingOrField + Units>(
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
    // A match of `len` edges starting at edge `start` (cyclic) touches
    // vertices `start, start+1, …, start+len` (mod n). Compute the
    // forward cyclic distance from `start` to `index` and accept if
    // it's ≤ len. This is correct regardless of whether `start + len`
    // wraps around the boundary — the previous formulation had a bug
    // exactly when `start + len == n`, because the `end <= n` branch
    // checked `index <= end == n` (which `index < n` would fail) and
    // missed the wrap vertex at `end % n == 0`.
    let cyclic_diff = (index + n - start % n) % n;
    cyclic_diff <= len
}

/// Whether two **edge-inclusive** cyclic arcs on a boundary of length
/// `n` share at least one edge.
///
/// Arc `A` = edges `[a, a+1, ..., a+l_a - 1]` (mod `n`);
/// arc `B` = edges `[b, b+1, ..., b+l_b - 1]` (mod `n`).
/// Both arcs are interpreted on a cycle of length `n`; either arc may
/// wrap. Empty arcs (`l_a == 0` or `l_b == 0`) never overlap.
///
/// Returns true iff the two arcs share at least one boundary edge.
/// Internal helper for [`GrowingPatch::get_matches_in_edge_range`].
pub(crate) fn cyclic_arcs_overlap(a: usize, l_a: usize, b: usize, l_b: usize, n: usize) -> bool {
    if l_a == 0 || l_b == 0 || n == 0 {
        return false;
    }
    // Two cyclic arcs overlap iff either's start is "in" the other.
    // "Start of B in A": forward cyclic distance from `a` to `b` is
    // strictly less than `l_a`. Symmetric for "start of A in B".
    let b_in_a = (b + n - a % n) % n < l_a;
    let a_in_b = (a + n - b % n) % n < l_b;
    b_in_a || a_in_b
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
///
/// # Caller contract
///
/// Inherits the precondition from [`angles::glue_raw_angles`]: `pm` must
/// describe a **real** match — `(pm.start_a, pm.len, pm.start_b)` must
/// satisfy the revcomp relation on the `pm.len - 1` interior angles.
/// Obtain `pm` from `GrowingPatch::get_all_matches` /
/// `GrowingPatch::get_matches_touching_vertex` /
/// `Rat::get_match` / `forward_match_length`; hand-constructed
/// `PatchMatch`es are not validated here and may produce geometrically
/// nonsensical glues that downstream self-intersection checks happen
/// to accept.
fn compute_glue_angles<T: IsRingOrField + Units>(
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
fn trace_polyline_from<T: IsRingOrField + Units>(
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

fn trace_boundary_positions<T: IsRingOrField + Units>(angles: &[i8]) -> Vec<T> {
    trace_polyline_from(T::zero(), 0, angles)
}

fn build_boundary_grid<T: IsRingOrField + Units>(
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
        for cell in UnitSquareGrid::seg_neighborhood_of(p1, p2) {
            grid.add(cell, i);
        }
    }
    (grid, seg_data, boundary_edge_ids, n)
}

fn register_segment<T: IsRingOrField + Units>(
    grid: &mut UnitSquareGrid,
    p1: T,
    p2: T,
    id: usize,
) {
    for cell in UnitSquareGrid::seg_neighborhood_of(p1, p2) {
        grid.add(cell, id);
    }
}

fn unregister_segment<T: IsRingOrField + Units>(
    grid: &mut UnitSquareGrid,
    seg_data: &[(T, T)],
    id: usize,
) {
    let (p1, p2) = seg_data[id];
    for cell in UnitSquareGrid::seg_neighborhood_of(p1, p2) {
        grid.remove(cell, id);
    }
}

fn check_segment_clear<T: IsRingOrField + Units>(
    grid: &UnitSquareGrid,
    seg_data: &[(T, T)],
    p1: T,
    p2: T,
    allowed_endpoint: Option<T>,
) -> bool {
    let is_allowed = allowed_endpoint == Some(p2);
    for cell in UnitSquareGrid::seg_neighborhood_of(p1, p2) {
        for &id in grid.get(cell) {
            let (x, y) = seg_data[id];
            if !is_allowed && (p2 == x || p2 == y) {
                return false;
            }
            if intersect_unit_segments(&(p1, p2), &(x, y)) {
                return false;
            }
        }
    }
    true
}

/// Walk consecutive segments of a polyline, checking each against the
/// boundary grid via [`check_segment_clear`]. The last segment's CCW
/// endpoint is permitted to coincide with `allowed_last_endpoint` (this
/// is where the new tile rejoins the existing boundary).
fn segments_all_clear<T: IsRingOrField + Units>(
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
fn dir_of_edge<T: IsRingOrField + Units>(from: T, to: T) -> i8 {
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

    fn ei(tile_id: usize, tile_offset: usize) -> EdgeInfo {
        EdgeInfo {
            tile_id,
            tile_offset,
        }
    }

    #[test]
    fn closed_vertex_type_canonicalises_to_lex_min_rotation() {
        // Four distinct rotations of the same ring should canonicalise
        // identically.
        let petals = [ei(2, 0), ei(0, 1), ei(1, 2), ei(0, 3)];
        let canonical = ClosedVertexType::from_cyclic(&petals);
        for shift in 0..petals.len() {
            let rotated: Vec<EdgeInfo> = (0..petals.len())
                .map(|i| petals[(shift + i) % petals.len()])
                .collect();
            assert_eq!(
                ClosedVertexType::from_cyclic(&rotated),
                canonical,
                "rotation by {shift} should canonicalise to the same VT"
            );
        }

        // The canonical form must start at the lex-min entry.
        let edges = canonical.edges();
        for k in 1..edges.len() {
            assert!(edges[0] <= edges[k]);
        }
    }

    #[test]
    fn closed_vertex_type_distinguishes_non_rotation_orderings() {
        let a = ClosedVertexType::from_cyclic(&[ei(0, 0), ei(0, 1), ei(0, 2)]);
        let b = ClosedVertexType::from_cyclic(&[ei(0, 0), ei(0, 2), ei(0, 1)]);
        assert_ne!(a, b, "different cyclic orderings must be distinct");
    }

    #[test]
    fn closed_vertex_type_from_open_via_closure() {
        let open = OpenVertexType {
            cw: ei(1, 5),
            inner: vec![ei(0, 0), ei(2, 3)],
            ccw: ei(1, 4),
        };
        let closed = ClosedVertexType::from_open_via_closure(&open);
        // Underlying raw ring is [cw, inner..., ccw] = [(1,5),(0,0),(2,3),(1,4)],
        // canonicalised to start at the lex-min entry (0,0).
        let expected = ClosedVertexType::from_cyclic(&[ei(0, 0), ei(2, 3), ei(1, 4), ei(1, 5)]);
        assert_eq!(closed, expected);
        assert_eq!(closed.len(), 4);
    }

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

    /// User-suggested hollow-ring construction: build a curving chain
    /// of hexagons by always gluing the latest hex's edge 1 to the
    /// new hex's edge 5 (= a 60° wedge angle, so the chain curves
    /// inward). The first 4 glues succeed and produce a 5-hex C
    /// around an empty hex-shaped center. The 5th glue (= closing
    /// into a 6-hex hollow ring around the empty center) would
    /// produce a non-simply-connected patch with a hole — which
    /// `GrowingPatch::add_tile` correctly rejects.
    ///
    /// At each step the latest hex's edge 1 sits at boundary position
    /// `boundary_len - 4` (= second of the latest hex's surviving
    /// edges in CCW order after the rotation applied by
    /// `glue_match_to_raw_boundary`).
    #[test]
    fn hollow_hex_ring_closure_rejected() {
        let mut gp = hex_patch();
        // First glue: hex_0's edge 1 → hex_1's edge 5 (start_a=1 on
        // the seed's notional boundary, len=1, start_b=0 so the
        // matched petal edge is start_b-1 = 5 mod 6).
        let first = PatchMatch {
            start_a: 1,
            len: 1,
            start_b: 0,
            tile_id: 0,
        };
        assert!(gp.add_tile(&first), "first glue should succeed");
        // Glues 2-4: continue the chain. start_a is tracked from
        // the post-glue boundary's "second surviving edge of latest
        // hex" = boundary_len - 4.
        for step in 2..=4 {
            let start_a = gp.boundary_len() - 4;
            let pm = PatchMatch {
                start_a,
                len: 1,
                start_b: 0,
                tile_id: 0,
            };
            assert!(
                gp.add_tile(&pm),
                "step {} glue (pm={:?}) should succeed",
                step,
                pm
            );
        }
        assert_eq!(
            gp.boundary_len(),
            22,
            "after 4 glues = 5 hexes in a C, boundary should be 22 edges"
        );

        // Step 5: would add hex_5 closing the chain into a 6-hex
        // ring AROUND AN EMPTY CENTER. The chain has curved enough
        // that hex_5's surviving edges would spatially coincide
        // with hex_0's exposed edges on the other side of the gap
        // (= the chain endpoints face each other across the empty
        // center). `check_segment_clear` rejects: the new tile's
        // segments would touch the existing boundary at non-
        // endpoint positions.
        let start_a = gp.boundary_len() - 4;
        let closing_pm = PatchMatch {
            start_a,
            len: 1,
            start_b: 0,
            tile_id: 0,
        };
        let ok = gp.add_tile(&closing_pm);
        assert!(
            !ok,
            "GrowingPatch::add_tile must refuse the closing glue \
             (= would build a 6-hex ring with a hex-shaped hole at \
             the center, which is non-simply-connected). \
             pm={:?}, current boundary_len={}",
            closing_pm,
            gp.boundary_len()
        );
        // After rejection the patch is unchanged.
        assert_eq!(
            gp.boundary_len(),
            22,
            "rejected glue must leave state unchanged"
        );
    }

    /// Build a 7-hex full corona (1 central + 6 ring tiles) via the
    /// user-suggested approach: glue the central as the FIRST chain
    /// step, then continue the same curving "edge 1 → edge 5"
    /// pattern. Returns the patch.
    ///
    /// Step 1 glues central to hex_0's edge 0 (`start_a = 0`,
    /// `len = 1`). After this glue the rotation puts hex_0's edge 1
    /// at boundary position 0, so step 2 also uses `start_a = 0`.
    /// From step 3 onward, the latest hex's edge 1 sits at
    /// `boundary_len - 4` per the rotation convention.
    fn seven_hex_full_corona() -> GrowingPatch<ZZ12> {
        // Build via three phases:
        //   1. Chain (4 glues): 5 hexes curving around an empty center.
        //   2. Fill (1 glue): central tile fills the inner concavity
        //      via mlen=5 (matching all 5 center-facing edges).
        //   3. Close (1 glue): 6th corona at the remaining wedge via
        //      mlen=3 (matching the 3 wedge-facing edges).
        let mut gp = hex_patch();
        let first = PatchMatch {
            start_a: 1,
            len: 1,
            start_b: 0,
            tile_id: 0,
        };
        assert!(gp.add_tile(&first), "chain glue 1");
        for _ in 2..=4 {
            let start_a = gp.boundary_len() - 4;
            let pm = PatchMatch {
                start_a,
                len: 1,
                start_b: 0,
                tile_id: 0,
            };
            assert!(gp.add_tile(&pm), "chain glue {:?}", pm);
        }
        // Brute-force pick: a match with mlen=5 fills central.
        let central = gp
            .get_all_matches()
            .into_iter()
            .find(|pm| {
                pm.len == 5 && {
                    let mut trial = gp.clone();
                    trial.add_tile(pm)
                }
            })
            .unwrap_or_else(|| panic!("no mlen=5 candidate to fill central"));
        assert!(gp.add_tile(&central), "central fill");
        // Close via mlen=3.
        let closer = gp
            .get_all_matches()
            .into_iter()
            .find(|pm| {
                pm.len == 3 && {
                    let mut trial = gp.clone();
                    trial.add_tile(pm) && trial.boundary_len() == 18
                }
            })
            .unwrap_or_else(|| panic!("no mlen=3 candidate closing the corona"));
        assert!(gp.add_tile(&closer), "ring closure");
        gp
    }

    /// User-flagged invariant (2): the 7-hex full corona's outer
    /// boundary has 18 edges and 6 junctions. Feeding that vt_seq
    /// to `construct_witness_from_vt_sequence` does NOT produce the
    /// 7-hex full corona — it produces a 7-hex CHAIN with boundary
    /// length 30 (= 6 adjacencies, no closure). The vt_seq encodes a
    /// chain of 6 junctions but doesn't enforce the closing
    /// adjacency that would form the ring.
    ///
    /// So neither "hollow ring" (= 6 hexes around empty center) nor
    /// "full corona" (= 7 hexes around central) is produced by
    /// minimal-witness reconstruction. The function returns a chain
    /// — a different geometric shape that also realizes 6 junctions.
    ///
    /// This documents the actual current behavior. It means closed
    /// SurroundedTile entries (which carry the corona's outer
    /// vt_seq, no central) **cannot** be reconstructed back to the
    /// corona via `construct_witness_from_vt_sequence`; doing so
    /// silently produces a chain. We don't currently reconstruct
    /// closed entries, so this is latent.
    /// Precondition test for [`GrowingPatch::construct_witness_from_vt_sequence`]:
    /// when the input vt_seq describes a patch that has FULLY INTERIOR
    /// tiles not captured by any junction's `inner` field, the function
    /// returns a geometrically invalid patch (= self-intersecting
    /// boundary, `Snake::try_from` rejects). This documents the
    /// precondition stated on the function: every realising-patch tile
    /// must appear in at least one junction (as cw / ccw / inner).
    ///
    /// Concrete case: 7-hex full corona has 6 outer junctions (= cw/ccw
    /// from the 2 outer hexes meeting at each corner; `inner` empty).
    /// The central hex is not in any junction's `inner`. Reconstruction
    /// from these 6 vtypes silently produces a self-intersecting chain.
    #[test]
    fn construct_witness_self_intersects_with_fully_interior_tile() {
        let gp = seven_hex_full_corona();
        assert_eq!(gp.boundary_len(), 18);
        let mut vt_seq: Vec<OpenVertexType> = Vec::new();
        for i in 0..gp.boundary_len() {
            if let Some(vt) = gp.junction_vertex_type_at(i) {
                vt_seq.push(vt);
            }
        }
        assert_eq!(vt_seq.len(), 6, "7-hex corona has 6 outer junctions");
        // None of the 6 outer junctions has the central in `inner`.
        for vt in &vt_seq {
            assert!(
                vt.inner.is_empty(),
                "outer-corner junctions have empty inner — central not captured"
            );
        }
        let mi = Arc::clone(gp.match_index());
        let (rebuilt, _) =
            GrowingPatch::construct_witness_from_vt_sequence(&vt_seq, mi).expect("returns Some");
        assert_eq!(
            rebuilt.boundary_len(),
            30,
            "reconstruction places 7 tiles in a spiral (wrong) instead of 6 \
             corona tiles around the central — the central's constraint is \
             missing from the vt_seq."
        );
        assert!(
            Snake::<ZZ12>::try_from(rebuilt.angles()).is_err(),
            "the spiral self-intersects — Snake rejects, but \
             construct_witness_from_vt_sequence does not validate."
        );
    }

    /// Pure unit tests for [`cyclic_range_contains`]. Computed via the
    /// brute reference "which vertices does a `len`-edge match anchored
    /// at `start` touch on a cyclic boundary of length `n`?":
    /// vertices `{start, start+1, …, start+len}` modulo `n` (i.e.
    /// `len + 1` vertices).
    ///
    /// Regression: the previous implementation had a wrap-around bug
    /// exactly when `start + len == n`. In that case the match's
    /// CCW-endpoint vertex is `n mod n = 0`, but the function's
    /// `end <= n` branch checked `index >= start && index <= end`
    /// which is false for `index = 0` whenever `start > 0`.
    /// This test pins all four interesting regimes (interior,
    /// CCW-endpoint-no-wrap, CW-endpoint, wrap-around) plus the
    /// `start + len == n` exact-fit boundary case.
    #[test]
    fn cyclic_range_contains_unit() {
        // (start, len, n) → set of vertex indices the match touches.
        fn brute(start: usize, len: usize, n: usize) -> std::collections::BTreeSet<usize> {
            if len == 0 || n == 0 {
                return std::collections::BTreeSet::new();
            }
            (0..=len).map(|i| (start + i) % n).collect()
        }

        // Pin the regression case directly:
        // start=25, len=1, n=26 should touch vertices {25, 0}.
        // Signature: cyclic_range_contains(start, len, index, n).
        assert!(
            cyclic_range_contains(25, 1, 0, 26),
            "regression: end-at-n-mod-n=0 wrap"
        );
        assert!(cyclic_range_contains(25, 1, 25, 26), "CW endpoint");

        // Exhaustive cross-check against brute over a moderate range.
        for n in [1, 2, 5, 13, 26] {
            for start in 0..n {
                for len in 0..=(n + 1) {
                    let want = brute(start, len, n);
                    for index in 0..n {
                        let got = cyclic_range_contains(start, len, index, n);
                        let expected = want.contains(&index);
                        assert_eq!(got, expected, "n={n} start={start} len={len} index={index}");
                    }
                }
            }
        }

        // Edge cases.
        assert!(!cyclic_range_contains(0, 0, 0, 10), "len=0 → false");
        assert!(!cyclic_range_contains(0, 5, 0, 0), "n=0 → false");
    }

    /// Pure unit test for [`cyclic_arcs_overlap`]. Exhaustively
    /// cross-checks against brute-force edge enumeration over small
    /// boundary sizes plus a handful of regression cases (empty arcs,
    /// zero-length boundary, full-cycle arcs, wraparound on both
    /// arcs).
    #[test]
    fn cyclic_arcs_overlap_unit() {
        fn brute_edges(start: usize, len: usize, n: usize) -> std::collections::BTreeSet<usize> {
            if len == 0 || n == 0 {
                return std::collections::BTreeSet::new();
            }
            (0..len).map(|i| (start + i) % n).collect()
        }
        fn brute_overlap(a: usize, l_a: usize, b: usize, l_b: usize, n: usize) -> bool {
            let arc_a = brute_edges(a, l_a, n);
            let arc_b = brute_edges(b, l_b, n);
            !arc_a.is_disjoint(&arc_b)
        }

        // Exhaustive cross-check over moderate sizes.
        for n in [1, 2, 5, 8, 13] {
            for a in 0..n {
                for l_a in 0..=(n + 1) {
                    for b in 0..n {
                        for l_b in 0..=(n + 1) {
                            let got = cyclic_arcs_overlap(a, l_a, b, l_b, n);
                            let want = brute_overlap(a, l_a, b, l_b, n);
                            assert_eq!(
                                got, want,
                                "mismatch: a={a} l_a={l_a} b={b} l_b={l_b} n={n}"
                            );
                        }
                    }
                }
            }
        }

        // Targeted edge cases.
        assert!(
            !cyclic_arcs_overlap(0, 0, 0, 5, 10),
            "empty arc never overlaps"
        );
        assert!(
            !cyclic_arcs_overlap(0, 5, 0, 0, 10),
            "empty arc never overlaps (other side)"
        );
        assert!(!cyclic_arcs_overlap(0, 5, 0, 5, 0), "n=0 → false");
        assert!(
            cyclic_arcs_overlap(0, 10, 5, 1, 10),
            "full-cycle A vs any non-empty B"
        );
        assert!(
            cyclic_arcs_overlap(7, 5, 1, 2, 10),
            "wraparound A vs interior B"
        );
        assert!(!cyclic_arcs_overlap(0, 3, 5, 3, 10), "disjoint interiors");
        assert!(cyclic_arcs_overlap(0, 3, 2, 3, 10), "edge 2 shared");
    }

    /// `get_matches_in_edge_range` agreement with the brute-force
    /// derivation (= filter `get_all_matches()` by edge-set
    /// intersection). Exhaustive across all start/end positions on
    /// real BFS-grown patches in hex and spectre tilesets.
    #[test]
    fn get_matches_in_edge_range_matches_brute_force() {
        for ts in [
            Arc::new(TileSet::new(vec![
                Rat::try_from(&tiles::hexagon::<ZZ12>()).unwrap()
            ])),
            Arc::new(TileSet::new(vec![
                Rat::try_from(&tiles::spectre::<ZZ12>()).unwrap()
            ])),
        ] {
            let mut gp = GrowingPatch::new(Arc::clone(&ts), 0);
            let first = gp.get_all_matches().into_iter().next().expect("seed match");
            assert!(gp.add_tile(&first), "seed add");
            let n = gp.boundary_len();
            assert!(n > 0);
            let all = gp.get_all_matches();
            for start in 0..n {
                for end in 0..n {
                    let range_len = (end + n - start) % n + 1;
                    let want: std::collections::BTreeSet<(usize, usize, usize, usize)> = all
                        .iter()
                        .filter(|pm| cyclic_arcs_overlap(start, range_len, pm.start_a, pm.len, n))
                        .map(|pm| (pm.start_a, pm.len, pm.start_b, pm.tile_id))
                        .collect();
                    let got: std::collections::BTreeSet<(usize, usize, usize, usize)> = gp
                        .get_matches_in_edge_range(start, end)
                        .into_iter()
                        .map(|pm| (pm.start_a, pm.len, pm.start_b, pm.tile_id))
                        .collect();
                    assert_eq!(
                        got, want,
                        "mismatch on n={n} start={start} end={end} range_len={range_len}"
                    );
                }
            }
        }
    }

    /// Identity test: a range covering the full boundary returns
    /// every match (= equivalent to `get_all_matches()`).
    #[test]
    fn get_matches_in_edge_range_full_boundary_equals_all() {
        let ts: Arc<TileSet<ZZ12>> =
            Arc::new(TileSet::new(vec![
                Rat::try_from(&tiles::spectre::<ZZ12>()).unwrap()
            ]));
        let mut gp = GrowingPatch::new(Arc::clone(&ts), 0);
        let first = gp.get_all_matches().into_iter().next().unwrap();
        assert!(gp.add_tile(&first));
        let n = gp.boundary_len();
        let mut all: Vec<_> = gp.get_all_matches();
        all.sort_by_key(|pm| (pm.start_a, pm.len, pm.start_b, pm.tile_id));
        for start in 0..n {
            let end = (start + n - 1) % n;
            let mut got = gp.get_matches_in_edge_range(start, end);
            got.sort_by_key(|pm| (pm.start_a, pm.len, pm.start_b, pm.tile_id));
            assert_eq!(got, all, "full-boundary range from start={start}");
        }
    }

    /// Verify that the *set* of `(OpenVertexType, OpenVertexType)`
    /// junction pairs on a patch boundary is invariant under
    /// `normalize()`. Used by the segbfs closure check, which skips
    /// `normalize()` on trial patches for speed.
    ///
    /// Rationale: `OpenVertexType` is built from `cw`, `inner`, and
    /// `ccw` fields, all `EdgeInfo` (= `tile_id` + `tile_offset`) —
    /// never `patch_tile_id`. And the set of consecutive pairs on a
    /// cyclic boundary is rotation-invariant. So normalize, which
    /// only rotates the boundary and renumbers `patch_tile_id`s,
    /// can't change the pair set.
    #[test]
    fn junction_pair_set_is_normalize_invariant() {
        for ts in [
            Arc::new(TileSet::new(vec![
                Rat::try_from(&tiles::hexagon::<ZZ12>()).unwrap()
            ])),
            Arc::new(TileSet::new(vec![
                Rat::try_from(&tiles::spectre::<ZZ12>()).unwrap()
            ])),
        ] {
            // Grow a patch a few tiles deep and check at each step.
            let mut gp = GrowingPatch::new(Arc::clone(&ts), 0);
            let first = gp.get_all_matches().into_iter().next().unwrap();
            assert!(gp.add_tile(&first));
            for _step in 0..4 {
                let pre_pairs = collect_pair_set(&gp);
                let mut normed = gp.clone();
                normed.normalize();
                let post_pairs = collect_pair_set(&normed);
                assert_eq!(
                    pre_pairs, post_pairs,
                    "junction pair set differs across normalize"
                );
                // Grow one more tile for the next iteration.
                if let Some(pm) = gp.get_all_matches().into_iter().next() {
                    if !gp.add_tile(&pm) {
                        break;
                    }
                } else {
                    break;
                }
            }
        }
    }

    fn collect_pair_set(
        patch: &GrowingPatch<ZZ12>,
    ) -> std::collections::BTreeSet<(OpenVertexType, OpenVertexType)> {
        let n = patch.boundary_len();
        let juncs: Vec<OpenVertexType> = (0..n)
            .filter_map(|i| patch.junction_vertex_type_at(i))
            .collect();
        let k = juncs.len();
        let mut out = std::collections::BTreeSet::new();
        if k < 2 {
            return out;
        }
        for j in 0..k {
            out.insert((juncs[j].clone(), juncs[(j + 1) % k].clone()));
        }
        out
    }

    /// Empty-boundary patch: `get_matches_in_edge_range` returns an
    /// empty vec without panicking.
    #[test]
    fn get_matches_in_edge_range_handles_empty_boundary() {
        // Take an end-state (Closed) patch via from_parts with empty
        // boundary. Direct construction is simpler.
        let ts: Arc<TileSet<ZZ12>> =
            Arc::new(TileSet::new(vec![
                Rat::try_from(&tiles::hexagon::<ZZ12>()).unwrap()
            ]));
        let gp = GrowingPatch::new(Arc::clone(&ts), 0);
        // The seed-state patch returns boundary_len = 0 since seed has
        // no boundary edges materialized yet. Confirm this.
        // (Note: seed state still has `cached_matches`, so
        // get_all_matches returns non-empty even with boundary_len=0.)
        if gp.boundary_len() == 0 {
            let got = gp.get_matches_in_edge_range(0, 0);
            assert!(got.is_empty(), "empty boundary → no matches");
        }
    }

    /// Regression test for the keystone-glue path in
    /// [`build_glued_edges`] and [`update_inner_chains`]: when
    /// `pm.len == m_tile` (= the petal's full perimeter is absorbed
    /// by the match), `seg_len_new == 0` and the petal contributes
    /// zero surviving boundary edges. The pre-fix code:
    ///
    /// - `build_glued_edges` unconditionally pushed one petal edge,
    ///   yielding `seg_len_old + 1` edges instead of `seg_len_old`.
    /// - `update_inner_chains` wrote `chain_cw` into
    ///   `new_inner[seg_len_old]`, which is one past the end of the
    ///   length-`seg_len_old` vector.
    ///
    /// This test calls the private helpers directly with crafted
    /// inputs that exercise the `seg_len_new == 0` path.
    /// Build 8 unit squares as 3x3-minus-top-left-corner. Glue order
    /// is hand-picked so every intermediate patch is simply connected.
    ///
    /// At each step we filter `get_all_matches()` for a match that
    /// adds the tile at a specific geometric position, by checking the
    /// resulting boundary edge count. This is brittle and only works
    /// for unit squares, but suffices for the fixture.
    fn square_grid_3x3_minus_top_left_corner() -> GrowingPatch<ZZ4> {
        // 8 unit squares forming a 3x3 minus the top-left corner
        // (X = present, . = missing):
        //   row 2: . X X
        //   row 1: X X X
        //   row 0: X X X
        //
        // The 7 glues below were extracted by greedy search (pick the
        // first match producing the right boundary length) and then
        // pinned. Each step keeps the cumulative patch simply
        // connected. Pinning the literals avoids the brute search.
        let glues = [
            // boundary 4 → 6: attach a strip-mate to the seed square.
            PatchMatch {
                start_a: 0,
                len: 1,
                start_b: 0,
                tile_id: 0,
            },
            // boundary 6 → 8: extend the row.
            PatchMatch {
                start_a: 0,
                len: 1,
                start_b: 1,
                tile_id: 0,
            },
            // boundary 8 → 10: extend again to make a 1×3 row.
            PatchMatch {
                start_a: 0,
                len: 1,
                start_b: 1,
                tile_id: 0,
            },
            // boundary 10 → 12: turn upward, starting the right column.
            PatchMatch {
                start_a: 0,
                len: 1,
                start_b: 1,
                tile_id: 0,
            },
            // boundary 12 → 12: continue upward.
            PatchMatch {
                start_a: 2,
                len: 2,
                start_b: 1,
                tile_id: 0,
            },
            // boundary 12 → 12: wrap left along the top.
            PatchMatch {
                start_a: 1,
                len: 2,
                start_b: 1,
                tile_id: 0,
            },
            // boundary 12 → 12: drop into the inner tile (1, 1).
            PatchMatch {
                start_a: 1,
                len: 2,
                start_b: 1,
                tile_id: 0,
            },
        ];
        let mut gp = square_patch();
        for (i, pm) in glues.iter().enumerate() {
            assert!(gp.add_tile(pm), "fixture glue {} failed: pm={:?}", i, pm);
        }
        gp
    }

    /// User-suggested scenario: 8 unit squares forming a 3x3 grid
    /// minus the top-left corner — a simply-connected patch with a
    /// concave notch where the missing tile would be. Extract the
    /// boundary's vt_seq (7 junctions: 6 "straight" tile-tile
    /// boundaries + 1 concave-notch corner) and feed it into
    /// `construct_witness_from_vt_sequence`.
    ///
    /// `construct_witness_from_vt_sequence` glues tiles one at a
    /// time around the seq, which on this input does NOT need the
    /// inner tile (1, 1) — the minimal witness for these 7 junctions
    /// is the 7-tile ring of corner+edge tiles around the notch. The
    /// reconstruction therefore succeeds without ever hitting a
    /// seg_len_new == 0 keystone glue. Pin: rebuilt.boundary_len ==
    /// original.boundary_len.
    #[test]
    fn reconstruct_3x3_minus_corner_from_vt_seq() {
        let gp = square_grid_3x3_minus_top_left_corner();
        assert_eq!(
            gp.boundary_len(),
            12,
            "fixture: 3x3-minus-corner has 12 boundary edges"
        );
        let n = gp.boundary_len();
        let mut corona_vt_seq: Vec<OpenVertexType> = Vec::new();
        for i in 0..n {
            if let Some(vt) = gp.junction_vertex_type_at(i) {
                corona_vt_seq.push(vt);
            }
        }
        assert_eq!(
            corona_vt_seq.len(),
            7,
            "fixture: 6 tile-tile straight junctions + 1 concave notch"
        );
        let mi = Arc::clone(gp.match_index());
        let (rebuilt, _junc_positions) =
            GrowingPatch::construct_witness_from_vt_sequence(&corona_vt_seq, mi)
                .expect("3x3-minus-corner vt_seq should reconstruct");
        assert_eq!(
            rebuilt.boundary_len(),
            gp.boundary_len(),
            "reconstructed boundary length should match original",
        );
    }

    #[test]
    fn build_glued_edges_keystone_len() {
        // Old boundary: 8 edges, all tile_id 0 (placeholders).
        let old_edges: Vec<EdgeInfo> = (0..8)
            .map(|i| EdgeInfo {
                tile_id: 0,
                tile_offset: i,
            })
            .collect();
        let old_ptids = vec![0usize; 8];
        // Keystone: petal of m_tile = 4 fully absorbed (pm.len = 4 = m_tile).
        let pm = PatchMatch {
            start_a: 2,
            len: 4,
            start_b: 0,
            tile_id: 1,
        };
        let m_tile = 4;
        let (new_edges, new_ptids) = build_glued_edges(&old_edges, &old_ptids, &pm, m_tile, 99);
        // new_len = seg_len_old + seg_len_new = (8-4) + (4-4) = 4 + 0 = 4.
        assert_eq!(new_edges.len(), 4, "keystone glue: new_len == seg_len_old");
        assert_eq!(new_ptids.len(), 4);
        // None of the new edges should belong to the petal (= no petal
        // edge survives the keystone glue).
        for e in &new_edges {
            assert_ne!(e.tile_id, pm.tile_id, "no surviving petal edges");
        }
    }

    #[test]
    fn update_inner_chains_keystone_no_oob() {
        // Same setup as build_glued_edges_keystone_len.
        let old_edges: Vec<EdgeInfo> = (0..8)
            .map(|i| EdgeInfo {
                tile_id: 0,
                tile_offset: i,
            })
            .collect();
        let old_inner: Vec<Vec<EdgeInfo>> = vec![Vec::new(); 8];
        let old_ptids = vec![0usize; 8];
        let pm = PatchMatch {
            start_a: 2,
            len: 4,
            start_b: 0,
            tile_id: 1,
        };
        let new_n = 4; // seg_len_old + seg_len_new = 4 + 0.
        let new_inner = update_inner_chains(&old_inner, &old_edges, &pm, new_n, &old_ptids);
        // The pre-fix code would have panicked here. Just verify
        // length and that the call succeeded.
        assert_eq!(new_inner.len(), new_n);
    }

    /// Verify [`angles::glue_raw_angles`] handles `mlen == m` (= the
    /// keystone case, `y_raw_len == 1`) by:
    /// - Returning a result of the correct length (`seg_len_old`).
    /// - Setting a non-`None` `a_yx` / `a_xy` (= the merged junction
    ///   angle written into `result[0]`).
    /// - Producing an angle that is a valid normalized turn (`|merged|
    ///   < hturn`).
    ///
    /// This is purely algebraic; the function makes no claim that the
    /// resulting boundary actually corresponds to a realizable simply
    /// connected patch.
    #[test]
    fn glue_raw_angles_keystone_returns_adjusted_result() {
        use crate::intgeom::angles::glue_raw_angles;
        // Self: 8 angles, with 4 consecutive angles forming a revcomp
        // pattern with a hypothetical 4-edge petal whose angles are all 1.
        // Petal angles = [1, 1, 1, 1]; revcomp(petal) reversed and
        // negated = [-1, -1, -1, -1]. So self_angles[3..=6] = [-1, -1, -1, -1].
        // Outside the match, fill arbitrarily; sum will not equal turn but
        // glue_raw_angles is a pure algebraic transform that doesn't care.
        let self_angles = vec![3, 3, 3, -1, -1, -1, -1, 3];
        let other_angles = vec![1, 1, 1, 1];
        // start_a = 3 (= start of match on self), mlen = 4 (= m, keystone),
        // start_b = 0 (= first surviving petal index; for mlen = m, none survive).
        let gr = glue_raw_angles::<ZZ12>(&self_angles, &other_angles, 3, 4, 0)
            .expect("glue should succeed on keystone");
        assert_eq!(
            gr.angles.len(),
            4,
            "keystone result length = seg_len_old = 8 - 4 = 4"
        );
        assert!(
            gr.a_yx.is_some() && gr.a_xy.is_some(),
            "keystone junction angle should be recorded"
        );
        // Pre-fix, `result[0]` was left as the raw old angle
        // (`self_angles[7] = 3`). Post-fix, it's set to the merged
        // junction angle = normalize(x_first + x_last + y - turn)
        // = normalize(self[7] + self[3] + other[0] - 12)
        // = normalize(3 + (-1) + 1 - 12) = normalize(-9) = 3
        // (since -9 mod 12 = 3 in [-6, 6]).
        assert_eq!(
            gr.angles[0], 3,
            "merged junction angle at result[0]: normalize(3 + (-1) + 1 - 12) = 3"
        );
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

    fn verify_edges_consistency<T: IsRingOrField + Units>(
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
    fn assert_minimal_witness_roundtrips_for<T: IsRingOrField + Units>(
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
    fn assert_witness_matches_brute_force<T: IsRingOrField + Units>(
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
    fn assert_junction_angle_sequence_valid<T: IsRingOrField + Units>(
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

    fn brute_force_recurse<T: IsRingOrField + Units>(
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

    fn brute_force_patches<T: IsRingOrField + Units>(
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
    fn classify_candidates<T: IsRingOrField + Units>(
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
    fn snapshot_growing<T: IsRingOrField + Units>(
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
            // Brute-side filter must NOT use `cyclic_range_contains`,
            // otherwise this cross-check is circular (a bug in
            // `cyclic_range_contains` would affect both sides
            // identically and pass). We instead use an explicit
            // "vertex `target` is in `{start, start+1, …, start+len}`
            // mod n" check via modular arithmetic — independent of the
            // function under test.
            let touching_brute: std::collections::BTreeSet<(usize, usize, usize, usize)> =
                brute_set
                    .iter()
                    .copied()
                    .filter(|(start_a, len, _, _)| {
                        let cyclic_diff = (target + n - *start_a % n) % n;
                        cyclic_diff <= *len
                    })
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
        // Five hexagons arranged as a cross: one central hex with four
        // petals on opposite-pair edges (= a 2-axis-symmetric shape,
        // 18-edge boundary). Glue sequence pinned to literal
        // PatchMatch values for reproducibility; targeted boundary
        // lengths after each step are 10, 14, 16, 18.
        let glues = [
            PatchMatch {
                start_a: 0,
                len: 1,
                start_b: 0,
                tile_id: 0,
            },
            PatchMatch {
                start_a: 1,
                len: 1,
                start_b: 1,
                tile_id: 0,
            },
            PatchMatch {
                start_a: 2,
                len: 2,
                start_b: 1,
                tile_id: 0,
            },
            PatchMatch {
                start_a: 9,
                len: 2,
                start_b: 1,
                tile_id: 0,
            },
        ];
        let mut gp = hex_patch();
        for (i, pm) in glues.iter().enumerate() {
            assert!(
                gp.add_tile(pm),
                "five_hex_cross glue {} failed: pm={:?}",
                i,
                pm
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
        // T-tetromino built incrementally as a chain of 4 unit
        // squares: pin each glue's PatchMatch directly. Boundary
        // length grows 4 → 6 → 8 → 10.
        let glues = [
            PatchMatch {
                start_a: 0,
                len: 1,
                start_b: 0,
                tile_id: 0,
            },
            PatchMatch {
                start_a: 0,
                len: 1,
                start_b: 1,
                tile_id: 0,
            },
            PatchMatch {
                start_a: 0,
                len: 1,
                start_b: 1,
                tile_id: 0,
            },
        ];
        let mut gp = square_patch();
        for (i, pm) in glues.iter().enumerate() {
            assert!(
                gp.add_tile(pm),
                "t_tetromino glue {} failed: pm={:?}",
                i,
                pm
            );
        }
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
