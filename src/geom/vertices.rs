//! Vertex and edge data types: how an edge identifies its tile of
//! origin, how junction vertices on a patch boundary are classified,
//! and which side of a junction a glue transition consumed.
//!
//! The data types are small (no behaviour beyond construction and
//! trivial accessors); the module also hosts the small set of
//! boundary-level helpers that operate purely on
//! `(angles, edges, tileset)` triples without needing a full
//! `GrowingPatch` -- `is_junction_at`, `compute_junctions`,
//! `compute_segments`. They live here so the partitioning of a
//! boundary into junctions / per-tile segments is in the same module
//! as the vertex / segment types it produces, and so
//! [`crate::geom::patch`] can stay focused on the live-patch state
//! machine.

use crate::cyclotomic::IsRing;
use crate::geom::matches::{EdgeRange, PatchSegment, Segment};
use crate::geom::rat::lex_min_rot;
use crate::geom::tileset::TileSet;

// ============================================================
// EdgeInfo: identifying one edge of one tile.
// ============================================================

/// Identifies a single tile edge by `(tile_id, tile_offset)`.
///
/// `tile_id` indexes into the patch's [`crate::geom::tileset::TileSet`];
/// `tile_offset` is the edge index within that tile's boundary
/// (`0..tile.len()`).
///
/// On a [`crate::geom::patch::GrowingPatch`] boundary, an `EdgeInfo`
/// at position `i` says "the boundary edge at position `i` is edge
/// `tile_offset` of an instance of tile `tile_id`". Different
/// boundary positions can share the same `EdgeInfo` if the same tile
/// shape appears multiple times in the patch (each tile *instance*
/// gets a separate `patch_tile_id`, tracked elsewhere in
/// `GrowingPatch`).
#[derive(
    Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, serde::Serialize, serde::Deserialize,
)]
pub struct EdgeInfo {
    pub tile_id: usize,
    pub tile_offset: usize,
}

// ============================================================
// OpenVertexType / ClosedVertexType: junction-vertex configurations.
// ============================================================

/// An **open** junction vertex: the arrangement of tiles meeting at a
/// boundary vertex that is *not* fully surrounded.
///
/// At a junction (where the boundary angle differs from the tile's
/// internal angle), multiple tile instances converge. `cw` is the
/// edge arriving from the CW direction, `ccw` is the edge departing
/// in the CCW direction, and `inner` lists any enclosed tile edges
/// between them.
///
/// "Open" means the vertex is on the patch boundary -- there is still
/// room for additional tiles to be glued in (i.e. the cumulative
/// angle around the vertex is strictly less than a full turn).
/// Equivalently, the boundary angle at the vertex is **not**
/// `±half-turn` -- a half-turn boundary angle would mean the
/// boundary doubles back on itself (a degenerate pinched vertex),
/// which is excluded by construction.
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
/// moment a closing glue ([`TransitionSide::Both`] -- both incident
/// boundary edges consumed by the same match) seals the focus
/// vertex: the closed petal ring is `[cw, inner[0], ..., inner[m-1],
/// ccw]`, with the new tile contributing implicitly between `ccw`
/// and `cw` cyclically. See [`Self::from_open_via_closure`].
///
/// # Identity caveat
///
/// Two distinct *new tiles* that close the same `OpenVertexType` and
/// produce the same cyclic edge ring compare equal under this type --
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
    /// lex-min cyclic rotation (computed via Booth's O(n) algorithm in
    /// [`lex_min_rot`]).
    pub fn from_cyclic(petals: &[EdgeInfo]) -> Self {
        let n = petals.len();
        let edges = if n <= 1 {
            petals.to_vec()
        } else {
            let offset = lex_min_rot(petals);
            (0..n).map(|i| petals[(offset + i) % n]).collect()
        };
        ClosedVertexType { edges }
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

// ============================================================
// CoarseJunction: equivalence-by-boundary-shape at a junction.
// ============================================================

/// A coarser variant of [`OpenVertexType`] that ignores the interior
/// tile arrangement at the junction and instead records the
/// **boundary angle** at the vertex.
///
/// At a junction vertex, the boundary angle (= what the patch stores
/// in its `angles` vec) is `full_turn - sum(incident interior tile
/// angles)`. It captures the boundary's local "shape" at the
/// junction -- what matters for future tile attachments -- without
/// depending on which specific tiles are buried in the interior.
///
/// Two junctions with the same `(cw_edge, ccw_edge, angle)` admit
/// the same outgoing tile attachments, even if their interior
/// populations differ. This makes `CoarseJunction` the right
/// equivalence for seg enumeration: it conflates configurations
/// whose only difference is invisible from outside the boundary.
#[derive(
    Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, serde::Serialize, serde::Deserialize,
)]
pub struct CoarseJunction {
    pub cw_edge: EdgeInfo,
    pub ccw_edge: EdgeInfo,
    pub angle: i8,
}

// ============================================================
// TransitionSide: which incident edges a glue consumed at a focus vertex.
// ============================================================

/// Which incident edge(s) of a focus junction vertex a transition
/// consumed.
///
/// Used by both VT and NT enumeration. At an open vertex with two
/// incident boundary edges (CW and CCW), a glue can:
///
/// - consume only the CW edge -> [`TransitionSide::Cw`],
/// - consume only the CCW edge -> [`TransitionSide::Ccw`],
/// - consume **both** incident edges, sealing the vertex ->
///   [`TransitionSide::Both`].
///
/// When the side is `Both`, the destination of the transition is
/// always the closed-vertex sentinel
/// (`analysis::vertextypes::CLOSED_ID`); the transition's
/// `tile_offset` is normalized to the **CW edge's** matching
/// tile-offset by convention.
///
/// For NT transitions (in `neighborhood`), only `Cw` and `Ccw` arise
/// -- NT enumeration steps the central tile along one direction at
/// a time and never produces a `Both` step.
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

// ============================================================
// Boundary-level helpers: junction detection + segment partition.
// ============================================================
//
// These work on raw `(angles, edges, tileset)` triples and don't
// need a live `GrowingPatch` -- they're shared by the patch state
// machine and by the witness-construction / raw-boundary code paths.
// Keeping them next to `EdgeInfo`, `OpenVertexType`, and `PatchSegment`
// puts the partitioning of a boundary into junctions / per-tile
// segments in the same module as the types it produces.

/// `true` if position `pos` on the boundary is a junction vertex --
/// i.e. the boundary angle there differs from the natural internal
/// angle of the tile occupying the boundary edge at that position.
/// Equivalently: at a non-junction position, the boundary curves
/// exactly the way the tile's own outgoing edge prescribes.
pub(crate) fn is_junction_at<T: IsRing>(
    angles: &[i8],
    edges: &[EdgeInfo],
    tileset: &TileSet<T>,
    pos: usize,
) -> bool {
    let ei = edges[pos];
    tileset.rat(ei.tile_id).seq()[ei.tile_offset] != angles[pos]
}

/// Indices of every junction vertex on the boundary, in increasing
/// order. See [`is_junction_at`] for the per-position predicate.
pub(crate) fn compute_junctions<T: IsRing>(
    angles: &[i8],
    edges: &[EdgeInfo],
    tileset: &TileSet<T>,
) -> Vec<usize> {
    (0..edges.len())
        .filter(|&i| is_junction_at(angles, edges, tileset, i))
        .collect()
}

/// Partition the boundary into [`PatchSegment`]s -- maximal contiguous
/// runs of edges from the same tile instance (no junction between them).
///
/// The segment-break condition is the canonical junction check
/// ([`is_junction_at`]): `angles[i] != tile.seq()[edges[i].tile_offset]`.
/// Every glue updates the boundary angle at the new junction (the
/// compute_glue_angles path produces an angle that for non-degenerate
/// tiles is provably different from either tile's natural internal
/// angle there), so this check reliably detects every tile-instance
/// boundary.
///
/// A simpler offset-arithmetic heuristic (same `tile_id` + contiguous
/// `tile_offset` from the segment start) is **not** used here. It
/// would over-segment on intra-tile wrap-arounds (a single tile
/// instance whose surviving boundary edges have offsets like
/// `2,3,4,5,0`) and could under-segment in pathological same-tile-id
/// glues whose offsets happen to line up. The angle-based check is
/// canonical and avoids both.
///
/// Cyclic-vs-linear caveat: when position 0 is not at a junction, a
/// single cyclic tile-instance run is split into two linear segments
/// at the array seam (`[0, first_junction)` and `[last_junction, n)`).
/// See [`PatchSegment`] for details.
pub(crate) fn compute_segments<T: IsRing>(
    angles: &[i8],
    edges: &[EdgeInfo],
    tileset: &TileSet<T>,
) -> Vec<PatchSegment> {
    let n = edges.len();
    if n == 0 {
        return vec![];
    }
    let mut segments: Vec<PatchSegment> = Vec::new();
    let mut seg_start = 0;
    for i in 1..=n {
        let break_here = i == n || is_junction_at(angles, edges, tileset, i);
        if break_here {
            let len = i - seg_start;
            segments.push(PatchSegment::new(
                EdgeRange::new(seg_start, len),
                Segment::new(
                    edges[seg_start].tile_id,
                    EdgeRange::new(edges[seg_start].tile_offset, len),
                ),
            ));
            seg_start = i;
        }
    }
    segments
}

/// Raw `(cw, inner, ccw)` triple extraction at a boundary position
/// (no junction check). Test-only; the live-patch path goes through
/// [`crate::geom::patch::GrowingPatch::junction_vertex_type_at`] which
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
