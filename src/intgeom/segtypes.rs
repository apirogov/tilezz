//! Open segment types and the BFS-derived segment type index.
//!
//! A **segment** on a patch open boundary is a maximal run of edges all
//! contributed by a single tile instance, bounded at each end by a
//! [`OpenVertexType`] junction. A **segment type** is the equivalence
//! class of such runs, keyed by the ordered pair of bounding VTs
//! `(vt_cw, vt_ccw)`. Two segments share a type iff their boundary
//! geometry is identical — same tile, same exposed tile offsets, same
//! neighboring tile arrangement at both endpoints.
//!
//! [`OpenSegmentTypeIndex`] is the registry of all segment types that
//! appear on the open boundary of any patch reachable by the BFS in
//! [`crate::intgeom::neighborhood::NeighborhoodIndex`]. It is built
//! once from a `NeighborhoodIndex` by enumerating every consecutive
//! junction pair on every reachable phase-1 ctx boundary, dedup'ing,
//! and assigning stable 0-based ids in canonical sort order.
//!
//! Phase-2 closed coronas are terminal in the BFS and never serve as a
//! ctx, so they do not contribute new segment types.

use std::sync::Arc;

use rustc_hash::{FxHashMap, FxHashSet};

use crate::cyclotomic::{IsComplex, IsRingOrField, Units};
use crate::intgeom::matchtypes::MatchTypeIndex;
use crate::intgeom::neighborhood::{NeighborhoodIndex, NeighborhoodType, NtEntry};
use crate::intgeom::patch::{GrowingPatch, OpenVertexType, PatchMatch};

/// An ordered pair of open vertex types defining one segment type.
/// `cw` is the segment's CW endpoint junction (the one reached by
/// walking CW into the segment from the interior); `ccw` is the CCW
/// endpoint. The exposed tile-offset run between the two endpoints is
/// uniquely determined by the pair.
pub type SegmentTypePair = (OpenVertexType, OpenVertexType);

/// Registry of every open segment type reachable on the open boundary
/// of any phase-1 BFS state in the source `NeighborhoodIndex`.
///
/// Ids are stable: built in canonical sort order over the dedup'd pair
/// set, 0-based.
#[derive(Clone, Debug)]
pub struct OpenSegmentTypeIndex {
    segments: Vec<SegmentTypePair>,
    by_pair: FxHashMap<SegmentTypePair, usize>,
}

impl OpenSegmentTypeIndex {
    /// Build the index from a `NeighborhoodIndex`. For each phase-1
    /// state, walks the **extended-match window** on the ctx (LHS:
    /// `j_cw → vt_seq → j_ccw`) and on the aug (RHS: `j_cw → cw_anchor'
    /// → ccw_anchor' → j_ccw`) and records every consecutive
    /// junction pair from both. Segs that appear only as "deep ctx"
    /// (= far from any state's match position) are intentionally
    /// excluded — they don't participate in any rewrite rule's LHS
    /// or RHS, so they have no role in the segment-level rewrite
    /// system.
    pub fn build_from_bfs<T: IsComplex + IsRingOrField + Units>(
        idx: &NeighborhoodIndex<T>,
    ) -> Self {
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        let mut seen: FxHashMap<SegmentTypePair, ()> = FxHashMap::default();

        for entry in idx.entries() {
            let NtEntry::Phase1(nt) = entry else {
                continue;
            };
            collect_rule_window_pairs(nt, &mi, &mut seen);
        }

        let mut segments: Vec<SegmentTypePair> = seen.into_keys().collect();
        segments.sort();
        let by_pair = segments
            .iter()
            .enumerate()
            .map(|(i, p)| (p.clone(), i))
            .collect();
        OpenSegmentTypeIndex { segments, by_pair }
    }

    /// Look up the segment id for a `(cw_vt, ccw_vt)` pair, if known.
    pub fn vertex_types_to_seg(
        &self,
        cw: &OpenVertexType,
        ccw: &OpenVertexType,
    ) -> Option<usize> {
        self.by_pair.get(&(cw.clone(), ccw.clone())).copied()
    }

    /// Look up the VT pair for a segment id, if in range.
    pub fn seg_to_vertex_types(&self, id: usize) -> Option<&SegmentTypePair> {
        self.segments.get(id)
    }

    pub fn len(&self) -> usize {
        self.segments.len()
    }

    pub fn is_empty(&self) -> bool {
        self.segments.is_empty()
    }

    pub fn iter(&self) -> impl Iterator<Item = (usize, &SegmentTypePair)> {
        self.segments.iter().enumerate()
    }
}

/// For one phase-1 state, reconstruct ctx, compute the LHS-extended
/// match window (`j_cw → vt_seq → j_ccw` on ctx) and the RHS-aug
/// window (`j_cw_aug → cw_anchor_aug → ccw_anchor_aug → j_ccw_aug`
/// on aug), and record every consecutive junction pair from both
/// into `seen`. Silently skips states whose geometry doesn't
/// produce a well-defined window.
fn collect_rule_window_pairs<T: IsComplex + IsRingOrField + Units>(
    nt: &NeighborhoodType,
    mi: &Arc<MatchTypeIndex<T>>,
    seen: &mut FxHashMap<SegmentTypePair, ()>,
) {
    let Some((ctx, junc_positions)) =
        GrowingPatch::construct_witness_from_vt_sequence(&nt.vt_seq, Arc::clone(mi))
    else {
        return;
    };
    let ctx_n = ctx.boundary_len();
    if ctx_n == 0 || junc_positions.is_empty() {
        return;
    }
    let first_junc = junc_positions[0];
    let last_junc = *junc_positions.last().unwrap();
    let cw_frontier_ctx = (first_junc + ctx_n - nt.cw_anchor_on_context) % ctx_n;
    let ccw_frontier_ctx = (last_junc + nt.ccw_anchor_on_context) % ctx_n;
    let Some(j_cw_ctx) = next_junction_cw(&ctx, cw_frontier_ctx) else {
        return;
    };
    let Some(j_ccw_ctx) = next_junction_ccw(&ctx, ccw_frontier_ctx) else {
        return;
    };
    let Some(lhs_chain) = walk_junctions_ccw(&ctx, j_cw_ctx, j_ccw_ctx) else {
        return;
    };
    insert_consecutive_pairs(&lhs_chain, seen);

    // Glue central to get aug, then collect the 3-segment RHS window.
    let central_n = mi.tileset().rat(nt.central_tile_id).seq().len();
    let match_len =
        (nt.cw_anchor_on_central + central_n - nt.ccw_anchor_on_central) % central_n;
    if match_len == 0 || match_len >= central_n {
        return;
    }
    let central_pm = PatchMatch {
        start_a: cw_frontier_ctx,
        len: match_len,
        start_b: nt.cw_anchor_on_central,
        tile_id: nt.central_tile_id,
    };
    let mut aug = ctx.clone();
    if !aug.add_tile(&central_pm) {
        return;
    }
    let aug_n = aug.boundary_len();
    let ctx_ptids: FxHashSet<usize> = ctx.patch_tile_ids().iter().copied().collect();
    let Some(central_ptid) = aug
        .patch_tile_ids()
        .iter()
        .find(|p| !ctx_ptids.contains(p))
        .copied()
    else {
        return;
    };
    let Some(cw_anchor_aug) = find_anchor(&aug, aug_n, central_ptid, AnchorSide::Cw) else {
        return;
    };
    let Some(ccw_anchor_aug) = find_anchor(&aug, aug_n, central_ptid, AnchorSide::Ccw) else {
        return;
    };
    let Some(j_cw_aug) = next_junction_cw(&aug, cw_anchor_aug) else {
        return;
    };
    let Some(j_ccw_aug) = next_junction_ccw(&aug, ccw_anchor_aug) else {
        return;
    };
    let rhs_chain: Vec<OpenVertexType> = [j_cw_aug, cw_anchor_aug, ccw_anchor_aug, j_ccw_aug]
        .iter()
        .filter_map(|&p| aug.junction_vertex_type_at(p))
        .collect();
    if rhs_chain.len() == 4 {
        insert_consecutive_pairs(&rhs_chain, seen);
    }
}

fn insert_consecutive_pairs(
    chain: &[OpenVertexType],
    seen: &mut FxHashMap<SegmentTypePair, ()>,
) {
    for w in chain.windows(2) {
        seen.insert((w[0].clone(), w[1].clone()), ());
    }
}

/// One rewrite rule: a single tile + attachment choice plus the
/// resulting RHS (= aug-side extended-match window). The LHS
/// (interface) is **not** stored — rules are grouped by interface
/// externally in [`SegSubstitutionTable::interface_to_rules`].
///
/// Geometric content of the rule:
/// - `(central_tile_id, cw_anchor_on_central)` — which tile is being
///   added and which of its edges is the CW anchor (= which edge of
///   the central matches up to the boundary on the CW side).
/// - `(cw_anchor_on_context, ccw_anchor_on_context)` — where on the
///   pad segments the CW/CCW frontier vertices fall, in edges from
///   the matched-side junction. Together with the (implicit) LHS
///   they pin the full match interval on the ctx.
/// - `rhs` = `[pad_cw', central_seg, pad_ccw']` — three seg ids on
///   the aug boundary after the central is glued.
///
/// Multiple rules can share an interface (= same LHS seg seq). The
/// set of all rules for one interface enumerates all possible tile
/// additions observed for that boundary configuration.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct SegRewriteRule {
    pub central_tile_id: usize,
    pub cw_anchor_on_central: usize,
    pub cw_anchor_on_context: usize,
    pub ccw_anchor_on_context: usize,
    pub rhs: [usize; 3],
}

/// Per-state rewrite rule + per-interface rule set + preimage map.
///
/// - `state_to_rule[i]` (0-based, state id is `i+1`): the `(interface,
///   rule)` pair for that phase-1 state. `None` for phase-2 states
///   and degenerate phase-1 states.
/// - `interface_to_rules[interface]` = all rewrite rules with that
///   LHS interface — the set of all possible tile additions for
///   that boundary configuration. Note: same `(tile, attachment)`
///   from different source states gets deduplicated; rules within a
///   set are distinct.
/// - `interface_to_state_ids[interface]` = phase-1 state ids whose
///   LHS matches the interface (preimage map used for ranking and
///   kind-distribution analysis). Multiple states share an interface
///   when their ctx geometry is identical up to the extended-match
///   window.
#[derive(Clone, Debug)]
pub struct SegSubstitutionTable {
    pub state_to_rule: Vec<Option<(Vec<usize>, SegRewriteRule)>>,
    pub interface_to_rules: FxHashMap<Vec<usize>, Vec<SegRewriteRule>>,
    pub interface_to_state_ids: FxHashMap<Vec<usize>, Vec<usize>>,
}

impl SegSubstitutionTable {
    /// Build the table from a `NeighborhoodIndex` and a pre-built
    /// segment type index. For each phase-1 state, computes the
    /// interface (LHS seg seq) and the rewrite rule, then groups
    /// rules and state ids by interface.
    pub fn build_from_bfs<T: IsComplex + IsRingOrField + Units>(
        idx: &NeighborhoodIndex<T>,
        seg_idx: &OpenSegmentTypeIndex,
    ) -> Self {
        let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));
        let n = idx.num_types();
        let mut state_to_rule: Vec<Option<(Vec<usize>, SegRewriteRule)>> = vec![None; n];
        let mut interface_to_state_ids: FxHashMap<Vec<usize>, Vec<usize>> = FxHashMap::default();
        let mut interface_to_rules: FxHashMap<Vec<usize>, Vec<SegRewriteRule>> =
            FxHashMap::default();
        let mut interface_rule_seen: FxHashMap<Vec<usize>, FxHashMap<SegRewriteRule, ()>> =
            FxHashMap::default();
        for (i, entry) in idx.entries().iter().enumerate() {
            let NtEntry::Phase1(nt) = entry else {
                continue;
            };
            let Some((interface, rule)) = compute_interface_and_rule(nt, seg_idx, &mi) else {
                continue;
            };
            interface_to_state_ids
                .entry(interface.clone())
                .or_default()
                .push(i + 1);
            let seen = interface_rule_seen.entry(interface.clone()).or_default();
            if !seen.contains_key(&rule) {
                seen.insert(rule.clone(), ());
                interface_to_rules
                    .entry(interface.clone())
                    .or_default()
                    .push(rule.clone());
            }
            state_to_rule[i] = Some((interface, rule));
        }
        SegSubstitutionTable {
            state_to_rule,
            interface_to_rules,
            interface_to_state_ids,
        }
    }
}

/// Seg-level kind, derived purely from the rewrite structure (no
/// reference to per-state BFS kinds).
///
/// - `Dead`: the seg never appears on any LHS — no rule can rewrite
///   a sequence containing it. It's a sink.
/// - `Cursed`: the seg appears on some LHS, but every interface
///   containing it has every rule producing an RHS that contains at
///   least one `Dead` or `Cursed` seg. So every rewrite chain
///   inevitably introduces a sink. Iteratively propagated to a
///   fixpoint.
/// - `Free`: the seg has at least one rewrite path that avoids
///   introducing Dead/Cursed segs.
/// - `Blessed`: TBD — no clean substitution-side definition yet
///   (would need a notion of "good terminal" at the seg level).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SegKind {
    Dead,
    Cursed,
    Blessed,
    Free,
}

/// Classify each segment type by the rewrite structure alone.
/// Returns a vec indexed by seg id, length `seg_idx.len()`.
pub fn classify_segments(
    seg_idx: &OpenSegmentTypeIndex,
    table: &SegSubstitutionTable,
) -> Vec<SegKind> {
    let n = seg_idx.len();
    let mut on_lhs = vec![false; n];
    let mut seg_to_interfaces: Vec<Vec<Vec<usize>>> = vec![Vec::new(); n];
    for interface in table.interface_to_rules.keys() {
        let mut seen_in_this_interface: FxHashSet<usize> = FxHashSet::default();
        for &seg in interface {
            on_lhs[seg] = true;
            if seen_in_this_interface.insert(seg) {
                seg_to_interfaces[seg].push(interface.clone());
            }
        }
    }

    // Dead: not on any LHS. (= no rule can be applied to rewrite it)
    let mut kinds = vec![SegKind::Free; n];
    for i in 0..n {
        if !on_lhs[i] {
            kinds[i] = SegKind::Dead;
        }
    }

    // Iterative propagation: a seg is Cursed if every interface
    // containing it has every rule producing an RHS with a
    // Dead/Cursed seg.
    loop {
        let mut changed = false;
        for i in 0..n {
            if kinds[i] != SegKind::Free {
                continue;
            }
            // Check every interface containing seg i.
            let all_trap = seg_to_interfaces[i].iter().all(|iface| {
                let rules = table.interface_to_rules.get(iface).expect("known interface");
                rules.iter().all(|r| {
                    r.rhs
                        .iter()
                        .any(|&s| matches!(kinds[s], SegKind::Dead | SegKind::Cursed))
                })
            });
            if all_trap && !seg_to_interfaces[i].is_empty() {
                kinds[i] = SegKind::Cursed;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }

    kinds
}

/// Walk the boundary CCW from `start` to `end`, collecting every
/// junction VT in order. Both endpoints must themselves be
/// junctions. If `start == end`, the walk traverses the full cycle
/// and the endpoint VT is pushed twice (once at the start, once at
/// the end). Returns `None` if either endpoint is not a junction or
/// the walk fails to reach `end` within `n` steps.
fn walk_junctions_ccw<T: IsComplex + IsRingOrField + Units>(
    patch: &GrowingPatch<T>,
    start: usize,
    end: usize,
) -> Option<Vec<OpenVertexType>> {
    let n = patch.boundary_len();
    if n == 0 {
        return None;
    }
    let start_vt = patch.junction_vertex_type_at(start)?;
    if !patch.is_junction(end) {
        return None;
    }
    let mut chain = vec![start_vt];
    let mut p = start;
    for _ in 0..n {
        p = (p + 1) % n;
        if let Some(vt) = patch.junction_vertex_type_at(p) {
            chain.push(vt);
            if p == end {
                return Some(chain);
            }
        }
    }
    None
}

/// Walk CW from `start` looking for the next junction (excluding
/// `start` itself). Returns the position, or `None` if no other
/// junction exists.
fn next_junction_cw<T: IsComplex + IsRingOrField + Units>(
    patch: &GrowingPatch<T>,
    start: usize,
) -> Option<usize> {
    let n = patch.boundary_len();
    if n == 0 {
        return None;
    }
    let mut p = (start + n - 1) % n;
    for _ in 0..n {
        if p == start {
            return None;
        }
        if patch.is_junction(p) {
            return Some(p);
        }
        p = (p + n - 1) % n;
    }
    None
}

/// Walk CCW from `start` looking for the next junction (excluding
/// `start` itself). Returns the position, or `None` if no other
/// junction exists.
fn next_junction_ccw<T: IsComplex + IsRingOrField + Units>(
    patch: &GrowingPatch<T>,
    start: usize,
) -> Option<usize> {
    let n = patch.boundary_len();
    if n == 0 {
        return None;
    }
    let mut p = (start + 1) % n;
    for _ in 0..n {
        if p == start {
            return None;
        }
        if patch.is_junction(p) {
            return Some(p);
        }
        p = (p + 1) % n;
    }
    None
}

fn pair_to_seg_ids(
    chain: &[OpenVertexType],
    seg_idx: &OpenSegmentTypeIndex,
) -> Option<Vec<usize>> {
    if chain.len() < 2 {
        return None;
    }
    let mut out = Vec::with_capacity(chain.len() - 1);
    for w in chain.windows(2) {
        out.push(seg_idx.vertex_types_to_seg(&w[0], &w[1])?);
    }
    Some(out)
}

fn compute_interface_and_rule<T: IsComplex + IsRingOrField + Units>(
    nt: &NeighborhoodType,
    seg_idx: &OpenSegmentTypeIndex,
    mi: &Arc<MatchTypeIndex<T>>,
) -> Option<(Vec<usize>, SegRewriteRule)> {
    // Reconstruct ctx.
    let (ctx, junc_positions) =
        GrowingPatch::construct_witness_from_vt_sequence(&nt.vt_seq, Arc::clone(mi))?;
    let ctx_n = ctx.boundary_len();
    if ctx_n == 0 || junc_positions.is_empty() {
        return None;
    }
    let first_junc = junc_positions[0];
    let last_junc = *junc_positions.last().unwrap();
    let cw_frontier_ctx = (first_junc + ctx_n - nt.cw_anchor_on_context) % ctx_n;
    let ccw_frontier_ctx = (last_junc + nt.ccw_anchor_on_context) % ctx_n;

    // -- LHS: walk ctx from j_cw to j_ccw CCW through the match interval. --
    let j_cw_ctx = next_junction_cw(&ctx, cw_frontier_ctx)?;
    let j_ccw_ctx = next_junction_ccw(&ctx, ccw_frontier_ctx)?;
    let lhs_chain = walk_junctions_ccw(&ctx, j_cw_ctx, j_ccw_ctx)?;
    let lhs_segs = pair_to_seg_ids(&lhs_chain, seg_idx)?;

    // -- Glue central → aug. --
    let central_n = mi.tileset().rat(nt.central_tile_id).seq().len();
    let match_len = (nt.cw_anchor_on_central + central_n - nt.ccw_anchor_on_central) % central_n;
    if match_len == 0 || match_len >= central_n {
        return None;
    }
    let central_pm = PatchMatch {
        start_a: cw_frontier_ctx,
        len: match_len,
        start_b: nt.cw_anchor_on_central,
        tile_id: nt.central_tile_id,
    };
    let mut aug = ctx.clone();
    if !aug.add_tile(&central_pm) {
        return None;
    }

    // -- RHS: walk aug from j_cw_aug to j_ccw_aug CCW through central_seg. --
    //
    // Locate the central tile's patch_tile_id on aug. Multiple tiles
    // can share the same `tile_id` (e.g., when the tileset has a
    // single tile shape), so we discriminate by *instance* via
    // `patch_tile_ids`. The central is the new tile, so its ptid is
    // the unique one present in aug but not in ctx.
    let aug_n = aug.boundary_len();
    let ctx_ptids: FxHashSet<usize> = ctx.patch_tile_ids().iter().copied().collect();
    let central_ptid = aug
        .patch_tile_ids()
        .iter()
        .find(|p| !ctx_ptids.contains(p))
        .copied()?;
    let cw_anchor_aug = find_anchor(&aug, aug_n, central_ptid, AnchorSide::Cw)?;
    let ccw_anchor_aug = find_anchor(&aug, aug_n, central_ptid, AnchorSide::Ccw)?;
    let j_cw_aug = next_junction_cw(&aug, cw_anchor_aug)?;
    let j_ccw_aug = next_junction_ccw(&aug, ccw_anchor_aug)?;
    let j_cw_vt = aug.junction_vertex_type_at(j_cw_aug)?;
    let cw_anchor_vt = aug.junction_vertex_type_at(cw_anchor_aug)?;
    let ccw_anchor_vt = aug.junction_vertex_type_at(ccw_anchor_aug)?;
    let j_ccw_vt = aug.junction_vertex_type_at(j_ccw_aug)?;
    let rhs = [
        seg_idx.vertex_types_to_seg(&j_cw_vt, &cw_anchor_vt)?,
        seg_idx.vertex_types_to_seg(&cw_anchor_vt, &ccw_anchor_vt)?,
        seg_idx.vertex_types_to_seg(&ccw_anchor_vt, &j_ccw_vt)?,
    ];

    Some((
        lhs_segs,
        SegRewriteRule {
            central_tile_id: nt.central_tile_id,
            cw_anchor_on_central: nt.cw_anchor_on_central,
            cw_anchor_on_context: nt.cw_anchor_on_context,
            ccw_anchor_on_context: nt.ccw_anchor_on_context,
            rhs,
        },
    ))
}

enum AnchorSide {
    Cw,  // CW anchor: ctx-edge CW, central-edge CCW
    Ccw, // CCW anchor: central-edge CW, ctx-edge CCW
}

fn find_anchor<T: IsComplex + IsRingOrField + Units>(
    aug: &GrowingPatch<T>,
    n: usize,
    central_ptid: usize,
    side: AnchorSide,
) -> Option<usize> {
    let ptids = aug.patch_tile_ids();
    for p in 0..n {
        let cw_ptid = ptids[(p + n - 1) % n];
        let ccw_ptid = ptids[p];
        match side {
            AnchorSide::Cw => {
                if cw_ptid != central_ptid && ccw_ptid == central_ptid {
                    return Some(p);
                }
            }
            AnchorSide::Ccw => {
                if cw_ptid == central_ptid && ccw_ptid != central_ptid {
                    return Some(p);
                }
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic::ZZ12;
    use crate::intgeom::neighborhood::NeighborhoodIndex;
    use crate::intgeom::rat::Rat;
    use crate::intgeom::tiles;
    use crate::intgeom::tileset::TileSet;
    use std::sync::OnceLock;

    fn hex_idx() -> &'static NeighborhoodIndex<ZZ12> {
        static IDX: OnceLock<NeighborhoodIndex<ZZ12>> = OnceLock::new();
        IDX.get_or_init(|| {
            let rat = Rat::try_from(&tiles::hexagon::<ZZ12>()).unwrap();
            NeighborhoodIndex::new(Arc::new(TileSet::new(vec![rat])))
        })
    }

    fn spectre_idx() -> &'static NeighborhoodIndex<ZZ12> {
        static IDX: OnceLock<NeighborhoodIndex<ZZ12>> = OnceLock::new();
        IDX.get_or_init(|| {
            let rat = Rat::try_from(&tiles::spectre::<ZZ12>()).unwrap();
            NeighborhoodIndex::new(Arc::new(TileSet::new(vec![rat])))
        })
    }

    #[test]
    fn seg_index_roundtrips_hex() {
        let seg_idx = OpenSegmentTypeIndex::build_from_bfs(hex_idx());
        assert!(seg_idx.len() > 0, "expected at least one segment type");
        for (id, pair) in seg_idx.iter() {
            let (cw, ccw) = pair;
            let looked_up = seg_idx
                .vertex_types_to_seg(cw, ccw)
                .expect("pair should be in index");
            assert_eq!(looked_up, id, "id roundtrip failed for {:?}", pair);
        }
    }

    #[test]
    fn seg_index_roundtrips_spectre() {
        let seg_idx = OpenSegmentTypeIndex::build_from_bfs(spectre_idx());
        assert!(seg_idx.len() > 0, "expected at least one segment type");
        for (id, pair) in seg_idx.iter() {
            let (cw, ccw) = pair;
            assert_eq!(
                seg_idx.vertex_types_to_seg(cw, ccw),
                Some(id),
                "id roundtrip failed for id {}",
                id
            );
        }
    }

    #[test]
    fn seg_counts_printed() {
        // Data-only sanity log. Doesn't pin counts (those will drift
        // as BFS structure changes); just records orders of magnitude.
        for (name, idx) in [
            ("hex", hex_idx()),
            ("spectre", spectre_idx()),
        ] {
            let seg_idx = OpenSegmentTypeIndex::build_from_bfs(idx);
            eprintln!("{name}: {} distinct segment types", seg_idx.len());
        }
    }

    #[test]
    fn substitution_table_builds_hex() {
        let seg_idx = OpenSegmentTypeIndex::build_from_bfs(hex_idx());
        let table = SegSubstitutionTable::build_from_bfs(hex_idx(), &seg_idx);
        let phase1_count = hex_idx()
            .entries()
            .iter()
            .filter(|e| matches!(e, NtEntry::Phase1(_)))
            .count();
        let rule_count = table.state_to_rule.iter().filter(|r| r.is_some()).count();
        assert_eq!(
            rule_count, phase1_count,
            "expected one rule per phase-1 state"
        );
        let total_unique_rules: usize =
            table.interface_to_rules.values().map(|v| v.len()).sum();
        eprintln!(
            "hex: {} phase-1 states, {} distinct interfaces, {} distinct (interface,rule) pairs",
            phase1_count,
            table.interface_to_state_ids.len(),
            total_unique_rules,
        );
    }

    #[test]
    fn substitution_table_builds_spectre() {
        let seg_idx = OpenSegmentTypeIndex::build_from_bfs(spectre_idx());
        let table = SegSubstitutionTable::build_from_bfs(spectre_idx(), &seg_idx);
        let phase1_count = spectre_idx()
            .entries()
            .iter()
            .filter(|e| matches!(e, NtEntry::Phase1(_)))
            .count();
        let rule_count = table.state_to_rule.iter().filter(|r| r.is_some()).count();
        assert_eq!(
            rule_count, phase1_count,
            "expected one rule per phase-1 state on spectre"
        );
        let total_unique_rules: usize =
            table.interface_to_rules.values().map(|v| v.len()).sum();
        eprintln!(
            "spectre: {} phase-1 states, {} distinct interfaces, {} distinct (interface,rule) pairs",
            phase1_count,
            table.interface_to_state_ids.len(),
            total_unique_rules,
        );
    }

    /// Print the seg-kind distribution per tileset. Diagnostic only
    /// — no hard count assertions. NOTE: the current scheme (Dead =
    /// not on any LHS, Cursed propagated, else Free) leaves all hex
    /// and spectre segs in Dead/Cursed because the rewrite system
    /// only captures phase-1 → phase-1 transitions. A full mapping
    /// requires phase-1 → phase-2 transitions too; that work is on
    /// hold pending a design rethink.
    /// Closure check: does the combined seg set from phase-1 (ctx +
    /// aug) and phase-2 boundaries cover every seg that arises from
    /// one further glue on any of those patches?
    ///
    /// If yes → the BFS catalog is seg-complete and we don't need a
    /// new algorithm. If no → the gap is real and a different
    /// enumeration is needed.
    #[test]
    #[ignore]
    fn neighborhood_bfs_seg_closure_check() {
        for (name, idx) in [("hex", hex_idx()), ("spectre", spectre_idx())] {
            let mi = Arc::new(MatchTypeIndex::new(Arc::clone(idx.tileset())));

            // Step 1: build S = union of all junction pairs visible
            // on phase-1 ctx, phase-1 aug, and phase-2 patch boundaries.
            let mut s: FxHashSet<SegmentTypePair> = FxHashSet::default();
            let mut all_patches: Vec<GrowingPatch<ZZ12>> = Vec::new();

            for entry in idx.entries() {
                match entry {
                    NtEntry::Phase1(nt) => {
                        let Some((ctx, junc_positions)) =
                            GrowingPatch::construct_witness_from_vt_sequence(
                                &nt.vt_seq,
                                Arc::clone(&mi),
                            )
                        else {
                            continue;
                        };
                        record_cyclic_pairs(&ctx, &mut s);
                        let ctx_n = ctx.boundary_len();
                        if ctx_n == 0 || junc_positions.is_empty() {
                            continue;
                        }
                        let first_junc = junc_positions[0];
                        let cw_frontier =
                            (first_junc + ctx_n - nt.cw_anchor_on_context) % ctx_n;
                        let central_n = mi.tileset().rat(nt.central_tile_id).seq().len();
                        let match_len = (nt.cw_anchor_on_central + central_n
                            - nt.ccw_anchor_on_central)
                            % central_n;
                        if match_len == 0 || match_len >= central_n {
                            continue;
                        }
                        let central_pm = PatchMatch {
                            start_a: cw_frontier,
                            len: match_len,
                            start_b: nt.cw_anchor_on_central,
                            tile_id: nt.central_tile_id,
                        };
                        let mut aug = ctx.clone();
                        if !aug.add_tile(&central_pm) {
                            continue;
                        }
                        record_cyclic_pairs(&aug, &mut s);
                        all_patches.push(aug);
                    }
                    NtEntry::Phase2(st) => {
                        let Some(patch) = st.to_patch(Arc::clone(&mi)) else {
                            continue;
                        };
                        record_cyclic_pairs(&patch, &mut s);
                        all_patches.push(patch);
                    }
                }
            }
            eprintln!(
                "{name}: phase-1+phase-2 union seg set |S| = {}, patches collected = {}",
                s.len(),
                all_patches.len()
            );

            // Step 2: for each collected patch, try every legal glue
            // and check if any boundary seg of the result isn't in S.
            let mut new_segs: FxHashSet<SegmentTypePair> = FxHashSet::default();
            let mut glues_tried = 0usize;
            let mut glues_succeeded = 0usize;
            for patch in &all_patches {
                let matches = patch.get_all_matches();
                for pm in matches {
                    glues_tried += 1;
                    let mut trial = patch.clone();
                    if !trial.add_tile(&pm) {
                        continue;
                    }
                    glues_succeeded += 1;
                    let n = trial.boundary_len();
                    let juncs: Vec<OpenVertexType> = (0..n)
                        .filter_map(|i| trial.junction_vertex_type_at(i))
                        .collect();
                    let k = juncs.len();
                    if k < 2 {
                        continue;
                    }
                    for j in 0..k {
                        let pair = (juncs[j].clone(), juncs[(j + 1) % k].clone());
                        if !s.contains(&pair) {
                            new_segs.insert(pair);
                        }
                    }
                }
            }
            eprintln!(
                "{name}: glues tried = {}, succeeded = {}, NEW segs not in S = {}",
                glues_tried,
                glues_succeeded,
                new_segs.len()
            );
        }
    }

    fn record_cyclic_pairs<T: IsComplex + IsRingOrField + Units>(
        patch: &GrowingPatch<T>,
        out: &mut FxHashSet<SegmentTypePair>,
    ) {
        let n = patch.boundary_len();
        let juncs: Vec<OpenVertexType> = (0..n)
            .filter_map(|i| patch.junction_vertex_type_at(i))
            .collect();
        let k = juncs.len();
        if k < 2 {
            return;
        }
        for j in 0..k {
            out.insert((juncs[j].clone(), juncs[(j + 1) % k].clone()));
        }
    }

    #[test]
    fn classify_segments_diagnostic() {
        for (name, idx) in [("hex", hex_idx()), ("spectre", spectre_idx())] {
            let seg_idx = OpenSegmentTypeIndex::build_from_bfs(idx);
            let table = SegSubstitutionTable::build_from_bfs(idx, &seg_idx);
            let seg_kinds = classify_segments(&seg_idx, &table);
            let dead = seg_kinds.iter().filter(|k| **k == SegKind::Dead).count();
            let cursed = seg_kinds.iter().filter(|k| **k == SegKind::Cursed).count();
            let blessed = seg_kinds.iter().filter(|k| **k == SegKind::Blessed).count();
            let free = seg_kinds.iter().filter(|k| **k == SegKind::Free).count();
            eprintln!(
                "{name}: segs={} dead={} cursed={} blessed={} free={}",
                seg_kinds.len(),
                dead,
                cursed,
                blessed,
                free
            );
        }
    }
}
