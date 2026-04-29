# New Vertex Types Refactoring Plan

## Goal
Replace `VertexType = Vec<(usize, usize)>` with `PatchVertexType { angle, cw, ccw }` throughout.
Simplify `GrowingPatch` to only maintain `edges: Vec<EdgeInfo>`.
Rewrite BFS vertex type enumeration using GrowingPatch with segment-based `recompute_matches()`.
Benchmark performance vs old approach.

## Key Insight
`PatchVertexType { angle: i8, cw: EdgeInfo, ccw: EdgeInfo }` captures exactly what matters for
vertex classification — the boundary gap angle plus the two adjacent edge identities. The inner
tiles between cw and ccw are irrelevant; we never distinguish vertices by that internal structure.

A witness `GrowingPatch` is stored alongside each discovered PatchVertexType for validation
and further computation. Closed (interior) vertices are not tracked — they're irrelevant for BFS.

## Files to Change
1. `src/intgeom/patch.rs` — core GrowingPatch refactor
2. `src/intgeom/vertextypes.rs` — rewrite BFS with PatchVertexType
3. `src/intgeom/boundary_vtype.rs` — DELETE entirely
4. `src/intgeom/mod.rs` — remove `boundary_vtype` module

## Phase 1: Simplify GrowingPatch (patch.rs)

### Remove
- `pub type VertexType = Vec<(usize, usize)>`
- `vertex_types: Vec<VertexType>` field from `PatchState::Growing`
- `vertex_types()` accessor
- `edge_types()` method (superseded by `edges()`)
- `canonicalize_vtx()` function
- `lex_min_rotation_clone()` function
- `synthetic_closed_vtypes_len3/len4()` test helpers
- All `old_local`/`new_local`/`closed_vertex_types` tracking in init_from_first_add and add_tile_growing
- `seed_vtx` construction in init_from_first_add

### Change
- `AddTileDiff` → `pub struct AddTileDiff;` (unit struct — add_tile returns `Option<AddTileDiff>`)
  - All callers just check `.is_some()` anyway. The old fields were only used by the old BFS.
- `from_parts(tileset, angles, vertex_types)` → `from_parts(tileset, angles, edges)`
- `restore_growing(angles, vertex_types, edges)` → `restore_growing(angles, edges)`
- `init_from_first_add` — simplify to only build `edges: Vec<EdgeInfo>`:
  - CCW junction edge = EdgeInfo { tile_id: pm.tile_id, tile_offset: (pm.start_b + mlen) % m }
  - Seed edges from ccw_pos+1 to cw_pos: EdgeInfo { tile_id: seed_id, tile_offset: (ccw_pos+i)%n }
  - CW junction edge = EdgeInfo { tile_id: seed_id, tile_offset: cw_pos }
  - New tile edges after CW junction: EdgeInfo { tile_id: pm.tile_id, tile_offset: (pm.start_b+mlen+k)%m }
  - New tile edges before CCW junction: already covered by CCW junction edge
- `add_tile_growing` — same simplification:
  - CCW junction: prepend new tile edge
  - Old edges from ccw_pos+1 to cw_pos: keep their EdgeInfo from existing edges array
  - CW junction: keep old EdgeInfo from existing edges[cw_pos]
  - New tile edges: EdgeInfo { tile_id: pm.tile_id, tile_offset: ... }

### New compute_seed_matches
Use MatchTypeIndex directly instead of MatchFinder::crossing:
```
for offset in 0..seed.len():
    for candidate in match_index.candidates_starting_at(seed_tile_id, offset):
        let (ns, len, ne) = seed.get_match((offset as i64, candidate.start_b as i64), tile)
        if len == 0 { continue }
        if !junction_gap_nonnegative(...) { continue }
        if let Ok(glued) = seed.try_glue_precomputed((ns, len, ne), tile, true):
            if Snake::try_from(glued.seq()).is_ok():
                push PatchMatch { start_a: ns as usize, len, start_b: ne as usize, tile_id: candidate.tile_b }
```
Plus single-edge candidates for all (ia, ib) pairs via is_single_edge_candidate.

### New recompute_matches (THE PERFORMANCE-CRITICAL CHANGE)
Replace CMI-based approach with segment-based lookup:
```
let rat = self.to_rat();
let segments = self.tile_segments();
let mut seen: BTreeSet<(usize, usize, usize, usize)> = BTreeSet::new(); // dedup by (start_a, len, start_b, tile_id)

for segment in segments:
    let tile_id = segment.tile_id;
    let tile = self.tileset.rat(tile_id);
    for local_k in 0..segment_length:
        let tile_offset = (segment.offset_start + local_k) % tile.len();
        let patch_pos = segment.patch_start + local_k;
        for candidate in self.match_index.candidates_starting_at(tile_id, tile_offset):
            let (ns, len, ne) = rat.get_match((patch_pos as i64, candidate.start_b as i64),
                                                self.tileset.rat(candidate.tile_b))
            if len == 0 { continue }
            if !junction_gap_nonnegative(...) { continue }
            if let Ok(glued) = rat.try_glue_precomputed((ns, len, ne), ..., true):
                if Snake::try_from(glued.seq()).is_ok():
                    let key = (ns as usize, len, ne as usize, candidate.tile_b)
                    if seen.insert(key):
                        push PatchMatch { start_a: ns as usize, len, start_b: ne as usize, tile_id: candidate.tile_b }

// Single-edge candidates at junction vertices
for (junc_idx, _) in self.junction_vertices():
    for tile_id in 0..self.tileset.num_tiles():
        let tile = self.tileset.rat(tile_id);
        for ib in 0..tile.len():
            if !is_single_edge_candidate(angles, junc_idx, tile.seq(), ib) { continue }
            let (ns, len, ne) = rat.get_match((junc_idx as i64, ib as i64), tile)
            if len != 1 { continue }
            if let Ok(glued) = rat.try_glue_precomputed((ns, len, ne), tile, true):
                if Snake::try_from(glued.seq()).is_ok():
                    push PatchMatch { ... }
```

### Tests to update in patch.rs
- **Remove**: tests that assert on `vertex_types()` contents or `AddTileDiff` fields
  - `first_add_produces_growing` — update to not check vertex_types
  - `segment_cyclic_invariant` — update to use edges
  - `junction_vertex_ids_nonempty_after_each_add` — update to use junction_vertices()
  - `hexagon_all_36_matches_produce_valid_bi_hexes` — keep but remove vertex_types checks
  - `square_all_16_matches_produce_valid_bi_squares` — same
  - `to_rat_matches_direct_glue_for_all_matches` — keep as-is
  - `mixed_hex_square_add_tile` — update to use edges
  - All brute-force tests — update to not reference VertexType
  - `*_systematic_vtypes`, `*_vertex_type_index`, `*_realizing_rat_grouping` — move to vertextypes.rs tests or remove
  - `gap_angle_verification` — update to use vertex_type_at().angle
  - `from_vertex_type_ids_roundtrip` — remove (depends on old VertexTypeIndex)
  - `edge_types_*_consistency` — replace with edges consistency tests
  - `verify_edge_type_consistency` — replace with verify_edges_consistency
  - `synthetic_closed_vtypes_*` — remove
- **Keep/upgrade**:
  - `edges_consistent_with_vertex_types` → `edges_self_consistent` (verify edges have correct tile_ids and offsets)
  - `seed_patch_has_no_boundary`, `seed_patch_has_matches`
  - All brute-force tests that just check add_tile succeeds and produces valid rats
  - `to_rat_matches_direct_glue_for_all_matches`

## Phase 2: Rewrite VertexTypeIndex (vertextypes.rs)

### New types
```rust
pub struct PatchVertexTypeInfo<T: IsComplex> {
    vtype: PatchVertexType,
    kind: VertexTypeKind,
    successors: Vec<usize>,
    predecessors: Vec<usize>,
    is_cursed: bool,
    witness: GrowingPatch<T>,
    witness_pos: usize,
}

pub struct PatchVertexTypeIndex<T: IsComplex> {
    tileset: Arc<TileSet<T>>,
    entries: Vec<PatchVertexTypeInfo<T>>,
    reverse: HashMap<PatchVertexType, usize>,
}
```

### BFS Algorithm
```
visited: BTreeSet<PatchVertexType>
queue: VecDeque<(GrowingPatch<T>, usize)>  // (patch, vertex_position)

// Seeds
for seed_id in 0..tileset.num_tiles():
    let seed = GrowingPatch::new(tileset.clone(), seed_id)
    for pm in seed.get_all_matches():
        let mut gp = seed.clone()
        if gp.add_tile(pm).is_some() && gp.is_growing():
            for pos in 0..gp.boundary_len():
                let pvt = gp.vertex_type_at(pos).unwrap()
                if visited.insert(pvt):
                    queue.push_back((gp.clone(), pos))

// BFS
while let Some((gp, pos)) = queue.pop_front():
    let pvt = gp.vertex_type_at(pos).unwrap()
    if kind_map.get(&pvt) == Some(Dead) { continue }

    let touching: Vec<&PatchMatch> = gp.get_matches_touching_vertex(pos).collect()
    if touching.is_empty():
        kind_map.insert(pvt, Dead)
        continue

    kind_map.entry(pvt).or_insert(Open)

    for pm in touching:
        let mut gp2 = gp.clone()
        if gp2.add_tile(pm).is_some() && gp2.is_growing():
            let n = gp.boundary_len()
            let junction_pos = if pm.start_a == pos { n - pm.len } else { 0 }
            let new_pvt = gp2.vertex_type_at(junction_pos).unwrap()
            transitions.push((pvt, new_pvt))
            for new_pos in 0..gp2.boundary_len():
                let nv = gp2.vertex_type_at(new_pos).unwrap()
                if visited.insert(nv):
                    queue.push_back((gp2.clone(), new_pos))
```

### Simplifications vs old code
- `compute_gap_angle` → just read `pvt.angle` (no tile offset summation)
- `reconstruct_minimal_rat` → `witness.to_rat()` (no tile stitching)
- `from_vertex_type_ids` → needs rethinking (can't reconstruct GrowingPatch from PatchVertexTypes alone)
- No closed type tracking at all
- `validate_closed_types` → remove

### Witness GrowingPatch storage
Each `PatchVertexTypeInfo` stores a `GrowingPatch<T>` that realizes this vertex type.
This is memory-intensive but correct. For large tilesets, could be changed to store just the
provenance chain, but for benchmarking purposes the GrowingPatch approach is fine.

### Tests
- Hex: expect same number of PatchVertexTypes as old BoundaryVertexTypes (42)
- Square: expect same as old (36)
- Mixed: expect same as old (274)
- Classification (open/dead/cursed) should match
- Successor/predecessor graph should be isomorphic to old

## Phase 3: Delete boundary_vtype.rs

Remove file and `pub mod boundary_vtype;` from mod.rs.
Its functionality is now in the new PatchVertexTypeIndex.

## Phase 4: Benchmark

Run new PatchVertexTypeIndex on hex, square, mixed tilesets.
Compare timing with old BoundaryVertexTypeIndex.
Expected improvement: no CMI construction per BFS layer, just table lookups.

## Commit Strategy
- Commit 1: Simplify GrowingPatch (remove vertex_types, simplify AddTileDiff, edges-only)
- Commit 2: New recompute_matches + compute_seed_matches using MatchTypeIndex
- Commit 3: Rewrite vertextypes.rs, delete boundary_vtype.rs, benchmark
