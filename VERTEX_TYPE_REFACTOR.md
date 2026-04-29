# Vertex Type Refactor Plan

## Goal

Replace patch-based BFS with vertex-type-driven BFS. Vertex types include inner tile chains,
enabling local transition admissibility on minimal witnesses without depending on full patch context.

## Key Concepts

- **VertexType** = `{cw: EdgeInfo, inner: Vec<EdgeInfo>, ccw: EdgeInfo}` — fully describes the local
  structure at a boundary vertex. The gap angle is derived from the tile geometry.
- **Minimal witness** = smallest GrowingPatch that realizes a VertexType. Constructed by gluing
  cw_tile + inner tiles + ccw_tile sequentially.
- **Inner chain** = sequence of tiles filling the gap between cw and ccw, going counterclockwise
  through the interior. Empty for seed vertex types (first glue of two tiles).
- **Closed** = artificial sink type. A transition covering both incident edges closes the vertex.
- **Cursed** = all transition paths eventually lead to dead types. Computed backwards from dead types.

## Transition Model

At vertex `pos` (incident edges at `pos-1` = cw, `pos` = ccw):

| Case              | Condition              | Result                            |
|-------------------|------------------------|-----------------------------------|
| Covers neither    | —                      | Skip (seed already has type)      |
| Covers ccw only   | pos in match, pos-1 not| CcwSide: inner += [ccw], new ccw  |
| Covers cw only    | pos-1 in match, pos not| CwSide: inner += [cw], new cw     |
| Covers both       | both in match          | Closed transition                 |

## Implementation Phases

### Phase 1: New Types
- Define `VertexType { cw, inner, ccw }` with Hash/Eq
- Define `TransitionSide { Cw, Ccw }` and `Transition { src_id, dst_id, side, tile_id, tile_offset }`
- Keep existing `PatchVertexType` for backward compat during migration

### Phase 2: Inner Chains in GrowingPatch
- Add `inner_chains: Vec<Vec<EdgeInfo>>` to Growing state
- Update `add_tile_growing`: compute inner chains at both junctions
  - Junction 0: `old_inner[ccw_pos] ++ [old_edges[(start_a+mlen-1)%n]]`
  - Junction seg_len_old: `old_inner[start_a] ++ [old_edges[start_a]]`
  - Non-junction surviving: reindexed from old
  - New tile-only positions: empty
- Update `init_from_first_add`: initialize inner chains for the 2 junctions
- Update `clone_for_mutation`: clone inner chains
- Add `vertex_type_at(pos) -> Option<VertexType>` method
- **Tests**: verify inner chain tracking for hex, square, mixed patches

### Phase 3: Minimal Witness Construction
- `fn construct_minimal_witness<V>(vtype: &VertexType, match_index: &Arc<MatchTypeIndex<V>>) -> Option<(GrowingPatch<V>, usize)>`
- Seed with cw.tile_id, glue inner tiles, glue ccw tile
- First glue: same-instance check (cw tile's own edge) → don't add to inner chain
- Subsequent glues: always add
- Assert final inner chain == vtype.inner
- **Tests**: roundtrip — extract VertexType from patch → construct witness → extract again → equal

### Phase 4: Seeding from MatchTypeIndex
- Enumerate all valid matches: for each tile pair, each valid (offset_a, offset_b)
- Glue → 2-tile patch → extract 2 junction vertex types (empty inner)
- Dedup by VertexType
- **Tests**: verify seed types for hex, square

### Phase 5: Transition Search
- New method on GrowingPatch: find transitions at vertex pos
  - Scan candidates_by_start at positions (pos-k)..=pos
  - Classify each match: CwSide, CcwSide, Closed
  - For each: compute new VertexType at the junction
- **Tests**: verify transition classification

### Phase 6: BFS Rewrite
- New VertexTypeIndex with HashMap<VertexType, usize> dedup
- Queue of unexplored type IDs
- For each type: get cached minimal witness, find transitions, glue, extract new types
- Cache minimal witnesses
- Record transitions

### Phase 7: Cursed Computation
- After BFS completes, compute cursed backwards:
  - Mark all dead types (no outgoing transitions) as cursed-source
  - BFS backwards: if all successors of a type are cursed → type is cursed
  - Use reverse adjacency from transitions

### Phase 8: vtype_enum Update
- Serialize VertexType (cw, inner count + entries, ccw)
- Serialize Transition (side, tile_id, tile_offset)
- Update validation to reconstruct from new format

## Edge Convention

- cw tile's "next edge" = `(cw.tile_offset + 1) % cw_tile_len`
- ccw tile's "receiving edge" = `(ccw.tile_offset - 1 + ccw_len) % ccw_len`
- These are the edges where consecutive tiles glue at the vertex

## Invariants

- BFS only explores minimal witnesses (invariant maintained by construction)
- Non-minimal patches only exist as intermediates during tile placement
- cw and ccw are always from different tile instances in BFS (ensured by 2-tile seeding)
- Inner chain always grows by exactly 1 per glue at a vertex (no same-instance issue in BFS)
