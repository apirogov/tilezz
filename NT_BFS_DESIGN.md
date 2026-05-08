# New NT BFS Design

## What is a Neighborhood Type (NT)?

An NT describes a local configuration: a **central tile** surrounded by a **context patch** of tiles glued to it.
The central tile has some edges covered by the context ("covered") and some exposed ("gap").
An NT is uniquely identified by `(central_tile_id, cw_anchor_on_central, vt_seq)`.

## Key Abstractions

### Context Patch
- A set of tiles arranged around the covered portion of the central tile.
- Represented by its **VT sequence**: an ordered list of junction VertexTypes between context tiles.
- The VT sequence does NOT include the two gap-boundary junctions (where context meets central tile).
- Reconstructed using `construct_witness_from_vt_sequence` (for non-empty vt_seq) or from the seed tile directly (for empty vt_seq).

### CW Junction (Fixed)
- The junction on the context boundary where context meets the central tile at the **CW end** of the coverage.
- This position on the context boundary is **fixed** — it never moves during BFS growth.
- The corresponding edge on the central tile (`cw_anchor_on_central`) is also fixed.

### Frontier (CCW end, moves during growth)
- The junction on the context boundary where context meets the central tile at the **CCW end** of the coverage.
- This is the **growth point** — new tiles are added here.
- As tiles are added, the frontier moves CCW on the context (= CW on the central tile), covering more of the central tile.
- The frontier junction angle must **monotonically decrease** with each growth step.

### Direction Convention
- The context extends **CCW** along the boundary from the fixed CW junction.
- From the central tile's POV, coverage extends **CW** from the CW anchor.
- The CW anchor on context is fixed. The CCW frontier moves outward.

## BFS State

```rust
struct NtBfsState {
    central_tile_id: usize,
    cw_anchor_on_central: usize,  // edge of central tile at CW junction (fixed)
    vt_seq: Vec<VertexType>,      // context junction VTs (grows during BFS)
    seed_tile_id: usize,          // initial context tile (for empty vt_seq reconstruction)
    cw_junc_on_seed: usize,       // CW junction position on seed tile boundary
    seed_match_len: usize,        // initial match length at CW junction
}
```

- `vt_seq` only grows (VTs appended). This means BFS exploration is monotonic — same state is unlikely to be reached via multiple paths.
- Dedup key: `(central_tile_id, cw_anchor_on_central, vt_seq)`. Track visited states as safety measure.

## Seeding

For each central tile T and each valid match with tile B (using MatchTypeIndex / compute_seed_matches):
- Match gives `(start_a on T, start_b on B, len)`.
- `cw_anchor_on_central = start_a` (edge of T at the CW end of coverage)
- `cw_junc_on_seed = start_b` (position on B's boundary at the CW junction)
- `seed_match_len = len`
- `vt_seq = []` (empty — no context junctions yet)
- If `len >= m` (all central tile edges covered): skip, this would be a trivially closed NT.

The seed tile B is the initial context. B's boundary IS the context boundary. The CW junction is at position `start_b` on B.

## Growth Step (Core Algorithm)

For each dequeued state:

### 1. Reconstruct context boundary
- **Empty vt_seq**: context = seed tile's full boundary (angles from `tileset.rat(seed_tile_id).seq()`). CW junction at `cw_junc_on_seed`.
- **Non-empty vt_seq**: call `construct_witness_from_vt_sequence(&vt_seq, match_index)`. Returns `(patch, junc_positions)`. CW junction is at or near `junc_positions[0]` (the first VT's position). The exact CW junction position needs to be determined — it's at the boundary between the CW fringe and the covered portion. One approach: the CW junction on context is the vertex where the first VT's `cw` edge is. From the returned patch, find the position where `edges[i] == vt_seq[0].cw`.

### 2. Glue central tile to context at CW junction
- Use `Rat::get_match` to find maximal match starting from:
  - Context side: `(cw_junc_pos_on_ctx, cw_anchor_on_central)` as the anchor point
  - Central tile side: `cw_anchor_on_central`
- The match extends CCW on context / CW on central tile, covering the already-covered portion.
- Apply `glue_raw_angles` + edge/inner_chain updates.
- Result: **augmented boundary** (context + central tile combined).

### 3. Find frontier on augmented boundary
- Walk CCW from CW junction on augmented boundary.
- The first junction vertex encountered is the **frontier** (gap-context boundary at CCW end).
- From edge info at this junction:
  - **CW edge** info → which edge of the central tile is at the frontier
  - **CCW edge** info → offset on the last context tile

### 4. Enumerate petal candidates at frontier
- For each tile T' in tileset, for each junction offset j (0..m'-1):
  - Compute `forward_match_length(augmented_angles, frontier_pos, T'_angles, j)`
  - If match > 0, try the glue (flower_petal_glue logic: compute new angles, edges, inner_chains)
  - Validate: `compute_glue_angles` gives valid result, `Snake::is_closed()`, geometry clear (no intersections)
  - Check: **monotonic angle decrease** — new junction angle at frontier must be strictly smaller (in absolute value) than old junction angle

### 5. Replay working petals on original context
For each petal that passed validation on the augmented boundary:
- Remember its `(tile_id, tile_offset)` — the junction edge of the petal.
- On the **original context boundary** (without central tile), find the frontier position.
- Apply the same petal at the frontier using `forward_match_length` + glue (flower_petal_glue logic).
- The match length on the context may differ from the augmented boundary (shorter, since no central tile edges), but the glue should still work.
- Extract the new junction VT from the resulting context boundary.
- **Append** this VT to `vt_seq`.

### 6. Closed check
- After gluing central tile (step 2), check if `get_match` covers ALL edges of the central tile.
- If yes → record a **closed transition** (NT_COMPLETED).
- This means the central tile is fully inside a patch with no exposed edges.

### 7. Enqueue new state
- New state has the updated `vt_seq` (with appended VT).
- All other fields unchanged (`cw_anchor_on_central` is fixed, seed info is fixed).
- Check dedup key `(central_tile_id, cw_anchor_on_central, vt_seq)` against visited set.
- If new, enqueue.

## Closed NTs
- A closed NT is one where the central tile has NO unmatched edges — it's completely surrounded by the context.
- These are the "completed" neighborhoods. They represent valid local configurations in a tiling.
- Closed NTs are important for the transition graph: they are the terminal states.

## Transition Graph & Classification
After BFS completes:
- Build transition graph: each growth step that produces a new NT is a transition.
- Terminal states (no outgoing transitions) that reach closed → **Blessed**.
- Terminal states that don't reach closed → **Dead**.
- States where all paths lead to Dead → **Undead**.
- States with mixed paths → **Free**.
- This is the same classification as the current implementation (`classify_all`).

## Key Functions Used (Existing)
- `construct_witness_from_vt_sequence(vt_seq, match_index)` — rebuild context from VT sequence
- `Rat::get_match((self_start, other_end), other)` — find maximal match between boundaries
- `forward_match_length(self_angles, self_start, other_angles, other_junction)` — local match at junction
- `glue_raw_angles<T>(self_angles, other_angles, start_a, mlen, start_b)` — combine boundaries
- `update_inner_chains(old_inner, old_edges, pm, new_n)` — update inner chain data after glue
- `flower_petal_glue(boundary, junc_pos, targets, tileset, first_step)` — add petal at junction
- `Snake::is_closed()` — validate boundary is closed loop
- `compute_glue_angles<T>(angles, pm, tileset)` — validate glue and get new angles

## Relationship to Existing Code
- The existing `neighborhood.rs` BFS stores full boundary state (angles, edges, inner_chains, gap_start, gap_len) and has bugs in gap tracking and deferred update mechanism.
- The new BFS uses VT sequences as compact state, reconstructing boundaries on-the-fly.
- The existing `construct_witness_from_vt_sequence` is used as-is — NO modification needed.
- The existing `NeighborhoodType` struct and validation can be reused for output.
- Classification logic (`classify_all`) is reused as-is.

## Implementation Order
1. Define `NtBfsState` and reconstruction logic
2. Implement seeding
3. Implement growth step (glue central tile → find frontier → enumerate petals → replay on context)
4. Implement closed detection
5. Implement dedup (visited set)
6. Wire up BFS loop
7. Implement transition recording
8. Reuse classification
9. Write tests: square, hex, spectre
10. Validate all NTs

## Edge Cases & Notes
- **Empty vt_seq seeds**: context is just the seed tile. CW junction at `cw_junc_on_seed`. The seed tile's full boundary is the context boundary.
- **Multiple VT sequences → same patch**: different VT sequences can produce the same context patch. These are treated as different NTs during BFS. Quotient out later if needed.
- **Monotonic angle decrease**: ensures termination. Each growth step strictly reduces the junction angle at the frontier.
- **The CW offset trick**: the offset from CW junction to frontier, measured going CW (the "opposite" direction from growth), remains stable across central tile gluing. This can be used to track frontier position, but walking CCW from CW junction is simpler.
- **Match at CW junction**: `get_match` extends maximally from the CW junction. This covers the entire already-covered portion of the central tile. After the match, the CCW end of the match on the augmented boundary is the frontier.
