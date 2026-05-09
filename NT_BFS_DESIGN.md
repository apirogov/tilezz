# New NT BFS Design

This document describes the replacement neighborhood-type (NT) enumeration
algorithm. The current `neighborhood.rs` BFS is a flawed prototype.
The new implementation rewrites `NeighborhoodIndex::new()` in-place,
keeping the public API (`NeighborhoodType`, `NtTransition`, `NtKind`,
`classify_all`, `validate`, `write_collection`, `parse_file`) unchanged.

## Goal

Enumerate open neighborhoods around a central tile by growing a context patch at
the moving CCW frontier. The BFS state is compact: it stores a vertex-type (VT)
sequence for the context and just enough anchor information to attach the
central tile.

An NT is an open local configuration:

- one central tile;
- a context patch glued to a contiguous covered portion of the central tile;
- a remaining unmatched central-tile gap;
- a context VT sequence ordered from CW to CCW from the context point of view.

Closed neighborhoods, where the central tile has no unmatched edge left, are not
stored as ordinary NT entries. They are represented by `NT_CLOSED_ID` in the
transition graph.

## File Layout

All code lives in `src/intgeom/neighborhood.rs`. The rewrite replaces the BFS
core (`NtBfsState`, `NtStateKey`, `NeighborhoodIndex::new`, and BFS-internal
helpers) while preserving the public data types and the non-BFS methods.

### What to keep unchanged

- `NeighborhoodType` struct, `gap_len()`, `validate()`
- `NtKind` enum
- `NtTransition` struct
- `NT_CLOSED_ID` constant
- `NeighborhoodIndex` struct fields and all public methods except `new()`
- `classify_all()`, `validate()` (the `NeighborhoodIndex` method),
  `write_collection()`, `parse_file()`
- Existing tests in `mod tests`

### What to keep (reusable free functions)

- `tile_boundary()` — produces a single-tile `RawBoundary`
- `valid_two_tile_glue()` — validates a two-tile glue via `compute_glue_angles`
- `in_consumed_range()` — tests if a position falls in the consumed segment
- `remap_surviving_position()` — remaps a surviving position through a glue
- `extract_covered_vertex_types()` — reads VTs from the covered portion

### What to rewrite

- `NtBfsState` — simplified, no `nt_id` field
- `NtStateKey` — single dedup key: `(central_tile_id, cw_anchor_on_central, vt_seq)`
- `reconstruct_context()` — direct anchor invariant, no scanning
- `attach_central()` — rewritten with anchor invariant
- `find_gap_frontier()` — simplified
- `NeighborhoodIndex::new()` — clean BFS loop

### What to remove

- `ContextBoundary` struct (inline the two fields)
- `CentralAttachment` struct (inline the fields)

## Direction And Edge Conventions

Boundary arrays use the existing project convention:

- `angles[i]` is the angle at vertex `i`;
- `edges[i]` is the edge leaving vertex `i` in the CCW traversal direction;
- `EdgeInfo { tile_id, tile_offset }` identifies that boundary edge on its
  source tile;
- `vertex_type_raw_from(edges, inner_chains, pos)` reads `cw = edges[pos - 1]`,
  `inner = inner_chains[pos]`, `ccw = edges[pos]`.

The context VT sequence is ordered CW to CCW from the context point of view.
`vt_seq.last()` is the active frontier VT when the sequence is non-empty.

The central tile is attached to the context by a fixed CW anchor. From the
context point of view, the covered central-tile boundary extends CCW away from
that anchor. From the central tile point of view, the remaining gap starts at
the stored central anchor edge.

## BFS State

```rust
#[derive(Clone)]
struct NtBfsState {
    central_tile_id: usize,

    // Fixed central-tile gap edge at the CW anchor. This is the `start_b`
    // value when gluing the central tile as the "other" boundary. Never changes.
    cw_anchor_on_central: usize,

    // Context VTs ordered CW -> CCW. Empty for seed states before the first
    // context frontier VT has been created.
    vt_seq: Vec<VertexType>,

    // Seed fields: used only to reconstruct empty vt_seq states.
    seed_tile_id: usize,
    cw_junc_on_seed: usize,
    seed_match_len: usize,

    // Fixed CW junction position on the reconstructed context boundary.
    // Invariant: after reconstruction, Rat::get_match at this position
    // must satisfy the anchor invariant (see "Context Reconstruction").
    cw_junc_on_ctx: usize,
}
```

### Deduplication key

```rust
(central_tile_id, cw_anchor_on_central, vt_seq)
```

This is a single key type. Empty `vt_seq` states with different seeds that
produce the same `(central_tile_id, cw_anchor_on_central, [])` coalesce into
one state. `cw_junc_on_ctx` is anchoring metadata, not part of the key.

## Seed Generation

For each pair of tile types (A, B), compute all maximal matches between them
using `GrowingPatch::new(tileset, A).get_all_matches()`.

For each valid match `(start_a, len, start_b, tile_id=B)`:

1. Validate: call `valid_two_tile_glue` on the seed angles + match.
2. Deduplicate oriented pairs: insert `(min(A,start_a,B,start_b), max(...))`
   into a `BTreeSet` to avoid processing the same glue twice.
3. For each of the two orientations:
   - **Orientation 1**: central=A, context seed=B
     - `central_tile_id = A`, `seed_tile_id = B`
     - `cw_anchor_on_central = start_a`, `cw_junc_on_seed = start_b`
     - `cw_junc_on_ctx = start_b`, `seed_match_len = len`
     - `vt_seq = []`
   - **Orientation 2**: central=B, context seed=A
     - `central_tile_id = B`, `seed_tile_id = A`
     - `cw_anchor_on_central = start_b`, `cw_junc_on_seed = start_a`
     - `cw_junc_on_ctx = start_a`, `seed_match_len = len`
     - `vt_seq = []`
4. Skip if `len >= central_tile_edge_count` (fully covered = already closed).
5. Build context boundary via `tile_boundary(tileset, context_tile_id)`.
6. Call `attach_central(context, cw_junc_on_ctx, central_tile_id,
   cw_anchor_on_central, tileset, first_step=true)`.
7. If attachment fails, skip.
8. Deduplicate via `(central_tile_id, cw_anchor_on_central, [])`.
9. Create `NeighborhoodType` entry with `id = entries.len() + 1`.
10. Enqueue state.

## Context Reconstruction

For a dequeued state:

### Case 1: `vt_seq` is empty

The context is the full boundary of `seed_tile_id`:

```rust
let context = tile_boundary(tileset, seed_tile_id);
// cw_junc_on_ctx comes directly from the state
```

### Case 2: `vt_seq` is non-empty

```rust
let (patch, _juncs) = GrowingPatch::construct_witness_from_vt_sequence(
    &vt_seq, match_index
)?;
let context = RawBoundary {
    angles: patch.angles().to_vec(),
    edges: patch.edges().to_vec(),
    inner_chains: patch.inner_chains().to_vec(),
    patch_tile_ids: patch.patch_tile_ids().to_vec(),
};
```

Then recover `cw_junc_on_ctx` using the **anchor invariant**:

```rust
let ctx_rat = Rat::from_slice_unchecked(&context.angles);
let central_rat = tileset.rat(central_tile_id);
let (ctx_match_start, covered_len, central_end) = ctx_rat.get_match(
    (cw_junc_on_ctx as i64, cw_anchor_on_central as i64),
    central_rat,
);
```

**Anchor invariant** (must hold, fail fast if not):

```
ctx_match_start.rem_euclid(ctx_n) == cw_junc_on_ctx
central_end.rem_euclid(central_n) == cw_anchor_on_central
covered_len > 0 && covered_len < central_n
```

If any part fails, return `None` — this is a bug, not a recoverable error.
Do NOT scan all positions as a fallback.

## Attaching The Central Tile

Given a context boundary and `cw_junc_on_ctx`:

1. Call `Rat::get_match` at `(cw_junc_on_ctx, cw_anchor_on_central)` to get
   `(ctx_match_start, covered_len, central_end)`.
2. Validate anchor invariant.
3. Build `PatchMatch { start_a: cw_junc_on_ctx, len: covered_len,
   start_b: cw_anchor_on_central, tile_id: central_tile_id }`.
4. Call `glue_match_to_raw_boundary(context, &pm, tileset, first_step, new_tile_id)`.
   - `first_step = true` only for seed states where `vt_seq` is empty.
   - `new_tile_id = 0` (we don't track ptids for NT boundaries).
5. Call `validate_raw_boundary` on the result.
6. Derive:
   - `gap_start = glue.old_survivor_len` (the CCW end of the old-context
     survivor segment, which is where the central gap begins)
   - `gap_len = central_n - covered_len`
7. Return the augmented boundary, covered_len, gap_start, gap_len.

## Finding The Frontier

The frontier is the CCW-most junction on the augmented boundary that borders
the central-tile gap.

Procedure:

1. Start at `gap_start` (the CW end of the gap).
2. Walk CCW along the gap (up to `gap_len` positions).
3. Return the first position that is a junction.

A position is a junction if `boundary.angles[pos] != tile_rat.seq()[edge.tile_offset]`
where `edge = boundary.edges[pos]`.

If no junction is found within the gap, return `(gap_start + gap_len) % n`
(the position just past the gap end). This means the frontier is at the CCW
end of the gap — there are no intermediate junctions.

The **context-side frontier** is:

```rust
let ctx_frontier = (cw_junc_on_ctx + covered_len) % ctx_n;
```

This is the replay point for growing the context without the central tile.

## Frontier Distance

To compare frontiers across different context boundaries, use an anchored
distance metric, not raw indices.

```
frontier_distance(ctx_n, cw_junc_on_ctx, ctx_frontier) =
    (cw_junc_on_ctx + ctx_n - ctx_frontier) % ctx_n
```

This is the CW distance from the fixed CW anchor to the frontier. It increases
when the frontier advances CCW.

## Petal Enumeration

At the augmented frontier, enumerate candidate petals `(tile_id, tile_offset)`
from the tileset.

**Pre-filter by forward match**: for each `(tile_id, tile_offset)`, compute
`forward_match_length(augmented_angles, augmented_frontier, tile_seq, tile_offset)`.
Skip candidates with `forward_match_length == 0`.

For each surviving candidate:

1. Glue the petal to the augmented boundary:
   ```rust
   let aug_glue = glue_tile_to_raw_boundary(
       &augmented, augmented_frontier, target, tileset, false, 0,
   );
   ```
   If glue fails, skip.
2. Validate the result: `validate_raw_boundary`.
3. Check frontier angle monotonicity for same-frontier refinements (see below).

### Frontier angle monotonicity

When a petal is glued at the same frontier (covered_len unchanged), the
frontier angle must be monotonically non-increasing using normalized signed
angles. This prevents revisiting the same frontier state.

For the **first** petal at a new frontier, there is no monotonicity check —
all candidates are valid. For subsequent same-frontier petals, the new frontier
angle must be `<=` the previous frontier angle.

## Closed Detection

After gluing the petal to the augmented boundary:

- If the augmented boundary has **no** edge with `tile_id == central_tile_id`,
  the central tile was fully consumed. Record a **closed transition**:
  ```rust
  NtTransition { src_id, dst_id: NT_CLOSED_ID, tile_id, tile_offset, ... }
  ```
  Do NOT replay on context. Do NOT enqueue a new state.

## Replay On The Context

For non-closed petals, replay the glue on the original context boundary
(without the central tile).

1. Glue at `ctx_frontier`:
   ```rust
   let replay_glue = glue_tile_to_raw_boundary(
       &context, ctx_frontier, target, tileset, first_step, 0,
   );
   ```
   - `first_step = true` only when `vt_seq` is empty (seed context).
2. Validate: `validate_raw_boundary`.
3. Remap `cw_junc_on_ctx`:
   ```rust
   let new_cw_junc = remap_surviving_position(
       cw_junc_on_ctx, ctx_frontier, replay_glue.match_len, ctx_n,
   );
   ```
   If `cw_junc_on_ctx` was consumed (returns `None`), the replay is invalid.
   Skip.
4. The new context boundary is `replay_glue.boundary`.

## Re-Attach Central To New Context

After replay, attach the central tile to the new context:

1. Call `attach_central(new_context, new_cw_junc, central_tile_id,
   cw_anchor_on_central, tileset, false)`.
2. If attachment returns `None` (fully covered), record a **closed transition**.
   Do NOT enqueue.
3. If attachment returns `Some(None)` (covered_len == central_n, which means
   closed), same as above — record closed transition.
4. If covered_len **decreased** from the source state's covered_len, skip.
   This should not happen for valid petals but is a safety check.
5. Derive the new augmented boundary, covered_len, gap_start, gap_len,
    and new frontier.

## Updating The VT Sequence

After re-attaching the central tile to the new context:

1. Compute the new frontier VT:
   ```rust
   let new_frontier_vt = vertex_type_raw_from(
       &new_augmented.edges, &new_augmented.inner_chains, new_augmented_frontier,
   );
   ```
2. Compute the new context-side frontier distance:
   ```
   new_ctx_n = new_context.angles.len()
   new_ctx_frontier = (new_cw_junc + new_covered_len) % new_ctx_n
   new_frontier_dist = frontier_distance(new_ctx_n, new_cw_junc, new_ctx_frontier)
   ```
3. Compare with the source state's frontier distance. **Important**: compute
   the source frontier distance on the source context, not the augmented
   boundary:
   ```
   src_ctx_frontier = (cw_junc_on_ctx + covered_len) % ctx_n
   src_frontier_dist = frontier_distance(ctx_n, cw_junc_on_ctx, src_ctx_frontier)
   ```
4. Update `vt_seq`:
   - If `vt_seq` is empty → **append** `new_frontier_vt`. This is the first
     context frontier VT.
   - If `new_frontier_dist == src_frontier_dist` (same frontier) → **replace**
     `vt_seq.last()` with `new_frontier_vt`.
   - If `new_frontier_dist != src_frontier_dist` (frontier advanced) →
     **append** `new_frontier_vt`.

## Recording NT Entries

Each non-closed BFS state (both seeds and grown states) corresponds to one
`NeighborhoodType` entry. The entry is derived fresh from the state by
reconstructing the context and attaching the central tile:

```rust
fn make_neighborhood_entry(
    id: usize,
    state: &NtBfsState,
    attachment: &AttachmentResult,  // augmented boundary, gap_start, etc.
) -> NeighborhoodType {
    NeighborhoodType {
        id,
        central_tile_id: state.central_tile_id,
        central_anchor_edge: state.cw_anchor_on_central,
        gap_len: u8::try_from(attachment.gap_len).unwrap(),
        context_vertex_types: state.vt_seq.clone(),
        angles: attachment.augmented.angles.clone(),
        edges: attachment.augmented.edges.clone(),
        inner_chains: attachment.augmented.inner_chains.clone(),
        gap_start: attachment.gap_start,
    }
}
```

For seed states, the entry is created during seed generation. For grown states,
the entry is created when the new state is first discovered (not yet in `seen`).

## Transition Graph

Each successful petal produces exactly one transition:

```rust
NtTransition {
    src_id: current_state_id,
    dst_id: existing_id or new_entry_id or NT_CLOSED_ID,
    tile_id: target.tile_id,
    tile_offset: target.tile_offset,
    match_start: augmented_frontier,
    match_len: aug_glue.match_len,
}
```

Classification (`Dead`, `Undead`, `Blessed`, `Free`) reuses the existing
`classify_all` logic.

## BFS Loop Pseudocode

```
fn new(tileset: Arc<TileSet<T>>) -> NeighborhoodIndex<T> {
    let match_index = Arc::new(MatchTypeIndex::new(tileset.clone()));
    let mut entries: Vec<NeighborhoodType> = Vec::new();
    let mut transitions: Vec<NtTransition> = Vec::new();
    let mut seen: FxHashMap<NtStateKey, usize> = FxHashMap::default();
    let mut queue: VecDeque<NtBfsState> = VecDeque::new();

    // --- Phase 1: Seed Generation ---
    for each pair of tiles (A, B) with valid maximal match (start_a, len, start_b):
        for each orientation (central, context):
            if len >= central_tile_n: continue  // already closed
            let state = NtBfsState { ... vt_seq: [], ... };
            let context = tile_boundary(tileset, context_tile_id);
            let attachment = attach_central(context, ..., first_step=true)?;
            let key = (central_tile_id, cw_anchor_on_central, []);
            if seen.contains(&key): continue
            let id = entries.len() + 1;
            entries.push(make_neighborhood_entry(id, &state, &attachment));
            seen.insert(key, id);
            queue.push_back(state);

    // --- Phase 2: BFS Growth ---
    while let Some(state) = queue.pop_front():
        let (context, cw_junc, covered_len) = reconstruct_and_attach(&state)?;
        let augmented = attach_central(context, cw_junc, ...)?;
        let aug_frontier = find_gap_frontier(&augmented, gap_start, gap_len);
        let ctx_frontier = (cw_junc + covered_len) % ctx_n;
        let src_frontier_dist = frontier_distance(ctx_n, cw_junc, ctx_frontier);

        for each (tile_id, tile_offset) with forward_match > 0 at aug_frontier:
            let aug_glue = glue_tile_to_raw_boundary(&augmented, aug_frontier, target, ...)?;
            validate_raw_boundary(&aug_glue.boundary)?

            // Closed detection
            if no central_tile_id edge in aug_glue.boundary:
                transitions.push(closed_transition);
                continue

            // Replay on context
            let replay = glue_tile_to_raw_boundary(&context, ctx_frontier, target, ...)?;
            validate_raw_boundary(&replay.boundary)?
            let new_cw_junc = remap_surviving_position(cw_junc, ctx_frontier, replay.match_len, ctx_n)?;
            let new_context = replay.boundary;

            // Re-attach central
            let new_attachment = attach_central(&new_context, new_cw_junc, ...)?;
            if new_attachment is closed:
                transitions.push(closed_transition);
                continue
            if new_attachment.covered_len < covered_len: continue

            // Compute new frontier VT and distance
            let new_frontier_vt = vertex_type_raw_from(...);
            let new_frontier_dist = frontier_distance(...);

            // Update vt_seq
            let mut new_vt_seq = state.vt_seq.clone();
            if new_vt_seq.is_empty():
                new_vt_seq.push(new_frontier_vt)
            elif new_frontier_dist == src_frontier_dist:
                *new_vt_seq.last_mut() = new_frontier_vt
            else:
                new_vt_seq.push(new_frontier_vt)

            // Dedup and enqueue
            let key = (central_tile_id, cw_anchor_on_central, new_vt_seq);
            let dst_id = seen.get(&key).cloned().unwrap_or_else(|| {
                let id = entries.len() + 1;
                let new_state = NtBfsState { vt_seq: new_vt_seq, cw_junc_on_ctx: new_cw_junc, ... };
                entries.push(make_neighborhood_entry(id, &new_state, &new_attachment));
                seen.insert(key, id);
                queue.push_back(new_state);
                id
            });
            transitions.push(NtTransition { src_id: state.nt_id, dst_id, ... });

    NeighborhoodIndex { tileset, entries, transitions }
}
```

## Required Tests

Existing tests must continue to pass:

- `square_has_valid_open_types` — square NT collection validates
- `square_roundtrip_collection` — write/parse roundtrip
- `classify_vector_matches_entry_count` — classify_all length matches entries
- `spectre_collection_validates` — spectre validates in release mode

New invariant tests to add:

- seed generation creates both orientations for each valid two-tile glue
- central glue anchor invariant holds for every dequeued state
- `cw_junc_on_ctx` survives and remaps correctly after every context replay
- same-frontier petal growth replaces the last VT
- frontier-advance petal growth appends a new VT
- empty-seed first petal appends the first VT
- closed transitions are recorded when central coverage reaches full tile

## Build Requirements

```bash
cargo fmt
cargo clippy -- -D warnings
cargo test --release
```

Fix all issues before committing.
