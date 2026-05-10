# NT BFS Implementation

## Goal

Enumerate all open neighborhood types (NTs) around a central tile. An NT is a local
configuration: a central tile plus a context patch glued to a contiguous covered
portion of the central tile, leaving an uncovered gap. The BFS grows the context
patch outward along the moving CCW frontier, producing a transition graph that
records how each NT evolves when a new tile (petal) is attached.

## Data Types

```rust
struct NeighborhoodType {
    central_tile_id: usize,
    cw_anchor_on_central: usize,
    cw_anchor_on_context: usize,
    vt_seq: Vec<VertexType>,
}

struct NtTransition {
    src_id: usize,
    dst_id: usize,
    tile_id: usize,
    tile_offset: usize,
}

enum NtKind { Dead, Undead, Blessed, Free }

const NT_CLOSED_ID: usize = 0;

struct NeighborhoodIndex<T> {
    tileset: Arc<TileSet<T>>,
    entries: Vec<NeighborhoodType>,
    transitions: Vec<NtTransition>,
}
```

An NT is uniquely identified by all four fields:
`(central_tile_id, cw_anchor_on_central, cw_anchor_on_context, vt_seq)`.

- `cw_anchor_on_central` — fixed edge on the central tile where the context attaches
- `cw_anchor_on_context` — position on the context boundary matching that anchor
- `vt_seq` — vertex type sequence describing the context patch frontier

All derived data (augmented boundary, gap_start, gap_len, frontier) is computed on
demand from these fields plus the tileset.

## Algorithm Pseudocode

```
fn new(tileset) -> NeighborhoodIndex:
    match_index = MatchTypeIndex::new(tileset)
    entries = []
    transitions = []
    seen = Map<NeighborhoodType, usize>
    queue = Deque<NeighborhoodType>

    // --- Phase 1: Seed Generation ---
    for each match_type in match_index (tile_a, start_a, tile_b, start_b, len):
        for each orientation:
            (central=tile_a, anchor=start_a, ctx_anchor=start_b)
            (central=tile_b, anchor=start_b, ctx_anchor=start_a)
        nt = NeighborhoodType { central_tile_id, cw_anchor_on_central,
                                cw_anchor_on_context, vt_seq=[] }
        if seen.contains(nt): skip
        seen.insert(nt, entries.len())
        entries.push(nt)
        queue.push_back(nt)

    // --- Phase 2: BFS Growth ---
    while let Some(state) = queue.pop_front():
        state_id = seen[state]
        (augmented, gap_start, gap_len, frontier) = reconstruct_and_attach(state, tileset)

        for each candidate petal (tile_id, tile_offset) at frontier:
            glued = glue_petal(augmented, frontier, petal)

            // Closed detection
            if central tile fully consumed in glued:
                transitions.push(NtTransition { src_id=state_id,
                    dst_id=NT_CLOSED_ID, tile_id, tile_offset })
                continue

            // Replay on context (without central tile)
            new_context = replay_on_context(state, frontier, petal, tileset)
            new_augmented, new_gap = re_attach_central(new_context, state, tileset)

            if new_augmented fully covers central:
                transitions.push(closed transition)
                continue

            // Update VT sequence
            new_vt = vertex_type_at_frontier(new_augmented, new_frontier)
            new_vt_seq = update_vt_seq(state.vt_seq, new_vt,
                                        old_frontier_dist, new_frontier_dist)

            new_nt = NeighborhoodType { state.central_tile_id,
                                        state.cw_anchor_on_central,
                                        new_cw_anchor_on_context,
                                        new_vt_seq }

            dst_id = seen.get(new_nt), or:
                new_id = entries.len()
                seen.insert(new_nt, new_id)
                entries.push(new_nt)
                queue.push_back(new_nt)
                dst_id = new_id

            transitions.push(NtTransition { src_id=state_id,
                dst_id, tile_id, tile_offset })

    return NeighborhoodIndex { tileset, entries, transitions }
```

## Incremental Implementation Plan

### Stage A: Minimal API shell
1. Rewrite `neighborhood.rs` with approved data types. `new()` panics,
   `validate`/`classify_all`/`write_collection`/`parse_file` are stubs.
   → `cargo check` passes

### Stage B: Seed generation (Phase 1) via MatchTypeIndex
2. In `new()`, create `MatchTypeIndex`, iterate all match types, create both
   orientations of `NeighborhoodType`, dedup via `FxHashMap`.
   → Test: seed count (square=13, hex=31)

### Stage C: BFS loop skeleton (Phase 2 framework)
3. Add queue + dequeue loop with empty body.
   → Test: dequeue count matches seed count

### Stage D: Reconstruction + attach
4. `reconstruct_context(nt, match_index)` — empty vt_seq returns seed tile
   boundary; non-empty uses `construct_witness_from_vt_sequence`.
   → Test: seed reconstruction roundtrip
5. `attach_central(context, nt, tileset)` — returns augmented boundary + gap info.
   → Test: known augmented boundary for square seed

### Stage E: Frontier + petal candidates
6. `find_gap_frontier(augmented, gap_start, gap_len)`.
   → Test: known frontier position
7. Enumerate petal candidates at frontier.
   → Test: candidate count at known frontier

### Stage F: Growth logic
8. Glue petal + closed detection.
   → Test: closed transitions appear
9. Replay on context + remap surviving position.
   → Test: new context is valid
10. Re-attach central to new context.
    → Test: new augmented boundary valid

### Stage G: VT sequence update + dedup
11. VT sequence update (append vs replace logic).
    → Test: correct vt_seq after growth
12. Dedup + enqueue — full BFS loop.
    → Test: BFS terminates, correct entry count for square

### Stage H: Remaining API
13. `classify_all()` — reuse existing classification logic.
    → Test: correct classification for square
14. `validate()` — reconstruct + check well-formedness.
    → Test: all square entries valid
15. `write_collection()` + `parse_file()`.
    → Test: roundtrip
16. Full integration — `cargo test --release` including spectre.
