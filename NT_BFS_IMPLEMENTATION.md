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

- `cw_anchor_on_central` — edge on the central tile where the context attaches
- `cw_anchor_on_context` — position on the reconstructed context boundary
- `vt_seq` — junction vertex types of the context patch boundary (always non-empty)

All derived data (augmented boundary, gap_start, gap_len, frontier) is computed on
demand by reconstructing the context from `vt_seq` via
`GrowingPatch::construct_witness_from_vt_sequence`, then attaching the central tile.

## Key Conventions

- Context patches are always **normalized** (lex-min-rot) before extracting VTs or
  anchors. This ensures positions are consistent between seed generation and
  reconstruction.
- `vt_seq` contains the junction vertex types of the context patch boundary — the
  VTs where context tiles meet each other. These do NOT include central tile edges.
- `construct_witness_from_vt_sequence` always has a non-empty input (empty vt_seq
  never occurs). The empty-vt-seed edge case is eliminated by design.
- Seeds are three-tile configurations: two-tile context + central tile added third.
  No two-tile (empty vt_seq) seeds exist.

## Algorithm Pseudocode

```
fn new(tileset) -> NeighborhoodIndex:
    match_index = MatchTypeIndex::new(tileset)
    entries = []
    transitions = []
    seen = Map<NeighborhoodType, usize>
    queue = Deque<NeighborhoodType>

    // --- Phase 1: Seed Generation ---
    // For each match type, build a normalized two-tile context patch,
    // then try adding every possible third tile as "central".
    for each match_type in match_index (tile_a, start_a, tile_b, start_b, len):
        patch = GrowingPatch::new(tileset, tile_a)
        patch.add_tile(PatchMatch { start_a, len, start_b, tile_id: tile_b })
        patch.normalize()
        vt_seq = extract junction vertex types from patch boundary

        for each third_tile_match in patch.get_all_matches():
            central_tile_id = third_tile_match.tile_id
            cw_anchor_on_central = third_tile_match.start_b
            cw_anchor_on_context = third_tile_match.start_a

            nt = NeighborhoodType { central_tile_id, cw_anchor_on_central,
                                    cw_anchor_on_context, vt_seq }
            if seen.contains(nt): skip
            seen.insert(nt, entries.len())
            entries.push(nt)
            queue.push_back(nt)

    // --- Phase 2: BFS Growth ---
    while let Some(state) = queue.pop_front():
        state_id = seen[state]

        // Reconstruct context from vt_seq, attach central tile
        context = construct_witness_from_vt_sequence(state.vt_seq, match_index)
        augmented = attach_central(context, state, tileset)

        for each candidate petal (tile_id, tile_offset) at frontier:
            glued = glue_petal(augmented, frontier, petal)

            // Closed detection
            if central tile fully consumed in glued:
                transitions.push(NtTransition { src_id=state_id,
                    dst_id=NT_CLOSED_ID, tile_id, tile_offset })
                continue

            // Replay on context (without central tile)
            new_context = replay_on_context(state, frontier, petal, tileset)
            new_augmented = re_attach_central(new_context, state, tileset)

            if new_augmented fully covers central:
                transitions.push(closed transition)
                continue

            // Update VT sequence and context anchor
            new_vt_seq = update_vt_seq(...)
            new_cw_anchor_on_context = remapped anchor on new context

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

## Validation

Every seed is validated by:
1. Reconstructing context from `vt_seq` via `construct_witness_from_vt_sequence`
2. Finding the match at `(cw_anchor_on_context, cw_anchor_on_central)`
3. Adding the central tile via `add_tile` — must succeed

## Incremental Implementation Plan

### Stage A: Minimal API shell ✓
- Approved data types, stub methods.

### Stage B+C: Seed generation + BFS skeleton ✓
- Three-tile seeds via normalized two-tile context + third-tile addition.
- Context VTs extracted from normalized two-tile patch boundary junctions.
- BFS dequeue loop (empty body for now).
- Validation test: all seeds reconstruct and re-attach.
- Counts: square=240, hex=1008.

### Stage D: BFS growth body
- Reconstruction via `construct_witness_from_vt_sequence`.
- Attach central tile, find gap frontier.
- Enumerate petal candidates.

### Stage E: VT sequence update + dedup
- Append vs replace logic for frontier VTs.
- Dedup + enqueue — full BFS loop.
- Test: BFS terminates, correct entry count for square.

### Stage F: Remaining API
- `classify_all()` — reuse existing classification logic.
- `validate()` — reconstruct + check well-formedness.
- `write_collection()` + `parse_file()`.
- Full integration — `cargo test --release` including spectre.
