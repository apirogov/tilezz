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
    central_tile_id: usize,      // tileset tile type of the central tile
    cw_anchor_on_central: usize, // edge on central tile where context attaches (CW side)
    cw_anchor_on_context: usize, // position on reconstructed context boundary (fixed, never changes)
    vt_seq: Vec<VertexType>,     // junction VTs of context boundary (context-context only)
}

struct NtTransition {
    src_id: usize,
    dst_id: usize,          // NT_CLOSED_ID = 0 for closed transitions
    tile_id: usize,         // tileset tile type of the petal
    tile_offset: usize,     // which petal edge faces the frontier junction
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

## Key Concepts

### Augmented boundary

The augmented boundary = context boundary + central tile glued at the anchor. It has
a fixed layout:

```
positions [0, gap_start)       → context surviving edges (seg_len_old)
positions [gap_start, aug_n)   → gap edges (central tile, uncovered)
```

where `seg_len_old = ctx_n - central_match_len`, `gap_len = central_len - central_match_len`,
and `gap_end = (gap_start + gap_len) % aug_n = 0` (always, since they sum to aug_n).

### Frontier

The frontier is the junction at the CCW end of the gap (position 0 on the augmented
boundary). This is always a junction between context and the central tile. The BFS
adds petals at this junction to extend the context coverage.

### vt_seq — context-only junctions

`vt_seq` contains the junction vertex types of the **context** boundary — the VTs
where context tiles meet each other. These do NOT include the central tile. The
frontier junction (which involves the central tile) is NOT in vt_seq.

The CW anchor junction is also NOT in vt_seq (it involves the central tile).
The anchor is tracked by `cw_anchor_on_context` which is a fixed property of the NT.

### Fixed anchor

`cw_anchor_on_central` and `cw_anchor_on_context` never change during BFS growth.
The witness reconstruction is deterministic and monotone (glues petals in vt_seq order),
so the reconstructed context always has the anchor at the same position.

### Patch tile IDs vs tileset tile IDs

- `edges[].tile_id` is the tileset tile type (not unique per instance)
- `patch_tile_ids[]` is the unique per-instance ID assigned by GrowingPatch
- For tracing segments (gap, petal), always use `patch_tile_ids[]`

## Key Conventions

- Context patches are always **normalized** (lex-min-rot) before extracting VTs or
  anchors. This ensures positions are consistent between seed generation and
  reconstruction.
- `construct_witness_from_vt_sequence` always has a non-empty input (empty vt_seq
  never occurs). The empty-vt-seed edge case is eliminated by design.
- Seeds are three-tile configurations: two-tile context + central tile added third.
  No two-tile (empty vt_seq) seeds exist.

## Algorithm: `explore_one(nt, match_index) → Vec<ExploreResult>`

```rust
enum ExploreOutcome { Closed, Open(NeighborhoodType) }
struct ExploreResult { tile_id: usize, tile_offset: usize, outcome: ExploreOutcome }
```

### Phase 1: Build augmented, find frontier, enumerate candidates

```
1. context = construct_witness_from_vt_sequence(nt.vt_seq, match_index)
2. central_match_len = forward_match_length(
       context.angles, nt.cw_anchor_on_context,
       central_seq, nt.cw_anchor_on_central)
3. central_pm = PatchMatch { start_a: nt.cw_anchor_on_context,
       len: central_match_len, start_b: nt.cw_anchor_on_central,
       tile_id: nt.central_tile_id }
4. ctx_n = context.boundary_len()
5. augmented = clone(context)
   augmented.add_tile(central_pm)
   aug_n = augmented.boundary_len()
   gap_start = ctx_n - central_match_len  // = seg_len_old
   gap_len = central_seq.len() - central_match_len
   central_ptid = augmented.patch_tile_ids()[gap_start]
6. frontier_pos = 0  // gap_end, always a junction
7. candidates = compute_candidates_covering_position(augmented, 0)

8. dist_to_frontier = gap_start  // default: no junction found
   for i in (0..gap_start).rev():
       if augmented.is_junction(i):
           dist_to_frontier = gap_start - i - 1
           break
   // dist_to_frontier = number of context edges between the last
   // context junction and the gap boundary
```

### Phase 2: For each petal candidate

```
for petal_pm in candidates:

    // --- Step 2a: Glue petal to augmented ---
    trial = clone(augmented)
    if trial.add_tile(petal_pm).is_none(): skip
    petal_ptid = trial.patch_tile_ids().iter().max().unwrap()

    // --- Step 2b: Closed detection ---
    remaining_gap_count = count positions in trial.patch_tile_ids()
                           where value == central_ptid
    if remaining_gap_count == 0:
        results.push(Closed { tile_id: petal_pm.tile_id,
                              tile_offset: petal_pm.start_b })
        continue

    // --- Step 2c: Find remaining gap and petal edge on trial ---
    // Remaining gap = contiguous block of central_ptid edges
    new_gap_start, new_gap_len = find_contiguous_block(trial, central_ptid)
    new_gap_end = (new_gap_start + new_gap_len) % trial_n

    // Petal edge at the frontier junction: the petal edge adjacent to
    // the remaining gap, on the CW side.
    // = trial.edges()[(new_gap_end + trial_n - 1) % trial_n]
    //   but must verify it belongs to petal_ptid
    pos_before_gap = (new_gap_end + trial_n - 1) % trial_n
    // Walk CW from new_gap_end to find first non-gap edge
    // That edge should be the petal edge at the junction
    petal_edge_info = find_petal_edge_adjacent_to_gap(trial, new_gap_end,
                                                       central_ptid, petal_ptid)

    // --- Step 2d: Abstract VT seq update ---
    new_vt_seq = nt.vt_seq.clone()

    if dist_to_frontier > 0:
        // Frontier was beyond last context junction → new VT entry
        new_junc_aug_pos = gap_start - dist_to_frontier
        cw_edge = augmented.edges()[(new_junc_aug_pos + aug_n - 1) % aug_n]
        // cw_edge is the context edge CW of the new junction position
        new_vt = VertexType { cw: cw_edge, inner: vec![], ccw: petal_edge_info }
        new_vt_seq.push(new_vt)
    else:
        // Frontier was at the last context junction → extend last VT
        last_vt = new_vt_seq.last_mut()
        last_vt.inner.push(last_vt.ccw)
        last_vt.ccw = petal_edge_info

    // --- Step 2e: Build new NT ---
    new_nt = NeighborhoodType {
        central_tile_id: nt.central_tile_id,         // unchanged
        cw_anchor_on_central: nt.cw_anchor_on_central, // unchanged
        cw_anchor_on_context: nt.cw_anchor_on_context, // unchanged (fixed)
        vt_seq: new_vt_seq,
    }
    results.push(Open(new_nt))
```

### Phase 2 detail: finding the petal edge adjacent to the gap

After gluing the petal, the trial boundary has three kinds of edges:
- Context edges (from the original context)
- Central tile edges (remaining gap, identified by `central_ptid`)
- Petal edges (identified by `petal_ptid`)

The remaining gap is a contiguous block of `central_ptid` edges. The petal edges
form a contiguous block adjacent to the gap (on the CW side). To find the petal
edge at the junction:

1. Find `new_gap_end = (new_gap_start + new_gap_len) % trial_n`
2. Walk CW from `new_gap_end`: position `(new_gap_end + trial_n - 1) % trial_n`
3. This position should have `patch_tile_id == petal_ptid`
4. Read `trial.edges()[pos]` → this is the petal EdgeInfo for the VT update

### Validation of new NTs

Every new NT produced by `explore_one` is validated by:
1. Reconstructing context from `new_vt_seq` via `construct_witness_from_vt_sequence`
2. Finding the match at `(cw_anchor_on_context, cw_anchor_on_central)`
3. Adding the central tile via `add_tile` — must succeed

This verifies that the abstract VT seq update produces a valid, reconstructible NT.

## Full BFS

```
fn new(tileset) -> NeighborhoodIndex:
    match_index = MatchTypeIndex::new(tileset)
    entries = []
    transitions = []
    seen = Map<NeighborhoodType, usize>
    queue = Deque<NeighborhoodType>

    // --- Seed Generation (Phase 1) ---
    for each match_type in match_index:
        patch = GrowingPatch::new(tileset, tile_a)
        patch.add_tile(PatchMatch { start_a, len, start_b, tile_id: tile_b })
        patch.normalize()
        vt_seq = extract junction vertex types from patch boundary

        for each third_tile_match in patch.get_all_matches():
            nt = NeighborhoodType { central_tile_id, cw_anchor_on_central,
                                    cw_anchor_on_context, vt_seq }
            if seen.contains(nt): skip
            seen.insert(nt, entries.len())
            entries.push(nt)
            queue.push_back(nt)

    // --- BFS Growth (Phase 2) ---
    while let Some(state) = queue.pop_front():
        state_id = seen[state]
        results = explore_one(state, match_index)

        for result in results:
            if result.outcome == Closed:
                transitions.push(NtTransition { src_id: state_id,
                    dst_id: NT_CLOSED_ID,
                    tile_id: result.tile_id,
                    tile_offset: result.tile_offset })
            else if Open(new_nt):
                dst_id = seen.get(new_nt), or:
                    new_id = entries.len()
                    seen.insert(new_nt, new_id)
                    entries.push(new_nt)
                    queue.push_back(new_nt)
                    dst_id = new_id

                transitions.push(NtTransition { src_id: state_id,
                    dst_id, tile_id: result.tile_id,
                    tile_offset: result.tile_offset })

    return NeighborhoodIndex { tileset, entries, transitions }
```

## Incremental Implementation Plan

### Stage A: Minimal API shell ✓
### Stage B+C: Seed generation + BFS skeleton ✓
  - Three-tile seeds, vt_seq extraction, BFS dequeue loop.
  - Counts: square=240, hex=1008.

### Stage D: Explore step infrastructure ✓
  - `attach_central`: reconstruct context, attach central, compute gap.
  - `find_gap_frontier`: find frontier junction at gap_end.
  - `explore_step`: combine both, enumerate candidates.
  - Tested on all square and hex seeds.

### Stage E: Full explore step (current)
  - E.1: `ExploreResult` / `ExploreOutcome` types
  - E.2: Closed detection helper
  - E.3: Find remaining gap + petal edge on trial
  - E.4: `dist_to_frontier` computation on augmented boundary
  - E.5: Abstract VT seq update (case A: extend last VT)
  - E.6: Abstract VT seq update (case B: create new VT)
  - E.7: `explore_one` combining everything
  - E.8: Validation test: all results reconstruct and re-attach
  - E.9: Wire into BFS loop with dedup

### Stage F: Remaining API
  - `classify_all()` — reuse existing classification logic.
  - `validate()` — reconstruct + check well-formedness.
  - `write_collection()` + `parse_file()`.
  - Full integration — `cargo test --release` including spectre.
