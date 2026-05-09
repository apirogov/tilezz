# New NT BFS Design

This document describes the replacement neighborhood-type (NT) enumeration
algorithm. The current `neighborhood.rs` BFS should be treated as a flawed
prototype, not as an implementation model.

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

## Why Replace The Current BFS

The existing BFS stores the full augmented boundary, `gap_start`, `gap_len`, and
deferred update state. That makes the implementation fragile and has already led
to fundamental bugs.

Known issues:

- `compute_new_edges_seed` and `compute_new_edges_growing` use the wrong offset
  convention for the newly glued tile after the first surviving edge. The correct
  offsets are `start_b, start_b + 1, ...`, not `start_b, start_b + len + 1, ...`.
- Gap tracking is spread across mutable boundary rotations and deferred updates.
  This is exactly the bookkeeping the new algorithm avoids.
- The old deferred-update mechanism models "same frontier, changed VT" as a
  side case. That is the wrong abstraction: same-frontier growth should update
  the active VT slot directly.
- Canonicalizing the augmented boundary risks quotienting states at the wrong
  level. The new BFS intentionally deduplicates by `(central_tile_id,
  cw_anchor_on_central, vt_seq)`.

Do not fix the old BFS incrementally unless explicitly choosing to abandon this
design. The intended implementation is a rewrite around anchored VT-sequence
states.

## Direction And Edge Conventions

Boundary arrays use the existing project convention:

- `angles[i]` is the angle at vertex `i`;
- `edges[i]` is the edge leaving vertex `i` in the CCW traversal direction;
- `EdgeInfo { tile_id, tile_offset }` identifies that boundary edge on its source
  tile;
- `vertex_type_raw_from(edges, inner_chains, pos)` reads `cw = edges[pos - 1]`,
  `inner = inner_chains[pos]`, `ccw = edges[pos]`.

The context VT sequence is ordered CW to CCW from the context point of view.
`vt_seq.last()` is the active frontier VT when the sequence is non-empty.

The central tile is attached to the context by a fixed CW anchor. From the
context point of view, the covered central-tile boundary extends CCW away from
that anchor. From the central tile point of view, the remaining gap starts at the
stored central anchor edge.

## BFS State

```rust
struct NtBfsState {
    central_tile_id: usize,

    // Fixed central-tile gap edge at the CW anchor. This is the `other_end` /
    // `start_b` value when gluing the central tile as the "other" boundary.
    cw_anchor_on_central: usize,

    // Context VTs ordered CW -> CCW. Empty only for the two-tile seed state
    // before the first context frontier VT has been created.
    vt_seq: Vec<VertexType>,

    // Used only to reconstruct empty vt_seq states.
    seed_tile_id: usize,
    cw_junc_on_seed: usize,
    seed_match_len: usize,

    // Fixed CW junction on the reconstructed context boundary. This is the
    // `self_start` value when gluing the central tile to the context.
    // Invariant: this index is in the coordinate system returned by the context
    // reconstruction for this state's vt_seq.
    cw_junc_on_ctx: usize,
}
```

Deduplication key:

```rust
(central_tile_id, cw_anchor_on_central, vt_seq)
```

`cw_junc_on_ctx` is not part of the semantic dedup key. It is anchoring metadata
needed to reconstruct and grow the state. The algorithm assumes the VT sequence
helper reconstructs the same minimal context witness coordinate system produced
by replay growth. Tests must check this invariant.

## Seed States

Each unique valid maximal glue of two tiles induces two BFS starting states.

Given a valid glue of tile `A` to tile `B` with match data:

```text
start_a on A, start_b on B, len
```

create:

- central `A`, context seed `B`:
  `central_tile_id = A`, `seed_tile_id = B`,
  `cw_anchor_on_central = start_a`, `cw_junc_on_seed = start_b`,
  `cw_junc_on_ctx = start_b`, `seed_match_len = len`, `vt_seq = []`.
- central `B`, context seed `A`:
  `central_tile_id = B`, `seed_tile_id = A`,
  `cw_anchor_on_central = start_b`, `cw_junc_on_seed = start_a`,
  `cw_junc_on_ctx = start_a`, `seed_match_len = len`, `vt_seq = []`.

There is no extra normalization. The two junctions of the glued pair can be read
in either direction: one reading treats `A` as the center, the other treats `B`
as the center.

Seed validity requires the actual two-tile glue to be valid:

- `glue_raw_angles` / `compute_glue_angles` succeeds;
- no degenerate `+/-hturn` junction is produced;
- `Snake::try_from(new_angles)` succeeds;
- `snake.is_closed()` is true.

If `len` covers the entire chosen central tile, the result is already closed and
is not enqueued as an open NT.

## Context Reconstruction

For a dequeued state:

- if `vt_seq.is_empty()`, the context is the full boundary of `seed_tile_id`,
  with `cw_junc_on_ctx = cw_junc_on_seed`;
- otherwise call
  `GrowingPatch::construct_witness_from_vt_sequence(&vt_seq, match_index)` and
  use the returned patch as the context.

Do not modify `construct_witness_from_vt_sequence` as part of the NT rewrite.
If its coordinate-system invariant proves false, add a thin anchored wrapper or
additional return data; do not reintroduce augmented-boundary canonicalization.

## Central Glue

Attach the central tile to the reconstructed context at the fixed CW anchor.

Use the context as `self` and the central tile as `other`:

```rust
let ctx_rat = Rat::from_slice_unchecked(ctx_angles);
let central_rat = tileset.rat(central_tile_id);
let (ctx_match_start, covered_len, central_gap_start) = ctx_rat.get_match(
    (cw_junc_on_ctx as i64, cw_anchor_on_central as i64),
    central_rat,
);
```

State invariant:

```text
ctx_match_start == cw_junc_on_ctx
central_gap_start == cw_anchor_on_central
```

If this invariant fails, the stored anchor does not describe the CW-most match
endpoint in the reconstructed context. That is a bug in state construction or in
the reconstruction coordinate invariant; do not silently scan for another anchor.

Build the augmented boundary by gluing the central tile to the context with this
match. The central unmatched gap is the surviving central-tile portion beginning
at `cw_anchor_on_central`. The open NT boundary and `gap_len` are derived from
this augmented boundary:

```text
gap_len = central_tile_edge_count - covered_len
```

For a dequeued state, `gap_len` must be positive. During growth, if a newly
replayed context reattaches with `gap_len == 0`, record a transition to
`NT_CLOSED_ID` instead of creating an open NT.

## Finding The Frontier

The active frontier is determined on the augmented boundary, not by assuming a
fixed raw index such as `0`.

Procedure:

1. Identify the fixed CW junction on the augmented boundary, i.e. the junction
   where the central unmatched gap starts.
2. Walk CCW along the unmatched central-tile gap.
3. Stop at the next actual junction. That junction is the CCW frontier.

An "actual junction" should use the same definition as the patch candidate code:
the boundary angle at a position differs from the source tile's own angle at the
corresponding `EdgeInfo` offset.

The equivalent context-side frontier before adding a petal is:

```text
ctx_frontier = (cw_junc_on_ctx + covered_len) mod ctx_boundary_len
```

This context-side position is the replay point for growing the context without
the central tile.

## Frontier Distance

Raw indices are unstable because every glue operation reindexes the boundary.
To decide whether the active frontier advanced, compare an anchored distance.

For a context boundary with fixed CW anchor `cw` and current context frontier
`frontier`, define:

```text
frontier_distance = (cw + context_boundary_len - frontier) mod context_boundary_len
```

This is the CW distance from the fixed CW anchor to the CCW-most frontier along
the unmatched side of the context.

After adding a petal and replaying it on the context, glue the central tile to
the new context again and recompute this distance.

- If the distance did not change, the central-tile match did not advance. The
  old frontier is still the active frontier.
- If the distance changed, the old frontier was closed/absorbed by the central
  tile match. The new CCW-most frontier is a new VT slot.

Do not compare raw boundary indices across different boundaries.

## Petal Enumeration

For the current augmented frontier, enumerate candidate petals by tile id and
tile junction offset.

A candidate is valid only if it touches the frontier. It does not need to match a
central-tile edge immediately. It may merely refine the same frontier junction.

For each candidate:

1. Compute the local forward match at the augmented frontier.
2. Glue the petal to the augmented boundary.
3. Reject invalid glue:
   - empty result;
   - degenerate `+/-hturn` junction;
   - `Snake::try_from` failure;
   - `!snake.is_closed()`.
4. Enforce the frontier angle monotonicity invariant for same-frontier
   refinements. The comparison must use the project's normalized signed angle
   convention; do not substitute an absolute-value comparison unless separately
   justified and tested.

The implementation should expose/reuse a single boundary-glue helper for this.
Do not copy the old `compute_new_edges_*` logic. The correct new-tile edge
offsets after a glue are:

```text
tile_junc, tile_junc + 1, tile_junc + 2, ...
```

modulo the new tile's edge count.

## Replay On The Context

A petal validated on the augmented boundary must be replayed on the original
context boundary, without the central tile.

Replay uses the same `(tile_id, tile_offset)` at the context-side frontier:

```text
ctx_frontier = (cw_junc_on_ctx + covered_len) mod ctx_boundary_len
```

The match length on the context may differ from the augmented match length. That
is expected. In particular, if the central tile would have closed the old
frontier on the augmented boundary, the context-only replay still has no central
tile and may expose a different local junction structure.

The replay step must return:

- the new context boundary (`angles`, `edges`, `inner_chains`);
- the remapped `cw_junc_on_ctx` for the new boundary;
- the empirically observed context frontier VT after reattaching the central
  tile to the new context, unless the central tile is closed.

The fixed CW anchor remaps through ordinary `glue_raw_angles` indexing. If the
context replay starts at `ctx_frontier` with match length `mlen`, then old context
positions that survive are reindexed from:

```text
ccw_pos = (ctx_frontier + mlen) mod old_context_len
new_pos = (old_pos + old_context_len - ccw_pos) mod old_context_len
```

The remapped `cw_junc_on_ctx` must land inside the surviving old-context segment.
If it does not, the replay consumed the fixed anchor and is invalid.

## Updating The VT Sequence

After replaying the petal on the context, reattach the central tile to the new
context and compute the new frontier distance.

Then update `vt_seq`:

- If `vt_seq` is empty, append the observed frontier VT. This is the first
  context frontier VT after the two-tile seed.
- Else if the frontier distance did not change, replace `vt_seq.last()` with the
  observed frontier VT. We are still at the same frontier junction; only its VT
  changed.
- Else append the observed frontier VT. The old frontier was closed/absorbed by
  the central tile match and the new frontier is a new VT slot.

This is the central distinction from the old BFS. Same-frontier growth is an
update of the active VT, not a deferred side branch.

If reattaching the central tile to the new context covers the full central tile,
record a closed transition and do not enqueue a new open state. No new frontier
VT is needed for a closed destination.

## Recording NT Entries

Each non-closed BFS state corresponds to one `NeighborhoodType` entry derived by
gluing the central tile to the reconstructed context.

Fields should be derived, not incrementally tracked:

- `central_tile_id` from the state;
- `central_anchor_edge = cw_anchor_on_central`;
- `gap_len = central_tile_edge_count - covered_len`;
- `angles`, `edges`, `inner_chains`, and `gap_start` from the augmented boundary;
- `context_vertex_types = vt_seq`.

This avoids carrying stale gap bookkeeping through growth.

## Transition Graph

Each successful petal produces exactly one transition:

- to another open NT state if the central tile remains open;
- to `NT_CLOSED_ID` if the central tile becomes fully covered.

The transition should record the petal `(tile_id, tile_offset)` and useful match
metadata for diagnostics. Classification (`Dead`, `Undead`, `Blessed`, `Free`)
can reuse the existing `classify_all` logic once transitions are correct.

## Required Tests

The rewrite needs tests that validate the algorithmic invariants, not merely
internal consistency of serialized output.

Required coverage:

- seed generation creates both orientations for each valid two-tile glue;
- central glue anchor invariant holds for every dequeued state;
- `cw_junc_on_ctx` survives and remaps correctly after every context replay;
- same-frontier petal growth replaces the last VT;
- frontier-advance petal growth appends a new VT;
- empty-seed first petal appends the first VT;
- closed transitions are recorded when central coverage reaches the full central
  tile edge count;
- all produced NT entries pass `NeighborhoodType::validate`;
- square, hexagon, and spectre enumerate and validate in release mode.

Run, at minimum:

```bash
cargo fmt
cargo clippy
cargo test --release
```

Fix all issues before committing.
