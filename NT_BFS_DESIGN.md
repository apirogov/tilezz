# NT BFS Design

This document describes the neighborhood-type (NT) enumeration in
`src/intgeom/neighborhood.rs`. It supersedes the original prototype design
and is kept in sync with the implementation.

## Concept

An open NT is a local configuration around a central tile:

- one central tile of a known type;
- a finite context patch that shares one contiguous boundary segment with
  the central tile;
- the remaining uncovered "gap" segment of the central tile (must be
  non-empty for the NT to be *open*);
- the context patch is described by a sequence of vertex types (VTs)
  along the *covered* segment.

Closed configurations (= the gap segment is empty) are not stored as
entries. They are referenced via the sentinel id `NT_CLOSED_ID = 0` in
transitions.

## Data Types

```rust
pub struct NeighborhoodType {
    pub central_tile_id: usize,
    // Tile offsets on the central tile at the two endpoints of the
    // matched segment (= the CW and CCW anchor *vertices*). The matched
    // edges on central are at offsets [cw_anchor_on_central - match_len,
    // cw_anchor_on_central - 1] going CW from cw_anchor_on_central,
    // because shared edges face opposite directions:
    //     CCW on context  <-->  CW on central.
    // So ccw_anchor_on_central = (cw_anchor_on_central + central_n
    //                             - match_len) % central_n.
    pub cw_anchor_on_central: usize,
    pub ccw_anchor_on_central: usize,
    // Frontier positions on context parameterised relative to vt_seq:
    //   cw_anchor_on_context  = CW dist from vt_seq[0].junction
    //                           to cw_anchor vertex on context.
    //   ccw_anchor_on_context = CCW dist from vt_seq.last().junction
    //                           to ccw_anchor vertex on context.
    // Both stay small (the vt_seq's endpoints are the junctions closest
    // to the anchors that lie *outside* the covered segment).
    pub cw_anchor_on_context: usize,
    pub ccw_anchor_on_context: usize,
    // Number of context tiles glued. Increments by exactly 1 per BFS
    // step. Gives the natural DAG layering and bounds for matching pairs
    // in the structural CW reconstruction.
    pub num_ctx_tiles: usize,
    // CCW-ordered junctions on the context boundary that lie inside the
    // covered segment (= "incident with the match"). vt_seq[0] is the
    // first junction reachable CCW from the CW anchor vertex; vt_seq
    // .last() is the last junction reachable CW from the CCW anchor
    // vertex.
    pub vt_seq: Vec<VertexType>,
}

pub const NT_CLOSED_ID: usize = 0;

pub enum NtKind { Dead, Undead, Blessed, Free }

pub struct NtTransition {
    pub src_id: usize,
    pub dst_id: usize,        // NT_CLOSED_ID == 0 for closed transitions
    pub side: TransitionSide, // Cw / Ccw — which frontier the petal grew
    pub tile_id: usize,       // petal tile id
    pub tile_offset: usize,   // petal_pm.start_b
}
```

An NT is uniquely identified by all fields of `NeighborhoodType` together.
IDs are 1-based: `entries[i]` has id `i + 1`. The sentinel id `0` only
appears as `NtTransition::dst_id` (never as `src_id`).

`match_len` is **derived**:
`(cw_anchor_on_central + central_n - ccw_anchor_on_central) % central_n`.

## Coordinate conventions

- Boundary arrays use the project convention: `angles[i]` is the angle at
  vertex `i`, `edges[i]` is the edge leaving vertex `i` going CCW.
- `vt_seq` is in CCW order: walking the context boundary CCW from the CW
  anchor visits `vt_seq[0]`, then `vt_seq[1]`, …, then `vt_seq.last()`
  before reaching the CCW anchor.
- "CW distance from vertex A to vertex B" means going CW (= backward in
  the array) from A, the number of edges until B. Equivalently, the CCW
  distance from B to A.

## BFS

The algorithm:

1. **Seed generation.** For each match-type `(tile_a, start_a, tile_b,
   start_b, len)` from `MatchTypeIndex`:
   - Build a two-tile patch by gluing `tile_b` to `tile_a` along the
     match.
   - Normalize the boundary (`lex_min_rot`).
   - For every `third_pm` returned by `patch.get_all_matches()`:
     - Try gluing the third tile to the two-tile patch on a clone; skip
       if the glue collides (non-convex tiles like spectre can fail this
       even when the edge match is valid).
     - Collect all junctions in the normalized patch.
     - Filter to those incident with the match (CCW distance from
       `third_pm.start_a` is at most `third_pm.len`) — that becomes
       `vt_seq`.
     - Compute the anchor fields with Convention 2 (see "Anchors" below).
     - Run `try_construct_nt_from_cw` (which also validates the abstract
       NT reconstructs); skip on failure.
     - Dedup against `seen`; if new, push to `entries` and `queue`.

2. **BFS growth.** Dequeue each state and call `explore_one` to enumerate
   one-petal successors. For each open successor, dedup against `seen`,
   enqueue if new, push a transition with the right side.

3. **Classification.** `classify_all` walks the transition graph and
   assigns `NtKind` per entry: `Dead` (no outgoing transitions),
   `Undead` (all outgoing reach Dead/Undead via the same), `Blessed`
   (every outgoing leads to Closed or another Blessed), `Free`
   otherwise.

## Anchors

The CW anchor *vertex* on context is at position
`(first_junc + ctx_n - cw_anchor_on_context) mod ctx_n` where
`first_junc` is the position of `vt_seq[0]`'s junction in the
reconstructed ctx.

The CCW anchor *vertex* on context is at position
`(last_junc + ccw_anchor_on_context) mod ctx_n` where `last_junc` is
`vt_seq.last()`'s position.

Match length on context = `(ccw_anchor_pos - cw_anchor_pos) mod ctx_n`.

On the central tile, `cw_anchor_on_central` is `start_b` from
`PatchMatch` (= the tile offset at the CW anchor vertex). The CCW anchor
vertex on the *central* is at tile offset `start_b - match_len mod
central_n`, which is what `ccw_anchor_on_central` stores.

## Helpers

- `build_attached_context(nt, mi)` reconstructs ctx, attaches central,
  finds both frontiers on the augmented patch, and returns enough
  metadata for one BFS step (CCW and/or CW).

- `try_construct_nt_from_cw(central, cw_ac, cw_ax, num_ctx_tiles,
  vt_seq, mi) -> Option<NeighborhoodType>`: reconstruct from CW-side
  anchor fields, derive the CCW-side fields via `forward_match_length`,
  fail if reconstruction or the central glue fails.

- `try_construct_nt_from_ccw(central, ccw_ac, ccw_ax, num_ctx_tiles,
  vt_seq, mi) -> Option<NeighborhoodType>`: symmetric helper using
  `backward_match_length`.

- `nt_is_valid(nt, mi)`: reconstructs from the CW side and verifies the
  stored CCW-side fields agree (= the NT is self-consistent).

## `explore_one`

Returns a vector of `ExploreOutcome { side, petal_pm, kind }` where
`kind` is `Closed` or `Open { nt, trial, central_ptid }`.

For each candidate petal at the CCW frontier:
- Glue on a clone of `aug` (skip on collision).
- If the trial has no central edges remaining, emit `Closed`.
- Otherwise compute the new vt_seq:
  - **CCW Case A** (`frontier_is_junction_in_ctx`): extend `vt_seq[-1]`
    — push the old ccw into `inner`, set `ccw = petal_edge`.
  - **CCW Case B**: append a new VT with `cw = last_covered_ctx_edge`,
    `ccw = petal_edge`, `inner = []`.

  Where
  `petal_edge.tile_offset = (start_b - i_anchor) mod petal_n` and
  `i_anchor = (frontier_pos_on_aug - start_a) mod aug_n` (= the index of
  the OLD CCW anchor vertex within the petal's match span).

- Call `try_construct_nt_from_cw` to finalize the new NT (fills CCW-side
  anchors and validates).

## BFS

The state graph is a DAG: every transition increments `num_ctx_tiles` by
exactly 1.

CW-direction exploration is **not implemented**. The CCW-only BFS reaches
the complete CCW-reachable state set, verified by a brute-force
completeness test (`spectre_ccw_only_is_complete`): for every NT in the
BFS set, every direct `get_match` candidate at the CCW frontier lands
either on a closed configuration or on an NT already in the set.

## Roadmap

- ✅ NT struct with CW and CCW anchors on both central and context, plus
  `num_ctx_tiles`.
- ✅ Helper `try_construct_nt_from_cw`.
- ✅ CCW-only BFS verified complete via brute-force cross-check.
- 🚧 CW-direction `try_step_cw`. Previous attempts produced spurious NTs
  whose reconstructed boundary didn't match the actual `trial.angles()`
  boundary — the Case A/B vt_seq mutation breaks down when the OLD CW
  anchor's junction status flips between old ctx and new ctx (= when the
  petal's angle contribution at that vertex is zero). Removed pending a
  clean re-derivation.
- 🚧 Performance check on square (109k entries, ~1M transitions). If
  classify_all becomes the bottleneck, switch its inner loop to a
  successor-set lookup instead of a per-iter filter over transitions.

## Required passing tests

- Public API: `NeighborhoodType`, `NtTransition`, `NtKind`,
  `classify_all`, `validate`, `write_collection`, `parse_file`.
- `square_seeds_validate`, `hex_seeds_validate`, `spectre_validates` —
  all entries reconstruct and accept the central tile.
- `square_roundtrip_collection`, `hex_roundtrip_collection` — write /
  parse round-trip equality on the produced collection.
- `bfs_produces_transitions` — id range invariants.
- `classify_all_returns_one_per_entry` — classification vector matches
  entry count and at least one Blessed entry exists.

## Build requirements

```
cargo fmt
cargo clippy -- -D warnings
cargo test --release
```
