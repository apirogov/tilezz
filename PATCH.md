# GrowingPatch: Patch Boundary with Segment/Vertex Tracking

## Overview

`GrowingPatch` is a data structure for incrementally building tile patches while
tracking which original tile type each boundary segment belongs to. This enables:

1. Mapping any boundary position back to a position on a seed tile
2. Determining signed match types for new gluing operations
3. Tracking vertex types (ordered lists of signed match IDs around each vertex)

The structure lives in `src/intgeom/patch.rs`.

## Lifecycle

`GrowingPatch` has three states via an internal `PatchState` enum:

1. **Seed**: Created from a tileset + seed tile ID. Has cached matches but no
   boundary, segments, or vertices. The user picks a match from `get_all_matches()`.
2. **Growing**: After the first `add_tile`, the patch has a proper cyclic boundary
   with segments and open vertices. Each subsequent `add_tile` grows the patch.
3. **Closed**: All vertices have closed (no open vertices remain). No more tiles
   can be added.

Constructor: `GrowingPatch::new(tileset, match_type_index, seed_tile_id)`

## Data Model

```rust
pub struct BoundarySegment {
    pub tile_id: usize,      // which tile TYPE in the tileset
    pub tile_start: usize,   // start position ON that tile (in seq() order)
    pub tile_len: usize,     // number of boundary edges
}

pub struct VertexInfo {
    pub matches: Vec<i32>,   // signed match type IDs at this vertex
    pub is_open: bool,       // true = gap remains, false = fully surrounded
}

pub struct PatchMatch {
    pub start_a: usize,      // start of matched range on patch boundary
    pub len: usize,          // match length
    pub start_b: usize,      // start on target tile (in seq() order)
    pub tile_id: usize,      // which target tile
}

pub struct AddTileDiff {
    pub closed_vertices: Vec<VertexInfo>,  // vertices that became interior
}
```

## Invariants (Growing state only)

- `segments.len() == vertices.len()` (cyclic: one vertex per segment)
- Vertex `i` is at the boundary between segment `i-1` and segment `i`
  (cyclically: vertex 0 between segment `n-1` and segment 0)
- `angles.len() == segments.iter().map(|s| s.tile_len).sum()`
- `angles` is in glue-order (not canonical). Use `to_rat()` to canonicalize.
- Only **open** vertices are stored. Closed vertices are emitted in `AddTileDiff`
  and then forgotten.
- Segments from the same tile type with contiguous tile positions are **merged**
  into a single segment. This means a bi-hexagon (same tile type on both sides)
  can produce 1 or 2 segments depending on whether tile positions are contiguous.

## Segment Construction (Position Walk)

When rebuilding segments after `add_tile`, we walk through each new boundary
position in order and determine which original tile each position came from:

1. For positions in the old portion (0..seg_len_old): map back to the old
   boundary position, find which old segment contains it, compute the tile
   position within that tile type.
2. For positions in the new tile portion (seg_len_old..end): tile positions
   are `(start_b + offset) % tile_len`.

Adjacent positions with the same `tile_id` and contiguous tile positions
(`tile_pos == (prev + 1) % tile_len_total`) are merged into one segment.

This guarantees correct segment ordering and correct `tile_start` values.

## Junction Vertices

After a tile is added, two junction vertices get signed match type IDs:
- **CW junction** (vertex 0, at position 0): gets `-signed_id`
- **CCW junction** (vertex `num_old_out`, at position `seg_len_old`): gets `+signed_id`

Where `num_old_out` is the number of old-derived segments in the output
(segments covering positions 0..seg_len_old).

If all segments merge into one (CW == CCW), vertex 0 gets both IDs.

## Signed Match Type Resolution

When adding a tile, the signed match type ID is resolved by:
1. Finding the boundary segment at position `pm.start_a`
2. Computing `orig_pos = (seg.tile_start + offset_in_seg) % tile_len`
3. Calling `seed_a.get_match((orig_pos, pm.start_b), seed_b)` on the original
   seed tiles
4. Looking up the resulting `MatchType` in `MatchTypeIndex::signed_id()`

This maps patch-level matches back to canonical seed-level match types.

## `to_rat() -> Rat<T>`

Creates a `Rat` from `angles` via `Rat::from_slice_unchecked`, which
canonicalizes internally. All segment/vertex metadata is lost.

## Example: Bi-Hexagon

### Step 1: Seed phase

```
gp = GrowingPatch::new(ts, mti, 0)  // hex tile_id=0
gp.get_all_matches()  // 36 matches (6×6)
```

### Step 2: First add_tile (start_a=0, start_b=0)

```
gp.add_tile(&matches[0])
```

Boundary (10 positions):
```
angles = [-2, 2, 2, 2, 2, -2, 2, 2, 2, 2]
```

Segments: same tile type on both sides, positions 1-5 then 0-4 are contiguous
(5→0 wraps mod 6), so they merge:
```
seg 0: (tile_id=0, tile_start=1, tile_len=10)
```

Vertex:
```
vtx 0: matches=[+m, -m]  (both junctions at the same vertex)
```

### Step 3: First add_tile (start_a=0, start_b=1)

Tile positions 1-5 then 1-5 are NOT contiguous (5→1 ≠ 0 mod 6), so:
```
seg 0: (tile_id=0, tile_start=1, tile_len=5)   // old portion
seg 1: (tile_id=0, tile_start=1, tile_len=5)   // new portion
```

Vertices:
```
vtx 0 (CW junction):  matches=[-m]
vtx 1 (CCW junction): matches=[+m]
```

## Tests

1. `seed_patch_*`: Verify seed phase has matches but no boundary
2. `first_add_*`: Verify first add_tile initializes correctly
3. `square_all_16_matches_produce_valid_bi_squares`: All 4×4 matches
4. `hexagon_all_36_matches_produce_valid_bi_hexes`: All 6×6 matches
5. `to_rat_matches_direct_glue_for_all_matches`: Glue correctness for every match
6. `segment_cyclic_invariant`: segments.len() == vertices.len() always
7. `mixed_hex_square_add_tile`: Cross-tile matches
8. `junction_vertex_ids_nonempty_after_each_add`: Multi-step growth
