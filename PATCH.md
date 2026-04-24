# GrowingPatch: Patch Boundary with Segment/Vertex Tracking

## Overview

`GrowingPatch` is a data structure for incrementally building tile patches while
tracking which original tile each boundary segment belongs to. This enables:

1. Mapping any boundary position back to a position on a seed tile
2. Determining signed match types for new gluing operations
3. Tracking vertex types (ordered lists of signed match IDs around each vertex)

The structure lives in `src/intgeom/patch.rs`.

## Data Model

```rust
pub struct BoundarySegment {
    pub tile_id: usize,      // which seed tile this segment came from
    pub tile_start: usize,   // start position ON the seed tile (in seq() order)
    pub tile_len: usize,     // number of boundary edges (= number of angles)
}

pub struct VertexInfo {
    pub matches: Vec<i32>,   // ordered signed match type IDs (CCW around vertex)
    pub is_open: bool,       // true = gap remains, false = fully surrounded
}

pub struct VertexType(pub Vec<i32>);  // standalone key for later canonicalization

pub struct PatchMatch {
    pub start_a: usize,      // start of matched range on patch boundary
    pub len: usize,          // match length
    pub start_b: usize,      // start on seed tile (in seq() order)
    pub tile_id: usize,      // which seed tile
}

pub struct AddTileDiff {
    pub affected_old: Vec<VertexInfo>,
    pub new_vertices: Vec<VertexInfo>,  // [cw_junction, closed..., ccw_junction]
}

pub struct GrowingPatch<T: IsComplex> {
    tileset: Arc<TileSet<T>>,
    match_type_index: Arc<MatchTypeIndex<T>>,
    angles: Vec<i8>,                    // boundary in glue-order (non-canonical)
    segments: Vec<BoundarySegment>,
    vertices: Vec<VertexInfo>,
    matches: Vec<PatchMatch>,           // cached valid matches
}
```

## Invariants

- `angles.len() == segments.iter().map(|s| s.tile_len).sum()`
- `vertices.len() == segments.len()` (one vertex at the start boundary of each segment)
- Vertex `i` sits between the end of segment `i-1` and the start of segment `i`
  (cyclically: vertex 0 between segment `n-1` and segment 0)
- Segment `k` covers `angles[sum(0..k of tile_len) .. sum(0..k+1 of tile_len)]`
- `angles` is in glue-order (not canonical). Use `to_rat()` to canonicalize.
- `matches` is recomputed after each `add_tile` call.

## Constructors

### `GrowingPatch::new(tileset, match_type_index)`

Empty patch: `angles = []`, `segments = []`, `vertices = []`, `matches = []`.

### `GrowingPatch::from_tile(tileset, match_type_index, tile_id)`

Single tile as full boundary:
- `angles = tile.seq().to_vec()`
- `segments = [BoundarySegment { tile_id, tile_start: 0, tile_len: angles.len() }]`
- `vertices = []` (no matches yet, no real vertices)
- `matches` computed via MatchFinder

Adding the first tile to an empty patch is equivalent.

## Match Finding

### `get_all_matches() -> &[PatchMatch]`

Returns the cached list of all valid matches between the current patch boundary
and all tiles in the seed tileset.

### `get_matches(tile_id: usize) -> impl Iterator<PatchMatch>`

Filters cached matches to those that add a specific seed tile.

### Match computation (internal)

1. Create a single-tile `TileSet` from `angles`
2. Use `MatchFinder::crossing(patch_ts, seed_ts)` to find all valid matches
3. Cache the results

Note: the `angles` are NOT canonical. The MatchFinder handles this correctly
since `TileSet::new` will canonicalize internally via `Rat::from_slice_unchecked`.
The match positions are relative to the non-canonical `angles` sequence.

## `add_tile(patch_match) -> AddTileDiff`

Given a `PatchMatch` from `get_all_matches()` or `get_matches()`:

### Step 1: Identify affected segments

Find which segments overlap the matched range `[start_a, start_a+len)` on the
patch boundary. The matched range may span multiple segment boundaries.

### Step 2: Compute glue

Use the same logic as `Rat::try_glue_with_match_info` to produce new `angles`:
- `x = angles[start_a+len .. start_a-1]` (patch unmatched, cyclically)
- `y = seed_seq[start_b .. start_b-1]` (tile unmatched, cyclically)
- Junction angles at the two transition points
- New `angles = x[..-1] ++ y[..-1]` with junction angle fixes

### Step 3: Update segments

For each old segment that overlaps the matched range:
- If fully inside the matched range: removed (becomes interior)
- If partially overlapped: trimmed (start or end adjusted)
- The unmatched portion keeps its original (tile_id, tile_start) mapping

Append a new segment for the unmatched part of the added tile:
- `tile_id = patch_match.tile_id`
- `tile_start = (start_b + len) % seed_len`  (unmatched starts after match end)
- `tile_len = seed_len - len`

### Step 4: Update vertices

Determine which vertices fall within the matched range:
- Vertices strictly interior to the match: closed (is_open = false)
- Vertices at the boundary of the match: updated (one side matched, one not)

Create/update the two junction vertices:
- CW junction: where patch unmatched meets new tile unmatched (at match end)
- CCW junction: where new tile unmatched meets patch unmatched (at match start)

For each junction vertex, determine the signed match type by looking up the
boundary segment at the junction and resolving via `MatchTypeIndex::signed_id`.

### Step 5: Recompute matches

Call the internal match computation to update the cached `matches`.

### Step 6: Return diff

```
AddTileDiff {
    affected_old: vertices that were modified or closed,
    new_vertices: [cw_junction, closed_vertices..., ccw_junction]
}
```

## Signed Match Type Resolution

When a new tile is attached at a junction, we need to determine the signed
match type ID. The match is between:

- The boundary segment at the junction (which came from some original tile at
  some position)
- The new tile

We look up the boundary segment to get `(tile_id, tile_start)`, then compute
the match on the original seed tiles and look up the signed ID in
`MatchTypeIndex`.

For junction vertices that are updates of existing vertices, we append the new
signed match ID to the existing `matches` list.

## `to_rat() -> Rat<T>`

Canonicalizes `angles` via `Rat::from_slice_unchecked`, losing all
segment/vertex metadata. Returns a standard `Rat<T>`.

## Example: Hexagon 3-patch

### Step 1: from_tile(0)

```
angles = [2, 2, 2, 2, 2, 2]
segments = [(tile_id=0, tile_start=0, tile_len=6)]
vertices = []
```

### Step 2: add_tile(match at start_a=0, len=1, start_b=0, tile_id=0)

New angles (bi-hexagon, 10 positions):
```
angles = [-2, 2, 2, 2, 2, -2, 2, 2, 2, 2]
```
(junction angles at positions 0 and 5 are normalize(2+2-6) = -2)

Segments:
```
seg 0: (tile_id=0, tile_start=1, tile_len=4)  // old tile unmatched: pos 1-4
seg 1: (tile_id=0, tile_start=1, tile_len=4)  // new tile unmatched: pos 1-4
```

Wait — need to check tile_start for new tile. The match on the new tile is at
start_b=0, len=1. The unmatched part starts at (start_b+len)%6 = 1, length 5.
But the boundary only gets 4 positions (5 unmatched minus 1 shared endpoint).
So tile_start=1, tile_len=4.

Vertices:
```
vtx 0 (between seg 1 and seg 0): matches=[+m], is_open=true
vtx 1 (between seg 0 and seg 1): matches=[-m], is_open=true
```

Where m is the signed match type ID for the first gluing.

### Step 3: add_tile(match at start_a=4, len=2, start_b=0, tile_id=0)

This match spans the junction at position 5 (which is vtx 1). On the seed tile,
the corresponding match is at (tile_start=5, len=1) (since the junction doesn't
exist on the original tile).

After gluing: tri-hexagon. Some vertices may close.

## Tests

1. **hexagon_add_tile**: Start with hex, add one tile, verify segment/vertex counts
2. **hexagon_three_close**: Add three hex tiles, verify a closed vertex appears
3. **square_add_tile**: Start with square, add one tile
4. **square_four_close**: Add four squares, verify closed vertex
5. **mixed_hex_square**: Tileset with hex+square, verify cross-tile matches work
6. **match_type_ids**: Verify signed match IDs are correct at each step
7. **diff_tracking**: Verify AddTileDiff contains correct old/new vertices
