# API design review — patch / tileset module

Scope: the public + `pub(crate)` API surface of the patch/tileset
subsystem of `intgeom`.

Files covered:
- `src/intgeom/patch.rs`
- `src/intgeom/tileset.rs`
- `src/intgeom/tiles.rs`
- `src/intgeom/matchtypes.rs`
- `src/intgeom/growing.rs`
- `src/intgeom/rat.rs` (the underlying `Rat` type is heavily used by everything above)

Findings are organized by severity. Each item includes context, the
concrete API impact, and a recommended fix.

A previous review of `patch.rs`'s internals (semantics, tests, refactor
opportunities) was completed in a prior pass; this review focuses on
**API surface, naming, scoping, documentation, and dead public surface**.

---

## Critical

### C1. Two distinct types both named `GrowingPatch`

Both modules define a type with the same name:

- `patch::GrowingPatch<T: IsComplex>` — the main public abstraction
  for incrementally grown, hole-free tilings (3,656-line `patch.rs`).
- `growing::GrowingPatch<T: HasPatchPos>` — an internal helper for
  Redelmeier polyomino enumeration (1,096-line `growing.rs`).

External code uses only the `patch.rs` one (5 files import it). The
`growing.rs` one has zero external callers. The name collision is
invisible at the call site because each call site imports just one,
but it surfaces in `cargo doc`, in cross-file code reading, and
whenever someone writes `use crate::intgeom::*`.

**Impact:** confusing for any reader doing module-level navigation;
hostile to crate-level documentation; ambiguous in error messages
that elide the path.

**Recommended fix:** rename `growing::GrowingPatch` to something that
reflects its role (`RedelmeierPatch`, `EnumerationPatch`, or
`PartialPatch`). It has no external callers, so the rename is purely
internal. ~5-minute change.

### C2. `growing.rs` has a near-zero exported surface but most items are `pub`

External use of `growing.rs`:
- `grow_redelmeier` — 1 external file (`bin/patch_bench`).
- `grow_redelmeier_free` — 1 external file (`bin/patch_bench`).

Everything else is `pub` but never referenced outside `growing.rs`:
- `GlueSite` (currently `pub(crate)`)
- `PatchPos` trait, `HasPatchPos` trait
- `Pos2`, `Pos4`, `Pos8` position-type structs
- `growing::GrowingPatch<T: HasPatchPos>` (see C1)
- `GrowingPatch::len`, `is_empty`, `from_seed`, `to_rat` methods
- `GrowStats` struct + its `Display` impl
- `grow_redelmeier_profiled`
- `make_free`

**Impact:** every `pub` item is part of the crate's public API, lands
in `cargo doc`, and constrains future internal refactors. A user
encountering this module via docs would have no way to tell the two
real entry points from the dozen incidental ones.

**Recommended fix:** downgrade everything except `grow_redelmeier`
and `grow_redelmeier_free` to `pub(crate)` (or private where the use
is single-file). Audit each entry: if there's a reason to keep it
`pub`, document the use case. ~15-minute change.

### C3. Pervasive lack of doc comments on public types

The following public types have **no documentation at all**:

- **`tileset.rs`**: `TileSet`, all five methods.
- **`matchtypes.rs`**: `MatchType` (incl. its non-obvious field
  semantics — see L1), `CandidateMatch`, `MatchFinder` + all 18 of
  its methods, `MatchTypeIndex` + all 9 of its methods, the two
  `pub(crate)` free functions.
- **`patch.rs`**: `EdgeInfo`, `PatchMatch`, `TransitionSide`, and many
  `GrowingPatch` methods (`new`, `is_growing`, `add_tile`,
  `boundary_len`, `angles`, `edges`, `inner_chains`, etc.).
- **`tiles.rs`**: 9 of 13 tile constructors (only `square`,
  `triangle`, `hexagon`, `spectre` have docs).
- **`growing.rs`**: every public item.
- **`rat.rs`**: many heavily-used methods.

For a published crate (`tilezz v0.0.3`) this is a significant gap.
Several of these types encode subtle invariants (start_a vs start_b
asymmetry, edge-offset conventions, rotation conventions) that are
not visible from the field names. The `OpenVertexType` doc-comment
pattern from the prior review pass is the right model: a paragraph
on what the type means, then explicit invariants where they exist.

**Recommended fix:** doc-pass across all public items.
Prioritize the heavily-used types (`TileSet`, `Rat` methods, `MatchTypeIndex`)
and the ones with non-obvious semantics (`MatchType`, `EdgeInfo`,
`PatchMatch`). Best done as a separate dedicated PR.

---

## High

### H1. `MatchFinder` exposes ~10 methods that nothing calls externally

External call sites (per the dependency survey):

| Used                              | Not used |
|-----------------------------------|----------|
| `new`                             | `set_a`, `set_b` |
| `crossing`                        | `rat_a`, `rat_b` |
| `valid_matches` (1 site)          | `apply_match` |
| `valid_matches_filtered` (1 site) | `shared_boundaries` |
| `valid_results_for_pairs`         | `valid_matches_with_rats` |
| `num_tiles_a` / `num_tiles_b`     | `all_valid_matches` |
|                                   | `valid_matches_for_pairs` |

That's roughly half the public surface of the struct unused.

**Impact:** each unused method is a maintenance burden — it must
compile, must be considered in API-evolution decisions, contributes
to module size, and clutters docs.

**Recommended fix:** delete the unused methods. If a future caller
needs one, the new call site will tell us its exact required shape.
Should also document the surviving ones (see C3, L2).

### H2. `tiles.rs` has 6 dead tetromino constructors

Defined but uncalled anywhere in production or tests:
`tetromino_O`, `tetromino_I`, `tetromino_S`, `tetromino_Z`,
`tetromino_J`, `tetromino_L`.

Only `tetromino_T` is used (one call site, in a `patch.rs` test).

**Impact:** small, but six dead public functions in a small file
(131 lines total) is a high dead-code density.

**Recommended fix:** delete the six unused constructors. If polyomino
test coverage is wanted as a public catalog, instead add a smoke test
that constructs each shape and verifies its angle sequence — that
keeps them alive with explicit purpose.

### H3. Internal `pub` and `pub(crate)` cluster in `patch.rs`

Despite the prior cleanup pass, these items are still visible beyond
their need:

| Item | Current visibility | External use | Recommended |
|------|---------------------|--------------|-------------|
| `update_inner_chains` | `pub` | 0 | private |
| `junction_angle_sequence` | `pub` | 0 | private |
| `compute_glue_angles` | `pub` | 0 | private |
| `is_junction_at` | `pub(crate)` | 0 | private |
| `vertex_type_raw_from` | `pub(crate)` `#[cfg(test)]` | 0 (tests in same file) | private+`#[cfg(test)]` |
| `forward_match_length` | `pub(crate)` | 1 (`neighborhood.rs`, 3 sites) | keep `pub(crate)`, add doc |
| `RawBoundary` | `pub(crate)` | 0 | private |
| `RawGlueResult` | `pub(crate)` | 0 | private |
| `glue_tile_to_raw_boundary` | `pub(crate)` | 0 | private |
| `glue_match_to_raw_boundary` | `pub(crate)` | 0 | private |
| `raw_is_junction` | `pub(crate)` `#[cfg(test)]` | 0 | private+`#[cfg(test)]` |
| `next_junction_on_raw_boundary` | `pub(crate)` `#[cfg(test)]` | 0 | private+`#[cfg(test)]` |

**Impact:** the `Raw*` cluster is internal infrastructure that the
witness-construction path leans on; nothing outside `patch.rs` uses
any of it. Keeping it `pub(crate)` advertises an internal contract to
the rest of the crate that doesn't exist.

**Recommended fix:** tighten each visibility to the minimum that
satisfies its actual callers.

### H4. `GrowingPatch` (patch.rs) has many `pub` methods with no external callers

These methods are `pub` but never invoked outside `patch.rs`:

`get_matches_for_tile`, `ensure_candidates_materialized`,
`next_tile_id`, `candidates_by_start`, `vertex_type_at`,
`junction_vertices`, `tile_segments`, `from_parts`,
`construct_minimal_witness`, `construct_witness_from_vt_sequence`,
`compute_candidates_covering_position`, `compute_all_candidates`.

These fall into three buckets:

1. **Intended public API** (likely): `from_parts`,
   `construct_minimal_witness`, `construct_witness_from_vt_sequence`,
   `compute_candidates_covering_position`, `compute_all_candidates`.
2. **Internal handles that leaked out**: `ensure_candidates_materialized`,
   `next_tile_id`, `candidates_by_start`.
3. **Test-only accessors**: `vertex_type_at`, `junction_vertices`,
   `tile_segments`.

**Impact:** category 2 and 3 should not be in the public API;
category 1 should be documented as such.

**Recommended fix:** triage each. Move category 2 to `pub(crate)` or
private. Move category 3 to `pub(crate)` or expose them only inside
the test module via `pub(super)`. Confirm category 1 is intended
public API and document.

---

## Medium

### M1. `Seed` / `Growing` state machine in `GrowingPatch` is leaky

`GrowingPatch` is conceptually two states (`Seed`, `Growing`) under
one type. Many accessors return defaults for `Seed` rather than
panicking or being unavailable: `boundary_len() → 0`,
`angles() → &[]`, `edges() → &[]`, etc. Callers must explicitly check
`is_growing()` to know which state they're dealing with.

**Impact:** type-level semantics are weak. A caller can hold a
`Seed`-state `GrowingPatch`, call `angles()`, get an empty slice,
and not realize until much later that they were treating an
unstarted seed as a real boundary.

**Recommended fix (optional, larger):** consider splitting into two
types: `SeedPatch` (just the seed metadata) and `GrowingPatch` (the
boundary state). The transition is
`SeedPatch::add_tile(self, pm) -> Option<GrowingPatch>`. This makes
the state machine a type-level guarantee and removes ~10 "return
default for Seed" branches.

Not urgent. Worth revisiting if a third state ever shows up (e.g.
"closed" patch), which would be the natural trigger for the proper
split.

### M2. `Rat` exposes 9 methods nothing uses externally

`Rat::cycle`, `cycled`, `slice`, `slice_from`, `is_convex`,
`is_canonical`, `match_length` (method form), `len`, `is_empty`.

These are likely small one-liners but each contributes to `Rat`'s
public surface area.

**Impact:** `Rat` is the central data type (18 external files use
it). Every public method is part of the crate's API surface; keeping
unused ones around makes future evolution harder and clutters docs.

**Recommended fix:** downgrade to `pub(crate)` where appropriate. If
a `pub` method is genuinely intended for external use, document it
and add an example test.

### M3. `TileSet::index_of` is tests-only

External use: only `matchtypes.rs` lines 1486 and 1488, both inside
the `mod tests` block.

**Recommended fix:** either gate with `#[cfg(test)]` and downgrade to
`pub(crate)`, or remove if its test use can be replaced with linear
search.

### M4. `PatchVertexInfo` is `pub` but unused externally

Returned by `vertex_type_at()` and `junction_vertices()`. Neither
method has external callers (see H4). The struct's `pub` status
serves no external caller today.

**Recommended fix:** downgrade `PatchVertexInfo` to `pub(crate)` (or
private) and the two methods to match.

### M5. `TileSegment` is `pub` but tests-only externally

`tile_segments()` is only invoked from `patch.rs`'s own tests. The
struct's `pub` status is unused.

**Recommended fix:** `pub(crate)` for both the struct and the
method, or keep `pub` if there's a planned external use case
(document it then).

### M6. `MatchType::apply` and several `MatchTypeIndex` methods have no external callers

Unused externally:
- `MatchType::apply`
- `MatchTypeIndex::tileset`
- `MatchTypeIndex::num_types`
- `MatchTypeIndex::signed_id`
- `MatchTypeIndex::apply`
- `MatchTypeIndex::all_candidates_for_tile`

`MatchTypeIndex::tileset` might be genuinely useful (an external
caller might need the underlying tileset from a stored index). Check
before dropping.

**Recommended fix:** `pub(crate)` or private after triage. Audit
each: removal is fine for the rest.

---

## Low

### L1. `MatchType` struct: field semantics undocumented

```rust
pub struct MatchType {
    pub tile_a: usize,
    pub start_a: usize,
    pub tile_b: usize,
    pub start_b: usize,
    pub len: usize,
}
```

The semantic asymmetry between `start_a` (first *matched* edge offset
on side A) and `start_b` (first *surviving* edge offset on side B,
just past the match end on B's side — see the rotation-convention
discussion in the prior `patch.rs` review) is non-obvious and not
documented anywhere.

**Recommended fix:** doc-comment on the struct explaining the
convention. Same comment can serve `PatchMatch` and `RawGlueResult`.

### L2. `MatchFinder::new` vs `crossing` distinction is undocumented

Both have the same return type. Naming hints at "self-match within
one tileset" vs "match between two tilesets" but isn't explicit.

**Recommended fix:** one-line doc on each.

### L3. `growing.rs` feature flag undocumented

The choice between `check_segments_cross_cyclotomic` (with feature
`cyclotomic_intersect`) and the `pos4`/`pos2`/`pos8`-based
intersection helpers (without it) is significant — different math
backends. Neither path documents what they do or why one would
prefer the feature on/off.

**Recommended fix:** module-level doc on `growing.rs` (or on the
relevant trait `HasPatchPos`) explaining the two backends, their
trade-offs, and when to enable the feature.

### L4. `growing::GrowStats::Display` impl is unused externally

The struct and its `Display` rendering of profile counters are `pub`
but never referenced outside `growing.rs`.

**Recommended fix:** `pub(crate)` for the struct + impl. Re-promote
if external profiling tooling is added.

---

## What's working well

A balanced review should call this out:

- **`Rat`** is the right central abstraction. Heavily used (18
  external files), most-touched API in the module. `try_glue`,
  `try_glue_precomputed`, `get_match` are clearly intentional.
- **`TileSet`** is well-scoped: a small read-only collection of
  `Rat`s with deterministic ordering, dedup, and chirality check.
  Five public methods, all used, all single-purpose.
- **`tiles`** entry points for `hexagon`, `square`, `spectre`,
  `triangle`, `penrose_p3_narrow`, `penrose_p3_wide` are clearly
  the intended public catalog. Just trim the dead tetrominos
  (H2).
- **`patch.rs`** (after the prior cleanup pass) has a well-defined
  `OpenVertexType` semantic, documented hole-free invariant,
  documented rotation convention, and a cross-checked geometric
  implementation.
- **`MatchTypeIndex`** has a coherent responsibility despite the
  undocumented surface — it's the central index for candidate
  enumeration, and its critical paths (`candidates_starting_at`,
  `new`, `get`) are heavily exercised and now cross-checked
  against brute force.

---

## Suggested order of attack

1. **C1** (rename second `GrowingPatch`) — minutes, isolates a name
   clash.
2. **C2** (tighten `growing.rs` visibilities) — small, large unused
   public surface gets dropped.
3. **H2** (delete dead tetrominos) — trivial.
4. **H1** (delete unused `MatchFinder` methods) — small and visible.
5. **H3** (tighten `patch.rs` internal visibilities) — surgical,
   guided by the dependency map.
6. **H4** (triage `GrowingPatch` public methods) — medium, requires
   decisions per item.
7. **M3, M4, M5, M6** (tighten visibilities of tests-only or
   externally-unused items) — bundled together, mechanical.
8. **M2** (`Rat` method audit) — small.
9. **C3 + L1, L2, L3, L4** (documentation pass) — larger
   investment; do as a dedicated PR after the visibility cleanup
   so the doc target is the actual public API, not the noise.
10. **M1** (Seed/Growing type split) — only if/when the state
    machine grows a third variant.

---

## Status table

| ID  | Severity | Type             | Status   | Notes |
|-----|----------|------------------|----------|-------|
| C1  | Critical | naming           | done     | renamed type to `RedelmeierPatch`; renamed module file `growing.rs` → `redelmeier.rs` |
| C2  | Critical | visibility       | done     | `redelmeier.rs`: `RedelmeierPatch` and `GrowStats` → `pub(crate)`; deleted unused `RedelmeierPatch::{len, is_empty}` and `grow_redelmeier_profiled`; `make_free` is now `#[cfg(test)]` private. `PatchPos`/`HasPatchPos`/`Pos2`/`Pos4`/`Pos8` kept `pub` (required by the bound on the public `grow_redelmeier*` entry functions) |
| C3  | Critical | docs             | pending  | doc-pass across patch/tileset/matchtypes/growing/rat |
| H1  | High     | dead code        | pending  | delete ~10 unused `MatchFinder` methods |
| H2  | High     | dead code        | pending  | delete 6 unused tetromino constructors |
| H3  | High     | visibility       | pending  | tighten internal `pub`/`pub(crate)` in `patch.rs` |
| H4  | High     | visibility       | pending  | triage `GrowingPatch` `pub` methods (3 buckets) |
| M1  | Medium   | architecture     | DEFERRED | split Seed/Growing types — wait for natural trigger |
| M2  | Medium   | visibility       | pending  | downgrade 9 unused `Rat` methods |
| M3  | Medium   | visibility       | pending  | `TileSet::index_of` → `pub(crate)`/`#[cfg(test)]` |
| M4  | Medium   | visibility       | pending  | `PatchVertexInfo` → `pub(crate)` |
| M5  | Medium   | visibility       | pending  | `TileSegment` → `pub(crate)` |
| M6  | Medium   | visibility       | pending  | tighten 5 `MatchTypeIndex` + `MatchType::apply` methods |
| L1  | Low      | docs             | pending  | doc-comment on `MatchType` field semantics |
| L2  | Low      | docs             | pending  | doc `MatchFinder::new` vs `crossing` |
| L3  | Low      | docs             | pending  | document `cyclotomic_intersect` feature flag |
| L4  | Low      | visibility       | done     | done as part of C2 |
