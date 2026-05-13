# Review: `src/intgeom/patch.rs`

A 3,188-line file mixing the `GrowingPatch` state machine, parallel raw-boundary
helpers, witness construction, candidate computation, and tests. The core
algorithm is sound, but the file has substantial duplication, a leaky state
machine, and several lazy or vacuous tests.

Findings are prioritized; resolution status is tracked at the bottom.

---

## High-severity

### H1. `PatchState::_Phantom` is dead and should be deleted (line 104)
`T` already appears in `match_index: Arc<MatchTypeIndex<T>>` (line 80), so no
`PhantomData` is needed. The variant is unreachable, yet every `match` over
`PatchState` must handle it (`_ => return …`, `_Phantom(_) => …`, `to_rat`
panics on it). Removing it eliminates ~12 dead branches and one `panic!`
(line 988).

### H2. Dead public API
- `AddTileDiff` (line 76): unit struct returned by `add_tile`; never inspected
  by any caller. Either remove or have it actually carry the diff.
- `Transition` (line 68): defined but unused. Real types are `NtTransition` /
  `ParsedTransition` elsewhere. Only `TransitionSide` is in use.
- `candidates_from_flat` (line 1739): no callers in workspace.
- `RawBoundary::normalize` (line 213): `#[allow(dead_code)]`, no callers.
- `junction_vertex_type_from` (free fn, line 129): superseded by
  `GrowingPatch::junction_vertex_type_at`; no other callers.

### H3. `add_tile_growing` is 220 lines with three manual rollbacks
The "destructure into locals with `mem::take`, then on failure rebuild via
`restore_growing`" pattern (called three times: 1383, 1401, 1476) is exactly
the failure mode that motivates a two-phase write:

```text
fn plan_glue(&self, pm) -> Option<NewGrowingState>   // pure, no mutation
fn commit_glue(&mut self, ns: NewGrowingState)       // infallible swap
```

That deletes `restore_growing`, the `mem::take` choreography, and the partial
grid mutation in the middle that must manually re-register segments on
failure.

### H4. `compute_glue_angles` enforced in the live path but not in the raw path
- `add_tile` / `init_from_first_add` route through `compute_glue_angles` (line
  1761), which rejects ±hturn glues.
- `glue_match_to_raw_boundary` (line 363) calls `angles::glue_raw_angles`
  directly, bypassing that check.

Witness construction can therefore accept configurations the live patch
rejects. Either both paths reject, or the asymmetry is intentional and must
be documented. No test exercises this.

### H5. `compute_candidates_covering_position_is_subset_of_all_candidates` is half a test
Only checks one direction (`covering ⊆ touching`). Passes vacuously if
`covering` is empty. The reverse direction — every match that genuinely
touches `target` is returned by `compute_candidates_covering_position` — is
the load-bearing invariant and it is not tested.

### H6. Failed `add_tile` has no "state-unchanged" test
`restore_growing` is exactly the kind of code that silently corrupts state
when someone adds a field and forgets one branch. No test snapshots the patch
state, attempts a failing `add_tile`, and asserts byte-equality afterward.

### H7. No self-intersection rejection test
`check_segment_clear` is the geometric safety net. No test constructs a patch
where the next tile would self-intersect and asserts `add_tile` returns
`None`. Without that, refactoring the grid logic could pass all current tests
while admitting overlapping tiles.

---

## Medium-severity

### M1. Duplicated junction-test predicate (5 copies)
`tileset.rat(edge.tile_id).seq()[edge.tile_offset] != angles[pos]` appears at
lines 138, 247, 1187, 1652, 1676. Extract into one helper.

### M2. Triple-duplicated edge-construction loop
The "copy surviving old edges, push first new edge, push remaining new edges"
pattern appears in:
- `glue_match_to_raw_boundary` (331–348)
- `init_from_first_add` (1291–1307)
- `add_tile_growing` (1517–1534)

### M3. `compute_candidates_covering_position` and `compute_all_candidates` share their junction block verbatim
Lines 702–735 and 777–808 are the same nested loop. Extract.

### M4. `get_matches_touching_vertex` has two near-identical branches
`Some(cbs)` / `None` differ only in source of `[Vec<PatchMatch>]`. Compute the
source first, then run one shared loop.

### M5. `flower_petal_glue` clones-and-reassigns RawBoundary in the loop
Destructures, builds a fresh `RawBoundary` inside the call, reassigns the
fields. Callers also clone all four vecs. Take `&mut RawBoundary` or shadow
per iteration.

### M6. The FIXME at line 403 is a design defect, not a TODO
The rotation convention forces every caller tracking positions to remap them
after each glue. Either fix or track in an issue.

### M7. `update_inner_chains` has an unused parameter (`_new_ptids`)
Delete or use it.

### M8. Brittle test fixtures
`t_tetromino` (3108) and `five_hex_cross` (2843) select matches by ordering
and structural predicates. If `get_all_matches` ordering ever changes, they
build a different patch silently.

### M9. Test family duplication
- `edges_self_consistent` already covers hex *and* square. `edges_bi_hex`,
  `edges_bi_square`, `edges_multi_hex` are strict subsets.
- `construct_minimal_witness_{hex,square,spectre,hex_with_inner}_roundtrip`,
  `construct_minimal_witness_{hex,square}_boundary_matches_brute_force`,
  `test_junction_angle_sequence_{hex,square}`: parameterize over fixture.
- `brute_force_zz4` / `brute_force_zz12`: make generic.

### M10. `seg_data` is append-only and indexed by global ID
`unregister_segment` clears the grid entry but `seg_data` keeps the stale
tuple forever (`seg_data[id]` lookup at 1825). Unbounded memory growth over
long lifetimes.

### M11. `dir_of_edge` panics on a non-unit edge
Should return `Option<i8>` or `debug_assert!`.

### M12. Tests that pass for the wrong reasons
- `glue_raw_angles_matches_rat_glue` (2741) sorts both sequences before
  comparison. Random permutation passes. Should compare canonical form.
- `construct_minimal_witness_*_boundary_matches_brute_force` (2617): same.

### M13. Vacuous tests
- `segment_cyclic_invariant` (1909): nested `if let` / `if .is_some()` guards
  can no-op away the second assertion.
- `mixed_hex_square_add_tile` (2027): only tests first match.
- `junction_vertex_ids_nonempty_after_each_add` (1929): iterates over
  matches computed against the seed, but mutates `gp` each loop — after the
  first iteration the `pm`s are stale.

### M14. Untested public/crate API
- `neighbor_junction_offsets` (1190)
- `vertex_type_at` (1132) — only used indirectly via `junction_vertices`
- `compute_segments` / `compute_junctions` — no unit tests on hand-crafted
  boundaries

---

## Low-severity

### L1. Two clone semantics
`GrowingPatch` derives `Clone`. `clone_for_mutation` differs only by clearing
caches. Name conveys mutability hint, not cache reset.

### L2. `from_parts` silently loses caller state
- `boundary_edge_ids` reset to `(0..n)` regardless of input.
- `grid` / `seg_data` rebuilt from scratch.

### L3. `cyclic_range_contains` (1749) has inclusive-on-both-ends semantics
Returns `true` for `index == start + len`. Useful for "match touches vertex"
semantics, but name suggests range-containment.

### L4. Missing doc comments on load-bearing helpers
- `junction_angle_sequence` (380): core of witness construction.
- `forward_match_length` (144): reversed indexing is opaque.
- `compute_glue_angles` (1761): ±hturn rejection invisible at call site.

### L5. `Seed::cached_matches` is wasted by `clone_for_mutation`
Computed eagerly on `new()`, discarded on `clone_for_mutation`. Make lazy or
keep through clones.

### L6. Add-tile growing has logical blocks worth extracting
- "trace boundary positions" (1442–1454)
- "segment clear" (1456–1469)

### L7. `Snake::try_from` validity check is in `init_from_first_add` (1281)
but not `add_tile_growing`. Resolve the asymmetry.

---

## Execution plan & status

The plan is to apply the cheap, mechanical fixes first, then revisit the
larger refactors after the file is shorter and the cruft is gone.

| ID  | Type       | Status   | Notes |
|-----|------------|----------|-------|
| H1  | dead code  | done     | `_Phantom` variant deleted |
| H2  | dead code  | done     | `AddTileDiff`/`Transition`/`candidates_from_flat`/`RawBoundary::normalize`/`junction_vertex_type_from` removed; `add_tile` now returns `bool` |
| M1  | duplicate  | done     | new `is_junction_at` helper; 4 inline copies collapsed |
| M7  | clean-up   | done     | dropped `_new_ptids` from `update_inner_chains` |
| M3  | duplicate  | done     | new `enumerate_junction_candidates_at` helper takes `keep`/`emit` |
| M4  | duplicate  | done     | `get_matches_touching_vertex` picks source up front, one shared loop |
| M9  | test churn | done     | generic `brute_force_patches`; shared helpers `assert_minimal_witness_roundtrips_for`, `assert_witness_matches_brute_force`, `assert_junction_angle_sequence_valid`; deleted 3 strict-subset edges tests |
| H5  | test       | done     | `compute_candidates_covering_position_matches_full_enumeration` asserts multiset equality with full enumeration (was: subset-only) |
| H6  | test       | done     | `add_tile_failure_leaves_state_unchanged` snapshots state+candidate classification, exercises early-bound and geometric rollback paths |
| H7  | test       | done     | `add_tile_rejects_geometrically_invalid_candidate` pins that some spectre candidates trigger the `check_segment_clear` rejection |
| H3  | refactor   | done     | hoisted paths 1, 2, 4 above `mem::take` as `debug_assert!` + release-mode `return false`; path 3 (geometric) remains the only `restore_growing` caller; latent path-4 bug fixed for free. Updated 2 tests that fed adversarial/stale pms |
| H4  | semantics  | done     | renamed `VertexType` → `OpenVertexType` (4 files, propagated to `OpenVertexTypeInfo`/`OpenVertexTypeIndex`); added documentation distinguishing open (boundary) vs closed (fully-surrounded) VTs; enforced the open-VT invariant in `construct_witness_from_vt_sequence_inner` — each glue step rejects with `None` if it produces a ±hturn boundary angle at the tracked junction. Low-level `glue_match_to_raw_boundary` stays permissive (allows future closed-VT demonstration constructors). |
| M2  | refactor   | done     | new `build_glued_edges` helper replaces three duplicated loops; `init_from_first_add` synthesizes the seed's old-edge list to share the same helper |
| X-check tests | new tests | done | four spectre-fixture tests: `add_tile_decision_agrees_with_snake_on_spectre` (incremental check vs Snake batch validator), `growing_patch_boundary_validates_as_snake_through_growth` (boundary is a valid Snake after each glue), `get_all_matches_matches_brute_force_on_spectre` (vs independent `(tile_b, ib, start_a)` enumeration), `get_matches_touching_vertex_matches_brute_force_on_spectre` |
| M11 | API change | done | `dir_of_edge` `panic!` → `unreachable!` with diagnostic; docstring on the upstream invariant |
| L3  | smell      | done | `cyclic_range_contains` doc explains inclusive-on-both-ends "vertex touched by match" semantics |
| L4  | smell      | done | doc comments on `forward_match_length`, `junction_angle_sequence`, `compute_glue_angles` |
| L7  | smell      | done | comment on `init_from_first_add`'s Snake-batch path explains why it differs from `add_tile_growing`'s incremental path (and points to the cross-check test) |
| M13 | test       | done | deleted `segment_cyclic_invariant` and `mixed_hex_square_add_tile` (vacuous subsets of `edges_self_consistent` / `edges_mixed_consistency`) |
| M14 | new tests  | done | `vertex_type_at_returns_consistent_info`, `neighbor_junction_offsets_returns_valid_offsets`, `tile_segments_partitions_boundary` (with explicit handling of the linear-vs-cyclic segment-start seam) |
| docs/refactor | refactor | done | `compute_segments` now uses the **canonical** `is_junction_at` check as the segment-break condition (the prior `tile_change` heuristic over-segmented on intra-tile wrap-arounds and could under-segment for pathological same-tile-id glues); `GrowingPatch` + `RawBoundary` doc the hole-free / edge-to-edge invariants; `TileSegment` docs the cyclic-vs-linear seam |
| L5 | latent bug | done | `clone_for_mutation` on a `Seed` previously discarded `cached_matches` (so `get_all_matches()` on the clone returned `[]`); fixed to preserve. Regression test `clone_for_mutation_preserves_seed_matches` |
| L2 | smell | done | `from_parts` docstring spells out the derived-state rebuild (grid, seg_data, boundary_edge_ids) so callers don't expect a faithful snapshot restore |
| L6 | smell | done | extracted `trace_polyline_from` (generic start + initial_dir polyline tracer; `trace_boundary_positions` now delegates) and `segments_all_clear` (multi-segment collision check with allowed-last-endpoint exception). `add_tile_growing`'s body is correspondingly shorter and more readable |
| M5  | refactor   | DEFERRED | RawBoundary ownership, after M3 |
| M6  | design     | DEFERRED | rotation convention — own ticket |
| M10 | refactor   | DEFERRED | seg_data GC — needs design |
| M11 | API change | DEFERRED | `dir_of_edge` signature |
| M14 | new tests  | DEFERRED | after refactors so we don't pin internals we'll move |
| M12 | test       | DEFERRED | strengthen comparisons after canonicalisation API exists |
| M13 | test       | DEFERRED | strengthen after de-duplication in M9 |
| M8  | test       | DEFERRED | re-spell fixtures |
| L1-L7 | smells   | DEFERRED | bundle into future polish PR |
