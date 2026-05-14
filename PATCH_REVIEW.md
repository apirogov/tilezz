# Review: vertex type classification (`vertextypes.rs`)

Scope: `src/intgeom/vertextypes.rs` — `OpenVertexTypeIndex`, the BFS
catalog of open vertex types reachable from a tileset's seed tiles,
and the `VTypeKind` (Dead / Undead / Blessed / Free) classification.

Findings are organized by severity. Each item includes context, the
concrete impact, and a recommended fix.

---

## Critical

### V1. Dead code block in BFS

Lines 286–300 of the BFS loop body:

```rust
{
    let (gp, _, _) = &witness_store[&vt];
    let all_matches: Vec<PatchMatch> = gp.get_all_matches();
    for pm in &all_matches {
        if cyclic_range_contains(pm.start_a, pm.len, pos, n)
            || cyclic_range_contains(pm.start_a, pm.len, (pos + n - 1) % n, n)
        {
            continue;
        }
        let mut gp2 = gp.clone();
        if !gp2.add_tile(pm) || !gp2.is_growing() {
            continue;
        }
    }
}
```

Iterates over every match that does **not** touch the focus vertex,
clones the patch, applies the glue, then does nothing observable —
no transitions recorded, no VTs added to the catalog.

Three possibilities, all bad:

- The author intended a sibling discovery loop (the touching-match
  loop above it at lines 271–283 does discover new VTs from
  successful glues) and forgot the body.
- It's a leftover from a refactor.
- It's a stale diagnostic.

**Recommended fix:** investigate via `git blame` what was originally
there. Either remove it or restore the intended behaviour.

### V2. Side / junction_pos derivation lacks correctness comments

Lines 236–257 contain dense modular arithmetic deriving the post-glue
junction position and the side (CW/CCW) of the transition:

```rust
let covers_ccw = cyclic_range_contains(pm.start_a, pm.len, pos, n);
let covers_cw = cyclic_range_contains(pm.start_a, pm.len, (pos + n - 1) % n, n);

let junction_pos = if pm.start_a == pos { n - pm.len } else { 0 };

let covers_both = covers_cw && covers_ccw;
let side = if covers_both || covers_cw {
    TransitionSide::Cw
} else {
    TransitionSide::Ccw
};

let edge_pos = if side == TransitionSide::Cw {
    (pos + n - 1) % n
} else {
    pos
};
let offset_in_match =
    (edge_pos as i64 - pm.start_a as i64).rem_euclid(n as i64) as usize;
let m = tileset.rat(pm.tile_id).len();
let tile_offset = (pm.start_b as i64 + pm.len as i64 - offset_in_match as i64)
    .rem_euclid(m as i64) as usize;
```

With patch.rs's rotation convention in mind, the math is correct, but
extremely hard to verify without re-deriving from scratch. The naming
also misleads: `covers_ccw` doesn't mean "covers the CCW side" — it
means "the match covers the focus vertex `pos`". Likewise
`covers_cw` is "covers vertex `pos - 1`".

There's a latent invariant being relied on: `side` is set to `Ccw` in
the fallthrough case where **neither** `covers_cw` nor `covers_ccw` is
true. That branch is currently unreachable because `touching =
get_matches_touching_vertex(pos)` filters to matches that touch
`pos`, but the code doesn't assert it — a future change to
`get_matches_touching_vertex` semantics could silently produce
spurious transitions with `side = Ccw, junction_pos = 0`.

**Recommended fix:** add a doc-comment explaining the rotation-
convention math, rename `covers_cw`/`covers_ccw`
(`match_covers_pos_minus_1` / `match_covers_pos`), and add
`debug_assert!(covers_cw || covers_ccw)` to pin the invariant.

---

## High

### V3. `compute_blessed` is `O(T·N²)` when it could be `O(T·N)`

```rust
let all_closed_or_blessed = transitions
    .iter()
    .filter(|t| t.src_id == id)
    .all(|t| t.is_closed() || blessed.get(&t.dst_id).copied().unwrap_or(false));
```

For every VT in every fixpoint iteration, this scans **all**
transitions to filter by `src_id`. With `T` transitions and `N` VTs
and up to `N` fixpoint iterations, total complexity is `O(T·N²)`.

`compute_cursed` already pre-buckets successors via `succ_sets`.
`compute_blessed` should do the same.

**Recommended fix:** pre-bucket transitions by `src_id` once (e.g.
`Vec<Vec<&TransitionInfo>>` indexed by `src_id - 1`), then
`compute_blessed` iterates over each VT's transitions directly.
Drops to `O(T·N)`.

For tilesets with thousands of VTs and tens of thousands of
transitions (spectre, mixed tilesets), the difference is significant.

### V4. `is_cursed` / `is_blessed` use `HashMap<usize, bool>` keyed by id

```rust
let is_cursed = compute_cursed(&entries, &transition_map, &succ_sets);
let is_blessed = compute_blessed(&entries, &transition_infos);
```

Both return `HashMap<usize, bool>` where the key is the 1-based id.
Every classification lookup is a hash; the value is just a bool.

**Recommended fix:** return `Vec<bool>` indexed by `id - 1`. Simple
swap, cleaner code, faster.

### V5. `OpenVertexTypeIndex::new` is one 240-line function

The constructor performs five logically distinct phases:

1. Seed phase (lines 164–192): iterate all tile-pair first-glues,
   extract initial VTs.
2. BFS phase (lines 197–301): pop VTs, explore touching matches,
   accumulate transitions, discover new VTs.
3. Reverse-map / id assignment (lines 303–308).
4. Successor / predecessor sets + transition list (lines 310–341).
5. Classification + final info entries (lines 343–390).

Each is independently testable but currently inlined with shared
mutable locals. The function is hard to read end-to-end.

**Recommended fix:** decompose into private helpers per phase. Each
becomes ~50 lines, individually documentable. Logging (V8) becomes
easier to isolate too.

### V6. `VTypeKind` semantics are not documented anywhere

Dead / Undead / Blessed / Free are domain-specific terms, not
standard. Reading the code reveals:

- **Dead**: VT has zero outgoing transitions (no glue touching this
  vertex succeeds, or none exist).
- **Undead** (≈ "cursed"): VT is reachable but trapped — every
  successor, transitively via fixpoint, is dead/undead. Glues exist
  but lead nowhere productive.
- **Blessed**: every outgoing transition either closes the junction
  (`dst_id == CLOSED_ID`) or leads to another blessed VT (fixpoint).
  Geometrically: every continuation eventually reaches a fully-
  surrounded state.
- **Free**: alive (not cursed) and not blessed. Has at least one
  continuation that's "in play".

A reader has to derive these from the fixpoint code.

**Recommended fix:** doc-comment on the enum + on each variant.
Similar treatment for `VTypeKind::segment_kind`'s lattice operations
(cursed = OR; dead = OR; blessed = AND).

### V7. `OpenVertexTypeIndex::new` has no doc-comment

A 240-line constructor with no docstring. A short paragraph
explaining "BFS over all VTs reachable from any tile's first-glue
junctions, classifying each as Dead/Undead/Blessed/Free; transitions
between VTs are recorded with `tile_id` / `tile_offset` describing
the glue that effects the transition" would orient anyone touching
this code.

---

## Medium

### V8. `eprintln!` progress logs are unconditional

Lines 165–169, 204–210, 233 all use `eprintln!` to report BFS
progress. For library code this is invasive — stderr noise on every
`OpenVertexTypeIndex::new()` call. Tests that build a VT index also
print these.

**Recommended fix:** gate behind a feature flag, expose a
`trace_progress: bool` setting on the constructor, or wire through a
`log` / `tracing` macro that the user can route.

### V9. `OpenVertexTypeInfo` has 3 `#[allow(dead_code)]` fields

```rust
#[allow(dead_code)] successors: Vec<usize>,
#[allow(dead_code)] predecessors: Vec<usize>,
#[allow(dead_code)] has_transitions: HasTransitions,
```

Computed during construction (from `succ_sets` / `pred_sets` /
`transition_map`) and stored on the info struct, but never read
outside `new`. Free memory + maintenance cost.

**Recommended fix:** if they're for future readers, comment why.
If truly dead, drop them. If they're for debugging,
`#[cfg(debug_assertions)]` gate.

### V10. `transition_map` is mutated mid-BFS

```rust
if t.is_empty() {
    transition_map.insert(vt, HasTransitions::No);
    continue;
}
transition_map
    .entry(vt.clone())
    .or_insert(HasTransitions::Yes);
```

The `transition_map` is populated incrementally as BFS visits each
VT, then later consumed by `compute_cursed` as a source-of-truth.
Mixing "populate during enumeration" with "consume after" muddies the
BFS loop's contract.

**Recommended fix:** BFS just builds `raw_transitions`; afterwards, a
separate pass over the catalog derives `transition_map`. Cleaner
phase separation. Pairs with V5.

### V11. `entry.or_insert_with` vs raw `insert` inconsistency

Seed phase:
```rust
witness_store
    .entry(vt.clone())
    .or_insert_with(|| (gp.clone(), pos, gp.angles()[pos]));
```

BFS phase:
```rust
witness_store
    .insert(nv.clone(), (gp2.clone(), new_pos, gp2.angles()[new_pos]));
```

Seed phase preserves the first-seen witness; BFS phase overwrites.
Since `visited.insert(nv.clone())` gates BFS-phase insertions (only
on first visit), the difference doesn't matter in practice. But the
asymmetric semantics is a footgun if someone later loosens the BFS
guard.

**Recommended fix:** use `or_insert_with` consistently.

### V12. Memory: full `GrowingPatch` per VT in `witness_store`

`witness_store` holds an `(GrowingPatch, usize, i8)` for each
discovered VT. The `GrowingPatch` carries the spatial grid (`grid`,
`seg_data`, `boundary_edge_ids`, `positions`) — a few KB per entry.
For tilesets with thousands of VTs (e.g. spectre at depth 3+), this
adds up.

The witness's spatial grid isn't used by the catalog after BFS
terminates.

**Recommended fix:** either drop the grid after BFS (re-build via
`from_parts` on demand), or store just `RawBoundary` data and rebuild
the patch lazily when `OpenVertexTypeInfo::witness()` is called. Not
urgent for current workloads; worth noting.

---

## Low

### V13. Test coverage is weak

Four tests in `mod tests`: three smoke tests (counts > 0, with
`eprintln!` reporting) and one witness-roundtrip test. None of them
asserts the classification logic is correct.

For a non-trivial fixpoint algorithm, the right tests are:

- **Specific VT counts** for known tilesets (hex, square, spectre):
  `assert_eq!(idx.num_types(), KNOWN_VALUE)`.
- **Specific classification outcomes**: e.g. "in the hex tileset, no
  VT should be Dead", "in the square tileset, the 4-square-corner
  VT is Blessed".
- **Transition correctness**: for every transition `(src, dst)`,
  gluing the recorded `(tile_id, tile_offset)` to `src`'s witness
  should produce a patch realizing `dst` (or closing the junction).
- **Fixpoint correctness**: `is_cursed` should equal "no path leads
  to a Free or Blessed VT"; `is_blessed` should equal "every path
  closes or stays blessed".

The transition-correctness check in particular would cross-validate
the `side` / `junction_pos` arithmetic (V2).

### V14. `for (i, _vt) in entries.iter().enumerate()` in fixpoint loops

`_vt` is unused. `for i in 0..entries.len()` directly. Both
`compute_cursed` and `compute_blessed` do this.

### V15. Public methods on `OpenVertexTypeIndex` are undocumented

`num_types`, `transitions`, `entries`, `range_by_cw`, `get_id`,
`get_type`, `get_info`, `tileset` — none have docs. Same for
`OpenVertexTypeInfo`'s accessors, `TransitionInfo`'s fields,
`CLOSED_ID`.

`CLOSED_ID = 0` is a sentinel meaning "this transition closes the
junction (no destination VT)". Worth a one-line doc.

### V16. `range_by_cw(tile_id) -> &[OpenVertexTypeInfo]` has an undocumented invariant

Uses `partition_point` on the entries' `cw.tile_id`. This works
because entries are sorted by VT (`(cw, inner, ccw)` lex order), and
`cw.tile_id` is the most-significant sort key. A future refactor
that changes entries' sort order would silently break this method.

**Recommended fix:** doc-comment noting "entries are sorted with
`cw.tile_id` monotonically non-decreasing; this method relies on
that invariant."

---

## What's working well

- The 4-kind classification is a finer-grained reachability analysis
  than typical, and the fixpoint convergence is clean once
  understood.
- `succ_sets` / `pred_sets` are built and stored efficiently as
  `BTreeSet<usize>`.
- `reverse: HashMap<OpenVertexType, usize>` for id lookup is the
  right structure.
- The witness consistency test (`hexagon_witness_consistency`) is
  the right shape of check.
- Deterministic ordering (`BTreeSet<OpenVertexType>` → sorted
  `entries`) makes the index reproducible.

---

## Recommended order of attack

1. **V1** — investigate the dead block. Remove or restore.
2. **V6 + V7** — doc-comments on `VTypeKind` variants + on
   `OpenVertexTypeIndex::new`.
3. **V2** — comment the side / junction_pos derivation +
   `debug_assert`.
4. **V5** — split `new` into phases.
5. **V3 + V4** — efficiency improvements to `compute_blessed`,
   replace `HashMap<usize, bool>` with `Vec<bool>`.
6. **V13** — stronger classification tests (especially the
   transition-correctness cross-check).
7. Lower-priority: **V8** (logging), **V9** (dead fields), **V10**
   (transition_map phase), **V11** (insert consistency), **V12**
   (memory), **V14 / V15 / V16** (style + docs).

---

## Status table

| ID  | Severity | Type      | Status   | Notes |
|-----|----------|-----------|----------|-------|
| V1  | Critical | dead code | pending  | investigate + remove or restore |
| V2  | Critical | docs+invariant | pending | comment math, rename booleans, add `debug_assert!` |
| V3  | High     | perf      | pending  | pre-bucket transitions in `compute_blessed` |
| V4  | High     | perf      | pending  | `HashMap<usize, bool>` → `Vec<bool>` |
| V5  | High     | refactor  | pending  | split `new` into 5 phase helpers |
| V6  | High     | docs      | pending  | document `VTypeKind` variants + `segment_kind` |
| V7  | High     | docs      | pending  | document `OpenVertexTypeIndex::new` algorithm |
| V8  | Medium   | observability | pending | gate `eprintln!` progress logs |
| V9  | Medium   | smell     | pending  | resolve 3 `#[allow(dead_code)]` fields |
| V10 | Medium   | refactor  | pending  | derive `transition_map` post-BFS, not mid |
| V11 | Medium   | smell     | pending  | use `or_insert_with` consistently |
| V12 | Medium   | memory    | pending  | drop grid in `witness_store` after BFS |
| V13 | Low      | tests     | pending  | classification + transition-correctness tests |
| V14 | Low      | style     | pending  | drop unused `_vt` in fixpoint loops |
| V15 | Low      | docs      | pending  | doc-comment public methods + `CLOSED_ID` |
| V16 | Low      | docs      | pending  | document `range_by_cw` sort-order invariant |
