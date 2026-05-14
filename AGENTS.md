# Agent Instructions

REVIEW YOUR OWN WORK BEFORE SHIPPING. The most common agent failure is
producing work the agent itself would criticize if asked, because the
self-review step was skipped. Before you say "done", re-read your diff
with the question "if asked to review this, what would I flag?" If you
would flag something, fix it now. Do not make the user discover
problems you could have caught. Apply the five sins below to your own
diff, not just to other people's code.

TALK WITH THE USER. Don't run through walls. Some questions cannot be
answered from the code: what the user wants done, why a past choice
was made, whether surprising behavior is intentional, which tradeoff
to pick. Don't hypothesize about intent. Stop, summarize what you
have, and ask. Same for destructive actions (delete, force-push,
drop, overwrite) that weren't explicitly authorized.

---

## The five sins

These produce most cleanup work. Look for each in your own diff before
shipping.

**1. Coincidentally-correct logic.** A predicate that works on the one
test input you tried but breaks under any variation. State in one
sentence what specific variation would break it. If you can't, you
don't understand it. Test on structurally different inputs.

**2. Circular tests.** A test that compares `f(x)` to `g(x)` where
both sides go through the same helper; any bug in the helper makes
both sides agree. The reference must be a different implementation:
brute force, handwritten table, anything but the function under test.

**3. Undocumented preconditions on public functions.** A `pub` function
that trusts its caller for something it doesn't say it trusts. If
yours trusts the caller for anything (validity, ordering, in-range,
satisfies-some-invariant), say so in the doc and add a `debug_assert!`
to enforce it in test builds.

**4. Dead scaffolding.** Code, variants, fields, helpers left from
earlier iterations. If your PR adds something nothing else uses,
delete it before committing. Never use `#[allow(dead_code)]` on
committed code. If you'll genuinely need it later, add it later.

**5. Code duplication.** Three near-identical loops, four copies of
the same predicate. Two: fine. Three: refactor. Before writing a new
helper, check whether one exists.

---

## Working empirically

- **Use data, not theory.** Suspect a bug? Write a 5-line diagnostic
  (print, one-off test, assertion) before the fix. Print actual
  values. If the diagnostic contradicts your hypothesis, update your
  model. Don't argue with the evidence.
- **Time-box and ask.** A handful of tool calls on the same approach
  with no progress means your model is wrong, not that you need
  another try. Stop. Tell the user what you found.
- **Bisect to the smallest reproducer.** "The build is broken" is
  harder to fix than "this specific input returns this specific
  wrong value."

---

## Code quality

- **Visibility.** Default to private or module-scoped. `pub` is a
  commitment to external callers; tightening it later requires
  breaking changes.
- **Naming.** If you find yourself qualifying a name in its doc ("this
  is only for boundary vertices"), the name is wrong.
- **Comments.** WHY, not WHAT. The convention being assumed, the
  invariant being relied on - not what the code obviously does.

---

# Project-specific (tilezz)

## Toolchain

`cargo fmt && cargo clippy && cargo test --lib` before committing.
Some tests are slow (~5 min); ask before running the full `cargo
test` repeatedly. Brute completeness tests are `#[ignore]`-gated and
need `cargo test --release ... -- --ignored`.

`PATCH_REVIEW.md` tracks open review items (V1, V2, ...) on
`patch.rs` / `vertextypes.rs`. Reference these in commit messages.

## Conventions

- Tiles and patch boundaries are **angle sequences** (turn angles at
  vertices), not edge directions. Index conventions are subtle;
  document them at every function boundary that uses them.
- **Junctions** are positions where the boundary angle differs from
  the underlying tile's interior angle (two tile instances meet).
- **Match anchors** use `(start_a on self, start_b on other)` where
  `start_b` is the END index on the other tile (matching runs in
  opposite directions). Shared by `rat::get_match`,
  `PatchMatch::start_b`, `glue_raw_angles`. Misread it and you get
  silently nonsensical glues.
- **Rotation after a glue**: new boundary's index 0 is the first
  surviving old edge immediately CCW of the match. See
  `glue_match_to_raw_boundary`'s doc.
- **`pub` surface is small by design.** Most helpers are `pub(crate)`
  or private. Don't promote without a stated reason.

## Receipts

Each was latent for a meaningful period. Each is a sin above.

- `7a95704` (sins 1+2): `cyclic_range_contains` wrap bug at `start +
  len == n`. Only cross-check filtered both sides with the buggy
  function. Missed 51 spectre VTs / 92 transitions.
- `c08a762` V2 (sin 1): `covers_both` used vertex-touch instead of
  edge-coverage. Worked on hex (3 x 120deg = 360deg); would fail on
  square (3 x 90deg = 270deg).
- `7a95704` V2d (sin 3): `glue_raw_angles` was `pub` and trusted
  callers blindly; a hand-constructed `PatchMatch` could produce a
  geometrically nonsensical glue that passed self-intersection.
- `06fafba`, `2576734`, `f18c0d5` (sin 4): rounds of removing Phantom
  variant, AddTileDiff, PetalOutcome, ExploreResult.
- `890e31d`, `f5d5157` (sin 5): junction predicate written 4x across
  files; three near-identical glue-edge loops; copy-pasted test
  fixtures.
- `be5d605` (naming): `VertexType` was ambiguous (interior or
  boundary? open or closed?), renamed to `OpenVertexType` with
  invariant enforcement.
