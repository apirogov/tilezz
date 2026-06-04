# OEIS contribution map for cyclotomic matchstick polygons

What sequences exist, which we **match**, which we can **improve**
(extend or correct), and which we can **add** — plus how far each is
reachable on this commodity machine.

A "matchstick polygon" = a closed self-avoiding polygonal walk with
unit-length sides on a cyclotomic lattice `Z[zeta_n]` (`ZZn`). The
reference family is **Hugo Pfoertner's 2018 sequences** (Al Zimmermann
"Snakes on a Plane" contest); each encodes the **free** count (up to
rotation *and* reflection). His enumerator is unpublished.

## The core idea: one free run -> a shelf of sequences

A single **free** enumeration of `ZZn` to perimeter `N` is a superset.
Every related sequence is a *filter / bucket* over its rats, computed
one rat at a time — no re-enumeration, no extra machinery. `rat_enum`
emits all of them into each dataset's `ro-crate-metadata.json` as
schema.org `variableMeasured`, and `tools/verify_counts.py` re-derives
and checks them:

| derived sequence | filter on the free rats |
|---|---|
| **free** | all rats (the base count) |
| **oneSided** | up to rotation only = `2*free - achiral` |
| **achiral** | equal to its own mirror image |
| **rotationSymmetric** | repetition factor > 1 |
| **symmetric** | achiral OR rotationSymmetric |
| **subring** | all turns even = the order-`n/2` sub-ring |
| **coset** | all turns odd (no straight segments; even perimeter) |

This is the maximum-effect/least-effort lever: **run free once per
ring as deep as the budget allows; everything else falls out.** In
particular the odd-turn ("bipartite", A316195/197/199) families need
no turn-restriction mechanism — they are the `coset` filter.

## Match (verified, reproduced exactly)

Free counts, multiple independent paths (prunes on/off, streaming,
subring-step cross-checks). Pinned in `opt_correctness_tests`.

| Ring | OEIS | verified terms |
|---|---|---|
| ZZ4 free (square, perim 2n) | A266549 | a(2..7) |
| ZZ6 free (triangular) | A284869 | a(3..14), 12 terms (deepest oracle) |
| ZZ8 free (perim 2n) | A316198 | a(2..6) |
| ZZ10 free | A316200 | a(4..10) *(a(11): OEIS overcounts; correction in Improve)* |
| ZZ12 free | A316192 | a(3..10) |

Now also matchable **by filtering** (no new code):

| OEIS | what | status |
|---|---|---|
| A316194 | ZZ4 symmetric | **verified exact thru perim 16** (= `symmetric` filter) |
| A316196 | ZZ6 symmetric | matchable (`symmetric` filter) |
| A316195 | ZZ10 coset (odd turns) | **verified** perim 6/8/10 = 2,1,18 |
| A316197 | ZZ14 coset | matchable (`coset` filter) — not yet checked deep |
| A316199 | ZZ18 coset | matchable (`coset` filter) — not yet checked deep |
| A057730 / A057729 | ZZ4 / ZZ6 holes-allowed | cross-check up to first hole (n=8 / n=11) |

## Improve (extend or correct published data)

- **CORRECT A316200(11): 19405 -> 9883.** Resolved this session by an
  independent clean-room exact-`Z[sqrt5]` enumerator (zero floating
  point) that reproduces a(4..10) and gives a(11)=9883; tilezz agrees
  via five paths. Full evidence + reproducible enumerators in
  `docs/oeis-A316200-correction/`. **Submittable.**
- **Extend the free sequences** past their published depth (all
  reachable; see Reach):
  - A316192 (ZZ12): a(11..14) = 89075, 597581, 4076855, 28499301.
  - A316198 (ZZ8): a(7) = 240549 (perim 14).
  - A316200 (ZZ10): a(12..14) once a(11) is accepted.
  - A284869 (ZZ6): OEIS lists to ~a(24); a day of compute reaches it.
- **Extend the coset/symmetric** sequences (A316195/197/199,
  A316194/196) past their published terms — same filtered free runs.
  Given A316200 was wrong at its deepest term, **re-checking the
  deepest published coset terms (A316197 ZZ14, A316199 ZZ18) is the
  best shot at another correction.**

## Add (new registrations, all by filtering)

None of these are in OEIS; each is a `variableMeasured` column of a
free run:

- **one-sided** for every ring (`2*free - achiral`) — uncatalogued.
- **odd sub-rings** ZZ3/5/7/9 (= `subring` filter of ZZ6/10/14/18) —
  no OEIS entry; also obtainable via `--step 2` directly.
- **symmetric / achiral / rotationSymmetric** for rings beyond ZZ4/6;
  and the achiral/rotational split (A323188/A323189-style).
- **full-ring free** for ZZ14, ZZ16, ZZ18, ZZ20, ZZ24 (no OEIS oracle
  today) — straightforward submissions with a methodology citation.

## Reach (measured; commodity machine, all cores)

Extrapolated from a `--mode bench --free --threads 0 --mod-prune
--closure-key-prune` sweep (count + timing per `n`). "10 min" and "1
day" are wall-clock for the **enumeration**; counts are cumulative
(<= perimeter `n`).

| Ring | branch | 10 min: n (~rats) | 1 day: n (~rats) | notes |
|---|--:|---|---|---|
| ZZ4  | 3  | ~28 (1.4e7) | ~34 (4e9)  | square; even perim only; tiny DAFSA |
| ZZ6  | 5  | ~20 (4e7)   | ~24 (8e9)  | triangular; **reaches A284869's frontier** |
| ZZ8  | 7  | 16  (5e6)   | ~20 (3e9)  | even perim only |
| ZZ10 | 9  | ~15 (1.5e7) | ~18 (4e9)  | |
| ZZ12 | 11 | ~14 (3e7)   | **~16 (1.6e9)** | north star (OEIS A316192) |
| ZZ14 | 13 | 12 (8e5)    | ~14 (2e8)  | cubic-root sign, slow per node |
| ZZ16 | 15 | 12 (1.7e6)  | ~14        | even perim only |
| ZZ18 | 17 | ~11 (1e6)   | ~14 (8e8)  | cubic-root sign |
| ZZ20 | 19 | ~13 (2e7)   | ~15 (6e8)  | |
| ZZ24 | 23 | ~11 (4e6)   | ~13 (4e8)  | |
| ZZ32 | 31 | ~8–9        | ~10        | not benchmarked; i128 sign + huge branching |
| ZZ60 | 59 | ~6–7        | ~8         | extreme branching; smallest reach |

**Memory / disk (the real ceiling at the deep end).** The DAFSA
*build* holds the automaton in RAM, and that scales with **states**
(sublinear in rats), not rats. Reference (ZZ12): n=14 -> 33M rats /
613K states / 2.67M edges / 9 MB gz; n=16 -> ~1.7e9 rats but only
~7M states -> ~28 GB build RAM, ~130 MB gz blocks. So:

- Sequences/counts via `--mode bench`: time-limited only (cheap RAM).
- Full **datasets** (DAFSA): the billion-rat "1 day" reaches need a
  32–64 GB box for the build; gz hosting is ~100s of MB. Use the
  streaming pipeline (`stream`/`merge`/`build`) to bound RAM during
  the sort; the build stage is the bottleneck.
- Per-lookup cost in the explorer is unaffected (lazy block fetch).

## Workflow per ring (the lever, concretely)

1. Build `rat_enum` from a tagged (`origin/main`) commit (reproducibility).
2. One deep free run: `rat_enum --ring R -n N --free --threads 0
   --mod-prune --closure-key-prune --mode dafsa-blocks -o ...`
   (or the streaming pipeline for big N).
3. The RO-Crate now carries all 7 sequences (`variableMeasured`);
   `tools/verify_counts.py` confirms them; `tools/count_by_length.py`
   prints the OEIS-style terms.
4. Cross-validate via a `--step` subset against a pinned ring
   (ZZ24-step2=ZZ12, ZZ20-step2=ZZ10, ...); for a submission, also the
   overflow-checks verification gate (MAINTENANCE.md §5).
5. Submit: free / coset / symmetric / one-sided extensions and the
   no-oracle full rings; correct A316200(11).

## Out of scope / caveats

- **ZZ22** (A316201 coset): no supported ring and no subring
  cross-check -> correctness would rest on internal consistency only.
- **ZZ32 / ZZ60**: reachable only to small `n`; no external oracle
  except `--step` subset cross-checks; least-anchored rings.
- **Open-walk** siblings (A306175/177-182) and **cell-count**
  polyomino sequences (A001168/A000105/A000988, A258206) enumerate
  different objects (non-closing walks; area-indexed) — would need an
  open-walk mode or area post-processing, not on this lever.
