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
schema.org `variableMeasured`, and `tools/count.py --verify` re-derives
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

OEIS provenance (re-verified against oeis.org, 2026-06): A266549 (ZZ4
free) is Luca Petrone, to a(20). A284869 (ZZ6 free) is Petrone for
a(1..14), Hugo Pfoertner for a(15), and **Walter Trump (Nov 2023) for
a(16..22), a(22)=374128188** -- so the disputed a(22) is Trump's, not
Petrone's. Every A3161xx sequence -- free ZZ8/ZZ10/ZZ12,
symmetric, and coset -- is Hugo Pfoertner (Jun/Jul 2018). The
holes-allowed siblings A057729/A057730 are N. J. A. Sloane. All carry
the OEIS `more` keyword, i.e. they are open for extension. Indexing
differs: A266549, A316194, A316195/197/199/201, A316198 and A057730 use
a(n) = perimeter 2n (square / coset / even-perimeter families); A284869,
A316192, A316196, A316200 and A057729 index by perimeter n directly.
Quote terms in each sequence's own convention.

Now also matchable **by filtering** (no new code):

| OEIS | what | status |
|---|---|---|
| A316194 | ZZ4 symmetric | **verified exact thru perim 16** (= `symmetric` filter) |
| A316196 | ZZ6 symmetric | matchable (`symmetric` filter) |
| A316195 | ZZ10 coset (odd turns) | **verified** perim 6/8/10 = 2,1,18 |
| A316197 | ZZ14 coset | verified a(1..7) — see Verification log |
| A316199 | ZZ18 coset | verified a(1..6) — see Verification log |
| A057730 / A057729 | ZZ4 / ZZ6 holes-allowed | cross-check up to first hole (n=8 / n=11) |

## Verification log (full screening against published OEIS, 2026-06)

A one-time term-by-term screen of every published sequence we touch
against tilezz (`tools/count.py` for the per-perimeter counts; the
`variableMeasured` symmetric / coset arrays for the filtered ones).
**Every published term we can reach agrees exactly, with exactly two
exceptions** -- the two disputes below. "Reached frontier?" = did we
verify through the deepest published term.

| OEIS | object | verified | published frontier | reached frontier? |
|---|---|---|---|---|
| A316192 | ZZ12 free | a(1..10) | a(10) | YES (full extent) |
| A316198 | ZZ8 free | a(1..6) | a(6) | YES (full extent) |
| A316194 | ZZ4 symmetric | a(1..8) | a(8) | YES (full extent) |
| A316196 | ZZ6 symmetric | a(1..15) | a(15) | YES (full extent) |
| A316195 | ZZ10 coset | a(1..7) | a(7) | YES (full extent) |
| A316199 | ZZ18 coset | a(1..6) | a(6) | YES (full extent) |
| A316197 | ZZ14 coset | a(1..7) | a(7) | YES (full extent) |
| A316200 | ZZ10 free | a(1..10) | a(11) | a(1..10) match; **a(11) DISPUTED** (tilezz 9883 vs 19405; see Improve / docs/oeis-A316200-correction) |
| A284869 | ZZ6 free | a(1..21) | a(22) | a(1..21) match; **a(22) DISPUTED** (tilezz 374128154 vs 374128188; see docs/oeis-A284869-zz6-a22) |
| A266549 | ZZ4 free | a(1..16) | a(20) | **NO** -- see exception below |

**The ZZ4 exception.** A266549 is the one sequence published well past
our reach: polyominoes / square-lattice self-avoiding polygons are
extremely well studied, and a(20) (perimeter 40) comes from specialized
finite-lattice / transfer-matrix methods. We verified a(1..16)
(perimeter 32, the depth we have enumerated) and match exactly; a(17..20)
are simply beyond our cheap reach, not a disagreement. ZZ4 is therefore
not an extension target -- we are behind the published frontier there.

For every other sequence we screened, we reached the **last published
term** and matched it (the two flagged disputes aside). That is the
clean base for the extension submissions in Improve/Add: where we match
the full published extent, the new terms beyond it rest on an enumerator
already confirmed against every known value.

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
  - A284869 (ZZ6): OEIS lists to a(22)=374128188 (Petrone); a day of
    compute (n~24) extends it.
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
- **symmetric / achiral / rotationSymmetric** for rings beyond ZZ4/6,
  and the achiral vs rotational split per family. (A323188/A323189 do an
  analogous mirror/point-symmetry split but for square-lattice
  self-avoiding *walks*, not closed polygons -- ours are the
  closed-polygon analogue, uncatalogued.)
- **full-ring free** for ZZ14, ZZ16, ZZ18, ZZ20, ZZ24 (no OEIS oracle
  today) — straightforward submissions with a methodology citation.

## Reach (measured 2026-06-05; 16-core / 61 GB box, all cores, full prunes)

Per-ring `--mode bench --free --threads 0 --reachability-prune
--closure-table-prune` timing sweeps: a measured anchor `(n, seconds)`
plus the measured per-`n` growth factor, extrapolated PESSIMISTICALLY
(factor held flat or worst-observed; reach rounded down). These are
ENUMERATION (count) times; a deployable dataset adds the build stage
(see Resources). Counts cumulative (<= perimeter `n`).

| ring | step | growth /+1 n | anchor | 10 min | ~1 hr | ~1 day | ~rats @1-day |
|---|---|---|---|---|---|---|---|
| ZZ4  | even  | x5.7 /+2     | n27=36s  | n28 | n30    | n34 (~9h)    | ~2e9 |
| ZZ6  | every | x3.65        | n18=42s  | n20 | n21    | n23 (~8h)    | ~1.9e9 |
| ZZ8  | even  | x22.9 /+2    | n16=83s  | n16 | n18    | n20 (~12h)   | ~3e9 |
| ZZ10 | every | x6 (rising)  | n14=41s  | n15 | n16    | n17-18       | 5e8 / 2.4e9 |
| ZZ12 | every | x7.0         | n13=85s  | n14 | n14-15 | n16 (~8h)    | 1.6e9 |
| ZZ14 | every | x16.8        | n12=53s  | n12 | n13    | n13-14 (~4h) | ~2e8 |
| ZZ16 | even  | >=x26 /+2    | n12=42s  | n12 | n12-14 | ~n14         | ~1e8 |
| ZZ18 | every | x9.0         | n10=14s  | n11 | n12    | n13 (~3h)    | ~1e8 |
| ZZ20 | every | x14.7        | n12=255s | n12 | n12-13 | n14 (~15h)   | ~1.7e9 |
| ZZ24 | every | x11.3        | n10=25s  | n11 | n12    | n13 (~10h)   | ~6e8 |

Even-only rings (ZZ4/8/16) gain new counts only on even `n`. ZZ16 is
the least-pinned: n14 alone takes >18 min, so its 1-day cell is a
rough placeholder (a dedicated overnight calibration is needed to
target it). ZZ6 1-day is n23 (~8h), NOT n24 (n24 ~28h). On the
rising/steep rings (ZZ10/14/20/24) treat the 1-day cell as an upper
bound -- bank one less `n` for certainty. ZZ32/ZZ60 not calibrated
(extreme branching; small-n only).

## Resources (measured 2026-06-05; constants from a stream->merge->build sweep, rings 4/6/8/10/12)

- **Streaming-build RAM ~= 20 MB + ~280 B/state** -> even ZZ12 n16
  (~7M states) ~= **~2 GB**. RAM is NOT the constraint. (A prior
  "~28 GB" estimate was wrong/conflated with the in-memory path.)
- **Final dataset (gz blocks) ~= ~5 B/edge** -> **~100-200 MB even
  for billion-rat sets**. Hosting is a non-issue.
- **Peak SCRATCH disk** -- the streaming `runs/` + `unique.bin` hold
  ~dihedral-order (~2n) duplicate copies of each canonical rat (the
  "no HashSet" design trades RAM for disk; the DFS canonical prune
  cannot fully dedup the dihedral images). These duplicate records are
  identical canonical bytes written adjacently in the sorted runs, so
  **the scratch is now gzipped** (`Compression::fast()` at all four
  scratch I/O points; final dataset unchanged, already gz). Measured
  ZZ12 n13: `runs/` 2069->345 MB (6.0x), `unique.bin` 2134->19 MB
  (112x), **peak 4203->364 MB (11.5x)**. Uncompressed this was ~450 B
  per free rat (ZZ12 n16 ~690 GB, larger than a 341 GB volume);
  gzipped, **ZZ12 n16 peak ~= ~60 GB (worst case ~115 GB) -- fits a
  341 GB volume with headroom.** The certificate hashes record content
  (not file bytes), so it is unaffected.

The in-memory `--mode dafsa-blocks` path instead holds all rats
(~rats x 50 B): fine to ~1e8 rats (<~5 GB RAM), but the billion-rat
reaches would need ~50-80 GB RAM -- use the streaming pipeline there.
Per-lookup cost in the explorer is unaffected (lazy block fetch).

## Workflow per ring (the lever, concretely)

1. Build `rat_enum` from a tagged (`origin/main`) commit (reproducibility).
2. One deep free run: `rat_enum --ring R -n N --free --threads 0
   --reachability-prune --closure-table-prune --mode dafsa-blocks -o ...`
   (or the streaming pipeline for big N).
3. The RO-Crate now carries all 7 sequences (`variableMeasured`);
   `tools/count.py --verify` confirms them; `tools/count.py --print`
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
