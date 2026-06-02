# OEIS cross-references for cyclotomic matchstick polygons

Tracks every OEIS sequence relevant to the polygon enumeration in
this crate, our verification status, and what we can contribute back.

The canonical reference for "matchstick polygons" — closed
self-avoiding polygonal walks with unit-length sides on a cyclotomic
lattice — is **Hugo Pfoertner's 2018 family of sequences**, all
submitted on 2018-07-07 around the [Al Zimmermann "Snakes on a Plane"
contest][snakes]. Hugo's source code is not public; the only artifact
is the OEIS data and illustrations at [randomwalk.de][rwd]. Every
sequence in the family encodes the **free** (also called **dihedral-canonical**) count:
distinct polygons up to rotation *and* reflection of the angle
sequence.

[snakes]: http://www.recmath.org/contest/Snakes/index.php
[rwd]: http://www.randomwalk.de/sequences/

## Verified matches

Each row's ✅ terms match exactly under multiple independent
enumeration paths (`--mode bench` with/without prunes,
`--mode stream`/`merge`/`build`, and where applicable a subring-step
cross-check like ZZ20 step=2 = ZZ10 step=1).

| Ring | OEIS | Description | OEIS terms (verified by us) | New terms we can contribute |
|------|------|-------------|------------------------------|-----------------------------|
| **ZZ4** free | [A266549] | square lattice, perim 2n | a(2..7) = 1, 1, 3, 6, 25, 86 ✅ | a(8+) trivially (ring is small) |
| **ZZ6** free (triangular) | [A284869] | perim n | a(3..14) = 1, 1, 1, 4, 5, 16, 37, 120, 344, 1175, 3807, 13224 ✅ | OEIS publishes a(3..24); deeper terms via further enumeration |
| **ZZ8** free | [A316198] | perim 2n | a(2..6) = 2, 6, 59, 695, 12198 ✅ | **a(7) = 240549** (new, perim 14) |
| **ZZ10** free | [A316200] | perim n | a(4..10) = 2, 2, 10, 15, 124, 352, 2378 ✅ | **a(12..14)** once a(11) is resolved (see below) |
| **ZZ12** free | [A316192] | perim n | a(3..10) = 1, 3, 4, 22, 69, 418, 2210, 14024 ✅ | **a(11..14) = 89075, 597581, 4076855, 28499301** (sum 1..14 = 33,279,563, matches the streaming-pipeline total) |

[A266549]: https://oeis.org/A266549
[A284869]: https://oeis.org/A284869
[A316192]: https://oeis.org/A316192
[A316198]: https://oeis.org/A316198
[A316200]: https://oeis.org/A316200

**Holes-allowed siblings.** Two OEIS sequences enumerate the same
objects as our ZZ4 / ZZ6 free matches but *additionally* include
polygons with holes (polyominoes with internal cavities):

| Ring | "Holes allowed" sibling | First divergence |
|------|--------------------------|------------------|
| ZZ4  | [A057730] vs A266549 | n = 8 (smallest hole: 8-cell ring around 1 cell, perim 16) |
| ZZ6  | [A057729] vs A284869 | n = 11 (smallest triangular hole) |

A057730 / A057729 are direct cross-checks for our enumeration at
small n, where holes can't fit yet. Verifying our match up to the
first divergence boundary is free; we already implicitly do for ZZ4
n ≤ 7 (matches A266549 = A057730 there). After divergence the
sequences diverge by *exactly* the count of holed polyominoes, so
either side gives a way to enumerate those if wanted.

[A057729]: https://oeis.org/A057729
[A057730]: https://oeis.org/A057730

## Discrepancy under investigation: A316200(11)

We get **9883**; Hugo's data says **19405**. Five independent paths
of ours agree on 9883:

- `--ring 10 --step 1 --free -n 11`
- `--ring 20 --step 2 --free -n 11` (same set, separate code path)
- single-threaded, no prunes
- 16-threaded, full prunes
- full streaming pipeline (stream → merge → build)

The OEIS value sits between our free count (9883) and our
rotation-canonical count (118 achiral + 2 × 9765 chiral pairs = 19648)
and doesn't match any clean orbit-counting formula — suggesting a
one-off computation or transcription error rather than a definitional
difference. The OEIS record is unchanged since revision 9
(2018-07-07); Hugo's enumerator is unpublished so we cannot audit it.

`oeis_a316200_zz10_pin` pins only **a(4..10)** so the test suite
doesn't lock in the disputed value. Resolution requires a third
independent source (hand enumeration of perim-11 ZZ10 polygons, or a
separate brute-force oracle) before submitting an OEIS correction.

## Rings with no OEIS oracle for full-ring counts

Candidates for new OEIS submissions. Our streaming pipeline can
produce the full per-length table out to whatever perim the compute
budget allows; submission needs just the sequence + a methodology
citation.

| Ring | Notes |
|------|-------|
| ZZ14 | 28-gon. Heptagonal `Z[zeta_14]`, cubic real subring `Z[2·cos(π/7)]`. A316197 is the bipartite-subset variant, not the full ring. |
| ZZ16 | 32-gon. Octagon-friendly tilings. |
| ZZ18 | 36-gon. Nonagonal `Z[zeta_18]`, cubic real subring `Z[2·cos(π/9)]`. Contains ZZ6 (6 \| 18); `--step 3` enumerates ZZ6 polygons through ZZ18's machinery, giving an external cross-check (`cross_validate_zz18_step3_matches_zz6`). A316199 is the bipartite-subset variant. |
| ZZ20 | 40-gon. ZZ20 step=2 collapses to ZZ10 step=1; step=4 to ZZ5. |
| ZZ24 | 48-gon. step=2 = ZZ12, step=3 = ZZ8. |
| ZZ32 | 64-gon. |
| ZZ60 | 120-gon. Contains ZZ4/ZZ12 substructure. |

ZZ14 and ZZ18 (added 2026-06) need their own sign-extraction
machinery because their real subrings are generated by roots of
irreducible cubics — no nested-sqrt closed form. We use
**adaptive-width bisection** in i128, terminating when
`min(|f(lo)|, |f(hi)|) > L · width` (provably guarantees `f` has
constant sign on the interval), with a `Ratio<BigInt>` Sturm-Tarski
fallback for inputs that exceed the i128 budget. The
`sturm_matches_bisection_*` fuzz tests pin equivalence between the
two paths. See `src/cyclotomic/sign.rs`.

**ZZ22 was dropped from scope**: no OEIS full-ring oracle (A316201
is bipartite-subset only) *and* no subring path that would give an
external cross-check, so the correctness story would rest entirely
on internal consistency.

## Bipartite-subset variants (A316195/197/199/201)

A separate object Hugo studied: matchstick polygons with the turn
set restricted to **odd multiples of π/k** only (no zero turn, no
even-multiple turns). The cumulative-direction parity flips at every
step, so the walks live on a bipartite sublattice — published as
distinct sequences:

| OEIS | Ring | Restricted turn set |
|------|------|---------------------|
| [A316195] | ZZ10 | ±π/5, ±3π/5 (4 turn angles vs ZZ10 full 9) |
| [A316197] | ZZ14 | ±π/7, ±3π/7, ±5π/7 (6 vs 13) |
| [A316199] | ZZ18 | ±π/9, ±3π/9, ±5π/9, ±7π/9 (8 vs 17) |
| [A316201] | ZZ22 | ±π/11, ±3π/11, ±5π/11, ±7π/11, ±9π/11 (10 vs 21) |

[A316195]: https://oeis.org/A316195
[A316197]: https://oeis.org/A316197
[A316199]: https://oeis.org/A316199
[A316201]: https://oeis.org/A316201

Matching any of these needs both the underlying ring AND a turn-set
restriction mechanism in `rat_enum`. Today's `--step N` keeps turns
that are *multiples of* N (the even-only subset for step=2); the
complement (odd-only) isn't expressible. Not on the roadmap.

## Other symmetry reductions (fixed, one-sided)

Hugo's whole family ("Verified matches" above) counts polygons up to
the full dihedral group of the lattice — rotations and reflections
both collapsed. Two coarser reductions exist:

- **Fixed** (translation only — every rotation and reflection
  counted as distinct)
- **One-sided** (rotation collapsed, mirror images distinct)

For SAPs *indexed by perimeter* (the objects our DFS enumerates),
the fixed counterpart exists in OEIS only for the two best-studied
lattices (ZZ4 and ZZ6); the one-sided counterpart isn't catalogued
anywhere we could find:

| Ring | Fixed (translation only) | One-sided (rot reduced) |
|------|--------------------------|-------------------------|
| ZZ4 (square) | [A002931] | **not in OEIS** |
| ZZ6 (triangular) | [A036418] | **not in OEIS** |
| ZZ8/10/12/14/16/18/20/24/32/60 | not in OEIS | **not in OEIS** |

The ratios A002931 / A266549 → 8 = \|D₄\| and A036418 / A284869 → 12
= \|D₆\| at large n, confirming these are the genuine fixed/dihedral
pairs over the same underlying object.

**Heads-up about polyominoes**: for polyominoes indexed by *cell
count* (a different object from our perimeter-indexed SAPs), the
full trichotomy IS catalogued — fixed [A001168], one-sided
[A000988], free [A000105]. Hooking our enumeration up to those
would require post-processing by cell count rather than perimeter,
which we don't currently do.

From our free output, both companions fall out for free if we
want to publish them:

- One-sided: from free count `D` (with `A` achiral entries),
  `one-sided = 2D − A`. Sample for ZZ10:

  | perim | D | A | one-sided = 2D − A |
  |------:|--:|--:|------:|
  |  6 |    10 |   5 |    15 |
  |  7 |    15 |   4 |    26 |
  |  8 |   124 |  34 |   214 |
  |  9 |   352 |  30 |   674 |
  | 10 |  2378 | 131 |  4625 |
  | 11 |  9883 | 118 | 19648 |

- Fixed: from the rotation-symmetry histogram (`--stats`), each
  polygon contributes `perim / rep_factor` translates. Reproduces
  A002931 (ZZ4) and A036418 (ZZ6) — useful third cross-check axis
  for both rings.

One further note:

- **[A346132]** lists ZZ12 step counts at which *no* closed walk
  exists. Trivially checkable: any n where ZZ12 free returns
  zero rats is in A346132. A "must-be-empty" pin, not a count pin.

[A001168]: https://oeis.org/A001168
[A000105]: https://oeis.org/A000105
[A000988]: https://oeis.org/A000988
[A002931]: https://oeis.org/A002931
[A036418]: https://oeis.org/A036418
[A346132]: https://oeis.org/A346132

## Not directly usable as oracles

OEIS sequences that came up during the crawl but enumerate
*different* objects:

- **Symmetric-only variants** (A316193, A316194, A316196): count
  only polygons with non-trivial symmetry. Derivable from our
  `--stats` rotational-symmetry histogram as a sub-count.
- **Open-walk siblings** (A306175, A306177-A306182): self-avoiding
  walks that don't have to close. Would need an `--mode open-walks`.
- **Cell-count polyomino enumerations** (A258206, A057779,
  A001168 / A000105 / A000988 family): indexed by cell count, not
  perimeter. Different indexing axis from our walks; useful only
  via post-processing by area.
- **Area / perimeter moments** (A056625, A056631, A056632, A057406,
  A056638): aggregate statistics; require post-processing.

## Test coverage

| Test | What it pins |
|------|-------------|
| `oeis_a266549_zz4_pin` | ZZ4 free a(2..7) (perim 4..14) |
| `oeis_a284869_zz6_pin` | ZZ6 free a(3..14) (12 terms — deepest oracle) |
| `oeis_a316198_zz8_pin` | ZZ8 free a(2..6) (perim 4..12); a(7)=240549 documented as our extension |
| `oeis_a316200_zz10_pin` | ZZ10 free a(4..10) **only** (a(11) skipped pending OEIS resolution) |
| `oeis_a316192_each_opt_combo` | ZZ12 free a(3..10) under every prune subset (~60–90 s) |

All live in `src/bin/rat_enum.rs` under `opt_correctness_tests`. Run
with `cargo test --release --lib --bin rat_enum -- oeis_`.
