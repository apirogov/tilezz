# OEIS cross-references for cyclotomic matchstick polygons

This document tracks every OEIS sequence we've identified as relevant to
the polygon enumeration in this crate, along with our verification status
and any data we can contribute back.

The full canonical reference for "matchstick polygons" — closed
self-avoiding polygonal walks with unit-length sides on a cyclotomic
lattice — is **Hugo Pfoertner's 2018 family of sequences**, all
submitted on 2018-07-07 around the time of the [Al Zimmermann "Snakes
on a Plane" programming contest][snakes]. Hugo's source code is not
public; the only artifact is the OEIS data plus illustration pages at
[randomwalk.de][rwd]. Each sequence in the family encodes the
**dihedral-canonical** count: distinct polygons up to both rotation and
reflection of the underlying angle sequence.

[snakes]: http://www.recmath.org/contest/Snakes/index.php
[rwd]: http://www.randomwalk.de/sequences/

## Verified matches

For each ring we implement, this table compares our enumeration counts
against Hugo's OEIS data, term-by-term. All cells marked ✅ match
exactly under multiple independent enumeration paths (`--mode bench`
with and without prunes, `--mode stream`/`merge`/`build`, and where
applicable a subring-step cross-check like ZZ20 step=2 = ZZ10 step=1).

| Ring | OEIS | Description | OEIS terms (verified by us) | New terms we can contribute |
|------|------|-------------|------------------------------|-----------------------------|
| **ZZ4** dihedral | [A266549] | square lattice, perim 2n | a(2..7) = 1, 1, 3, 6, 25, 86 ✅ | a(8+) trivially (ring is small) |
| **ZZ8** dihedral | [A316198] | perim 2n | a(2..6) = 2, 6, 59, 695, 12198 ✅ | **a(7) = 240549** (new, perim 14) |
| **ZZ10** dihedral | [A316200] | perim n | a(4..10) = 2, 2, 10, 15, 124, 352, 2378 ✅ | **a(12..14)** once a(11) is resolved (see below) |
| **ZZ12** dihedral | [A316192] | perim n | a(3..10) = 1, 3, 4, 22, 69, 418, 2210, 14024 ✅ (existing test) | **a(11) = 89075, a(12) = 597581, a(13) = 4076855, a(14) = 28499301** (sum 1..14 = 33,279,563, matches the streaming-pipeline total) |

[A266549]: https://oeis.org/A266549
[A316192]: https://oeis.org/A316192
[A316198]: https://oeis.org/A316198
[A316200]: https://oeis.org/A316200

## Discrepancy under investigation: A316200(11)

We say a(11) = **9883**. Hugo's data says a(11) = **19405**.

Evidence we are correct:

- `--ring 10 --step 1 --dihedral -n 11` → 9883
- `--ring 20 --step 2 --dihedral -n 11` (mathematically the same set, completely different code path) → 9883
- single-threaded, no prunes → 9883
- 16-threaded, full prunes → 9883
- full streaming pipeline (stream → merge → build) → 9883

Five independent paths agree. The OEIS value sits between our dihedral
count (9883) and our rotation-canonical count (118 achiral + 2 × 9765
chiral pairs = 19648), and matches no clean orbit-counting formula —
which suggests a one-off transcription or computation error rather than
a definitional difference. The OEIS record is unchanged since revision
9, 2018-07-07. Hugo's enumerator is unpublished so we cannot audit it.

Resolution path: cross-check against an **independent third source**
(hand enumeration of perim-11 ZZ10 polygons, or a separate brute-force
oracle) before submitting a correction to OEIS.

The test `oeis_a316200_zz10_pin` currently pins **only a(4..10)** —
the range where we agree — so test infrastructure doesn't lock in the
disputed value.

## Rings we have, OEIS doesn't

These are clean candidates for new OEIS submissions if you want to
publish:

| Ring | Notes |
|------|-------|
| ZZ16 | 32-gon angle subset (octagon-friendly tilings) |
| ZZ20 | 40-gon. Note ZZ20 step=2 collapses to ZZ10 step=1; step=4 to ZZ5 |
| ZZ24 | 48-gon. ZZ24 step=2 = ZZ12, step=3 = ZZ8 |
| ZZ32 | 64-gon |
| ZZ60 | 120-gon. Contains ZZ4/ZZ12 substructure |

For each, our streaming pipeline can produce the full per-length table
out to whatever perim we have compute budget for; submission requires
just the sequence of counts plus a citation of the methodology.

## Rings missing from our codebase, OEIS has oracles

These would unlock external cross-validation if we implement the rings
(see the existing `define_integral_zz!` macro and the discussion in the
ring-implementation notes for what's involved per ring):

| Ring | OEIS | Description | OEIS terms |
|------|------|-------------|------------|
| ZZ6 (triangular) | [A284869] | perim n | a(3..24), 22 terms |
| ZZ14 | [A316197] | perim 2n | a(3..7), 5 terms |
| ZZ18 | [A316199] | perim 2n | a(3..6), 4 terms |
| ZZ22 | [A316201] | perim 2n | a(3..6), 4 terms |

[A284869]: https://oeis.org/A284869
[A316197]: https://oeis.org/A316197
[A316199]: https://oeis.org/A316199
[A316201]: https://oeis.org/A316201

ZZ6 is the most useful addition by external-oracle coverage: A284869
goes out to a(24), giving us a deep independent check that none of the
other rings have. The rest are short pins but still independent.

## "Non-dihedrally-reduced" (one-sided / rotation-canonical) variants

By default our DFS produces rotation-canonical output (each polygon
represented once per cyclic-rotation class). The `--dihedral` flag adds
the reflection collapse on top, mapping each chiral pair to a single
rep. The OEIS family above is all dihedral.

OEIS has **no sequence** for the "rotation-canonical, chiral pairs
distinct" intermediate that any of our rings produce by default.

For each dihedral count `D` (with `A` achiral entries), the
rotation-canonical (one-sided) count is `O = A + 2 * (D - A) = 2D - A`.
Both `D` and `A` come straight out of `--mode bench --stats`. Sample
for ZZ10 dihedral:

| perim | dihedral D | achiral A | one-sided O = 2D − A |
|-------|----------:|---------:|---------------------:|
|  6    |        10 |        5 |                   15 |
|  7    |        15 |        4 |                   26 |
|  8    |       124 |       34 |                  214 |
|  9    |       352 |       30 |                  674 |
| 10    |      2378 |      131 |                 4625 |
| 11    |      9883 |      118 |                19648 |

If we want to publish one-sided variants for any of our rings, the data
falls out of every dihedral run for free via the histogram. Two further
points to note before submitting:

- **No-symmetry-reduced sequences exist for ZZ4 (A002931) and ZZ6
  (A036418)** — those count every (rotation, reflection) pair as
  distinct, which is even less reduced than one-sided. We can derive
  those totals from our dihedral output via the rotation-symmetry
  histogram in `--stats` (each polygon contributes `perim / rep_factor`
  fixed copies), giving a third cross-check axis for ZZ4 specifically.

- **A346132** is a follow-up by Hugo listing the step counts at which
  *no* closed walk exists on the ZZ12 lattice. Cross-checkable from our
  enumeration trivially: any n where ZZ12 dihedral returns zero rats is
  in A346132. Not a count-pin but a "must-be-empty" pin.

## Combinatorial objects we are NOT trying to enumerate

For completeness, these are OEIS sequences that came up during the
crawl but enumerate *different* objects than our matchstick polygons,
so they are not directly usable as oracles for our DFS:

- **Symmetric-only variants** (A316193, A316194, A316196): count only
  polygons with non-trivial symmetry. Our `--stats` rotational-symmetry
  histogram already breaks this out, so we could check derived totals,
  but it's a sub-count.
- **Open-walk siblings** (A306175, A306177-A306182): self-avoiding
  walks that don't have to close. Different DFS — we'd need to add an
  `--mode open-walks` if we wanted to enumerate them.
- **Polyhexes / polyominoes** (A258206, A057729, A057779): hex / square
  cells with shared edges, not polygons drawn by walks. Different
  combinatorial object.
- **Area / perimeter moments** (A056625, A056631, A056632, A057406,
  A056638): aggregate statistics over polygon families. Could be
  computed from our enumeration but require post-processing.

## Test coverage

| Test | Location | What it pins |
|------|----------|-------------|
| `oeis_a316192_each_opt_combo` | `src/bin/rat_enum.rs` `opt_correctness_tests` | ZZ12 dihedral a(3..10) under every prune subset (~60–90 s) |
| `oeis_a266549_zz4_pin` | same | ZZ4 dihedral a(2..7) (perim 4..14) |
| `oeis_a316198_zz8_pin` | same | ZZ8 dihedral a(2..6) (perim 4..12), with a(7)=240549 documented as our extension |
| `oeis_a316200_zz10_pin` | same | ZZ10 dihedral a(4..10) **only** (a(11) skipped pending OEIS resolution) |

Run with `cargo test --release --lib --bin rat_enum -- oeis_`.
