# A284869 a(22): tilezz gets 374128154, OEIS lists 374128188 (UNRESOLVED)

**Status: OPEN / under investigation. This is NOT (yet) a correction claim.**
Unlike the A316200 case (where the published value was an outlier and an
independent clean-room enumerator confirmed tilezz), here the published value
is recent and from a highly reliable source, and we have *not* established
which side is right. This document is the running map of evidence so it can
become either a **bug post-mortem** (if tilezz is undercounting) or a
**correction case** (if the published term is wrong).

## The discrepancy

A284869 -- "Number of n-step 2-dimensional closed self-avoiding paths on the
triangular lattice, reduced for symmetry" -- is the **free** (rotation +
reflection) count of simple closed unit-edge polygons on Z[zeta_6], i.e. what
tilezz enumerates as "ZZ6 free". `a(n)` indexes perimeter `n` directly.

| | a(22) (perimeter 22) |
|---|---|
| tilezz ZZ6 free | **374,128,154** |
| OEIS A284869 (published) | **374,128,188** |
| difference | **-34** (tilezz is lower) |

tilezz also yields a(23) = 1,390,909,413, but since a(22) is in dispute that
term is **not** trustworthy and is not being proposed as an extension.

## Provenance of the published term (sets the prior)

From the OEIS entry's extension history:

- a(1..14): Luca Petrone (2017)
- a(15): Hugo Pfoertner (Jun 2018)
- **a(16)..a(22): Walter Trump (Nov 2023)**, with a methodology PDF
  (`/A284869/a284869.pdf`, "Self-avoiding closed walks on a triangular
  lattice") and the comment "a(n) is the number of simply connected
  polyiamonds with perimeter n."

Walter Trump is an exceptionally careful enumerator with documented methods.
So the a-priori expectation is that **the published a(22) is right and tilezz
undercounts by 34** -- the opposite prior to the A316200 dispute. We treat it
that way until evidence says otherwise.

## What matches exactly

tilezz reproduces every published term from a(3) through **a(21)** -- 19
consecutive terms -- exactly. The divergence is at a(22) alone:

```
perim   published a(n)        tilezz ZZ6        ZZ12 step-2 (independent)
 15            45,645            45,645            45,645
 16           161,705           161,705           161,705
 17           575,325           575,325           575,325
 18         2,074,088         2,074,088         2,074,088
 19         7,521,818         7,521,818         7,521,818
 20        27,502,445        27,502,445        27,502,445
 21       101,134,999       101,134,999       101,134,999
 22       374,128,188       374,128,154       374,128,154   <-- diverge (-34)
```

A bug that corrupts only the single deepest term while reproducing 19 prior
terms exactly (including a(21) = 101,134,999) is unusual -- most enumeration /
canonicalization / arithmetic bugs perturb many lengths and grow with n. That
argues *against* a systematic tilezz bug, but does not prove anything.

## Evidence chain so far

| # | Check | Method | Result |
|---|-------|--------|--------|
| 1 | OEIS published value | live fetch (oeis.org JSON) | a(22) = **374,128,188** |
| 2 | tilezz ZZ6 free, full pipeline | bench (in-memory HashSet) == stream->merge cert == build n_sequences == count_by_length; cumulative <=23 = 1,904,072,327 = sum of per-length | a(22) = **374,128,154**, internally consistent |
| 3 | Independent ring cross-check | **ZZ12 step-2** = ZZ6 sub-ring; different ring arithmetic + different mod-prune moduli + different closure-key tables | a(22) = **374,128,154** (agrees with native ZZ6) |
| 4 | Duplicate-storage failure mode | a DAFSA stores each string at most once; intermediate dihedral-image dups collapse at merge and are idempotent in the build | **excluded** -- the gap is a true distinct-count difference, not a storage artifact |
| 5 | Non-canonical-duplicate failure mode | `tools/verify_canonical.py` independently recomputes each rat's dihedral-canonical CCW form; a wrong canonical would also have disagreed between the two rings in #3 | **excluded** (overcount-via-bad-canonicalization) |
| 6 | Build provenance | binary built at commit 89d34a0f, whose `src/`, `build.rs`, `Cargo.*` are byte-identical to the v0.1.0 tag (dd4fd14) | runs are v0.1.0-equivalent enumeration code |
| 7 | closure-key-prune soundness | ZZ6 n22 with **mod-prune only** (closure-key disabled; ~3.9 h, ~1.5e10 records emitted then deduped) | a(22) = **374,128,154** -- closure-key prune **EXONERATED** (count unchanged without it) |
| 8 | mod-prune cross-validation | native ZZ6 uses ZZ6 reachability moduli; ZZ12 step-2 (#3) uses ZZ12 moduli -- *different* mod-prune tables | both give **374,128,154**; an unsound mod-prune would have to drop the same 34 under two different moduli sets |

Note on #3: the two ring representations share the DFS *core* and the prune
*algorithm* (only the ring-specific constants/tables differ), so they are not
fully independent of a logic-level bug. But the same core enumerated ZZ4 to
perimeter 32 and matched A266549 a(1..16) **exactly** (435,646,127 polygons),
which exercises the core at large scale on another ring.

## Failure-mode status

| Hypothesis | Direction | Status |
|---|---|---|
| Duplicate strings stored | overcount | impossible by construction (DAFSA = set) |
| Non-canonical duplicate (bad canonicalization) | overcount | excluded (#5, and two rings agree) |
| Pipeline drops records (stream/merge/gzip) | undercount | unlikely: ZZ4 n32 (435M) matched OEIS exactly through the same pipeline; ZZ12 n13 stream==in-memory verified earlier |
| closure-key-prune unsound | undercount | **excluded** (#7: mod-prune-only gives the same 374,128,154) |
| mod-prune unsound | undercount | strongly constrained (#8: ZZ6 and ZZ12 moduli differ yet agree); a clean no-mod-prune run is ~infeasible at n22 |
| DFS core / geometry bug | either | strongly constrained: a(3..21) exact, ZZ4 n32 exact |
| Published a(22) wrong | (Trump high) | open; would require strong independent confirmation given the source |
| Different objects (polyiamond vs self-avoiding walk; vertex-pinch) | either | the smallest vertex-pinch ("bowtie", perimeter 6) is excluded by both (a(6) = 4 matches), so this is not a simple definitional gap |

## Where we stand

**Three** tilezz configurations now agree at 374,128,154 -- native ZZ6 (both
prunes), native ZZ6 (closure-key off, #7), and ZZ12 step-2 (different ring,
different prune tables, #3) -- versus the published 374,128,188 (Walter Trump,
2023), differing by 34 at perimeter 22 only. Every tilezz-side failure mode is
now excluded or strongly constrained: storage dups (impossible by
construction), non-canonical dups (excluded, #5), closure-key prune
(exonerated, #7), mod-prune (two different moduli agree, #8), DFS core (ZZ4 n32
= 435,646,127 matches A266549 exactly). Internal consistency holds (per-length
sum = n_sequences). No tilezz mechanism for the -34 has been found.

This is starting to resemble the A316200 situation (tilezz consistently right,
published term off) -- but the prior is reversed: the published a(22) is recent
and from a meticulous, documented source (Walter Trump, 2023), so the bar to
contradict it is high. We do NOT assert a correction. The decisive next step is
a genuinely independent enumeration sharing no code with tilezz -- the analogue
of the clean-room check that settled A316200 -- reaching perimeter 22, and/or
contacting Walter Trump with the discrepancy. Until one direction is
established independently, this stays OPEN.

## Next steps

1. Finish the mod-prune-only run (#7) and record the result here.
2. If closure-key is exonerated: assess mod-prune (cost permitting) and/or a
   clean-room enumerator anchored at a smaller perimeter, and cross-reference
   A057729 (the holes-allowed sibling) where it still coincides.
3. Do not submit anything to OEIS until one direction is established
   independently.

## Reproduce

tilezz ZZ6 free to perimeter 23 (streaming pipeline), then read per-length
terms (build rat_enum from a clean checkout; see the dataset's
`tools/reproduce.sh` and the v0.1.1 provenance guard):

```sh
rat_enum --ring 6 -n 23 --free --threads 0 --mod-prune --closure-key-prune --mode stream -o zz6
rat_enum --ring 6 -n 23 --free --threads 0 --mod-prune --closure-key-prune --mode merge -o zz6
rat_enum --ring 6 -n 23 --free --threads 0 --mod-prune --closure-key-prune --mode build --oeis-a-number A284869 -o zz6
python3 zz6/dafsa/tools/count_by_length.py zz6/dafsa   # a(22) line = 374128154
```

Independent ring cross-check (ZZ12 step-2 = ZZ6):

```sh
rat_enum --ring 12 --step 2 -n 22 --free --threads 0 --mod-prune --closure-key-prune --mode stream -o chk
# ... merge, build ... ; count_by_length a(22) line = 374128154
```

## Artifacts in this directory

- `zz6_native_bylen.txt` -- per-perimeter counts, native ZZ6 free, n<=23.
- `zz12step2_bylen.txt` -- per-perimeter counts, ZZ12 step-2, n<=22.
- `zz6_modpruneonly_bylen.txt` -- per-perimeter counts, ZZ6 with closure-key
  prune disabled, n<=22 (same a(22) = 374128154; #7 above).

## Timings (16-core / 58 GB box, this run)

- ZZ6 n23 free, full prunes: stream 9366 s + merge 1545 s + build 2040 s
  = **~3 h 36 m** (1.9e9 cumulative records; gzipped scratch peaked ~18 GB).
- ZZ12 step-2 n22 cross-check, full prunes: stream 2816 s + merge 373 s +
  build 532 s = **~1 h 2 m**.
- ZZ6 n22 mod-prune-only (closure-key off): **~3.9 h** (14170 s); ~1.5e10
  records emitted (the closure-key prune normally suppresses ~10x the
  emission) then merged/deduped to the same a(22) = 374,128,154.
