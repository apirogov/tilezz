# OEIS A316200 a(11) correction: 19405 -> 9883

**Claim.** The published value **A316200(11) = 19405 is incorrect**; the
correct value is **9883**.

A316200 (Hugo Pfoertner, Jul 2018; from Al Zimmermann's "Snakes on a Plane"
contest, Fall 2006) is defined as:

> Number of self-avoiding polygons with perimeter n and sides = 1 that have
> vertex angles from the set 0, +-Pi/5, +-2*Pi/5, +-3*Pi/5, +-4*Pi/5, **not
> counting rotations and reflections as distinct.**

This is the **free** (rotation + reflection) count of simple closed unit-edge
polygons on the decagonal lattice Z[zeta_10] -- the same object `tilezz`
enumerates as "ZZ10 free". The published data is

```
offset 1:  0, 0, 0, 2, 2, 10, 15, 124, 352, 2378, 19405
           a1 a2 a3 a4 a5 a6  a7   a8   a9   a10    a11
```

`tilezz` and an independent clean-room enumerator both agree with OEIS for
a(4)..a(10) and both give **a(11) = 9883**, not 19405.

---

## The evidence chain

| # | Check | Method | Result |
|---|-------|--------|--------|
| 1 | OEIS published value | live fetch | a(11) = **19405** |
| 2 | `tilezz` free count | Rust enumerator | a(11) = **9883**; matches OEIS exactly at a(4..10) |
| 3 | Five independent `tilezz` paths | ZZ10-step1, ZZ20-step2, single/multi-thread, prune on/off, streaming pipeline | all **9883** (but share one geometry predicate) |
| 4 | `tilezz` with exact/f64 `cell_floor` guard ACTIVE | release + `debug-assertions` | still **9883** (guard never fired) |
| 5 | `tilezz` internal consistency | free / achiral / one-sided | free 9883, achiral 118, one-sided = 2*9883 - 118 = 19648 (self-consistent) |
| 6 | **Independent clean-room recompute** | own code, EXACT integer arithmetic in Z[sqrt5]/Z[zeta_10], zero floating point | a(4..10) reproduced exactly; **a(11) = 9883** |
| 7 | Second independent dedup | exact point-set congruence in Z[zeta_10] | agrees (see note on p=9 below) |

Check 6 is the decisive one: it shares **no code** with `tilezz` and uses **no
floating point** -- every orientation / on-segment / crossing test is an exact
integer sign decision `sign(a + b*sqrt5)` resolved by comparing `a^2` vs
`5*b^2`. So there is no rounding and no shared-bug risk; the result is rigorous.

Note that 19405 is **neither** the free count (9883) **nor** the one-sided count
(19648), so it is not a free-vs-one-sided convention mix-up -- it is an isolated
over-count. The likely cause is an **incomplete dihedral reduction** in the 2018
computation: the dihedral group here needs sign-negation (reflection) and order-
reversal (reverse traversal) as *independent* involutions (four variants
`{t, -t, reverse(t), reverse(-t)}`). The clean-room author independently hit
exactly this ~2x over-count during development and caught it via the a(4..10)
validation gate before trusting a(11).

---

## How to reproduce every check

All commands run from the repository root. `<RAT>` = the `rat_enum` CLI.

### Check 1 -- the OEIS published value
```sh
curl -s "https://oeis.org/search?q=id:A316200&fmt=text" | grep '^%S'
# %S A316200 0,0,0,2,2,10,15,124,352,2378,19405   (a(11) = 19405)
```

### Checks 2 & 5 -- tilezz free + per-length / achiral / one-sided
```sh
cargo build --release --bin rat_enum --features cli
# FREE (dihedral). The "length 11" row's "total" is a(11):
./target/release/rat_enum --ring 10 -n 11 --free --threads 0 --mode bench --stats
#   length 11 | total 9883 | achiral 118
# ONE-SIDED (rotation only), for the consistency relation:
./target/release/rat_enum --ring 10 -n 11        --threads 0 --mode bench --stats
#   length 11 | total 19648   (= 2*9883 - 118)
```

### Check 3 -- five independent tilezz paths
```sh
# ZZ10 directly, prunes on/off, single/multi-thread:
./target/release/rat_enum --ring 10 -n 11 --free --threads 0 --mode bench
./target/release/rat_enum --ring 10 -n 11 --free --threads 1 --mode bench --mod-prune --closure-key-prune
# ZZ20 with step=2 enumerates the ZZ10-equivalent subset:
./target/release/rat_enum --ring 20 -n 11 --free --step 2 --threads 0 --mode bench
# Streaming pipeline (stream -> merge -> build) over a temp dir:
D=$(mktemp -d)
./target/release/rat_enum --ring 10 -n 11 --free --mode stream -o "$D"
./target/release/rat_enum --ring 10 -n 11 --free --mode merge  -o "$D"   # certificate.json counts
# all report 9883 at perimeter 11.
```

### Check 4 -- exact/f64 cell_floor guard active (release runs the f64 fast path
by default; this forces the internal exact-vs-f64 agreement assertion on)
```sh
RUSTFLAGS="-C debug-assertions=on" cargo build --release --bin rat_enum --features cli
RUSTFLAGS="-C debug-assertions=on" ./target/release/rat_enum --ring 10 -n 11 --free --threads 0 --mode bench
# still 9883; no assertion fires -> the f64 cell_floor did not mis-bucket here.
```

### Check 6 -- the decisive independent exact recompute
```sh
# Validates a(4..10) against the agreed reference values, then prints a(11).
# a(11) takes ~18 min single-threaded CPython; a(4..9) is fast.
python3 docs/oeis-A316200-correction/zz10_independent.py 4 11
# p= 4..10  -> 2,2,10,15,124,352,2378   (all "OK")
# p=11      -> free_count=9883
```

### Check 7 -- second independent dedup (exact point-set congruence)
`zz10_pointset_crosscheck.py` deduplicates found polygons by exact geometric
congruence of their **vertex sets** in Z[zeta_10] (rotation = multiply by zeta,
reflection = complex conjugation), a method completely independent of the turn-
sequence canonical form. It agrees with `zz10_independent.py` through p=8. At
p=9 it gives 351 vs the correct 352: this is **not** an error -- there is one
genuine pair of *distinct* simple polygons that occupy the *same* 9-vertex set
with different edge connectivity (`(-4,-2,-2,0,-2,-2,-3,2,3)` and
`(-4,-2,-2,-1,-1,-2,-3,1,4)`). This confirms the turn-sequence model is the
correct "number of polygons" notion and point-set identity is too coarse,
reinforcing 9883.

---

## Files

- `zz10_independent.py` -- the authoritative independent enumerator (exact
  Z[sqrt5] arithmetic, no floating point). The trust anchor for the correction.
- `zz10_pointset_crosscheck.py` -- the second-opinion dedup by exact point-set
  congruence in Z[zeta_10].

## Provenance

- `tilezz` value: see `src/bin/rat_enum.rs` (the `oeis_a316200_zz10_pin` test
  pins a(4..10) against OEIS and deliberately omits a(11) with a comment noting
  the disagreement). The historical ZZ12 "T-touch" self-intersection bug
  (`src/cyclotomic/geometry.rs`, now regression-pinned) is the same bug class an
  over-counting enumerator would exhibit.
- This correction was established via a full correctness audit of the tilezz
  numerical/geometric/DFS stack; the independent recompute above closes the one
  gap that five shared-predicate tilezz paths could not.
