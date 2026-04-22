# Sign-Only Geometry Refactor Plan

## Background

The Redelmeier patch enumeration algorithm requires segment-crossing checks during growth.
Two implementations exist:

1. **Pos4/SR fast path** — uses hand-coded `[i32; 4]` positions with `(i64, i64)` SR type
   representing `a + b*sqrt(3)` with no denominator. Pure integer arithmetic.
2. **Generic cyclotomic path** — uses `geometry::intersect<ZZ12>` operating on `ZZ::Real = Z12`
   which stores `Ratio<i64>` coefficients. Generic but ~46x slower.

### Profile Comparison (hexagon patches, size 8)

**Pos4 path: 290ms**

| %    | Function                                    |
|------|---------------------------------------------|
| 54.2 | `check_new_segments_cross_existing_pos4`    |
| 9.1  | `lex_min_rot`                               |
| 5.1  | Vec allocation                              |
| 2.8  | HashMap insert                              |

**Cyclotomic path: 13.5s**

| %    | Function                                    |
|------|---------------------------------------------|
| 6.1  | `Ratio::reduce` (GCD normalization)         |
| 5.8  | `zz_partial_signum_2_sym`                   |
| 5.1  | `zz6_mul` (Ratio multiply)                  |
| 4.6  | `Ratio::from_f64`                           |
| 3.9  | `Ratio::add`                                |
| 3.9  | `Ratio::mul`                                |
| 3.7  | `check_segments_cross_cyclotomic`           |
| 3.6  | `geometry::intersect`                       |
| 3.6  | `Ratio::sub`                                |
| 2.7  | `linalg::wedge`                             |
| 2.0  | `linalg::is_ccw`                            |

The top 7 hotspots are all `Ratio<i64>` operations. Together they account for ~31% of
total samples, and they propagate through every downstream function (`wedge`, `is_ccw`,
`intersect`, `signum`).

## Why Scaled i64 Coefficients Don't Work for General Ring Types

Z12 values have the form `(a + b*sqrt(3)) / l` where `l = scaling_fac = 2`.

The coefficients stored in `Z12` are numerators with implicit denominator 2.
When you multiply two such values:

```
(a + b*sqrt(3))/2  *  (c + d*sqrt(3))/2  =  (ac + 3bd + (ad+bc)*sqrt(3)) / 4
```

The denominator **doubles** from 2 to 4. Crucially, the numerator isn't always divisible
by the scaling factor. Counterexample:

```
sqrt(3)/2 * sqrt(3)/2 = 3/4
```

Here `zz6_mul` produces numerator `3` for the first coefficient — not divisible by `l=2`.
So you cannot normalize back to denominator 2 with integer-only coefficients.
The denominator must grow, which is exactly what `Ratio<i64>` handles.

**This is why the manual `GaussInt<i64>`-based ZZ12 did not help**: even though `ZZ12`
itself used integer coefficients, `ZZ12::Real = Z12` is **always** the macro-generated
`GaussInt<Ratio<i64>>` type. The entire geometry pipeline (`wedge` -> `dot` -> `is_ccw`
-> `is_between` -> `signum`) operates on `Z12`, not `ZZ12`. So Ratio overhead was
unavoidable.

## Key Insight: Sign-Only Geometry

The geometry pipeline only needs **signs** of Real values:

- `is_ccw` checks `wedge().is_positive()`
- `is_between` checks `wedge().is_zero() && dot().is_positive()`
- `is_colinear` checks `(a*re + b*im + c).is_zero()`

Since all denominators in the ring representation are positive,
**`sign(numerator/denom) = sign(numerator)`**.

This means we can work with **raw i64 numerators** and skip Ratio entirely.
No division, no GCD, no normalization.

The multiplication functions (`zz6_mul`, `zz8_mul`, `zz12_mul`, etc.) are **already
generic** over `T: IntRing + FromPrimitive` — they work with `i64` directly.

## Ring Classification

### Category A: Integer root squares — exact integer sign

| Ring  | Root squares           | # roots | Sign function          | `signum_sum_sqrt_expr` works? |
|-------|------------------------|---------|------------------------|-------------------------------|
| ZZ4   | `[1]`                  | 1       | `signum_1_sym`         | trivial (sign of single int)  |
| ZZ6   | `[1, 3]`               | 2       | `signum_2_sym`         | yes: `sign(a + b*sqrt(3))`    |
| ZZ8   | `[1, 2]`               | 2       | `signum_2_sym`         | yes: `sign(a + b*sqrt(2))`    |
| ZZ12  | `[1, 3]` (same as ZZ6) | 2       | `signum_2_sym`         | yes                           |
| ZZ24  | `[1, 2, 3, 6]`         | 4       | `signum_4_sym`         | yes: k=1, l=m*n (2*3=6)       |

For these rings, the existing `signum_sum_sqrt_expr_{2,4}` functions work directly with
`i64` arguments. The sign computation is **exact** and uses only integer comparisons
(with squaring to avoid irrationals).

### Category B: Nested radical roots — f64 sign, but still no Ratio

| Ring  | Root squares include              | # roots | Current sign function | Proposed strategy              |
|-------|-----------------------------------|---------|-----------------------|--------------------------------|
| ZZ10  | `2(5-sqrt(5))`                    | 4       | f64 fallback          | i64 mul -> f64 sign            |
| ZZ16  | `2+sqrt(2)`, `2(2+sqrt(2))`       | 4       | `signum_4_sym`*       | i64 mul -> f64 sign            |
| ZZ20  | same as ZZ10                      | 4       | f64 fallback          | i64 mul -> f64 sign            |
| ZZ30  | `2(5-sqrt(5))`, etc.              | 8       | f64 fallback          | i64 mul -> f64 sign            |
| ZZ32  | `2+sqrt(2+sqrt(2))`, etc.         | 8       | f64 fallback          | i64 mul -> f64 sign            |
| ZZ60  | all of the above                  | 8       | f64 fallback          | i64 mul -> f64 sign            |

*Note: ZZ16 currently uses `signum_4_sym` with `Ratio<i64>` approximations of irrational
root squares (e.g., `from_f64(2+sqrt(2))`). This is **already inexact** — the Ratio
is an approximation of an irrational. Switching to f64 sign from i64 coefficients would
be equally accurate and simpler.

For Category B, the sign check converts i64 coefficients -> f64 -> signum. This is
identical in accuracy to the current approach but avoids all Ratio overhead.

### Summary

| Category | Ratio overhead eliminated | Sign accuracy     |
|----------|--------------------------|-------------------|
| A        | ALL (mul + sign)         | exact             |
| B        | ALL (mul only, sign f64) | same as current   |

## How Pos4/SR Already Does This

The Pos4/SR approach is a working proof of concept for ZZ12:

- `SR(i64, i64)` represents `a + b*sqrt(3)` with no denominator
- Multiplication: `SR(a*c + 3*b*d, a*d + b*c)` — stays in integers
- Sign check (`sign_sqrt3`): compares `a^2` vs `3*b^2` via squaring — avoids irrationals
- `Pos4 = [i32; 4]` encodes scaled real/imaginary parts directly

This is exactly the integer-numerator approach, but hardcoded for ZZ12 only.

## Proposed Design

### Trait for sign-only real arithmetic

```rust
/// Integer-only representation of a cyclotomic real value for sign-only operations.
/// Stores raw i64 numerators without the common denominator.
trait SignReal: Clone + Copy + Sized {
    fn zero() -> Self;
    fn add(&self, other: &Self) -> Self;
    fn sub(&self, other: &Self) -> Self;
    fn mul(&self, other: &Self) -> Self;
    fn neg(&self) -> Self;
    fn is_zero(&self) -> bool;
    fn is_positive(&self) -> bool;
}
```

### Generic geometry functions

Rewrite `wedge`, `dot`, `is_ccw`, `is_between`, `intersect` to be generic over `SignReal`:

```rust
fn wedge_fast<R: SignReal>(re: impl Fn(&Pos) -> R, im: impl Fn(&Pos) -> R,
                            p1: &Pos, p2: &Pos) -> R {
    re(p1).mul(&im(p2)).sub(&im(p1).mul(&re(p2)))
}

fn is_ccw_fast<R: SignReal>(re: impl Fn(&Pos) -> R, im: impl Fn(&Pos) -> R,
                             p: &Pos, a: &Pos, b: &Pos) -> bool {
    let v = sub_pos(a, p);
    let w = sub_pos(b, p);
    wedge_fast(re, im, &v, &w).is_positive()
}

fn intersect_fast<R: SignReal>(re: impl Fn(&Pos) -> R, im: impl Fn(&Pos) -> R,
                                a: &Pos, b: &Pos, c: &Pos, d: &Pos) -> bool {
    // same logic as geometry::intersect, using SignReal operations
}
```

### Per-ring SignReal implementations

Each ring implements `SignReal` on a fixed-size i64 array:

```rust
// ZZ12: 2 roots, i64 coefficients for a + b*sqrt(3)
#[derive(Clone, Copy)]
struct Z12Sign([i64; 2]);

impl SignReal for Z12Sign {
    fn mul(&self, other: &Self) -> Self {
        let [a, b] = self.0;
        let [c, d] = other.0;
        Z12Sign([a*c + 3*b*d, a*d + b*c])  // zz6_mul::<i64>()
    }

    fn is_positive(&self) -> bool {
        let [a, b] = self.0;
        signum_sum_sqrt_expr_2(a, 1i64, b, 3i64) > 0
    }
    // ...
}

// ZZ24: 4 roots, i64 coefficients for a + b*sqrt(2) + c*sqrt(3) + d*sqrt(6)
#[derive(Clone, Copy)]
struct Z24Sign([i64; 4]);

impl SignReal for Z24Sign {
    fn mul(&self, other: &Self) -> Self {
        // use zz24_mul::<i64>()
    }
    fn is_positive(&self) -> bool {
        // use signum_sum_sqrt_expr_4 with integer root squares
    }
}

// ZZ16: 4 roots with nested radicals -> f64 sign fallback
#[derive(Clone, Copy)]
struct Z16Sign([i64; 4]);

impl SignReal for Z16Sign {
    fn is_positive(&self) -> bool {
        // compute value as f64 from i64 coefficients and root values
        let val: f64 = self.0.iter().zip(ROOT_VALUES).map(|(c, r)| *c as f64 * r).sum();
        val > 0.0
    }
}
```

### Position storage

The growing algorithm builds positions by starting at origin and adding unit vectors.
Unit coefficients are already integers (from `ccw_unit_coeffs`). Track positions as
i64 arrays to avoid Ratio entirely:

```rust
// Generic position for any ring with N symbolic roots
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
struct IntPos<const N: usize> {
    // [re_0, re_1, ..., re_{N-1}, im_0, im_1, ..., im_{N-1}]
    coeffs: [i64; 2 * N],
}

impl<const N: usize> IntPos<N> {
    fn add_unit(&mut self, direction: i8, unit_coeffs: &[[i64; 2]]) {
        for i in 0..N {
            self.coeffs[i] += unit_coeffs[direction as usize][0]; // real part
            self.coeffs[N + i] += unit_coeffs[direction as usize][1]; // imag part
        }
    }
}
```

### Integration with existing code

The `SignReal` implementations and generic geometry functions would live in a new module
(e.g., `src/cyclotomic/sign_geom.rs`). The growing algorithm would use `intersect_fast`
instead of `geometry::intersect` when the feature is enabled.

The existing generic `geometry::intersect` would remain unchanged for non-performance-
critical use cases (tests, one-off computations, etc.).

## What This Eliminates

| Hotspot (current profile)    | After refactor                     |
|------------------------------|------------------------------------|
| `Ratio::reduce` (6.1%)      | eliminated - i64 arithmetic        |
| `Ratio::from_f64` (4.6%)    | eliminated - constants are i64     |
| `Ratio::add/mul/sub` (~11%) | eliminated - i64 ops               |
| `Ratio::cmp` (1.5%)         | eliminated - no Ord needed         |
| `zz_partial_signum_2_sym`    | much faster - direct i64 sign      |
| `zz6_mul` (Ratio multiply)  | much faster - i64 multiply         |

Expected result: near-Pos4 performance for Category A rings, significant improvement
for Category B rings.

## Overflow Considerations

i64 overflow must be considered for large patch sizes or rings with many roots.

**ZZ12 estimation** (scaling_fac = 2, 2 roots):
- Initial numerators: ~2 (from unit coefficients)
- After k multiplications: ~O(k^2) growth (each mul adds products of 2 terms)
- For size-8 patches (~7 multiplications in path): numerators ~O(100)
- Wedge/dot products involve 2 multiplications: ~O(10,000)
- Well within i64 range (9.2e18)

**ZZ60 estimation** (scaling_fac = 16, 8 roots):
- Much faster growth due to 8-root multiplication formula
- May need i128 for safety, or early rescaling
- Should benchmark to determine practical limits

Mitigation options:
- Use `i128` for safety (2x memory but still no Ratio/GCD)
- Add overflow checks with `checked_mul` in debug builds
- Limit patch sizes per ring based on empirical overflow testing

## Implementation Steps

1. **Create `SignReal` trait and generic geometry functions** in `sign_geom.rs`
2. **Implement `SignReal` for Category A rings** (ZZ4, ZZ6, ZZ8, ZZ12, ZZ24)
   using exact integer sign functions
3. **Implement `SignReal` for Category B rings** (ZZ10, ZZ16, ZZ20, ZZ30, ZZ32, ZZ60)
   using f64 sign fallback from i64 coefficients
4. **Create `IntPos<N>` generic position type** with unit-based construction
5. **Integrate with growing algorithm** — use `intersect_fast` + `IntPos` instead of
   `geometry::intersect` + Ratio-based positions
6. **Benchmark all ring types** and compare against current implementation
7. **Remove Pos4/SR specialized code** once generic path is fast enough

## Open Design Questions

- **Position type in the growing algorithm**: Should `IntPos<N>` replace the current
  `ZZ12` positions entirely, or should it be a parallel representation used only for
  geometry?
- **Macro vs manual implementation**: Should `SignReal` implementations be generated by
  macro (like `impl_symnum!`) or written manually (only ~10 types)?
- **Long-term refactor**: Should the entire ring type hierarchy be made parametric over
  the scalar type (`Ratio<i64>` vs `i64`)? This would be the cleanest design but is a
  massive refactor touching `symnum.rs`, `types.rs`, and all downstream code.
