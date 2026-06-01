//! Generic engine for integer-basis cyclotomic rings.
//!
//! Every ring `ZZ_n` whose minimal polynomial of `zeta = e^(2*pi*i/n)` over
//! `Z` is the cyclotomic polynomial `Phi_n(x)` is stored as an integer
//! coefficient vector against the power basis `{1, zeta, ..., zeta^(phi(n)-1)}`.
//! The pieces that vary per ring are:
//!
//! * `n` -- the order of `zeta`,
//! * `phi(n)` -- the dimension of the basis (= storage width),
//! * the *reduction rule* for `zeta^phi(n)` back into the basis (the
//!   coefficient list of `-Phi_n(x)` up to the leading term), and
//! * a small set of *precomputed projections* (conjugation matrix, Cartesian
//!   embedding, real/imag-part decomposition into the ring's real subring,
//!   ...) that the macro-driven per-ring impls feed into the generic
//!   helpers.
//!
//! This module is the *engine* for that representation: it exposes pure
//! free functions parametric over `const PHI` (and where relevant a second
//! `const K = phi(n)/2` for the real subring), all operating on bare
//! `[i64; PHI]` storage and per-ring constant tables. Types are not
//! defined here -- they live in `rings.rs`, one per `define_integral_zz!`
//! invocation.
//!
//! All helpers are `#[inline]` and `const PHI`-generic so LLVM can fully
//! unroll the inner loops when `PHI` is known at the call site (it always
//! is, per-ring).

use num_complex::Complex64;

// ----------------
// Public trait for borrowing the raw integer-basis coefficients.

/// Borrowed-slice view of a ring element's integer-basis coefficients.
///
/// Each ring's inherent `int_coeffs()` returns `[i64; PHI]` by value,
/// which is per-type and so awkward in generic code. This trait
/// provides a uniform `&[i64]` accessor that callers can use without
/// committing to a fixed PHI; the slice has length PHI for the
/// underlying ring.
pub trait IntCoeffsSlice {
    fn int_coeffs_slice(&self) -> &[i64];
}

// ----------------
// Trivial linear ops.

/// Coefficient-wise addition.
#[inline]
pub fn add_basis<const PHI: usize>(x: &[i64; PHI], y: &[i64; PHI]) -> [i64; PHI] {
    let mut out = [0i64; PHI];
    let mut i = 0;
    while i < PHI {
        out[i] = x[i] + y[i];
        i += 1;
    }
    out
}

/// Coefficient-wise subtraction.
#[inline]
pub fn sub_basis<const PHI: usize>(x: &[i64; PHI], y: &[i64; PHI]) -> [i64; PHI] {
    let mut out = [0i64; PHI];
    let mut i = 0;
    while i < PHI {
        out[i] = x[i] - y[i];
        i += 1;
    }
    out
}

/// Coefficient-wise negation.
#[inline]
pub fn neg_basis<const PHI: usize>(x: &[i64; PHI]) -> [i64; PHI] {
    let mut out = [0i64; PHI];
    let mut i = 0;
    while i < PHI {
        out[i] = -x[i];
        i += 1;
    }
    out
}

// ----------------
// Multiplication, mod the cyclotomic polynomial.

/// Multiply two elements in the power basis `{1, zeta, ..., zeta^(PHI-1)}`,
/// reducing modulo `Phi_n(x) = 0`.
///
/// `reduction[k]` is the coefficient of `zeta^k` in the basis expansion of
/// `zeta^PHI`. For example, for ZZ12 with `Phi_12(x) = x^4 - x^2 + 1`,
/// `zeta^4 = zeta^2 - 1`, so `reduction = [-1, 0, 1, 0]`.
///
/// Strategy: compute the unreduced polynomial product into a scratch buffer
/// of length `2*PHI - 1`, then fold every power `>= PHI` into the lower
/// half by repeated application of the reduction rule. Because `reduction`
/// is applied at most `PHI - 1` times, the cost is `O(PHI^2)` total --
/// the same complexity as ZZ12's hand-written `Mul`.
///
/// LLVM unrolls all inner loops when `PHI` is a compile-time constant.
///
/// # Example
///
/// A tiny 2-element ring (think `Z[i]` with the basis `{1, i}` and the
/// reduction `i^2 = -1`, so `reduction = [-1, 0]`):
///
/// ```
/// use tilezz::cyclotomic::integral_basis::mul_basis;
/// let one_i = [0i64, 1];     // = i
/// let reduction = [-1i64, 0]; // i^2 = -1
/// // (1 + i) * (1 + i) = 1 + 2i + i^2 = 2i
/// let x = [1i64, 1];
/// let y = [1i64, 1];
/// assert_eq!(mul_basis::<2>(&x, &y, &reduction), [0, 2]);
/// // i * i = -1
/// assert_eq!(mul_basis::<2>(&one_i, &one_i, &reduction), [-1, 0]);
/// ```
#[inline]
pub fn mul_basis<const PHI: usize>(
    x: &[i64; PHI],
    y: &[i64; PHI],
    reduction: &[i64; PHI],
) -> [i64; PHI] {
    // Allocate the scratch buffer on the stack at size `2 * PHI`. We only
    // populate the first `2*PHI - 1` entries; the trailing slot is unused
    // and never read. Using `2 * PHI` here (rather than `2 * PHI - 1`) lets
    // us avoid `generic_const_exprs` -- on stable Rust we cannot name
    // `[i64; 2*PHI - 1]` directly from a generic.
    let mut scratch = [0i64; 32];
    debug_assert!(
        2 * PHI <= scratch.len(),
        "mul_basis: PHI={} too large for 32-slot scratch",
        PHI,
    );

    // Polynomial product.
    let mut i = 0;
    while i < PHI {
        let xi = x[i];
        if xi != 0 {
            let mut j = 0;
            while j < PHI {
                scratch[i + j] += xi * y[j];
                j += 1;
            }
        }
        i += 1;
    }

    // Reduce powers `>= PHI` down into `0..PHI`. We work from the highest
    // power downward, so each fold is a single application of the rule.
    // After folding `zeta^k` (for `k > PHI`), the contribution lives in
    // `scratch[k - PHI ..= k - 1]` (i.e. shifted into the next-lower power
    // window), which then gets folded in turn.
    //
    // Concretely: `zeta^k = zeta^(k - PHI) * zeta^PHI
    //                     = zeta^(k - PHI) * sum_j reduction[j] * zeta^j
    //                     = sum_j reduction[j] * zeta^(k - PHI + j)`.
    let mut k = 2 * PHI - 2;
    while k >= PHI {
        let coeff = scratch[k];
        if coeff != 0 {
            scratch[k] = 0;
            let base = k - PHI;
            let mut j = 0;
            while j < PHI {
                scratch[base + j] += coeff * reduction[j];
                j += 1;
            }
        }
        k -= 1;
    }

    // Extract the low `PHI` slots.
    let mut out = [0i64; PHI];
    let mut i = 0;
    while i < PHI {
        out[i] = scratch[i];
        i += 1;
    }
    out
}

// ----------------
// Conjugation (precomputed matrix).

/// Apply a precomputed conjugation matrix: returns `conj(z) = M * z` where
/// `M[i][j]` is the coefficient of `zeta^i` in the basis expansion of
/// `conj(zeta^j)`.
///
/// The matrix is derivable from `reduction` via [`derive_conj_matrix`].
/// ZZ12 specializes this with the inline formula
/// `conj((a, b, c, d)) = (a + c, b, -c, -b - d)` (see `rings.rs`).
#[inline]
pub fn conj_basis<const PHI: usize>(x: &[i64; PHI], conj_matrix: &[[i64; PHI]; PHI]) -> [i64; PHI] {
    let mut out = [0i64; PHI];
    let mut i = 0;
    while i < PHI {
        let mut acc: i64 = 0;
        let mut j = 0;
        while j < PHI {
            acc += conj_matrix[i][j] * x[j];
            j += 1;
        }
        out[i] = acc;
        i += 1;
    }
    out
}

// ----------------
// Cartesian projection.

/// Project an integral-basis element to Cartesian `Complex64` via a
/// precomputed per-basis-element complex value table.
///
/// `cartesian[k]` is the `Complex64` value of `zeta^k`. The result is
/// `sum_k x[k] * cartesian[k]`.
///
/// This is the default `complex64` path; rings that need a faster
/// projection (e.g. ZZ12 with its inline `Re`/`Im` formulas) provide a
/// hand-rolled `complex64_fn` instead.
#[inline]
pub fn complex64_basis<const PHI: usize>(
    x: &[i64; PHI],
    cartesian: &[Complex64; PHI],
) -> Complex64 {
    let mut acc = Complex64::new(0.0, 0.0);
    let mut i = 0;
    while i < PHI {
        let xi = x[i] as f64;
        acc += Complex64::new(xi * cartesian[i].re, xi * cartesian[i].im);
        i += 1;
    }
    acc
}

// ----------------
// Re/Im sign extraction via real-subring projection.
//
// For a cyclotomic ring `ZZ_n` with `phi(n) = 2K`, the real subring has
// dimension `K = phi(n)/2`. Re(zeta^k) and Im(zeta^k) both live in this
// real subring, and each can be written as a vector of `K` integer
// coefficients against some symbolic real-subring basis (the same basis
// the per-ring `signum_sum_sqrt_expr_*` function consumes).
//
// `re_decomp[k]` is the K-vector of coefficients of Re(zeta^k); likewise
// `im_decomp[k]`. The denominator implicit in these decompositions (e.g.
// the `/2` in ZZ12's `Re(zeta) = sqrt(3)/2`) is positive and shared across
// all entries, so it doesn't affect the sign and is omitted from the
// generic shape.

/// Sign of `Re(z)` via projection to the ring's real subring, then per-ring
/// `real_sign`.
///
/// `re_decomp[k]` is the K-vector of real-subring coefficients for
/// `Re(zeta^k)`. The accumulated K-vector for `Re(z)` is
/// `sum_k x[k] * re_decomp[k]`, fed to `real_sign`.
///
/// ZZ12 specializes this with the inline `sign_m_plus_n_sqrt3(2a + c, b)`
/// formula (see `rings.rs`).
#[inline]
pub fn re_sign_basis<const PHI: usize, const K: usize>(
    x: &[i64; PHI],
    re_decomp: &[[i64; K]; PHI],
    real_sign: fn(&[i64; K]) -> i8,
) -> i8 {
    let mut acc = [0i64; K];
    let mut i = 0;
    while i < PHI {
        let xi = x[i];
        if xi != 0 {
            let mut j = 0;
            while j < K {
                acc[j] += xi * re_decomp[i][j];
                j += 1;
            }
        }
        i += 1;
    }
    real_sign(&acc)
}

/// Sign of `Im(z)` via projection to the ring's real subring, then per-ring
/// `real_sign`. Same shape as [`re_sign_basis`] but indexes `im_decomp`.
#[inline]
pub fn im_sign_basis<const PHI: usize, const K: usize>(
    x: &[i64; PHI],
    im_decomp: &[[i64; K]; PHI],
    real_sign: fn(&[i64; K]) -> i8,
) -> i8 {
    let mut acc = [0i64; K];
    let mut i = 0;
    while i < PHI {
        let xi = x[i];
        if xi != 0 {
            let mut j = 0;
            while j < K {
                acc[j] += xi * im_decomp[i][j];
                j += 1;
            }
        }
        i += 1;
    }
    real_sign(&acc)
}

// ----------------
// Unit-segment intersection. This is the one canonical algorithm; every
// `IntersectUnitSegments` impl in the crate routes through here.

/// Decide whether two unit-length segments intersect, in exact ring
/// arithmetic over the integer basis. Touching only at shared endpoints
/// (`a == c` etc.) is NOT counted as intersection.
///
/// # Math
///
/// Given `s1 = (a, b)`, `s2 = (c, d)`, all four endpoints distinct:
///
/// ```text
///   uA    = b - a,  uB = d - c,  delta = c - a    (uA, uB are units)
///   v_z   = conj(uA) * delta     mul #1
///   k_z   = conj(uA) * uB        mul #2
///   w_z   = conj(delta) * uB     mul #3 (only when needed)
/// ```
///
/// Then `Im(v_z) = wedge(uA, delta)` is the sidedness of `c` vs line ab,
/// `Im(v_z + k_z) = wedge(uA, d - a)` of `d` vs ab, `Im(w_z) = wedge(delta, uB)`
/// of `a` vs cd, and `Im(w_z - k_z) = wedge(c - b, uB)` of `b` vs cd.
///
/// Case analysis:
///
/// * `sign(K) == 0`: lines parallel. If `sign(V) != 0` the lines are
///   distinct parallels (no intersection). Otherwise both segments
///   are colinear; check overlap via `T = Re(v_z) = dot(uA, delta)`
///   against the unit-length endpoints (`-1 < T < 1` when `uA == uB`,
///   `0 < T < 2` when `uA == -uB`).
///
/// * `sign(K) != 0`: lines cross at some point. The four sidedness
///   signs `sign(V)`, `sign(V + K)`, `sign(W)`, `sign(W - K)` classify:
///
///   - **Proper crossing**: all four nonzero, both pairs disagree.
///   - **T-touch**: exactly one sidedness sign is zero (= that endpoint
///     lies on the other segment's *line*). Intersection iff the
///     endpoint also lies strictly between the other segment's
///     endpoints, decided by the corresponding real part:
///     `Re(v_z) / Re(v_z + k_z) / Re(w_z) / Re(w_z - k_z)`, each
///     compared against 0 and +/- 1 (the unit segment's endpoints).
///   - **Same side**: both members of either pair are strictly
///     same-signed: no crossing.
///
/// # Re/Im basis tables
///
/// `im_decomp` (= K-vector of `Im(zeta^k)`) drives [`im_sign_basis`] for
/// the four sidedness signs. `re_decomp` (= K-vector of `Re(zeta^k)`)
/// drives [`project_re`] for the dot-product position checks in the
/// colinear and T-touch sub-cases. `one_in_real_basis` is the K-vector
/// of the real constant `+1`, used to compare positions against the
/// unit-segment endpoints; we accept it as a parameter rather than
/// re-deriving it from `re_decomp[0]` so the call site doesn't have to
/// know the basis layout.
#[inline]
#[allow(clippy::too_many_arguments)]
pub fn intersect_unit_segments_basis<const PHI: usize, const K: usize>(
    s1: &([i64; PHI], [i64; PHI]),
    s2: &([i64; PHI], [i64; PHI]),
    reduction: &[i64; PHI],
    conj_matrix: &[[i64; PHI]; PHI],
    re_decomp: &[[i64; K]; PHI],
    im_decomp: &[[i64; K]; PHI],
    one_in_real_basis: &[i64; K],
    real_sign: fn(&[i64; K]) -> i8,
) -> bool {
    let (a, b) = (s1.0, s1.1);
    let (c, d) = (s2.0, s2.1);

    // Touching endpoints do not count as a proper intersection.
    if a == c || a == d || b == c || b == d {
        return false;
    }

    let u_a = sub_basis::<PHI>(&b, &a);
    let u_b = sub_basis::<PHI>(&d, &c);
    let delta = sub_basis::<PHI>(&c, &a);

    let conj_u_a = conj_basis::<PHI>(&u_a, conj_matrix);

    // Mul #1: v_z = conj(uA) * delta.
    let v_z = mul_basis::<PHI>(&conj_u_a, &delta, reduction);
    let sign_v = im_sign_basis::<PHI, K>(&v_z, im_decomp, real_sign);

    // Mul #2: k_z = conj(uA) * uB.
    let k_z = mul_basis::<PHI>(&conj_u_a, &u_b, reduction);
    let sign_k = im_sign_basis::<PHI, K>(&k_z, im_decomp, real_sign);

    if sign_k == 0 {
        // uA parallel to uB. If V != 0 the lines are distinct parallels.
        if sign_v != 0 {
            return false;
        }
        // Colinear. T = Re(v_z) = dot(uA, delta) lives in the real subring;
        // its K-vector is sum_k v_z[k] * re_decomp[k].
        let t = project_re::<PHI, K>(&v_z, re_decomp);
        // k_z is +/-1 (a real unit). Decide which by checking
        // sign(k_z[0] - 1) etc, but it is cheaper to read the integer-basis
        // coeff of `1`, i.e. k_z[0], plus the others being zero -- those
        // are guaranteed when sign_k == 0 in the unit-segment setting.
        // For robustness across rings we use the projected real part of
        // k_z (a single integer comparison), which is +1 or -1 in either
        // case.
        let k_real = project_re::<PHI, K>(&k_z, re_decomp);
        let sign_k_real = real_sign(&k_real);
        if sign_k_real > 0 {
            // uA == uB: interior overlap iff -1 < T < 1.
            let t_plus_1 = add_kvec::<K>(&t, one_in_real_basis);
            let t_minus_1 = sub_kvec::<K>(&t, one_in_real_basis);
            return real_sign(&t_plus_1) > 0 && real_sign(&t_minus_1) < 0;
        } else {
            // uA == -uB: interior overlap iff 0 < T < 2.
            let two = scale_kvec::<K>(one_in_real_basis, 2);
            let t_minus_2 = sub_kvec::<K>(&t, &two);
            return real_sign(&t) > 0 && real_sign(&t_minus_2) < 0;
        }
    }

    // Lines not parallel. Process the two sidedness pairs in turn.
    // For each pair we (a) cover the T-touch sub-case when one
    // sidedness sign is zero, then (b) reject if both are strictly
    // same-signed. The remaining case -- both strictly differ -- falls
    // through to the next pair (or to the final return). See the
    // function-level docstring for the algorithm structure.

    // `t` is a real quantity in the K-vector representation. Strictly
    // interior to a unit segment iff `0 < t < 1` or, in the
    // opposite-sign convention, `-1 < t < 0`.
    let in_open_0_1 = |t: &[i64; K]| -> bool {
        let t_minus_1 = sub_kvec::<K>(t, one_in_real_basis);
        real_sign(t) > 0 && real_sign(&t_minus_1) < 0
    };
    let in_open_neg1_0 = |t: &[i64; K]| -> bool {
        let t_plus_1 = add_kvec::<K>(t, one_in_real_basis);
        real_sign(t) < 0 && real_sign(&t_plus_1) > 0
    };

    // -- Pair 1: c, d vs line ab.
    let v_plus_k = add_basis::<PHI>(&v_z, &k_z);
    let sign_v_plus_k = im_sign_basis::<PHI, K>(&v_plus_k, im_decomp, real_sign);
    // T-touch: c on line ab, position along ab = Re(v_z), in (0, 1).
    if sign_v == 0 {
        return in_open_0_1(&project_re::<PHI, K>(&v_z, re_decomp));
    }
    // T-touch: d on line ab, position along ab = Re(v_z + k_z), in (0, 1).
    if sign_v_plus_k == 0 {
        return in_open_0_1(&project_re::<PHI, K>(&v_plus_k, re_decomp));
    }
    // c, d strictly on the same side: no crossing.
    if (sign_v > 0) == (sign_v_plus_k > 0) {
        return false;
    }

    // -- Pair 2: a, b vs line cd. (Only reached if pair 1 disagrees.)
    let conj_delta = conj_basis::<PHI>(&delta, conj_matrix);
    let w_z = mul_basis::<PHI>(&conj_delta, &u_b, reduction);
    let w_minus_k = sub_basis::<PHI>(&w_z, &k_z);
    let sign_w = im_sign_basis::<PHI, K>(&w_z, im_decomp, real_sign);
    let sign_w_minus_k = im_sign_basis::<PHI, K>(&w_minus_k, im_decomp, real_sign);
    // T-touch: a on line cd. Position along cd = dot(uB, a - c)
    // = -dot(uB, delta) = -Re(w_z). In (0, 1) iff Re(w_z) in (-1, 0).
    if sign_w == 0 {
        return in_open_neg1_0(&project_re::<PHI, K>(&w_z, re_decomp));
    }
    // T-touch: b on line cd. Position along cd = dot(uB, b - c)
    // = Re(k_z) - Re(w_z) = -Re(w_z - k_z). In (0, 1) iff
    // Re(w_z - k_z) in (-1, 0).
    if sign_w_minus_k == 0 {
        return in_open_neg1_0(&project_re::<PHI, K>(&w_minus_k, re_decomp));
    }
    // Proper crossing iff a, b strictly on opposite sides of cd.
    (sign_w > 0) != (sign_w_minus_k > 0)
}

#[inline]
fn project_re<const PHI: usize, const K: usize>(
    x: &[i64; PHI],
    re_decomp: &[[i64; K]; PHI],
) -> [i64; K] {
    let mut acc = [0i64; K];
    let mut i = 0;
    while i < PHI {
        let xi = x[i];
        if xi != 0 {
            let mut j = 0;
            while j < K {
                acc[j] += xi * re_decomp[i][j];
                j += 1;
            }
        }
        i += 1;
    }
    acc
}

#[inline]
fn add_kvec<const K: usize>(x: &[i64; K], y: &[i64; K]) -> [i64; K] {
    let mut out = [0i64; K];
    let mut i = 0;
    while i < K {
        out[i] = x[i] + y[i];
        i += 1;
    }
    out
}

#[inline]
fn sub_kvec<const K: usize>(x: &[i64; K], y: &[i64; K]) -> [i64; K] {
    let mut out = [0i64; K];
    let mut i = 0;
    while i < K {
        out[i] = x[i] - y[i];
        i += 1;
    }
    out
}

#[inline]
fn scale_kvec<const K: usize>(x: &[i64; K], s: i64) -> [i64; K] {
    let mut out = [0i64; K];
    let mut i = 0;
    while i < K {
        out[i] = x[i] * s;
        i += 1;
    }
    out
}

// ----------------
// Derived per-ring tables (intended for use at startup, behind `OnceLock`).

/// Build the conjugation matrix `M` such that `conj(z) = M * z` (column-major:
/// column `j` is `conj(zeta^j)` expanded into the integer basis).
///
/// `n` is the order of `zeta`; conjugation is `conj(zeta^j) = zeta^(n - j)`.
/// We compute `zeta^k` for `k = 0..n` by iterated multiplication by `zeta`
/// (using [`mul_basis`] with the given `reduction`), then read out
/// `conj_matrix[i][j] = (zeta^(n - j))[i]`.
///
/// # Example
///
/// For the 2-element `Z[i]` ring with `n = 4`, basis `{1, i}`, reduction
/// `[-1, 0]`:
///
/// ```
/// use tilezz::cyclotomic::integral_basis::derive_conj_matrix;
/// let m = derive_conj_matrix::<2>(4, &[-1, 0]);
/// // conj(1) = 1            -> column 0 = [1, 0]
/// // conj(zeta) = zeta^3 = -i = (0, -1) (since zeta = i, zeta^2 = -1, zeta^3 = -i)
/// //   -> column 1 = [0, -1]
/// assert_eq!(m, [[1, 0], [0, -1]]);
/// ```
pub fn derive_conj_matrix<const PHI: usize>(n: usize, reduction: &[i64; PHI]) -> [[i64; PHI]; PHI] {
    let units = derive_units_lookup::<PHI>(n, reduction);
    let mut out = [[0i64; PHI]; PHI];
    let mut j = 0;
    while j < PHI {
        // conj(zeta^j) = zeta^(n - j), with zeta^0 used when j == 0.
        let idx = if j == 0 { 0 } else { n - j };
        let col = &units[idx];
        let mut i = 0;
        while i < PHI {
            out[i][j] = col[i];
            i += 1;
        }
        j += 1;
    }
    out
}

/// Compute the unit lookup table: `zeta^k` in the integer basis for
/// `k = 0..=n`. Index `0` is `1`, index `n` is `1` again (full turn).
///
/// Implementation: start at `zeta^0 = (1, 0, ..., 0)` and repeatedly multiply
/// by `zeta = (0, 1, 0, ..., 0)`, reducing each step via [`mul_basis`].
///
/// For each per-ring `OnceLock` unit table, the macro-generated impl wraps a
/// call to this helper and indexes by `angle.rem_euclid(n)`.
pub fn derive_units_lookup<const PHI: usize>(n: usize, reduction: &[i64; PHI]) -> Vec<[i64; PHI]> {
    assert!(PHI >= 1, "PHI must be >= 1");
    let mut out: Vec<[i64; PHI]> = Vec::with_capacity(n + 1);

    let mut one = [0i64; PHI];
    one[0] = 1;
    out.push(one);

    if PHI == 1 {
        // Degenerate case: ZZ_1 / ZZ_2. zeta is +/-1, captured entirely by
        // the reduction rule `reduction[0] = +/-1`. zeta^k cycles with
        // period 1 or 2.
        let mut k = 1;
        while k <= n {
            let prev = out[k - 1];
            // zeta * prev: only the zeta^0 slot got prev[0], pushed up to
            // zeta^1 which immediately reduces.
            let mut next = [0i64; PHI];
            next[0] = prev[0] * reduction[0];
            out.push(next);
            k += 1;
        }
        return out;
    }

    let mut zeta = [0i64; PHI];
    zeta[1] = 1;

    let mut k = 1;
    while k <= n {
        let prev = out[k - 1];
        let next = mul_basis::<PHI>(&prev, &zeta, reduction);
        out.push(next);
        k += 1;
    }
    out
}

// ----------------
// Macro skeleton.

/// Define an integer-basis cyclotomic ring `ZZ_n`.
///
/// Emits the full per-ring impl bag (Add/Sub/Mul/Neg, SymNum, Conj, Units,
/// ReImSign, WithinRadius, IntersectUnitSegments, ...) by routing each impl
/// to the generic `integral_basis::*` helpers parameterized on the per-ring
/// constants listed below.
///
/// # Parameters
///
/// * `name` -- ring type name (`ZZ12`, `ZZ24`, ...).
/// * `n` -- order of `zeta`.
/// * `phi` -- `phi(n)`, the storage width.
/// * `real_dim` -- `phi(n)/2`, the real-subring dimension `K`.
/// * `reduction` -- `[i64; PHI]`, expansion of `zeta^PHI` in the integer basis.
/// * `re_decomp` -- `[[i64; K]; PHI]`, K-vector of `Re(zeta^k)` per `k`.
/// * `im_decomp` -- `[[i64; K]; PHI]`, K-vector of `Im(zeta^k)` per `k`.
/// * `cartesian` -- `[Complex64; PHI]`, Cartesian value of `zeta^k` per `k`.
/// * `one_in_real_basis` -- `[i64; K]`, the K-vector for the real-subring
///   element `1` (used by the intersect-unit-segments colinear branch).
/// * `display_fn` -- `fn(&[i64; PHI], &mut std::fmt::Formatter<'_>) -> std::fmt::Result`,
///   the ring's Display impl. Defined per-ring because the symbolic shape
///   of `Re`/`Im` differs.
/// * `complex64_fn` -- `fn(&[i64; PHI]) -> num_complex::Complex64`, the
///   Cartesian projection used by `SymNum::complex64`. Defined per-ring
///   so the ring can supply a hand-rolled inline `Re/Im` expansion
///   instead of the generic `complex64_basis` loop -- `complex64` is on
///   the `cell_of` grid-bucketing hot path. Rings without a hot-path
///   need can delegate to `integral_basis::complex64_basis::<PHI>` and
///   pass the per-ring `CARTESIAN` table.
/// * `has` -- list of `HasZZk*Impl` subring-containment markers, e.g.
///   `[HasZZ4Impl, HasZZ6Impl, HasZZ12Impl]`.
#[macro_export]
macro_rules! define_integral_zz {
    (
        name: $name:ident,
        n: $n:expr,
        phi: $phi:expr,
        real_dim: $k:expr,
        reduction: $reduction:expr,
        re_decomp: $re_decomp:expr,
        im_decomp: $im_decomp:expr,
        cartesian: $cartesian:expr,
        one_in_real_basis: $one_real:expr,
        display_fn: $display_fn:path,
        complex64_fn: $complex64_fn:path,
        has: [$($has:ident),* $(,)?] $(,)?
    ) => {
        #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
        pub struct $name {
            coeffs: [i64; $phi],
        }

        impl $name {
            /// The order of `zeta` for this ring.
            pub const N: usize = $n;
            /// The storage width, `phi(n)`.
            pub const PHI: usize = $phi;
            /// The real-subring dimension, `phi(n) / 2`.
            pub const REAL_DIM: usize = $k;
            /// The reduction rule for `zeta^PHI` in the integer basis.
            pub const REDUCTION: [i64; $phi] = $reduction;
            #[doc(hidden)]
            pub const RE_DECOMP: [[i64; $k]; $phi] = $re_decomp;
            #[doc(hidden)]
            pub const IM_DECOMP: [[i64; $k]; $phi] = $im_decomp;
            #[doc(hidden)]
            pub const CARTESIAN: [num_complex::Complex64; $phi] = $cartesian;
            #[doc(hidden)]
            pub const ONE_IN_REAL_BASIS: [i64; $k] = $one_real;

            #[inline]
            pub const fn from_int_coeffs(coeffs: [i64; $phi]) -> Self {
                Self { coeffs }
            }

            /// Raw access to the integer-basis coefficients.
            #[inline]
            pub const fn int_coeffs(&self) -> [i64; $phi] {
                self.coeffs
            }
        }

        // Per-ring `OnceLock`-backed caches (conjugation matrix, units table).
        impl $name {
            #[doc(hidden)]
            #[inline]
            pub fn __conj_matrix() -> &'static [[i64; $phi]; $phi] {
                static M: std::sync::OnceLock<[[i64; $phi]; $phi]> =
                    std::sync::OnceLock::new();
                M.get_or_init(|| {
                    $crate::cyclotomic::integral_basis::derive_conj_matrix::<$phi>(
                        $n, &Self::REDUCTION,
                    )
                })
            }

            #[doc(hidden)]
            #[inline]
            pub fn __units_table() -> &'static Vec<[i64; $phi]> {
                static U: std::sync::OnceLock<Vec<[i64; $phi]>> =
                    std::sync::OnceLock::new();
                U.get_or_init(|| {
                    $crate::cyclotomic::integral_basis::derive_units_lookup::<$phi>(
                        $n, &Self::REDUCTION,
                    )
                })
            }
        }

        // ---- Display ----
        impl std::fmt::Display for $name {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                $display_fn(&self.coeffs, f)
            }
        }

        // ---- Neg / Add / Sub / Mul ----
        impl std::ops::Neg for $name {
            type Output = Self;
            #[inline]
            fn neg(self) -> Self {
                Self {
                    coeffs: $crate::cyclotomic::integral_basis::neg_basis::<$phi>(&self.coeffs),
                }
            }
        }

        impl std::ops::Add<$name> for $name {
            type Output = Self;
            #[inline]
            fn add(self, other: Self) -> Self {
                Self {
                    coeffs: $crate::cyclotomic::integral_basis::add_basis::<$phi>(
                        &self.coeffs,
                        &other.coeffs,
                    ),
                }
            }
        }

        impl std::ops::Sub<$name> for $name {
            type Output = Self;
            #[inline]
            fn sub(self, other: Self) -> Self {
                Self {
                    coeffs: $crate::cyclotomic::integral_basis::sub_basis::<$phi>(
                        &self.coeffs,
                        &other.coeffs,
                    ),
                }
            }
        }

        // NOTE: `impl Mul<$name> for $name` is NOT emitted by this macro.
        // `Mul` is the single hottest hot-path; we make every ring opt into
        // either the generic-via-basis path (via `impl_integral_mul_via_basis!`)
        // or a hand-rolled per-ring fast-path body. See `rings.rs` for ZZ12's
        // hand-rolled 4x4 Mul.

        // ---- Zero / One ----
        impl num_traits::Zero for $name {
            #[inline]
            fn zero() -> Self {
                Self { coeffs: [0i64; $phi] }
            }
            #[inline]
            fn is_zero(&self) -> bool {
                self.coeffs == [0i64; $phi]
            }
        }

        impl num_traits::One for $name {
            #[inline]
            fn one() -> Self {
                let mut c = [0i64; $phi];
                c[0] = 1;
                Self { coeffs: c }
            }
            #[inline]
            fn is_one(&self) -> bool {
                if self.coeffs[0] != 1 {
                    return false;
                }
                let mut i = 1;
                while i < $phi {
                    if self.coeffs[i] != 0 {
                        return false;
                    }
                    i += 1;
                }
                true
            }
        }

        impl $crate::cyclotomic::IntRing for $name {}

        // ---- Pow ----
        impl num_traits::Pow<u8> for $name {
            type Output = Self;
            fn pow(self, other: u8) -> Self {
                use $crate::cyclotomic::SymNum;
                self.zz_pow(other)
            }
        }
        impl num_traits::Pow<i8> for $name {
            type Output = Self;
            fn pow(self, other: i8) -> Self {
                use $crate::cyclotomic::SymNum;
                assert!(other >= 0, "Negative powers are not supported!");
                self.zz_pow(other as u8)
            }
        }

        // ---- InnerIntType + From<i64> + From<(i64, i64)> if HasZZ4 ----
        impl $crate::cyclotomic::numtraits::InnerIntType for $name {
            type IntType = i64;
        }

        impl From<i64> for $name {
            #[inline]
            fn from(value: i64) -> Self {
                let mut c = [0i64; $phi];
                c[0] = value;
                Self { coeffs: c }
            }
        }

        // ---- SymNum ----
        impl $crate::cyclotomic::SymNum for $name {
            type Scalar = i64;

            #[inline]
            fn zz_coeffs(&self) -> &[Self::Scalar] {
                &self.coeffs
            }

            #[inline]
            fn zz_coeffs_mut(&mut self) -> &mut [Self::Scalar] {
                &mut self.coeffs
            }

            #[inline]
            fn turn() -> i8 {
                $n as i8
            }

            #[inline]
            fn complex64(&self) -> num_complex::Complex64 {
                $complex64_fn(&self.coeffs)
            }
        }

        // NOTE: `impl Conj for $name` is NOT emitted by this macro.
        // Conj is on the hot intersect_unit_segments path; the macro lets
        // each ring choose between the generic-via-basis path (via
        // `impl_integral_conj_via_basis!`) and a hand-rolled per-ring inline
        // body. See `rings.rs` for ZZ12's hand-rolled conj
        // (`(a + c, b, -c, -b - d)`).

        // ---- Ccw ----
        impl $crate::cyclotomic::Ccw for $name {
            #[inline]
            fn ccw() -> Self {
                <Self as $crate::cyclotomic::Units>::unit(1)
            }
            #[inline]
            fn is_ccw(&self) -> bool {
                *self == <Self as $crate::cyclotomic::Ccw>::ccw()
            }
        }

        // ---- ZZComplex ----
        impl $crate::cyclotomic::ZZComplex for $name {
            fn is_real(&self) -> bool {
                <Self as $crate::cyclotomic::ReImSign>::im_sign(self) == 0
            }
            fn is_imag(&self) -> bool {
                <Self as $crate::cyclotomic::ReImSign>::re_sign(self) == 0
            }
        }

        // ---- IntCoeffsSlice ----
        // Borrowed-slice view of the integer-basis coefficients. Useful
        // for generic code that needs the raw coords (e.g. modular
        // pruning) without committing to a per-ring fixed PHI.
        impl $crate::cyclotomic::integral_basis::IntCoeffsSlice for $name {
            #[inline]
            fn int_coeffs_slice(&self) -> &[i64] {
                &self.coeffs
            }
        }

        // NOTE: `impl Units for $name` is NOT emitted by this macro.
        // `Units::unit` shows up in `rat_enum`'s tight loop; per-ring
        // hand-rolled `static UNIT_TABLE: [[i64; PHI]; N] = [...];` is
        // faster than the OnceLock + Vec indirection. The default macro
        // `impl_integral_units_via_basis!` is available for rings that
        // don't yet have a hand-rolled table.

        // ---- Catch-all "is a cyclotomic ring" marker impl ----
        impl $crate::cyclotomic::IsRing for $name {}

        // ---- Subring-containment markers ----
        $(
            impl $crate::cyclotomic::traits::$has for $name {}
        )*
    };
}

/// Default `Units` impl for an integer-basis ring: cache `derive_units_lookup`
/// in a `OnceLock<Vec<[i64; PHI]>>` and index by `angle.rem_euclid(n)`.
///
/// Rings whose `unit(k)` is called in a tight loop can provide their own
/// `impl Units for $name` with a `static UNIT_TABLE: [[i64; PHI]; n] = [...]`
/// instead, which avoids the OnceLock check and the heap-resident Vec.
#[macro_export]
macro_rules! impl_integral_units_via_basis {
    ($name:ident, $n:expr) => {
        impl $crate::cyclotomic::Units for $name {
            #[inline]
            fn unit(angle: i8) -> Self {
                let tab = <$name>::__units_table();
                let idx = angle.rem_euclid($n as i8) as usize;
                Self::from_int_coeffs(tab[idx])
            }
        }
    };
}

/// Default `Mul<Self>` impl for an integer-basis ring: route through
/// [`mul_basis`] with the per-ring `REDUCTION` table.
///
/// Rings that have a faster hand-rolled `Mul` can opt out of this macro and
/// provide their own `impl Mul<$name> for $name { ... }` body (e.g. ZZ12
/// uses the inline 4x4 polynomial expansion in `rings.rs`).
#[macro_export]
macro_rules! impl_integral_mul_via_basis {
    ($name:ident, $phi:expr) => {
        impl std::ops::Mul<$name> for $name {
            type Output = Self;
            #[inline]
            fn mul(self, other: Self) -> Self {
                Self::from_int_coeffs($crate::cyclotomic::integral_basis::mul_basis::<$phi>(
                    &self.int_coeffs(),
                    &other.int_coeffs(),
                    &<$name>::REDUCTION,
                ))
            }
        }
    };
}

/// Default `Conj` impl for an integer-basis ring: route through
/// [`conj_basis`] with the cached per-ring conjugation matrix.
///
/// Rings that have a cheaper inline conjugation (e.g. ZZ12's
/// `(a + c, b, -c, -b - d)`) can opt out of this macro and provide their
/// own `impl Conj for $name`.
#[macro_export]
macro_rules! impl_integral_conj_via_basis {
    ($name:ident, $phi:expr) => {
        impl $crate::cyclotomic::Conj for $name {
            #[inline]
            fn conj(&self) -> Self {
                Self::from_int_coeffs($crate::cyclotomic::integral_basis::conj_basis::<$phi>(
                    &self.int_coeffs(),
                    <$name>::__conj_matrix(),
                ))
            }
        }
    };
}

/// Default `ReImSign` impl for an integer-basis ring: route through
/// `re_sign_basis` / `im_sign_basis` against the per-ring `RE_DECOMP` /
/// `IM_DECOMP` tables and the supplied `real_sign_fn`.
///
/// Rings that have a cheaper hand-rolled sign extraction (e.g. ZZ12 with
/// `sign_m_plus_n_sqrt3` reading the (m, n) components directly) can opt out
/// of this macro and provide their own `impl ReImSign for $name { ... }`.
#[macro_export]
macro_rules! impl_integral_re_im_sign_via_basis {
    ($name:ident, $phi:expr, $k:expr, $real_sign:path) => {
        impl $crate::cyclotomic::ReImSign for $name {
            #[inline]
            fn re_sign(&self) -> i8 {
                $crate::cyclotomic::integral_basis::re_sign_basis::<$phi, $k>(
                    &<$name>::int_coeffs(self),
                    &<$name>::RE_DECOMP,
                    $real_sign,
                )
            }
            #[inline]
            fn im_sign(&self) -> i8 {
                $crate::cyclotomic::integral_basis::im_sign_basis::<$phi, $k>(
                    &<$name>::int_coeffs(self),
                    &<$name>::IM_DECOMP,
                    $real_sign,
                )
            }
        }
    };
}

/// Default `IntersectUnitSegments` impl for an integer-basis ring: route
/// through `intersect_unit_segments_basis` with the per-ring constant bundle.
///
/// Rings that have a faster hand-rolled fast path (e.g. ZZ12's inline 3-mul
/// path that exploits the `[i64; 4]` layout and ZZ12-specific
/// `sign_m_plus_n_sqrt3`) can opt out of this macro and provide their own
/// `impl IntersectUnitSegments for $name { ... }`.
#[macro_export]
macro_rules! impl_integral_intersect_unit_segments_via_basis {
    ($name:ident, $phi:expr, $k:expr, $real_sign:path) => {
        impl $crate::cyclotomic::IntersectUnitSegments for $name {
            #[inline]
            fn intersect_unit_segments(s1: &($name, $name), s2: &($name, $name)) -> bool {
                $crate::cyclotomic::integral_basis::intersect_unit_segments_basis::<$phi, $k>(
                    &(s1.0.int_coeffs(), s1.1.int_coeffs()),
                    &(s2.0.int_coeffs(), s2.1.int_coeffs()),
                    &<$name>::REDUCTION,
                    <$name>::__conj_matrix(),
                    &<$name>::RE_DECOMP,
                    &<$name>::IM_DECOMP,
                    &<$name>::ONE_IN_REAL_BASIS,
                    $real_sign,
                )
            }
        }
    };
}

/// Default `WithinRadius` impl for an integer-basis ring: compare
/// `|z|^2` to `radius^2` exactly via ring arithmetic plus one `re_sign`
/// call. No f64.
///
/// `|z|^2 = z * conj(z)` is real (lives in the real subring), so
/// `(radius^2 - z*conj(z)).re_sign() >= 0` is the exact decision.
///
/// Rings with a tighter pure-i64 fast path (e.g. ZZ12 inlining its
/// `{1, sqrt(3)}` norm-squared formula) can opt out of this macro and
/// provide their own `impl WithinRadius for $name`.
#[macro_export]
macro_rules! impl_integral_within_radius_via_norm_sq {
    ($name:ident) => {
        impl $crate::cyclotomic::WithinRadius for $name {
            #[inline]
            fn within_radius(&self, radius: i64) -> bool {
                use $crate::cyclotomic::{Conj, ReImSign};
                let norm_sq = *self * self.conj();
                let r_sq = <$name as From<i64>>::from(radius * radius);
                (r_sq - norm_sq).re_sign() >= 0
            }
        }
    };
}

// ----------------
// Generic test suite for integer-basis cyclotomic rings.
//
// `zz_integral_ring_tests!(name: SomeRing)` emits a `mod` of ring-axiom and
// cyclotomic-structure tests for the given ring type. Invoked once per
// ring in `rings.rs`.
//
// The suite covers:
//   * Algebra: associativity / commutativity / distributivity / identity /
//     inverse for `+` and `*`, plus `-`, plus zero-absorbing.
//   * Cyclotomic structure: `unit(0) = 1`, `unit(N/2) = -1`, the unit cycle
//     group property `unit(i)*unit(j) = unit((i+j) mod N)`, `unit(1)^k = unit(k)`.
//   * Conjugation: involution, additive / multiplicative homomorphism,
//     reality of `x * conj(x)`, conj fixes `1`.
//   * Cartesian projection: `complex64` is additive / multiplicative within
//     epsilon, `unit(k).complex64() ~= exp(2*pi*i*k/N)`.
//   * `ReImSign`: agreement with `complex64()` for non-boundary values.
//   * `WithinRadius`: agreement with `|z| <= r`.
//   * `IntersectUnitSegments`: matches the generic `intersect` and is
//     symmetric in its arguments.

/// Emit a `mod` of generic algebraic / cyclotomic-structure tests for the
/// given integer-basis cyclotomic ring type.
///
/// # Shape
///
/// ```ignore
/// zz_integral_ring_tests!(name: ZZ12);
/// ```
///
/// emits roughly:
///
/// ```ignore
/// #[cfg(test)]
/// mod zz12_generic_ring_tests {
///     // ... ~25 #[test] functions for algebra, conjugation, complex64
///     // projection, signs, within-radius, intersect-unit-segments ...
/// }
/// ```
///
/// The macro is hygienic across multiple instantiations -- the generated
/// `mod` is named `<ring_lower>_generic_ring_tests` (e.g. `zz12_generic_ring_tests`).
///
/// # Requirements on `name`
///
/// The ring type must implement the usual cyclotomic-ring trait bundle:
/// `SymNum`, `Units`, `Ccw`, `Conj`, `ReImSign`, `WithinRadius`,
/// `IntersectUnitSegments`, plus `Add`/`Sub`/`Mul`/`Neg`, `Zero`/`One`,
/// and `From<i64>`.
#[macro_export]
macro_rules! zz_integral_ring_tests {
    (name: $ring:ident) => {
        ::paste::paste! {
            #[cfg(test)]
            mod [<$ring:lower _generic_ring_tests>] {
                #[allow(unused_imports)]
                use super::*;
                use $crate::cyclotomic::{
                    Ccw, Conj, IntersectUnitSegments, ReImSign, SymNum, Units, WithinRadius,
                };
                use $crate::cyclotomic::geometry::intersect;
                use num_traits::{One, Zero};

                type R = $ring;

                /// Small grid of "interesting" ring elements: zero, the
                /// integer constants {-2, -1, 1, 2}, every `unit(k)`, and a
                /// handful of integer combinations of small units.
                fn sample_elems() -> Vec<R> {
                    let n = R::turn();
                    let mut out: Vec<R> = vec![
                        R::zero(),
                        R::one(),
                        -R::one(),
                        R::from(2i64),
                        R::from(-2i64),
                    ];
                    let mut k = 0i8;
                    while k < n {
                        out.push(<R as Units>::unit(k));
                        k += 1;
                    }
                    // A few mixed elements.
                    out.push(<R as Units>::unit(1) + <R as Units>::unit(2));
                    out.push(<R as Units>::unit(0) - <R as Units>::unit(3));
                    out.push(<R as Units>::unit(1) + R::from(3i64));
                    out.push(R::from(2i64) * <R as Units>::unit(1) - <R as Units>::unit(2));
                    out
                }

                /// Smaller sample for O(n^3) tests (associativity etc.).
                fn small_sample_elems() -> Vec<R> {
                    vec![
                        R::zero(),
                        R::one(),
                        -R::one(),
                        <R as Units>::unit(1),
                        <R as Units>::unit(2),
                        <R as Units>::unit(1) + <R as Units>::unit(2),
                        R::from(2i64) - <R as Units>::unit(1),
                    ]
                }

                // ---- Algebra ----

                #[test]
                fn add_identity() {
                    let zero = R::zero();
                    for x in sample_elems() {
                        assert_eq!(x + zero, x, "x + 0 != x for x = {x}");
                        assert_eq!(zero + x, x, "0 + x != x for x = {x}");
                    }
                }

                #[test]
                fn add_inverse() {
                    let zero = R::zero();
                    for x in sample_elems() {
                        assert_eq!(x + (-x), zero, "x + (-x) != 0 for x = {x}");
                    }
                }

                #[test]
                fn add_commutativity() {
                    let xs = sample_elems();
                    for x in &xs {
                        for y in &xs {
                            assert_eq!(*x + *y, *y + *x, "x + y != y + x for x={x}, y={y}");
                        }
                    }
                }

                #[test]
                fn add_associativity() {
                    let xs = small_sample_elems();
                    for x in &xs {
                        for y in &xs {
                            for z in &xs {
                                assert_eq!(
                                    (*x + *y) + *z,
                                    *x + (*y + *z),
                                    "associativity of + failed at x={x}, y={y}, z={z}",
                                );
                            }
                        }
                    }
                }

                #[test]
                fn mul_identity() {
                    let one = R::one();
                    for x in sample_elems() {
                        assert_eq!(x * one, x, "x * 1 != x for x = {x}");
                        assert_eq!(one * x, x, "1 * x != x for x = {x}");
                    }
                }

                #[test]
                fn mul_commutativity() {
                    let xs = sample_elems();
                    for x in &xs {
                        for y in &xs {
                            assert_eq!(*x * *y, *y * *x, "x * y != y * x for x={x}, y={y}");
                        }
                    }
                }

                #[test]
                fn mul_associativity() {
                    let xs = small_sample_elems();
                    for x in &xs {
                        for y in &xs {
                            for z in &xs {
                                assert_eq!(
                                    (*x * *y) * *z,
                                    *x * (*y * *z),
                                    "associativity of * failed at x={x}, y={y}, z={z}",
                                );
                            }
                        }
                    }
                }

                #[test]
                fn distributivity() {
                    let xs = small_sample_elems();
                    for x in &xs {
                        for y in &xs {
                            for z in &xs {
                                assert_eq!(
                                    *x * (*y + *z),
                                    *x * *y + *x * *z,
                                    "x*(y+z) != x*y + x*z for x={x}, y={y}, z={z}",
                                );
                            }
                        }
                    }
                }

                #[test]
                fn zero_absorbing() {
                    let zero = R::zero();
                    for x in sample_elems() {
                        assert_eq!(x * zero, zero, "x * 0 != 0 for x = {x}");
                        assert_eq!(zero * x, zero, "0 * x != 0 for x = {x}");
                    }
                }

                #[test]
                fn sub_eq_add_neg() {
                    let xs = sample_elems();
                    for x in &xs {
                        for y in &xs {
                            assert_eq!(
                                *x - *y,
                                *x + (-(*y)),
                                "x - y != x + (-y) for x={x}, y={y}",
                            );
                        }
                    }
                }

                // ---- Cyclotomic structure ----

                #[test]
                fn unit_zero_is_one() {
                    assert_eq!(<R as Units>::unit(0), R::one());
                    let n = R::turn();
                    assert_eq!(<R as Units>::unit(n), R::one());
                    if n % 2 == 0 {
                        assert_eq!(<R as Units>::unit(n / 2), -R::one());
                    }
                }

                #[test]
                fn unit_cycle_group() {
                    let n = R::turn();
                    for i in 0..n {
                        for j in 0..n {
                            let lhs = <R as Units>::unit(i) * <R as Units>::unit(j);
                            let rhs = <R as Units>::unit(
                                ((i as i16 + j as i16).rem_euclid(n as i16)) as i8,
                            );
                            assert_eq!(lhs, rhs, "unit({i})*unit({j}) != unit((i+j) mod N)");
                        }
                    }
                }

                #[test]
                fn powers_of_zeta() {
                    let n = R::turn();
                    let zeta = R::ccw();
                    let mut acc = R::one();
                    let mut k = 0i8;
                    while k < n {
                        assert_eq!(
                            acc,
                            <R as Units>::unit(k),
                            "unit(1)^{k} != unit({k})",
                        );
                        acc = acc * zeta;
                        k += 1;
                    }
                    // After a full turn, we're back at one.
                    assert_eq!(acc, R::one(), "unit(1)^N != 1");
                }

                // ---- Conjugation ----

                #[test]
                fn conj_involution() {
                    for x in sample_elems() {
                        assert_eq!(x.conj().conj(), x, "conj(conj(x)) != x for x = {x}");
                    }
                }

                #[test]
                fn conj_distributive_over_add() {
                    let xs = sample_elems();
                    for x in &xs {
                        for y in &xs {
                            assert_eq!(
                                (*x + *y).conj(),
                                x.conj() + y.conj(),
                                "conj(x+y) != conj(x)+conj(y) for x={x}, y={y}",
                            );
                        }
                    }
                }

                #[test]
                fn conj_distributive_over_mul() {
                    let xs = small_sample_elems();
                    for x in &xs {
                        for y in &xs {
                            assert_eq!(
                                (*x * *y).conj(),
                                x.conj() * y.conj(),
                                "conj(x*y) != conj(x)*conj(y) for x={x}, y={y}",
                            );
                        }
                    }
                }

                #[test]
                fn conj_fixes_one() {
                    let one = <R as Units>::unit(0);
                    assert_eq!(one.conj(), one);
                }

                #[test]
                fn x_times_conj_x_is_real() {
                    for x in sample_elems() {
                        let prod = x * x.conj();
                        assert_eq!(
                            prod.im_sign(),
                            0,
                            "x*conj(x) is not real for x = {x}: prod = {prod}",
                        );
                    }
                }

                // ---- Cartesian projection (complex64) ----

                const COMPLEX_EPS: f64 = 1e-9;

                fn approx_eq(a: num_complex::Complex64, b: num_complex::Complex64, eps: f64) -> bool {
                    (a.re - b.re).abs() < eps && (a.im - b.im).abs() < eps
                }

                #[test]
                fn complex64_linear() {
                    let xs = sample_elems();
                    for x in &xs {
                        for y in &xs {
                            let lhs = (*x + *y).complex64();
                            let rhs = x.complex64() + y.complex64();
                            assert!(
                                approx_eq(lhs, rhs, COMPLEX_EPS),
                                "(x+y).complex64() != x.complex64() + y.complex64() for x={x}, y={y}: {lhs} vs {rhs}",
                            );
                        }
                    }
                }

                #[test]
                fn complex64_mul_matches() {
                    let xs = small_sample_elems();
                    for x in &xs {
                        for y in &xs {
                            let lhs = (*x * *y).complex64();
                            let rhs = x.complex64() * y.complex64();
                            // Multiplicative slack: integer combinations of small units.
                            assert!(
                                approx_eq(lhs, rhs, 1e-8),
                                "(x*y).complex64() != x.complex64() * y.complex64() for x={x}, y={y}: {lhs} vs {rhs}",
                            );
                        }
                    }
                }

                #[test]
                fn complex64_unit_matches_exp() {
                    let n = R::turn();
                    for k in 0..n {
                        let u = <R as Units>::unit(k);
                        let c = u.complex64();
                        let angle = 2.0 * std::f64::consts::PI * (k as f64) / (n as f64);
                        let expected = num_complex::Complex64::new(angle.cos(), angle.sin());
                        assert!(
                            approx_eq(c, expected, COMPLEX_EPS),
                            "unit({k}).complex64() != exp(2*pi*i*{k}/{n}): {c} vs {expected}",
                        );
                    }
                }

                // ---- ReImSign ----

                #[test]
                fn re_im_sign_agrees_with_complex64() {
                    for x in sample_elems() {
                        let c = x.complex64();
                        let re_sign = x.re_sign();
                        let im_sign = x.im_sign();
                        if c.re.abs() > 1e-9 {
                            let expected = if c.re > 0.0 { 1 } else { -1 };
                            assert_eq!(
                                re_sign, expected,
                                "re_sign disagrees with complex64 for x = {x}: re={}, sign={re_sign}",
                                c.re,
                            );
                        }
                        if c.im.abs() > 1e-9 {
                            let expected = if c.im > 0.0 { 1 } else { -1 };
                            assert_eq!(
                                im_sign, expected,
                                "im_sign disagrees with complex64 for x = {x}: im={}, sign={im_sign}",
                                c.im,
                            );
                        }
                    }
                }

                // ---- WithinRadius ----

                #[test]
                fn within_radius_agrees_with_norm() {
                    let xs = sample_elems();
                    let radii: [i64; 5] = [0, 1, 2, 3, 5];
                    for x in &xs {
                        let n = x.complex64().norm();
                        for &r in &radii {
                            let expected_inside = n <= (r as f64) + 1e-9;
                            let expected_outside = n >= (r as f64) + 1e-9;
                            let actual = x.within_radius(r);
                            // Skip near-boundary cases where float epsilon
                            // could legitimately decide either way.
                            if (n - r as f64).abs() < 1e-6 {
                                continue;
                            }
                            if expected_inside {
                                assert!(
                                    actual,
                                    "within_radius({r}) should be true for x = {x} (|x| = {n})",
                                );
                            }
                            if expected_outside && n > (r as f64) + 1e-6 {
                                assert!(
                                    !actual,
                                    "within_radius({r}) should be false for x = {x} (|x| = {n})",
                                );
                            }
                        }
                    }
                }

                // ---- IntersectUnitSegments ----

                #[test]
                fn intersect_unit_segments_matches_generic() {
                    let n = R::turn();
                    let zero = R::zero();
                    // Anchor points to translate the second segment by.
                    let anchors: Vec<R> = vec![
                        zero,
                        <R as Units>::unit(0),
                        <R as Units>::unit(1),
                        <R as Units>::unit(2),
                        -<R as Units>::unit(0),
                    ];
                    // For each (a, b) pair with a < b, the segment is
                    // (origin, origin + unit(a)) vs (p + unit(0), p + unit(0) + unit(b)).
                    // We use a sparse grid: every 2 angles, to keep test runtime
                    // moderate for rings with large N.
                    let step = (n / 6).max(1);
                    for ai in (0..n).step_by(step as usize) {
                        for bi in (0..n).step_by(step as usize) {
                            for p in &anchors {
                                let s1 = (zero, <R as Units>::unit(ai));
                                let s2 = (
                                    *p,
                                    *p + <R as Units>::unit(bi),
                                );
                                let specialized = <R as IntersectUnitSegments>::intersect_unit_segments(&s1, &s2);
                                let generic = intersect::<R>(&s1, &s2);
                                assert_eq!(
                                    specialized, generic,
                                    "intersect_unit_segments disagrees with generic intersect: s1=({}, {}), s2=({}, {}); specialized={specialized}, generic={generic}",
                                    s1.0, s1.1, s2.0, s2.1,
                                );
                            }
                        }
                    }
                }

                #[test]
                fn intersect_unit_segments_symmetric() {
                    let n = R::turn();
                    let zero = R::zero();
                    let anchors: Vec<R> = vec![
                        zero,
                        <R as Units>::unit(0),
                        <R as Units>::unit(1),
                        <R as Units>::unit(2),
                    ];
                    let step = (n / 6).max(1);
                    for ai in (0..n).step_by(step as usize) {
                        for bi in (0..n).step_by(step as usize) {
                            for p in &anchors {
                                let s1 = (zero, <R as Units>::unit(ai));
                                let s2 = (
                                    *p,
                                    *p + <R as Units>::unit(bi),
                                );
                                let ab = <R as IntersectUnitSegments>::intersect_unit_segments(&s1, &s2);
                                let ba = <R as IntersectUnitSegments>::intersect_unit_segments(&s2, &s1);
                                assert_eq!(ab, ba, "intersect_unit_segments not symmetric: s1=({}, {}), s2=({}, {})", s1.0, s1.1, s2.0, s2.1);
                            }
                        }
                    }
                }

                // ---- ZZComplex predicates ----

                #[test]
                fn is_real_imag_complex_predicates() {
                    use $crate::cyclotomic::ZZComplex;

                    assert!(R::zero().is_real(), "0 should be real");
                    assert!(R::zero().is_imag(), "0 should be imag");
                    assert!(!R::zero().is_complex(), "0 should not be (strictly) complex");

                    assert!(R::one().is_real(), "1 should be real");
                    assert!((-R::one()).is_real(), "-1 should be real");
                    assert!(!R::one().is_imag(), "1 should not be imag");
                    assert!(!R::one().is_complex(), "1 should not be (strictly) complex");

                    let n = R::turn();
                    let qturn = n / 4;
                    if n % 4 == 0 {
                        // i = unit(qturn) for any ring containing i.
                        let i = <R as Units>::unit(qturn);
                        assert!(!i.is_real(), "i should not be real");
                        assert!(i.is_imag(), "i should be imag");
                        assert!(!i.is_complex(), "i should not be (strictly) complex");
                        assert!((-i).is_imag(), "-i should be imag");
                    }

                    // unit(1) is a non-axis root of unity for N >= 8.
                    if n >= 8 && n % 4 == 0 {
                        let p = <R as Units>::unit(1);
                        assert!(!p.is_real(), "unit(1) should not be real for N={n}");
                        assert!(!p.is_imag(), "unit(1) should not be imag for N={n}");
                        assert!(p.is_complex(), "unit(1) should be strictly complex for N={n}");
                    }
                }

                // ---- Pow ----

                #[test]
                fn pow_at_hturn_is_neg_one() {
                    use num_traits::Pow;
                    let n = R::turn();
                    let hturn = (n / 2) as u8;
                    let zeta = R::ccw();
                    assert_eq!(zeta.pow(hturn), -R::one(), "zeta^(N/2) should equal -1");
                    let turn = n as u8;
                    assert_eq!(zeta.pow(turn), R::one(), "zeta^N should equal 1");
                }

                #[test]
                #[should_panic]
                fn pow_negative_panics() {
                    use num_traits::Pow;
                    let _ = R::one().pow(-1i8);
                }

                // ---- Sample complex64 values ----

                #[test]
                fn complex64_sample_values() {
                    use num_complex::Complex64;
                    use num_traits::{One, Zero};

                    assert_eq!(R::zero().complex64(), Complex64::zero());
                    assert_eq!(R::one().complex64(), Complex64::one());
                    assert_eq!((-R::one()).complex64(), -Complex64::one());
                    let two = R::one() + R::one();
                    assert!(
                        approx_eq(two.complex64(), Complex64::new(2.0, 0.0), COMPLEX_EPS),
                        "(1+1).complex64() != 2",
                    );
                }

                // ---- Hashable round-trip ----

                #[test]
                fn hashable_round_trip() {
                    use std::collections::HashSet;
                    let mut s: HashSet<R> = HashSet::new();
                    s.insert(R::zero());
                    s.insert(R::one());
                    s.insert(R::ccw());
                    assert!(s.contains(&R::ccw()));
                    assert!(s.contains(&(R::ccw() + R::zero())));
                    assert!(!s.contains(&(R::ccw() + R::one())));
                }
            }
        }
    };
}

// ----------------
// Tests

#[cfg(test)]
mod tests {
    use super::*;

    // A tiny "ZZ_4 as Z[i]" smoke test: basis {1, i}, reduction i^2 = -1.
    // Not a real ring of interest, just a fast sanity check of the engine.
    const RED_I: [i64; 2] = [-1, 0];

    #[test]
    fn add_sub_neg() {
        let x = [3i64, 4];
        let y = [1i64, 2];
        assert_eq!(add_basis::<2>(&x, &y), [4, 6]);
        assert_eq!(sub_basis::<2>(&x, &y), [2, 2]);
        assert_eq!(neg_basis::<2>(&x), [-3, -4]);
    }

    #[test]
    fn mul_smoke_z_i() {
        // (a + bi)(c + di) = (ac - bd) + (ad + bc)i
        let x = [3i64, 4]; // 3 + 4i
        let y = [1i64, 2]; // 1 + 2i
                           // (3 + 4i)(1 + 2i) = 3 + 6i + 4i + 8i^2 = 3 + 10i - 8 = -5 + 10i
        assert_eq!(mul_basis::<2>(&x, &y, &RED_I), [-5, 10]);

        // i * i = -1
        let i = [0i64, 1];
        assert_eq!(mul_basis::<2>(&i, &i, &RED_I), [-1, 0]);
    }

    #[test]
    fn derive_conj_matrix_z_i() {
        // n = 4, basis {1, i}, conj(1) = 1, conj(i) = -i.
        let m = derive_conj_matrix::<2>(4, &RED_I);
        assert_eq!(m, [[1, 0], [0, -1]]);

        // Apply: conj(3 + 4i) = 3 - 4i.
        let z = [3i64, 4];
        assert_eq!(conj_basis::<2>(&z, &m), [3, -4]);
    }

    #[test]
    fn derive_units_lookup_z_i() {
        // n = 4: zeta^0 = 1, zeta^1 = i, zeta^2 = -1, zeta^3 = -i, zeta^4 = 1.
        let u = derive_units_lookup::<2>(4, &RED_I);
        assert_eq!(u.len(), 5);
        assert_eq!(u[0], [1, 0]);
        assert_eq!(u[1], [0, 1]);
        assert_eq!(u[2], [-1, 0]);
        assert_eq!(u[3], [0, -1]);
        assert_eq!(u[4], [1, 0]);
    }

    #[test]
    fn mul_smoke_zz12_pattern() {
        // Use ZZ12's reduction rule: Phi_12 = x^4 - x^2 + 1, so
        // zeta^4 = zeta^2 - 1, i.e. reduction = [-1, 0, 1, 0].
        const RED_ZZ12: [i64; 4] = [-1, 0, 1, 0];

        // zeta * zeta = zeta^2.
        let zeta = [0i64, 1, 0, 0];
        assert_eq!(mul_basis::<4>(&zeta, &zeta, &RED_ZZ12), [0, 0, 1, 0]);

        // zeta^2 * zeta^2 = zeta^4 = zeta^2 - 1.
        let zeta_sq = [0i64, 0, 1, 0];
        assert_eq!(mul_basis::<4>(&zeta_sq, &zeta_sq, &RED_ZZ12), [-1, 0, 1, 0],);

        // zeta^3 * zeta^3 = zeta^6 = zeta * zeta^5 = ... = -1.
        let zeta_cu = [0i64, 0, 0, 1];
        assert_eq!(mul_basis::<4>(&zeta_cu, &zeta_cu, &RED_ZZ12), [-1, 0, 0, 0]);
    }

    #[test]
    fn re_im_sign_z_i() {
        // For Z[i]: Re(1) = 1, Re(i) = 0, Im(1) = 0, Im(i) = 1.
        // Real subring is just Z, so K = 1 and real_sign is plain sign of
        // the lone integer coefficient.
        const RE: [[i64; 1]; 2] = [[1], [0]];
        const IM: [[i64; 1]; 2] = [[0], [1]];
        fn sgn(x: &[i64; 1]) -> i8 {
            x[0].signum() as i8
        }
        let z = [3i64, -4];
        assert_eq!(re_sign_basis::<2, 1>(&z, &RE, sgn), 1);
        assert_eq!(im_sign_basis::<2, 1>(&z, &IM, sgn), -1);
    }

    #[test]
    fn complex64_z_i() {
        use num_complex::Complex64;
        let cart = [Complex64::new(1.0, 0.0), Complex64::new(0.0, 1.0)];
        let z = [3i64, -4];
        let c = complex64_basis::<2>(&z, &cart);
        assert!((c.re - 3.0).abs() < 1e-12);
        assert!((c.im + 4.0).abs() < 1e-12);
    }
}
