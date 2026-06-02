//! Canonical-form predicates and finalizers for the rat enumeration.
//!
//! The DFS produces every cyclic rotation of every closed polygon by
//! default; canonical-form bookkeeping collapses each equivalence
//! class to one representative. Two equivalence flavors are
//! supported:
//!
//! - **Rotation-canonical**: one rep per cyclic rotation class.
//!   Prefix prune is sound + complete (lex-min check on the walk
//!   prefix is monotone: a rotation that already lost can never
//!   recover). `canonicalize` is the identity.
//!
//! - **Dihedral-canonical**: one rep per cyclic + reflection class.
//!   Prefix prune is sound but **incomplete** -- reflection involves
//!   the unknown tail of the walk, so we can only check against
//!   complement (negated) rotations at the prefix level. The closure-
//!   time mapper [`dihedral_canonical`] picks up the slack by mapping
//!   every chirality-normalized rotation to its lex-min dihedral
//!   form, so both members of every chiral pair hash to the same
//!   bucket.
//!
//! # Why the dihedral prefix prune is heuristic
//!
//! For a closed walk `c_1..c_N`, the dihedral group has 2N elements
//! (N rotations + N reflections). For a *prefix* `c_1..c_d`:
//!
//! * Rotation k's prefix is `c_{k+1}..c_{k+d}` -- determined by
//!   `c_1..c_d` for positions where `k+i < d`. The comparable region
//!   is non-empty for every k > 0, and a verdict there is permanent
//!   (no future angles can flip it), so the prefix prune is exact.
//! * Reflection's prefix is `-c_N, -c_{N-1}, ..., -c_{N-d+1}` --
//!   every position lies in the unknown future of the current walk.
//!   There is no comparable region until N is known.
//!
//! The closure-time `dihedral_canonical` mapping handles what the
//! prefix prune can't see: it maps every closed canonical rotation
//! to the lex-min over rotations AND reversed rotations, so both
//! members of every chiral pair collapse to one rep in the result
//! HashSet.

/// Pair of canonical-check and output-mapping functions that
/// parameterise the DFS. Both rotation-canonical and dihedral-
/// canonical enumeration share the same core walk; they differ only
/// in which pair of functions they supply.
///
/// * `is_canonical` -- prefix prune applied *before* `Snake::add`.
///   Returns `false` when the extended walk cannot be the lex-min
///   rotation (or dihedral image) of any closure it could grow into.
/// * `canonicalize` -- applied to the chirality-normalised canonical
///   rotation at closure, producing the key inserted into the result
///   `HashSet`. For rotation-canonical this is the identity; for
///   dihedral-canonical it picks the lex-min over rotations and
///   reversed-rotations.
#[derive(Clone, Copy)]
pub struct CanonicalOps {
    pub is_canonical: fn(&[i8], i8) -> bool,
    pub canonicalize: fn(&[i8]) -> Vec<i8>,
}

/// Bind `prefix ++ [new]` so that index `i` (for `i < prefix.len() + 1`)
/// resolves to the corresponding walk angle. Out-of-range indices are
/// not handled -- callers must keep accesses within `0..d` where
/// `d = prefix.len() + 1`.
fn walk_get(prefix: &[i8], new: i8, i: usize) -> i8 {
    if i < prefix.len() {
        prefix[i]
    } else {
        // Caller must guarantee i == prefix.len(). debug_assert documents
        // the boundary; release-mode behavior is still correct (returns
        // `new`) but a violation indicates a bug in the caller.
        debug_assert_eq!(i, prefix.len(), "walk_get: index past extended prefix");
        new
    }
}

/// For the extended walk `prefix ++ [new]` of length `d`, return
/// `false` if some rotation `k > 0` makes the rotation lex-smaller
/// than the identity within the fully-known comparable region. The
/// caller pairs this with a complementary loop for the dihedral case
/// (see [`is_dihedral_canonical_extended`]).
fn rotation_lex_min_violated(prefix: &[i8], new: i8) -> bool {
    let d = prefix.len() + 1;
    for k in 1..d {
        for i in 0..(d - k) {
            match walk_get(prefix, new, k + i).cmp(&walk_get(prefix, new, i)) {
                std::cmp::Ordering::Less => return true,
                std::cmp::Ordering::Greater => break,
                std::cmp::Ordering::Equal => continue,
            }
        }
    }
    false
}

/// Canonical-rotation prune for the open walk `prefix` virtually
/// extended by one more angle `new`. Lets the caller skip walks
/// that cannot be the lex-min cyclic rotation of any closure they
/// could grow into, *before* paying the cyclotomic intersect cost
/// of `Snake::add`.
///
/// Returns `false` when there's a rotation index `k > 0` such that
/// the rotation `[prefix++new][k..]` is strictly lex-less than
/// `[prefix++new][0..]` within the wrap-free comparable region; that
/// decision is permanent (no future angles can flip it) so the
/// eventual closed polygon's canonical rotation cannot start at
/// position 0 -- pruning is safe.
pub fn is_canonical_extended(prefix: &[i8], new: i8) -> bool {
    !rotation_lex_min_violated(prefix, new)
}

pub fn canonical_identity(seq: &[i8]) -> Vec<i8> {
    seq.to_vec()
}

/// Heuristic dihedral prefix prune. Like [`is_canonical_extended`]
/// but also compares the walk prefix against complement rotations
/// (negated angle sequences). Returns `false` when any complement
/// rotation of the eventual closure has a lex-smaller prefix than
/// the identity walk, indicating the walk is suboptimal under the
/// dihedral group.
///
/// This is sound (never prunes a walk that could produce the
/// dihedral-min output) but incomplete (does not eliminate all
/// chiral duplicates -- see module-level comment above). At depth 0
/// it prunes all positive first angles (`-a_0 < a_0`), roughly
/// halving the root branching factor.
///
/// # Index bounds
///
/// Complement rotation `k` compares `-a_{k - i}` with `a_i` at
/// positions `i = 0..=k`. Both sides are decided by the known walk
/// only when `0 <= k - i < d` and `0 <= i < d`; for the inner-loop
/// range `i = 0..=k` to stay within the prefix we need `k < d`.
/// Hence the outer loop is `for k in 0..d` (NOT `0..=d`) -- accessing
/// `walk_get(prefix, new, d)` would silently return `new` again for
/// a position whose value is genuinely undetermined, which would
/// over-prune.
pub fn is_dihedral_canonical_extended(prefix: &[i8], new: i8) -> bool {
    if rotation_lex_min_violated(prefix, new) {
        return false;
    }

    let d = prefix.len() + 1;
    for k in 0..d {
        for i in 0..=k {
            match (-walk_get(prefix, new, k - i)).cmp(&walk_get(prefix, new, i)) {
                std::cmp::Ordering::Less => return false,
                std::cmp::Ordering::Greater => break,
                std::cmp::Ordering::Equal => continue,
            }
        }
    }

    true
}

/// Correctness layer for the dihedral DFS. Computes the lex-minimum
/// over all cyclic rotations and reversed rotations of `seq`.
///
/// The inputs are already chirality-normalized to CW by the DFS
/// before this function runs. Under that normalization, two mirror
/// images of the same polygon shape (an enantiomer pair) end up as
/// CW walks of the original and the mirror. Algebraically these CW
/// walks differ by SEQUENCE REVERSAL (one polygon walked CW from
/// vertex V0, vs its mirror walked CW from V0' = mirror of V0,
/// gives sequences related by reverse). Plain reversal -- not
/// revcomp -- is therefore the right relation here.
///
/// Both members of an enantiomer pair map to the same value, so the
/// HashSet deduplicates them. Without this step, the output would
/// contain both members of every non-achiral chiral pair.
pub fn dihedral_canonical(seq: &[i8]) -> Vec<i8> {
    let n = seq.len();
    let mut best: Vec<i8> = seq.to_vec();
    let mut rot: Vec<i8> = seq.to_vec();
    for _ in 0..n {
        if rot < best {
            best = rot.clone();
        }
        let mut rev = rot.clone();
        rev.reverse();
        if rev < best {
            best = rev.clone();
        }
        rot.rotate_left(1);
    }
    best
}

/// Build the [`CanonicalOps`] pair from the `dihedral` CLI flag.
pub fn make_ops(dihedral: bool) -> CanonicalOps {
    if dihedral {
        CanonicalOps {
            is_canonical: is_dihedral_canonical_extended,
            canonicalize: dihedral_canonical,
        }
    } else {
        CanonicalOps {
            is_canonical: is_canonical_extended,
            canonicalize: canonical_identity,
        }
    }
}
