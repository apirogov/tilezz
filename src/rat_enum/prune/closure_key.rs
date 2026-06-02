//! Closure-key prune.
//!
//! Pre-pass: DFS-enumerate every simple open snake up to length `L`
//! and store the set of `(endpoint, facing)` pairs they reach -- the
//! "closure keys" `K_{<=L}`. During the main DFS, when
//! `remaining_after <= L`, the candidate prefix can close iff its
//! required suffix's `(endpoint, facing)` is in `K_{<=L}`. Otherwise
//! no length-`<=remaining_after` continuation can complete it.
//!
//! Strictly stronger than any modular projection (it sees exact
//! lattice + facing info, not just a quotient). The cost is memory +
//! pre-pass time + a ring multiplication per hot-path check (to
//! compute the target suffix endpoint via
//! `target = -unit(-facing) * disp`).

use crate::cyclotomic::{IsRing, Units, ZZ10, ZZ12, ZZ16, ZZ20, ZZ24, ZZ32, ZZ4, ZZ60, ZZ8};
use crate::geom::snake::Snake;

/// Set of closure keys indexed by ring-coefficient vector + facing.
pub struct ClosureKeyPrune {
    /// Maximum tabulated suffix length. The prune fires only when
    /// `remaining_after <= max_l`.
    pub max_l: usize,
    /// `K_{<=L} = {(endpoint coords, facing) : reachable by some
    /// simple open snake of length 0..=L}`. Length 0 (the empty
    /// snake, key `(0, 0)`) is included so already-closed candidates
    /// aren't incorrectly rejected.
    pub keys: rustc_hash::FxHashSet<(Vec<i64>, i8)>,
}

fn collect_closure_keys_dfs<ZZ: IsRing>(
    snake: &mut Snake<ZZ>,
    max_l: usize,
    keys: &mut rustc_hash::FxHashSet<(Vec<i64>, i8)>,
) {
    if snake.angles().len() >= max_l {
        return;
    }
    let turn = ZZ::turn();
    for direction in ((-ZZ::hturn() + 1)..ZZ::hturn()).rev() {
        if !snake.add(direction) {
            continue;
        }
        // Normalize facing to 0..n-1: `Snake::direction` returns
        // truncating-mod (can be negative for net-CW walks); the
        // hot-path lookup uses `rem_euclid` (always 0..n-1), so we
        // canonicalize here to make the hash keys compatible.
        let facing = snake.direction().rem_euclid(turn);
        keys.insert((snake.offset().int_coeffs_slice().to_vec(), facing));
        collect_closure_keys_dfs::<ZZ>(snake, max_l, keys);
        snake.pop();
    }
}

pub fn collect_closure_keys<ZZ: IsRing>(
    max_l: usize,
) -> rustc_hash::FxHashSet<(Vec<i64>, i8)> {
    let mut keys = rustc_hash::FxHashSet::default();
    // Empty snake's key: closure is already achieved.
    let phi = <ZZ as Units>::unit(0).int_coeffs_slice().len();
    keys.insert((vec![0i64; phi], 0));
    let mut snake: Snake<ZZ> = Snake::new();
    collect_closure_keys_dfs::<ZZ>(&mut snake, max_l, &mut keys);
    keys
}

/// Ring-dispatch wrapper for [`collect_closure_keys`].
pub fn collect_closure_keys_for_ring(
    ring: u8,
    max_l: usize,
) -> rustc_hash::FxHashSet<(Vec<i64>, i8)> {
    match ring {
        4 => collect_closure_keys::<ZZ4>(max_l),
        8 => collect_closure_keys::<ZZ8>(max_l),
        10 => collect_closure_keys::<ZZ10>(max_l),
        12 => collect_closure_keys::<ZZ12>(max_l),
        16 => collect_closure_keys::<ZZ16>(max_l),
        20 => collect_closure_keys::<ZZ20>(max_l),
        24 => collect_closure_keys::<ZZ24>(max_l),
        32 => collect_closure_keys::<ZZ32>(max_l),
        60 => collect_closure_keys::<ZZ60>(max_l),
        _ => panic!("invalid ring selected"),
    }
}
