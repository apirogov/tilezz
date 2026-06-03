//! Multi-tile patch enumeration by incremental boundary glue.
//!
//! [`enum_patches`] is the public entry point: given an arbitrary
//! [`TileSet`] and a maximum size, it enumerates every distinct
//! patch of up to `max_size` tile copies (drawn from any combination
//! of tiles in the set) and returns them grouped by size.
//!
//! # Algorithm
//!
//! Plain layer-BFS over the boundary `Rat`. At size 1 every tile is
//! its own patch. At each subsequent layer the surviving patches
//! become the A-side of a [`MatchFinder`] crossing against the
//! original tileset as B-side; the union of valid glue products is
//! the next layer.
//!
//! Snake-intersection validation and canonical `Rat` dedup are
//! handled inside `MatchFinder`. No bespoke geometry, no per-tile
//! special-casing — works the same for a single seed (e.g. hex,
//! square) or a mixed tileset (e.g. square+hex, tetrominoes).
//!
//! # Memory
//!
//! Layer-BFS holds every patch of the previous layer in memory
//! simultaneously to form the A-side `TileSet`. For tilesets that
//! explode (spectre at size ≥ 5, large mixed sets) this can be
//! significant — but the structure stays simple and the per-patch
//! cost (one `Rat` clone) is minimal.

use std::collections::BTreeMap;
use std::sync::Arc;

use rustc_hash::FxHashSet;

use crate::analysis::matchfinder::{BpSeed, MatchFinder};
use crate::cyclotomic::IsRing;
use crate::geom::rat::Rat;
use crate::geom::tileset::TileSet;

/// Enumerate every distinct **onesided** patch of up to `max_size`
/// tile copies drawn from `tileset`. Patches are keyed by their
/// canonical boundary `Rat`; rotations of the same shape are
/// counted once, but chiral reflections count as two.
///
/// Returns a map `size -> {patches of that size}`. Empty `tileset`
/// or `max_size == 0` yield an empty result.
pub fn enum_patches<T>(
    tileset: Arc<TileSet<T>>,
    max_size: usize,
) -> BTreeMap<usize, FxHashSet<Rat<T>>>
where
    T: IsRing,
{
    let mut results: BTreeMap<usize, FxHashSet<Rat<T>>> = BTreeMap::new();
    if max_size == 0 || tileset.num_tiles() == 0 {
        return results;
    }

    // Size 1: every tile is a patch of size 1.
    results.insert(1, tileset.rats().iter().cloned().collect());

    // B-side bit-parallel state, built once and shared across layers.
    let bp_seed = BpSeed::new(Arc::clone(&tileset));

    for k in 2..=max_size {
        let prev = match results.get(&(k - 1)) {
            Some(s) if !s.is_empty() => s,
            _ => break,
        };
        let patches_vec: Vec<Rat<T>> = prev.iter().cloned().collect();
        let patch_ts = Arc::new(TileSet::new(patches_vec));
        let mf = MatchFinder::crossing_with_seed(Arc::clone(&patch_ts), bp_seed.clone());

        let pairs: Vec<(usize, usize)> = (0..mf.num_tiles_a())
            .flat_map(|i| (0..mf.num_tiles_b()).map(move |j| (i, j)))
            .collect();

        let layer: FxHashSet<Rat<T>> = mf.valid_results_for_pairs(&pairs).into_iter().collect();
        if layer.is_empty() {
            break;
        }
        results.insert(k, layer);
    }

    results
}

/// Like [`enum_patches`] but identifies chiral pairs: every patch
/// is stored as `min(rat, rat.reflected())`, so a shape and its
/// mirror image collapse to a single entry. The result is the set
/// of **free** patches.
pub fn enum_patches_free<T>(
    tileset: Arc<TileSet<T>>,
    max_size: usize,
) -> BTreeMap<usize, FxHashSet<Rat<T>>>
where
    T: IsRing,
{
    enum_patches(tileset, max_size)
        .into_iter()
        .map(|(k, set)| {
            let free: FxHashSet<Rat<T>> = set
                .into_iter()
                .map(|r| std::cmp::min(r.clone(), r.reflected()))
                .collect();
            (k, free)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic::{ZZ4, ZZ12};
    use crate::geom::tileset;

    /// Independent brute-force reference: for each size, glue every
    /// `(patch, tile, ia, ib)` combination via `try_glue` and collect
    /// the resulting Rats. No optimization, no canonicalization
    /// shortcuts — just the most direct enumeration possible.
    fn brute_force_grow<T>(
        tileset: &Arc<TileSet<T>>,
        max_size: usize,
    ) -> BTreeMap<usize, FxHashSet<Rat<T>>>
    where
        T: IsRing,
    {
        let mut results: BTreeMap<usize, FxHashSet<Rat<T>>> = BTreeMap::new();
        if max_size == 0 || tileset.num_tiles() == 0 {
            return results;
        }
        results.insert(1, tileset.rats().iter().cloned().collect());
        for k in 2..=max_size {
            let prev: Vec<Rat<T>> = results[&(k - 1)].iter().cloned().collect();
            let mut next: FxHashSet<Rat<T>> = FxHashSet::default();
            for patch in &prev {
                for tile_idx in 0..tileset.num_tiles() {
                    let tile = tileset.rat(tile_idx);
                    for ia in 0..patch.len() {
                        for ib in 0..tile.len() {
                            if let Ok(glued) = patch.try_glue((ia as i64, ib as i64), tile) {
                                next.insert(glued);
                            }
                        }
                    }
                }
            }
            results.insert(k, next);
        }
        results
    }

    fn assert_grow_matches_brute_force<T>(label: &str, tileset: Arc<TileSet<T>>, max_size: usize)
    where
        T: IsRing,
    {
        let fast = enum_patches(Arc::clone(&tileset), max_size);
        let brute = brute_force_grow(&tileset, max_size);
        for k in 1..=max_size {
            let fast_set = fast.get(&k).cloned().unwrap_or_default();
            let brute_set = brute.get(&k).cloned().unwrap_or_default();
            assert_eq!(
                fast_set.len(),
                brute_set.len(),
                "{label} size {k}: enum_patches={} brute_force={}",
                fast_set.len(),
                brute_set.len()
            );
            for r in &brute_set {
                assert!(
                    fast_set.contains(r),
                    "{label} size {k}: brute force has rat {r:?} that enum_patches missed"
                );
            }
            for r in &fast_set {
                assert!(
                    brute_set.contains(r),
                    "{label} size {k}: enum_patches has rat {r:?} that brute force missed"
                );
            }
        }
    }

    /// OEIS A000105: number of free polyominoes (no holes) with n cells.
    /// Indices 0..13 (we use 1..max_size).
    const FREE_POLYOMINOES_NO_HOLES: &[usize] = &[
        1, 1, 1, 2, 5, 12, 35, 107, 363, 1248, 4460, 16094, 58937, 217117,
    ];

    /// Free polyomino counts via `enum_patches_free` must match OEIS
    /// A000105 — this is the strongest external sanity check we have.
    #[test]
    fn square_zz4_free_polyominoes_match_oeis() {
        let max_size = 9;
        let ts = tileset::square::<ZZ4>();
        let free = enum_patches_free(ts, max_size);
        for (k, &expected) in FREE_POLYOMINOES_NO_HOLES
            .iter()
            .enumerate()
            .take(max_size + 1)
            .skip(1)
        {
            let got = free.get(&k).map(|s| s.len()).unwrap_or(0);
            assert_eq!(got, expected, "size {k}: got {got}, expected {expected}");
        }
    }

    /// Hexagon-only (ZZ12): cross-check enum_patches vs brute-force
    /// `try_glue` up to size 4. Verifies the BP-accelerated path
    /// agrees with the most direct enumeration on a non-square ring.
    #[test]
    fn hex_zz12_matches_brute_force_size4() {
        assert_grow_matches_brute_force("hex", tileset::hex::<ZZ12>(), 4);
    }

    /// Spectre (ZZ12, k=14 boundary): cross-check at size 3.
    /// Larger boundary, sparser matches; exercises the BP single-edge
    /// path's interaction with non-trivial tile shapes.
    #[test]
    fn spectre_zz12_matches_brute_force_size3() {
        assert_grow_matches_brute_force("spectre", tileset::spectre::<ZZ12>(), 3);
    }

    /// **Multi-tile** cross-check: mixed square + hexagon tileset
    /// over ZZ12. This is the case the deleted Redelmeier code never
    /// handled — verifies enum_patches generalizes correctly.
    #[test]
    fn mixed_square_hex_zz12_matches_brute_force_size3() {
        assert_grow_matches_brute_force("mixed", tileset::mixed::<ZZ12>(), 3);
    }

    /// **Multi-tile** stress: all seven tetrominoes (ZZ4) up to
    /// size 3. Confirms the multi-tile path handles 7 distinct tile
    /// shapes correctly.
    #[test]
    fn tetrominoes_zz4_matches_brute_force_size3() {
        assert_grow_matches_brute_force("tetrominoes", tileset::tetrominoes::<ZZ4>(), 3);
    }

    /// Empty / degenerate inputs.
    #[test]
    fn max_size_zero_returns_empty() {
        let result = enum_patches(tileset::square::<ZZ4>(), 0);
        assert!(result.is_empty());
    }

    /// At size 1 the result is exactly the input tileset.
    #[test]
    fn size_one_equals_input_tileset() {
        let ts = tileset::mixed::<ZZ12>();
        let result = enum_patches(Arc::clone(&ts), 1);
        let size_1 = result.get(&1).expect("size 1 must be present");
        assert_eq!(size_1.len(), ts.num_tiles());
        for r in ts.rats() {
            assert!(size_1.contains(r));
        }
    }
}
