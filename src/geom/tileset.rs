//! A `TileSet` is a small immutable collection of tile shapes (`Rat`s)
//! used as the alphabet for patch construction.
//!
//! - All tiles must have the same chirality (asserted at construction).
//! - Tiles are sorted into a deterministic canonical order and deduped,
//!   so logically-equal tilesets compare equal regardless of input
//!   order.
//! - Tile indices `0..num_tiles()` are the stable identifiers used by
//!   the rest of the patch/match machinery (`tile_id` fields on
//!   `EdgeInfo`, `PatchMatch`, etc.).

use std::sync::Arc;

use crate::cyclotomic::{HasZZ4, HasZZ6, HasZZ10, HasZZ12, IsRing};
use crate::geom::rat::Rat;
use crate::geom::tiles;

/// An immutable, ordered, deduplicated collection of tile shapes.
///
/// All tiles in a `TileSet` share the same chirality. Tiles are stored
/// in their `Rat`-canonical sort order; `TileSet::rat(i)` returns the
/// tile at index `i` in that order. The index is stable for the
/// lifetime of the `TileSet` and is what `tile_id` fields in
/// `PatchMatch` / `EdgeInfo` / etc. refer to.
pub struct TileSet<T: IsRing> {
    rats: Vec<Rat<T>>,
}

impl<T: IsRing> TileSet<T> {
    /// Build a `TileSet` from the given tile shapes.
    ///
    /// Sorts and deduplicates the input. Panics if `rats` is empty or
    /// if the tiles have mixed chirality (some CW and some CCW).
    pub fn new(mut rats: Vec<Rat<T>>) -> Self {
        assert!(!rats.is_empty(), "need at least one tile");

        let chirality = rats[0].chirality();
        for (i, r) in rats.iter().enumerate() {
            assert_eq!(
                r.chirality(),
                chirality,
                "mixed chirality: tile {i} has chirality {} but expected {}",
                r.chirality(),
                chirality,
            );
        }

        rats.sort();
        rats.dedup();

        TileSet { rats }
    }

    /// Number of distinct tile shapes in the set.
    pub fn num_tiles(&self) -> usize {
        self.rats.len()
    }

    /// Tile at index `i` (`0..num_tiles()`).
    pub fn rat(&self, i: usize) -> &Rat<T> {
        &self.rats[i]
    }

    /// All tiles in their sorted/deduped order.
    pub fn rats(&self) -> &[Rat<T>] {
        &self.rats
    }

    /// Look up the index of a specific tile, if present. Uses a binary
    /// search over the canonical-sort order.
    #[allow(dead_code)]
    pub(crate) fn index_of(&self, rat: &Rat<T>) -> Option<usize> {
        self.rats.binary_search(rat).ok()
    }
}

// ---- Stock tileset constructors ----
//
// Common tile families packaged as ready-to-use `Arc<TileSet<T>>`,
// generic over a compatible cyclotomic ring. Each helper picks tile
// shapes from [`crate::geom::tiles`] and wraps them in a `TileSet`.
// Use these to avoid the same `Rat::from_unchecked(&tiles::foo()) →
// TileSet::new(vec![..]) → Arc::new(..)` boilerplate at every call site.

/// Singleton hexagon tileset over a ZZ6-compatible ring.
pub fn hex<T: HasZZ6>() -> Arc<TileSet<T>> {
    Arc::new(TileSet::new(vec![Rat::from_unchecked(
        &tiles::hexagon::<T>(),
    )]))
}

/// Singleton square tileset over a ZZ4-compatible ring.
pub fn square<T: HasZZ4>() -> Arc<TileSet<T>> {
    Arc::new(TileSet::new(vec![Rat::from_unchecked(
        &tiles::square::<T>(),
    )]))
}

/// Square + hexagon tileset over a ring supporting both ZZ4 and ZZ6
/// (the smallest concrete choice is `ZZ12`).
pub fn mixed<T: HasZZ4 + HasZZ6>() -> Arc<TileSet<T>> {
    Arc::new(TileSet::new(vec![
        Rat::from_unchecked(&tiles::square::<T>()),
        Rat::from_unchecked(&tiles::hexagon::<T>()),
    ]))
}

/// All seven standard tetrominoes (O, I, T, S, Z, J, L) as a single
/// tileset over a ZZ4-compatible ring.
pub fn tetrominoes<T: HasZZ4>() -> Arc<TileSet<T>> {
    Arc::new(TileSet::new(vec![
        Rat::from_unchecked(&tiles::tetromino_O::<T>()),
        Rat::from_unchecked(&tiles::tetromino_I::<T>()),
        Rat::from_unchecked(&tiles::tetromino_T::<T>()),
        Rat::from_unchecked(&tiles::tetromino_S::<T>()),
        Rat::from_unchecked(&tiles::tetromino_Z::<T>()),
        Rat::from_unchecked(&tiles::tetromino_J::<T>()),
        Rat::from_unchecked(&tiles::tetromino_L::<T>()),
    ]))
}

/// Singleton spectre tileset over a ZZ12-compatible ring.
pub fn spectre<T: HasZZ12>() -> Arc<TileSet<T>> {
    Arc::new(TileSet::new(vec![Rat::from_unchecked(
        &tiles::spectre::<T>(),
    )]))
}

/// Penrose P3 (narrow + wide rhomb) tileset over a ZZ10-compatible
/// ring.
pub fn penrose<T: HasZZ10>() -> Arc<TileSet<T>> {
    Arc::new(TileSet::new(vec![
        Rat::from_unchecked(&tiles::penrose_p3_narrow::<T>()),
        Rat::from_unchecked(&tiles::penrose_p3_wide::<T>()),
    ]))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic::ZZ12;
    use crate::geom::tiles::{hexagon, spectre};

    #[test]
    #[should_panic(expected = "mixed chirality")]
    fn mixed_chirality_panics() {
        let ccw = Rat::<ZZ12>::from_slice_unchecked(&[1, 1, 1, 1]);
        let cw = ccw.reversed();
        let _ts = TileSet::new(vec![ccw, cw]);
    }

    #[test]
    fn dedup_removes_duplicates() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&hexagon());
        let ts = TileSet::new(vec![r.clone(), r.clone(), r.clone()]);
        assert_eq!(ts.num_tiles(), 1);
    }

    #[test]
    fn sorts_canonically() {
        let r1 = Rat::<ZZ12>::from_slice_unchecked(&[4, 4, 4]);
        let r2 = Rat::<ZZ12>::from_slice_unchecked(&[1, 1, 1, 1]);
        let ts = TileSet::new(vec![r1, r2]);
        assert_eq!(ts.num_tiles(), 2);
        assert_eq!(ts.rat(0).len(), 4);
        assert_eq!(ts.rat(1).len(), 3);
    }

    #[test]
    fn index_of_finds_tile() {
        let r1: Rat<ZZ12> = Rat::from_unchecked(&hexagon());
        let r2: Rat<ZZ12> = Rat::from_unchecked(&spectre());
        let ts = TileSet::new(vec![r1.clone(), r2.clone()]);
        assert!(ts.index_of(&r1).is_some());
        assert!(ts.index_of(&r2).is_some());
    }

    #[test]
    fn index_of_missing_returns_none() {
        let r1: Rat<ZZ12> = Rat::from_unchecked(&hexagon());
        let r2: Rat<ZZ12> = Rat::from_unchecked(&spectre());
        let ts = TileSet::new(vec![r1]);
        assert!(ts.index_of(&r2).is_none());
    }

    #[test]
    #[should_panic(expected = "need at least one tile")]
    fn empty_panics() {
        let _: TileSet<ZZ12> = TileSet::new(vec![]);
    }

    #[test]
    fn preserves_unique_tiles() {
        let r1: Rat<ZZ12> = Rat::from_unchecked(&hexagon());
        let r2: Rat<ZZ12> = Rat::from_unchecked(&spectre());
        let ts = TileSet::new(vec![r1, r2]);
        assert_eq!(ts.num_tiles(), 2);
    }

    #[test]
    fn chirality_check_passes_for_consistent_tiles() {
        let r1 = Rat::<ZZ12>::from_slice_unchecked(&[0, 1, 0, 1, 0, 1, 0, 1]);
        let r2 = Rat::<ZZ12>::from_slice_unchecked(&[4, 4, 4]);
        let ts = TileSet::new(vec![r1, r2]);
        assert_eq!(ts.num_tiles(), 2);
    }
}
