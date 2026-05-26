use rustc_hash::FxHashMap;

use crate::cyclotomic::IsRing;

/// Precomputed unions of the two 5-cell `+` crosses centered on cells
/// that are `(dx, dy)` apart, indexed by `(dx + 1) * 3 + (dy + 1)` for
/// `(dx, dy) in {-1, 0, 1}^2`. All 9 offsets relative to `cell(p1)`,
/// sorted and deduplicated.
///
/// `cell_of` uses `floor`, so a unit-length segment is guaranteed to
/// have `(dx, dy)` in this range -- no boundary cases.
static NEIGHBOR_OFFSETS: [&[(i64, i64)]; 9] = [
    // (dx=-1, dy=-1): p2 down-left of p1
    &[(-2, -1), (-1, -2), (-1, -1), (-1, 0), (0, -1), (0, 0), (0, 1), (1, 0)],
    // (dx=-1, dy=0): p2 left of p1
    &[(-2, 0), (-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 0), (0, 1), (1, 0)],
    // (dx=-1, dy=1): p2 up-left of p1
    &[(-2, 1), (-1, 0), (-1, 1), (-1, 2), (0, -1), (0, 0), (0, 1), (1, 0)],
    // (dx=0, dy=-1): p2 below p1
    &[(-1, -1), (-1, 0), (0, -2), (0, -1), (0, 0), (0, 1), (1, -1), (1, 0)],
    // (dx=0, dy=0): same cell, single 5-cell cross
    &[(-1, 0), (0, -1), (0, 0), (0, 1), (1, 0)],
    // (dx=0, dy=1): p2 above p1
    &[(-1, 0), (-1, 1), (0, -1), (0, 0), (0, 1), (0, 2), (1, 0), (1, 1)],
    // (dx=1, dy=-1): p2 down-right of p1
    &[(-1, 0), (0, -1), (0, 0), (0, 1), (1, -2), (1, -1), (1, 0), (2, -1)],
    // (dx=1, dy=0): p2 right of p1
    &[(-1, 0), (0, -1), (0, 0), (0, 1), (1, -1), (1, 0), (1, 1), (2, 0)],
    // (dx=1, dy=1): p2 up-right of p1
    &[(-1, 0), (0, -1), (0, 0), (0, 1), (1, 0), (1, 1), (1, 2), (2, 1)],
];

/// Spatial hash from integer unit-square cells to a list of `usize` IDs
/// stored at each cell.
///
/// Cells are the half-open unit squares `[cx, cx+1) x [cy, cy+1)` of
/// the Cartesian plane, indexed by their lower-left integer corner.
/// What an "ID" means is up to the caller -- it's usually an index
/// into the caller's own data:
///
/// - [`crate::geom::snake::Snake`] stores **point indices**: each
///   boundary vertex is filed at one cell (its `cell_of(point)`).
///   Lookups find which other boundary vertices share that cell or a
///   nearby one.
/// - [`crate::geom::patch::GrowingPatch`] stores **segment IDs**:
///   each boundary edge is filed at *every* cell its
///   [`Self::seg_neighborhood_of`] reaches (up to 8 cells per
///   segment). Lookups during incremental glue checks ask "is this
///   candidate segment within striking distance of an existing
///   segment?" by reading the neighborhood cells.
///
/// The grid is a multimap: multiple IDs can share a cell, and the same
/// ID can appear in many cells (the segment use case above). Inserts
/// (`add`) and removals (`remove`) are O(1) per cell; cell lookups
/// (`get`, `get_cells`) are O(1) per cell plus the size of the cell's
/// payload.
#[derive(Debug, Clone)]
pub struct UnitSquareGrid {
    cells: FxHashMap<(i64, i64), Vec<usize>>,
}

impl Default for UnitSquareGrid {
    fn default() -> Self {
        Self::new()
    }
}

impl UnitSquareGrid {
    /// Compute the half-open cell `[cx, cx+1) x [cy, cy+1)` containing
    /// `zz`, i.e. `(floor(Re(zz)), floor(Im(zz)))`. Delegates to the
    /// ring's `CellFloor::cell_floor` impl.
    ///
    /// # Speed and exactness
    ///
    /// The default `cell_floor` is the fast f64 floor of
    /// `complex64()`. In **debug builds** it `debug_assert!`s the
    /// result against the per-ring `cell_floor_exact` (an
    /// integer-arithmetic sign-verification path that's exact at all
    /// points, including the measure-zero half-integer grid lines).
    /// In **release**, no verification -- the f64 path stands alone.
    ///
    /// This is safe in practice because canonical ring elements with
    /// small integer coefficients essentially never land exactly on a
    /// half-integer boundary; the f64 floor matches `cell_floor_exact`
    /// on every snake / patch state the codebase actually constructs.
    /// Debug builds catch any boundary-case divergence in tests.
    ///
    /// Every blessed ring has an exact path available (`HasZZ4` rings
    /// via cell-corner construction + the `rect_signs` sign primitive;
    /// `ZZ10` via per-axis pentagonal sign helpers, since `i` is not
    /// in the ring and corner construction isn't available). Callers
    /// that need bit-exact behaviour in release can invoke
    /// `zz.cell_floor_exact()` directly.
    ///
    /// # Why floor (not round)
    ///
    /// A unit-length segment never spans more than 2 cells along an
    /// axis. With round-half-away-from-zero, points at exactly
    /// +/-0.5 round in opposite directions, so a unit segment
    /// crossing both sides of zero would give a cell-offset of 2 --
    /// forcing the neighborhood lookup table in
    /// [`Self::seg_neighborhood_of`] to enumerate more cases than the
    /// conceptual 9. Floor partitions space into `[a, a+1)` cells, so
    /// the unit shift always changes the floor by 0 or +/-1.
    pub fn cell_of<T: IsRing>(zz: T) -> (i64, i64) {
        zz.cell_floor()
    }

    /// Cells covering the neighborhood of a **unit-length** segment from
    /// `p1` to `p2`, with no intermediate allocations and no runtime sort
    /// or dedup.
    ///
    /// Both endpoints must be at most one unit apart, so the cell offset
    /// `(dx, dy) = cell(p2) - cell(p1)` always lies in `{-1, 0, 1}^2` --
    /// 9 cases. For each, the union of the two 5-cell `+` crosses (one
    /// centered on each endpoint's cell) is a fixed list of at most 8
    /// distinct cells, precomputed and sorted at compile time. The
    /// runtime lookup is a single array index on the packed offset
    /// `(dx + 1) * 3 + (dy + 1)`.
    ///
    /// Yields absolute cell coordinates `(cell(p1).0 + ox, cell(p1).1 + oy)`.
    /// Panics in debug if `(dx, dy)` escapes `{-1, 0, 1}^2` (which would
    /// indicate non-unit input or a `cell_of` rounding bug).
    #[inline]
    pub fn seg_neighborhood_of<T: IsRing>(
        p1: T,
        p2: T,
    ) -> impl Iterator<Item = (i64, i64)> {
        let (cx, cy) = Self::cell_of(p1);
        let (cx2, cy2) = Self::cell_of(p2);
        let dx = cx2 - cx;
        let dy = cy2 - cy;
        debug_assert!(
            dx.abs() <= 1 && dy.abs() <= 1,
            "seg_neighborhood_of: endpoints not unit-distance ({dx}, {dy})"
        );
        let idx = ((dx + 1) * 3 + (dy + 1)) as usize;
        NEIGHBOR_OFFSETS[idx]
            .iter()
            .map(move |&(ox, oy)| (cx + ox, cy + oy))
    }

    /// Empty grid; no cells, no entries.
    pub fn new() -> Self {
        Self {
            cells: FxHashMap::default(),
        }
    }

    /// Total number of stored IDs (summed across all cells). An ID
    /// filed in `k` cells contributes `k` to the count -- this is
    /// `cell-entries` count, not "distinct IDs".
    pub fn len(&self) -> usize {
        self.cells.values().map(|v| v.len()).sum()
    }

    pub fn is_empty(&self) -> bool {
        self.cells.values().all(|v| v.is_empty())
    }

    /// File `value` at `cell`. Multiple values at the same cell are
    /// kept in insertion order (a `Vec`), and the same value can be
    /// filed at the same cell multiple times -- the grid is a
    /// multiset per cell, not a set.
    pub fn add(&mut self, cell: (i64, i64), value: usize) {
        self.cells.entry(cell).or_default().push(value);
    }

    /// Remove one occurrence of `value` from `cell`. No-op if `cell`
    /// has no entries or doesn't contain `value`. When the last
    /// entry of a cell is removed, the cell itself is dropped from
    /// the underlying map (so subsequent `get(cell)` returns the
    /// empty slice via the `unwrap_or(&[])` fallback rather than
    /// via an empty `Vec`).
    pub fn remove(&mut self, cell: (i64, i64), value: usize) {
        if let Some(vals) = self.cells.get_mut(&cell) {
            if let Some(pos) = vals.iter().position(|&v| v == value) {
                vals.swap_remove(pos);
            }
            if vals.is_empty() {
                self.cells.remove(&cell);
            }
        }
    }

    /// Return the IDs stored at a single cell, in insertion order.
    /// Empty slice if the cell has no entries (never been `add`ed to,
    /// or all its entries have been `remove`d).
    ///
    /// Borrow-only; no allocation. For querying a fixed list of
    /// cells in one call (the typical "what's in this neighborhood?"
    /// pattern), prefer [`Self::get_cells`] which fuses the lookups
    /// and returns an owned `Vec`.
    pub fn get(&self, cell: (i64, i64)) -> &[usize] {
        self.cells.get(&cell).map(Vec::as_slice).unwrap_or(&[])
    }

    /// Return the concatenation of all IDs stored across the given
    /// cells, in the order: outer = `cells` iteration order, inner =
    /// insertion order within each cell.
    ///
    /// This is the natural pair with [`Self::seg_neighborhood_of`]:
    /// given a candidate unit-length segment, ask the grid for every
    /// ID stored in any of the (up to) 8 cells the segment could
    /// possibly interact with -- one allocation, one pass.
    ///
    /// IDs that live in multiple of the queried cells appear once per
    /// containing cell. Callers that need uniqueness (e.g. unique
    /// segment IDs from a segment-spreading grid like
    /// `GrowingPatch`'s) must dedup the returned `Vec` themselves.
    pub fn get_cells(&self, cells: &[(i64, i64)]) -> Vec<usize> {
        let total_entries: usize = cells
            .iter()
            .map(|&cell| self.cells.get(&cell).map(Vec::len).unwrap_or(0))
            .sum();

        let mut result = Vec::with_capacity(total_entries);
        for &cell in cells {
            if let Some(vals) = self.cells.get(&cell) {
                result.extend(vals);
            }
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic::{SymNum, Units, ZZ12};
    use num_traits::{One, Zero};

    #[test]
    fn test_cell_of() {
        // `cell_of` uses floor: cells are `[a, a+1) x [b, b+1)`.
        // Unit cartesian coordinates for zeta^k:
        //   k=0:  (1.000,  0.000)
        //   k=1:  (0.866,  0.500)
        //   k=2:  (0.500,  0.866)
        //   k=3:  (0.000,  1.000)
        //   k=4: (-0.500,  0.866)
        //   k=5: (-0.866,  0.500)
        //   k=6: (-1.000,  0.000)
        //   k=7: (-0.866, -0.500)
        //   k=8: (-0.500, -0.866)
        //   k=9:  (0.000, -1.000)
        //   k=10: (0.500, -0.866)
        //   k=11: (0.866, -0.500)
        assert_eq!(UnitSquareGrid::cell_of(ZZ12::zero()), (0, 0));
        assert_eq!(UnitSquareGrid::cell_of(<ZZ12 as Units>::unit(0)), (1, 0));
        assert_eq!(UnitSquareGrid::cell_of(<ZZ12 as Units>::unit(1)), (0, 0));
        assert_eq!(UnitSquareGrid::cell_of(<ZZ12 as Units>::unit(2)), (0, 0));
        assert_eq!(UnitSquareGrid::cell_of(<ZZ12 as Units>::unit(3)), (0, 1));
        assert_eq!(UnitSquareGrid::cell_of(<ZZ12 as Units>::unit(4)), (-1, 0));
        assert_eq!(UnitSquareGrid::cell_of(<ZZ12 as Units>::unit(5)), (-1, 0));
        assert_eq!(UnitSquareGrid::cell_of(<ZZ12 as Units>::unit(6)), (-1, 0));
        assert_eq!(UnitSquareGrid::cell_of(<ZZ12 as Units>::unit(7)), (-1, -1));
        assert_eq!(UnitSquareGrid::cell_of(<ZZ12 as Units>::unit(8)), (-1, -1));
        assert_eq!(UnitSquareGrid::cell_of(<ZZ12 as Units>::unit(9)), (0, -1));
        assert_eq!(UnitSquareGrid::cell_of(<ZZ12 as Units>::unit(10)), (0, -1));
        assert_eq!(UnitSquareGrid::cell_of(<ZZ12 as Units>::unit(11)), (0, -1));
    }

    #[test]
    fn test_seg_neighborhood_of() {
        // `seg_neighborhood_of` returns the sorted-deduplicated union of
        // the two 5-cell `+` crosses around the endpoint cells, as an
        // iterator. We collect to compare against expected lists.
        let cells_of = |p1: ZZ12, p2: ZZ12| -> Vec<(i64, i64)> {
            UnitSquareGrid::seg_neighborhood_of(p1, p2).collect()
        };

        let tmp1 =
            (ZZ12::one().scale(2) - <ZZ12 as Units>::unit(2) - <ZZ12 as Units>::unit(-1).scale(2))
                * <ZZ12 as Units>::unit(1);
        let tmp2 =
            (ZZ12::one().scale(3) - <ZZ12 as Units>::unit(2) - <ZZ12 as Units>::unit(-1).scale(3))
                * <ZZ12 as Units>::unit(2);
        let p1 = (tmp1 - tmp2) * <ZZ12 as Units>::unit(-2);
        let p2 = p1 + <ZZ12 as Units>::unit(2);
        // p1 in cell (0, 0), p2 in cell (-1, -1) under floor semantics.
        assert_eq!(
            cells_of(p1, p2),
            vec![
                (-2, -1),
                (-1, -2),
                (-1, -1),
                (-1, 0),
                (0, -1),
                (0, 0),
                (0, 1),
                (1, 0),
            ]
        );

        // unit(0) = (1.0, 0.0); cell(0, 0) and cell(1, 0).
        assert_eq!(
            cells_of(ZZ12::zero(), <ZZ12 as Units>::unit(0)),
            vec![
                (-1, 0),
                (0, -1),
                (0, 0),
                (0, 1),
                (1, -1),
                (1, 0),
                (1, 1),
                (2, 0),
            ]
        );

        // unit(3) = (0.0, 1.0); cell(0, 0) and cell(0, 1).
        assert_eq!(
            cells_of(ZZ12::zero(), <ZZ12 as Units>::unit(3)),
            vec![
                (-1, 0),
                (-1, 1),
                (0, -1),
                (0, 0),
                (0, 1),
                (0, 2),
                (1, 0),
                (1, 1),
            ]
        );

        // unit(1) = (sqrt(3)/2, 1/2) ~= (0.866, 0.5); both endpoints in cell (0, 0).
        assert_eq!(
            cells_of(ZZ12::zero(), <ZZ12 as Units>::unit(1)),
            vec![(-1, 0), (0, -1), (0, 0), (0, 1), (1, 0)]
        );
    }

    /// Cross-check the precomputed `seg_neighborhood_of` table against the
    /// reference union of two 5-cell `+` crosses for every `(p1, dir)`
    /// pair on a coarse integer-coordinate grid spanning both axis-
    /// aligned and half-integer boundary positions. Any regression in
    /// the lookup table or in `cell_of`'s `{-1, 0, 1}^2` invariant
    /// would surface here.
    #[test]
    fn test_seg_neighborhood_of_matches_reference() {
        let mut diffs = 0usize;
        for ax in -3..=3i64 {
            for ay in -3..=3i64 {
                for bx in -3..=3i64 {
                    for by in -3..=3i64 {
                        for dir in 0..12i8 {
                            // p1 = ax + bx*sqrt(3)/2 + i*(ay + by*sqrt(3)/2)
                            // -- sweeps both integer and half-integer cell
                            // boundary positions.
                            let p1: ZZ12 = ZZ12::from((ax, ay))
                                + <ZZ12 as Units>::unit(2).scale(bx)
                                + <ZZ12 as Units>::unit(2).scale(by)
                                    * <ZZ12 as Units>::unit(3);
                            let p2 = p1 + <ZZ12 as Units>::unit(dir);
                            let mut from_table: Vec<(i64, i64)> =
                                UnitSquareGrid::seg_neighborhood_of(p1, p2).collect();
                            from_table.sort_unstable();
                            from_table.dedup();
                            // Reference: union of the two 5-cell `+` crosses
                            // centered on each endpoint's cell.
                            let cross = |p: ZZ12| -> [(i64, i64); 5] {
                                let (x, y) = UnitSquareGrid::cell_of(p);
                                [
                                    (x - 1, y),
                                    (x, y - 1),
                                    (x, y),
                                    (x, y + 1),
                                    (x + 1, y),
                                ]
                            };
                            let mut reference: Vec<(i64, i64)> = Vec::with_capacity(10);
                            reference.extend_from_slice(&cross(p1));
                            reference.extend_from_slice(&cross(p2));
                            reference.sort_unstable();
                            reference.dedup();
                            if from_table != reference {
                                if diffs < 5 {
                                    let c1 = UnitSquareGrid::cell_of(p1);
                                    let c2 = UnitSquareGrid::cell_of(p2);
                                    eprintln!(
                                        "DIFF at (ax,ay,bx,by,dir)=({ax},{ay},{bx},{by},{dir}): \
                                         cell(p1)={c1:?}, cell(p2)={c2:?}"
                                    );
                                    eprintln!("  from_table: {from_table:?}");
                                    eprintln!("  reference:  {reference:?}");
                                }
                                diffs += 1;
                            }
                        }
                    }
                }
            }
        }
        assert_eq!(diffs, 0, "{diffs} (p1, dir) cases disagree with the reference");
    }

    #[test]
    fn test_get_cells() {
        let mut grid = UnitSquareGrid::new();
        grid.add((0, 0), 0);
        grid.add((0, 0), 1);
        grid.add((1, 0), 2);
        grid.add((1, 0), 0);

        let empty: &[usize] = &[];
        assert_eq!(grid.get((0, 0)), &[0, 1]);
        assert_eq!(grid.get((1, 0)), &[2, 0]);
        assert_eq!(grid.get((2, 0)), empty);
        assert_eq!(grid.get_cells(&[(0, 0)]), &[0, 1]);
        assert_eq!(grid.get_cells(&[(1, 0)]), &[2, 0]);
        assert_eq!(grid.get_cells(&[(2, 0)]), empty);
        assert_eq!(grid.get_cells(&[(2, 0), (1, 0), (0, 0)]), &[2, 0, 0, 1]);
    }

    #[test]
    fn test_remove_basic() {
        let mut grid = UnitSquareGrid::new();
        grid.add((0, 0), 10);
        grid.add((0, 0), 20);
        grid.add((0, 0), 30);
        grid.add((1, 1), 40);

        assert_eq!(grid.len(), 4);

        grid.remove((0, 0), 20);
        assert_eq!(grid.get((0, 0)), &[10, 30]);
        assert_eq!(grid.len(), 3);

        grid.remove((0, 0), 10);
        assert_eq!(grid.get((0, 0)), &[30]);
        assert_eq!(grid.len(), 2);

        grid.remove((0, 0), 30);
        let empty: &[usize] = &[];
        assert_eq!(grid.get((0, 0)), empty);
        assert!(!grid.cells.contains_key(&(0, 0)));
        assert_eq!(grid.len(), 1);

        grid.remove((1, 1), 40);
        assert!(!grid.cells.contains_key(&(1, 1)));
        assert!(grid.is_empty());
    }

    #[test]
    fn test_remove_nonexistent() {
        let mut grid = UnitSquareGrid::new();
        grid.add((0, 0), 1);
        grid.remove((0, 0), 999);
        assert_eq!(grid.get((0, 0)), &[1]);
        assert_eq!(grid.len(), 1);

        grid.remove((5, 5), 1);
        assert_eq!(grid.get((0, 0)), &[1]);
    }

    #[test]
    fn test_remove_preserves_other_cells() {
        let mut grid = UnitSquareGrid::new();
        grid.add((0, 0), 1);
        grid.add((0, 0), 2);
        grid.add((1, 0), 3);
        grid.add((1, 1), 4);

        grid.remove((0, 0), 1);
        assert_eq!(grid.get((0, 0)), &[2]);
        assert_eq!(grid.get((1, 0)), &[3]);
        assert_eq!(grid.get((1, 1)), &[4]);
    }

    #[test]
    fn test_add_remove_roundtrip() {
        let mut grid = UnitSquareGrid::new();
        let cells = [(0, 0), (1, 0), (0, 1), (-1, 0), (0, -1), (2, 2), (-2, -2)];
        for (i, &cell) in cells.iter().enumerate() {
            grid.add(cell, i);
        }
        assert_eq!(grid.len(), 7);

        for (i, &cell) in cells.iter().enumerate() {
            grid.remove(cell, i);
        }
        assert!(grid.is_empty());
        assert!(grid.cells.is_empty());
    }

    #[test]
    fn test_multiple_values_same_cell() {
        let mut grid = UnitSquareGrid::new();
        for i in 0..10 {
            grid.add((3, 7), i);
        }
        assert_eq!(grid.len(), 10);
        assert_eq!(grid.get((3, 7)).len(), 10);

        for i in (0..10).rev() {
            grid.remove((3, 7), i);
        }
        assert!(grid.is_empty());
        assert!(!grid.cells.contains_key(&(3, 7)));
    }

    #[test]
    fn test_clone_independent() {
        let mut grid = UnitSquareGrid::new();
        grid.add((0, 0), 1);
        grid.add((1, 1), 2);

        let mut clone = grid.clone();
        clone.remove((0, 0), 1);
        clone.add((2, 2), 3);

        let empty: &[usize] = &[];
        assert_eq!(grid.get((0, 0)), &[1]);
        assert_eq!(grid.get((2, 2)), empty);
        assert_eq!(clone.get((0, 0)), empty);
        assert_eq!(clone.get((1, 1)), &[2]);
        assert_eq!(clone.get((2, 2)), &[3]);
    }
}
