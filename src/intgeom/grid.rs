//! Utilities for managing information attached to a unit square lattice.

use std::{collections::BTreeMap, ops::Bound};

use crate::cyclotomic::SymNum;

#[derive(Debug, Clone)]
pub struct UnitSquareGrid {
    num_points: usize,
    cols_rows: BTreeMap<i64, BTreeMap<i64, Vec<usize>>>,
}

impl Default for UnitSquareGrid {
    fn default() -> Self {
        Self::new()
    }
}

impl UnitSquareGrid {
    /// Return coordinates of closest unit square lattice cell
    /// that a complex integer falls in. The cells are centered
    /// at points that are pairs of two integers in the complex plane.
    pub fn cell_of<T: SymNum>(zz: T) -> (i64, i64) {
        let c = zz.complex64();
        (c.re.round() as i64, c.im.round() as i64)
    }

    /// Given a complex integer point, return the 5 unit square lattice cells
    /// that make up the neighborhood of the cell of the point.
    pub fn neighborhood_of<T: SymNum>(zz: T) -> [(i64, i64); 5] {
        let (x, y) = Self::cell_of(zz);
        [(x - 1, y), (x, y - 1), (x, y), (x, y + 1), (x + 1, y)]
    }

    /// Given endpoints of a unit length line segment,
    /// returns the 5 or 8 distinct unit square lattice cells that:
    /// 1. the line segment is guaranteed to be fully contained inside, and
    /// 2. all intersecting unit length segments also have at least one point in it.
    pub fn seg_neighborhood_of<T: SymNum>(p1: T, p2: T) -> Vec<(i64, i64)> {
        let mut result = Vec::with_capacity(8);
        result.extend_from_slice(&Self::neighborhood_of(p1));
        result.extend_from_slice(&Self::neighborhood_of(p2));
        result.sort_unstable();
        result.dedup();
        result
    }

    // --------

    pub fn new() -> Self {
        Self {
            num_points: 0,
            cols_rows: BTreeMap::new(),
        }
    }

    /// Add a new index to the grid cell containing the given point
    pub fn add(&mut self, (x, y): (i64, i64), value: usize) {
        self.num_points += 1;
        self.cols_rows
            .entry(x)
            .or_default()
            .entry(y)
            .or_default()
            .push(value);
    }

    /// Get all values attached to the given grid cell.
    ///
    /// Returns a reference to the indices stored in this cell, or an empty slice.
    pub fn get(&self, (x, y): (i64, i64)) -> &[usize] {
        self.cols_rows
            .get(&x)
            .and_then(|col| col.get(&y))
            .map(Vec::as_slice)
            .unwrap_or(&[])
    }

    /// Returns indices of all points in the given unit square lattice
    /// that are located inside one of the given cells.
    ///
    /// Note that this does neither sort nor deduplicate points.
    pub fn get_cells(&self, cells: &[(i64, i64)]) -> Vec<usize> {
        let total_entries: usize = cells
            .iter()
            .map(|&(x, y)| {
                self.cols_rows
                    .get(&x)
                    .and_then(|col| col.get(&y))
                    .map(Vec::len)
                    .unwrap_or(0)
            })
            .sum();

        let mut pt_indices = Vec::with_capacity(total_entries);
        for &(x, y) in cells {
            if let Some(col) = self.cols_rows.get(&x) {
                if let Some(vals) = col.get(&y) {
                    pt_indices.extend(vals);
                }
            }
        }
        pt_indices
    }

    /// Returns indices of all points in all cells in the given range
    /// of grid coordinates.
    pub fn get_range(
        &self,
        x_range: (Bound<i64>, Bound<i64>),
        y_range: (Bound<i64>, Bound<i64>),
    ) -> Vec<usize> {
        let mut pt_indices = Vec::new();
        for (_, col) in self.cols_rows.range(x_range) {
            for (_, cell) in col.range(y_range) {
                pt_indices.reserve(cell.len());
                pt_indices.extend(cell);
            }
        }
        pt_indices
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic::{SymNum, Units, ZZ12};
    use num_traits::{One, Zero};

    #[test]
    fn test_cell_of() {
        assert_eq!(UnitSquareGrid::cell_of(ZZ12::zero()), (0, 0));
        assert_eq!(UnitSquareGrid::cell_of(<ZZ12 as Units>::unit(0)), (1, 0));
        assert_eq!(UnitSquareGrid::cell_of(<ZZ12 as Units>::unit(1)), (1, 1));
        assert_eq!(UnitSquareGrid::cell_of(<ZZ12 as Units>::unit(2)), (1, 1));
        assert_eq!(UnitSquareGrid::cell_of(<ZZ12 as Units>::unit(3)), (0, 1));
        assert_eq!(UnitSquareGrid::cell_of(<ZZ12 as Units>::unit(4)), (-1, 1));
        assert_eq!(UnitSquareGrid::cell_of(<ZZ12 as Units>::unit(5)), (-1, 1));
        assert_eq!(UnitSquareGrid::cell_of(<ZZ12 as Units>::unit(6)), (-1, 0));
        assert_eq!(UnitSquareGrid::cell_of(<ZZ12 as Units>::unit(7)), (-1, -1));
        assert_eq!(UnitSquareGrid::cell_of(<ZZ12 as Units>::unit(8)), (-1, -1));
        assert_eq!(UnitSquareGrid::cell_of(<ZZ12 as Units>::unit(9)), (0, -1));
        assert_eq!(UnitSquareGrid::cell_of(<ZZ12 as Units>::unit(10)), (1, -1));
        assert_eq!(UnitSquareGrid::cell_of(<ZZ12 as Units>::unit(11)), (1, -1));
    }

    #[test]
    fn test_seg_neighborhood_of() {
        // segment fully contained in one cell
        let tmp1 =
            (ZZ12::one().scale(2) - <ZZ12 as Units>::unit(2) - <ZZ12 as Units>::unit(-1).scale(2))
                * <ZZ12 as Units>::unit(1);
        let tmp2 =
            (ZZ12::one().scale(3) - <ZZ12 as Units>::unit(2) - <ZZ12 as Units>::unit(-1).scale(3))
                * <ZZ12 as Units>::unit(2);
        let p1 = (tmp1 - tmp2) * <ZZ12 as Units>::unit(-2);
        let p2 = p1 + <ZZ12 as Units>::unit(2);
        let n0 = UnitSquareGrid::seg_neighborhood_of(p1, p2);
        let n0_exp = &[(-1, 0), (0, -1), (0, 0), (0, 1), (1, 0)];
        assert_eq!(n0, n0_exp);

        // segment touches 2 cells (horizontal)
        let n1 = UnitSquareGrid::seg_neighborhood_of(ZZ12::zero(), <ZZ12 as Units>::unit(0));
        let n1_exp = &[
            (-1, 0),
            (0, -1),
            (0, 0),
            (0, 1),
            (1, -1),
            (1, 0),
            (1, 1),
            (2, 0),
        ];
        assert_eq!(n1, n1_exp);

        // segment touches 2 cells (vertical)
        let n2 = UnitSquareGrid::seg_neighborhood_of(ZZ12::zero(), <ZZ12 as Units>::unit(3));
        let n2_exp = &[
            (-1, 0),
            (-1, 1),
            (0, -1),
            (0, 0),
            (0, 1),
            (0, 2),
            (1, 0),
            (1, 1),
        ];
        assert_eq!(n2, n2_exp);

        // segment touches 2 cells (diagonal)
        let n3 = UnitSquareGrid::seg_neighborhood_of(ZZ12::zero(), <ZZ12 as Units>::unit(1));
        let n3_exp = &[
            (-1, 0),
            (0, -1),
            (0, 0),
            (0, 1),
            (1, 0),
            (1, 1),
            (1, 2),
            (2, 1),
        ];
        assert_eq!(n3, n3_exp);
    }

    #[test]
    fn test_indices_from_cells() {
        let mut grid = UnitSquareGrid::new();
        grid.add((0, 0), 0);
        grid.add((0, 0), 1);
        grid.add((1, 0), 2);
        grid.add((1, 0), 0);

        // sanity-check
        assert_eq!(grid.get((0, 0)), &[0, 1]);
        assert_eq!(grid.get((1, 0)), &[2, 0]);
        assert_eq!(grid.get((2, 0)), &[]);
        assert_eq!(grid.get_cells(&[(0, 0)]), &[0, 1]);
        assert_eq!(grid.get_cells(&[(1, 0)]), &[2, 0]);
        assert_eq!(UnitSquareGrid::get_cells(&grid, &[(2, 0)]), &[]);

        // combined results are simply concatenated
        assert_eq!(
            UnitSquareGrid::get_cells(&grid, &[(2, 0), (1, 0), (0, 0)]),
            &[2, 0, 0, 1]
        );
    }
}
