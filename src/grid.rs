/// Unit square lattice utils
use super::zzbase::ZZNum;
use std::{collections::BTreeMap, ops::Bound};

#[derive(Debug, Clone)]
pub struct UnitSquareGrid {
    num_points: usize,
    cols_rows: BTreeMap<i64, BTreeMap<i64, Vec<usize>>>,
}

impl UnitSquareGrid {
    /// Return coordinates of closest unit square lattice cell
    /// that a complex integer falls in. The cells are centered
    /// at points that are pairs of two integers in the complex plane.
    pub fn cell_of<T: ZZNum>(zz: T) -> (i64, i64) {
        let c = zz.complex();
        (c.re.round() as i64, c.im.round() as i64)
    }

    /// Given a complex integer point, return the 5 unit square lattice cells
    /// that make up the neighborhood of the cell of the point.
    pub fn neighborhood_of<T: ZZNum>(zz: T) -> Vec<(i64, i64)> {
        let c @ (x, y) = Self::cell_of(zz);
        let r = (x + 1, y);
        let l = (x - 1, y);
        let u = (x, y + 1);
        let d = (x, y - 1);
        vec![l, d, c, u, r] // sorted by first x, then y
    }

    /// Given endpoints of a unit length line segment,
    /// returns the 5 or 8 distinct unit square lattice cells that:
    /// 1. the line segment is guaranteed to be fully contained inside, and
    /// 2. all intersecting unit length segments also have at least one point in it.
    pub fn seg_neighborhood_of<T: ZZNum>(p1: T, p2: T) -> Vec<(i64, i64)> {
        let mut result = Self::neighborhood_of(p1);
        result.extend(Self::neighborhood_of(p2));
        result.sort();
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
    pub fn get(&self, (x, y): (i64, i64)) -> Vec<usize> {
        match self.cols_rows.get(&x) {
            None => vec![],
            Some(col) => match col.get(&y) {
                None => vec![],
                Some(vals) => vals.clone(),
            },
        }
    }

    /// Returns indices of all points in the given unit square lattice
    /// that are located inside one of the given cells.
    ///
    /// Note that this does neither sort not deduplicate points.
    pub fn get_cells(&self, cells: &[(i64, i64)]) -> Vec<usize> {
        let mut pt_indices: Vec<usize> = Vec::new();
        for cell in cells {
            pt_indices.extend(self.get(*cell));
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
        let mut pt_indices: Vec<usize> = Vec::new();
        for (_, col) in self.cols_rows.range(x_range) {
            for (_, cell) in col.range(y_range) {
                pt_indices.extend(cell);
            }
        }
        pt_indices
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::zz::ZZ12;
    use crate::zzbase::ZZBase;
    use num_traits::{One, Zero};

    #[test]
    fn test_cell_of() {
        assert_eq!(UnitSquareGrid::cell_of(ZZ12::zero()), (0, 0));
        assert_eq!(UnitSquareGrid::cell_of(ZZ12::unit(0)), (1, 0));
        assert_eq!(UnitSquareGrid::cell_of(ZZ12::unit(1)), (1, 1));
        assert_eq!(UnitSquareGrid::cell_of(ZZ12::unit(2)), (1, 1));
        assert_eq!(UnitSquareGrid::cell_of(ZZ12::unit(3)), (0, 1));
        assert_eq!(UnitSquareGrid::cell_of(ZZ12::unit(4)), (-1, 1));
        assert_eq!(UnitSquareGrid::cell_of(ZZ12::unit(5)), (-1, 1));
        assert_eq!(UnitSquareGrid::cell_of(ZZ12::unit(6)), (-1, 0));
        assert_eq!(UnitSquareGrid::cell_of(ZZ12::unit(7)), (-1, -1));
        assert_eq!(UnitSquareGrid::cell_of(ZZ12::unit(8)), (-1, -1));
        assert_eq!(UnitSquareGrid::cell_of(ZZ12::unit(9)), (0, -1));
        assert_eq!(UnitSquareGrid::cell_of(ZZ12::unit(10)), (1, -1));
        assert_eq!(UnitSquareGrid::cell_of(ZZ12::unit(11)), (1, -1));
    }

    #[test]
    fn test_seg_neighborhood_of() {
        // segment fully contained in one cell
        let tmp1 = (ZZ12::one().scale(2) - ZZ12::unit(2) - ZZ12::unit(-1).scale(2)) * ZZ12::unit(1);
        let tmp2 = (ZZ12::one().scale(3) - ZZ12::unit(2) - ZZ12::unit(-1).scale(3)) * ZZ12::unit(2);
        let p1 = (tmp1 - tmp2) * ZZ12::unit(-2);
        let p2 = p1 + ZZ12::unit(2);
        let n0 = UnitSquareGrid::seg_neighborhood_of(p1, p2);
        let n0_exp = &[(-1, 0), (0, -1), (0, 0), (0, 1), (1, 0)];
        assert_eq!(n0, n0_exp);

        // segment touches 2 cells (horizontal)
        let n1 = UnitSquareGrid::seg_neighborhood_of(ZZ12::zero(), ZZ12::unit(0));
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
        let n2 = UnitSquareGrid::seg_neighborhood_of(ZZ12::zero(), ZZ12::unit(3));
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
        let n3 = UnitSquareGrid::seg_neighborhood_of(ZZ12::zero(), ZZ12::unit(1));
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
