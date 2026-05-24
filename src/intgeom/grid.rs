use rustc_hash::FxHashMap;

use crate::cyclotomic::SymNum;

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
    pub fn cell_of<T: SymNum>(zz: T) -> (i64, i64) {
        let c = zz.complex64();
        (c.re.round() as i64, c.im.round() as i64)
    }

    pub fn neighborhood_of<T: SymNum>(zz: T) -> [(i64, i64); 5] {
        let (x, y) = Self::cell_of(zz);
        [(x - 1, y), (x, y - 1), (x, y), (x, y + 1), (x + 1, y)]
    }

    pub fn seg_neighborhood_of<T: SymNum>(p1: T, p2: T) -> Vec<(i64, i64)> {
        let mut result = Vec::with_capacity(8);
        result.extend_from_slice(&Self::neighborhood_of(p1));
        result.extend_from_slice(&Self::neighborhood_of(p2));
        result.sort_unstable();
        result.dedup();
        result
    }

    pub fn new() -> Self {
        Self {
            cells: FxHashMap::default(),
        }
    }

    pub fn len(&self) -> usize {
        self.cells.values().map(|v| v.len()).sum()
    }

    pub fn is_empty(&self) -> bool {
        self.cells.values().all(|v| v.is_empty())
    }

    pub fn add(&mut self, cell: (i64, i64), value: usize) {
        self.cells.entry(cell).or_default().push(value);
    }

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

    pub fn get(&self, cell: (i64, i64)) -> &[usize] {
        self.cells.get(&cell).map(Vec::as_slice).unwrap_or(&[])
    }

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
        let tmp1 =
            (ZZ12::one().scale(2) - <ZZ12 as Units>::unit(2) - <ZZ12 as Units>::unit(-1).scale(2))
                * <ZZ12 as Units>::unit(1);
        let tmp2 =
            (ZZ12::one().scale(3) - <ZZ12 as Units>::unit(2) - <ZZ12 as Units>::unit(-1).scale(3))
                * <ZZ12 as Units>::unit(2);
        let p1 = (tmp1 - tmp2) * <ZZ12 as Units>::unit(-2);
        let p2 = p1 + <ZZ12 as Units>::unit(2);
        let n0 = UnitSquareGrid::seg_neighborhood_of(p1, p2);
        assert_eq!(n0, &[(-1, 0), (0, -1), (0, 0), (0, 1), (1, 0)]);

        let n1 = UnitSquareGrid::seg_neighborhood_of(ZZ12::zero(), <ZZ12 as Units>::unit(0));
        assert_eq!(
            n1,
            &[
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

        let n2 = UnitSquareGrid::seg_neighborhood_of(ZZ12::zero(), <ZZ12 as Units>::unit(3));
        assert_eq!(
            n2,
            &[
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

        let n3 = UnitSquareGrid::seg_neighborhood_of(ZZ12::zero(), <ZZ12 as Units>::unit(1));
        assert_eq!(
            n3,
            &[
                (-1, 0),
                (0, -1),
                (0, 0),
                (0, 1),
                (1, 0),
                (1, 1),
                (1, 2),
                (2, 1),
            ]
        );
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
