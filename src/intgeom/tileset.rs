use crate::cyclotomic::{IsComplex, IsRingOrField, Units};
use crate::intgeom::rat::Rat;

pub struct TileSet<T: IsComplex> {
    rats: Vec<Rat<T>>,
}

impl<T: IsComplex + IsRingOrField + Units> TileSet<T> {
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

    pub fn num_tiles(&self) -> usize {
        self.rats.len()
    }

    pub fn rat(&self, i: usize) -> &Rat<T> {
        &self.rats[i]
    }

    pub fn rats(&self) -> &[Rat<T>] {
        &self.rats
    }

    #[allow(dead_code)]
    pub(crate) fn index_of(&self, rat: &Rat<T>) -> Option<usize> {
        self.rats.binary_search(rat).ok()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic::ZZ12;
    use crate::intgeom::tiles::{hexagon, spectre};

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
