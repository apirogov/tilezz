use std::collections::{BTreeMap, BTreeSet};
use std::marker::PhantomData;

use crate::cyclotomic::{IsComplex, IsRingOrField, Units};
use crate::intgeom::angles::normalize_angle;
use crate::intgeom::rat::{lex_min_rot, Rat};
use crate::intgeom::tileset::TileSet;

pub(crate) struct GlueSite {
    pub(crate) counter: usize,
    pub(crate) norm_start: usize,
    pub(crate) match_len: usize,
    pub(crate) norm_end: usize,
}

pub struct GrowingPatch<T: IsComplex> {
    angles: Vec<i8>,
    counters: Vec<usize>,
    next_counter: usize,
    _phantom: PhantomData<T>,
}

impl<T: IsComplex> GrowingPatch<T> {
    pub fn len(&self) -> usize {
        self.angles.len()
    }

    pub fn is_empty(&self) -> bool {
        self.angles.is_empty()
    }
}

impl<T: IsComplex + IsRingOrField + Units> GrowingPatch<T> {
    pub fn from_seed(seed: &Rat<T>) -> Self {
        let canonical = seed.clone().canonical();
        let seq = canonical.seq();
        let n = seq.len();
        GrowingPatch {
            angles: seq.to_vec(),
            counters: (0..n).collect(),
            next_counter: n,
            _phantom: PhantomData,
        }
    }

    pub fn to_rat(&self) -> Rat<T> {
        Rat::from_slice_unchecked(&self.angles)
    }

    pub(crate) fn enumerate_sites(&self, seed: &Rat<T>, min_counter: usize) -> Vec<GlueSite> {
        let patch_rat = self.to_rat();
        let ts = TileSet::new(vec![patch_rat, seed.clone()]);
        let glues = ts.valid_glues(0, 1);

        let mut sites: Vec<GlueSite> = Vec::new();
        for glue in glues {
            let ns = glue.start_a as usize;
            if self.counters[ns] < min_counter {
                continue;
            }
            sites.push(GlueSite {
                counter: self.counters[ns],
                norm_start: ns,
                match_len: glue.match_len,
                norm_end: glue.end_b as usize,
            });
        }
        sites.sort_by_key(|s| s.counter);
        sites
    }

    fn apply_site(&self, site: &GlueSite, seed: &Rat<T>) -> GrowingPatch<T> {
        let seed_seq = seed.seq();
        let m = seed_seq.len();
        let n = self.angles.len();
        let ml = site.match_len;
        let ns = site.norm_start;
        let ne = site.norm_end;

        let x_len = n - ml + 1;
        let y_len = m - ml + 1;

        let x: Vec<i8> = (0..x_len).map(|i| self.angles[(ns + ml + i) % n]).collect();
        let y: Vec<i8> = (0..y_len).map(|i| seed_seq[(ne + i) % m]).collect();

        let mut new_angles: Vec<i8> = Vec::with_capacity(x_len + y_len - 2);
        new_angles.extend_from_slice(&x[..x_len - 1]);
        new_angles.extend_from_slice(&y[..y_len - 1]);

        let a_yx = normalize_angle::<T>(x[0] + y[y_len - 1] - T::hturn());
        let a_xy = normalize_angle::<T>(y[0] + x[x_len - 1] - T::hturn());
        new_angles[0] = a_yx;
        new_angles[x_len - 1] = a_xy;

        let mut new_counters: Vec<usize> = Vec::with_capacity(new_angles.len());
        let mut nc = self.next_counter;
        for i in 0..(x_len - 1) {
            let old_pos = (ns + ml + i) % n;
            new_counters.push(self.counters[old_pos]);
        }
        for _ in 0..(y_len - 1) {
            new_counters.push(nc);
            nc += 1;
        }

        let offset = lex_min_rot(&new_angles);
        new_angles.rotate_left(offset);
        new_counters.rotate_left(offset);

        GrowingPatch {
            angles: new_angles,
            counters: new_counters,
            next_counter: nc,
            _phantom: PhantomData,
        }
    }
}

pub fn grow_redelmeier<T>(seed: &Rat<T>, max_size: usize) -> BTreeMap<usize, BTreeSet<Rat<T>>>
where
    T: IsComplex + IsRingOrField + Units,
{
    let mut results: BTreeMap<usize, BTreeSet<Rat<T>>> = BTreeMap::new();
    let initial = GrowingPatch::from_seed(seed);
    results.entry(1).or_default().insert(initial.to_rat());
    grow_recursive(&initial, seed, 0, 1, max_size, &mut results);
    results
}

fn grow_recursive<T>(
    patch: &GrowingPatch<T>,
    seed: &Rat<T>,
    min_counter: usize,
    current_size: usize,
    max_size: usize,
    results: &mut BTreeMap<usize, BTreeSet<Rat<T>>>,
) where
    T: IsComplex + IsRingOrField + Units,
{
    if current_size >= max_size {
        return;
    }
    let sites = patch.enumerate_sites(seed, min_counter);
    for site in &sites {
        let new_patch = patch.apply_site(site, seed);
        let rat = new_patch.to_rat();
        let new_size = current_size + 1;
        if results.entry(new_size).or_default().insert(rat) {
            grow_recursive(&new_patch, seed, site.counter, new_size, max_size, results);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic::{ZZ12, ZZ4};
    use crate::intgeom::tiles;

    #[test]
    fn hexagon_matches_old_approach_size4() {
        let seed: Rat<ZZ12> = Rat::from_unchecked(&tiles::hexagon());
        let new_results = grow_redelmeier(&seed, 4);

        let mut old_results: BTreeMap<usize, BTreeSet<Rat<ZZ12>>> = BTreeMap::new();
        old_results.insert(1, std::iter::once(seed.clone()).collect());

        for k in 2..=4 {
            let prev: Vec<Rat<ZZ12>> = old_results[&(k - 1)].iter().cloned().collect();
            let count_a = prev.len();
            let mut all_tiles = prev;
            all_tiles.push(seed.clone());
            let ts = crate::intgeom::tileset::TileSet::new(all_tiles);
            let pairs: Vec<(usize, usize)> = (0..count_a).map(|i| (i, count_a)).collect();
            let (results, _) = ts.valid_rats_for_pairs(&pairs);
            old_results.insert(k, results);
        }

        for k in 1..=4 {
            let old_count = old_results.get(&k).map(|s| s.len()).unwrap_or(0);
            let new_count = new_results.get(&k).map(|s| s.len()).unwrap_or(0);
            assert_eq!(
                old_count, new_count,
                "size {k}: old={old_count} new={new_count}"
            );
            if let (Some(old_set), Some(new_set)) = (old_results.get(&k), new_results.get(&k)) {
                assert_eq!(old_set, new_set, "size {k}: sets differ");
            }
        }
    }

    #[test]
    fn spectre_matches_old_approach_size3() {
        let seed: Rat<ZZ12> = Rat::from_unchecked(&tiles::spectre());
        let new_results = grow_redelmeier(&seed, 3);

        let mut old_results: BTreeMap<usize, BTreeSet<Rat<ZZ12>>> = BTreeMap::new();
        old_results.insert(1, std::iter::once(seed.clone()).collect());

        for k in 2..=3 {
            let prev: Vec<Rat<ZZ12>> = old_results[&(k - 1)].iter().cloned().collect();
            let count_a = prev.len();
            let mut all_tiles = prev;
            all_tiles.push(seed.clone());
            let ts = crate::intgeom::tileset::TileSet::new(all_tiles);
            let pairs: Vec<(usize, usize)> = (0..count_a).map(|i| (i, count_a)).collect();
            let (results, _) = ts.valid_rats_for_pairs(&pairs);
            old_results.insert(k, results);
        }

        for k in 1..=3 {
            let old_count = old_results.get(&k).map(|s| s.len()).unwrap_or(0);
            let new_count = new_results.get(&k).map(|s| s.len()).unwrap_or(0);
            assert_eq!(
                old_count, new_count,
                "size {k}: old={old_count} new={new_count}"
            );
            if let (Some(old_set), Some(new_set)) = (old_results.get(&k), new_results.get(&k)) {
                assert_eq!(old_set, new_set, "size {k}: sets differ");
            }
        }
    }

    #[test]
    fn square_zz4_matches_old_approach_size6() {
        let seed: Rat<ZZ4> = Rat::from_unchecked(&tiles::square());
        let new_results = grow_redelmeier(&seed, 6);

        let mut old_results: BTreeMap<usize, BTreeSet<Rat<ZZ4>>> = BTreeMap::new();
        old_results.insert(1, std::iter::once(seed.clone()).collect());

        for k in 2..=6 {
            let prev: Vec<Rat<ZZ4>> = old_results[&(k - 1)].iter().cloned().collect();
            let count_a = prev.len();
            let mut all_tiles = prev;
            all_tiles.push(seed.clone());
            let ts = crate::intgeom::tileset::TileSet::new(all_tiles);
            let pairs: Vec<(usize, usize)> = (0..count_a).map(|i| (i, count_a)).collect();
            let (results, _) = ts.valid_rats_for_pairs(&pairs);
            old_results.insert(k, results);
        }

        for k in 1..=6 {
            let old_count = old_results.get(&k).map(|s| s.len()).unwrap_or(0);
            let new_count = new_results.get(&k).map(|s| s.len()).unwrap_or(0);
            assert_eq!(
                old_count, new_count,
                "size {k}: old={old_count} new={new_count}"
            );
            if let (Some(old_set), Some(new_set)) = (old_results.get(&k), new_results.get(&k)) {
                assert_eq!(old_set, new_set, "size {k}: sets differ");
            }
        }
    }
}
