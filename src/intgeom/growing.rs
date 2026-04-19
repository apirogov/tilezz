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

#[derive(Default)]
pub struct GrowStats {
    pub enumerate_calls: usize,
    pub enumerate_ns: u64,
    pub to_rat_ns: u64,
    pub tileset_new_ns: u64,
    pub valid_glues_ns: u64,
    pub apply_ns: u64,
    pub apply_to_rat_ns: u64,
}

impl std::fmt::Display for GrowStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let total = self.enumerate_ns + self.apply_ns;
        write!(
            f,
            "enumerate: {} calls, {:.2}s total ({:.0}% to_rat, {:.0}% tileset_new, {:.0}% valid_glues) | apply: {:.2}s ({:.0}% to_rat)",
            self.enumerate_calls,
            self.enumerate_ns as f64 / 1e9,
            if self.enumerate_ns > 0 { self.to_rat_ns as f64 / self.enumerate_ns as f64 * 100.0 } else { 0.0 },
            if self.enumerate_ns > 0 { self.tileset_new_ns as f64 / self.enumerate_ns as f64 * 100.0 } else { 0.0 },
            if self.enumerate_ns > 0 { self.valid_glues_ns as f64 / self.enumerate_ns as f64 * 100.0 } else { 0.0 },
            self.apply_ns as f64 / 1e9,
            if self.apply_ns > 0 { self.apply_to_rat_ns as f64 / self.apply_ns as f64 * 100.0 } else { 0.0 },
        )
    }
}

pub fn grow_redelmeier<T>(seed: &Rat<T>, max_size: usize) -> BTreeMap<usize, BTreeSet<Rat<T>>>
where
    T: IsComplex + IsRingOrField + Units,
{
    grow_redelmeier_profiled(seed, max_size).0
}

pub fn grow_redelmeier_profiled<T>(
    seed: &Rat<T>,
    max_size: usize,
) -> (BTreeMap<usize, BTreeSet<Rat<T>>>, GrowStats)
where
    T: IsComplex + IsRingOrField + Units,
{
    use std::time::Instant;
    let mut stats = GrowStats::default();
    let mut results: BTreeMap<usize, BTreeSet<Rat<T>>> = BTreeMap::new();
    let initial = GrowingPatch::from_seed(seed);
    results.entry(1).or_default().insert(initial.to_rat());
    grow_recursive_profiled(&initial, seed, 0, 1, max_size, &mut results, &mut stats);
    (results, stats)
}

fn grow_recursive_profiled<T>(
    patch: &GrowingPatch<T>,
    seed: &Rat<T>,
    min_counter: usize,
    current_size: usize,
    max_size: usize,
    results: &mut BTreeMap<usize, BTreeSet<Rat<T>>>,
    stats: &mut GrowStats,
) where
    T: IsComplex + IsRingOrField + Units,
{
    use std::time::Instant;
    if current_size >= max_size {
        return;
    }

    let sites = {
        stats.enumerate_calls += 1;
        let t0 = Instant::now();

        let t1 = Instant::now();
        let patch_rat = patch.to_rat();
        stats.to_rat_ns += t1.elapsed().as_nanos() as u64;

        let t1 = Instant::now();
        let ts = TileSet::new(vec![patch_rat, seed.clone()]);
        stats.tileset_new_ns += t1.elapsed().as_nanos() as u64;

        let t1 = Instant::now();
        let glues = ts.valid_glues(0, 1);
        stats.valid_glues_ns += t1.elapsed().as_nanos() as u64;

        let mut sites: Vec<GlueSite> = Vec::new();
        for glue in glues {
            let ns = glue.start_a as usize;
            if patch.counters[ns] < min_counter {
                continue;
            }
            sites.push(GlueSite {
                counter: patch.counters[ns],
                norm_start: ns,
                match_len: glue.match_len,
                norm_end: glue.end_b as usize,
            });
        }
        sites.sort_by_key(|s| s.counter);
        stats.enumerate_ns += t0.elapsed().as_nanos() as u64;
        sites
    };

    for site in &sites {
        let t0 = Instant::now();
        let new_patch = patch.apply_site(site, seed);
        let t1 = Instant::now();
        let rat = new_patch.to_rat();
        stats.apply_to_rat_ns += t1.elapsed().as_nanos() as u64;
        stats.apply_ns += t0.elapsed().as_nanos() as u64;

        let new_size = current_size + 1;
        if results.entry(new_size).or_default().insert(rat) {
            grow_recursive_profiled(
                &new_patch,
                seed,
                site.counter,
                new_size,
                max_size,
                results,
                stats,
            );
        }
    }
}

pub fn make_free<T>(onesided: &BTreeSet<Rat<T>>) -> BTreeSet<Rat<T>>
where
    T: IsComplex + IsRingOrField + Units,
{
    onesided
        .iter()
        .map(|r| std::cmp::min(r.clone(), r.reflected()))
        .collect()
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

    fn brute_force_grow<T: IsComplex + IsRingOrField + Units>(
        seed: &Rat<T>,
        max_size: usize,
    ) -> BTreeMap<usize, BTreeSet<Rat<T>>> {
        let mut results: BTreeMap<usize, BTreeSet<Rat<T>>> = BTreeMap::new();
        results.insert(1, std::iter::once(seed.clone()).collect());

        for k in 2..=max_size {
            let prev: Vec<Rat<T>> = results[&(k - 1)].iter().cloned().collect();
            let mut next: BTreeSet<Rat<T>> = BTreeSet::new();
            for patch in &prev {
                for ia in 0..patch.len() {
                    for ib in 0..seed.len() {
                        if let Ok(glued) = patch.try_glue((ia as i64, ib as i64), seed) {
                            next.insert(glued);
                        }
                    }
                }
            }
            results.insert(k, next);
        }
        results
    }

    #[test]
    fn square_zz4_brute_force_vs_tileset_size7() {
        let seed: Rat<ZZ4> = Rat::from_unchecked(&tiles::square());
        let bf = brute_force_grow(&seed, 7);

        eprintln!("=== Brute force (try_glue) ===");
        for k in 1..=7 {
            let count = bf.get(&k).map(|s| s.len()).unwrap_or(0);
            eprintln!("  size {k}: {count}");
        }

        let redel = grow_redelmeier(&seed, 7);
        eprintln!("=== Redelmeier (TileSet) ===");
        for k in 1..=7 {
            let count = redel.get(&k).map(|s| s.len()).unwrap_or(0);
            eprintln!("  size {k}: {count}");
        }

        for k in 1..=7 {
            let bf_count = bf.get(&k).map(|s| s.len()).unwrap_or(0);
            let redel_count = redel.get(&k).map(|s| s.len()).unwrap_or(0);
            assert_eq!(
                bf_count, redel_count,
                "size {k}: brute_force={bf_count} redelmeier={redel_count}"
            );
            if let (Some(bf_set), Some(redel_set)) = (bf.get(&k), redel.get(&k)) {
                let missing: Vec<_> = bf_set.difference(redel_set).collect();
                if !missing.is_empty() {
                    eprintln!(
                        "  size {k}: {} patches in brute_force but not redelmeier:",
                        missing.len()
                    );
                    for m in &missing {
                        eprintln!(
                            "    {}",
                            m.seq()
                                .iter()
                                .map(|a| a.to_string())
                                .collect::<Vec<_>>()
                                .join(", ")
                        );
                    }
                }
                assert_eq!(bf_set, redel_set, "size {k}: sets differ");
            }
        }
    }

    const FREE_POLYOMINOES_NO_HOLES: &[usize] = &[
        1, 1, 1, 2, 5, 12, 35, 107, 363, 1248, 4460, 16094, 58937, 217117,
    ];

    #[test]
    fn square_zz4_free_polyominoes_match_oeis() {
        let max_size = 9;
        let seed: Rat<ZZ4> = Rat::from_unchecked(&tiles::square());
        let onesided = grow_redelmeier(&seed, max_size);

        for k in 1..=max_size {
            let onesided_set = onesided.get(&k).unwrap();
            let free_set = make_free(onesided_set);
            let free_count = free_set.len();
            let expected = FREE_POLYOMINOES_NO_HOLES[k];
            assert_eq!(
                free_count, expected,
                "size {k}: got {free_count}, expected {expected} free polyominoes (no holes)"
            );
        }
    }
}
