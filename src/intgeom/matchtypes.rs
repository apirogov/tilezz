use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

use crate::cyclotomic::{IsComplex, IsRingOrField, Units};
use crate::intgeom::rat::Rat;
use crate::intgeom::snake::Snake;
use crate::intgeom::tileset::TileSet;
use crate::stringmatch::repetition_factor;
use crate::stringmatch::CyclicMatchIndex;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct MatchType {
    pub tile_a: usize,
    pub start_a: usize,
    pub tile_b: usize,
    pub start_b: usize,
    pub len: usize,
}

impl MatchType {
    pub fn apply<T: IsComplex + IsRingOrField + Units>(
        &self,
        rats_a: &[Rat<T>],
        rats_b: &[Rat<T>],
    ) -> Rat<T> {
        rats_a[self.tile_a]
            .try_glue_precomputed(
                (self.start_a as i64, self.len, self.start_b as i64),
                &rats_b[self.tile_b],
                true,
            )
            .expect("match was pre-validated")
    }

    fn involution(self, len_a: usize, len_b: usize) -> Self {
        MatchType {
            tile_a: self.tile_b,
            start_a: (self.start_b + len_b - self.len) % len_b,
            tile_b: self.tile_a,
            start_b: (self.start_a + self.len) % len_a,
            len: self.len,
        }
    }

    fn normalize_by_periods(self, periods: &[usize]) -> Self {
        MatchType {
            tile_a: self.tile_a,
            start_a: self.start_a % periods[self.tile_a],
            tile_b: self.tile_b,
            start_b: self.start_b % periods[self.tile_b],
            len: self.len,
        }
    }
}

struct TypeEntry {
    first: MatchType,
    second: MatchType,
}

pub struct MatchFinder<T: IsComplex> {
    set_a: Arc<TileSet<T>>,
    set_b: Arc<TileSet<T>>,
    offset_b: usize,
    cmi: CyclicMatchIndex,
}

impl<T: IsComplex + IsRingOrField + Units> MatchFinder<T> {
    pub fn new(tileset: Arc<TileSet<T>>) -> Self {
        let sequences: Vec<Vec<i8>> = tileset.rats().iter().map(|r| r.seq().to_vec()).collect();
        let cmi = CyclicMatchIndex::new(&sequences);
        MatchFinder {
            set_a: Arc::clone(&tileset),
            set_b: tileset,
            offset_b: 0,
            cmi,
        }
    }

    pub fn crossing(a: Arc<TileSet<T>>, b: Arc<TileSet<T>>) -> Self {
        let sequences: Vec<Vec<i8>> = a
            .rats()
            .iter()
            .chain(b.rats().iter())
            .map(|r| r.seq().to_vec())
            .collect();
        let cmi = CyclicMatchIndex::new(&sequences);
        let offset_b = a.num_tiles();
        MatchFinder {
            set_a: a,
            set_b: b,
            offset_b,
            cmi,
        }
    }

    pub fn set_a(&self) -> &Arc<TileSet<T>> {
        &self.set_a
    }

    pub fn set_b(&self) -> &Arc<TileSet<T>> {
        &self.set_b
    }

    pub fn rat_a(&self, i: usize) -> &Rat<T> {
        self.set_a.rat(i)
    }

    pub fn rat_b(&self, j: usize) -> &Rat<T> {
        self.set_b.rat(j)
    }

    pub fn num_tiles_a(&self) -> usize {
        self.set_a.num_tiles()
    }

    pub fn num_tiles_b(&self) -> usize {
        self.set_b.num_tiles()
    }

    pub fn apply_match(&self, m: &MatchType) -> Rat<T> {
        m.apply(self.set_a.rats(), self.set_b.rats())
    }

    pub fn shared_boundaries(&self, i: usize, j: usize) -> Vec<crate::stringmatch::CyclicMatch> {
        self.cmi.maximal_rc_matches(i, self.offset_b + j)
    }

    pub fn valid_matches(&self, i: usize, j: usize) -> Vec<MatchType> {
        let mut m = Self::validated_matches(self.candidates_for_pair(i, j));
        Self::sort_by_interval(&mut m);
        m
    }

    pub fn all_valid_matches(&self) -> Vec<MatchType> {
        self.valid_matches_for_pairs(&self.all_pairs())
    }

    pub fn valid_matches_for_pairs(&self, pairs: &[(usize, usize)]) -> Vec<MatchType> {
        let mut m = Self::validated_matches(self.collect_candidates(pairs));
        Self::sort_global(&mut m);
        m
    }

    pub fn valid_results_for_pairs(&self, pairs: &[(usize, usize)]) -> BTreeSet<Rat<T>> {
        Self::validated_results(self.collect_candidates(pairs))
    }

    fn candidates_for_pair(&self, i: usize, j: usize) -> BTreeMap<Rat<T>, Vec<MatchType>> {
        let a = self.set_a.rat(i);
        let b = self.set_b.rat(j);
        let cmi_j = self.offset_b + j;
        let n_a = a.len();
        let n_b = b.len();

        let mut groups: BTreeMap<Rat<T>, Vec<MatchType>> = BTreeMap::new();
        if n_a == 0 || n_b == 0 {
            return groups;
        }

        let cmi_matches = self.cmi.maximal_rc_matches(i, cmi_j);
        for m in &cmi_matches {
            let (ns, len, ne) = a.get_match((m.pos_a as i64, m.pos_b as i64), b);
            if len <= 1 {
                continue;
            }
            if !junction_gap_nonnegative(a.seq(), ns as usize, len, b.seq(), ne as usize) {
                continue;
            }
            if let Ok(glued) = a.try_glue_precomputed((ns, len, ne), b, true) {
                groups.entry(glued).or_default().push(MatchType {
                    tile_a: i,
                    start_a: ns as usize,
                    tile_b: j,
                    start_b: ne as usize,
                    len,
                });
            }
        }

        let seq_a = a.seq();
        let seq_b = b.seq();
        for ia in 0..n_a {
            for ib in 0..n_b {
                if !is_single_edge_candidate(seq_a, ia, seq_b, ib) {
                    continue;
                }
                let (ns, len, ne) = a.get_match((ia as i64, ib as i64), b);
                if len != 1 {
                    continue;
                }
                if let Ok(glued) = a.try_glue_precomputed((ns, len, ne), b, true) {
                    groups.entry(glued).or_default().push(MatchType {
                        tile_a: i,
                        start_a: ns as usize,
                        tile_b: j,
                        start_b: ne as usize,
                        len,
                    });
                }
            }
        }

        groups
    }

    fn all_pairs(&self) -> Vec<(usize, usize)> {
        let na = self.num_tiles_a();
        let nb = self.num_tiles_b();
        (0..na).flat_map(|i| (0..nb).map(move |j| (i, j))).collect()
    }

    fn collect_candidates(&self, pairs: &[(usize, usize)]) -> BTreeMap<Rat<T>, Vec<MatchType>> {
        let mut all: BTreeMap<Rat<T>, Vec<MatchType>> = BTreeMap::new();
        for &(i, j) in pairs {
            for (rat, matches) in self.candidates_for_pair(i, j) {
                all.entry(rat).or_default().extend(matches);
            }
        }
        all
    }

    fn validated_matches(groups: BTreeMap<Rat<T>, Vec<MatchType>>) -> Vec<MatchType> {
        groups
            .into_iter()
            .filter(|(rat, _)| Snake::<T>::try_from(rat.seq()).is_ok())
            .flat_map(|(_, matches)| matches)
            .collect()
    }

    fn validated_results(groups: BTreeMap<Rat<T>, Vec<MatchType>>) -> BTreeSet<Rat<T>> {
        groups
            .into_iter()
            .filter(|(rat, _)| Snake::<T>::try_from(rat.seq()).is_ok())
            .map(|(rat, _)| rat)
            .collect()
    }

    fn sort_by_interval(results: &mut [MatchType]) {
        results.sort_by(|a, b| {
            a.start_a
                .cmp(&b.start_a)
                .then_with(|| a.start_b.cmp(&b.start_b))
        });
    }

    fn sort_global(results: &mut [MatchType]) {
        results.sort_by(|a, b| {
            a.tile_a
                .cmp(&b.tile_a)
                .then_with(|| a.tile_b.cmp(&b.tile_b))
                .then_with(|| a.start_a.cmp(&b.start_a))
                .then_with(|| a.start_b.cmp(&b.start_b))
        });
    }
}

pub struct MatchTypeIndex<T: IsComplex> {
    tileset: Arc<TileSet<T>>,
    entries: Vec<TypeEntry>,
    periods: Vec<usize>,
}

impl<T: IsComplex + IsRingOrField + Units> MatchTypeIndex<T> {
    pub fn new(tileset: Arc<TileSet<T>>) -> Self {
        let periods: Vec<usize> = tileset
            .rats()
            .iter()
            .map(|r| {
                let seq = r.seq();
                seq.len() / repetition_factor(seq)
            })
            .collect();

        let mf = MatchFinder::new(Arc::clone(&tileset));
        let all_matches = mf.all_valid_matches();
        let rats = tileset.rats();

        let mut seen: BTreeSet<(MatchType, MatchType)> = BTreeSet::new();

        for mt in &all_matches {
            let len_a = rats[mt.tile_a].len();
            let len_b = rats[mt.tile_b].len();

            let fwd = mt.normalize_by_periods(&periods);
            let bwd = mt.involution(len_a, len_b).normalize_by_periods(&periods);

            let (first, second) = if fwd <= bwd { (fwd, bwd) } else { (bwd, fwd) };

            seen.insert((first, second));
        }

        let entries: Vec<TypeEntry> = seen
            .into_iter()
            .map(|(first, second)| TypeEntry { first, second })
            .collect();

        MatchTypeIndex {
            tileset,
            entries,
            periods,
        }
    }

    pub fn tileset(&self) -> &Arc<TileSet<T>> {
        &self.tileset
    }

    pub fn num_types(&self) -> usize {
        self.entries.len()
    }

    pub fn get(&self, unsigned_id: usize) -> MatchType {
        assert!(
            unsigned_id >= 1 && unsigned_id <= self.entries.len(),
            "match type id {unsigned_id} out of range [1, {}]",
            self.entries.len()
        );
        self.entries[unsigned_id - 1].first
    }

    pub fn period(&self, tile_id: usize) -> usize {
        self.periods[tile_id]
    }

    pub fn signed_id(&self, mt: &MatchType) -> Option<i32> {
        let norm = mt.normalize_by_periods(&self.periods);
        for (idx, entry) in self.entries.iter().enumerate() {
            if norm == entry.first {
                return Some((idx + 1) as i32);
            }
            if norm == entry.second {
                return Some(-((idx + 1) as i32));
            }
        }
        None
    }

    pub fn apply(&self, signed_id: i32) -> Rat<T> {
        let unsigned = signed_id.unsigned_abs() as usize;
        assert!(
            unsigned >= 1 && unsigned <= self.entries.len(),
            "signed match type id {signed_id} out of range"
        );
        let entry = &self.entries[unsigned - 1];
        let rats = self.tileset.rats();
        let mt = if signed_id > 0 {
            entry.first
        } else {
            entry.second
        };
        mt.apply(rats, rats)
    }
}

fn is_single_edge_candidate(a: &[i8], ia: usize, b: &[i8], ib: usize) -> bool {
    let na = a.len();
    let nb = b.len();
    let left = a[ia] as i32 + b[ib] as i32;
    let right = a[(ia + 1) % na] as i32 + b[(ib + nb - 1) % nb] as i32;
    left > 0 && right > 0
}

fn junction_gap_nonnegative(a: &[i8], ns: usize, mlen: usize, b: &[i8], ne: usize) -> bool {
    let na = a.len();
    let nb = b.len();
    let left = a[(ns + mlen) % na] as i32 + b[(ne + nb - mlen) % nb] as i32;
    let right = a[ns] as i32 + b[ne] as i32;
    left >= 0 && right >= 0
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic::ZZ12;
    use crate::intgeom::tiles::{hexagon, spectre, square};
    use std::collections::{BTreeSet, HashSet};

    const MYSTIC_ZZ12: &[i8] = &[
        0, 2, -3, 2, 3, -2, 3, -2, 3, 2, -3, 2, 0, 2, -3, 2, 3, 2, -3, 2,
    ];

    fn brute_force_valid_matches<T: IsComplex + IsRingOrField + Units>(
        a: &Rat<T>,
        b: &Rat<T>,
        i: usize,
        j: usize,
    ) -> Vec<MatchType> {
        let mut seen: HashSet<(usize, usize, usize)> = HashSet::new();
        let mut results: Vec<MatchType> = Vec::new();
        for ia in 0..a.len() {
            for ib in 0..b.len() {
                let (ns, len, ne) = a.get_match((ia as i64, ib as i64), b);
                if len == 0 || !seen.insert((ns as usize, len, ne as usize)) {
                    continue;
                }
                if a.try_glue_precomputed((ns, len, ne), b, false).is_ok() {
                    results.push(MatchType {
                        tile_a: i,
                        start_a: ns as usize,
                        tile_b: j,
                        start_b: ne as usize,
                        len,
                    });
                }
            }
        }
        results.sort_by(|x, y| {
            x.start_a
                .cmp(&y.start_a)
                .then_with(|| x.start_b.cmp(&y.start_b))
        });
        results
    }

    fn assert_match_sets_match<T: IsComplex + IsRingOrField + Units>(
        mf_matches: &[MatchType],
        bf_matches: &[MatchType],
        rats_a: &[Rat<T>],
        rats_b: &[Rat<T>],
        label: &str,
    ) {
        assert_eq!(
            mf_matches.len(),
            bf_matches.len(),
            "{label}: count mismatch: MatchFinder={}, brute_force={}",
            mf_matches.len(),
            bf_matches.len(),
        );
        for (idx, (mf_m, bf_m)) in mf_matches.iter().zip(bf_matches.iter()).enumerate() {
            assert_eq!(
                (mf_m.start_a, mf_m.len, mf_m.start_b),
                (bf_m.start_a, bf_m.len, bf_m.start_b),
                "{label}: interval mismatch at index {idx}",
            );
            assert_eq!(
                mf_m.apply(rats_a, rats_b),
                bf_m.apply(rats_a, rats_b),
                "{label}: result mismatch at ({}, {})",
                mf_m.start_a,
                mf_m.start_b,
            );
        }
    }

    fn self_match_ts(rats: Vec<Rat<ZZ12>>) -> MatchFinder<ZZ12> {
        let ts = Arc::new(TileSet::new(rats));
        MatchFinder::new(ts)
    }

    #[test]
    fn spectre_self_finds_mystic() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&spectre());
        let mf = self_match_ts(vec![r]);

        let matches = mf.valid_matches(0, 0);
        assert!(
            matches
                .iter()
                .any(|m| mf.apply_match(m).seq() == MYSTIC_ZZ12),
            "mystic not found in {} matches: {:?}",
            matches.len(),
            matches,
        );
    }

    #[test]
    fn hexagon_self_single_edge() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&hexagon());
        let mf = self_match_ts(vec![r]);

        let matches = mf.valid_matches(0, 0);
        assert!(
            !matches.is_empty(),
            "hexagon should have valid self-matches"
        );

        let expected = &[-2, 2, 2, 2, 2, -2, 2, 2, 2, 2];
        assert!(
            matches.iter().any(|m| mf.apply_match(m).seq() == expected),
            "hex+hex match not found in {} matches",
            matches.len(),
        );
    }

    #[test]
    fn self_intersecting_filtered() {
        let r1 = Rat::<ZZ12>::from_slice_unchecked(&[0, 0, 3, 0, 3, 3, -3, -3, 3, 3, 0, 3]);
        let r2 = Rat::<ZZ12>::from_slice_unchecked(&[0, 0, 3, 3, 0, 0, 3, 3]);

        assert!(r1.try_glue((8, 0), &r2).is_err());

        let mf = self_match_ts(vec![r1, r2]);
        let matches = mf.valid_matches(0, 1);
        for m in &matches {
            assert!(
                !(m.start_a == 8 && m.start_b == 0),
                "self-intersecting glue should be filtered"
            );
        }
    }

    #[test]
    fn vertex_touching_rejected() {
        let r3 = Rat::<ZZ12>::try_from(
            &Snake::try_from(&[
                0, 0, 3, 2, 4, -3, -3, 3, 2, 4, -3, -3, 3, 0, 2, 4, -3, -3, 3, 2, 4, -3, 3, -3,
            ])
            .unwrap(),
        )
        .unwrap();
        let r4 = Rat::<ZZ12>::from_slice_unchecked(&[0, 0, 3, 3, 0, 0, 3, 3]);

        assert!(r3.try_glue((7, 0), &r4).is_err());

        let mf = self_match_ts(vec![r3, r4]);
        let matches = mf.valid_matches(0, 1);
        for m in &matches {
            assert!(
                !(m.start_a == 7 && m.start_b == 0),
                "vertex-touching glue should be filtered"
            );
        }
    }

    #[test]
    fn deterministic_order() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&hexagon());
        let mf = self_match_ts(vec![r]);
        let matches = mf.valid_matches(0, 0);
        for w in matches.windows(2) {
            assert!(
                w[0].start_a <= w[1].start_a
                    && (w[0].start_a < w[1].start_a || w[0].start_b <= w[1].start_b),
                "matches not in deterministic order"
            );
        }
    }

    #[test]
    fn multi_tile_all_valid_matches() {
        let sq = Rat::<ZZ12>::from_slice_unchecked(&[0, 1, 0, 1, 0, 1, 0, 1]);
        let tri = Rat::<ZZ12>::from_slice_unchecked(&[4, 4, 4]);

        let mf = self_match_ts(vec![sq, tri]);
        let all = mf.all_valid_matches();
        for w in all.windows(2) {
            assert!(
                w[0].tile_a < w[1].tile_a
                    || (w[0].tile_a == w[1].tile_a && w[0].tile_b < w[1].tile_b)
                    || (w[0].tile_a == w[1].tile_a
                        && w[0].tile_b == w[1].tile_b
                        && (w[0].start_a < w[1].start_a
                            || (w[0].start_a == w[1].start_a && w[0].start_b <= w[1].start_b))),
                "all_valid_matches not sorted"
            );
        }
    }

    #[test]
    fn spectre_foldback_cases_accepted() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&spectre());
        let cases: Vec<(usize, usize)> = vec![(0, 4), (0, 8), (6, 4), (6, 8)];
        for (ia, ib) in cases {
            assert!(
                r.try_glue((ia as i64, ib as i64), &r).is_ok(),
                "spectre foldback ({ia},{ib}) should be accepted by try_glue",
            );
        }
    }

    #[test]
    fn hexagon_pair_exhaustive() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&hexagon());
        let mf = self_match_ts(vec![r.clone()]);

        let mf_matches = mf.valid_matches(0, 0);
        let bf_matches = brute_force_valid_matches(&r, &r, 0, 0);
        assert_match_sets_match(
            &mf_matches,
            &bf_matches,
            mf.set_a().rats(),
            mf.set_b().rats(),
            "hexagon self",
        );

        for m in &mf_matches {
            assert_eq!(m.len, 1, "hexagon pair match should be single-edge");
            assert_eq!(
                mf.apply_match(m).len(),
                10,
                "hex+hex result should be 10-gon"
            );
        }

        let canonical = mf.apply_match(&mf_matches[0]);
        for (idx, m) in mf_matches.iter().enumerate() {
            assert_eq!(
                mf.apply_match(m),
                canonical,
                "hex+hex match #{idx} should yield same shape",
            );
        }

        assert_eq!(
            mf_matches.len(),
            36,
            "hexagon pair should have 36 valid matches"
        );
    }

    #[test]
    fn hexamino_pair_exhaustive() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&hexagon());
        let hexamino = r.glue((0, 0), &r);

        let mf = self_match_ts(vec![hexamino.clone()]);
        let mf_matches = mf.valid_matches(0, 0);
        let bf_matches = brute_force_valid_matches(&hexamino, &hexamino, 0, 0);
        assert_match_sets_match(
            &mf_matches,
            &bf_matches,
            mf.set_a().rats(),
            mf.set_b().rats(),
            "hexamino self",
        );

        assert!(
            !mf_matches.is_empty(),
            "hexamino pair should have valid matches"
        );

        let multi_count = mf_matches.iter().filter(|m| m.len > 1).count();
        let single_count = mf_matches.iter().filter(|m| m.len == 1).count();
        assert!(
            multi_count > 0,
            "hexamino pair should have multi-edge matches",
        );
        assert!(
            single_count > 0,
            "hexamino pair should have single-edge matches",
        );
    }

    #[test]
    fn spectre_pair_exhaustive() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&spectre());
        let mf = self_match_ts(vec![r.clone()]);

        let mf_matches = mf.valid_matches(0, 0);
        let bf_matches = brute_force_valid_matches(&r, &r, 0, 0);
        assert_match_sets_match(
            &mf_matches,
            &bf_matches,
            mf.set_a().rats(),
            mf.set_b().rats(),
            "spectre self",
        );

        assert!(
            mf_matches
                .iter()
                .any(|m| mf.apply_match(m).seq() == MYSTIC_ZZ12),
            "spectre pair should find mystic, got {} matches",
            mf_matches.len(),
        );

        let multi_count = mf_matches.iter().filter(|m| m.len > 1).count();
        let single_count = mf_matches.iter().filter(|m| m.len == 1).count();
        assert!(
            multi_count > 0,
            "spectre pair should have multi-edge matches"
        );
        assert!(
            single_count > 0,
            "spectre pair should have single-edge matches"
        );
    }

    #[test]
    fn mixed_shapes_exhaustive() {
        let sq = Rat::<ZZ12>::from_slice_unchecked(&[0, 1, 0, 1, 0, 1, 0, 1]);
        let tri = Rat::<ZZ12>::from_slice_unchecked(&[4, 4, 4]);
        let mf = self_match_ts(vec![sq, tri]);

        for i in 0..2 {
            for j in 0..2 {
                let mf_matches = mf.valid_matches(i, j);
                let bf_matches = brute_force_valid_matches(mf.rat_a(i), mf.rat_b(j), i, j);
                assert_match_sets_match(
                    &mf_matches,
                    &bf_matches,
                    mf.set_a().rats(),
                    mf.set_b().rats(),
                    &format!("mixed shapes ({i},{j})"),
                );
            }
        }
    }

    #[test]
    fn single_edge_candidate_unit() {
        let hex: &[i8] = &[2, 2, 2, 2, 2, 2];
        assert!(
            is_single_edge_candidate(hex, 0, hex, 0),
            "hex (0,0): 2+2=4 > 0 on both sides"
        );
        assert!(
            is_single_edge_candidate(hex, 0, hex, 3),
            "hex (0,3): 2+2=4 > 0 on both sides"
        );

        let hexamino: &[i8] = &[-2, 2, 2, 2, 2, -2, 2, 2, 2, 2];
        assert!(
            !is_single_edge_candidate(hexamino, 0, hexamino, 0),
            "hexamino (0,0): left=-2+(-2)=-4 < 0, overflow"
        );
        assert!(
            !is_single_edge_candidate(hexamino, 1, hexamino, 1),
            "hexamino (1,1): right=2+(-2)=0, collinear extension"
        );
        assert!(
            is_single_edge_candidate(hexamino, 1, hexamino, 2),
            "hexamino (1,2): left=2+2=4, right=2+2=4, both gaps"
        );

        let spectre: &[i8] = &[3, 2, 0, 2, -3, 2, 3, 2, -3, 2, 3, -2, 3, -2];
        assert!(
            !is_single_edge_candidate(spectre, 3, spectre, 4),
            "spectre (3,4): left=2+(-3)=-1 < 0, overflow"
        );
        assert!(
            is_single_edge_candidate(spectre, 0, spectre, 1),
            "spectre (0,1): left=3+2=5, right=2+3=5, both gaps"
        );
    }

    #[test]
    fn interior_angle_heuristic_sound() {
        let hex = Rat::<ZZ12>::from_unchecked(&hexagon());
        let spec = Rat::<ZZ12>::from_unchecked(&spectre());
        let hexamino = hex.glue((0, 0), &hex);

        let tiles: Vec<(&str, &Rat<ZZ12>)> =
            vec![("hex", &hex), ("spec", &spec), ("hexamino", &hexamino)];

        for (label, rat) in &tiles {
            let n = rat.len();
            let seq = rat.seq();
            for ia in 0..n {
                for ib in 0..n {
                    let (_, len, _) = rat.get_match((ia as i64, ib as i64), rat);
                    if len != 1 {
                        continue;
                    }
                    if rat.try_glue((ia as i64, ib as i64), rat).is_ok() {
                        assert!(
                            is_single_edge_candidate(seq, ia, seq, ib),
                            "{label}: try_glue accepted single-edge ({ia},{ib}) but heuristic rejected",
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn empty_tile_valid_matches() {
        let hex = Rat::<ZZ12>::from_unchecked(&hexagon());
        let mf = self_match_ts(vec![hex]);
        assert!(!mf.valid_matches(0, 0).is_empty());
    }

    #[test]
    fn shared_boundaries_returns_cmi_matches() {
        let h: Snake<ZZ12> = hexagon();
        let r: Rat<ZZ12> = Rat::from_unchecked(&h);
        let mf = self_match_ts(vec![r]);

        let boundaries = mf.shared_boundaries(0, 0);
        assert!(
            boundaries.is_empty(),
            "hexagon self should have no RC matches: {:?}",
            boundaries,
        );

        let s: Snake<ZZ12> = spectre();
        let spec: Rat<ZZ12> = Rat::from_unchecked(&s);
        let mf2 = self_match_ts(vec![spec]);
        let self_bounds = mf2.shared_boundaries(0, 0);
        assert!(
            !self_bounds.is_empty(),
            "spectre self should have RC matches",
        );
        let max_len = self_bounds.iter().map(|m| m.len).max().unwrap();
        assert_eq!(max_len, 3, "spectre max RC self-match should be 3");
    }

    #[test]
    fn heuristic_sound_on_size4_hexagon_patches() {
        let hex: Rat<ZZ12> = Rat::from_unchecked(&hexagon());

        let ts1 = Arc::new(TileSet::new(vec![hex.clone()]));
        let mf1 = MatchFinder::new(ts1);
        let size2: Vec<Rat<ZZ12>> = mf1
            .all_valid_matches()
            .iter()
            .map(|m| mf1.apply_match(m))
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();

        let ts2 = Arc::new(TileSet::new(size2.clone()));
        let mf2 = MatchFinder::new(ts2);
        let size4: Vec<Rat<ZZ12>> = mf2
            .all_valid_matches()
            .iter()
            .map(|m| mf2.apply_match(m))
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();

        let mut false_rejects = 0;

        for a in &size4 {
            for b in &size4 {
                let seq_a = a.seq();
                let seq_b = b.seq();
                for ia in 0..a.len() {
                    for ib in 0..b.len() {
                        let (_, len, _) = a.get_match((ia as i64, ib as i64), b);
                        if len != 1 {
                            continue;
                        }
                        if a.try_glue((ia as i64, ib as i64), b).is_err() {
                            continue;
                        }
                        if !is_single_edge_candidate(seq_a, ia, seq_b, ib) {
                            false_rejects += 1;
                        }
                    }
                }
            }
        }

        assert_eq!(
            false_rejects, 0,
            "heuristic falsely rejected {false_rejects} valid single-edge glues"
        );
    }

    #[test]
    fn size8_hexagon_valid_matches_matches_brute_force() {
        let hex: Rat<ZZ12> = Rat::from_unchecked(&hexagon());

        let ts1 = Arc::new(TileSet::new(vec![hex.clone()]));
        let mf1 = MatchFinder::new(ts1);
        let size2: Vec<Rat<ZZ12>> = mf1
            .all_valid_matches()
            .iter()
            .map(|m| mf1.apply_match(m))
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();

        let ts2 = Arc::new(TileSet::new(size2.clone()));
        let mf2 = MatchFinder::new(ts2);
        let size4: Vec<Rat<ZZ12>> = mf2
            .all_valid_matches()
            .iter()
            .map(|m| mf2.apply_match(m))
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();

        let ts4 = Arc::new(TileSet::new(size4.clone()));
        let mf4 = MatchFinder::new(ts4);

        let mut index_results: BTreeSet<Rat<ZZ12>> = BTreeSet::new();
        for i in 0..mf4.num_tiles_a() {
            for j in 0..mf4.num_tiles_b() {
                for m in mf4.valid_matches(i, j) {
                    index_results.insert(mf4.apply_match(&m));
                }
            }
        }

        let mut bf_results: BTreeSet<Rat<ZZ12>> = BTreeSet::new();
        for a in &size4 {
            for b in &size4 {
                for ia in 0..a.len() {
                    for ib in 0..b.len() {
                        if let Ok(glued) = a.try_glue((ia as i64, ib as i64), b) {
                            bf_results.insert(glued);
                        }
                    }
                }
            }
        }

        eprintln!(
            "index: {}, brute_force: {}",
            index_results.len(),
            bf_results.len(),
        );

        let missing: Vec<_> = bf_results.difference(&index_results).collect();
        assert_eq!(missing.len(), 0, "index missed {} patches", missing.len(),);
    }

    #[test]
    fn crossing_constructor_distinct_sets() {
        let hex: Rat<ZZ12> = Rat::from_unchecked(&hexagon());
        let spec: Rat<ZZ12> = Rat::from_unchecked(&spectre());

        let ts_a = Arc::new(TileSet::new(vec![hex.clone()]));
        let ts_b = Arc::new(TileSet::new(vec![spec.clone()]));
        let mf = MatchFinder::crossing(ts_a, ts_b);

        assert_eq!(mf.num_tiles_a(), 1);
        assert_eq!(mf.num_tiles_b(), 1);

        let matches = mf.valid_matches(0, 0);
        assert!(
            !matches.is_empty(),
            "hex x spectre crossing should have valid matches"
        );

        for m in &matches {
            assert_eq!(m.tile_a, 0);
            assert_eq!(m.tile_b, 0);
        }
    }

    #[test]
    fn crossing_self_matches_new_equivalent() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&hexagon());

        let ts = Arc::new(TileSet::new(vec![r.clone()]));
        let mf_self = MatchFinder::new(Arc::clone(&ts));
        let mf_cross = MatchFinder::crossing(Arc::clone(&ts), ts);

        let self_matches = mf_self.valid_matches(0, 0);
        let cross_matches = mf_cross.valid_matches(0, 0);

        assert_eq!(
            self_matches.len(),
            cross_matches.len(),
            "self and crossing should produce same number of matches"
        );

        for (s, c) in self_matches.iter().zip(cross_matches.iter()) {
            assert_eq!(
                (s.start_a, s.len, s.start_b),
                (c.start_a, c.len, c.start_b),
                "intervals should match"
            );
        }
    }

    #[test]
    fn crossing_multi_patch_against_seed() {
        let hex: Rat<ZZ12> = Rat::from_unchecked(&hexagon());
        let hex2 = hex.glue((0, 0), &hex);

        let ts_patches = Arc::new(TileSet::new(vec![hex.clone(), hex2]));
        let ts_seed = Arc::new(TileSet::new(vec![hex]));
        let mf = MatchFinder::crossing(ts_patches, ts_seed);

        assert_eq!(mf.num_tiles_a(), 2);
        assert_eq!(mf.num_tiles_b(), 1);

        let all = mf.all_valid_matches();
        assert!(
            !all.is_empty(),
            "multi-patch x seed should have valid matches"
        );

        for m in &all {
            assert!(m.tile_a < 2, "tile_a should index into patches");
            assert_eq!(m.tile_b, 0, "tile_b should index into seed");
        }
    }

    // --- MatchTypeIndex tests ---

    #[test]
    fn hexagon_self_one_type() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&hexagon());
        let ts = Arc::new(TileSet::new(vec![r]));
        let idx = MatchTypeIndex::new(ts);

        assert_eq!(idx.period(0), 1, "hexagon period should be 1");
        assert_eq!(idx.num_types(), 1, "hexagon self should have 1 match type");

        let mt = idx.get(1);
        assert_eq!(mt.tile_a, 0);
        assert_eq!(mt.tile_b, 0);
        assert_eq!(mt.start_a, 0);
        assert_eq!(mt.start_b, 0);
        assert_eq!(mt.len, 1);
    }

    #[test]
    fn square_self_one_type() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&square());
        let ts = Arc::new(TileSet::new(vec![r]));
        let idx = MatchTypeIndex::new(ts);

        assert_eq!(idx.period(0), 1, "square period should be 1");
        assert_eq!(idx.num_types(), 1, "square self should have 1 match type");

        let mt = idx.get(1);
        assert_eq!((mt.tile_a, mt.tile_b), (0, 0));
        assert_eq!((mt.start_a, mt.start_b), (0, 0));
        assert_eq!(mt.len, 1);
    }

    #[test]
    fn spectre_self_many_types() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&spectre());
        let ts = Arc::new(TileSet::new(vec![r]));
        let idx = MatchTypeIndex::new(ts);

        assert_eq!(idx.period(0), 14, "spectre period should be 14");
        assert!(
            idx.num_types() > 10,
            "spectre self should have many match types, got {}",
            idx.num_types()
        );

        for id in 1..=idx.num_types() {
            let mt = idx.get(id);
            assert!(mt.start_a < 14, "start_a should be < period");
            assert!(mt.start_b < 14, "start_b should be < period");
        }
    }

    #[test]
    fn signed_id_forward_and_backward() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&spectre());
        let ts = Arc::new(TileSet::new(vec![r]));
        let idx = MatchTypeIndex::new(ts);

        let mut found_asymmetric = false;
        for id in 1..=idx.num_types() {
            let mt = idx.get(id);
            if mt.start_a == mt.start_b {
                continue;
            }
            let len_a = idx.tileset().rat(mt.tile_a).len();
            let len_b = idx.tileset().rat(mt.tile_b).len();
            let inv = mt.involution(len_a, len_b);
            let sid_fwd = idx.signed_id(&mt).unwrap();
            let sid_inv = idx.signed_id(&inv).unwrap();
            assert_eq!(sid_fwd.abs(), sid_inv.abs(), "same unsigned id");
            if inv == mt {
                assert_eq!(sid_fwd, sid_inv, "self-involutive: both same sign");
                assert!(sid_fwd > 0, "should be positive");
            } else {
                assert_eq!(sid_fwd, -sid_inv, "asymmetric: opposite signs");
                assert!(sid_fwd > 0, "first should be positive");
                found_asymmetric = true;
            }
        }
        assert!(
            found_asymmetric,
            "should find at least one asymmetric match type"
        );
    }

    #[test]
    fn rotational_equivalence_same_signed_id() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&hexagon());
        let ts = Arc::new(TileSet::new(vec![r]));
        let idx = MatchTypeIndex::new(ts);

        let mut seen_ids: HashSet<i32> = HashSet::new();
        for ia in 0..6 {
            for ib in 0..6 {
                let mt = MatchType {
                    tile_a: 0,
                    start_a: ia,
                    tile_b: 0,
                    start_b: ib,
                    len: 1,
                };
                if let Some(id) = idx.signed_id(&mt) {
                    seen_ids.insert(id);
                }
            }
        }
        assert_eq!(
            seen_ids.len(),
            1,
            "all hex self-matches should produce same signed id"
        );
    }

    #[test]
    fn apply_produces_valid_rat() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&spectre());
        let ts = Arc::new(TileSet::new(vec![r]));
        let idx = MatchTypeIndex::new(ts);

        for id in 1..=idx.num_types() {
            let result_pos = idx.apply(id as i32);
            assert!(
                Snake::<ZZ12>::try_from(result_pos.seq()).is_ok(),
                "apply(+{id}) should produce valid snake"
            );
            let result_neg = idx.apply(-(id as i32));
            assert!(
                Snake::<ZZ12>::try_from(result_neg.seq()).is_ok(),
                "apply(-{id}) should produce valid snake"
            );
            assert_eq!(
                result_pos, result_neg,
                "forward and backward apply should produce same rat"
            );
        }
    }

    #[test]
    fn bi_hex_bi_square_cross_matches() {
        let hex: Rat<ZZ12> = Rat::from_unchecked(&hexagon());
        let sq: Rat<ZZ12> = Rat::from_unchecked(&square());
        let bi_hex = hex.glue((0, 0), &hex);
        let bi_sq = sq.glue((0, 0), &sq);

        let ts = Arc::new(TileSet::new(vec![bi_hex.clone(), bi_sq.clone()]));
        let idx = MatchTypeIndex::new(Arc::clone(&ts));

        assert_eq!(ts.num_tiles(), 2);
        assert_eq!(idx.period(0), 5, "bi-hex period should be 5");
        assert_eq!(idx.period(1), 3, "bi-sq period should be 3");

        let n = idx.num_types();
        assert!(n > 0, "should have match types");

        let mut has_self_hex = false;
        let mut has_self_sq = false;
        let mut has_cross = false;
        for id in 1..=n {
            let mt = idx.get(id);
            if mt.tile_a == 0 && mt.tile_b == 0 {
                has_self_hex = true;
            }
            if mt.tile_a == 1 && mt.tile_b == 1 {
                has_self_sq = true;
            }
            if mt.tile_a != mt.tile_b {
                has_cross = true;
            }
        }
        assert!(has_self_hex, "should have bi-hex self-match types");
        assert!(has_self_sq, "should have bi-sq self-match types");
        assert!(has_cross, "should have cross-match types");

        for id in 1..=n {
            let mt = idx.get(id);
            assert!(
                mt.start_a < idx.period(mt.tile_a),
                "start_a {} should be < period {}",
                mt.start_a,
                idx.period(mt.tile_a)
            );
            assert!(
                mt.start_b < idx.period(mt.tile_b),
                "start_b {} should be < period {}",
                mt.start_b,
                idx.period(mt.tile_b)
            );
        }
    }

    #[test]
    fn signed_id_roundtrip_bi_hex_bi_square() {
        let hex: Rat<ZZ12> = Rat::from_unchecked(&hexagon());
        let sq: Rat<ZZ12> = Rat::from_unchecked(&square());
        let bi_hex = hex.glue((0, 0), &hex);
        let bi_sq = sq.glue((0, 0), &sq);

        let ts = Arc::new(TileSet::new(vec![bi_hex, bi_sq]));
        let idx = MatchTypeIndex::new(Arc::clone(&ts));

        let mf = MatchFinder::new(Arc::clone(&ts));
        let all = mf.all_valid_matches();

        let mut signed_ids: HashSet<i32> = HashSet::new();
        for mt in &all {
            let sid = idx
                .signed_id(mt)
                .unwrap_or_else(|| panic!("match {:?} not found in index", mt));
            signed_ids.insert(sid);

            let result_mf = mf.apply_match(mt);
            let result_idx = idx.apply(sid);
            assert_eq!(
                result_mf, result_idx,
                "MatchFinder and MatchTypeIndex should produce same rat for {:?}",
                mt
            );
        }

        assert!(
            signed_ids.len() >= idx.num_types(),
            "should have at least one signed id per type"
        );

        for &sid in &signed_ids {
            let unsigned = sid.unsigned_abs() as usize;
            assert!(unsigned >= 1 && unsigned <= idx.num_types());
        }
    }

    #[test]
    fn get_match_involution_semantics() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&spectre());
        let n = r.len();

        let mf = self_match_ts(vec![r.clone()]);
        let matches = mf.valid_matches(0, 0);

        let mut checked = 0;
        for mt in &matches {
            if mt.start_a == mt.start_b {
                continue;
            }

            let fwd = (mt.start_a as i64, mt.len, mt.start_b as i64);
            let rat_fwd = r
                .try_glue_precomputed(fwd, &r, true)
                .expect("forward glue should work");

            let ns = mt.start_a as i64;
            let ne = mt.start_b as i64;
            let len = mt.len as i64;

            let candidates: Vec<(&str, (i64, i64))> = vec![
                ("(ne, ns)", (ne, ns)),
                ("(ne, ns+len)", (ne, ns + len)),
                ("(ne, ns-len)", (ne, ns - len)),
                ("(ne-len, ns)", (ne - len, ns)),
                ("(ne-len, ns+len)", (ne - len, ns + len)),
                ("(ne+len, ns)", (ne + len, ns)),
                ("(ne+len, ns-len)", (ne + len, ns - len)),
            ];

            let mut found = false;
            for (label, (s, e)) in &candidates {
                let (ns2, len2, ne2) = r.get_match((*s, *e), &r);
                if len2 != mt.len {
                    continue;
                }
                if let Ok(rat_rev) = r.try_glue_precomputed((ns2, len2, ne2), &r, true) {
                    if rat_rev == rat_fwd {
                        eprintln!(
                            "  MATCH {} works: fwd=({},{},{}) rev=({},{},{}) rat_fwd==rat_rev",
                            label, ns, mt.len, ne, ns2, len2, ne2
                        );
                        found = true;

                        assert_eq!(
                            ns2,
                            (ne - len).rem_euclid(n as i64),
                            "reverse ext_start should be (ne - len) % n for {}",
                            label
                        );
                        assert_eq!(
                            ne2,
                            (ns + len).rem_euclid(n as i64),
                            "reverse ext_end should be (ns + len) % n for {}",
                            label
                        );
                    }
                }
            }

            assert!(
                found,
                "no reverse candidate produced same rat for fwd=({ns},{len},{ne})"
            );

            checked += 1;
            if checked >= 5 {
                break;
            }
        }
        assert!(
            checked > 0,
            "need at least one match with start_a != start_b"
        );
    }

    #[test]
    fn hexagon_self_apply_both_signs_same_rat() {
        let hex: Rat<ZZ12> = Rat::from_unchecked(&hexagon());
        let ts = Arc::new(TileSet::new(vec![hex]));
        let idx = MatchTypeIndex::new(ts);

        assert_eq!(idx.num_types(), 1, "hex self should have exactly 1 type");

        let r_pos = idx.apply(1);
        let r_neg = idx.apply(-1);
        assert!(
            Snake::<ZZ12>::try_from(r_pos.seq()).is_ok(),
            "apply(+1) should produce valid snake"
        );
        assert!(
            Snake::<ZZ12>::try_from(r_neg.seq()).is_ok(),
            "apply(-1) should produce valid snake"
        );
        assert_eq!(
            r_pos, r_neg,
            "apply(+1) and apply(-1) should produce same rat"
        );
    }

    #[test]
    fn hex_square_cross_apply_both_signs_same_rat() {
        let hex: Rat<ZZ12> = Rat::from_unchecked(&hexagon());
        let sq: Rat<ZZ12> = Rat::from_unchecked(&square());
        let ts = Arc::new(TileSet::new(vec![hex, sq]));
        let idx = MatchTypeIndex::new(ts);

        let cross_ids: Vec<usize> = (1..=idx.num_types())
            .filter(|id| {
                let mt = idx.get(*id);
                mt.tile_a != mt.tile_b
            })
            .collect();

        assert!(
            !cross_ids.is_empty(),
            "hex+square should have cross-match types"
        );

        for id in cross_ids {
            let r_pos = idx.apply(id as i32);
            let r_neg = idx.apply(-(id as i32));
            assert!(
                Snake::<ZZ12>::try_from(r_pos.seq()).is_ok(),
                "apply(+{id}) should produce valid snake"
            );
            assert!(
                Snake::<ZZ12>::try_from(r_neg.seq()).is_ok(),
                "apply(-{id}) should produce valid snake"
            );
            assert_eq!(
                r_pos, r_neg,
                "apply(+{id}) and apply(-{id}) should produce same rat"
            );
        }
    }

    #[test]
    fn hex_bisq_match_type_count() {
        let hex: Rat<ZZ12> = Rat::from_unchecked(&hexagon());
        let sq: Rat<ZZ12> = Rat::from_unchecked(&square());
        let bi_sq = sq.glue((0, 0), &sq);
        let ts = Arc::new(TileSet::new(vec![hex, bi_sq]));
        let idx = MatchTypeIndex::new(Arc::clone(&ts));

        let hex_idx = ts
            .index_of(&Rat::<ZZ12>::from_unchecked(&hexagon()))
            .unwrap();
        let bisq_idx = ts.index_of(&sq.glue((0, 0), &sq)).unwrap();
        assert_ne!(hex_idx, bisq_idx, "tiles should be distinct");

        let hex_period = idx.period(hex_idx);
        let bisq_period = idx.period(bisq_idx);
        assert_eq!(hex_period, 1, "hex period should be 1");
        assert_eq!(bisq_period, 3, "bisq period should be 3");

        let mut hex_self = 0usize;
        let mut bisq_self = 0usize;
        let mut cross = 0usize;
        for id in 1..=idx.num_types() {
            let mt = idx.get(id);
            if mt.tile_a == hex_idx && mt.tile_b == hex_idx {
                hex_self += 1;
            } else if mt.tile_a == bisq_idx && mt.tile_b == bisq_idx {
                bisq_self += 1;
            } else {
                cross += 1;
            }

            let r_pos = idx.apply(id as i32);
            let r_neg = idx.apply(-(id as i32));
            assert!(
                Snake::<ZZ12>::try_from(r_pos.seq()).is_ok(),
                "apply(+{id}) should produce valid snake"
            );
            assert!(
                Snake::<ZZ12>::try_from(r_neg.seq()).is_ok(),
                "apply(-{id}) should produce valid snake"
            );
            assert_eq!(
                r_pos, r_neg,
                "apply(+{id}) and apply(-{id}) should produce same rat"
            );
        }

        assert_eq!(hex_self, 1, "should have exactly 1 hex self-match type");
        assert_eq!(cross, 3, "should have exactly 3 cross-match types");
        assert_eq!(bisq_self, 6, "should have exactly 6 bisq self-match types");
    }
}
