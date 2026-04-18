use std::collections::HashSet;
use std::ops::Range;

use crate::cyclotomic::{IsComplex, IsRingOrField, Units};
use crate::intgeom::rat::Rat;
use crate::stringmatch::CyclicMatchIndex;

pub struct GlueOp<T: IsComplex> {
    pub tile_a: usize,
    pub tile_b: usize,
    pub start_a: i64,
    pub end_b: i64,
    pub match_len: usize,
    pub result: Rat<T>,
}

impl<T: IsComplex + IsRingOrField + Units> std::fmt::Debug for GlueOp<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "GlueOp(a={}, b={}, start_a={}, end_b={}, len={}, result={})",
            self.tile_a, self.tile_b, self.start_a, self.end_b, self.match_len, self.result,
        )
    }
}

pub struct TileSet<T: IsComplex> {
    rats: Vec<Rat<T>>,
    cmi: CyclicMatchIndex,
}

impl<T: IsComplex + IsRingOrField + Units> TileSet<T> {
    pub fn new(rats: Vec<Rat<T>>) -> Self {
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

        let sequences: Vec<Vec<i8>> = rats.iter().map(|r| r.seq().to_vec()).collect();
        let cmi = CyclicMatchIndex::new(&sequences);

        TileSet { rats, cmi }
    }

    pub fn num_tiles(&self) -> usize {
        self.rats.len()
    }

    pub fn rat(&self, i: usize) -> &Rat<T> {
        &self.rats[i]
    }

    pub fn shared_boundaries(&self, i: usize, j: usize) -> Vec<crate::stringmatch::CyclicMatch> {
        self.cmi.maximal_rc_matches(i, j)
    }

    fn get_new_match(
        a: &Rat<T>,
        b: &Rat<T>,
        ia: i64,
        ib: i64,
        seen: &mut HashSet<(i64, usize, i64)>,
    ) -> Option<(i64, usize, i64)> {
        let (ns, len, ne) = a.get_match((ia, ib), b);
        if len > 0 && seen.insert((ns, len, ne)) {
            Some((ns, len, ne))
        } else {
            None
        }
    }

    fn build_glue_op(
        a: &Rat<T>,
        b: &Rat<T>,
        ia: i64,
        ib: i64,
        i: usize,
        j: usize,
        (ns, len, ne): (i64, usize, i64),
    ) -> Option<GlueOp<T>> {
        a.try_glue((ia, ib), b).ok().map(|glued| GlueOp {
            tile_a: i,
            tile_b: j,
            start_a: ns,
            end_b: ne,
            match_len: len,
            result: glued,
        })
    }

    fn try_add_glue(
        a: &Rat<T>,
        b: &Rat<T>,
        ia: i64,
        ib: i64,
        i: usize,
        j: usize,
        seen: &mut HashSet<(i64, usize, i64)>,
    ) -> Option<GlueOp<T>> {
        let match_info = Self::get_new_match(a, b, ia, ib, seen)?;
        Self::build_glue_op(a, b, ia, ib, i, j, match_info)
    }

    fn sort_by_interval(results: &mut [GlueOp<T>]) {
        results.sort_by(|a, b| {
            a.start_a
                .cmp(&b.start_a)
                .then_with(|| a.end_b.cmp(&b.end_b))
        });
    }

    /// Find all valid glue operations between tiles `i` and `j`.
    ///
    /// Uses a two-phase approach:
    ///
    /// **Phase 1 (multi-edge):** Seeds from the CMI index find maximal
    /// reverse-complement matches. Each CMI match position is also shifted
    /// by one edge to catch boundary extensions where `Rat::get_match`
    /// extends beyond the CMI match (since `Rat` starts unconditionally at
    /// length 1 while CMI requires angle equality at every position).
    ///
    /// **Phase 2 (single-edge):** Enumerates all position pairs and applies
    /// the interior-angle overflow heuristic (`is_single_edge_candidate`)
    /// to reject locally inconsistent candidates before the global
    /// self-intersection check via `try_glue`.
    ///
    /// Results are sorted by `(start_a, end_b)` in deterministic order.
    pub fn valid_glues(&self, i: usize, j: usize) -> Vec<GlueOp<T>> {
        let a = &self.rats[i];
        let b = &self.rats[j];
        let n_a = a.len();
        let n_b = b.len();

        if n_a == 0 || n_b == 0 {
            return vec![];
        }

        let mut seen: HashSet<(i64, usize, i64)> = HashSet::new();
        let mut results: Vec<GlueOp<T>> = Vec::new();

        let cmi_matches = self.cmi.maximal_rc_matches(i, j);
        for m in &cmi_matches {
            let pa = m.pos_a as i64;
            let pb = m.pos_b as i64;
            results.extend(Self::try_add_glue(a, b, pa, pb, i, j, &mut seen));
            results.extend(Self::try_add_glue(
                a,
                b,
                (pa - 1).rem_euclid(n_a as i64),
                (pb + 1).rem_euclid(n_b as i64),
                i,
                j,
                &mut seen,
            ));
        }

        let seq_a = a.seq();
        let seq_b = b.seq();
        for ia in 0..n_a {
            for ib in 0..n_b {
                if !is_single_edge_candidate(seq_a, ia, seq_b, ib) {
                    continue;
                }
                results.extend(Self::try_add_glue(
                    a, b, ia as i64, ib as i64, i, j, &mut seen,
                ));
            }
        }

        Self::sort_by_interval(&mut results);
        results
    }

    pub fn all_valid_glues(&self) -> Vec<GlueOp<T>> {
        let mut results = Vec::new();
        for i in 0..self.rats.len() {
            for j in i..self.rats.len() {
                results.extend(self.valid_glues(i, j));
            }
        }
        results.sort_by(|a, b| {
            a.tile_a
                .cmp(&b.tile_a)
                .then_with(|| a.tile_b.cmp(&b.tile_b))
                .then_with(|| a.start_a.cmp(&b.start_a))
                .then_with(|| a.end_b.cmp(&b.end_b))
        });
        results
    }

    pub fn valid_glues_containing(
        &self,
        i: usize,
        range: Range<usize>,
        j: usize,
    ) -> Vec<GlueOp<T>> {
        let n_a = self.rats[i].len();
        let n_b = self.rats[j].len();

        if n_a == 0 || n_b == 0 || range.is_empty() || range.end > n_a {
            return vec![];
        }

        let mut results: Vec<GlueOp<T>> = self
            .valid_glues(i, j)
            .into_iter()
            .filter(|g| match_covers_range(g.start_a, g.match_len, &range, n_a))
            .collect();
        Self::sort_by_interval(&mut results);
        results
    }
}

fn match_covers_range(norm_start: i64, match_len: usize, range: &Range<usize>, n: usize) -> bool {
    if match_len < range.end - range.start {
        return false;
    }
    for v in range.clone() {
        let d = ((v as i64 - norm_start).rem_euclid(n as i64)) as usize;
        if d >= match_len {
            return false;
        }
    }
    true
}

/// Check whether a single-edge glue at positions `(ia, ib)` is a locally
/// consistent candidate based on interior angle sums at both junction vertices.
///
/// When two tiles share a single edge, two junction vertices are formed
/// at the endpoints. At each junction, four edges meet — two from each tile.
/// The tiles' interior angles at the junction sum to:
///
/// > (hturn − exterior_a) + (hturn − exterior_b) = 2·hturn − (exterior_a + exterior_b)
///
/// * **Overflow** (`exterior_a + exterior_b < 0`): the interior angles sum to
///   more than a full circle, so the neighbor edges physically overlap → the
///   glue creates a self-intersection.
/// * **Collinear** (`exterior_a + exterior_b == 0`): the interior angles sum
///   to exactly a full circle, so the neighbor edges are collinear → this is a
///   multi-edge extension already handled by the CMI index.
/// * **Gap** (`exterior_a + exterior_b > 0`): the interior angles sum to less
///   than a full circle, so there is a gap between the neighbor edges → the
///   candidate is locally consistent.
///
/// Returns `true` only when **both** junctions have a positive gap,
/// meaning the candidate is a true single-edge match that needs only
/// the global self-intersection check (`try_glue`).
fn is_single_edge_candidate(a: &[i8], ia: usize, b: &[i8], ib: usize) -> bool {
    let na = a.len();
    let nb = b.len();
    let left = a[ia] as i32 + b[ib] as i32;
    let right = a[(ia + 1) % na] as i32 + b[(ib + nb - 1) % nb] as i32;
    left > 0 && right > 0
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic::ZZ12;
    use crate::intgeom::rat::Rat;
    use crate::intgeom::snake::Snake;
    use crate::intgeom::tiles::{hexagon, spectre};
    use std::collections::HashSet;

    const MYSTIC_ZZ12: &[i8] = &[
        0, 2, -3, 2, 3, -2, 3, -2, 3, 2, -3, 2, 0, 2, -3, 2, 3, 2, -3, 2,
    ];

    #[test]
    fn spectre_self_finds_mystic() {
        let s: Snake<ZZ12> = spectre();
        let r: Rat<ZZ12> = Rat::from_unchecked(&s);
        let ts = TileSet::new(vec![r.clone()]);

        let glues = ts.valid_glues(0, 0);
        assert!(
            glues.iter().any(|g| g.result.seq() == MYSTIC_ZZ12),
            "mystic not found in {} glues: {:?}",
            glues.len(),
            glues,
        );
    }

    #[test]
    fn hexagon_self_single_edge() {
        let h: Snake<ZZ12> = hexagon();
        let r: Rat<ZZ12> = Rat::from_unchecked(&h);
        let ts = TileSet::new(vec![r.clone()]);

        let glues = ts.valid_glues(0, 0);
        assert!(!glues.is_empty(), "hexagon should have valid self-glues");

        let expected = &[-2, 2, 2, 2, 2, -2, 2, 2, 2, 2];
        assert!(
            glues.iter().any(|g| g.result.seq() == expected),
            "hex+hex glue not found in {} glues",
            glues.len(),
        );
    }

    #[test]
    fn self_intersecting_filtered() {
        let r1 = Rat::<ZZ12>::from_slice_unchecked(&[0, 0, 3, 0, 3, 3, -3, -3, 3, 3, 0, 3]);
        let r2 = Rat::<ZZ12>::from_slice_unchecked(&[0, 0, 3, 3, 0, 0, 3, 3]);
        let ts = TileSet::new(vec![r1.clone(), r2.clone()]);

        assert!(r1.try_glue((8, 0), &r2).is_err());

        let glues = ts.valid_glues(0, 1);
        for g in &glues {
            assert!(
                !(g.start_a == 8 && g.end_b == 0),
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
        let ts = TileSet::new(vec![r3.clone(), r4.clone()]);

        assert!(r3.try_glue((7, 0), &r4).is_err());

        let glues = ts.valid_glues(0, 1);
        for g in &glues {
            assert!(
                !(g.start_a == 7 && g.end_b == 0),
                "vertex-touching glue should be filtered"
            );
        }
    }

    #[test]
    #[should_panic(expected = "mixed chirality")]
    fn mixed_chirality_panics() {
        let ccw = Rat::<ZZ12>::from_slice_unchecked(&[1, 1, 1, 1]);
        let cw = ccw.reversed();
        let _ts = TileSet::new(vec![ccw, cw]);
    }

    #[test]
    fn deterministic_order() {
        let h: Snake<ZZ12> = hexagon();
        let r: Rat<ZZ12> = Rat::from_unchecked(&h);
        let ts = TileSet::new(vec![r]);
        let glues = ts.valid_glues(0, 0);
        for w in glues.windows(2) {
            assert!(
                w[0].start_a <= w[1].start_a
                    && (w[0].start_a < w[1].start_a || w[0].end_b <= w[1].end_b),
                "glues not in deterministic order"
            );
        }
    }

    #[test]
    fn multi_tile_all_valid_glues() {
        let sq = Rat::<ZZ12>::from_slice_unchecked(&[0, 1, 0, 1, 0, 1, 0, 1]);
        let tri = Rat::<ZZ12>::from_slice_unchecked(&[4, 4, 4]);
        let ts = TileSet::new(vec![sq, tri]);

        let all = ts.all_valid_glues();
        for g in &all {
            assert!(g.tile_a <= g.tile_b);
        }
        for w in all.windows(2) {
            assert!(
                w[0].tile_a < w[1].tile_a
                    || (w[0].tile_a == w[1].tile_a && w[0].tile_b < w[1].tile_b)
                    || (w[0].tile_a == w[1].tile_a
                        && w[0].tile_b == w[1].tile_b
                        && (w[0].start_a < w[1].start_a
                            || (w[0].start_a == w[1].start_a && w[0].end_b <= w[1].end_b))),
                "all_valid_glues not sorted"
            );
        }
    }

    #[test]
    fn containing_query_spectre() {
        let s: Snake<ZZ12> = spectre();
        let r: Rat<ZZ12> = Rat::from_unchecked(&s);
        let ts = TileSet::new(vec![r]);

        let glues = ts.valid_glues_containing(0, 2..3, 0);
        assert!(
            glues.iter().any(|g| g.result.seq() == MYSTIC_ZZ12),
            "containing query should find mystic, got {} glues: {:?}",
            glues.len(),
            glues,
        );
    }

    fn brute_force_valid_glues<T: IsComplex + IsRingOrField + Units>(
        a: &Rat<T>,
        b: &Rat<T>,
        i: usize,
        j: usize,
    ) -> Vec<GlueOp<T>> {
        let mut seen: HashSet<(i64, usize, i64)> = HashSet::new();
        let mut results: Vec<GlueOp<T>> = Vec::new();
        for ia in 0..a.len() {
            for ib in 0..b.len() {
                let (ns, len, ne) = a.get_match((ia as i64, ib as i64), b);
                if len == 0 || !seen.insert((ns, len, ne)) {
                    continue;
                }
                if let Ok(glued) = a.try_glue((ia as i64, ib as i64), b) {
                    results.push(GlueOp {
                        tile_a: i,
                        tile_b: j,
                        start_a: ns,
                        end_b: ne,
                        match_len: len,
                        result: glued,
                    });
                }
            }
        }
        results.sort_by(|x, y| {
            x.start_a
                .cmp(&y.start_a)
                .then_with(|| x.end_b.cmp(&y.end_b))
        });
        results
    }

    fn assert_glue_sets_match<T: IsComplex + IsRingOrField + Units>(
        ts_glues: &[GlueOp<T>],
        bf_glues: &[GlueOp<T>],
        label: &str,
    ) {
        assert_eq!(
            ts_glues.len(),
            bf_glues.len(),
            "{label}: count mismatch: TileSet={}, brute_force={}",
            ts_glues.len(),
            bf_glues.len(),
        );
        for (idx, (ts_g, bf_g)) in ts_glues.iter().zip(bf_glues.iter()).enumerate() {
            assert_eq!(
                (ts_g.start_a, ts_g.match_len, ts_g.end_b),
                (bf_g.start_a, bf_g.match_len, bf_g.end_b),
                "{label}: interval mismatch at index {idx}",
            );
            assert_eq!(
                ts_g.result, bf_g.result,
                "{label}: result mismatch at ({}, {})",
                ts_g.start_a, ts_g.end_b,
            );
        }
    }

    #[test]
    fn spectre_foldback_cases_accepted() {
        let s: Snake<ZZ12> = spectre();
        let r: Rat<ZZ12> = Rat::from_unchecked(&s);
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
        let h: Snake<ZZ12> = hexagon();
        let r: Rat<ZZ12> = Rat::from_unchecked(&h);
        let ts = TileSet::new(vec![r.clone(), r.clone()]);

        let ts_glues = ts.valid_glues(0, 1);
        let bf_glues = brute_force_valid_glues(&r, &r, 0, 1);
        assert_glue_sets_match(&ts_glues, &bf_glues, "hexagon pair");

        for g in &ts_glues {
            assert_eq!(g.match_len, 1, "hexagon pair glue should be single-edge");
            assert_eq!(g.result.len(), 10, "hex+hex result should be 10-gon");
        }

        let canonical = &ts_glues[0].result;
        for (idx, g) in ts_glues.iter().enumerate() {
            assert_eq!(
                g.result, *canonical,
                "hex+hex glue #{idx} should yield same shape",
            );
        }

        assert_eq!(
            ts_glues.len(),
            36,
            "hexagon pair should have 36 valid glues"
        );
    }

    #[test]
    fn hexamino_pair_exhaustive() {
        let h: Snake<ZZ12> = hexagon();
        let r: Rat<ZZ12> = Rat::from_unchecked(&h);
        let hexamino = r.glue((0, 0), &r);

        let ts = TileSet::new(vec![hexamino.clone(), hexamino.clone()]);
        let ts_glues = ts.valid_glues(0, 1);
        let bf_glues = brute_force_valid_glues(&hexamino, &hexamino, 0, 1);
        assert_glue_sets_match(&ts_glues, &bf_glues, "hexamino pair");

        assert!(
            !ts_glues.is_empty(),
            "hexamino pair should have valid glues"
        );

        let multi_count = ts_glues.iter().filter(|g| g.match_len > 1).count();
        let single_count = ts_glues.iter().filter(|g| g.match_len == 1).count();
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
        let s: Snake<ZZ12> = spectre();
        let r: Rat<ZZ12> = Rat::from_unchecked(&s);
        let ts = TileSet::new(vec![r.clone(), r.clone()]);

        let ts_glues = ts.valid_glues(0, 1);
        let bf_glues = brute_force_valid_glues(&r, &r, 0, 1);
        assert_glue_sets_match(&ts_glues, &bf_glues, "spectre pair");

        assert!(
            ts_glues.iter().any(|g| g.result.seq() == MYSTIC_ZZ12),
            "spectre pair should find mystic, got {} glues",
            ts_glues.len(),
        );

        let multi_count = ts_glues.iter().filter(|g| g.match_len > 1).count();
        let single_count = ts_glues.iter().filter(|g| g.match_len == 1).count();
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
        let ts = TileSet::new(vec![sq.clone(), tri.clone()]);

        for i in 0..2 {
            for j in i..2 {
                let ts_glues = ts.valid_glues(i, j);
                let bf_glues = brute_force_valid_glues(ts.rat(i), ts.rat(j), i, j);
                assert_glue_sets_match(&ts_glues, &bf_glues, &format!("mixed shapes ({i},{j})"));
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
    fn containing_query_exhaustive() {
        let spec = Rat::<ZZ12>::from_unchecked(&spectre());
        let hexamino =
            Rat::<ZZ12>::from_unchecked(&hexagon()).glue((0, 0), &Rat::from_unchecked(&hexagon()));

        for (label, rat) in [("spec", &spec), ("hexamino", &hexamino)] {
            let ts = TileSet::new(vec![rat.clone()]);
            let all = ts.valid_glues(0, 0);
            let n = rat.len();

            for edge in 0..n {
                let containing = ts.valid_glues_containing(0, edge..edge + 1, 0);
                let expected_intervals: Vec<_> = all
                    .iter()
                    .filter(|g| match_covers_range(g.start_a, g.match_len, &(edge..edge + 1), n))
                    .map(|g| (g.start_a, g.match_len, g.end_b))
                    .collect();
                let got_intervals: Vec<_> = containing
                    .iter()
                    .map(|g| (g.start_a, g.match_len, g.end_b))
                    .collect();
                assert_eq!(
                    got_intervals, expected_intervals,
                    "{label}: edge {edge}: containing query mismatch",
                );
            }
        }
    }
}
