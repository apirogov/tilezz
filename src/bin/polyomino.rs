use std::collections::BTreeMap;
use std::time::Instant;

use rustc_hash::FxHashSet;
use std::collections::BTreeSet;

use clap::Parser;
use tilezz::intgeom::angles::normalize_angle;
use tilezz::intgeom::rat::{lex_min_rot, Rat};
use tilezz::intgeom::tiles;

const HTURN: i8 = 2;
const TURN: i32 = 4;

#[derive(Parser)]
#[command(
    name = "polyomino",
    about = "Specialized free hole-free polyomino enumerator"
)]
struct Args {
    #[arg(long, default_value_t = 12)]
    max_size: usize,

    #[arg(long, default_value_t = false)]
    profile: bool,
}

struct Site {
    counter: usize,
    ns: usize,
    ml: usize,
    ne: usize,
}

struct Patch {
    angles: Vec<i8>,
    counters: Vec<usize>,
    next_counter: usize,
    positions: Vec<(i32, i32)>,
    visited: FxHashSet<(i32, i32)>,
    rat: Rat<tilezz::cyclotomic::ZZ4>,
}

fn dir_step(dir: i32) -> (i32, i32) {
    match dir.rem_euclid(TURN) {
        0 => (1, 0),
        1 => (0, 1),
        2 => (-1, 0),
        _ => (0, -1),
    }
}

fn trace_positions(start: (i32, i32), start_dir: i32, angles: &[i8]) -> Vec<(i32, i32)> {
    let mut pos = vec![start];
    let mut dir = start_dir;
    for &a in angles {
        dir = (dir + a as i32).rem_euclid(TURN);
        let (dx, dy) = dir_step(dir);
        let last = *pos.last().unwrap();
        pos.push((last.0 + dx, last.1 + dy));
    }
    pos
}

impl Patch {
    fn from_seed(seed: &Rat<tilezz::cyclotomic::ZZ4>) -> Self {
        let canonical = seed.clone().canonical();
        let seq = canonical.seq();
        let n = seq.len();
        let points = trace_positions((0, 0), 0, seq);
        let positions = points[..n].to_vec();
        let visited: FxHashSet<(i32, i32)> = points[..=n].iter().copied().collect();
        let rat = Rat::from_slice_unchecked(&seq);
        Patch {
            angles: seq.to_vec(),
            counters: (0..n).collect(),
            next_counter: n,
            positions,
            visited,
            rat,
        }
    }

    fn to_rat(&self) -> Rat<tilezz::cyclotomic::ZZ4> {
        self.rat.clone()
    }

    fn enumerate_sites(&self, seed_seq: &[i8], min_counter: usize) -> Vec<Site> {
        let n = self.angles.len();
        let m = seed_seq.len();
        let mut seen: BTreeSet<(usize, usize, usize)> = BTreeSet::new();
        let mut sites: Vec<Site> = Vec::new();

        for ia in 0..n {
            if self.counters[ia] < min_counter {
                continue;
            }
            for ib in 0..m {
                let (ns, ml, ne) = self.compute_match(ia, seed_seq, ib);
                if ml == 0 || !seen.insert((ns, ml, ne)) {
                    continue;
                }
                if ml >= 2 {
                    let gap_left =
                        self.angles[(ns + ml) % n] as i32 + seed_seq[(ne + m - ml) % m] as i32;
                    let gap_right = self.angles[ns] as i32 + seed_seq[ne] as i32;
                    if gap_left < 0 || gap_right < 0 {
                        continue;
                    }
                } else {
                    let left = self.angles[ia] as i32 + seed_seq[ib] as i32;
                    let right =
                        self.angles[(ia + 1) % n] as i32 + seed_seq[(ib + m - 1) % m] as i32;
                    if left <= 0 || right <= 0 {
                        continue;
                    }
                }
                sites.push(Site {
                    counter: self.counters[ns],
                    ns,
                    ml,
                    ne,
                });
            }
        }

        sites.sort_by_key(|s| s.counter);
        sites
    }

    fn compute_match(&self, ia: usize, seed_seq: &[i8], ib: usize) -> (usize, usize, usize) {
        let n = self.angles.len();
        let m = seed_seq.len();
        let x: Vec<i8> = (0..n).map(|i| self.angles[(ia + i) % n]).collect();
        let seed_slice: Vec<i8> = (0..m).map(|i| seed_seq[(ib + 1 + i) % m]).collect();
        let y: Vec<i8> = seed_slice.iter().rev().map(|a| -a).collect();

        let min_len = x.len().min(y.len());
        if min_len < 2 {
            return (ia, 1, ib);
        }

        let mut len = min_len;
        for i in 1..min_len {
            if x[i] != y[i] {
                len = i;
                break;
            }
        }
        if x[0] != y[0] {
            return (
                ((ia as i64).rem_euclid(n as i64)) as usize,
                len,
                ((ib as i64).rem_euclid(m as i64)) as usize,
            );
        }

        let remaining = min_len - len;
        for i in 1..remaining {
            if x[x.len() - i] != y[y.len() - i] {
                return (
                    ((ia as i64 - i as i64).rem_euclid(n as i64)) as usize,
                    len + i,
                    ((ib as i64 + i as i64).rem_euclid(m as i64)) as usize,
                );
            }
        }
        (
            ((ia as i64 - remaining as i64).rem_euclid(n as i64)) as usize,
            len + remaining,
            ((ib as i64 + remaining as i64).rem_euclid(m as i64)) as usize,
        )
    }

    fn try_apply(&self, site: &Site, seed_seq: &[i8]) -> Option<Patch> {
        let m = seed_seq.len();
        let n = self.angles.len();
        let ml = site.ml;
        let ns = site.ns;
        let ne = site.ne;

        let x_len = n - ml + 1;
        let y_len = m - ml + 1;

        let x: Vec<i8> = (0..x_len).map(|i| self.angles[(ns + ml + i) % n]).collect();
        let y: Vec<i8> = (0..y_len).map(|i| seed_seq[(ne + i) % m]).collect();

        let mut new_angles: Vec<i8> = Vec::with_capacity(x_len + y_len - 2);
        new_angles.extend_from_slice(&x[..x_len - 1]);
        new_angles.extend_from_slice(&y[..y_len - 1]);

        let a_yx = normalize_angle::<tilezz::cyclotomic::ZZ4>(x[0] + y[y_len - 1] - HTURN);
        let a_xy = normalize_angle::<tilezz::cyclotomic::ZZ4>(y[0] + x[x_len - 1] - HTURN);
        if a_yx.abs() == HTURN || a_xy.abs() == HTURN {
            return None;
        }
        new_angles[0] = a_yx;
        new_angles[x_len - 1] = a_xy;

        let junction_a = self.positions[(ns + ml) % n];
        let junction_b = self.positions[ns];

        let prev_idx = (ns + n - 1) % n;
        let in_dir_b = dir_between(self.positions[prev_idx], self.positions[ns]);

        let out_dir_b = (in_dir_b + a_xy as i32).rem_euclid(TURN);
        let mut trace_angles = vec![0i8];
        trace_angles.extend_from_slice(&y[1..y_len - 1]);
        let new_points = trace_positions(junction_b, out_dir_b, &trace_angles);

        let num_new = new_points.len();
        for &p in &new_points[1..num_new] {
            if p != junction_a && self.visited.contains(&p) {
                return None;
            }
        }

        let mut new_positions: Vec<(i32, i32)> = Vec::with_capacity(new_angles.len());
        for i in 0..(x_len - 1) {
            new_positions.push(self.positions[(ns + ml + i) % n]);
        }
        for &p in &new_points[..num_new - 1] {
            new_positions.push(p);
        }

        let mut new_visited = self.visited.clone();
        for &p in &new_points[..num_new] {
            new_visited.insert(p);
        }

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
        new_positions.rotate_left(offset);

        let rat = Rat::from_canonical_angles_unchecked(new_angles.clone());

        Some(Patch {
            angles: new_angles,
            counters: new_counters,
            next_counter: nc,
            positions: new_positions,
            visited: new_visited,
            rat,
        })
    }
}

fn dir_between(a: (i32, i32), b: (i32, i32)) -> i32 {
    let dx = b.0 - a.0;
    let dy = b.1 - a.1;
    match (dx, dy) {
        (1, 0) => 0,
        (0, 1) => 1,
        (-1, 0) => 2,
        (0, -1) => 3,
        _ => panic!("non-unit step from ({a:?}) to ({b:?})"),
    }
}

#[allow(dead_code)]
fn make_free(
    _onesided: &FxHashSet<Rat<tilezz::cyclotomic::ZZ4>>,
) -> FxHashSet<Rat<tilezz::cyclotomic::ZZ4>> {
    // Note: this function is now dead code since we store free polyominoes directly
    // Keeping it for API compatibility - would need the same transformation logic
    FxHashSet::default()
}

fn main() {
    let args = Args::parse();

    #[cfg(feature = "pprof")]
    let guard = pprof::ProfilerGuardBuilder::default()
        .frequency(1000)
        .blocklist(&["libc", "libgcc", "pthread", "vdso"])
        .build()
        .unwrap();

    let seed: Rat<tilezz::cyclotomic::ZZ4> = Rat::from_unchecked(&tiles::square());
    let seed_seq = seed.seq().to_vec();

    let t0 = Instant::now();
    let mut results: BTreeMap<usize, FxHashSet<Rat<tilezz::cyclotomic::ZZ4>>> = BTreeMap::new();
    let initial = Patch::from_seed(&seed);
    let initial_rat = initial.to_rat();
    let initial_free = std::cmp::min(initial_rat.clone(), initial_rat.reflected());
    results.entry(1).or_default().insert(initial_free);

    let mut profile_enumerate = 0u64;
    let mut profile_apply = 0u64;
    let mut enumerate_calls = 0usize;

    grow_recursive(
        &initial,
        &seed_seq,
        0,
        1,
        args.max_size,
        &mut results,
        &mut profile_enumerate,
        &mut profile_apply,
        &mut enumerate_calls,
    );
    let elapsed = t0.elapsed();

    // Results already in free form - just print them
    println!("=== Free hole-free polyominoes (incremental) ===");
    let mut total = 0;
    for size in 1..=args.max_size {
        let patches = results.get(&size);
        println!("  size {size:>3}: {:>6} patches", patches.map(|s| s.len()).unwrap_or(0));
        total += patches.map(|s| s.len()).unwrap_or(0);
    }
    println!("  total:   {total:>6} patches");
    println!("  time:    {elapsed:.2?}");

    if args.profile {
        eprintln!(
            "\nProfile: {} enumerate calls, {:.2}s enumerate, {:.2}s apply",
            enumerate_calls,
            profile_enumerate as f64 / 1e9,
            profile_apply as f64 / 1e9,
        );
    }

    #[cfg(feature = "pprof")]
    {
        if let Ok(report) = guard.report().build() {
            let flamegraph_file = "flamegraph.svg";
            let mut file = std::fs::File::create(flamegraph_file).unwrap();
            report.flamegraph(&mut file).unwrap();
            eprintln!("\nFlame graph written to {}", flamegraph_file);
        }
    }
}

fn grow_recursive(
    patch: &Patch,
    seed_seq: &[i8],
    min_counter: usize,
    current_size: usize,
    max_size: usize,
    results: &mut BTreeMap<usize, FxHashSet<Rat<tilezz::cyclotomic::ZZ4>>>,
    profile_enumerate: &mut u64,
    profile_apply: &mut u64,
    enumerate_calls: &mut usize,
) {
    if current_size >= max_size {
        return;
    }

    let t0 = Instant::now();
    *enumerate_calls += 1;
    let sites = patch.enumerate_sites(seed_seq, min_counter);
    *profile_enumerate += t0.elapsed().as_nanos() as u64;

    for site in &sites {
        let t0 = Instant::now();
        let new_patch = match patch.try_apply(site, seed_seq) {
            Some(p) => p,
            None => continue,
        };
        let rat = new_patch.to_rat();
        *profile_apply += t0.elapsed().as_nanos() as u64;

        // Store free (canonical) form directly to cut search space in half
        let free_form = std::cmp::min(rat.clone(), rat.reflected());

        let new_size = current_size + 1;
        if results.entry(new_size).or_default().insert(free_form) {
            grow_recursive(
                &new_patch,
                seed_seq,
                site.counter,
                new_size,
                max_size,
                results,
                profile_enumerate,
                profile_apply,
                enumerate_calls,
            );
        }
    }
}

#[allow(dead_code)]
fn enumerate_onesided(max_size: usize) -> BTreeMap<usize, FxHashSet<Rat<tilezz::cyclotomic::ZZ4>>> {
    let seed: Rat<tilezz::cyclotomic::ZZ4> = Rat::from_unchecked(&tiles::square());
    let seed_seq = seed.seq().to_vec();
    let initial = Patch::from_seed(&seed);
    let mut results: BTreeMap<usize, FxHashSet<Rat<tilezz::cyclotomic::ZZ4>>> = BTreeMap::new();
    results.entry(1).or_default().insert(initial.to_rat());
    let mut pe = 0u64;
    let mut pa = 0u64;
    let mut ec = 0usize;
    grow_recursive(
        &initial,
        &seed_seq,
        0,
        1,
        max_size,
        &mut results,
        &mut pe,
        &mut pa,
        &mut ec,
    );
    results
}

#[cfg(test)]
mod tests {
    use super::*;
    use tilezz::intgeom::rat::Rat;
    use tilezz::intgeom::tiles;

    fn square_seed() -> (Rat<tilezz::cyclotomic::ZZ4>, Vec<i8>) {
        let seed: Rat<tilezz::cyclotomic::ZZ4> = Rat::from_unchecked(&tiles::square());
        let seq = seed.seq().to_vec();
        (seed, seq)
    }

    #[test]
    fn test_dir_step() {
        assert_eq!(dir_step(0), (1, 0));
        assert_eq!(dir_step(1), (0, 1));
        assert_eq!(dir_step(2), (-1, 0));
        assert_eq!(dir_step(3), (0, -1));
        assert_eq!(dir_step(-1), (0, -1));
        assert_eq!(dir_step(4), (1, 0));
    }

    #[test]
    fn test_dir_between() {
        assert_eq!(dir_between((0, 0), (1, 0)), 0);
        assert_eq!(dir_between((0, 0), (0, 1)), 1);
        assert_eq!(dir_between((0, 0), (-1, 0)), 2);
        assert_eq!(dir_between((0, 0), (0, -1)), 3);
    }

    #[test]
    fn test_trace_positions_square() {
        let pts = trace_positions((0, 0), 0, &[1, 1, 1, 1]);
        assert_eq!(pts.len(), 5);
        assert_eq!(pts[0], (0, 0));
        assert_eq!(pts[1], (0, 1));
        assert_eq!(pts[2], (-1, 1));
        assert_eq!(pts[3], (-1, 0));
        assert_eq!(pts[4], (0, 0));
    }

    #[test]
    fn test_trace_positions_straight() {
        let pts = trace_positions((0, 0), 0, &[0, 0]);
        assert_eq!(pts, vec![(0, 0), (1, 0), (2, 0)]);
    }

    #[test]
    fn test_trace_positions_empty() {
        let pts = trace_positions((5, 5), 2, &[]);
        assert_eq!(pts, vec![(5, 5)]);
    }

    #[test]
    fn test_from_seed_square() {
        let (seed, _) = square_seed();
        let patch = Patch::from_seed(&seed);
        assert_eq!(patch.angles.len(), 4);
        assert_eq!(patch.positions.len(), 4);
        assert_eq!(patch.counters, vec![0, 1, 2, 3]);
        assert_eq!(patch.next_counter, 4);
        assert_eq!(patch.visited.len(), 4);
        assert!(patch.visited.contains(&(0, 0)));
        assert!(patch.visited.contains(&(0, 1)));
        assert!(patch.visited.contains(&(-1, 1)));
        assert!(patch.visited.contains(&(-1, 0)));
    }

    #[test]
    fn test_compute_match_matches_general_code() {
        let (seed, seed_seq) = square_seed();
        let patch = Patch::from_seed(&seed);
        let n = patch.angles.len();
        let m = seed_seq.len();

        for ia in 0..n {
            for ib in 0..m {
                let (ns, ml, ne) = patch.compute_match(ia, &seed_seq, ib);

                let (general_ns, general_ml, general_ne) =
                    seed.get_match((ia as i64, ib as i64), &seed);

                assert_eq!(
                    (ns, ml, ne),
                    (general_ns as usize, general_ml, general_ne as usize),
                    "ia={ia}, ib={ib}: polyomino ({ns},{ml},{ne}) != general ({general_ns},{general_ml},{general_ne})"
                );
            }
        }
    }

    #[test]
    fn test_compute_match_revcomp_correct() {
        let (seed, seed_seq) = square_seed();
        let patch = Patch::from_seed(&seed);
        let n = patch.angles.len();
        let m = seed_seq.len();

        for ia in 0..n {
            for ib in 0..m {
                let (ns, ml, ne) = patch.compute_match(ia, &seed_seq, ib);
                if ml < 2 {
                    continue;
                }

                for k in 0..ml {
                    let patch_angle = patch.angles[(ns + k) % n];
                    let seed_angle = seed_seq[(ne + m - k) % m];
                    assert_eq!(
                        patch_angle, -seed_angle,
                        "match at ia={ia}, ib={ib}, ne={ne}, ml={ml}, k={k}: patch angle {patch_angle} != neg seed angle {}",
                        -seed_angle
                    );
                }
            }
        }
    }

    #[test]
    fn test_enumerate_sites_square_finds_sites() {
        let (seed, seed_seq) = square_seed();
        let patch = Patch::from_seed(&seed);
        let sites = patch.enumerate_sites(&seed_seq, 0);
        assert!(!sites.is_empty(), "single square should have glue sites");
    }

    #[test]
    fn test_try_apply_produces_valid_rat() {
        let (seed, seed_seq) = square_seed();
        let patch = Patch::from_seed(&seed);
        let sites = patch.enumerate_sites(&seed_seq, 0);

        let mut valid_count = 0;
        for site in &sites {
            if let Some(new_patch) = patch.try_apply(site, &seed_seq) {
                let rat = new_patch.to_rat();
                assert_eq!(rat.len(), 6, "2 glued squares should produce a 6-gon");
                assert!(
                    rat.seq().iter().all(|&a| a.abs() <= HTURN),
                    "angles should be normalized: {:?}",
                    rat.seq()
                );
                assert_eq!(new_patch.positions.len(), 6);
                assert_eq!(new_patch.counters.len(), 6);
                assert_eq!(new_patch.visited.len(), 6);
                valid_count += 1;
            }
        }
        assert!(
            valid_count > 0,
            "at least one site should apply successfully"
        );
    }

    #[test]
    fn test_size2_matches_general_glue() {
        let seed: Rat<tilezz::cyclotomic::ZZ4> = Rat::from_unchecked(&tiles::square());
        let general_results: FxHashSet<Rat<tilezz::cyclotomic::ZZ4>> = (0..seed.len())
            .filter_map(|ia| seed.try_glue((ia as i64, 0), &seed).ok())
            .collect();

        let (_, seed_seq) = square_seed();
        let patch = Patch::from_seed(&seed);
        let sites = patch.enumerate_sites(&seed_seq, 0);
        let incremental_results: FxHashSet<Rat<tilezz::cyclotomic::ZZ4>> = sites
            .iter()
            .filter_map(|site| patch.try_apply(site, &seed_seq))
            .map(|p| p.to_rat())
            .collect();

        assert_eq!(
            incremental_results, general_results,
            "size-2 results differ: got {:?}, expected {:?}",
            incremental_results, general_results
        );
    }

    #[test]
    fn test_free_matches_oeis() {
        use tilezz::intgeom::rat::Rat;
        use tilezz::intgeom::tiles;
        let seed: Rat<tilezz::cyclotomic::ZZ4> = Rat::from_unchecked(&tiles::square());
        let max_size = 12;
        let results = enumerate_onesided(max_size);
        for k in 1..=max_size {
            let count = results.get(&k).map(|s| s.len()).unwrap_or(0);
            eprintln!("size {k}: {count}");
        }

        // Our implementation now stores free polyominoes directly
        // (not one-sided like before)
        // These should match OEIS A000105 (free polyominoes, one less than one-sided)
        let expected = [
            (1, 1),
            (2, 1),
            (3, 2),
            (4, 5),
            (5, 12),
            (6, 35),
            (7, 107),
            (8, 363),
            (9, 1248),
            (10, 4460),
            (11, 16094),
            (12, 58937),
        ];
        for (size, expected_count) in expected.iter() {
            let actual = results.get(size).map(|s| s.len()).unwrap_or(0);
            assert_eq!(
                actual,
                *expected_count,
                "size {}: expected {}, got {}",
                size,
                expected_count,
                actual
            );
        }
    }

    // Note: make_free is now unused since we store free form directly
    // Keeping it for API compatibility and potential future use

    #[test]
    fn test_visited_grows_correctly() {
        let (seed, seed_seq) = square_seed();
        let patch = Patch::from_seed(&seed);
        let sites = patch.enumerate_sites(&seed_seq, 0);

        for site in &sites {
            if let Some(new_patch) = patch.try_apply(site, &seed_seq) {
                for &p in &patch.visited {
                    assert!(
                        new_patch.visited.contains(&p),
                        "visited should be superset: missing {p:?}"
                    );
                }
                assert!(
                    new_patch.visited.len() > patch.visited.len(),
                    "visited should grow after a successful apply"
                );
                return;
            }
        }
        panic!("no valid apply found");
    }

    #[test]
    fn test_positions_form_cycle() {
        let (seed, seed_seq) = square_seed();
        let patch = Patch::from_seed(&seed);
        let sites = patch.enumerate_sites(&seed_seq, 0);

        for site in &sites {
            if let Some(new_patch) = patch.try_apply(site, &seed_seq) {
                let pts = trace_positions(
                    new_patch.positions[0],
                    if new_patch.positions.len() > 1 {
                        dir_between(new_patch.positions[0], new_patch.positions[1])
                    } else {
                        0
                    },
                    &new_patch.angles,
                );
                let n = new_patch.angles.len();
                assert_eq!(
                    pts.len(),
                    n + 1,
                    "tracing {} angles should give {} points",
                    n,
                    n + 1
                );
                assert_eq!(
                    *pts.first().unwrap(),
                    *pts.last().unwrap(),
                    "positions should form a closed cycle"
                );
                for i in 0..n {
                    assert_eq!(
                        new_patch.positions[i], pts[i],
                        "position {i} mismatch with traced points"
                    );
                }
            }
        }
    }
}
