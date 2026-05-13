use std::collections::BTreeMap;
use std::marker::PhantomData;
use std::sync::OnceLock;

use rustc_hash::FxHashSet;

use crate::cyclotomic::gaussint::GaussInt;
#[cfg(feature = "cyclotomic_intersect")]
use crate::cyclotomic::geometry::intersect;
use crate::cyclotomic::{IsComplex, IsRingOrField, SymNum, Units};
use crate::intgeom::angles::normalize_angle;
use crate::intgeom::rat::{lex_min_rot, Rat};
use crate::intgeom::snake::Snake;

pub(crate) struct GlueSite {
    pub(crate) counter: usize,
    pub(crate) norm_start: usize,
    pub(crate) match_len: usize,
    pub(crate) norm_end: usize,
}

pub trait PatchPos: Copy + Clone + Eq + std::hash::Hash + 'static {
    fn zero() -> Self;
    fn add_unit(&self, dir: i8) -> Self;
    fn sub_pos(&self, other: &Self) -> Self;
    fn direction_of(&self) -> i8;
    fn turn() -> i8;
}

fn scale_ratio(r: &num_rational::Ratio<i64>, scaling_fac: i64) -> i32 {
    (*r.numer() * (scaling_fac / *r.denom())) as i32
}

fn build_scaled_table<T, const N: usize>(scaling_fac: i64) -> Vec<[i32; N]>
where
    T: Units + SymNum<Scalar = GaussInt<num_rational::Ratio<i64>>>,
{
    (0..T::turn())
        .map(|d| {
            let unit = T::unit(d);
            let mut arr = [0i32; N];
            let mut idx = 0;
            for gi in unit.zz_coeffs() {
                arr[idx] = scale_ratio(&gi.real, scaling_fac);
                arr[idx + 1] = scale_ratio(&gi.imag, scaling_fac);
                idx += 2;
            }
            arr
        })
        .collect()
}

macro_rules! define_pos {
    ($name:ident, $n:expr) => {
        #[derive(Copy, Clone, PartialEq, Eq, std::hash::Hash, Debug, Default)]
        pub struct $name([i32; $n]);
    };
}

define_pos!(Pos2, 2);
define_pos!(Pos4, 4);
define_pos!(Pos8, 8);

macro_rules! impl_patch_pos {
    ($pos:ty, $n:expr, $zz:ty, $sf:expr) => {
        impl $pos {
            fn scaled_unit_table() -> &'static Vec<[i32; $n]> {
                static TABLE: OnceLock<Vec<[i32; $n]>> = OnceLock::new();
                TABLE.get_or_init(|| build_scaled_table::<$zz, $n>($sf))
            }
        }

        impl PatchPos for $pos {
            fn zero() -> Self {
                Self([0i32; $n])
            }

            fn add_unit(&self, dir: i8) -> Self {
                let tab = Self::scaled_unit_table();
                let mut r = self.0;
                let step = &tab[dir.rem_euclid(<$zz as SymNum>::turn()) as usize];
                for i in 0..$n {
                    r[i] += step[i];
                }
                Self(r)
            }

            fn sub_pos(&self, other: &Self) -> Self {
                let mut r = self.0;
                for i in 0..$n {
                    r[i] -= other.0[i];
                }
                Self(r)
            }

            fn direction_of(&self) -> i8 {
                let tab = Self::scaled_unit_table();
                for d in 0..tab.len() {
                    if self.0 == tab[d] {
                        return d as i8;
                    }
                }
                panic!("direction_of: not a unit step");
            }

            fn turn() -> i8 {
                <$zz as SymNum>::turn()
            }
        }
    };
}

impl_patch_pos!(Pos2, 2, crate::cyclotomic::ZZ4, 1);
impl_patch_pos!(Pos4, 4, crate::cyclotomic::ZZ12, 2);
impl_patch_pos!(Pos8, 8, crate::cyclotomic::ZZ10, 8);

#[cfg(not(feature = "cyclotomic_intersect"))]
#[derive(Copy, Clone, Debug)]
struct SR(i64, i64);

#[cfg(not(feature = "cyclotomic_intersect"))]
impl SR {
    fn mul(self, o: Self) -> Self {
        SR(self.0 * o.0 + 3 * self.1 * o.1, self.0 * o.1 + self.1 * o.0)
    }
    fn add(self, o: Self) -> Self {
        SR(self.0 + o.0, self.1 + o.1)
    }
    fn is_zero(self) -> bool {
        self.0 == 0 && self.1 == 0
    }
    fn is_positive(self) -> bool {
        sign_sqrt3(self.0, self.1) == std::cmp::Ordering::Greater
    }
}

#[cfg(not(feature = "cyclotomic_intersect"))]
fn sign_sqrt3(a: i64, b: i64) -> std::cmp::Ordering {
    if b == 0 {
        return a.cmp(&0);
    }
    if a == 0 {
        return b.cmp(&0);
    }
    let sa = a.cmp(&0);
    let sb = b.cmp(&0);
    if sa == sb {
        return sa;
    }
    let aa = a.unsigned_abs();
    let bb = b.unsigned_abs();
    match (aa * aa).cmp(&(3 * bb * bb)) {
        std::cmp::Ordering::Greater => sa,
        std::cmp::Ordering::Less => sb,
        std::cmp::Ordering::Equal => std::cmp::Ordering::Equal,
    }
}

#[cfg(not(feature = "cyclotomic_intersect"))]
fn re4(p: Pos4) -> SR {
    SR(i64::from(p.0[0]), i64::from(p.0[2]))
}

#[cfg(not(feature = "cyclotomic_intersect"))]
fn im4(p: Pos4) -> SR {
    SR(i64::from(p.0[1]), i64::from(p.0[3]))
}

#[cfg(not(feature = "cyclotomic_intersect"))]
fn wedge4(p1: Pos4, p2: Pos4) -> SR {
    re4(p1).mul(im4(p2)).sub(im4(p1).mul(re4(p2)))
}

#[cfg(not(feature = "cyclotomic_intersect"))]
fn dot4(p1: Pos4, p2: Pos4) -> SR {
    re4(p1).mul(re4(p2)).add(im4(p1).mul(im4(p2)))
}

#[cfg(not(feature = "cyclotomic_intersect"))]
impl SR {
    fn sub(self, o: Self) -> Self {
        SR(self.0 - o.0, self.1 - o.1)
    }
}

#[cfg(not(feature = "cyclotomic_intersect"))]
fn is_ccw4(p: Pos4, a: Pos4, b: Pos4) -> bool {
    wedge4(a.sub_pos(&p), b.sub_pos(&p)).is_positive()
}

#[cfg(not(feature = "cyclotomic_intersect"))]
fn is_between4(p: Pos4, a: Pos4, b: Pos4) -> bool {
    let v = a.sub_pos(&p);
    let w = p.sub_pos(&b);
    wedge4(v, w).is_zero() && dot4(v, w).is_positive()
}

#[cfg(not(feature = "cyclotomic_intersect"))]
fn intersect4(a: Pos4, b: Pos4, c: Pos4, d: Pos4) -> bool {
    if a == c || a == d || b == c || b == d {
        return false;
    }
    let l_ab_0 = im4(a.sub_pos(&b));
    let l_ab_1 = re4(b.sub_pos(&a));
    let l_ab_2 = wedge4(a, b);
    let colinear_c = l_ab_0
        .mul(re4(c))
        .add(l_ab_1.mul(im4(c)))
        .add(l_ab_2)
        .is_zero();
    let colinear_d = l_ab_0
        .mul(re4(d))
        .add(l_ab_1.mul(im4(d)))
        .add(l_ab_2)
        .is_zero();
    if colinear_c && colinear_d {
        is_between4(c, a, b) || is_between4(d, a, b)
    } else {
        is_ccw4(a, c, d) != is_ccw4(b, c, d) && is_ccw4(a, b, c) != is_ccw4(a, b, d)
    }
}

#[cfg(not(feature = "cyclotomic_intersect"))]
fn check_new_segments_cross_existing_pos4(
    new_points: &[Pos4],
    parent_positions: &[Pos4],
    n: usize,
    ns: usize,
    ml: usize,
) -> bool {
    let num_new = new_points.len();
    if num_new < 2 {
        return false;
    }
    let x_len = n - ml + 1;
    for ni in 0..(num_new - 1) {
        let na = new_points[ni];
        let nb = new_points[ni + 1];
        for xi in 0..(x_len - 1) {
            let ei = (ns + ml + xi) % n;
            let ea = parent_positions[ei];
            let eb = parent_positions[(ei + 1) % n];
            if intersect4(na, nb, ea, eb) {
                return true;
            }
        }
    }
    false
}

pub trait HasPatchPos: IsComplex + IsRingOrField + Units {
    type Pos: PatchPos;
    fn check_segment_crossings(
        new_points: &[Self::Pos],
        parent_positions: &[Self::Pos],
        n: usize,
        ns: usize,
        ml: usize,
    ) -> bool;
    fn needs_snake_validation() -> bool {
        true
    }
}

impl HasPatchPos for crate::cyclotomic::ZZ4 {
    type Pos = Pos2;
    fn check_segment_crossings(
        _new_points: &[Self::Pos],
        _parent_positions: &[Self::Pos],
        _n: usize,
        _ns: usize,
        _ml: usize,
    ) -> bool {
        false
    }
    fn needs_snake_validation() -> bool {
        false
    }
}

impl HasPatchPos for crate::cyclotomic::ZZ10 {
    type Pos = Pos8;
    fn check_segment_crossings(
        _new_points: &[Self::Pos],
        _parent_positions: &[Self::Pos],
        _n: usize,
        _ns: usize,
        _ml: usize,
    ) -> bool {
        false
    }
}

#[cfg(feature = "cyclotomic_intersect")]
impl PatchPos for crate::cyclotomic::ZZ12 {
    fn zero() -> Self {
        <Self as num_traits::Zero>::zero()
    }
    fn add_unit(&self, dir: i8) -> Self {
        *self + Self::unit(dir)
    }
    fn sub_pos(&self, other: &Self) -> Self {
        *self - *other
    }
    fn direction_of(&self) -> i8 {
        for d in 0..<Self as SymNum>::turn() {
            if *self == Self::unit(d) {
                return d;
            }
        }
        panic!("direction_of: not a unit step");
    }
    fn turn() -> i8 {
        <Self as SymNum>::turn()
    }
}

#[cfg(feature = "cyclotomic_intersect")]
fn check_segments_cross_cyclotomic(
    new_points: &[crate::cyclotomic::ZZ12],
    parent_positions: &[crate::cyclotomic::ZZ12],
    n: usize,
    ns: usize,
    ml: usize,
) -> bool {
    let num_new = new_points.len();
    if num_new < 2 {
        return false;
    }
    let x_len = n - ml + 1;
    for ni in 0..(num_new - 1) {
        let na = new_points[ni];
        let nb = new_points[ni + 1];
        for xi in 0..(x_len - 1) {
            let ei = (ns + ml + xi) % n;
            let ea = parent_positions[ei];
            let eb = parent_positions[(ei + 1) % n];
            if intersect(&(na, nb), &(ea, eb)) {
                return true;
            }
        }
    }
    false
}

#[cfg(feature = "cyclotomic_intersect")]
impl HasPatchPos for crate::cyclotomic::ZZ12 {
    type Pos = Self;
    fn check_segment_crossings(
        new_points: &[Self::Pos],
        parent_positions: &[Self::Pos],
        n: usize,
        ns: usize,
        ml: usize,
    ) -> bool {
        check_segments_cross_cyclotomic(new_points, parent_positions, n, ns, ml)
    }
    fn needs_snake_validation() -> bool {
        false
    }
}

#[cfg(not(feature = "cyclotomic_intersect"))]
impl HasPatchPos for crate::cyclotomic::ZZ12 {
    type Pos = Pos4;
    fn check_segment_crossings(
        new_points: &[Self::Pos],
        parent_positions: &[Self::Pos],
        n: usize,
        ns: usize,
        ml: usize,
    ) -> bool {
        check_new_segments_cross_existing_pos4(new_points, parent_positions, n, ns, ml)
    }
    fn needs_snake_validation() -> bool {
        false
    }
}

fn trace_positions_from<P: PatchPos>(start: P, start_dir: i8, angles: &[i8]) -> Vec<P> {
    let mut positions = vec![start];
    let mut dir = start_dir;
    for &a in angles {
        dir = (dir as i64 + a as i64).rem_euclid(P::turn() as i64) as i8;
        let last = *positions.last().unwrap();
        positions.push(last.add_unit(dir));
    }
    positions
}

pub struct RedelmeierPatch<T: HasPatchPos> {
    angles: Vec<i8>,
    counters: Vec<usize>,
    next_counter: usize,
    cached_rat: Option<Rat<T>>,
    positions: Vec<T::Pos>,
    visited: FxHashSet<T::Pos>,
    _phantom: PhantomData<T>,
}

impl<T: HasPatchPos> RedelmeierPatch<T> {
    pub fn len(&self) -> usize {
        self.angles.len()
    }

    pub fn is_empty(&self) -> bool {
        self.angles.len() == 0
    }
}

impl<T: HasPatchPos> RedelmeierPatch<T> {
    pub fn from_seed(seed: &Rat<T>) -> Self {
        let canonical = seed.clone().canonical();
        let seq = canonical.seq();
        let n = seq.len();
        let rat = Rat::from_slice_unchecked(seq);
        let all_positions = trace_positions_from(T::Pos::zero(), 0, seq);
        let positions = all_positions[..n].to_vec();
        let visited = all_positions.into_iter().collect();
        RedelmeierPatch {
            angles: seq.to_vec(),
            counters: (0..n).collect(),
            next_counter: n,
            cached_rat: Some(rat),
            positions,
            visited,
            _phantom: PhantomData,
        }
    }

    pub fn to_rat(&self) -> Rat<T> {
        self.cached_rat
            .clone()
            .unwrap_or_else(|| Rat::from_slice_unchecked(&self.angles))
    }

    fn compute_match_inline(&self, ia: usize, seed_seq: &[i8], ib: usize) -> (usize, usize, usize) {
        let n = self.angles.len();
        let m = seed_seq.len();
        let min_len = n.min(m);

        if min_len < 2 {
            return (ia, 1, ib);
        }

        let mut len = min_len;
        for i in 1..min_len {
            if self.angles[(ia + i) % n] != -seed_seq[(ib + m - i) % m] {
                len = i;
                break;
            }
        }

        if self.angles[ia] != -seed_seq[ib] {
            return (ia, len, ib);
        }

        let remaining = min_len - len;
        for i in 1..remaining {
            if self.angles[(ia + n - i) % n] != -seed_seq[(ib + i) % m] {
                return (
                    (ia as i64 - i as i64).rem_euclid(n as i64) as usize,
                    len + i,
                    (ib as i64 + i as i64).rem_euclid(m as i64) as usize,
                );
            }
        }
        (
            (ia as i64 - remaining as i64).rem_euclid(n as i64) as usize,
            len + remaining,
            (ib as i64 + remaining as i64).rem_euclid(m as i64) as usize,
        )
    }

    fn enumerate_sites_inline(&self, seed_seq: &[i8], min_counter: usize) -> Vec<GlueSite> {
        let n = self.angles.len();
        let m = seed_seq.len();
        let mut seen: FxHashSet<(usize, usize, usize)> = FxHashSet::default();
        let mut sites: Vec<GlueSite> = Vec::new();

        for ia in 0..n {
            if self.counters[ia] < min_counter {
                continue;
            }
            for ib in 0..m {
                let (ns, ml, ne) = self.compute_match_inline(ia, seed_seq, ib);
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
                sites.push(GlueSite {
                    counter: self.counters[ns],
                    norm_start: ns,
                    match_len: ml,
                    norm_end: ne,
                });
            }
        }

        sites.sort_by_key(|s| s.counter);
        sites
    }

    fn apply_site(&self, site: &GlueSite, seed: &Rat<T>) -> Option<RedelmeierPatch<T>> {
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
        if a_yx.abs() == T::hturn() || a_xy.abs() == T::hturn() {
            return None;
        }
        new_angles[0] = a_yx;
        new_angles[x_len - 1] = a_xy;

        let junction_a = self.positions[(ns + ml) % n];
        let junction_b = self.positions[ns];

        let prev_idx = (ns + n - 1) % n;
        let step_in = self.positions[ns].sub_pos(&self.positions[prev_idx]);
        let in_dir = step_in.direction_of();
        let out_dir: i8 = (in_dir as i64 + a_xy as i64).rem_euclid(T::turn() as i64) as i8;

        let mut trace_angles: Vec<i8> = Vec::with_capacity(y_len - 1);
        trace_angles.push(0i8);
        if y_len > 2 {
            trace_angles.extend_from_slice(&y[1..y_len - 1]);
        }
        let new_points = trace_positions_from(junction_b, out_dir, &trace_angles);

        let num_new = new_points.len();
        for p in &new_points[1..num_new] {
            if *p != junction_a && self.visited.contains(p) {
                return None;
            }
        }

        if T::check_segment_crossings(&new_points, &self.positions, n, ns, ml) {
            return None;
        }

        let mut new_positions: Vec<T::Pos> = Vec::with_capacity(new_angles.len());
        for i in 0..(x_len - 1) {
            new_positions.push(self.positions[(ns + ml + i) % n]);
        }
        for p in &new_points[..num_new - 1] {
            new_positions.push(*p);
        }

        let mut new_visited = self.visited.clone();
        for p in &new_points {
            new_visited.insert(*p);
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

        let cached_rat = Rat::from_canonical_angles_unchecked(new_angles.clone());

        Some(RedelmeierPatch {
            angles: new_angles,
            counters: new_counters,
            next_counter: nc,
            cached_rat: Some(cached_rat),
            positions: new_positions,
            visited: new_visited,
            _phantom: PhantomData,
        })
    }
}

#[derive(Default)]
pub struct GrowStats {
    pub enumerate_calls: usize,
    pub enumerate_ns: u64,
    pub apply_ns: u64,
}

impl std::fmt::Display for GrowStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "enumerate: {} calls, {:.2}s | apply: {:.2}s",
            self.enumerate_calls,
            self.enumerate_ns as f64 / 1e9,
            self.apply_ns as f64 / 1e9,
        )
    }
}

struct GrowState<T: HasPatchPos> {
    results: BTreeMap<usize, FxHashSet<Rat<T>>>,
    stats: GrowStats,
}

pub fn grow_redelmeier<T>(seed: &Rat<T>, max_size: usize) -> BTreeMap<usize, FxHashSet<Rat<T>>>
where
    T: HasPatchPos,
{
    grow_redelmeier_inner(seed, max_size, false).0
}

pub fn grow_redelmeier_free<T>(seed: &Rat<T>, max_size: usize) -> BTreeMap<usize, FxHashSet<Rat<T>>>
where
    T: HasPatchPos,
{
    grow_redelmeier_inner(seed, max_size, true).0
}

pub fn grow_redelmeier_profiled<T>(
    seed: &Rat<T>,
    max_size: usize,
) -> (BTreeMap<usize, FxHashSet<Rat<T>>>, GrowStats)
where
    T: HasPatchPos,
{
    grow_redelmeier_inner(seed, max_size, false)
}

fn grow_redelmeier_inner<T>(
    seed: &Rat<T>,
    max_size: usize,
    free: bool,
) -> (BTreeMap<usize, FxHashSet<Rat<T>>>, GrowStats)
where
    T: HasPatchPos,
{
    let mut state = GrowState {
        results: BTreeMap::new(),
        stats: GrowStats::default(),
    };
    let initial = RedelmeierPatch::from_seed(seed);
    let initial_rat = initial.to_rat();
    let stored = if free {
        std::cmp::min(initial_rat.clone(), initial_rat.reflected())
    } else {
        initial_rat
    };
    state.results.entry(1).or_default().insert(stored);
    grow_recursive_inner(&initial, seed, 0, 1, max_size, free, &mut state);
    (state.results, state.stats)
}

fn grow_recursive_inner<T>(
    patch: &RedelmeierPatch<T>,
    seed: &Rat<T>,
    min_counter: usize,
    current_size: usize,
    max_size: usize,
    free: bool,
    state: &mut GrowState<T>,
) where
    T: HasPatchPos,
{
    if current_size >= max_size {
        return;
    }

    let seed_seq = seed.seq();

    let sites = {
        state.stats.enumerate_calls += 1;
        let t0 = std::time::Instant::now();
        let sites = patch.enumerate_sites_inline(seed_seq, min_counter);
        state.stats.enumerate_ns += t0.elapsed().as_nanos() as u64;
        sites
    };

    for site in &sites {
        let t0 = std::time::Instant::now();
        let new_patch = match patch.apply_site(site, seed) {
            Some(p) => p,
            None => {
                state.stats.apply_ns += t0.elapsed().as_nanos() as u64;
                continue;
            }
        };
        let rat = new_patch.to_rat();

        if T::needs_snake_validation() && Snake::<T>::try_from(rat.seq()).is_err() {
            state.stats.apply_ns += t0.elapsed().as_nanos() as u64;
            continue;
        }

        state.stats.apply_ns += t0.elapsed().as_nanos() as u64;

        let new_size = current_size + 1;
        let stored = if free {
            std::cmp::min(rat.clone(), rat.reflected())
        } else {
            rat
        };
        if state.results.entry(new_size).or_default().insert(stored) {
            grow_recursive_inner(
                &new_patch,
                seed,
                site.counter,
                new_size,
                max_size,
                free,
                state,
            );
        }
    }
}

pub fn make_free<T>(onesided: &FxHashSet<Rat<T>>) -> FxHashSet<Rat<T>>
where
    T: HasPatchPos,
{
    onesided
        .iter()
        .map(|r| std::cmp::min(r.clone(), r.reflected()))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[cfg(not(feature = "cyclotomic_intersect"))]
    use crate::cyclotomic::geometry::intersect;
    #[cfg(not(feature = "cyclotomic_intersect"))]
    use crate::cyclotomic::OneImag;
    use crate::cyclotomic::{ZZ12, ZZ4};
    use crate::intgeom::matchtypes::MatchFinder;
    use crate::intgeom::tiles;
    use crate::intgeom::tileset::TileSet;
    #[cfg(not(feature = "cyclotomic_intersect"))]
    use num_traits::{One, Zero};
    use std::sync::Arc;

    #[cfg(not(feature = "cyclotomic_intersect"))]
    fn zz12_to_pos4(t: ZZ12) -> Pos4 {
        let mut arr = [0i32; 4];
        let mut idx = 0;
        for gi in t.zz_coeffs() {
            arr[idx] = scale_ratio(&gi.real, 2);
            arr[idx + 1] = scale_ratio(&gi.imag, 2);
            idx += 2;
        }
        Pos4(arr)
    }

    #[cfg(not(feature = "cyclotomic_intersect"))]
    fn build_patch_positions(angles: &[i8]) -> Vec<ZZ12> {
        let mut positions = vec![ZZ12::zero()];
        let mut dir: i8 = 0;
        for &a in angles {
            dir = (dir as i64 + a as i64).rem_euclid(12) as i8;
            let last = *positions.last().unwrap();
            positions.push(last + ZZ12::unit(dir));
        }
        positions
    }

    #[cfg(not(feature = "cyclotomic_intersect"))]
    #[test]
    fn segment_intersect_matches_cyclotomic() {
        let cases: Vec<(ZZ12, ZZ12, ZZ12, ZZ12)> = vec![
            (ZZ12::zero(), ZZ12::one(), ZZ12::from(2), ZZ12::from(3)),
            (ZZ12::zero(), ZZ12::one(), ZZ12::one(), ZZ12::from(2)),
            (ZZ12::zero(), ZZ12::from(2), ZZ12::one(), ZZ12::from(3)),
            (
                ZZ12::zero(),
                ZZ12::one_i(),
                ZZ12::one(),
                ZZ12::one() + ZZ12::one_i(),
            ),
            (
                ZZ12::zero(),
                ZZ12::one() + ZZ12::one_i(),
                ZZ12::one(),
                ZZ12::one_i(),
            ),
        ];
        for (a, b, c, d) in &cases {
            let expected = intersect(&(*a, *b), &(*c, *d));
            let got = intersect4(
                zz12_to_pos4(*a),
                zz12_to_pos4(*b),
                zz12_to_pos4(*c),
                zz12_to_pos4(*d),
            );
            assert_eq!(
                got, expected,
                "intersect({a:?}, {b:?}, {c:?}, {d:?}): expected {expected}, got {got}"
            );
        }
    }

    #[cfg(not(feature = "cyclotomic_intersect"))]
    #[test]
    fn segment_intersect_fuzz_against_cyclotomic() {
        let hex_angles: &[i8] = &[2, 2, 2, 2, 2, 2];
        let spectre_angles: &[i8] = &[3, 2, 0, 2, -3, 2, 3, 2, -3, 2, 3, -2, 3, -2];
        for angles in [hex_angles, spectre_angles] {
            let positions = build_patch_positions(angles);
            let n = positions.len() - 1;
            for i in 0..n {
                for k in 0..n {
                    let a = positions[i];
                    let b = positions[(i + 1) % (n + 1)];
                    let c = positions[k];
                    let d = positions[(k + 1) % (n + 1)];
                    let expected = intersect(&(a, b), &(c, d));
                    let got = intersect4(
                        zz12_to_pos4(a),
                        zz12_to_pos4(b),
                        zz12_to_pos4(c),
                        zz12_to_pos4(d),
                    );
                    assert_eq!(
                        got, expected,
                        "angles={angles:?} seg({i},{}) vs seg({k},{}): expected {expected}, got {got}",
                        (i+1) % (n+1), (k+1) % (n+1)
                    );
                }
            }
        }
    }

    #[test]
    fn hexagon_matches_old_approach_size4() {
        let seed: Rat<ZZ12> = Rat::from_unchecked(&tiles::hexagon());
        let new_results = grow_redelmeier(&seed, 4);

        let mut old_results: BTreeMap<usize, FxHashSet<Rat<ZZ12>>> = BTreeMap::new();
        old_results.insert(1, std::iter::once(seed.clone()).collect());

        for k in 2..=4 {
            let prev: Vec<Rat<ZZ12>> = old_results[&(k - 1)].iter().cloned().collect();
            let patch_ts = Arc::new(TileSet::new(prev));
            let seed_ts = Arc::new(TileSet::new(vec![seed.clone()]));
            let mf = MatchFinder::crossing(patch_ts, seed_ts);
            let pairs: Vec<(usize, usize)> = (0..mf.num_tiles_a())
                .flat_map(|i| (0..mf.num_tiles_b()).map(move |j| (i, j)))
                .collect();
            let results = mf.valid_results_for_pairs(&pairs);
            old_results.insert(k, results.into_iter().collect());
        }

        for k in 1..=4 {
            let old_count = old_results.get(&k).map(|s| s.len()).unwrap_or(0);
            let new_count = new_results.get(&k).map(|s| s.len()).unwrap_or(0);
            assert_eq!(
                old_count, new_count,
                "size {k}: old={old_count} new={new_count}"
            );
            if let (Some(old_set), Some(new_set)) = (old_results.get(&k), new_results.get(&k)) {
                assert_eq!(
                    old_set.len(),
                    new_set.len(),
                    "size {k}: sets differ in length"
                );
                for r in old_set.iter() {
                    assert!(new_set.contains(r), "size {k}: missing {r}");
                }
            }
        }
    }

    #[test]
    fn spectre_matches_old_approach_size3() {
        let seed: Rat<ZZ12> = Rat::from_unchecked(&tiles::spectre());
        let new_results = grow_redelmeier(&seed, 3);

        let mut old_results: BTreeMap<usize, FxHashSet<Rat<ZZ12>>> = BTreeMap::new();
        old_results.insert(1, std::iter::once(seed.clone()).collect());

        for k in 2..=3 {
            let prev: Vec<Rat<ZZ12>> = old_results[&(k - 1)].iter().cloned().collect();
            let patch_ts = Arc::new(TileSet::new(prev));
            let seed_ts = Arc::new(TileSet::new(vec![seed.clone()]));
            let mf = MatchFinder::crossing(patch_ts, seed_ts);
            let pairs: Vec<(usize, usize)> = (0..mf.num_tiles_a())
                .flat_map(|i| (0..mf.num_tiles_b()).map(move |j| (i, j)))
                .collect();
            let results = mf.valid_results_for_pairs(&pairs);
            old_results.insert(k, results.into_iter().collect());
        }

        for k in 1..=3 {
            let old_count = old_results.get(&k).map(|s| s.len()).unwrap_or(0);
            let new_count = new_results.get(&k).map(|s| s.len()).unwrap_or(0);
            assert_eq!(
                old_count, new_count,
                "size {k}: old={old_count} new={new_count}"
            );
            if let (Some(old_set), Some(new_set)) = (old_results.get(&k), new_results.get(&k)) {
                assert_eq!(
                    old_set.len(),
                    new_set.len(),
                    "size {k}: sets differ in length"
                );
                for r in old_set.iter() {
                    assert!(new_set.contains(r), "size {k}: missing {r}");
                }
            }
        }
    }

    #[test]
    fn square_zz4_matches_old_approach_size6() {
        let seed: Rat<ZZ4> = Rat::from_unchecked(&tiles::square());
        let new_results = grow_redelmeier(&seed, 6);

        let mut old_results: BTreeMap<usize, FxHashSet<Rat<ZZ4>>> = BTreeMap::new();
        old_results.insert(1, std::iter::once(seed.clone()).collect());

        for k in 2..=6 {
            let prev: Vec<Rat<ZZ4>> = old_results[&(k - 1)].iter().cloned().collect();
            let patch_ts = Arc::new(TileSet::new(prev));
            let seed_ts = Arc::new(TileSet::new(vec![seed.clone()]));
            let mf = MatchFinder::crossing(patch_ts, seed_ts);
            let pairs: Vec<(usize, usize)> = (0..mf.num_tiles_a())
                .flat_map(|i| (0..mf.num_tiles_b()).map(move |j| (i, j)))
                .collect();
            let results = mf.valid_results_for_pairs(&pairs);
            old_results.insert(k, results.into_iter().collect());
        }

        for k in 1..=6 {
            let old_count = old_results.get(&k).map(|s| s.len()).unwrap_or(0);
            let new_count = new_results.get(&k).map(|s| s.len()).unwrap_or(0);
            assert_eq!(
                old_count, new_count,
                "size {k}: old={old_count} new={new_count}"
            );
            if let (Some(old_set), Some(new_set)) = (old_results.get(&k), new_results.get(&k)) {
                assert_eq!(
                    old_set.len(),
                    new_set.len(),
                    "size {k}: sets differ in length"
                );
                for r in old_set.iter() {
                    assert!(new_set.contains(r), "size {k}: missing {r}");
                }
            }
        }
    }

    fn brute_force_grow<T: HasPatchPos>(
        seed: &Rat<T>,
        max_size: usize,
    ) -> BTreeMap<usize, FxHashSet<Rat<T>>> {
        let mut results: BTreeMap<usize, FxHashSet<Rat<T>>> = BTreeMap::new();
        results.insert(1, std::iter::once(seed.clone()).collect());

        for k in 2..=max_size {
            let prev: Vec<Rat<T>> = results[&(k - 1)].iter().cloned().collect();
            let mut next: FxHashSet<Rat<T>> = FxHashSet::default();
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
        eprintln!("=== Redelmeier (MatchFinder) ===");
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
                for r in bf_set.iter() {
                    assert!(
                        redel_set.contains(r),
                        "size {k}: brute_force patch missing from redelmeier: {:?}",
                        r.seq()
                    );
                }
                for r in redel_set.iter() {
                    assert!(
                        bf_set.contains(r),
                        "size {k}: redelmeier patch missing from brute_force: {:?}",
                        r.seq()
                    );
                }
            }
        }
    }

    const FREE_POLYOMINOES_NO_HOLES: &[usize] = &[
        1, 1, 1, 2, 5, 12, 35, 107, 363, 1248, 4460, 16094, 58937, 217117,
    ];

    #[test]
    fn square_zz4_free_polyominoes_match_oeis() {
        let max_size: usize = 9;
        let seed: Rat<ZZ4> = Rat::from_unchecked(&tiles::square());
        let onesided = grow_redelmeier(&seed, max_size);

        for (k, &expected) in FREE_POLYOMINOES_NO_HOLES
            .iter()
            .enumerate()
            .take(max_size + 1)
            .skip(1)
        {
            let onesided_set = onesided.get(&k).unwrap();
            let free_set = make_free(onesided_set);
            let free_count = free_set.len();
            assert_eq!(
                free_count, expected,
                "size {k}: got {free_count}, expected {expected} free polyominoes (no holes)"
            );
        }
    }

    #[test]
    fn grow_redelmeier_free_matches_oeis() {
        let max_size: usize = 9;
        let seed: Rat<ZZ4> = Rat::from_unchecked(&tiles::square());
        let free_results = grow_redelmeier_free(&seed, max_size);

        for (k, &expected) in FREE_POLYOMINOES_NO_HOLES
            .iter()
            .enumerate()
            .take(max_size + 1)
            .skip(1)
        {
            let free_count = free_results.get(&k).map(|s| s.len()).unwrap_or(0);
            assert_eq!(
                free_count, expected,
                "size {k}: got {free_count}, expected {expected} free polyominoes (no holes)"
            );
        }
    }
}
