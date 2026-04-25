use std::collections::{BTreeMap, VecDeque};
use std::sync::Arc;

use crate::cyclotomic::{IsComplex, IsRingOrField, Units};
use crate::intgeom::matchtypes::{MatchFinder, MatchType, MatchTypeIndex};
use crate::intgeom::rat::Rat;
use crate::intgeom::tileset::TileSet;

type BoundaryMap = Vec<(usize, usize)>;

struct VertexTypeEntry<T: IsComplex> {
    matches: Vec<i32>,
    rat: Rat<T>,
    vertex_pos: usize,
    boundary_map: BoundaryMap,
    is_closed: bool,
    is_terminal: bool,
}

pub struct VertexTypeIndex<T: IsComplex> {
    seed_tileset: Arc<TileSet<T>>,
    match_type_index: Arc<MatchTypeIndex<T>>,
    entries: Vec<VertexTypeEntry<T>>,
}

impl<T: IsComplex + IsRingOrField + Units> VertexTypeIndex<T> {
    pub fn new(seed: Arc<TileSet<T>>) -> Self {
        let mti = Arc::new(MatchTypeIndex::new(Arc::clone(&seed)));

        let mut entries: Vec<VertexTypeEntry<T>> = Vec::new();
        let mut seq_index: BTreeMap<Vec<i32>, usize> = BTreeMap::new();
        let mut queue: VecDeque<usize> = VecDeque::new();

        for id in 1..=mti.num_types() {
            let rat = mti.apply(id as i32);
            let first = mti.get(id);
            let n_a = seed.rat(first.tile_a).len();
            let n_b = seed.rat(first.tile_b).len();
            let mlen = first.len;
            let result_len = n_a + n_b - 2 * mlen;

            let mut bm = Vec::with_capacity(result_len);
            for i in 0..n_a - mlen {
                bm.push((first.tile_a, (first.start_a + mlen + i) % n_a));
            }
            for j in 0..n_b - mlen {
                bm.push((first.tile_b, (first.start_b + j) % n_b));
            }

            let vpos_plus = n_a - mlen;
            let matches_plus = vec![id as i32];
            let idx_plus = entries.len();
            entries.push(VertexTypeEntry {
                matches: matches_plus.clone(),
                rat: rat.clone(),
                vertex_pos: vpos_plus,
                boundary_map: bm.clone(),
                is_closed: false,
                is_terminal: false,
            });
            seq_index.insert(matches_plus, idx_plus);
            queue.push_back(idx_plus);

            let vpos_minus = 0usize;
            let matches_minus = vec![-(id as i32)];
            let idx_minus = entries.len();
            entries.push(VertexTypeEntry {
                matches: matches_minus.clone(),
                rat: rat.clone(),
                vertex_pos: vpos_minus,
                boundary_map: bm,
                is_closed: false,
                is_terminal: false,
            });
            seq_index.insert(matches_minus, idx_minus);
            queue.push_back(idx_minus);
        }

        let mut vt_index = VertexTypeIndex {
            seed_tileset: seed,
            match_type_index: mti,
            entries,
        };

        while let Some(vt_idx) = queue.pop_front() {
            let is_closed = vt_index.entries[vt_idx].is_closed;
            let is_terminal = vt_index.entries[vt_idx].is_terminal;
            if is_closed || is_terminal {
                continue;
            }

            let extensions = vt_index.compute_extensions(vt_idx);

            if extensions.is_empty() {
                vt_index.entries[vt_idx].is_terminal = true;
                continue;
            }

            for ext in extensions {
                if seq_index.contains_key(&ext.new_matches) {
                    continue;
                }

                let new_idx = vt_index.entries.len();
                seq_index.insert(ext.new_matches.clone(), new_idx);
                vt_index.entries.push(VertexTypeEntry {
                    matches: ext.new_matches,
                    rat: ext.rat,
                    vertex_pos: ext.vertex_pos,
                    boundary_map: ext.boundary_map,
                    is_closed: ext.is_closed,
                    is_terminal: false,
                });

                if !ext.is_closed {
                    queue.push_back(new_idx);
                }

                if vt_index.entries.len() > 50000 {
                    eprintln!(
                        "WARNING: vertex type limit reached at {} types, stopping",
                        vt_index.entries.len()
                    );
                    return vt_index;
                }
            }
        }

        vt_index
    }

    fn compute_extensions(&self, vt_idx: usize) -> Vec<Extension<T>> {
        let vt = &self.entries[vt_idx];
        let patch_ts = Arc::new(TileSet::new(vec![vt.rat.clone()]));
        let mf = MatchFinder::crossing(Arc::clone(&patch_ts), Arc::clone(&self.seed_tileset));

        let n = vt.rat.len();
        let mut extensions = Vec::new();

        for j in 0..mf.num_tiles_b() {
            for mt in mf.valid_matches(0, j) {
                if !in_cyclic_range(vt.vertex_pos, mt.start_a, mt.len, n) {
                    continue;
                }

                let signed_id = match self.resolve_signed_id(&vt.boundary_map, &mt) {
                    Some(id) => id,
                    None => continue,
                };

                let new_rat = mf.apply_match(&mt);
                let new_bm = update_boundary_map(
                    &vt.boundary_map,
                    &vt.rat,
                    &mt,
                    self.seed_tileset.rat(j),
                    j,
                );

                let n_patch = vt.rat.len();
                let _n_tile = self.seed_tileset.rat(j).len();
                let mlen = mt.len;

                let (is_closed, new_vpos) = if mt.start_a == vt.vertex_pos {
                    (false, n_patch - mlen)
                } else if (mt.start_a + mt.len - 1) % n_patch == vt.vertex_pos {
                    (false, 0)
                } else {
                    (true, 0)
                };

                let mut new_matches = vt.matches.clone();
                new_matches.push(signed_id);

                if is_closed {
                    let rotation = lex_min_rotation(&new_matches);
                    new_matches.rotate_left(rotation);
                    if new_vpos != 0 {
                        unreachable!("closed vertex should have vpos 0");
                    }
                }

                extensions.push(Extension {
                    new_matches,
                    rat: new_rat,
                    vertex_pos: new_vpos,
                    boundary_map: new_bm,
                    is_closed,
                });
            }
        }

        extensions
    }

    fn resolve_signed_id(&self, boundary_map: &BoundaryMap, mt: &MatchType) -> Option<i32> {
        let (orig_tile, orig_pos) = boundary_map[mt.start_a];

        let seed_a = self.seed_tileset.rat(orig_tile);
        let seed_b = self.seed_tileset.rat(mt.tile_b);
        let (ns, seed_len, ne) = seed_a.get_match((orig_pos as i64, mt.start_b as i64), seed_b);

        if seed_len == 0 {
            return None;
        }

        let orig_mt = MatchType {
            tile_a: orig_tile,
            start_a: ns as usize,
            tile_b: mt.tile_b,
            start_b: ne as usize,
            len: seed_len,
        };

        self.match_type_index.signed_id(&orig_mt)
    }

    pub fn num_types(&self) -> usize {
        self.entries.len()
    }

    pub fn num_open(&self) -> usize {
        self.entries.iter().filter(|e| !e.is_closed).count()
    }

    pub fn num_closed(&self) -> usize {
        self.entries.iter().filter(|e| e.is_closed).count()
    }

    pub fn num_terminal(&self) -> usize {
        self.entries.iter().filter(|e| e.is_terminal).count()
    }

    #[cfg(test)]
    fn entries(&self) -> &[VertexTypeEntry<T>] {
        &self.entries
    }
}

struct Extension<T: IsComplex> {
    new_matches: Vec<i32>,
    rat: Rat<T>,
    vertex_pos: usize,
    boundary_map: BoundaryMap,
    is_closed: bool,
}

fn update_boundary_map<T: IsComplex + IsRingOrField + Units>(
    old_bm: &BoundaryMap,
    old_rat: &Rat<T>,
    mt: &MatchType,
    new_tile: &Rat<T>,
    new_tile_id: usize,
) -> BoundaryMap {
    let n_patch = old_rat.len();
    let n_tile = new_tile.len();
    let mlen = mt.len;

    let mut bm = Vec::with_capacity(n_patch + n_tile - 2 * mlen);
    for i in 0..n_patch - mlen {
        bm.push(old_bm[(mt.start_a + mlen + i) % n_patch]);
    }
    for j in 0..n_tile - mlen {
        bm.push((new_tile_id, (mt.start_b + j) % n_tile));
    }

    bm
}

fn in_cyclic_range(v: usize, start: usize, len: usize, n: usize) -> bool {
    if len == 0 {
        return false;
    }
    (0..len).any(|i| (start + i) % n == v)
}

fn lex_min_rotation(seq: &[i32]) -> usize {
    if seq.len() <= 1 {
        return 0;
    }
    let n = seq.len();
    let mut best = 0;
    for i in 1..n {
        for j in 0..n {
            let a = seq[(best + j) % n];
            let b = seq[(i + j) % n];
            if b < a {
                best = i;
                break;
            }
            if b > a {
                break;
            }
        }
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic::ZZ12;
    use crate::intgeom::tiles::{hexagon, square};

    #[test]
    fn lex_min_rotation_trivial() {
        assert_eq!(lex_min_rotation(&[1, 2, 3]), 0);
        assert_eq!(lex_min_rotation(&[3, 1, 2]), 1);
        assert_eq!(lex_min_rotation(&[2, 3, 1]), 2);
        assert_eq!(lex_min_rotation(&[5]), 0);
        assert_eq!(lex_min_rotation(&[]), 0);
    }

    #[test]
    fn lex_min_rotation_uniform() {
        assert_eq!(lex_min_rotation(&[1, 1, 1]), 0);
    }

    #[test]
    fn in_cyclic_range_basic() {
        assert!(in_cyclic_range(2, 1, 3, 6));
        assert!(in_cyclic_range(1, 1, 3, 6));
        assert!(!in_cyclic_range(4, 1, 3, 6));
    }

    #[test]
    fn in_cyclic_range_wrapping() {
        assert!(in_cyclic_range(0, 5, 3, 6));
        assert!(in_cyclic_range(5, 5, 3, 6));
        assert!(!in_cyclic_range(3, 5, 3, 6));
    }

    #[test]
    fn hexagon_initial_vertex_types() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&hexagon());
        let ts = Arc::new(TileSet::new(vec![r]));
        let mti = MatchTypeIndex::new(Arc::clone(&ts));
        let vti = VertexTypeIndex::new(ts);

        let expected_initial = 2 * mti.num_types();
        assert!(
            vti.num_types() >= expected_initial,
            "hexagon should have at least {} vertex types (2x {} match types), got {}",
            expected_initial,
            mti.num_types(),
            vti.num_types()
        );
        eprintln!(
            "hexagon: {} types total, {} open, {} closed, {} terminal",
            vti.num_types(),
            vti.num_open(),
            vti.num_closed(),
            vti.num_terminal()
        );
    }

    #[test]
    fn square_initial_vertex_types() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&square());
        let ts = Arc::new(TileSet::new(vec![r]));
        let mti = MatchTypeIndex::new(Arc::clone(&ts));
        let vti = VertexTypeIndex::new(ts);

        let expected_initial = 2 * mti.num_types();
        assert!(
            vti.num_types() >= expected_initial,
            "square should have at least {} vertex types, got {}",
            expected_initial,
            vti.num_types()
        );
    }

    #[test]
    fn hexagon_vertex_type_details() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&hexagon());
        let ts = Arc::new(TileSet::new(vec![r]));
        let vti = VertexTypeIndex::new(ts);

        eprintln!(
            "hexagon: {} types total, {} open, {} closed, {} terminal",
            vti.num_types(),
            vti.num_open(),
            vti.num_closed(),
            vti.num_terminal()
        );
        for entry in vti.entries() {
            eprintln!(
                "  matches={:?} vpos={} closed={} terminal={} rat_len={}",
                entry.matches,
                entry.vertex_pos,
                entry.is_closed,
                entry.is_terminal,
                entry.rat.len()
            );
        }
    }

    #[test]
    fn hexagon_debug_extension_matches() {
        let r: Rat<ZZ12> = Rat::from_unchecked(&hexagon());
        let ts = Arc::new(TileSet::new(vec![r.clone()]));
        let mti = MatchTypeIndex::new(Arc::clone(&ts));

        let first = mti.get(1);
        let rat = mti.apply(1);
        eprintln!("match type 1: {:?}", first);
        eprintln!("result rat seq: {:?}", rat.seq());

        let patch_ts = Arc::new(TileSet::new(vec![rat.clone()]));
        let mf = MatchFinder::crossing(Arc::clone(&patch_ts), Arc::clone(&ts));

        let n = rat.len();
        let vpos = 5usize;

        let n_a = 6;
        let mlen = first.len;
        let mut bm = Vec::with_capacity(rat.len());
        for i in 0..n_a - mlen {
            bm.push((first.tile_a, (first.start_a + mlen + i) % n_a));
        }
        for j in 0..6 - mlen {
            bm.push((first.tile_b, (first.start_b + j) % 6));
        }

        for j in 0..mf.num_tiles_b() {
            for mt in mf.valid_matches(0, j) {
                if !in_cyclic_range(vpos, mt.start_a, mt.len, n) {
                    continue;
                }
                let (orig_tile, orig_pos) = bm[mt.start_a];
                let seed_a = ts.rat(orig_tile);
                let seed_b = ts.rat(mt.tile_b);
                let (ns, seed_len, ne) =
                    seed_a.get_match((orig_pos as i64, mt.start_b as i64), seed_b);
                let orig_mt = MatchType {
                    tile_a: orig_tile,
                    start_a: ns as usize,
                    tile_b: mt.tile_b,
                    start_b: ne as usize,
                    len: seed_len,
                };
                let sid = mti.signed_id(&orig_mt);

                let new_rat = mf.apply_match(&mt);
                let is_closed = !(mt.start_a == vpos || (mt.start_a + mt.len - 1) % n == vpos);
                let new_vpos = if mt.start_a == vpos { n - mt.len } else { 0 };

                eprintln!(
                    "  patch_match=({},{},{}) seed_match=({},{},{}) sid={:?} new_rat_len={} closed={} new_vpos={}",
                    mt.start_a, mt.start_b, mt.len,
                    ns, ne, seed_len,
                    sid, new_rat.len(), is_closed, new_vpos
                );
            }
        }
    }
}
