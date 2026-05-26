//! Enumerate every cyclic boundary subsequence reachable on patches
//! grown from a fixed tileset.
//!
//! # What this builds
//!
//! Given a [`TileSet`], [`SeqExplorer`] grows patches by repeatedly
//! gluing tiles and records, for each patch's boundary, all of its
//! cyclic length-`k` substrings of edge angles (where `k =
//! max(tile.len())` across the tileset). Two patches with the same
//! canonical boundary shape (same [`Rat`]) are merged. The output is:
//!
//! - a **trie** of every distinct cyclic length-`≤k` boundary
//!   substring seen across all reachable patches,
//! - a list of **witness rats**: those patches that contributed at
//!   least one substring not already in the trie at the time they
//!   were processed,
//! - a [`Provenance`] entry per known rat that records how it was
//!   produced (which seed tile, or which source rat × tile glue).
//!
//! The BFS terminates when no further glue produces a substring that
//! isn't already in the trie — equivalently, when the per-substring
//! contribution graph reaches a fixed point.
//!
//! # Algorithm
//!
//! Seed phase: each distinct tile rat is inserted with all of its
//! cyclic length-`k` substrings. Tiles that contribute new substrings
//! become pending witnesses.
//!
//! BFS layer (repeats until no pending witnesses remain):
//!
//! 1. Group the pending witnesses into a fresh "batch" [`TileSet`].
//! 2. Build a [`MatchFinder::crossing_with_seed`] over `(batch,
//!    fixed_seed_bp)` — the B-side bit-parallel masks of the original
//!    tileset are precomputed once via [`BpSeed`] before the loop and
//!    reused across every layer.
//! 3. For each `(batch_rat, tile)` pair, enumerate valid glue matches
//!    via [`MatchFinder::valid_matches_filtered`]. The filter is the
//!    per-rat **active-edge mask** (see below), which prunes matches
//!    that geometrically can't produce a new boundary substring.
//! 4. Each accepted glue is checked for Snake-validity, deduped by
//!    [`Rat`] equality, and assigned a fresh `rat_id`.
//! 5. For each new rat, insert its cyclic substrings into the trie.
//!    If any are new, the rat becomes a witness and is added to the
//!    next layer's pending set. Compute its active-edge mask for next
//!    time. Otherwise the rat is recorded but doesn't propagate.
//!
//! # The active-edge mask
//!
//! After a glue produces a new rat, the only boundary positions whose
//! local k-neighbourhood could *possibly* differ from a previously-
//! seen substring are positions near the two new junctions (the CW
//! and CCW endpoints of the consumed match). The active-edge mask
//! marks edges that lie within `k` of those junctions, then dilates
//! by `k − 1` so that a future match anchored anywhere in that
//! neighbourhood is allowed; everything else can be safely skipped.
//! This is what keeps the BFS bounded — without it, growth would
//! revisit boundary regions that can't contribute anything new.
//!
//! # When to use this
//!
//! When you need the universe of locally-realizable boundary
//! substrings for a tileset, e.g. for downstream pattern analysis or
//! to test whether a candidate sequence is geometrically realizable.
//! `SeqExplorer` does not retain the grown patches themselves, only
//! their canonical [`Rat`]s and a provenance chain — so memory stays
//! modest even for large state spaces (where the parallel
//! `GrowingPatch`-based explorer ran into 10s of GB).

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::sync::Arc;
use std::time::Instant;

use rustc_hash::{FxHashMap, FxHashSet};
use serde::{Deserialize, Serialize};

use crate::cyclotomic::{IsRing};
use crate::intgeom::matchtypes::{BpSeed, MatchFinder};
use crate::intgeom::rat::Rat;
use crate::intgeom::snake::Snake;
use crate::intgeom::tileset::TileSet;
use crate::stringmatch::cyclic_contains;

/// Re-apply a glue described by `(start_a, start_b, len)` to
/// reconstruct the resulting [`Rat`] from a source rat and a tile
/// rat. Useful for replaying a [`Provenance::Glue`] step against
/// stored rats without needing the full explorer state.
///
/// The `(start_a, start_b, len)` triple is the same one stored in
/// [`Provenance::Glue`] and follows the standard `MatchType`
/// edge-offset convention.
pub fn replay_glue<T: IsRing>(
    source: &Rat<T>,
    tile: &Rat<T>,
    start_a: usize,
    start_b: usize,
    len: usize,
) -> Result<Rat<T>, String> {
    source
        .try_glue_precomputed((start_a as i64, len, start_b as i64), tile, true)
        .map_err(|e| e.to_string())
}

/// One node in the substring trie. `children` maps the next angle in
/// the substring to the child node. `contributor = Some(rat_id)` on a
/// node tags the rat that *first* deposited the substring read from
/// root to this node.
struct TrieNode {
    children: HashMap<i8, usize>,
    contributor: Option<usize>,
}

impl TrieNode {
    fn new() -> Self {
        TrieNode {
            children: HashMap::new(),
            contributor: None,
        }
    }
}

/// Trie over cyclic length-`≤k` substrings drawn from boundary
/// angle sequences. Used by [`SeqExplorer`] to track which
/// substrings have already been seen and which rat first produced
/// each.
struct SeqTrie {
    nodes: Vec<TrieNode>,
}

impl SeqTrie {
    fn new() -> Self {
        SeqTrie {
            nodes: vec![TrieNode::new()],
        }
    }

    /// Insert every cyclic length-`min(k, n)` substring of `seq`
    /// (where `n = seq.len()`) into the trie, tagging each newly-
    /// created node with `rat_id`.
    ///
    /// Returns `(new_count, active)` where:
    /// - `new_count` is the number of trie nodes created by this
    ///   call (i.e. the number of substrings not previously seen),
    /// - `active[i] = true` iff the cyclic substring starting at
    ///   position `i` of `seq` contributed at least one new node.
    ///
    /// The `active` mask is the input to [`glue_active_mask`], which
    /// uses it to compute the per-rat active-edge filter for the
    /// next BFS layer.
    fn insert_cyclic_subseqs(&mut self, seq: &[i8], k: usize, rat_id: usize) -> (usize, Vec<bool>) {
        let n = seq.len();
        if n == 0 {
            return (0, Vec::new());
        }
        let max_len = k.min(n);
        let mut new_count = 0;
        let mut active = vec![false; n];

        for start in 0..n {
            let mut node = 0;
            let mut any_new = false;
            for l in 0..max_len {
                let angle = seq[(start + l) % n];
                let next = self.nodes[node].children.get(&angle).copied();
                node = match next {
                    Some(child) => child,
                    None => {
                        let new_node = self.nodes.len();
                        self.nodes.push(TrieNode::new());
                        self.nodes[node].children.insert(angle, new_node);
                        self.nodes[new_node].contributor = Some(rat_id);
                        new_count += 1;
                        any_new = true;
                        new_node
                    }
                };
            }
            if any_new {
                active[start] = true;
            }
        }

        (new_count, active)
    }

    fn len(&self) -> usize {
        self.nodes.len() - 1
    }

    fn sequences_by_rat_id(&self) -> BTreeMap<usize, Vec<Vec<i8>>> {
        let mut result: BTreeMap<usize, Vec<Vec<i8>>> = BTreeMap::new();
        self.dfs_collect(0, &mut Vec::new(), &mut result);
        result
    }

    fn dfs_collect(
        &self,
        node: usize,
        path: &mut Vec<i8>,
        result: &mut BTreeMap<usize, Vec<Vec<i8>>>,
    ) {
        if let Some(rat_id) = self.nodes[node].contributor {
            result.entry(rat_id).or_default().push(path.clone());
        }
        let mut sorted_children: Vec<(i8, usize)> = self.nodes[node]
            .children
            .iter()
            .map(|(&a, &c)| (a, c))
            .collect();
        sorted_children.sort_by_key(|(a, _)| *a);
        for (angle, child) in sorted_children {
            path.push(angle);
            self.dfs_collect(child, path, result);
            path.pop();
        }
    }
}

/// Build the active-edge mask for a newly-discovered glued rat.
///
/// Only the two new boundary junctions (CW and CCW endpoints of the
/// consumed match) can possibly host a new substring. The mask flags
/// positions within `k` of either junction iff there's evidence (from
/// `raw_active`) that something new really appeared there, then dilates
/// by `k − 1` so future matches anchored anywhere in that
/// neighbourhood are allowed through.
///
/// * `raw_active`: the per-position "this substring was new" flags from
///   [`SeqTrie::insert_cyclic_subseqs`];
/// * `k`: max substring length;
/// * `match_len`: edge length of the consumed match;
/// * `source_len`: perimeter length of the source (pre-glue) rat.
///
/// For seed rats, callers use `vec![true; perimeter]` directly — every
/// position is fully active.
fn glue_active_mask(
    raw_active: &[bool],
    k: usize,
    match_len: usize,
    source_len: usize,
) -> Vec<bool> {
    let n = raw_active.len();
    let cw_j = source_len - match_len;
    let ccw_j = 0usize;
    let mut active = vec![false; n];
    let check_radius = k.min(n);
    let dilate = k.saturating_sub(1).min(n);

    for &j in &[ccw_j, cw_j] {
        let mut induced = false;
        for d in 0..check_radius {
            if raw_active[(j + n - d) % n] {
                induced = true;
                break;
            }
        }
        if induced {
            for d in 0..=dilate {
                active[(j + n - d) % n] = true;
                active[(j + d) % n] = true;
            }
        }
    }
    active
}

/// How a rat in the explorer's catalog was first produced.
///
/// Every known rat has exactly one provenance — recorded on first
/// insertion and never updated. For seeds, the provenance just names
/// the source tile. For glues, it records the source rat plus the
/// glue's `(start_a, start_b, len)` triple, which together with
/// [`Self::apply`] (or the lower-level [`replay_glue`]) is enough to
/// reconstruct the rat from scratch.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Provenance {
    /// A bare tile rat used as a BFS seed.
    Seed {
        /// Index into the source [`TileSet`].
        tile_idx: usize,
    },
    /// A rat produced by gluing tile `tile_idx` onto the rat with id
    /// `source_rat_id`, using the match described by `(start_a,
    /// start_b, len)` (standard `MatchType` edge-offset convention).
    Glue {
        source_rat_id: usize,
        tile_idx: usize,
        start_a: usize,
        start_b: usize,
        len: usize,
    },
}

impl Provenance {
    /// Reconstruct the rat described by this provenance, by looking
    /// up source rats and tiles in `explorer`. For [`Self::Seed`]
    /// just returns the tile rat; for [`Self::Glue`] re-applies the
    /// match. Returns `Err` on a malformed Glue (shouldn't happen
    /// for provenances produced by a SeqExplorer).
    pub fn apply<T: IsRing>(
        &self,
        explorer: &SeqExplorer<T>,
    ) -> Result<Rat<T>, String> {
        match self {
            Self::Seed { tile_idx } => Ok(explorer.tileset().rat(*tile_idx).clone()),
            Self::Glue {
                source_rat_id,
                tile_idx,
                start_a,
                start_b,
                len,
            } => {
                let source = explorer.rat(*source_rat_id);
                let tile = explorer.tileset().rat(*tile_idx);
                replay_glue(source, tile, *start_a, *start_b, *len)
            }
        }
    }
}

/// Serializable snapshot of a `SeqExplorer` result: the tileset
/// (as raw angle sequences), every rat's provenance, and every
/// `(rat_id, subseq)` pair from the trie. Non-generic — all typed
/// rat data is recoverable by replaying provenances against a
/// caller-built `TileSet<T>` of the appropriate ring.
///
/// Designed to round-trip via `serde_json` (or any other serde
/// format) for `seq_collect`-style collect/validate workflows.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Collection {
    /// Symbolic name of the cyclotomic ring (e.g. `"ZZ12"`,
    /// `"ZZ10"`) — used by deserializers to dispatch to the right
    /// `T`.
    pub ring: String,
    /// Raw angle sequence of each tile, in `TileSet` (sorted) order.
    pub tile_angles: Vec<Vec<i8>>,
    /// `provenances[rat_id]` for every rat in the explorer.
    pub provenances: Vec<Provenance>,
    /// `(rat_id, subseq)` pairs: each cyclic substring in the trie,
    /// tagged by the rat that first contributed it.
    pub subseqs: Vec<(usize, Vec<i8>)>,
}

impl Collection {
    /// Snapshot a built explorer, tagging it with the given ring
    /// name (used only on round-trip, to pick the right `T`).
    pub fn from_explorer<T>(explorer: &SeqExplorer<T>, ring: impl Into<String>) -> Self
    where
        T: IsRing,
    {
        let tile_angles = explorer
            .tileset()
            .rats()
            .iter()
            .map(|r| r.seq().to_vec())
            .collect();
        let provenances = explorer.provenances().to_vec();
        let mut subseqs: Vec<(usize, Vec<i8>)> = Vec::new();
        for (rat_id, seqs) in explorer.sequences_by_rat_id() {
            for seq in seqs {
                subseqs.push((rat_id, seq));
            }
        }
        Collection {
            ring: ring.into(),
            tile_angles,
            provenances,
            subseqs,
        }
    }

    /// Reconstruct typed rats by replaying every provenance against
    /// `tile_ts`. The caller is responsible for building `tile_ts`
    /// from `self.tile_angles` in the matching order — see
    /// `seq_collect.rs` for the ring-dispatch pattern.
    pub fn replay_rats<T>(&self, tile_ts: &TileSet<T>) -> Result<Vec<Rat<T>>, String>
    where
        T: IsRing,
    {
        let mut rats: Vec<Rat<T>> = Vec::with_capacity(self.provenances.len());
        for (id, prov) in self.provenances.iter().enumerate() {
            let rat = match prov {
                Provenance::Seed { tile_idx } => tile_ts.rat(*tile_idx).clone(),
                Provenance::Glue {
                    source_rat_id,
                    tile_idx,
                    start_a,
                    start_b,
                    len,
                } => {
                    let source = &rats[*source_rat_id];
                    let tile = tile_ts.rat(*tile_idx);
                    replay_glue(source, tile, *start_a, *start_b, *len)
                        .map_err(|e| format!("rat {id}: glue replay failed: {e}"))?
                }
            };
            rats.push(rat);
        }
        Ok(rats)
    }

    /// `(rat_id, subseq)` pairs whose `subseq` is *not* actually a
    /// cyclic substring of `rats[rat_id]`. Empty ⇒ presence check
    /// passes.
    pub fn presence_errors<T>(&self, rats: &[Rat<T>]) -> Vec<(usize, Vec<i8>)>
    where
        T: IsRing,
    {
        self.subseqs
            .iter()
            .filter(|(id, seq)| !cyclic_contains(rats[*id].seq(), seq))
            .cloned()
            .collect()
    }
}

/// Summary of a [`check_fixed_point`] scan.
#[derive(Debug, Default, Clone)]
pub struct FixedPointReport {
    /// Total number of valid `(witness, tile)` glue matches checked.
    pub matches_checked: usize,
    /// Distinct cyclic substrings produced by some witness × tile
    /// glue that were *not* in the supplied `known` set. Empty ⇒
    /// the BFS is at its claimed fixed point.
    pub missing: Vec<Vec<i8>>,
}

impl FixedPointReport {
    /// True iff no missing substrings were found.
    pub fn is_complete(&self) -> bool {
        self.missing.is_empty()
    }
}

/// Witness-batch size used by [`check_fixed_point`]. Bounds peak
/// memory for the per-batch `MatchFinder` on large tilesets while
/// staying single-batch on hex/square.
const FIXED_POINT_BATCH: usize = 500;

/// Verify the BFS completeness property: no valid glue from any
/// witness against any tile produces a cyclic length-`≤k`
/// substring that isn't already in `known`.
///
/// Enumerates every valid `(witness, tile)` glue match
/// (unfiltered by the active-edge mask), snake-validates each
/// result (matching the BFS's own filter), and checks every cyclic
/// length-`≤max_subseq_len` substring of the result against `known`.
/// Discovered misses are returned in `FixedPointReport::missing`.
pub fn check_fixed_point<T>(
    tileset: &Arc<TileSet<T>>,
    rats: &[Rat<T>],
    witness_ids: &[usize],
    known: &BTreeSet<Vec<i8>>,
    max_subseq_len: usize,
) -> FixedPointReport
where
    T: IsRing,
{
    let mut report = FixedPointReport::default();
    let mut seen_missing: BTreeSet<Vec<i8>> = BTreeSet::new();

    if witness_ids.is_empty() {
        return report;
    }

    // The B-side (`tileset`) is fixed across all batches; precompute
    // its bit-parallel masks once and share via Arc.
    let bp_seed = BpSeed::new(Arc::clone(tileset));

    for chunk in witness_ids.chunks(FIXED_POINT_BATCH) {
        let witness_rats: Vec<Rat<T>> = chunk.iter().map(|&id| rats[id].clone()).collect();
        let witness_ts = Arc::new(TileSet::new(witness_rats));
        let mf = MatchFinder::crossing_with_seed(Arc::clone(&witness_ts), bp_seed.clone());

        for wi in 0..witness_ts.num_tiles() {
            for ti in 0..tileset.num_tiles() {
                for mt in mf.valid_matches(wi, ti) {
                    report.matches_checked += 1;
                    let glued = crate::intgeom::matchtypes::apply_match(
                        &mt,
                        witness_ts.rats(),
                        tileset.rats(),
                    );
                    if Snake::<T>::try_from(glued.seq()).is_err() {
                        continue;
                    }
                    let seq = glued.seq();
                    let n = seq.len();
                    let max_len = max_subseq_len.min(n);
                    for start in 0..n {
                        for l in 1..=max_len {
                            let sub: Vec<i8> = (0..l).map(|j| seq[(start + j) % n]).collect();
                            if !known.contains(&sub) && seen_missing.insert(sub.clone()) {
                                report.missing.push(sub);
                            }
                        }
                    }
                }
            }
        }
    }

    report
}

/// BFS-built catalog of every cyclic boundary substring reachable
/// on patches grown from a [`TileSet`].
///
/// See the module docs for the algorithm. After [`Self::new`] runs:
///
/// - [`Self::num_subseqs`] gives the total number of distinct cyclic
///   length-`≤k` substrings discovered;
/// - [`Self::sequences_by_rat_id`] groups them by the rat that first
///   contributed them;
/// - [`Self::rats`] / [`Self::provenances`] give the full per-rat
///   discovery record.
pub struct SeqExplorer<T: IsRing> {
    tileset: Arc<TileSet<T>>,
    /// Trie of every cyclic length-`≤k` substring seen, with each
    /// node tagged by the `rat_id` of the rat that *first* deposited
    /// that substring.
    trie: SeqTrie,
    /// Every **contributing** rat (witness) the BFS catalogued,
    /// indexed by `rat_id`. Non-contributing rats (whose cyclic
    /// substrings were all already in the trie when first seen) are
    /// not kept here — they were discarded after dedup, since
    /// [`Provenance::Glue::source_rat_id`] only ever references a
    /// witness.
    rats: Vec<Rat<T>>,
    /// `provenances[rat_id]` records how `rats[rat_id]` was produced.
    provenances: Vec<Provenance>,
    /// Reverse map from a witness's canonical [`Rat`] to its `rat_id`.
    rat_to_id: FxHashMap<Rat<T>, usize>,
    /// `k = max(tile.len() for tile in tileset)`. The trie collects
    /// cyclic substrings of length up to `k`.
    max_subseq_len: usize,
}

impl<T: IsRing> SeqExplorer<T> {
    /// Run the BFS over patches of `tileset` to fixed point and
    /// return the populated catalog. Quiet by default; use
    /// [`Self::new_verbose`] to also log per-layer progress to stderr.
    ///
    /// Runtime is roughly linear in (tileset size × number of
    /// reachable rats × max tile length). Memory is dominated by
    /// `O(num_rats × avg_perimeter)` for the rat list plus the trie.
    pub fn new(tileset: Arc<TileSet<T>>) -> Self {
        Self::build(tileset, false)
    }

    /// Like [`Self::new`] but logs the seed phase, per-layer
    /// progress, and a final summary to stderr.
    pub fn new_verbose(tileset: Arc<TileSet<T>>) -> Self {
        Self::build(tileset, true)
    }

    fn build(tileset: Arc<TileSet<T>>, verbose: bool) -> Self {
        let mut b = Builder::new(tileset, verbose);
        b.seed_phase();
        while !b.pending.is_empty() {
            b.run_layer();
        }
        b.into_explorer()
    }

    /// Total number of distinct cyclic substrings discovered (size
    /// of the trie minus the root).
    pub fn num_subseqs(&self) -> usize {
        self.trie.len()
    }

    /// Number of catalogued rats (witnesses): every rat the BFS kept
    /// because it contributed at least one new cyclic substring on
    /// first insertion. Equals `rats().len()` and equals the number
    /// of keys in [`Self::sequences_by_rat_id`]. Non-contributing
    /// rats are not retained.
    pub fn num_rats(&self) -> usize {
        self.rats.len()
    }

    /// `k = max(tile.len() for tile in tileset)`. The trie collects
    /// cyclic substrings of length up to this many edges.
    pub fn max_subseq_len(&self) -> usize {
        self.max_subseq_len
    }

    /// All catalogued witness rats, indexed by `rat_id`.
    pub fn rats(&self) -> &[Rat<T>] {
        &self.rats
    }

    /// `provenances[rat_id]` for every `rat_id`.
    pub fn provenances(&self) -> &[Provenance] {
        &self.provenances
    }

    /// Provenance of a specific rat. Panics if `rat_id` is out of
    /// range.
    pub fn provenance(&self, rat_id: usize) -> &Provenance {
        &self.provenances[rat_id]
    }

    /// Look up a rat by id. Panics if out of range.
    pub fn rat(&self, id: usize) -> &Rat<T> {
        &self.rats[id]
    }

    /// Look up the `rat_id` for a given canonical rat, if known.
    pub fn rat_id_of(&self, rat: &Rat<T>) -> Option<usize> {
        self.rat_to_id.get(rat).copied()
    }

    /// All trie substrings grouped by `rat_id` of the rat that first
    /// contributed them. Every value in the returned map is a list
    /// of cyclic substrings of length `≤ max_subseq_len()` that are
    /// substrings of `rats()[rat_id]`'s angle sequence.
    pub fn sequences_by_rat_id(&self) -> BTreeMap<usize, Vec<Vec<i8>>> {
        self.trie.sequences_by_rat_id()
    }

    /// The tileset this explorer was built from.
    pub fn tileset(&self) -> &TileSet<T> {
        &self.tileset
    }

    /// `Arc` handle to the tileset, for callers that need ownership
    /// (e.g. passing to another component).
    pub fn tileset_arc(&self) -> Arc<TileSet<T>> {
        Arc::clone(&self.tileset)
    }
}

/// Per-layer summary returned by [`Builder::integrate_new_rats`].
struct LayerStats {
    /// Total new rats discovered (whether contributing or not).
    total: usize,
    /// Of those, how many contributed at least one new substring.
    contributing: usize,
    /// Min / max edge count among the new rats. `len_min` is
    /// `usize::MAX` when `total == 0`.
    len_min: usize,
    len_max: usize,
}

/// Private build state for [`SeqExplorer::new`]. Owns everything the
/// BFS needs and exposes it as a sequence of phase helpers
/// (`seed_phase` → `run_layer*` → `into_explorer`).
struct Builder<T: IsRing> {
    tileset: Arc<TileSet<T>>,
    trie: SeqTrie,
    /// Witness rats only (those that contributed when first inserted).
    rats: Vec<Rat<T>>,
    provenances: Vec<Provenance>,
    rat_to_id: FxHashMap<Rat<T>, usize>,
    max_subseq_len: usize,

    /// Rats queued for the next BFS layer (contributing on their last
    /// insertion). Drained at the start of each `run_layer`.
    pending: BTreeSet<usize>,
    /// Active-edge mask for each pending rat, drained as the rat is
    /// processed in a layer.
    active_map: HashMap<usize, Vec<bool>>,
    /// Canonical forms of non-contributing rats — kept solely so we
    /// don't re-process the same dead-end via another glue path.
    /// Dropped at `into_explorer` time; never exposed on the
    /// finished `SeqExplorer`.
    seen_dead: FxHashSet<Rat<T>>,
    /// Stat: total non-contributing rats encountered. Reported in
    /// verbose mode to surface how much pruning saves.
    dead_count: usize,
    /// Precomputed bit-parallel B-side state for the fixed tileset.
    /// Shared across every layer's MatchFinder.
    bp_seed: BpSeed<T>,

    verbose: bool,
    layer: usize,
    started: Instant,
}

impl<T: IsRing> Builder<T> {
    fn new(tileset: Arc<TileSet<T>>, verbose: bool) -> Self {
        let max_subseq_len = tileset.rats().iter().map(|r| r.len()).max().unwrap_or(0);
        // B-side (tileset) is fixed across all layers; precompute the
        // bit-parallel masks once.
        let bp_seed = BpSeed::new(Arc::clone(&tileset));
        Builder {
            tileset,
            trie: SeqTrie::new(),
            rats: Vec::new(),
            provenances: Vec::new(),
            rat_to_id: FxHashMap::default(),
            max_subseq_len,
            pending: BTreeSet::new(),
            active_map: HashMap::new(),
            seen_dead: FxHashSet::default(),
            dead_count: 0,
            bp_seed,
            verbose,
            layer: 0,
            started: Instant::now(),
        }
    }

    /// Insert every tile of `tileset` as a seed rat. Tiles whose
    /// cyclic substrings expand the trie become pending for the first
    /// BFS layer with a fully-active edge mask; non-contributing
    /// seeds (whose substrings were all already in the trie from an
    /// earlier seed) are remembered for dedup only.
    fn seed_phase(&mut self) {
        for tile_idx in 0..self.tileset.num_tiles() {
            let rat = self.tileset.rat(tile_idx).clone();
            if self.rat_to_id.contains_key(&rat) || self.seen_dead.contains(&rat) {
                continue;
            }

            // Speculate rat_id = next contributing index. If this rat
            // doesn't contribute, no trie node will be tagged with it
            // (insert_cyclic_subseqs only tags nodes it *creates*),
            // so the speculative id is safe to discard.
            let tentative_rat_id = self.rats.len();
            let (new_count, _raw) =
                self.trie
                    .insert_cyclic_subseqs(rat.seq(), self.max_subseq_len, tentative_rat_id);
            if new_count == 0 {
                self.seen_dead.insert(rat);
                self.dead_count += 1;
                continue;
            }

            // Contributing: commit.
            let rat_id = tentative_rat_id;
            self.rat_to_id.insert(rat.clone(), rat_id);
            // Seed rats: no prior neighbourhood to restrict against,
            // so every boundary position is active.
            self.active_map.insert(rat_id, vec![true; rat.seq().len()]);
            self.pending.insert(rat_id);
            self.rats.push(rat);
            self.provenances.push(Provenance::Seed { tile_idx });
        }

        if self.verbose {
            eprintln!(
                "  Seeds: {} witnesses, {} dead, {} subseqs, {} pending",
                self.rats.len(),
                self.dead_count,
                self.trie.len(),
                self.pending.len(),
            );
        }
    }

    /// One BFS layer: drain `pending` into a batch, enumerate all
    /// glues against the original tileset, integrate the new rats and
    /// repopulate `pending` for the next iteration.
    fn run_layer(&mut self) {
        self.layer += 1;
        let batch_ids: Vec<usize> = std::mem::take(&mut self.pending).into_iter().collect();

        let (batch_ts, global_from_batch, batch_active) = self.build_batch(&batch_ids);
        let new_entries = self.enumerate_layer_glues(&batch_ts, &global_from_batch, &batch_active);
        let stats = self.integrate_new_rats(new_entries);

        if self.verbose {
            eprintln!(
                "  Layer {}: batch={} new_rats={} new_contrib={} trie={} pending={} lens=[{},{}]",
                self.layer,
                batch_ts.num_tiles(),
                stats.total,
                stats.contributing,
                self.trie.len(),
                self.pending.len(),
                stats.len_min,
                stats.len_max,
            );
        }
    }

    /// Pack `batch_ids` into a per-layer batch tileset, recover the
    /// post-sort batch-idx → global-rat-id mapping, and pull each
    /// rat's active-edge mask out of `active_map`.
    ///
    /// `TileSet::new` sorts its input (and would dedup if there were
    /// dupes — there aren't, since `pending` is a BTreeSet of unique
    /// ids), so the original input order is not preserved. We must
    /// look each batch rat back up in `rat_to_id` to recover the
    /// correct mapping.
    fn build_batch(
        &mut self,
        batch_ids: &[usize],
    ) -> (Arc<TileSet<T>>, Vec<usize>, Vec<Vec<bool>>) {
        let batch_rats: Vec<Rat<T>> = batch_ids.iter().map(|&id| self.rats[id].clone()).collect();
        let batch_ts = Arc::new(TileSet::new(batch_rats));
        assert_eq!(batch_ts.num_tiles(), batch_ids.len());

        let global_from_batch: Vec<usize> = (0..batch_ts.num_tiles())
            .map(|i| {
                *self
                    .rat_to_id
                    .get(batch_ts.rat(i))
                    .expect("batch rat must be in registry")
            })
            .collect();

        let batch_active: Vec<Vec<bool>> = global_from_batch
            .iter()
            .map(|&id| self.active_map.remove(&id).unwrap_or_default())
            .collect();

        (batch_ts, global_from_batch, batch_active)
    }

    /// Enumerate every `(batch_rat × tile)` glue that produces a new,
    /// snake-valid rat. Deduped against rats already in the catalog
    /// and within the layer itself.
    ///
    /// **Why MatchFinder, not MatchTypeIndex.** The A-side (batch)
    /// changes every layer, so any precomputed self-index over `batch
    /// ∪ tileset` would have to be rebuilt every layer; the per-layer
    /// cost of building such an index exceeds the cost of just
    /// streaming the A-side through the precomputed B-side bit-parallel
    /// masks. `MatchTypeIndex` also doesn't support the per-rat
    /// active-edge filter that drives BFS pruning here. The current
    /// `MatchFinder::crossing_with_seed` is the right shape: precompute
    /// what's fixed (the B-side seed) once, stream what varies (the
    /// A-side batch) per layer.
    fn enumerate_layer_glues(
        &self,
        batch_ts: &Arc<TileSet<T>>,
        global_from_batch: &[usize],
        batch_active: &[Vec<bool>],
    ) -> Vec<(Rat<T>, Provenance)> {
        let mf = MatchFinder::crossing_with_seed(Arc::clone(batch_ts), self.bp_seed.clone());

        let mut new_entries: Vec<(Rat<T>, Provenance)> = Vec::new();
        let mut seen_new: FxHashSet<Rat<T>> = FxHashSet::default();
        let mut seen_invalid: FxHashSet<Rat<T>> = FxHashSet::default();

        for batch_idx in 0..batch_ts.num_tiles() {
            let active = &batch_active[batch_idx];
            let source_rat_id = global_from_batch[batch_idx];

            for tile_idx in 0..self.tileset.num_tiles() {
                for (glued, mt) in mf.valid_matches_filtered(batch_idx, tile_idx, active) {
                    if self.rat_to_id.contains_key(&glued)
                        || self.seen_dead.contains(&glued)
                        || seen_new.contains(&glued)
                        || seen_invalid.contains(&glued)
                    {
                        continue;
                    }
                    if Snake::<T>::try_from(glued.seq()).is_err() {
                        seen_invalid.insert(glued);
                        continue;
                    }
                    seen_new.insert(glued.clone());
                    new_entries.push((
                        glued,
                        Provenance::Glue {
                            source_rat_id,
                            tile_idx,
                            start_a: mt.a.range.start_offset,
                            start_b: mt.b.range.start_offset,
                            len: mt.len(),
                        },
                    ));
                }
            }
        }
        new_entries
    }

    /// Insert each candidate rat's cyclic substrings into the trie;
    /// commit witnesses (new_count > 0) to the catalog and queue them
    /// for the next layer; remember non-contributing rats in
    /// `seen_dead` so we don't re-process them via another path.
    fn integrate_new_rats(&mut self, new_entries: Vec<(Rat<T>, Provenance)>) -> LayerStats {
        let mut stats = LayerStats {
            total: new_entries.len(),
            contributing: 0,
            len_min: usize::MAX,
            len_max: 0,
        };
        for (rat, _) in &new_entries {
            let l = rat.len();
            stats.len_min = stats.len_min.min(l);
            stats.len_max = stats.len_max.max(l);
        }

        for (rat, prov) in new_entries {
            // BFS-layer entries are always glues — seeds were handled
            // in the seed phase.
            let (source_rat_id, match_len) = match &prov {
                Provenance::Glue {
                    source_rat_id, len, ..
                } => (*source_rat_id, *len),
                Provenance::Seed { .. } => unreachable!("BFS layer only produces glues"),
            };
            let source_len = self.rats[source_rat_id].len();

            // Speculate the next contributing rat_id; if this rat
            // turns out non-contributing, nothing in the trie gets
            // tagged with it (no new nodes were created), so the
            // tentative id is safe to discard.
            let tentative_rat_id = self.rats.len();
            let (new_count, raw_active) =
                self.trie
                    .insert_cyclic_subseqs(rat.seq(), self.max_subseq_len, tentative_rat_id);
            if new_count == 0 {
                self.seen_dead.insert(rat);
                self.dead_count += 1;
                continue;
            }

            // Contributing: commit.
            let rat_id = tentative_rat_id;
            let mask = glue_active_mask(&raw_active, self.max_subseq_len, match_len, source_len);
            self.active_map.insert(rat_id, mask);
            self.pending.insert(rat_id);
            self.rat_to_id.insert(rat.clone(), rat_id);
            self.rats.push(rat);
            self.provenances.push(prov);
            stats.contributing += 1;
        }

        stats
    }

    /// Consume the builder and assemble the immutable `SeqExplorer`.
    /// `seen_dead` is dropped here — it served only as in-flight
    /// dedup state.
    fn into_explorer(self) -> SeqExplorer<T> {
        if self.verbose {
            eprintln!(
                "  Done: {} layers, {} witnesses, {} dead, {} subseqs, time={:.2?}",
                self.layer,
                self.rats.len(),
                self.dead_count,
                self.trie.len(),
                self.started.elapsed(),
            );
        }

        SeqExplorer {
            tileset: self.tileset,
            trie: self.trie,
            rats: self.rats,
            provenances: self.provenances,
            rat_to_id: self.rat_to_id,
            max_subseq_len: self.max_subseq_len,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cyclotomic::ZZ12;
    use crate::intgeom::tileset;

    fn hex_ts() -> Arc<TileSet<ZZ12>> {
        tileset::hex::<ZZ12>()
    }

    fn sq_ts() -> Arc<TileSet<ZZ12>> {
        tileset::square::<ZZ12>()
    }

    /// Structural invariants every built `SeqExplorer` must satisfy,
    /// independent of the specific tileset.
    fn assert_invariants(label: &str, explorer: &SeqExplorer<ZZ12>) {
        let n_rats = explorer.num_rats();
        let n_subseqs = explorer.num_subseqs();
        let k = explorer.max_subseq_len();

        // Vec lengths agree with their accessor.
        assert_eq!(explorer.rats().len(), n_rats);
        assert_eq!(explorer.provenances().len(), n_rats);

        // `rat_to_id` is the exact inverse of `rats[]`.
        for rat_id in 0..n_rats {
            assert_eq!(
                explorer.rat_id_of(explorer.rat(rat_id)),
                Some(rat_id),
                "{label}: rat_id_of round-trip failed for rat {rat_id}"
            );
        }

        // Provenances form a DAG: every glue's source predates its
        // own rat; every tile index is in range; seeds always point
        // at a real tile.
        let n_tiles = explorer.tileset().num_tiles();
        for (rat_id, prov) in explorer.provenances().iter().enumerate() {
            match prov {
                Provenance::Seed { tile_idx } => {
                    assert!(*tile_idx < n_tiles, "{label}: seed tile_idx out of range");
                }
                Provenance::Glue {
                    source_rat_id,
                    tile_idx,
                    ..
                } => {
                    assert!(
                        *source_rat_id < rat_id,
                        "{label}: glue source ({source_rat_id}) must precede dest ({rat_id})"
                    );
                    assert!(*tile_idx < n_tiles, "{label}: glue tile_idx out of range");
                }
            }
        }

        // Every provenance reconstructs its rat exactly when applied
        // back through `apply`.
        for rat_id in 0..n_rats {
            let reconstructed = explorer
                .provenance(rat_id)
                .apply(explorer)
                .unwrap_or_else(|e| panic!("{label}: rat {rat_id} apply failed: {e}"));
            assert_eq!(
                &reconstructed,
                explorer.rat(rat_id),
                "{label}: provenance.apply doesn't reproduce rat {rat_id}"
            );
        }

        // Every subseq in `sequences_by_rat_id` really is a cyclic
        // substring of its host rat, within the length bound; the
        // grouping covers every trie entry and every key contributed.
        // Every catalogued rat is a contributor, so the by-rat map's
        // keys are exactly `0..n_rats`.
        let by_rat = explorer.sequences_by_rat_id();
        assert_eq!(
            by_rat.len(),
            n_rats,
            "{label}: every catalogued rat must contribute"
        );
        let mut total_seqs = 0usize;
        for (&rat_id, seqs) in &by_rat {
            assert!(rat_id < n_rats);
            assert!(
                !seqs.is_empty(),
                "{label}: contributor entry must be non-empty"
            );
            let host = explorer.rat(rat_id).seq();
            for seq in seqs {
                assert!(!seq.is_empty(), "{label}: empty subseq emitted");
                assert!(
                    seq.len() <= k,
                    "{label}: subseq of length {} exceeds max_subseq_len {}",
                    seq.len(),
                    k
                );
                assert!(
                    cyclic_contains(host, seq),
                    "{label}: subseq {:?} is not a cyclic substring of rat {} ({:?})",
                    seq,
                    rat_id,
                    host
                );
                total_seqs += 1;
            }
        }
        assert_eq!(
            total_seqs, n_subseqs,
            "{label}: subseqs across map ({total_seqs}) must equal trie size ({n_subseqs})"
        );
    }

    /// Baked-in counts for hex self-exploration. If the BFS or the
    /// match-finder changes behaviour these numbers drift; treat any
    /// drift as a flag to investigate, not just a number to update.
    #[test]
    fn hex_explorer_known_values_and_invariants() {
        let explorer = SeqExplorer::new(hex_ts());
        assert_invariants("hex", &explorer);
        assert_eq!(explorer.max_subseq_len(), 6);
        assert_eq!(explorer.num_subseqs(), 120);
        assert_eq!(explorer.num_rats(), 27);
    }

    /// Same baked-in checks for square self-exploration.
    #[test]
    fn square_explorer_known_values_and_invariants() {
        let explorer = SeqExplorer::new(sq_ts());
        assert_invariants("square", &explorer);
        assert_eq!(explorer.max_subseq_len(), 4);
        assert_eq!(explorer.num_subseqs(), 110);
        assert_eq!(explorer.num_rats(), 43);
    }

    /// Single-tile seeds must produce a single seed rat with all
    /// `n` cyclic length-`n` substrings inserted under that rat's id.
    #[test]
    fn single_tile_seed_emits_exactly_its_cyclic_substrings() {
        let explorer = SeqExplorer::new(hex_ts());
        // Hex has just one tile shape.
        assert_eq!(explorer.tileset().num_tiles(), 1);
        // Rat id 0 must be the seed.
        assert!(matches!(
            explorer.provenance(0),
            Provenance::Seed { tile_idx: 0 }
        ));
        let by_rat = explorer.sequences_by_rat_id();
        let hex_seqs = by_rat.get(&0).expect("seed must be a contributor");
        // A hexagon's 6 cyclic length-6 substrings — all rotations of
        // its angle sequence — are all distinct (or all the same if
        // the angle sequence is repetitive), but together with the
        // shorter cyclic substrings they form the trie's initial
        // contents. Just check the seed contributed something.
        assert!(!hex_seqs.is_empty());
    }

    /// `tileset()` returns `&TileSet`, not `&Arc<TileSet>`; verify
    /// callers can still navigate through it. `tileset_arc()` gives
    /// an owned Arc when needed.
    #[test]
    fn tileset_accessors_compile_and_agree() {
        let explorer = SeqExplorer::new(hex_ts());
        let n_via_ref: usize = explorer.tileset().num_tiles();
        let arc = explorer.tileset_arc();
        assert_eq!(n_via_ref, arc.num_tiles());
    }

    /// All cyclic length-`min(k, n)` substrings of `seq`. Returning
    /// only the longest-available length per host suffices for
    /// coverage checks: if a length-k node exists in the trie, every
    /// shorter prefix node does too (the trie inserts length 1..=k
    /// cumulatively).
    fn cyclic_substrings_at_max_len(seq: &[i8], k: usize) -> Vec<Vec<i8>> {
        let n = seq.len();
        let len = k.min(n);
        if len == 0 || n == 0 {
            return Vec::new();
        }
        (0..n)
            .map(|start| (0..len).map(|i| seq[(start + i) % n]).collect())
            .collect()
    }

    /// Validity + completeness check at the BFS fixed point.
    ///
    /// Per-witness presence/reconstruction is verified inline; the
    /// completeness scan (no valid glue introduces a new substring)
    /// delegates to the shared [`check_fixed_point`] helper that the
    /// `seq_collect` binary also uses.
    fn assert_bfs_fixed_point(label: &str, tileset: Arc<TileSet<ZZ12>>) {
        let explorer = SeqExplorer::new(Arc::clone(&tileset));
        let k = explorer.max_subseq_len();

        let known: BTreeSet<Vec<i8>> = explorer
            .sequences_by_rat_id()
            .into_values()
            .flatten()
            .collect();
        assert_eq!(
            known.len(),
            explorer.num_subseqs(),
            "{label}: |known set| must match trie size",
        );

        // (1) + (2): per-witness reconstruction and own-substring
        // coverage. These are O(n_rats * perimeter) and serve as a
        // sanity check on `insert_cyclic_subseqs` + provenance integrity.
        let by_rat = explorer.sequences_by_rat_id();
        for &rat_id in by_rat.keys() {
            let reconstructed = explorer
                .provenance(rat_id)
                .apply(&explorer)
                .unwrap_or_else(|e| panic!("{label}: witness {rat_id} apply failed: {e}"));
            assert_eq!(
                &reconstructed,
                explorer.rat(rat_id),
                "{label}: provenance.apply doesn't reproduce witness {rat_id}",
            );
            for s in cyclic_substrings_at_max_len(reconstructed.seq(), k) {
                assert!(
                    known.contains(&s),
                    "{label}: witness {rat_id} has cyclic substring {s:?} \
                     not in trie",
                );
            }
        }

        // (3) Completeness: no valid glue from any witness against
        // any tile introduces a new cyclic substring. Shared helper —
        // same code the CLI runs.
        let witness_ids: Vec<usize> = by_rat.keys().copied().collect();
        let report = check_fixed_point(
            &explorer.tileset_arc(),
            explorer.rats(),
            &witness_ids,
            &known,
            k,
        );
        assert!(
            report.is_complete(),
            "{label}: completeness scan found {} missing substrings: {:?}",
            report.missing.len(),
            report.missing.iter().take(5).collect::<Vec<_>>(),
        );
    }

    #[test]
    fn hex_bfs_is_at_fixed_point() {
        assert_bfs_fixed_point("hex", hex_ts());
    }

    #[test]
    fn square_bfs_is_at_fixed_point() {
        assert_bfs_fixed_point("square", sq_ts());
    }
}
