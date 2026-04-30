# SeqExplorer Speed Analysis

## Baseline (release mode)

| Tileset | SeqExplorer | PatchSeqExplorer | Notes |
|---|---|---|---|
| hex | 73ms | 467ms | 120 subseqs |
| square | 92ms | 264ms | 110 subseqs |
| mixed | 795s | OOM at layer 6 | 29,066 subseqs |

## Mixed Layer-by-Layer (SeqExplorer baseline)

| Layer | Batch | New Rats | New Wit | Trie | Pending | Skip/Total |
|---|---|---|---|---|---|---|
| 1 | 2 | 3 | 3 | 79 | 3 | 0/100 |
| 2 | 3 | 18 | 18 | 439 | 18 | 0/228 |
| 3 | 18 | 141 | 138 | 1,724 | 138 | 0/1,750 |
| 4 | 138 | 1,209 | 779 | 5,009 | 779 | 448/15,580 |
| 5 | 779 | 9,008 | 2,617 | 11,215 | 2,617 | 19,824/99,804 |
| 6 | 2,617 | 33,659 | 5,000 | 19,312 | 5,000 | 137,524/372,084 |
| 7 | 5,000 | 64,903 | 5,252 | 25,859 | 5,252 | 373,978/762,540 |
| 8 | 5,252 | 63,867 | 2,555 | 28,596 | 2,555 | 485,790/830,570 |
| 9 | 2,555 | 28,839 | 458 | 29,054 | 458 | 268,212/414,268 |
| 10 | 458 | 4,703 | 12 | 29,066 | 12 | 52,330/75,374 |
| 11 | 12 | 102 | 0 | 29,066 | 0 | 1,448/1,952 |
| **Total** | | 206,454 | 16,834 | 29,066 | | 1,339,554/2,574,250 |

**52% skip rate** — active mask filter prunes about half the candidates.

## Bottleneck Breakdown

### Per-layer costs (dominant at layers 6-8):

1. **CyclicMatchIndex construction** (suffix array over batch + tileset): O(S log S) where S = sum of doubled+RC sequences
2. **Match iteration + glue** (~400K candidates per layer): each `try_glue_precomputed` costs O(n) for concatenation + O(n) for Booth's normalization
3. **Rat dedup** (`BTreeMap` lookup): O(n × log N) per check, where N = total known rats
4. **Trie insertion** (`insert_cyclic_subseqs`): O(n × k) per new rat, called for ALL new rats (even non-witnesses)
5. **seen_new dedup** (`BTreeSet`): O(n × log M) per check, where M = new rats in current layer

### Cost profiling estimate for layer 7 (batch=5000):
- ~762K total candidates, ~388K pass active filter
- ~65K glue operations at O(n=30) each = ~2M operations
- ~65K BTreeMap lookups at O(30 × 18) each = ~35M comparisons
- ~65K BTreeSet inserts at O(30 × 16) each = ~31M comparisons
- ~65K trie insertions at O(30 × 6) each = ~12M operations

**Dedup dominates** — BTreeMap + BTreeSet together account for ~66M comparisons vs ~14M for glue+trie.

## Optimization Opportunities

### A. HashMap Dedup (IMPLEMENTED)

**Problem**: `rat_to_id: BTreeMap<Rat<T>, usize>` and `seen_new: BTreeSet<Rat<T>>` do O(n) comparison per step × O(log N) steps. For 200K rats with n≈30, that's ~18 × 30 = 540 comparisons per lookup.

**Fix**: Replace both with `HashMap`/`HashSet` using Rat's existing `Hash` impl (hashes the doubled angle sequence).

**Estimated impact**: 5-10x speedup on dedup, which dominates late layers.

### B. Pre-Dedup Before Glue (IMPLEMENTED)

**Problem**: `valid_matches()` returns all matches for a (batch_idx, tile_idx) pair. Multiple matches can produce the same glued rat (same boundary, different match location). Each calls `try_glue_precomputed` (O(n)) even if the result is a duplicate.

**Fix**: Within a (batch_idx, tile_idx) pair, collect matches into a HashSet keyed by a cheap fingerprint (start_a, len) before computing glue. This avoids redundant glue calls for matches that would produce the same boundary.

**Estimated impact**: Avoids ~10-30% of glue calls that produce duplicates within the same pair.

### C. Pruned Trie Insertion (IMPLEMENTED)

**Problem**: `insert_cyclic_subseqs` is called for ALL new rats (including non-witnesses that contribute no new subseqs). At layer 7, 65K rats get full trie insertion even though only ~5K become witnesses.

Actually, trie insertion already returns early if no new subseqs are found. The real issue is that non-witness rats still traverse the trie (O(n × k) each) even though they contribute nothing new.

**Fix**: Only insert subseqs at positions near active junctions. Positions far from junctions can't contribute new subseqs that weren't already found by shorter rats. Use the junction positions from the provenance to limit trie insertion to a window around junctions.

**Estimated impact**: For n=30 rats with 2 junctions and k=6, insertion drops from ~180 operations to ~36 per rat (only 2 windows of k positions). ~5x reduction for non-witness rats.

### D. Rayon Parallelism (NOT YET IMPLEMENTED)

**Fix**: Parallelize the match iteration loop with `rayon`.

**Estimated impact**: 4-8x on layers.

### E. Try-Rollback for PatchSeqExplorer (NOT YET IMPLEMENTED)

**Fix**: Implement `add_tile_in_place` with rollback, eliminating clone_for_mutation for failed attempts.

### F. Incremental Candidate Update for PatchSeqExplorer (NOT YET IMPLEMENTED)

**Fix**: After add_tile, only recompute candidates near junctions instead of the full boundary.

**Estimated impact**: 10-25x on add_tile.

### G. Hybrid Approach (NOT YET IMPLEMENTED)

Use SeqExplorer's MatchFinder::crossing() for batch matching, construct GrowingPatches only for extensions that produce new subseqs.

## Results

| Optimization | Mixed Time | Speedup |
|---|---|---|
| Baseline | 795s | 1.0x |
| A: FxHashMap dedup (rat_to_id + seen_new) | 789s | 1.01x |
| B: Avoid redundant glue (valid_matches_with_rats) | 783s | 1.02x |
| Active-filtered CMI scan (position-aware matching) | 464s | 1.71x |
| + FxHashSet in valid_matches_filtered | 406s | 1.96x |
| + Deferred Snake validation (per-layer invalid cache) | 391s | 2.03x |

| Tileset | Baseline | Optimized | Speedup |
|---|---|---|---|
| hex | 73ms | 33ms | 2.2x |
| square | 92ms | 37ms | 2.5x |
| mixed | 795s | 391s | 2.0x |

## Profiling Results

Manual timing instrumentation (release mode, mixed tileset, before Snake deferral):
- **CMI suffix array scan**: 311s (75%)
- **Snake validation**: 100s (24%) — 304K calls, 53K fail (17%)
- **Glue + single-edge**: 2s (<1%)
- **Trie + matchfind + other**: 2s (<1%)

Post-Snake-deferral profile (estimated):
- **CMI scan**: ~311s (80%)
- **Snake validation**: ~80s (20%) — reduced from 304K to ~210K calls
- **Other**: ~2s

## Remaining Bottleneck

CMI suffix array scan at active positions accounts for ~80% of remaining time.
The scan walks SA neighbors for each active position in each batch rat.
For larger tilesets this will dominate further. Options to explore:
- Direct string matching at active positions (KMP / rolling hash) instead of full SA
- Precomputed per-tile match templates applied as filters
- Batch position scanning (coalesce active positions across batch rats)

**Conclusion**: 100% of time is in match filtering. The CMI suffix array scan was the real bottleneck, not dedup or glue.

## Key Finding

The bottleneck was `find_matches_with_filter` scanning the suffix array at ALL positions for every (batch_idx, tile_idx) pair, even though the active mask means only ~10-20% of positions can produce useful matches.

**Fix**: Added `maximal_rc_matches_at_positions(i, j, positions)` that only scans SA at positions near active mask entries. Also replaced `BTreeSet` with `FxHashSet` for the dedup inside `valid_matches_filtered`.

**Result**: Mixed went from 795s → 406s (1.96x). Hex: 73ms → 46ms. Square: 92ms → 48ms.
