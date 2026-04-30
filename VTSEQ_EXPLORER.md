# PatchSeqExplorer

Seq explorer reimplemented with GrowingPatch witnesses instead of Rats.

Same algorithm as `seq_explorer` (enumerate all cyclic angle subsequences of
length up to k), but uses GrowingPatch internally for precomputed candidates
and edge provenance tracking. Dedup by Rat (angle sequence only).

## What changed from seq_explorer

- **Witness type**: `GrowingPatch<T>` instead of `Rat<T>`
- **Dedup**: `BTreeMap<Rat<T>, usize>` — one GrowingPatch per unique boundary Rat
- **Extension**: `gp.get_all_matches()` + `gp.add_tile()` instead of `MatchFinder::crossing()` + `mt.apply()`
- **Everything else identical**: trie structure, active mask, BFS layers, provenance

## Data Structures

```
SeqTrie (i8 keys)               — identical to seq_explorer
PatchSeqExplorer<T> {
    witnesses: Vec<GrowingPatch<T>>,
    rat_to_id: BTreeMap<Rat<T>, usize>,
    provenances: Vec<Provenance>,
    witness_indices: Vec<usize>,
    active_map: HashMap<usize, Vec<bool>>,
}
```

## Algorithm

1. **Seed**: Insert single-tile angle sequences. Create 2-tile GrowingPatches via `add_tile()`, insert subsequences, track active masks.
2. **BFS layers**: For each pending GrowingPatch:
   - Get matches from precomputed `candidates_by_start`
   - Filter by active mask
   - `add_tile()` → new GrowingPatch
   - Dedup via `to_rat()` → `rat_to_id`
   - Insert cyclic subsequences, track active masks
   - Patches with new subseqs → next pending set
3. **Terminate** when no pending patches remain.

## Files

- `src/intgeom/vtseq_explorer.rs` — core module
- `src/bin/vtseq_collect.rs` — CLI binary (planned)
