# Stringmatch Module Roadmap

Efficient maximal match enumeration over integer sequences using suffix array data structures.

## Architecture

```
Layer 2 (Domain): cyclic doubling, reverse complement, length capping
Layer 1 (Core):   src/stringmatch/ — linear string matching
```

The core `MatchIndex` is agnostic to domain semantics. Domain concerns (cyclicity, reverse complement) are handled by transforming inputs before indexing and post-processing query results.

## Phases

### Phase 1: Suffix Array Construction (`sais.rs`) — DONE

- [x] SA-IS algorithm for integer sequences (small alphabet, u8 after remapping)
- [x] Comprehensive tests: empty, single element, repeated chars, known examples, random vs naive

### Phase 2: LCP Array (`lcp.rs`) — DONE

- [x] Kasai's algorithm: O(n) LCP array construction from SA + original string
- [x] Tests: known examples, verify LCP[rank[i]] = longest common prefix of suffix i with its predecessor in SA order

### Phase 3: Range Minimum Query (`rmq.rs`) — DONE

- [x] Sparse table: O(n log n) build, O(1) query for minimum in any LCP subarray
- [x] Enables O(1) longest common prefix between arbitrary suffix pairs
- [x] Tests: correctness against brute-force, boundary queries

### Phase 4: MatchIndex Builder (`mod.rs`) — DONE

- [x] Concatenate input strings with unique sentinels
- [x] Build SA, LCP, document array (string ID + position), sparse table RMQ
- [x] Tests: verify document array, concatenation integrity

### Phase 5: Maximal Match Enumeration (`maximal.rs`) — DONE

- [x] SA scan with running LCP minimum per position
- [x] For pair (i, j): enumerate all maximal common substrings
  - Right-maximality: automatic (LCP = exact match length, can't extend)
  - Left-maximality: check preceding characters differ or at boundary
- [x] Self-match handling (i == j): normalized deduplication
- [x] Tests: hand-computed examples, verify against brute-force for small inputs, random stress

### Phase 6: Integration & Polish — DONE

- [x] End-to-end tests with realistic tile-matching inputs (spectre, hexagon, tetrominos, penrose rhombs)
- [x] Match content verification (RC relationship checked for every match)
- [x] Naive brute-force cross-validation for all e2e test cases
- [x] Public API: `CyclicMatchIndex` and `CyclicMatch` re-exported from `stringmatch` module
- [x] Doc comments on all public types and methods
- [ ] Benchmark against brute-force to validate performance (optional, scale is small)

## Key Invariants

- All algorithms target linear or near-linear time
- Total concatenated string size is assumed small (<100K), so constant factors are secondary to clarity
- Alphabet is i8, remapped to u32 — SA-IS exploits the bounded alphabet
- MatchIndex is built once, queried many times for different pairs
