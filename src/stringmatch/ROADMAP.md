# Stringmatch Module Roadmap

Efficient maximal match enumeration over integer sequences using suffix array data structures.

## Architecture

```
Layer 2 (Domain): cyclic doubling, reverse complement, length capping
Layer 1 (Core):   src/stringmatch/ — linear string matching
```

The core `MatchIndex` is agnostic to domain semantics. Domain concerns (cyclicity, reverse complement) are handled by transforming inputs before indexing and post-processing query results.

## Phases

### Phase 1: Suffix Array Construction (`sais.rs`) — **CURRENT**

- [x] SA-IS algorithm for integer sequences (small alphabet, u8 after remapping)
- [x] Comprehensive tests: empty, single element, repeated chars, known examples, random vs naive

### Phase 2: LCP Array (`lcp.rs`)

- [ ] Kasai's algorithm: O(n) LCP array construction from SA + original string
- [ ] Tests: known examples, verify LCP[rank[i]] = longest common prefix of suffix i with its predecessor in SA order

### Phase 3: Range Minimum Query (`rmq.rs`)

- [ ] Sparse table: O(n log n) build, O(1) query for minimum in any LCP subarray
- [ ] Enables O(1) longest common prefix between arbitrary suffix pairs
- [ ] Tests: correctness against brute-force, boundary queries

### Phase 4: MatchIndex Builder (`mod.rs`)

- [ ] Concatenate input strings with unique sentinels
- [ ] Build SA, LCP, document array (string ID + position), sparse table RMQ
- [ ] Tests: verify document array, concatenation integrity

### Phase 5: Maximal Match Enumeration (`maximal.rs`)

- [ ] LCP interval traversal (stack-based bottom-up suffix tree walk)
- [ ] For pair (i, j): enumerate all maximal common substrings
  - Right-maximality: automatic by construction (different children in LCP interval)
  - Left-maximality: check preceding characters differ or at boundary
- [ ] Tests: hand-computed examples, verify against brute-force for small inputs

### Phase 6: Integration & Polish

- [ ] End-to-end tests with realistic tile-matching inputs
- [ ] Benchmark against brute-force to validate performance
- [ ] Doc examples, final API ergonomics

## Key Invariants

- All algorithms target linear or near-linear time
- Total concatenated string size is assumed small (<100K), so constant factors are secondary to clarity
- Alphabet is i8, remapped to u8 — SA-IS exploits the small bounded alphabet
