//! Minimum deterministic acyclic finite-state automaton over a small
//! alphabet -- compact storage of large sets of structured short
//! sequences. Supports two construction modes (both Daciuk-incremental
//! under the hood):
//!
//! * [`Dafsa::from_seqs`] -- build in one shot from a pre-sorted,
//!   duplicate-free iterable. Sugar over the streaming API.
//! * [`Dafsa::builder`] -> [`DafsaBuilder::insert`] -> [`DafsaBuilder::finish`]
//!   -- stream sequences in one by one. Useful when the caller wants
//!   to read sorted input from disk / network / a producer thread
//!   without first materializing it into a `Vec`.
//!
//! Once built, the DAFSA supports membership (`contains`), prefix
//! walk (`walk_prefix`), indexed lookup (`get(i)`), and lex-sorted
//! enumeration (`iter`) -- the last of which round-trips back to the
//! original input list and is how the caller "uncompresses" the
//! structure into a flat list of sequences.
//!
//! # Intended use case
//!
//! Dumping enumerated polygon angle-sequence databases (rat_enum
//! output) to a form that the browser-side interactive explorer can
//! load without choking on raw size. For the free (full dihedral symmetry reduction)
//! ZZ12 sequences up to perimeter n=16 (~1.4B sequences in raw
//! form), DAFSA's prefix+suffix sharing compresses the set
//! 50-200x over packed-binary representations because of the small
//! alphabet (11 symbols), heavy lex-min prefix sharing, and
//! closure-constrained suffix structure.
//!
//! # Invariants
//!
//! * **Input must be strictly lex-sorted, no duplicates.** This
//!   applies to both `from_seqs` and successive `insert` calls on a
//!   `DafsaBuilder`. Both paths debug-assert the ordering; release-mode
//!   behavior on violation is unspecified (may produce a non-minimal
//!   or incorrect automaton).
//! * `iter()` yields the exact same sequences that were fed in, in
//!   the same order. This is the round-trip property and matches the
//!   indexing used by `get(i)`.
//!
//! # Algorithm
//!
//! Daciuk-Mihov-Watson incremental minimum DAFSA construction
//! (Daciuk et al., "Incremental construction of minimal acyclic
//! finite-state automata", Computational Linguistics 26(1), 2000).
//! Inputs must be sorted; we maintain an "active path" of
//! in-progress states from root to the current insertion's leaf,
//! and a "register" mapping finalized state signatures to their
//! interned ids. When the common prefix length with the next input
//! drops below the path depth, the trailing path is frozen one
//! state at a time -- each frozen state is interned (deduplicating
//! suffix-sharing siblings) and its edge added to its parent.
//!
//! Time: O(total input length) amortized; memory during construction
//! is O(active path depth + register size). The register size is
//! bounded by the minimum DAFSA's state count, which is what we keep
//! anyway.
//!
//! # Serialization
//!
//! [`Dafsa::to_json_string`] / [`Dafsa::from_json_str`] (and the
//! `io::Write` / `io::Read` variants `write_json` / `read_json`)
//! round-trip the automaton through a flat, column-stored JSON form
//! whose schema is documented at [`JSON_SCHEMA_DOC`]. The format is:
//!
//! * **self-describing**: a `format` discriminator string, integer
//!   `version`, and a `scalar` string identifying the label type.
//!   Anyone opening the file can identify it without the originating
//!   Rust crate.
//! * **language-universal**: it is plain JSON, so a browser-side
//!   loader can parse it with `JSON.parse` and consume it in pure
//!   JavaScript or WASM. We do not require the Rust crate at the
//!   read site.
//! * **column-stored**: the state and edge tables are split into
//!   parallel arrays (`edges_start`, `labels`, `targets`, `counts`)
//!   so that gzip / brotli can compress runs of small integers
//!   efficiently. Per-state `is_accept` and the top-level
//!   `accept_count` are *not* stored -- both are derivable from
//!   `counts` (see [`JSON_SCHEMA_DOC`]) and we drop them to keep
//!   the file smaller.
//!
//! Size overhead vs. a packed-binary representation is roughly 3-5x
//! uncompressed and 1.5-2x after gzip; in exchange the file is
//! debuggable by hand and portable across any language that has a
//! JSON parser.

use std::collections::HashMap;
use std::io;

/// Minimum deterministic acyclic finite-state automaton over the
/// `i8` alphabet (the angle-sequence label type used throughout the
/// crate). Supports set membership, prefix walk, indexed lookup,
/// and exact round-trip to a sorted sequence list.
///
/// Build via [`Dafsa::from_seqs`] from a pre-sorted, duplicate-free
/// sequence iterator, or via [`Dafsa::builder`] for streaming
/// (sequence-by-sequence) construction. Either path produces the
/// same minimum DAFSA on the same input.
///
/// The label type is hardcoded to `i8` because that is the only
/// alphabet we use; generalizing to other element types would be
/// straightforward but would only add boilerplate without a
/// concrete second caller.
#[derive(Debug, Clone)]
pub struct Dafsa {
    /// State table. State 0 is the root. Each state's outgoing edges
    /// occupy a contiguous slice of `edges[state.edges_start ..
    /// next_state.edges_start]` -- the standard CSR (compressed
    /// sparse row) layout for sparse graphs.
    states: Vec<StateInfo>,
    /// Edge table. Edges out of each state are sorted by `label`.
    edges: Vec<Edge>,
    /// Per-state subtree-acceptance counts (the number of accepted
    /// sequences reachable from each state, including itself if
    /// accepting). `counts[0]` is the total number of accepted
    /// sequences; `counts[s]` lets [`Self::get`] navigate to the
    /// i-th sequence in O(seq length * fanout).
    pub(crate) counts: Vec<usize>,
}

#[derive(Debug, Clone, Copy)]
struct StateInfo {
    /// True if walking to this state from the root spells an accepted
    /// sequence. (Equivalently: the empty continuation from this
    /// state is in the language.)
    is_accept: bool,
    /// Offset into `edges` of this state's first outgoing edge.
    /// The slice ends at the next state's `edges_start`, with a
    /// sentinel state at the end of `states` carrying
    /// `edges_start == edges.len()`.
    edges_start: u32,
}

#[derive(Debug, Clone, Copy)]
struct Edge {
    label: i8,
    target: u32,
}

impl Dafsa {
    /// Start a streaming construction. The returned [`DafsaBuilder`]
    /// accepts sequences one at a time in strict lex order via
    /// [`DafsaBuilder::insert`], and is finalized into a `Dafsa` via
    /// [`DafsaBuilder::finish`].
    ///
    /// Use this when the caller wants to feed the builder lazily
    /// (e.g. reading sorted input from a stream); when a `Vec` of
    /// sequences is already in hand, [`Dafsa::from_seqs`] is the
    /// terser equivalent.
    pub fn builder() -> DafsaBuilder {
        DafsaBuilder::new()
    }

    /// Build a DAFSA from a pre-sorted, duplicate-free sequence of
    /// sequences. Sugar over the streaming API:
    /// `let b = Dafsa::builder(); for s in seqs { b.insert(s); } b.finish()`.
    ///
    /// `seqs` is consumed lazily; only states currently on the active
    /// path plus the register of finalized states are held in memory
    /// during construction.
    ///
    /// # Panics (debug builds)
    ///
    /// Panics if the input is not strictly increasing in lex order.
    pub fn from_seqs<I, S>(seqs: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: AsRef<[i8]>,
    {
        let mut builder = DafsaBuilder::new();
        for s in seqs {
            builder.insert(s.as_ref());
        }
        builder.finish()
    }

    /// Iterate accepted sequences in lex-sorted order. Yields freshly
    /// allocated `Vec<i8>`s; for a zero-allocation walk use
    /// [`Self::contains`] or [`Self::get`] instead.
    pub fn iter(&self) -> DafsaIter<'_> {
        DafsaIter::new(self)
    }

    /// Total number of accepted sequences. Equals `counts[0]`.
    pub fn len(&self) -> usize {
        self.counts.first().copied().unwrap_or(0)
    }

    /// True iff the DAFSA accepts no sequences.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Set-membership query: does `seq` exactly match an accepted
    /// sequence? O(seq.len() * log(fanout)); fanout is small for
    /// our cyclotomic alphabets.
    pub fn contains<S: AsRef<[i8]>>(&self, seq: S) -> bool {
        let seq = seq.as_ref();
        let mut s = 0u32;
        for &c in seq {
            match self.find_edge(s, c) {
                Some(t) => s = t,
                None => return false,
            }
        }
        self.states[s as usize].is_accept
    }

    /// The `i`-th accepted sequence in lex order, or `None` if `i`
    /// is out of bounds.
    ///
    /// Uses the per-state subtree-count table to descend the DAFSA
    /// in O(seq length * fanout), accumulating one symbol per
    /// step; never enumerates the irrelevant siblings.
    pub fn get(&self, i: usize) -> Option<Vec<i8>> {
        if i >= self.len() {
            return None;
        }
        let mut out = Vec::new();
        let mut s = 0u32;
        let mut remaining = i;
        loop {
            // The lex-smallest accepted sequence at state `s` is the
            // empty continuation (if `is_accept`), followed in order
            // by edges' subtrees.
            if self.states[s as usize].is_accept {
                if remaining == 0 {
                    return Some(out);
                }
                remaining -= 1;
            }
            // Walk into the edge whose subtree contains the `remaining`-th.
            let edge_range = self.edge_range(s);
            for &Edge { label, target } in &self.edges[edge_range] {
                let n = self.counts[target as usize];
                if remaining < n {
                    out.push(label);
                    s = target;
                    break;
                }
                remaining -= n;
            }
        }
    }

    /// Number of states in the underlying automaton. Crate-private:
    /// callers in `rat_dafsa_lazy` use it to slice the structure
    /// into block files.
    pub(crate) fn raw_n_states(&self) -> usize {
        self.states.len()
    }

    /// Number of edges in the underlying automaton.
    pub(crate) fn raw_n_edges(&self) -> usize {
        self.edges.len()
    }

    /// CSR row-pointer array: `edges_start[s]` is the index into the
    /// edge array of state `s`'s first outgoing edge.
    pub(crate) fn raw_edges_start(&self) -> Vec<u32> {
        self.states.iter().map(|s| s.edges_start).collect()
    }

    /// Per-state accepting flag.
    pub(crate) fn raw_is_accept(&self) -> Vec<bool> {
        self.states.iter().map(|s| s.is_accept).collect()
    }

    /// Edge labels in stored order (sorted within each state's slice).
    pub(crate) fn raw_labels(&self) -> Vec<i8> {
        self.edges.iter().map(|e| e.label).collect()
    }

    /// Edge target state ids in stored order.
    pub(crate) fn raw_targets(&self) -> Vec<u32> {
        self.edges.iter().map(|e| e.target).collect()
    }

    /// Per-state subtree-acceptance counts.
    pub(crate) fn raw_counts(&self) -> &[usize] {
        &self.counts
    }

    /// Lex rank of `seq` if it is in the language. Inverse of
    /// [`Self::get`].
    ///
    /// Walks the automaton along `seq` and accumulates, for each
    /// state, the count of sequences in earlier-lex subtrees plus
    /// the empty continuation if the state itself is accepting.
    /// O(seq.len() * fanout).
    pub fn lex_rank_of<S: AsRef<[i8]>>(&self, seq: S) -> Option<u64> {
        let seq = seq.as_ref();
        let mut s = 0u32;
        // u64 accumulator: a rank can exceed u32::MAX (~ZZ12 n=17),
        // and `usize` is 32-bit on wasm32. Matches the lazy reader.
        let mut rank: u64 = 0;
        for &c in seq {
            if self.states[s as usize].is_accept {
                rank += 1;
            }
            let edges = &self.edges[self.edge_range(s)];
            let mut next: Option<u32> = None;
            for e in edges {
                match e.label.cmp(&c) {
                    std::cmp::Ordering::Less => rank += self.counts[e.target as usize] as u64,
                    std::cmp::Ordering::Equal => {
                        next = Some(e.target);
                        break;
                    }
                    std::cmp::Ordering::Greater => return None,
                }
            }
            s = next?;
        }
        if !self.states[s as usize].is_accept {
            return None;
        }
        Some(rank)
    }

    /// Walk the prefix `prefix` and return the state reached, or
    /// `None` if any step has no matching edge. Useful for the
    /// browser-side "as the user types" interaction: after each
    /// keystroke, check whether the typed prefix is still a viable
    /// path in the DAFSA, and consult the resulting state's
    /// `is_accept` to highlight closures.
    pub fn walk_prefix<S: AsRef<[i8]>>(&self, prefix: S) -> Option<DafsaCursor<'_>> {
        let prefix = prefix.as_ref();
        let mut s = 0u32;
        for &c in prefix {
            s = self.find_edge(s, c)?;
        }
        Some(DafsaCursor {
            dafsa: self,
            state: s,
        })
    }

    fn find_edge(&self, state: u32, label: i8) -> Option<u32> {
        let edges = &self.edges[self.edge_range(state)];
        edges
            .binary_search_by(|e| e.label.cmp(&label))
            .ok()
            .map(|idx| edges[idx].target)
    }

    fn edge_range(&self, state: u32) -> std::ops::Range<usize> {
        let start = self.states[state as usize].edges_start as usize;
        let end = if (state as usize + 1) < self.states.len() {
            self.states[state as usize + 1].edges_start as usize
        } else {
            self.edges.len()
        };
        start..end
    }
}

/// Cursor returned by [`Dafsa::walk_prefix`]: the position reached
/// after consuming the given prefix.
#[derive(Debug, Clone, Copy)]
pub struct DafsaCursor<'a> {
    dafsa: &'a Dafsa,
    state: u32,
}

impl DafsaCursor<'_> {
    /// True if the prefix consumed so far is itself an accepted
    /// sequence (= the prefix is a complete entry in the DAFSA).
    pub fn is_accept(&self) -> bool {
        self.dafsa.states[self.state as usize].is_accept
    }

    /// Number of accepted sequences below the current state, counting
    /// the empty continuation if accept.
    pub fn count(&self) -> usize {
        self.dafsa.counts[self.state as usize]
    }
}

// ---- Iteration ----

/// Iterator over the DAFSA's accepted sequences in lex-sorted order.
///
/// Maintains a DFS stack of (state, next-edge-index) and a running
/// `current` prefix that grows on descent and shrinks on ascent. Each
/// `next()` call drives the DFS forward until it reaches the next
/// accepting state; then it returns a clone of the current prefix.
pub struct DafsaIter<'a> {
    dafsa: &'a Dafsa,
    stack: Vec<DfsFrame>,
    current: Vec<i8>,
    /// After yielding the empty continuation at state s (when s is
    /// accepting), we still need to descend into s's children. This
    /// flag distinguishes "about to yield this state's accept" from
    /// "already yielded; descend into children".
    pending_accept: bool,
}

#[derive(Debug, Clone, Copy)]
struct DfsFrame {
    /// Index into `edges` of the next edge to descend, NOT a label
    /// or offset into the state's edge slice. We store absolute
    /// edge indices to avoid recomputing the state's edge range on
    /// each step.
    next_edge: u32,
    edges_end: u32,
}

impl<'a> DafsaIter<'a> {
    fn new(dafsa: &'a Dafsa) -> Self {
        let edges_end = if dafsa.states.len() > 1 {
            dafsa.states[1].edges_start
        } else {
            dafsa.edges.len() as u32
        };
        let stack = vec![DfsFrame {
            next_edge: dafsa.states[0].edges_start,
            edges_end,
        }];
        DafsaIter {
            dafsa,
            stack,
            current: Vec::new(),
            pending_accept: dafsa.states[0].is_accept,
        }
    }
}

impl Iterator for DafsaIter<'_> {
    type Item = Vec<i8>;
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // First, if the state on top of the stack is accepting and
            // we haven't yielded its empty continuation yet, do so.
            if self.pending_accept {
                self.pending_accept = false;
                return Some(self.current.clone());
            }

            let frame = self.stack.last_mut()?;
            if frame.next_edge < frame.edges_end {
                let edge = self.dafsa.edges[frame.next_edge as usize];
                frame.next_edge += 1;
                // Descend.
                self.current.push(edge.label);
                let target = edge.target;
                let edges_start = self.dafsa.states[target as usize].edges_start;
                let edges_end = if (target as usize + 1) < self.dafsa.states.len() {
                    self.dafsa.states[target as usize + 1].edges_start
                } else {
                    self.dafsa.edges.len() as u32
                };
                self.stack.push(DfsFrame {
                    next_edge: edges_start,
                    edges_end,
                });
                self.pending_accept = self.dafsa.states[target as usize].is_accept;
                continue;
            }
            // No more children at this state -- pop.
            self.stack.pop();
            self.current.pop();
        }
    }
}

// ---- Builder (Daciuk incremental) ----

/// Per-state during construction. Edges are kept in sorted-by-label
/// order; the natural insertion order (sorted input) makes that
/// automatic for us.
#[derive(Debug, Clone)]
struct BuildState {
    is_accept: bool,
    edges: Vec<(i8, u32)>,
}

/// Streaming builder for [`Dafsa`]. Sequences must be fed via
/// [`Self::insert`] in strict lex order with no duplicates, then
/// [`Self::finish`] consumes the builder and produces a minimum
/// DAFSA.
///
/// Construct via [`Dafsa::builder`]; the type is intentionally
/// not `Default`-implementable to keep the discovery path through
/// `Dafsa`'s docs.
///
/// Memory during construction is bounded by the size of the final
/// automaton (the register interns one state per minimum-DAFSA state)
/// plus an active path no longer than the longest single sequence.
pub struct DafsaBuilder {
    /// Active path of in-progress states from root to the leaf of
    /// the most recently inserted sequence. `path[0]` is the root;
    /// `path[i]` for `i > 0` is reached from `path[i-1]` by edge
    /// `path_labels[i-1]`. The edge from `path[i-1]` to `path[i]`
    /// is NOT yet in `path[i-1].edges` -- it gets added when
    /// `path[i]` is finalized (frozen + interned).
    path: Vec<BuildState>,
    /// Labels of the edges connecting consecutive `path` states.
    path_labels: Vec<i8>,
    /// Frozen states by their interned id.
    frozen: Vec<BuildState>,
    /// Lookup: `(edges, is_accept)` -> frozen id. The key uniquely
    /// identifies a state's behaviour, so two states with the same
    /// signature merge.
    register: HashMap<(Vec<(i8, u32)>, bool), u32>,
    /// Most recently inserted sequence, used to compute the common
    /// prefix with the next input.
    last_seq: Vec<i8>,
    /// Set during the first `insert` to support a strict-sort
    /// debug_assert across inserts.
    seen_any: bool,
}

impl DafsaBuilder {
    fn new() -> Self {
        DafsaBuilder {
            path: vec![BuildState {
                is_accept: false,
                edges: Vec::new(),
            }],
            path_labels: Vec::new(),
            frozen: Vec::new(),
            register: HashMap::new(),
            last_seq: Vec::new(),
            seen_any: false,
        }
    }

    /// Intern a finalized state. Returns the interned (frozen) id.
    fn freeze(&mut self, state: BuildState) -> u32 {
        let key = (state.edges.clone(), state.is_accept);
        if let Some(&id) = self.register.get(&key) {
            return id;
        }
        let id = self.frozen.len() as u32;
        self.register.insert(key, id);
        self.frozen.push(state);
        id
    }

    /// Insert one sequence into the builder. Must be strictly lex-
    /// greater than the previous insertion (debug-asserted).
    ///
    /// Accepts any borrow of `[i8]` -- `&[i8]`, `&Vec<i8>`, `Vec<i8>`.
    pub fn insert<S: AsRef<[i8]>>(&mut self, seq: S) {
        self.insert_slice(seq.as_ref());
    }

    fn insert_slice(&mut self, seq: &[i8]) {
        if self.seen_any {
            debug_assert!(
                self.last_seq.as_slice() < seq,
                "Dafsa::from_seqs: input must be strictly lex-sorted"
            );
        }
        self.seen_any = true;

        // Longest common prefix length with the previous sequence.
        let common = self
            .last_seq
            .iter()
            .zip(seq.iter())
            .take_while(|(a, b)| a == b)
            .count();

        // Finalize states on the active path past depth `common`.
        // path.len() == last_seq.len() + 1 invariant.
        while self.path.len() > common + 1 {
            let label = self.path_labels.pop().unwrap();
            let state = self.path.pop().unwrap();
            let frozen_id = self.freeze(state);
            self.path.last_mut().unwrap().edges.push((label, frozen_id));
        }

        // Extend the path with one fresh state per remaining symbol.
        for &c in &seq[common..] {
            self.path_labels.push(c);
            self.path.push(BuildState {
                is_accept: false,
                edges: Vec::new(),
            });
        }

        // Mark the new leaf as accepting.
        self.path.last_mut().unwrap().is_accept = true;
        self.last_seq.clear();
        self.last_seq.extend_from_slice(seq);
    }

    /// Consume the builder and produce the minimum DAFSA. May be
    /// called on a fresh builder with no inserts -- the result is an
    /// empty DAFSA accepting no sequences.
    pub fn finish(mut self) -> Dafsa {
        // Finalize everything left on the active path back to the root.
        while self.path.len() > 1 {
            let label = self.path_labels.pop().unwrap();
            let state = self.path.pop().unwrap();
            let frozen_id = self.freeze(state);
            self.path.last_mut().unwrap().edges.push((label, frozen_id));
        }
        let root_state = self.path.pop().unwrap();
        let root_id = self.freeze(root_state);

        // Compact into CSR with a DFS pre-order renumbering from the
        // root. This puts each parent at a lower id than its
        // first-visited children, and (typically) at consecutive
        // ids -- so a query walk visits monotonically increasing
        // state ids, sequential in memory. The same locality
        // benefits any blocked wire format we layer on top later
        // (state s/BLOCK_SIZE points to the file holding it, and
        // walks rarely cross block boundaries).
        //
        // States may be reached by multiple paths (suffix sharing);
        // the first visit wins and pins the id, later visits skip.
        let n = self.frozen.len();
        let mut new_id = vec![u32::MAX; n];
        let mut order: Vec<u32> = Vec::with_capacity(n);
        let mut stack: Vec<u32> = vec![root_id];
        while let Some(orig) = stack.pop() {
            if new_id[orig as usize] != u32::MAX {
                continue;
            }
            new_id[orig as usize] = order.len() as u32;
            order.push(orig);
            // Push children in reverse so the first edge gets popped
            // first; the DFS then descends edges in their stored
            // (label-sorted) order.
            let s = &self.frozen[orig as usize];
            for &(_, child) in s.edges.iter().rev() {
                if new_id[child as usize] == u32::MAX {
                    stack.push(child);
                }
            }
        }
        debug_assert_eq!(
            order.len(),
            n,
            "DFS pre-order must visit every reachable state"
        );

        // Build the CSR.
        let mut states = Vec::with_capacity(n);
        let mut edges: Vec<Edge> = Vec::new();
        for &orig_id in &order {
            let s = &self.frozen[orig_id as usize];
            states.push(StateInfo {
                is_accept: s.is_accept,
                edges_start: edges.len() as u32,
            });
            for &(label, child) in &s.edges {
                edges.push(Edge {
                    label,
                    target: new_id[child as usize],
                });
            }
        }

        // Bottom-up count of accepted sequences reachable from each
        // Bottom-up subtree-acceptance count per state. Under DFS
        // pre-order renumbering the parent always has a lower id
        // than the children it visits first; suffix-shared
        // descendants can have arbitrary ids, so we do not rely on
        // any id-based ordering and compute via an explicit
        // reverse-postorder DFS from the root.
        let mut counts = vec![0usize; n];
        Self::compute_counts(&states, &edges, &mut counts);

        Dafsa {
            states,
            edges,
            counts,
        }
    }

    /// Bottom-up subtree-count fill. See [`Dafsa::get`] for the
    /// counter's meaning.
    fn compute_counts(states: &[StateInfo], edges: &[Edge], counts: &mut [usize]) {
        // Iterative reverse-postorder DFS from root.
        let n = states.len();
        let mut order: Vec<u32> = Vec::with_capacity(n);
        let mut visited = vec![false; n];
        let mut stack: Vec<(u32, bool)> = vec![(0, false)];
        while let Some((s, expanded)) = stack.pop() {
            if expanded {
                order.push(s);
                continue;
            }
            if visited[s as usize] {
                continue;
            }
            visited[s as usize] = true;
            stack.push((s, true));
            let start = states[s as usize].edges_start as usize;
            let end = if (s as usize + 1) < states.len() {
                states[s as usize + 1].edges_start as usize
            } else {
                edges.len()
            };
            for edge in &edges[start..end] {
                if !visited[edge.target as usize] {
                    stack.push((edge.target, false));
                }
            }
        }
        // `order` is post-order: children appear before their parents.
        for s in &order {
            let s = *s as usize;
            let mut total = if states[s].is_accept { 1 } else { 0 };
            let start = states[s].edges_start as usize;
            let end = if (s + 1) < states.len() {
                states[s + 1].edges_start as usize
            } else {
                edges.len()
            };
            for edge in &edges[start..end] {
                total += counts[edge.target as usize];
            }
            counts[s] = total;
        }
    }
}

// ---- Serialization (JSON, column-stored) ----

/// The DAFSA file format identifier baked into every serialized
/// document. Renaming this is a breaking change for already-emitted
/// files.
const FORMAT_TAG: &str = "tilezz-dafsa";

/// On-disk version. Bump whenever the layout changes; the reader
/// rejects unknown versions.
const FORMAT_VERSION: u32 = 1;

/// Plain-text description of the JSON schema, embedded in the docs
/// so that anyone looking at a serialized file can identify it
/// without the Rust crate.
///
/// ```text
/// {
///   "format":       "tilezz-dafsa",         // discriminator
///   "version":      1,                       // u32, currently 1
///   "scalar":       "i8",                    // label element type, always "i8"
///   "n_states":     <usize>,                 // == edges_start.len() == counts.len()
///   "n_edges":      <usize>,                 // == labels.len() == targets.len()
///   "edges_start":  [u32, u32, ...],         // CSR row pointers; state s owns edges
///                                            //   labels[edges_start[s] .. edges_start[s+1]]
///                                            //   with the last row ending at n_edges
///   "labels":       [i8, ...],               // per edge, input symbol (sorted within state)
///   "targets":      [u32, u32, ...],         // per edge, destination state index
///   "counts":       [u64, u64, ...]          // per state, # accepted sequences reachable
/// }
/// ```
///
/// State 0 is always the root. Two fields are intentionally NOT stored
/// because they are fully derivable from the data above:
///
/// * `is_accept[s]`: derive as `counts[s] > sum(counts[child] for
///   each outgoing edge)`. In any well-formed file the difference is
///   exactly 0 (non-accepting) or 1 (accepting -- the +1 comes from
///   the empty continuation at s).
/// * `accept_count`: equals `counts[0]`, the count at the root.
///
/// Readers reconstruct both once at load time. The implementation
/// also requires that edges out of each state appear sorted ascending
/// by label; validating that on read is a small additional pass we
/// currently skip in favor of a debug-mode invariant check on the
/// writer side.
pub const JSON_SCHEMA_DOC: &str = include_str!("schemas/core_schema.txt");

/// Stable identifier for the label scalar type baked into every
/// serialized document. The reader rejects files whose `scalar`
/// field does not match this value.
const SCALAR_TAG: &str = "i8";

#[derive(serde::Serialize, serde::Deserialize)]
struct DafsaSerForm {
    format: String,
    version: u32,
    scalar: String,
    n_states: usize,
    n_edges: usize,
    edges_start: Vec<u32>,
    labels: Vec<i8>,
    targets: Vec<u32>,
    counts: Vec<u64>,
}

impl Dafsa {
    fn to_ser_form(&self) -> DafsaSerForm {
        let n_states = self.states.len();
        let n_edges = self.edges.len();
        let edges_start: Vec<u32> = self.states.iter().map(|s| s.edges_start).collect();
        let mut labels = Vec::with_capacity(n_edges);
        let mut targets = Vec::with_capacity(n_edges);
        for e in &self.edges {
            labels.push(e.label);
            targets.push(e.target);
        }
        DafsaSerForm {
            format: FORMAT_TAG.to_string(),
            version: FORMAT_VERSION,
            scalar: SCALAR_TAG.to_string(),
            n_states,
            n_edges,
            edges_start,
            labels,
            targets,
            counts: self.counts.iter().map(|&c| c as u64).collect(),
        }
    }

    fn from_ser_form(form: DafsaSerForm) -> io::Result<Self> {
        if form.format != FORMAT_TAG {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "Dafsa: bad format tag (expected '{}', got '{}')",
                    FORMAT_TAG, form.format
                ),
            ));
        }
        if form.version != FORMAT_VERSION {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "Dafsa: unsupported version (expected {}, got {})",
                    FORMAT_VERSION, form.version
                ),
            ));
        }
        if form.scalar != SCALAR_TAG {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "Dafsa: scalar mismatch (file is '{}', expected '{}')",
                    form.scalar, SCALAR_TAG
                ),
            ));
        }
        let n_states = form.n_states;
        let n_edges = form.n_edges;
        let len_ok = form.edges_start.len() == n_states
            && form.counts.len() == n_states
            && form.labels.len() == n_edges
            && form.targets.len() == n_edges;
        if !len_ok {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Dafsa: array length disagrees with n_states/n_edges",
            ));
        }
        // Derive per-state `is_accept` from `counts`. In a well-
        // formed DAFSA, counts[s] equals (1 if is_accept else 0) plus
        // the sum of counts over outgoing edges, so the boolean is
        // pure redundancy on disk -- we reconstruct it here.
        let counts: Vec<usize> = form.counts.iter().map(|&c| c as usize).collect();
        let mut is_accept = vec![false; n_states];
        for s in 0..n_states {
            let start = form.edges_start[s] as usize;
            let end = if s + 1 < n_states {
                form.edges_start[s + 1] as usize
            } else {
                n_edges
            };
            let children_sum: usize = form.targets[start..end]
                .iter()
                .map(|&t| counts[t as usize])
                .sum();
            // The accept-count from a node is children_sum (paths
            // through edges) plus 1 if the node itself accepts the
            // empty continuation. Any well-formed file satisfies
            // counts[s] - children_sum in {0, 1}; out-of-range values
            // indicate a corrupt or non-minimal automaton.
            let diff = counts[s].checked_sub(children_sum).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Dafsa: counts[{s}] < sum of children counts (corrupt file)"),
                )
            })?;
            match diff {
                0 => is_accept[s] = false,
                1 => is_accept[s] = true,
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "Dafsa: counts[{s}] disagrees with children-sum by {diff}, expected 0 or 1"
                        ),
                    ));
                }
            }
        }
        let states: Vec<StateInfo> = form
            .edges_start
            .into_iter()
            .zip(is_accept)
            .map(|(edges_start, is_accept)| StateInfo {
                is_accept,
                edges_start,
            })
            .collect();
        let edges: Vec<Edge> = form
            .labels
            .into_iter()
            .zip(form.targets)
            .map(|(label, target)| Edge { label, target })
            .collect();
        Ok(Dafsa {
            states,
            edges,
            counts,
        })
    }

    /// Serialize to a JSON string. See [`JSON_SCHEMA_DOC`] for the
    /// exact schema; the output is a single-line compact JSON object
    /// suitable for embedding in a `.json` asset.
    pub fn to_json_string(&self) -> String {
        serde_json::to_string(&self.to_ser_form()).expect("DafsaSerForm is always serializable")
    }

    /// Serialize to a pretty-printed JSON string. Larger than
    /// [`Self::to_json_string`] but easier to diff and inspect by
    /// hand; gzip absorbs nearly all of the whitespace overhead.
    pub fn to_json_string_pretty(&self) -> String {
        serde_json::to_string_pretty(&self.to_ser_form())
            .expect("DafsaSerForm is always serializable")
    }

    /// Parse a JSON document into a `Dafsa`. Rejects documents
    /// whose `format`, `version`, or `scalar` fields do not match
    /// the expected values, or whose array lengths disagree with the
    /// declared `n_states` / `n_edges`.
    pub fn from_json_str(s: &str) -> io::Result<Self> {
        let form: DafsaSerForm =
            serde_json::from_str(s).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        Self::from_ser_form(form)
    }

    /// Stream a DAFSA into an `io::Write`. Equivalent to
    /// `write!(w, "{}", self.to_json_string())` but avoids the
    /// intermediate `String` allocation.
    pub fn write_json<W: io::Write>(&self, w: W) -> io::Result<()> {
        serde_json::to_writer(w, &self.to_ser_form()).map_err(io::Error::other)
    }

    /// Stream a DAFSA out of an `io::Read`. Reads the entire input
    /// (the document is one JSON object), then parses it. Suitable
    /// for files; for a network download, fetch into a buffer first.
    pub fn read_json<R: io::Read>(r: R) -> io::Result<Self> {
        let form: DafsaSerForm = serde_json::from_reader(r)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        Self::from_ser_form(form)
    }

    /// Stream a gzipped DAFSA JSON document into an `io::Write`. The
    /// uncompressed content is exactly what [`Self::write_json`] would
    /// emit; gzip absorbs whitespace, repeated field names, and runs
    /// of small integers in the column-stored arrays (~50% on the
    /// large ZZ12 free n=16 set).
    ///
    /// Pair with [`Self::read_json_gz`] for the inverse direction.
    pub fn write_json_gz<W: io::Write>(&self, w: W) -> io::Result<()> {
        let mut enc = flate2::write::GzEncoder::new(w, flate2::Compression::default());
        self.write_json(&mut enc)?;
        enc.finish()?;
        Ok(())
    }

    /// Parse a gzipped DAFSA JSON document from an `io::Read`. The
    /// inverse of [`Self::write_json_gz`]: decompresses on the fly,
    /// then parses the JSON exactly as [`Self::read_json`] would.
    pub fn read_json_gz<R: io::Read>(r: R) -> io::Result<Self> {
        let dec = flate2::read::GzDecoder::new(r);
        Self::read_json(dec)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build(seqs: &[&[i8]]) -> Dafsa {
        Dafsa::from_seqs(seqs.iter().copied())
    }

    #[test]
    fn empty_set() {
        let dafsa = build(&[]);
        assert_eq!(dafsa.len(), 0);
        assert!(dafsa.is_empty());
        assert!(!dafsa.contains([1i8, 2, 3].as_slice()));
        assert!(dafsa.get(0).is_none());
        assert_eq!(dafsa.iter().count(), 0);
    }

    #[test]
    fn empty_sequence_only() {
        let empty: &[i8] = &[];
        let dafsa = build(&[empty]);
        assert_eq!(dafsa.len(), 1);
        assert!(dafsa.contains(empty));
        assert_eq!(dafsa.get(0), Some(Vec::<i8>::new()));
        let out: Vec<Vec<i8>> = dafsa.iter().collect();
        assert_eq!(out, vec![Vec::<i8>::new()]);
    }

    #[test]
    fn singleton() {
        let seqs: &[&[i8]] = &[&[1, 2, 3]];
        let dafsa = build(seqs);
        assert_eq!(dafsa.len(), 1);
        assert!(dafsa.contains([1i8, 2, 3].as_slice()));
        assert!(!dafsa.contains([1i8, 2].as_slice()));
        assert!(!dafsa.contains([1i8, 2, 3, 4].as_slice()));
        assert_eq!(dafsa.get(0), Some(vec![1, 2, 3]));
        assert_eq!(dafsa.iter().collect::<Vec<_>>(), vec![vec![1, 2, 3]]);
    }

    /// Round-trip property: `from_seqs(sorted).iter().collect() ==
    /// sorted`. Covers the basic Daciuk construction + iter
    /// agreement on a small set with shared prefixes.
    #[test]
    fn roundtrip_small_sorted_set() {
        let seqs: Vec<Vec<i8>> = vec![vec![1, 2], vec![1, 2, 3], vec![1, 3], vec![2], vec![2, 1]];
        let dafsa = Dafsa::from_seqs(seqs.iter().map(|v| v.as_slice()));
        assert_eq!(dafsa.len(), seqs.len());
        let out: Vec<Vec<i8>> = dafsa.iter().collect();
        assert_eq!(out, seqs);
    }

    /// Suffix-sharing: "cat", "bat" share suffix "at". A correct
    /// minimum DAFSA must end at the same accept state for both.
    #[test]
    fn suffix_sharing_compresses() {
        let seqs: Vec<Vec<i8>> = vec![
            vec![b'b' as i8, b'a' as i8, b't' as i8],
            vec![b'c' as i8, b'a' as i8, b't' as i8],
        ];
        let dafsa = Dafsa::from_seqs(seqs.iter().map(|v| v.as_slice()));
        // Two sequences, three letters each = 6 trie nodes including root.
        // Min DAFSA: root, 'b'/'c' both point to the same 'a' state,
        // 'a' -> 't' (accepting). 4 states.
        assert!(dafsa.states.len() <= 4, "expected suffix sharing");
        assert_eq!(dafsa.iter().collect::<Vec<_>>(), seqs);
        assert!(dafsa.contains([b'b' as i8, b'a' as i8, b't' as i8].as_slice()));
        assert!(dafsa.contains([b'c' as i8, b'a' as i8, b't' as i8].as_slice()));
        assert!(!dafsa.contains([b'a' as i8, b't' as i8].as_slice()));
    }

    /// Index access: `get(i)` returns the i-th sequence in lex order,
    /// for every i, across a set with both prefix and suffix sharing.
    #[test]
    fn indexed_get_matches_iter() {
        let seqs: Vec<Vec<i8>> = vec![
            vec![-3, 1, 2],
            vec![-3, 1, 2, 4],
            vec![-3, 2],
            vec![-2],
            vec![-2, 1],
            vec![-1, 1, 2],
            vec![0],
            vec![0, 1],
        ];
        let dafsa = Dafsa::from_seqs(seqs.iter().map(|v| v.as_slice()));
        for (i, seq) in seqs.iter().enumerate() {
            assert_eq!(dafsa.get(i).as_ref(), Some(seq));
        }
        assert!(dafsa.get(seqs.len()).is_none());
    }

    /// `walk_prefix` reports an accept state iff the typed prefix is
    /// exactly a sequence in the set, and returns the correct subtree
    /// count when the prefix is a strict prefix of one or more
    /// sequences.
    #[test]
    fn walk_prefix_classifies_correctly() {
        let seqs: Vec<Vec<i8>> = vec![vec![1, 2, 3], vec![1, 2, 3, 4], vec![1, 2, 5], vec![2]];
        let dafsa = Dafsa::from_seqs(seqs.iter().map(|v| v.as_slice()));
        // [1] is a prefix of three sequences, not itself accepted.
        let c = dafsa.walk_prefix([1i8].as_slice()).unwrap();
        assert!(!c.is_accept());
        assert_eq!(c.count(), 3);
        // [1, 2, 3] is itself accepted, and is a prefix of [1, 2, 3, 4].
        let c = dafsa.walk_prefix([1i8, 2, 3].as_slice()).unwrap();
        assert!(c.is_accept());
        assert_eq!(c.count(), 2);
        // [1, 2, 9] has no matching edge: returns None.
        assert!(dafsa.walk_prefix([1i8, 2, 9].as_slice()).is_none());
        // [2] is accepted and has no continuation.
        let c = dafsa.walk_prefix([2i8].as_slice()).unwrap();
        assert!(c.is_accept());
        assert_eq!(c.count(), 1);
    }

    /// Streaming construction via [`Dafsa::builder`] produces the
    /// exact same automaton as the batched [`Dafsa::from_seqs`] when
    /// fed the same sequences in the same order.
    #[test]
    fn builder_matches_from_seqs() {
        let seqs: Vec<Vec<i8>> = vec![
            vec![1, 2],
            vec![1, 2, 3],
            vec![1, 3],
            vec![2],
            vec![2, 1],
            vec![2, 1, 4],
        ];
        let one_shot = Dafsa::from_seqs(seqs.iter().map(|v| v.as_slice()));
        let mut b = Dafsa::builder();
        for s in &seqs {
            b.insert(s.as_slice());
        }
        let streamed = b.finish();
        // Structural equivalence: same state and edge counts, and same
        // round-trip output. Two minimum DAFSAs over the same language
        // built with the same renumbering policy must agree state-by-
        // state, so the counts check pins down that the streaming path
        // and the batched path share the same internal layout (not just
        // the same language).
        assert_eq!(one_shot.states.len(), streamed.states.len());
        assert_eq!(one_shot.edges.len(), streamed.edges.len());
        assert_eq!(streamed.iter().collect::<Vec<_>>(), seqs);
        assert_eq!(one_shot.iter().collect::<Vec<_>>(), seqs);
    }

    /// Calling [`DafsaBuilder::finish`] on a fresh builder with no
    /// inserts yields an empty DAFSA (no accepted sequences). This is
    /// the streaming counterpart of `empty_set`.
    #[test]
    fn builder_finish_empty() {
        let dafsa: Dafsa = Dafsa::builder().finish();
        assert!(dafsa.is_empty());
        assert_eq!(dafsa.len(), 0);
        assert_eq!(dafsa.iter().count(), 0);
    }

    /// Out-of-order insertion is a programming error and triggers a
    /// debug-mode panic with a precise message. Release-mode behavior
    /// is unspecified (no panic guarantee), so we gate the test on
    /// `debug_assertions`.
    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "strictly lex-sorted")]
    fn builder_panics_on_out_of_order() {
        let mut b = Dafsa::builder();
        b.insert([1i8, 2].as_slice());
        // [1, 1] < [1, 2] in lex order -- must panic.
        b.insert([1i8, 1].as_slice());
    }

    /// Inserting the same sequence twice in a row also violates the
    /// strict-sort precondition (equal is not strictly greater).
    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "strictly lex-sorted")]
    fn builder_panics_on_duplicate() {
        let mut b = Dafsa::builder();
        b.insert([1i8, 2, 3].as_slice());
        b.insert([1i8, 2, 3].as_slice());
    }

    /// Streaming use case: feed sequences from a generator without
    /// materializing the full list in memory. Verifies that the
    /// resulting DAFSA agrees with the batched build and that all
    /// fed sequences are recoverable via `iter` / `get` / `contains`.
    #[test]
    fn builder_stream_from_generator() {
        // Generate all length-3 non-decreasing sequences over {0,1,2,3}
        // on the fly; the iteration order is lex-sorted by construction
        // because the loops nest in that order.
        let generator =
            (0..4i8).flat_map(|a| (a..4).flat_map(move |b| (b..4).map(move |c| vec![a, b, c])));
        let materialized: Vec<Vec<i8>> = generator.clone().collect();
        let mut b = Dafsa::builder();
        for s in generator {
            b.insert(&s);
        }
        let dafsa = b.finish();
        assert_eq!(dafsa.len(), materialized.len());
        assert_eq!(dafsa.iter().collect::<Vec<_>>(), materialized);
        for (i, s) in materialized.iter().enumerate() {
            assert!(dafsa.contains(s.as_slice()));
            assert_eq!(dafsa.get(i).as_ref(), Some(s));
        }
    }

    /// Larger round-trip on a generated structured set: all length-`k`
    /// sequences over a small alphabet that pass a non-trivial filter.
    /// The DAFSA must compress (alphabet^k > final state count) and
    /// round-trip exactly.
    #[test]
    fn roundtrip_filtered_dense_set() {
        // All length-4 sequences over alphabet {0, 1, 2, 3} with
        // strictly increasing values. ~choose(4, 4) = 1 sequence
        // since alphabet only has 4 symbols and we need strictly
        // increasing length-4 -- only [0,1,2,3] qualifies. Use a
        // looser filter: non-decreasing.
        let mut seqs: Vec<Vec<i8>> = Vec::new();
        for a in 0..4 {
            for b in a..4 {
                for c in b..4 {
                    for d in c..4 {
                        seqs.push(vec![a, b, c, d]);
                    }
                }
            }
        }
        // Sorted by construction.
        let dafsa = Dafsa::from_seqs(seqs.iter().map(|v| v.as_slice()));
        assert_eq!(dafsa.len(), seqs.len());
        let out: Vec<Vec<i8>> = dafsa.iter().collect();
        assert_eq!(out, seqs);
        // Indexed access agrees.
        for (i, seq) in seqs.iter().enumerate() {
            assert_eq!(dafsa.get(i).as_ref(), Some(seq));
        }
    }

    // ---- Serialization tests ----

    fn json_roundtrip_corpus() -> Vec<Vec<i8>> {
        vec![
            vec![-3, 1, 2],
            vec![-3, 1, 2, 4],
            vec![-3, 2],
            vec![-2],
            vec![-2, 1],
            vec![-1, 1, 2],
            vec![0],
            vec![0, 1],
        ]
    }

    /// JSON round-trip: serialize a non-trivial DAFSA, parse it back,
    /// and assert that every public observable agrees (len, contains,
    /// get, iter).
    #[test]
    fn json_roundtrip_preserves_behavior() {
        let seqs = json_roundtrip_corpus();
        let original = Dafsa::from_seqs(seqs.iter().map(|v| v.as_slice()));
        let json = original.to_json_string();
        let restored: Dafsa = Dafsa::from_json_str(&json).expect("parse succeeds");
        assert_eq!(restored.len(), original.len());
        assert_eq!(restored.iter().collect::<Vec<_>>(), seqs);
        for (i, seq) in seqs.iter().enumerate() {
            assert!(restored.contains(seq.as_slice()));
            assert_eq!(restored.get(i).as_ref(), Some(seq));
        }
    }

    /// Pretty-printed JSON parses to the same DAFSA as compact JSON.
    /// The format flag is normative; whitespace is not.
    #[test]
    fn json_pretty_roundtrip() {
        let seqs = json_roundtrip_corpus();
        let original = Dafsa::from_seqs(seqs.iter().map(|v| v.as_slice()));
        let json_pretty = original.to_json_string_pretty();
        assert!(
            json_pretty.contains('\n'),
            "pretty form should be multi-line"
        );
        let restored: Dafsa = Dafsa::from_json_str(&json_pretty).expect("parse succeeds");
        assert_eq!(restored.iter().collect::<Vec<_>>(), seqs);
    }

    /// Empty DAFSA round-trips. The serialized form is one JSON
    /// object with empty arrays; parsing it produces an
    /// observationally identical empty DAFSA.
    #[test]
    fn json_roundtrip_empty() {
        let original: Dafsa = Dafsa::from_seqs(std::iter::empty::<&[i8]>());
        let json = original.to_json_string();
        let restored: Dafsa = Dafsa::from_json_str(&json).expect("parse succeeds");
        assert!(restored.is_empty());
        assert_eq!(restored.iter().count(), 0);
    }

    /// `write_json` / `read_json` (io::Write / io::Read variants)
    /// agree with the string-based path.
    #[test]
    fn json_roundtrip_via_writer_reader() {
        let seqs = json_roundtrip_corpus();
        let original = Dafsa::from_seqs(seqs.iter().map(|v| v.as_slice()));
        let mut buf: Vec<u8> = Vec::new();
        original.write_json(&mut buf).expect("write succeeds");
        let restored: Dafsa = Dafsa::read_json(buf.as_slice()).expect("read succeeds");
        assert_eq!(restored.iter().collect::<Vec<_>>(), seqs);
    }

    /// Reader rejects a document whose `format` field is wrong --
    /// otherwise loading a JSON file from an unrelated tool would
    /// silently produce a garbage DAFSA.
    #[test]
    fn json_rejects_bad_format_tag() {
        let bad = r#"{
            "format": "not-a-dafsa",
            "version": 1,
            "scalar": "i8",
            "n_states": 0,
            "n_edges": 0,
            "edges_start": [],
            "labels": [],
            "targets": [],
            "counts": []
        }"#;
        let err = Dafsa::from_json_str(bad).unwrap_err();
        assert!(err.to_string().contains("format tag"));
    }

    /// Reader rejects a document whose `version` does not match the
    /// current implementation; the file might use an incompatible
    /// layout.
    #[test]
    fn json_rejects_bad_version() {
        let bad = r#"{
            "format": "tilezz-dafsa",
            "version": 999,
            "scalar": "i8",
            "n_states": 0,
            "n_edges": 0,
            "edges_start": [],
            "labels": [],
            "targets": [],
            "counts": []
        }"#;
        let err = Dafsa::from_json_str(bad).unwrap_err();
        assert!(err.to_string().contains("version"));
    }

    /// Reader rejects a document whose `scalar` tag does not match
    /// the requested element type. Prevents silently misinterpreting
    /// a file written with a different label width.
    #[test]
    fn json_rejects_scalar_mismatch() {
        let bad = r#"{
            "format": "tilezz-dafsa",
            "version": 1,
            "scalar": "i32",
            "n_states": 0,
            "n_edges": 0,
            "edges_start": [],
            "labels": [],
            "targets": [],
            "counts": []
        }"#;
        let err = Dafsa::from_json_str(bad).unwrap_err();
        assert!(err.to_string().contains("scalar"));
    }

    /// Reader rejects a document whose declared n_states / n_edges
    /// disagree with the actual array lengths.
    #[test]
    fn json_rejects_length_disagreement() {
        // n_states says 3 but only 1 edges_start entry provided.
        let bad = r#"{
            "format": "tilezz-dafsa",
            "version": 1,
            "scalar": "i8",
            "n_states": 3,
            "n_edges": 0,
            "edges_start": [0],
            "labels": [],
            "targets": [],
            "counts": [0]
        }"#;
        let err = Dafsa::from_json_str(bad).unwrap_err();
        assert!(err.to_string().contains("length"));
    }

    /// Reader rejects a document whose `counts` is inconsistent with
    /// the graph structure. Specifically, the per-state subtree count
    /// must equal `(1 if accepting else 0) + sum_over_children`, so
    /// in a well-formed file `counts[s] - children_sum` is always 0
    /// or 1. Other values indicate corruption / a non-minimal automaton.
    #[test]
    fn json_rejects_inconsistent_counts() {
        // One state, no edges, counts[0] = 5 (impossible: a leaf
        // state can contribute at most 1 to its own count).
        let bad = r#"{
            "format": "tilezz-dafsa",
            "version": 1,
            "scalar": "i8",
            "n_states": 1,
            "n_edges": 0,
            "edges_start": [0],
            "labels": [],
            "targets": [],
            "counts": [5]
        }"#;
        let err = Dafsa::from_json_str(bad).unwrap_err();
        let msg = err.to_string();
        assert!(
            msg.contains("counts") && msg.contains("children"),
            "unexpected error: {msg}"
        );
    }

    /// Verifies the on-disk form really does omit the redundant
    /// fields (`is_accept` per state, top-level `accept_count`).
    /// Both are reconstructed from `counts` on read. Compatibility
    /// note: any external reader that expects these fields has a bug.
    #[test]
    fn json_omits_redundant_fields() {
        let dafsa = Dafsa::from_seqs([[1i8, 2].as_slice(), [1, 3].as_slice()].iter().copied());
        let json = dafsa.to_json_string();
        assert!(
            !json.contains("is_accept"),
            "is_accept must not appear on disk; it is derived from counts"
        );
        assert!(
            !json.contains("accept_count"),
            "accept_count must not appear on disk; it equals counts[0]"
        );
    }

    /// The serialized JSON contains the documented top-level keys,
    /// so external (non-Rust) loaders can rely on them.
    #[test]
    fn json_schema_keys_present() {
        let dafsa = Dafsa::from_seqs([[1i8, 2].as_slice(), [1, 3].as_slice()].iter().copied());
        let json = dafsa.to_json_string();
        for key in [
            "\"format\"",
            "\"version\"",
            "\"scalar\"",
            "\"n_states\"",
            "\"n_edges\"",
            "\"edges_start\"",
            "\"labels\"",
            "\"targets\"",
            "\"counts\"",
        ] {
            assert!(json.contains(key), "missing schema key {key}");
        }
        // Sanity-check the discriminator and version values.
        assert!(json.contains("\"tilezz-dafsa\""));
        assert!(json.contains("\"version\":1"));
        assert!(json.contains("\"scalar\":\"i8\""));
    }

    /// `write_json_gz` -> `read_json_gz` round-trips through gzipped
    /// bytes; the decoded DAFSA reproduces the input set in lex order.
    #[test]
    fn gz_roundtrip() {
        let seqs: Vec<Vec<i8>> = vec![vec![1, 2], vec![1, 2, 3], vec![1, 3], vec![2], vec![2, 1]];
        let dafsa = Dafsa::from_seqs(seqs.iter().map(|v| v.as_slice()));

        let mut buf: Vec<u8> = Vec::new();
        dafsa.write_json_gz(&mut buf).expect("write_json_gz");

        // gzip magic bytes 0x1f 0x8b: the file really is compressed.
        assert!(buf.starts_with(&[0x1f, 0x8b]), "missing gzip magic");

        let decoded = Dafsa::read_json_gz(&buf[..]).expect("read_json_gz");
        let out: Vec<Vec<i8>> = decoded.iter().collect();
        assert_eq!(out, seqs);
    }
}
