//! `LazyRatDafsa`: a `RatDafsa`-shaped reader that fetches a
//! "tilezz-rat-dafsa-blocks" asset one block at a time, on demand.
//!
//! Use this when shipping a large rat set to a browser or any other
//! consumer that wants the DAFSA's query semantics without
//! materializing the whole automaton up front. The wire format
//! (one manifest file plus N gzipped block files) is described in
//! the bundled schema [`JSON_SCHEMA_DOC`].
//!
//! The reader is parameterised over a `fetch_block: Fn(usize) ->
//! io::Result<Vec<u8>>` callback. In tests we point it at a
//! filesystem directory; in production a WASM port points it at
//! `fetch(URL + template.format(id))`. The same `LazyRatDafsa`
//! type works for both.
//!
//! # API parity with [`RatDafsa`]
//!
//! [`LazyRatDafsa`] exposes the same query surface as
//! [`crate::stringmatch::RatDafsa`] minus `iter`: `len`,
//! `contains`, `get`, `index_of`. `iter` is intentionally absent --
//! it would force a full sweep of every block, defeating the
//! purpose. Consumers that want every rat should ship the
//! single-file `tilezz-rat-dafsa` JSON form instead.
//!
//! Length-prefix encoding is hidden from the caller exactly as in
//! [`RatDafsa`]: queries take and return raw rats; the length byte
//! lives only in the on-disk stored sequences.

use std::cell::RefCell;
use std::collections::HashMap;
use std::io;

use serde::{Deserialize, Serialize};

use super::rat::RatDafsa;

pub const JSON_SCHEMA_DOC: &str = include_str!("blocks_schema.txt");

const MANIFEST_FORMAT: &str = "tilezz-rat-dafsa-blocks";
const MANIFEST_VERSION: u32 = 1;
const BLOCK_MAGIC: &[u8; 4] = b"TRB1";
const SCALAR_TAG: &str = "i8";
const BLOCK_FORMAT_TAG: &str = "tilezz-rat-block";
const BLOCK_FORMAT_VERSION: u32 = 1;
const DEFAULT_BLOCK_FILENAME_TEMPLATE: &str = "block_{:06}.bin";

/// Bytes per state record on disk. See the schema for the layout.
const STATE_RECORD_BYTES: usize = 16;
/// Bytes per edge record on disk. See the schema for the layout.
const EDGE_RECORD_BYTES: usize = 8;
/// Bytes in the block header (magic + 3 u32 counters).
const BLOCK_HEADER_BYTES: usize = 4 + 4 + 4 + 4;

/// The manifest file produced by [`RatDafsa::write_blocks`] and
/// parsed by [`LazyRatDafsa::open`]. See [`JSON_SCHEMA_DOC`] for the
/// authoritative wire description.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BlockManifest {
    pub format: String,
    pub version: u32,
    pub scalar: String,
    pub block_format: String,
    pub block_version: u32,
    pub block_size: u32,
    pub n_states: u32,
    pub n_edges: u32,
    pub n_sequences: u32,
    pub n_blocks: u32,
    pub root_state: u32,
    pub block_filename_template: String,
    /// Maximum rat length covered by the asset (= the largest length
    /// byte appearing as a root-edge label, i.e. the maximum of the
    /// length-prefix encoding). A consumer can short-circuit queries
    /// for longer inputs without walking the DAFSA at all. Optional
    /// for backwards compatibility with manifests emitted before this
    /// field was added; readers should treat absence as "unknown".
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_indexed_length: Option<u32>,
}

impl BlockManifest {
    /// Render the canonical filename of a block file under this
    /// manifest's template. Currently we only honour the `{:06}`
    /// substitution; other format specs are not supported in v1 and
    /// the reader will reject manifests whose template does not
    /// match that pattern.
    pub fn block_filename(&self, block_id: u32) -> String {
        // Tiny templating: replace the first "{...}" occurrence with
        // the zero-padded id. The only template shape we support in
        // v1 is `block_{:06}.bin`, but allowing any width keeps the
        // schema flexible without pulling in a full format parser.
        if let Some(start) = self.block_filename_template.find('{') {
            if let Some(end_rel) = self.block_filename_template[start..].find('}') {
                let end = start + end_rel;
                let spec = &self.block_filename_template[start + 1..end];
                let width: usize = spec
                    .strip_prefix(':')
                    .and_then(|s| s.strip_prefix('0'))
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(6);
                let rendered = format!("{:0>width$}", block_id, width = width);
                let mut s = String::new();
                s.push_str(&self.block_filename_template[..start]);
                s.push_str(&rendered);
                s.push_str(&self.block_filename_template[end + 1..]);
                return s;
            }
        }
        // No substitution slot; return the literal name and let the
        // caller decide whether that's an error.
        self.block_filename_template.clone()
    }
}

/// A decoded block held in the in-memory cache. Owns the state and
/// edge arrays for the contiguous state-id range it represents.
#[derive(Debug, Clone)]
struct Block {
    first_state_id: u32,
    n_states: u32,
    states: Vec<StateRec>,
    edges: Vec<EdgeRec>,
}

#[derive(Debug, Clone, Copy)]
struct StateRec {
    edges_offset: u32,
    count: u64,
    is_accept: bool,
}

#[derive(Debug, Clone, Copy)]
struct EdgeRec {
    label: i8,
    target: u32,
}

impl Block {
    /// Decode a gzipped block file produced by [`RatDafsa::write_blocks`].
    fn from_gz_bytes(bytes: &[u8]) -> io::Result<Self> {
        let mut buf = Vec::new();
        let mut dec = flate2::read::GzDecoder::new(bytes);
        io::Read::read_to_end(&mut dec, &mut buf)?;
        Self::from_bytes(&buf)
    }

    fn from_bytes(buf: &[u8]) -> io::Result<Self> {
        if buf.len() < BLOCK_HEADER_BYTES {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "block: truncated header",
            ));
        }
        let magic = &buf[..4];
        if magic != BLOCK_MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "block: bad magic (expected {:?}, got {:?})",
                    std::str::from_utf8(BLOCK_MAGIC).unwrap(),
                    String::from_utf8_lossy(magic)
                ),
            ));
        }
        let first_state_id = read_u32(buf, 4);
        let n_states = read_u32(buf, 8);
        let n_edges = read_u32(buf, 12);

        let states_start = BLOCK_HEADER_BYTES;
        let states_end = states_start + n_states as usize * STATE_RECORD_BYTES;
        let edges_start = states_end;
        let edges_end = edges_start + n_edges as usize * EDGE_RECORD_BYTES;
        if buf.len() < edges_end {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "block: truncated body",
            ));
        }

        let mut states = Vec::with_capacity(n_states as usize);
        for i in 0..n_states as usize {
            let o = states_start + i * STATE_RECORD_BYTES;
            let edges_offset = read_u32(buf, o);
            let count = read_u64(buf, o + 4);
            let is_accept = buf[o + 12] != 0;
            states.push(StateRec {
                edges_offset,
                count,
                is_accept,
            });
        }
        let mut edges = Vec::with_capacity(n_edges as usize);
        for i in 0..n_edges as usize {
            let o = edges_start + i * EDGE_RECORD_BYTES;
            let label = buf[o] as i8;
            let target = read_u32(buf, o + 4);
            edges.push(EdgeRec { label, target });
        }
        Ok(Block {
            first_state_id,
            n_states,
            states,
            edges,
        })
    }

    /// Edges of the i-th state in this block, as a slice into
    /// `self.edges` (CSR-style; the last state runs to `n_edges`).
    fn edges_for(&self, local_index: usize) -> &[EdgeRec] {
        let start = self.states[local_index].edges_offset as usize;
        let end = if local_index + 1 < self.n_states as usize {
            self.states[local_index + 1].edges_offset as usize
        } else {
            self.edges.len()
        };
        &self.edges[start..end]
    }
}

fn read_u32(buf: &[u8], off: usize) -> u32 {
    u32::from_le_bytes([buf[off], buf[off + 1], buf[off + 2], buf[off + 3]])
}
fn read_u64(buf: &[u8], off: usize) -> u64 {
    u64::from_le_bytes([
        buf[off],
        buf[off + 1],
        buf[off + 2],
        buf[off + 3],
        buf[off + 4],
        buf[off + 5],
        buf[off + 6],
        buf[off + 7],
    ])
}
fn write_u32(buf: &mut Vec<u8>, v: u32) {
    buf.extend_from_slice(&v.to_le_bytes());
}
fn write_u64(buf: &mut Vec<u8>, v: u64) {
    buf.extend_from_slice(&v.to_le_bytes());
}

/// `RatDafsa`-shaped query API that fetches blocks lazily via a
/// caller-supplied callback. See module docs.
pub struct LazyRatDafsa<F>
where
    F: Fn(u32) -> io::Result<Vec<u8>>,
{
    manifest: BlockManifest,
    fetch_block: F,
    /// Interior mutability: queries are `&self` (so the same reader
    /// can be shared / borrowed during a walk) but they populate the
    /// block cache as a side effect. `RefCell` is fine because all
    /// of `LazyRatDafsa`'s public methods take `&self` and there is
    /// no concurrent re-entrancy: a `contains` call cannot recurse
    /// into another `contains` on the same reader.
    cache: RefCell<HashMap<u32, Block>>,
}

impl<F> LazyRatDafsa<F>
where
    F: Fn(u32) -> io::Result<Vec<u8>>,
{
    /// Parse `manifest_json` (the contents of `block_index.json`)
    /// and wrap it with the provided fetch callback. No blocks are
    /// fetched yet; the first walk pulls them in lazily.
    pub fn open(manifest_json: &str, fetch_block: F) -> io::Result<Self> {
        let manifest: BlockManifest = serde_json::from_str(manifest_json)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        if manifest.format != MANIFEST_FORMAT {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "LazyRatDafsa: bad manifest format (expected '{}', got '{}')",
                    MANIFEST_FORMAT, manifest.format
                ),
            ));
        }
        if manifest.version != MANIFEST_VERSION {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "LazyRatDafsa: unsupported manifest version (expected {}, got {})",
                    MANIFEST_VERSION, manifest.version
                ),
            ));
        }
        if manifest.scalar != SCALAR_TAG {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "LazyRatDafsa: scalar mismatch (file is '{}', expected '{}')",
                    manifest.scalar, SCALAR_TAG
                ),
            ));
        }
        if manifest.block_format != BLOCK_FORMAT_TAG {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "LazyRatDafsa: bad block format tag (expected '{}', got '{}')",
                    BLOCK_FORMAT_TAG, manifest.block_format
                ),
            ));
        }
        if manifest.block_version != BLOCK_FORMAT_VERSION {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "LazyRatDafsa: unsupported block format version (expected {}, got {})",
                    BLOCK_FORMAT_VERSION, manifest.block_version
                ),
            ));
        }
        if manifest.block_size == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "LazyRatDafsa: block_size must be >= 1",
            ));
        }
        Ok(Self {
            manifest,
            fetch_block,
            cache: RefCell::new(HashMap::new()),
        })
    }

    /// Number of rats stored.
    pub fn len(&self) -> usize {
        self.manifest.n_sequences as usize
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Set membership.
    pub fn contains<S: AsRef<[i8]>>(&self, rat: S) -> bool {
        let prefixed = prefix(rat.as_ref());
        self.walk(&prefixed)
            .ok()
            .is_some_and(|(rec, _)| rec.is_accept)
    }

    /// Assigned index of `rat` under `(length, lex)` ordering, or
    /// `None` if not in the set.
    pub fn index_of<S: AsRef<[i8]>>(&self, rat: S) -> Option<u32> {
        let prefixed = prefix(rat.as_ref());
        // The inner DAFSA's lex rank of the prefixed sequence equals
        // the rat's (length, lex) rank.
        let mut state = self.manifest.root_state;
        let mut rank: u32 = 0;
        for &symbol in &prefixed {
            let (rec, edges) = self.lookup_state(state).ok()?;
            if rec.is_accept {
                rank = rank.checked_add(1)?;
            }
            let mut found = None;
            for e in edges {
                match e.label.cmp(&symbol) {
                    std::cmp::Ordering::Less => {
                        let (sub_rec, _) = self.lookup_state(e.target).ok()?;
                        rank = rank.checked_add(sub_rec.count as u32)?;
                    }
                    std::cmp::Ordering::Equal => {
                        found = Some(e.target);
                        break;
                    }
                    std::cmp::Ordering::Greater => return None,
                }
            }
            state = found?;
        }
        let (rec, _) = self.lookup_state(state).ok()?;
        if !rec.is_accept {
            return None;
        }
        Some(rank)
    }

    /// The rat at assigned index `i` in `(length, lex)` order, or
    /// `None` if `i >= len()`.
    pub fn get(&self, i: usize) -> Option<Vec<i8>> {
        if i >= self.len() {
            return None;
        }
        let mut out: Vec<i8> = Vec::new();
        let mut state = self.manifest.root_state;
        let mut remaining = i as u64;
        loop {
            let (rec, edges) = self.lookup_state(state).ok()?;
            if rec.is_accept {
                if remaining == 0 {
                    // `out` currently holds the length-prefixed
                    // sequence; strip the length byte before
                    // returning.
                    if !out.is_empty() {
                        out.remove(0);
                    }
                    return Some(out);
                }
                remaining -= 1;
            }
            let mut descended = false;
            for e in edges {
                let (sub_rec, _) = self.lookup_state(e.target).ok()?;
                let n = sub_rec.count;
                if remaining < n {
                    out.push(e.label);
                    state = e.target;
                    descended = true;
                    break;
                }
                remaining -= n;
            }
            if !descended {
                // Well-formed DAFSA with `i < len()` always descends;
                // defensive fallback returns None.
                return None;
            }
        }
    }

    /// Materialise the full set of rats and rebuild a single-file
    /// [`RatDafsa`]. Walks every block in the asset (DFS over the
    /// underlying length-prefixed DAFSA), so it loads everything into
    /// memory -- this is the eager "decompress the lazy view" path,
    /// intended for restoring the JSON-form artifact from a blocked
    /// asset or for batch processing. For browser query workloads
    /// keep using the lazy `contains` / `get` / `index_of` methods.
    pub fn to_rat_dafsa(&self) -> io::Result<RatDafsa> {
        // Collect rats in the natural DFS (lex on prefixed = length,
        // lex on raw) order. `RatDafsa::from_rats` then re-sorts and
        // dedups, but the input is already in the right shape so
        // that step is a fast no-op for this input.
        let mut rats: Vec<Vec<i8>> = Vec::with_capacity(self.len());
        let mut current: Vec<i8> = Vec::new();
        let mut stack: Vec<(Vec<EdgeRec>, usize)> = Vec::new();

        let (root_rec, root_edges) = self.lookup_state(self.manifest.root_state)?;
        if root_rec.is_accept {
            // Defensive: a root accept implies an empty prefixed
            // sequence, which `RatDafsa::from_rats` cannot legally
            // represent (a rat needs at least its length byte). Skip
            // rather than emit a bogus zero-length rat.
        }
        stack.push((root_edges, 0));

        while let Some(frame) = stack.last_mut() {
            if frame.1 < frame.0.len() {
                let edge = frame.0[frame.1];
                frame.1 += 1;
                current.push(edge.label);
                let (rec, child_edges) = self.lookup_state(edge.target)?;
                if rec.is_accept {
                    // `current` is the full length-prefixed sequence;
                    // strip the leading length byte to recover the rat.
                    let mut rat = current.clone();
                    rat.remove(0);
                    rats.push(rat);
                }
                stack.push((child_edges, 0));
            } else {
                stack.pop();
                current.pop();
            }
        }

        debug_assert_eq!(
            rats.len(),
            self.len(),
            "to_rat_dafsa: DFS yielded {} rats but manifest claims {}",
            rats.len(),
            self.len()
        );
        Ok(RatDafsa::from_rats(rats))
    }

    /// Walk the entire prefixed sequence. Returns the terminal
    /// state's record + edges (whether or not it accepts). Caller
    /// inspects `rec.is_accept` to decide membership. Returns
    /// `Err(())` if the walk falls off the language at any step.
    fn walk(&self, prefixed: &[i8]) -> Result<(StateRec, Vec<EdgeRec>), ()> {
        let mut state = self.manifest.root_state;
        for &symbol in prefixed {
            let (_, edges) = self.lookup_state(state).map_err(|_| ())?;
            let mut next: Option<u32> = None;
            for e in &edges {
                match e.label.cmp(&symbol) {
                    std::cmp::Ordering::Less => continue,
                    std::cmp::Ordering::Equal => {
                        next = Some(e.target);
                        break;
                    }
                    std::cmp::Ordering::Greater => return Err(()),
                }
            }
            state = next.ok_or(())?;
        }
        self.lookup_state(state).map_err(|_| ())
    }

    /// Fetch the state record and a snapshot of its edges for
    /// `state_id`. Each call independently borrows the block cache
    /// for the duration of the copy, so callers can chain
    /// `lookup_state` calls without RefCell-borrow conflicts.
    fn lookup_state(&self, state_id: u32) -> io::Result<(StateRec, Vec<EdgeRec>)> {
        let block_id = state_id / self.manifest.block_size;
        let local = (state_id % self.manifest.block_size) as usize;
        self.ensure_block_loaded(block_id)?;
        let cache = self.cache.borrow();
        let block = cache
            .get(&block_id)
            .expect("ensure_block_loaded just populated it");
        let rec = block.states[local];
        let edges = block.edges_for(local).to_vec();
        Ok((rec, edges))
    }

    fn ensure_block_loaded(&self, block_id: u32) -> io::Result<()> {
        if self.cache.borrow().contains_key(&block_id) {
            return Ok(());
        }
        let bytes = (self.fetch_block)(block_id)?;
        let block = Block::from_gz_bytes(&bytes)?;
        let expected_first = block_id * self.manifest.block_size;
        if block.first_state_id != expected_first {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "block {block_id} declares first_state_id {} but should be {expected_first}",
                    block.first_state_id
                ),
            ));
        }
        self.cache.borrow_mut().insert(block_id, block);
        Ok(())
    }
}

// ---- Async (future-based) reader ----

/// Mirror of [`LazyRatDafsa`] for environments where block fetching is
/// inherently asynchronous (browsers / WASM / async Rust runtimes).
/// The struct is *not* generic over the fetcher type: each query
/// (`contains` / `get` / `index_of`) takes the fetcher as an
/// argument, so callers can construct one inline (e.g. capturing the
/// asset URL prefix from outer scope) without having to store an
/// owned closure inside the reader. Block cache and validation
/// semantics match the sync version exactly.
///
/// Borrows on the block cache never cross `.await`: each lookup pulls
/// out an owned `(StateRec, Vec<EdgeRec>)` copy of the relevant
/// state's data and drops the borrow before the next await point, so
/// the type is sound under single-threaded WASM concurrency where
/// multiple in-flight queries interleave at await boundaries.
pub struct LazyRatDafsaAsync {
    manifest: BlockManifest,
    cache: RefCell<HashMap<u32, Block>>,
}

impl LazyRatDafsaAsync {
    /// Parse the JSON manifest and construct an empty-cache reader.
    /// No blocks are fetched until the first query runs.
    pub fn open(manifest_json: &str) -> io::Result<Self> {
        let manifest: BlockManifest = serde_json::from_str(manifest_json)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        if manifest.format != MANIFEST_FORMAT {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "LazyRatDafsaAsync: bad manifest format (expected '{}', got '{}')",
                    MANIFEST_FORMAT, manifest.format
                ),
            ));
        }
        if manifest.version != MANIFEST_VERSION {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "LazyRatDafsaAsync: unsupported manifest version (expected {}, got {})",
                    MANIFEST_VERSION, manifest.version
                ),
            ));
        }
        if manifest.scalar != SCALAR_TAG {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "LazyRatDafsaAsync: scalar mismatch (file is '{}', expected '{}')",
                    manifest.scalar, SCALAR_TAG
                ),
            ));
        }
        if manifest.block_format != BLOCK_FORMAT_TAG {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "LazyRatDafsaAsync: bad block format tag (expected '{}', got '{}')",
                    BLOCK_FORMAT_TAG, manifest.block_format
                ),
            ));
        }
        if manifest.block_version != BLOCK_FORMAT_VERSION {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "LazyRatDafsaAsync: unsupported block format version (expected {}, got {})",
                    BLOCK_FORMAT_VERSION, manifest.block_version
                ),
            ));
        }
        if manifest.block_size == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "LazyRatDafsaAsync: block_size must be >= 1",
            ));
        }
        Ok(Self {
            manifest,
            cache: RefCell::new(HashMap::new()),
        })
    }

    pub fn manifest(&self) -> &BlockManifest {
        &self.manifest
    }

    pub fn len(&self) -> usize {
        self.manifest.n_sequences as usize
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Set membership: is `rat` in the stored set?
    pub async fn contains<S, F, Fut>(&self, rat: S, fetch: &F) -> bool
    where
        S: AsRef<[i8]>,
        F: Fn(u32) -> Fut,
        Fut: core::future::Future<Output = io::Result<Vec<u8>>>,
    {
        let prefixed = prefix(rat.as_ref());
        match self.walk(&prefixed, fetch).await {
            Ok((rec, _)) => rec.is_accept,
            Err(_) => false,
        }
    }

    /// Assigned (length, lex)-rank index of `rat`, or `None` if it
    /// is not in the set or is longer than the asset's
    /// `max_indexed_length` (when that hint is present).
    pub async fn index_of<S, F, Fut>(&self, rat: S, fetch: &F) -> Option<u32>
    where
        S: AsRef<[i8]>,
        F: Fn(u32) -> Fut,
        Fut: core::future::Future<Output = io::Result<Vec<u8>>>,
    {
        let rat_slice = rat.as_ref();
        if let Some(max) = self.manifest.max_indexed_length {
            if rat_slice.len() as u32 > max {
                return None;
            }
        }
        let prefixed = prefix(rat_slice);
        let mut state = self.manifest.root_state;
        let mut rank: u32 = 0;
        for &symbol in &prefixed {
            let (rec, edges) = self.lookup_state(state, fetch).await.ok()?;
            if rec.is_accept {
                rank = rank.checked_add(1)?;
            }
            let mut next: Option<u32> = None;
            for e in edges {
                match e.label.cmp(&symbol) {
                    std::cmp::Ordering::Less => {
                        let (sub_rec, _) = self.lookup_state(e.target, fetch).await.ok()?;
                        rank = rank.checked_add(sub_rec.count as u32)?;
                    }
                    std::cmp::Ordering::Equal => {
                        next = Some(e.target);
                        break;
                    }
                    std::cmp::Ordering::Greater => return None,
                }
            }
            state = next?;
        }
        let (rec, _) = self.lookup_state(state, fetch).await.ok()?;
        if !rec.is_accept {
            return None;
        }
        Some(rank)
    }

    /// Rat at assigned index `i` in (length, lex) order.
    pub async fn get<F, Fut>(&self, i: usize, fetch: &F) -> Option<Vec<i8>>
    where
        F: Fn(u32) -> Fut,
        Fut: core::future::Future<Output = io::Result<Vec<u8>>>,
    {
        if i >= self.len() {
            return None;
        }
        let mut out: Vec<i8> = Vec::new();
        let mut state = self.manifest.root_state;
        let mut remaining = i as u64;
        loop {
            let (rec, edges) = self.lookup_state(state, fetch).await.ok()?;
            if rec.is_accept {
                if remaining == 0 {
                    if !out.is_empty() {
                        out.remove(0);
                    }
                    return Some(out);
                }
                remaining -= 1;
            }
            let mut descended = false;
            for e in edges {
                let (sub_rec, _) = self.lookup_state(e.target, fetch).await.ok()?;
                let n = sub_rec.count;
                if remaining < n {
                    out.push(e.label);
                    state = e.target;
                    descended = true;
                    break;
                }
                remaining -= n;
            }
            if !descended {
                return None;
            }
        }
    }

    /// Walk the entire prefixed sequence; return the terminal state's
    /// record + edges (whether or not it accepts). Used by
    /// [`Self::contains`].
    async fn walk<F, Fut>(&self, prefixed: &[i8], fetch: &F) -> Result<(StateRec, Vec<EdgeRec>), ()>
    where
        F: Fn(u32) -> Fut,
        Fut: core::future::Future<Output = io::Result<Vec<u8>>>,
    {
        let mut state = self.manifest.root_state;
        for &symbol in prefixed {
            let (_, edges) = self.lookup_state(state, fetch).await.map_err(|_| ())?;
            let mut next: Option<u32> = None;
            for e in &edges {
                match e.label.cmp(&symbol) {
                    std::cmp::Ordering::Less => continue,
                    std::cmp::Ordering::Equal => {
                        next = Some(e.target);
                        break;
                    }
                    std::cmp::Ordering::Greater => return Err(()),
                }
            }
            state = next.ok_or(())?;
        }
        self.lookup_state(state, fetch).await.map_err(|_| ())
    }

    async fn lookup_state<F, Fut>(
        &self,
        state_id: u32,
        fetch: &F,
    ) -> io::Result<(StateRec, Vec<EdgeRec>)>
    where
        F: Fn(u32) -> Fut,
        Fut: core::future::Future<Output = io::Result<Vec<u8>>>,
    {
        let block_id = state_id / self.manifest.block_size;
        let local = (state_id % self.manifest.block_size) as usize;
        self.ensure_block_loaded(block_id, fetch).await?;
        let cache = self.cache.borrow();
        let block = cache
            .get(&block_id)
            .expect("ensure_block_loaded just populated it");
        let rec = block.states[local];
        let edges = block.edges_for(local).to_vec();
        Ok((rec, edges))
    }

    async fn ensure_block_loaded<F, Fut>(&self, block_id: u32, fetch: &F) -> io::Result<()>
    where
        F: Fn(u32) -> Fut,
        Fut: core::future::Future<Output = io::Result<Vec<u8>>>,
    {
        // Quick check; we drop the borrow before awaiting so the
        // RefCell doesn't bridge an await point.
        if self.cache.borrow().contains_key(&block_id) {
            return Ok(());
        }
        let bytes = fetch(block_id).await?;
        let block = Block::from_gz_bytes(&bytes)?;
        let expected_first = block_id * self.manifest.block_size;
        if block.first_state_id != expected_first {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "block {block_id} declares first_state_id {} but should be {expected_first}",
                    block.first_state_id
                ),
            ));
        }
        self.cache.borrow_mut().insert(block_id, block);
        Ok(())
    }
}

// ---- Writer ----

impl RatDafsa {
    /// Slice this `RatDafsa` into a block_index + per-block files
    /// suitable for serving from a static-asset host. See
    /// [`JSON_SCHEMA_DOC`] for the on-disk layout and
    /// [`LazyRatDafsa`] for the corresponding reader.
    ///
    /// `dir` must already exist. `block_size` is the number of
    /// states per file (except the last); typical values are
    /// 1024..8192.
    pub fn write_blocks(&self, dir: &std::path::Path, block_size: u32) -> io::Result<()> {
        assert!(block_size >= 1, "block_size must be >= 1");
        let dafsa = self.inner();
        let n_states = dafsa.raw_n_states() as u32;
        let n_edges = dafsa.raw_n_edges() as u32;
        let counts = dafsa.raw_counts();
        let n_sequences = counts.first().copied().unwrap_or(0) as u32;
        let n_blocks = n_states.div_ceil(block_size);

        let edges_start = dafsa.raw_edges_start();
        let labels = dafsa.raw_labels();
        let targets = dafsa.raw_targets();
        let is_accept = dafsa.raw_is_accept();

        for block_id in 0..n_blocks {
            let first = block_id * block_size;
            let last = ((block_id + 1) * block_size).min(n_states);
            let n_in_block = last - first;

            // Slice the edges owned by this block's states.
            let edge_lo = edges_start[first as usize] as usize;
            let edge_hi = if last as usize == n_states as usize {
                n_edges as usize
            } else {
                edges_start[last as usize] as usize
            };
            let n_in_edges = (edge_hi - edge_lo) as u32;

            let mut bytes = Vec::with_capacity(
                BLOCK_HEADER_BYTES
                    + n_in_block as usize * STATE_RECORD_BYTES
                    + n_in_edges as usize * EDGE_RECORD_BYTES,
            );
            bytes.extend_from_slice(BLOCK_MAGIC);
            write_u32(&mut bytes, first);
            write_u32(&mut bytes, n_in_block);
            write_u32(&mut bytes, n_in_edges);

            // State records: edges_offset is block-local (subtract
            // edge_lo from the global edges_start).
            for local in 0..n_in_block as usize {
                let global = first as usize + local;
                let edges_offset_local = (edges_start[global] as usize - edge_lo) as u32;
                write_u32(&mut bytes, edges_offset_local);
                write_u64(&mut bytes, counts[global] as u64);
                bytes.push(is_accept[global] as u8);
                bytes.push(0);
                bytes.push(0);
                bytes.push(0);
            }
            // Edge records.
            for i in edge_lo..edge_hi {
                bytes.push(labels[i] as u8);
                bytes.push(0);
                bytes.push(0);
                bytes.push(0);
                write_u32(&mut bytes, targets[i]);
            }

            // Gzip and write.
            let filename = render_template(DEFAULT_BLOCK_FILENAME_TEMPLATE, block_id);
            let path = dir.join(filename);
            let file = std::fs::File::create(&path)?;
            let mut enc = flate2::write::GzEncoder::new(file, flate2::Compression::default());
            io::Write::write_all(&mut enc, &bytes)?;
            enc.finish()?;
        }

        // Write the manifest last so the reader never sees an index
        // pointing at a block that hasn't been flushed yet.
        // Maximum rat length covered = the largest label on a root
        // edge (each root edge's label is the length-prefix byte of
        // a stored rat). Lets consumers short-circuit oversized
        // queries without walking.
        let root_first = edges_start[0] as usize;
        let root_last = if n_states > 1 {
            edges_start[1] as usize
        } else {
            n_edges as usize
        };
        let max_indexed_length = labels[root_first..root_last]
            .iter()
            .copied()
            .max()
            .map(|v| v as u32);

        let manifest = BlockManifest {
            format: MANIFEST_FORMAT.to_string(),
            version: MANIFEST_VERSION,
            scalar: SCALAR_TAG.to_string(),
            block_format: BLOCK_FORMAT_TAG.to_string(),
            block_version: BLOCK_FORMAT_VERSION,
            block_size,
            n_states,
            n_edges,
            n_sequences,
            n_blocks,
            root_state: 0,
            block_filename_template: DEFAULT_BLOCK_FILENAME_TEMPLATE.to_string(),
            max_indexed_length,
        };
        let manifest_path = dir.join("block_index.json");
        let manifest_file = std::fs::File::create(&manifest_path)?;
        serde_json::to_writer_pretty(manifest_file, &manifest)
            .map_err(|e| io::Error::other(format!("write manifest: {e}")))?;
        Ok(())
    }
}

fn render_template(template: &str, block_id: u32) -> String {
    if let Some(start) = template.find('{') {
        if let Some(end_rel) = template[start..].find('}') {
            let end = start + end_rel;
            let spec = &template[start + 1..end];
            let width: usize = spec
                .strip_prefix(':')
                .and_then(|s| s.strip_prefix('0'))
                .and_then(|s| s.parse().ok())
                .unwrap_or(6);
            let rendered = format!("{:0>width$}", block_id, width = width);
            let mut s = String::new();
            s.push_str(&template[..start]);
            s.push_str(&rendered);
            s.push_str(&template[end + 1..]);
            return s;
        }
    }
    template.to_string()
}

fn prefix(rat: &[i8]) -> Vec<i8> {
    assert!(
        rat.len() <= i8::MAX as usize,
        "LazyRatDafsa: rat length {} exceeds the single-byte length-prefix limit ({})",
        rat.len(),
        i8::MAX as usize
    );
    let mut v = Vec::with_capacity(rat.len() + 1);
    v.push(rat.len() as i8);
    v.extend_from_slice(rat);
    v
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    /// Build a small `RatDafsa`, write it as blocks to a tempdir,
    /// open as `LazyRatDafsa` via a filesystem-backed callback, and
    /// assert that `len`, `contains`, `get`, and `index_of` all
    /// agree with the source `RatDafsa` across every entry.
    #[test]
    fn blocks_roundtrip_matches_ratdafsa() {
        // Mix of lengths + shared prefixes so the underlying DAFSA
        // actually has interesting structure to slice.
        let input: Vec<Vec<i8>> = vec![
            vec![-2],
            vec![-1],
            vec![0],
            vec![1],
            vec![2],
            vec![-2, 1, 1],
            vec![-1, 0, 1],
            vec![0, 1, 2],
            vec![1, 2, 3],
            vec![1, 2, 3, 4],
            vec![1, 2, 4],
            vec![1, 3, 5],
            vec![2, 2, 2],
        ];
        let rd = RatDafsa::from_rats(input.iter().map(|v| v.as_slice()));

        let dir = tempdir();
        rd.write_blocks(&dir, 4).expect("write_blocks");
        // Manifest + at least one block file got written.
        let manifest_path = dir.join("block_index.json");
        assert!(manifest_path.exists(), "manifest must be created");
        let manifest_text = fs::read_to_string(&manifest_path).unwrap();
        let manifest: BlockManifest = serde_json::from_str(&manifest_text).unwrap();
        assert!(manifest.n_blocks >= 1);
        // A few blocks ought to be created for non-trivial input.
        assert!(
            manifest.n_blocks >= 2,
            "expected >=2 blocks for {} states, got {}",
            manifest.n_states,
            manifest.n_blocks
        );

        // Open through a filesystem fetch callback.
        let dir_for_cb = dir.clone();
        let manifest_clone = manifest.clone();
        let fetch = move |block_id: u32| -> io::Result<Vec<u8>> {
            let name = manifest_clone.block_filename(block_id);
            fs::read(dir_for_cb.join(name))
        };
        let lazy = LazyRatDafsa::open(&manifest_text, fetch).expect("open");

        // len agrees.
        assert_eq!(lazy.len(), rd.len());

        // get(i) agrees for every i.
        for i in 0..rd.len() {
            assert_eq!(lazy.get(i), rd.get(i), "get({i}) mismatch");
        }
        assert_eq!(lazy.get(rd.len()), None, "out-of-range");

        // contains + index_of agree for every entry and for a
        // negative probe.
        for i in 0..rd.len() {
            let rat = rd.get(i).unwrap();
            assert!(lazy.contains(rat.as_slice()), "missing rat #{i}");
            assert_eq!(lazy.index_of(rat.as_slice()), Some(i as u32));
        }
        assert!(!lazy.contains([42i8, 42, 42].as_slice()));
        assert_eq!(lazy.index_of([42i8, 42, 42].as_slice()), None);
    }

    /// With `block_size = 1` every state lives in its own block, so
    /// any non-trivial walk forces fetches across multiple block
    /// boundaries. Guards the cross-block edge-traversal path.
    #[test]
    fn cross_block_walks_are_correct() {
        let input: Vec<Vec<i8>> = vec![
            vec![-2, 1, 1],
            vec![1, 2, 3, 4],
            vec![1, 2, 4],
            vec![0, 0, 0],
        ];
        let rd = RatDafsa::from_rats(input.iter().map(|v| v.as_slice()));

        let dir = tempdir();
        rd.write_blocks(&dir, 1).expect("write_blocks");
        let manifest_text = fs::read_to_string(dir.join("block_index.json")).unwrap();
        let manifest: BlockManifest = serde_json::from_str(&manifest_text).unwrap();
        assert_eq!(manifest.block_size, 1);
        assert_eq!(manifest.n_blocks, manifest.n_states);

        let dir_for_cb = dir.clone();
        let manifest_clone = manifest.clone();
        let fetch = move |id: u32| fs::read(dir_for_cb.join(manifest_clone.block_filename(id)));
        let lazy = LazyRatDafsa::open(&manifest_text, fetch).expect("open");

        for i in 0..rd.len() {
            let rat = rd.get(i).unwrap();
            assert!(lazy.contains(rat.as_slice()));
            assert_eq!(lazy.get(i), Some(rat.clone()));
            assert_eq!(lazy.index_of(rat.as_slice()), Some(i as u32));
        }
    }

    /// One query should not require fetching every block. Counts
    /// fetches and asserts that a `contains` of an n=4-symbol
    /// sequence touches significantly fewer blocks than the total
    /// for a non-trivial DAFSA with `block_size = 2`.
    #[test]
    fn fetch_callback_is_called_minimally() {
        let input: Vec<Vec<i8>> = vec![
            vec![1, 2, 3, 4],
            vec![1, 2, 3, 5],
            vec![1, 2, 4, 5],
            vec![1, 3, 4, 5],
            vec![2, 3, 4, 5],
            vec![-1, 0, 1, 2],
            vec![-2, -1, 0, 1],
            vec![0, 0, 0, 0],
        ];
        let rd = RatDafsa::from_rats(input.iter().map(|v| v.as_slice()));
        let dir = tempdir();
        rd.write_blocks(&dir, 2).expect("write_blocks");
        let manifest_text = fs::read_to_string(dir.join("block_index.json")).unwrap();
        let manifest: BlockManifest = serde_json::from_str(&manifest_text).unwrap();
        let total_blocks = manifest.n_blocks;

        let dir_for_cb = dir.clone();
        let manifest_clone = manifest.clone();
        let counter = std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0));
        let counter_for_cb = counter.clone();
        let fetch = move |id: u32| {
            counter_for_cb.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            fs::read(dir_for_cb.join(manifest_clone.block_filename(id)))
        };
        let lazy = LazyRatDafsa::open(&manifest_text, fetch).expect("open");

        // Single contains: walk a 4-symbol rat (5 with length
        // prefix). Touches at most ~6 distinct blocks even in the
        // worst case where every step changes block.
        assert!(lazy.contains([1i8, 2, 3, 4].as_slice()));
        let fetches = counter.load(std::sync::atomic::Ordering::Relaxed);
        assert!(
            fetches < total_blocks as usize,
            "single query fetched all {total_blocks} blocks, no locality"
        );

        // After the cached blocks are warm, a repeat query is free.
        counter.store(0, std::sync::atomic::Ordering::Relaxed);
        assert!(lazy.contains([1i8, 2, 3, 4].as_slice()));
        assert_eq!(
            counter.load(std::sync::atomic::Ordering::Relaxed),
            0,
            "warm repeat should hit cache only"
        );
    }

    /// `LazyRatDafsa::to_rat_dafsa` reverses `write_blocks`:
    /// build a RatDafsa, write blocks, open lazily, restore -- the
    /// restored RatDafsa must yield the same sequences in the same
    /// order as the original.
    #[test]
    fn to_rat_dafsa_round_trips() {
        let input: Vec<Vec<i8>> = vec![
            vec![-3, 2, 1],
            vec![-1, 0, 1],
            vec![0, 0, 0, 0],
            vec![1, 2, 3, 4, 5],
            vec![1, 2, 4, 5],
            vec![1],
            vec![2, 3],
            vec![3, 3, 3],
            vec![-2, -1, 0, 1],
        ];
        let original = RatDafsa::from_rats(input.iter().map(|v| v.as_slice()));

        let dir = tempdir();
        original.write_blocks(&dir, 3).expect("write_blocks");
        let manifest_text = fs::read_to_string(dir.join("block_index.json")).unwrap();
        let manifest: BlockManifest = serde_json::from_str(&manifest_text).unwrap();
        // Non-trivial blocking: more than one block, smaller-than-set.
        assert!(manifest.n_blocks >= 2);

        let dir_for_cb = dir.clone();
        let manifest_clone = manifest.clone();
        let fetch = move |id: u32| fs::read(dir_for_cb.join(manifest_clone.block_filename(id)));
        let lazy = LazyRatDafsa::open(&manifest_text, fetch).expect("open");

        let restored = lazy.to_rat_dafsa().expect("to_rat_dafsa");
        assert_eq!(restored.len(), original.len());

        // Same sequences, same order (both yield in (length, lex)
        // order). Compare via `iter` and via random-access `get`.
        let original_seqs: Vec<Vec<i8>> = original.iter().collect();
        let restored_seqs: Vec<Vec<i8>> = restored.iter().collect();
        assert_eq!(original_seqs, restored_seqs, "iter() mismatch");

        for i in 0..original.len() {
            assert_eq!(original.get(i), restored.get(i), "get({i}) mismatch");
        }

        // The restored RatDafsa's gzipped JSON form is interchangeable
        // with the original's. We compare observable behaviour rather
        // than byte-equal output because the JSON state numbering can
        // legitimately differ.
        let mut a = Vec::new();
        let mut b = Vec::new();
        original.write_json_gz(&mut a).unwrap();
        restored.write_json_gz(&mut b).unwrap();
        // Both files must round-trip back to the same query results.
        let a2 = RatDafsa::read_json_gz(&a[..]).unwrap();
        let b2 = RatDafsa::read_json_gz(&b[..]).unwrap();
        assert_eq!(a2.iter().collect::<Vec<_>>(), b2.iter().collect::<Vec<_>>());
    }

    /// Manifests with wrong format / version / scalar are rejected.
    #[test]
    fn manifest_validation() {
        let valid = BlockManifest {
            format: MANIFEST_FORMAT.to_string(),
            version: MANIFEST_VERSION,
            scalar: SCALAR_TAG.to_string(),
            block_format: BLOCK_FORMAT_TAG.to_string(),
            block_version: BLOCK_FORMAT_VERSION,
            block_size: 2,
            n_states: 3,
            n_edges: 2,
            n_sequences: 1,
            n_blocks: 2,
            root_state: 0,
            block_filename_template: DEFAULT_BLOCK_FILENAME_TEMPLATE.to_string(),
            max_indexed_length: None,
        };
        let unused = |_: u32| -> io::Result<Vec<u8>> {
            Err(io::Error::other(
                "fetch should never be called in this test",
            ))
        };

        // Bad format tag.
        let mut bad = valid.clone();
        bad.format = "tilezz-something-else".to_string();
        let json = serde_json::to_string(&bad).unwrap();
        assert!(LazyRatDafsa::open(&json, unused).is_err());

        // Bad version.
        let mut bad = valid.clone();
        bad.version = 999;
        let json = serde_json::to_string(&bad).unwrap();
        assert!(LazyRatDafsa::open(&json, unused).is_err());

        // Wrong scalar.
        let mut bad = valid.clone();
        bad.scalar = "u16".to_string();
        let json = serde_json::to_string(&bad).unwrap();
        assert!(LazyRatDafsa::open(&json, unused).is_err());

        // Wrong block_format.
        let mut bad = valid.clone();
        bad.block_format = "something-else".to_string();
        let json = serde_json::to_string(&bad).unwrap();
        assert!(LazyRatDafsa::open(&json, unused).is_err());

        // block_size == 0.
        let mut bad = valid.clone();
        bad.block_size = 0;
        let json = serde_json::to_string(&bad).unwrap();
        assert!(LazyRatDafsa::open(&json, unused).is_err());

        // The valid one must open.
        let json = serde_json::to_string(&valid).unwrap();
        assert!(LazyRatDafsa::open(&json, unused).is_ok());
    }

    // ---- Test infrastructure ----

    /// Create a fresh temp directory under the system tmp dir.
    /// Returned path is also cleaned up by the test runner when the
    /// process exits; we keep it simple (no `Drop`) because our
    /// runs are short-lived.
    fn tempdir() -> std::path::PathBuf {
        let base = std::env::temp_dir();
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let p = base.join(format!("tilezz-lazyratdafsa-test-{nanos}"));
        fs::create_dir_all(&p).expect("tempdir mkdir");
        p
    }
}
