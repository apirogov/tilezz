//! `LazyRatDafsa`: a `RatDafsa`-shaped reader that fetches a
//! "tilezz-rat-dafsa-blocks" asset one block at a time, on demand.
//!
//! Use this when shipping a large rat set to a browser or any other
//! consumer that wants the DAFSA's query semantics without
//! materializing the whole automaton up front. The wire format
//! (one manifest file plus N gzipped block files) is described in
//! the bundled schema [`JSON_SCHEMA_DOC`].
//!
//! The reader is parameterised over a `fetch_block: Fn(u32) ->
//! io::Result<Vec<u8>>` callback that receives a block index into
//! the manifest's `blocks` array. In tests we point it at a
//! filesystem directory; in production a WASM port resolves the
//! index to a URL via `manifest.block_filename(&manifest.blocks[i])`
//! and calls `fetch(URL_PREFIX + filename)`. The same
//! `LazyRatDafsa` type works for both.
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
use sha2::{Digest, Sha256};

use super::rat::RatDafsa;

pub const JSON_SCHEMA_DOC: &str = include_str!("schemas/blocks_schema.txt");

const MANIFEST_FORMAT: &str = "tilezz-rat-dafsa-blocks";
const MANIFEST_VERSION: u32 = 1;
const BLOCK_MAGIC: &[u8; 4] = b"TRB1";
const SCALAR_TAG: &str = "i8";
const BLOCK_FORMAT_TAG: &str = "tilezz-rat-block";
const BLOCK_FORMAT_VERSION: u32 = 1;

/// Default target on-disk block file size for the writer. Real
/// blocks may be a bit larger (we close a block once the buffer
/// crosses this threshold, after the next state's records are
/// appended) and the final block can be smaller. 1 MiB is the
/// sweet spot for HTTP-served assets: roughly a hundred files per
/// gigabyte of dataset, each large enough that per-file overhead
/// (gzip header, HTTP round-trip) is well amortised.
pub const DEFAULT_TARGET_BLOCK_BYTES: u32 = 1 << 20;

/// Bytes per state record on disk. See the schema for the layout.
const STATE_RECORD_BYTES: usize = 16;
/// Bytes per edge record on disk. See the schema for the layout.
const EDGE_RECORD_BYTES: usize = 8;
/// Bytes in the block header (magic + 3 u32 counters).
const BLOCK_HEADER_BYTES: usize = 4 + 4 + 4 + 4;

/// The manifest file produced by [`RatDafsa::write_blocks`] and
/// parsed by [`LazyRatDafsa::open`]. See [`JSON_SCHEMA_DOC`] for the
/// authoritative wire description.
///
/// Key design properties:
///
/// * **Root state in the manifest, not in a block file.** State 0
///   (the DAFSA root) is the one state whose record necessarily
///   changes when extending the asset to deeper perim (it gains
///   length-prefix edges for the new length classes). Keeping it
///   out of the block files means every block file's bytes are
///   determined purely by its inner-DAFSA content -- which is
///   structurally stable across forward extensions.
///
/// * **Content-addressed block filenames.** Each block file's
///   name is the SHA-256 of its gzipped bytes. Two assets that
///   contain the same block (= same content range) reference the
///   same URL. CDN caches keyed on URL stay valid forever; an
///   asset extension that re-cuts the last partial block just
///   produces a new SHA-256 / new URL, leaving every other URL
///   untouched.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BlockManifest {
    pub format: String,
    pub version: u32,
    pub scalar: String,
    pub block_format: String,
    pub block_version: u32,
    /// Writer's target on-disk byte budget per block. Hint for the
    /// build pipeline; readers don't enforce it. Real blocks vary
    /// (the writer closes a block when its serialised size crosses
    /// the threshold; the final block may be smaller).
    pub target_block_bytes: u32,
    /// Total DAFSA state count, including the root. Block files
    /// cover states 1..n_states; state 0 lives in `root`.
    pub n_states: u32,
    /// Total DAFSA edges, including the root's edges.
    pub n_edges: u32,
    /// Number of accepted sequences (rats).
    pub n_sequences: u32,
    /// Maximum rat length covered by the asset (= the largest
    /// label appearing on a root-edge, i.e. the maximum length
    /// prefix). A consumer can short-circuit queries for longer
    /// inputs without walking the DAFSA at all.
    pub max_indexed_length: u32,
    /// Root state record (state 0). Lives in the manifest rather
    /// than in a block file so adding length classes doesn't
    /// reshape any `blocks/*.bin`.
    pub root: RootState,
    /// Block index. Each entry names the block's first state ID
    /// and its SHA-256 content hash (which is also its filename
    /// stem in `blocks/`). Sorted strictly ascending by
    /// `first_state`. Block N covers states
    /// `[blocks[N].first_state, blocks[N+1].first_state)`; the
    /// last block covers up to `n_states`.
    pub blocks: Vec<BlockEntry>,
    /// Optional URL prefix for fetching block files. When `Some`,
    /// readers resolve each block as `{block_base_url}{sha256}.bin`
    /// instead of `{asset_dir}/blocks/{sha256}.bin`. Lets the
    /// (small) manifest + RO-Crate live on one host (e.g. GitHub
    /// Pages) while the (large) block files live somewhere with
    /// fewer size constraints -- a GitHub Release, Cloudflare R2,
    /// IPFS, etc. The reader's integrity check (verify SHA-256
    /// against the manifest) is unchanged: the trust anchor stays
    /// the manifest URL, the block-content provenance follows from
    /// the content hash.
    ///
    /// Skipped in serialised output when `None` so existing
    /// single-host assets stay byte-identical.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub block_base_url: Option<String>,
}

/// State 0 (root) record, inlined into the manifest. Same shape
/// as a [`StateRec`] but rendered as JSON so a consumer that only
/// fetched the manifest can already enumerate the language's
/// length classes (via [`RootState::edges`]).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RootState {
    /// Number of accepted sequences reachable from the root.
    /// Equal to `n_sequences`.
    pub count: u64,
    /// Whether the root itself is an accepting state. For a
    /// `tilezz-rat-dafsa` (length-prefixed sequences) this is
    /// always `false` since every rat has at least the length
    /// byte; recorded explicitly for forward compatibility.
    pub is_accept: bool,
    /// Root's outgoing edges. For a `tilezz-rat-dafsa`, each
    /// edge's `label` is the length byte of one length class and
    /// `target` is the entry state of that class's sub-DAFSA.
    pub edges: Vec<RootEdge>,
}

/// One root edge in the manifest. Same shape as an [`EdgeRec`]
/// but rendered as JSON.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RootEdge {
    pub label: i8,
    pub target: u32,
}

/// One entry in the manifest's `blocks` array.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BlockEntry {
    /// Smallest state ID this block contains. Strictly greater
    /// than the previous entry's `first_state` (sorted).
    pub first_state: u32,
    /// Lowercase hex SHA-256 of the gzipped block file. Doubles
    /// as the filename stem: the block file lives at
    /// `blocks/<sha256>.bin`.
    pub sha256: String,
    /// Gzipped block file size in bytes. Lets clients show
    /// progress / preallocate buffers before fetching.
    pub size: u32,
}

impl BlockManifest {
    /// Canonical filename of one block file, relative to the asset
    /// root. The SHA-256 stem is the same one used as the cache /
    /// CDN key, so two assets that contain the same block file
    /// share the same URL.
    ///
    /// Returns the relative `blocks/<sha256>.bin` form regardless
    /// of whether the manifest carries a `block_base_url` -- it's
    /// the on-disk filename. For the actual fetch target (which
    /// honours `block_base_url`) see [`Self::block_url`].
    pub fn block_filename(&self, entry: &BlockEntry) -> String {
        format!("blocks/{}.bin", entry.sha256)
    }

    /// Where to fetch `entry`'s block file from. With
    /// `block_base_url` set, returns `{block_base_url}{sha256}.bin`
    /// (no `blocks/` prefix -- the URL points directly at the
    /// content-addressed object on its host of choice). Without it,
    /// returns the relative `blocks/<sha256>.bin`, which the caller
    /// resolves against the asset directory.
    pub fn block_url(&self, entry: &BlockEntry) -> String {
        match self.block_base_url.as_deref() {
            Some(base) => format!("{}{}.bin", base, entry.sha256),
            None => self.block_filename(entry),
        }
    }

    /// Find which block holds `state_id`. Returns the index into
    /// `self.blocks`, or `None` if `state_id == 0` (root, in the
    /// manifest) or `state_id >= n_states`. O(log blocks.len()).
    pub fn block_index_for_state(&self, state_id: u32) -> Option<usize> {
        if state_id == 0 || state_id >= self.n_states {
            return None;
        }
        // Largest i with blocks[i].first_state <= state_id. Since
        // `blocks` is sorted and `first_state` starts at 1, the
        // search returns at least 0 for any state_id >= 1.
        let pos = self.blocks.partition_point(|b| b.first_state <= state_id);
        // pos is the count of blocks whose first_state <= state_id;
        // we want the largest such index, which is pos - 1.
        if pos == 0 {
            // No block has first_state <= state_id. Manifest is
            // malformed (or state_id is < blocks[0].first_state,
            // which means state_id == 0 -- already filtered above).
            None
        } else {
            Some(pos - 1)
        }
    }

    /// The state range a given block covers: `[first_state,
    /// end_state)`. Convenience for readers.
    pub fn block_state_range(&self, block_index: usize) -> Option<(u32, u32)> {
        let entry = self.blocks.get(block_index)?;
        let end = self
            .blocks
            .get(block_index + 1)
            .map(|next| next.first_state)
            .unwrap_or(self.n_states);
        Some((entry.first_state, end))
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
        validate_block_index(&manifest)?;
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
        let mut state = 0u32;
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
        let mut state = 0u32;
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

        let (root_rec, root_edges) = self.lookup_state(0u32)?;
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
        let mut state = 0u32;
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
        // State 0 (root) lives in the manifest; serve it without
        // touching the block cache.
        if state_id == 0 {
            return Ok((
                root_as_state_rec(&self.manifest.root),
                root_edges(&self.manifest.root),
            ));
        }
        let block_index = self
            .manifest
            .block_index_for_state(state_id)
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "state_id {state_id} not within asset (n_states = {})",
                        self.manifest.n_states
                    ),
                )
            })?;
        let entry = &self.manifest.blocks[block_index];
        let local = (state_id - entry.first_state) as usize;
        self.ensure_block_loaded(block_index)?;
        let cache = self.cache.borrow();
        let block = cache
            .get(&(block_index as u32))
            .expect("ensure_block_loaded just populated it");
        let rec = block.states[local];
        let edges = block.edges_for(local).to_vec();
        Ok((rec, edges))
    }

    fn ensure_block_loaded(&self, block_index: usize) -> io::Result<()> {
        let key = block_index as u32;
        if self.cache.borrow().contains_key(&key) {
            return Ok(());
        }
        let entry = &self.manifest.blocks[block_index];
        // The fetcher receives the BLOCK INDEX (not a u32 ID derived
        // from state numbering as before). Callers translate it to a
        // URL via the manifest's `block_filename(entry)`.
        let bytes = (self.fetch_block)(block_index as u32)?;
        verify_block_hash(&bytes, &entry.sha256, block_index)?;
        let block = Block::from_gz_bytes(&bytes)?;
        if block.first_state_id != entry.first_state {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "block {block_index} declares first_state_id {} but manifest says {}",
                    block.first_state_id, entry.first_state
                ),
            ));
        }
        self.cache.borrow_mut().insert(key, block);
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
        validate_block_index(&manifest)?;
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
        if rat_slice.len() as u32 > self.manifest.max_indexed_length {
            return None;
        }
        let prefixed = prefix(rat_slice);
        let mut state = 0u32;
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
        let mut state = 0u32;
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
        let mut state = 0u32;
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
        // State 0 (root) served from the manifest -- no I/O.
        if state_id == 0 {
            return Ok((
                root_as_state_rec(&self.manifest.root),
                root_edges(&self.manifest.root),
            ));
        }
        let block_index = self
            .manifest
            .block_index_for_state(state_id)
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "state_id {state_id} not within asset (n_states = {})",
                        self.manifest.n_states
                    ),
                )
            })?;
        let entry = self.manifest.blocks[block_index].clone();
        let local = (state_id - entry.first_state) as usize;
        self.ensure_block_loaded(block_index, fetch).await?;
        let cache = self.cache.borrow();
        let block = cache
            .get(&(block_index as u32))
            .expect("ensure_block_loaded just populated it");
        let rec = block.states[local];
        let edges = block.edges_for(local).to_vec();
        Ok((rec, edges))
    }

    async fn ensure_block_loaded<F, Fut>(&self, block_index: usize, fetch: &F) -> io::Result<()>
    where
        F: Fn(u32) -> Fut,
        Fut: core::future::Future<Output = io::Result<Vec<u8>>>,
    {
        let key = block_index as u32;
        // Quick check; we drop the borrow before awaiting so the
        // RefCell doesn't bridge an await point.
        if self.cache.borrow().contains_key(&key) {
            return Ok(());
        }
        let entry = self.manifest.blocks[block_index].clone();
        let bytes = fetch(block_index as u32).await?;
        verify_block_hash(&bytes, &entry.sha256, block_index)?;
        let block = Block::from_gz_bytes(&bytes)?;
        if block.first_state_id != entry.first_state {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "block {block_index} declares first_state_id {} but manifest says {}",
                    block.first_state_id, entry.first_state
                ),
            ));
        }
        self.cache.borrow_mut().insert(key, block);
        Ok(())
    }
}

// ---- Helpers shared by sync + async readers ----

/// Convert a manifest [`RootState`] (parsed from JSON) into the
/// in-memory [`StateRec`] used by the lookup paths.
fn root_as_state_rec(root: &RootState) -> StateRec {
    StateRec {
        // Edges_offset is meaningless for the root: its edges
        // live in the manifest, not in a block file's edges array.
        edges_offset: 0,
        count: root.count,
        is_accept: root.is_accept,
    }
}

/// Convert a manifest's [`RootEdge`] list into the in-memory
/// [`EdgeRec`] list used by the lookup paths.
fn root_edges(root: &RootState) -> Vec<EdgeRec> {
    root.edges
        .iter()
        .map(|e| EdgeRec {
            label: e.label,
            target: e.target,
        })
        .collect()
}

/// SHA-256-hex of `bytes`, compared to `want`. Used to integrity-
/// check each fetched block file against the hash baked into the
/// manifest.
fn verify_block_hash(bytes: &[u8], want: &str, block_index: usize) -> io::Result<()> {
    let mut hasher = Sha256::new();
    hasher.update(bytes);
    let got = format!("{:x}", hasher.finalize());
    if got != want {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("block {block_index} sha256 mismatch (got {got}, manifest wants {want})"),
        ));
    }
    Ok(())
}

/// Manifest-time validation of the `blocks` index. Catches the
/// malformed cases that would otherwise cause confusing reader
/// errors deep into a query: empty index, out-of-order
/// `first_state`, gaps, n_states inconsistent with last block.
fn validate_block_index(manifest: &BlockManifest) -> io::Result<()> {
    if manifest.blocks.is_empty() {
        if manifest.n_states > 1 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "manifest has n_states={} but no blocks (only root in manifest)",
                    manifest.n_states
                ),
            ));
        }
        return Ok(());
    }
    if manifest.blocks[0].first_state != 1 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "first block must start at state 1 (root is state 0); got {}",
                manifest.blocks[0].first_state
            ),
        ));
    }
    for w in manifest.blocks.windows(2) {
        if w[1].first_state <= w[0].first_state {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "block index not strictly ascending by first_state",
            ));
        }
    }
    let last = &manifest.blocks[manifest.blocks.len() - 1];
    if last.first_state >= manifest.n_states {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "last block first_state={} >= n_states={}",
                last.first_state, manifest.n_states
            ),
        ));
    }
    Ok(())
}

// ---- Writer ----

impl RatDafsa {
    /// Slice this `RatDafsa` into a block_index + per-block files
    /// suitable for serving from a static-asset host. See
    /// [`JSON_SCHEMA_DOC`] for the on-disk layout and
    /// [`LazyRatDafsa`] for the corresponding reader.
    ///
    /// `dir` must already exist. `target_block_bytes` is the
    /// writer's target on-disk size per block: we close a block
    /// once its uncompressed serialised size crosses this threshold,
    /// then gzip + name it by content hash. The last block may be
    /// smaller. Typical value: `DEFAULT_TARGET_BLOCK_BYTES` (1 MiB).
    ///
    /// Asset layout produced:
    ///
    /// ```text
    /// dir/
    ///   block_index.json   (manifest with root + sha256 block index)
    ///   blocks/
    ///     <sha256>.bin     (one per block; filename = SHA-256 of file)
    /// ```
    pub fn write_blocks(&self, dir: &std::path::Path, target_block_bytes: u32) -> io::Result<()> {
        assert!(target_block_bytes >= 1, "target_block_bytes must be >= 1");
        let dafsa = self.inner();
        let n_states = dafsa.raw_n_states() as u32;
        let n_edges = dafsa.raw_n_edges() as u32;
        let counts = dafsa.raw_counts();
        let n_sequences = counts.first().copied().unwrap_or(0) as u32;

        let edges_start = dafsa.raw_edges_start();
        let labels = dafsa.raw_labels();
        let targets = dafsa.raw_targets();
        let is_accept = dafsa.raw_is_accept();

        // Build the root state record for the manifest. Root is
        // state 0; its outgoing edges live at edges_start[0]..edges_start[1].
        let root_edge_lo = edges_start[0] as usize;
        let root_edge_hi = if n_states > 1 {
            edges_start[1] as usize
        } else {
            n_edges as usize
        };
        let root_edges_list: Vec<RootEdge> = (root_edge_lo..root_edge_hi)
            .map(|i| RootEdge {
                label: labels[i],
                target: targets[i],
            })
            .collect();
        let max_indexed_length = root_edges_list
            .iter()
            .map(|e| e.label as u32)
            .max()
            .unwrap_or(0);
        let root = RootState {
            count: counts.first().copied().unwrap_or(0) as u64,
            is_accept: *is_accept.first().unwrap_or(&false),
            edges: root_edges_list,
        };

        let blocks_dir = dir.join("blocks");
        std::fs::create_dir_all(&blocks_dir)?;

        // Walk states 1..n_states (state 0 lives in the manifest)
        // packing them into byte-budgeted blocks. The header and
        // record sizes are known constants so we can predict the
        // exact serialised cost of adding the next state without
        // re-serialising.
        let mut block_index: Vec<BlockEntry> = Vec::new();
        let mut buf_states: Vec<(u32, u64, bool)> = Vec::new();
        let mut buf_edges: Vec<(i8, u32)> = Vec::new();
        let mut block_first_state: u32 = 1;
        let mut block_edge_lo: usize = root_edge_hi;

        for state in 1..n_states {
            let global = state as usize;
            let next_edge_offset = if global + 1 == n_states as usize {
                n_edges as usize
            } else {
                edges_start[global + 1] as usize
            };
            let state_edge_lo = edges_start[global] as usize;
            let state_edge_count = next_edge_offset - state_edge_lo;
            let edges_offset_local = (state_edge_lo - block_edge_lo) as u32;
            // Predict the new size if we append this state.
            let next_bytes = BLOCK_HEADER_BYTES
                + (buf_states.len() + 1) * STATE_RECORD_BYTES
                + (buf_edges.len() + state_edge_count) * EDGE_RECORD_BYTES;
            // If the predicted size crosses the threshold AND we
            // already have at least one state in the buffer, flush
            // first. (Guarantees every block holds >= 1 state, even
            // if one state's edges blow past the threshold alone.)
            if !buf_states.is_empty() && next_bytes as u32 > target_block_bytes {
                flush_block(
                    &blocks_dir,
                    block_first_state,
                    &buf_states,
                    &buf_edges,
                    &mut block_index,
                )?;
                buf_states.clear();
                buf_edges.clear();
                block_first_state = state;
                block_edge_lo = state_edge_lo;
            }
            // Append. edges_offset is block-local (relative to the
            // current block's edges_lo).
            let local_offset = if buf_states.is_empty() {
                0
            } else {
                edges_offset_local
            };
            buf_states.push((local_offset, counts[global] as u64, is_accept[global]));
            for i in state_edge_lo..next_edge_offset {
                buf_edges.push((labels[i], targets[i]));
            }
        }
        if !buf_states.is_empty() {
            flush_block(
                &blocks_dir,
                block_first_state,
                &buf_states,
                &buf_edges,
                &mut block_index,
            )?;
        }

        // Write the manifest last so a reader never sees an index
        // pointing at a block that hasn't been flushed yet.
        let manifest = BlockManifest {
            format: MANIFEST_FORMAT.to_string(),
            version: MANIFEST_VERSION,
            scalar: SCALAR_TAG.to_string(),
            block_format: BLOCK_FORMAT_TAG.to_string(),
            block_version: BLOCK_FORMAT_VERSION,
            target_block_bytes,
            n_states,
            n_edges,
            n_sequences,
            max_indexed_length,
            root,
            blocks: block_index,
            // Default: blocks live next to the manifest under
            // `blocks/`. A publisher who later moves block files to
            // an external host (release CDN, R2, IPFS) rewrites this
            // field in the emitted manifest before deployment.
            block_base_url: None,
        };
        let manifest_path = dir.join("block_index.json");
        let manifest_file = std::fs::File::create(&manifest_path)?;
        serde_json::to_writer_pretty(manifest_file, &manifest)
            .map_err(|e| io::Error::other(format!("write manifest: {e}")))?;
        Ok(())
    }
}

/// Serialise + gzip + hash + write one block, appending its
/// [`BlockEntry`] to `index`.
fn flush_block(
    blocks_dir: &std::path::Path,
    first_state: u32,
    states: &[(u32, u64, bool)],
    edges: &[(i8, u32)],
    index: &mut Vec<BlockEntry>,
) -> io::Result<()> {
    // Build the uncompressed block bytes.
    let n_states = states.len() as u32;
    let n_edges = edges.len() as u32;
    let mut bytes = Vec::with_capacity(
        BLOCK_HEADER_BYTES + states.len() * STATE_RECORD_BYTES + edges.len() * EDGE_RECORD_BYTES,
    );
    bytes.extend_from_slice(BLOCK_MAGIC);
    write_u32(&mut bytes, first_state);
    write_u32(&mut bytes, n_states);
    write_u32(&mut bytes, n_edges);
    for &(off, count, accept) in states {
        write_u32(&mut bytes, off);
        write_u64(&mut bytes, count);
        bytes.push(accept as u8);
        bytes.push(0);
        bytes.push(0);
        bytes.push(0);
    }
    for &(label, target) in edges {
        bytes.push(label as u8);
        bytes.push(0);
        bytes.push(0);
        bytes.push(0);
        write_u32(&mut bytes, target);
    }
    // Gzip.
    let mut gz: Vec<u8> = Vec::new();
    {
        let mut enc = flate2::write::GzEncoder::new(&mut gz, flate2::Compression::default());
        io::Write::write_all(&mut enc, &bytes)?;
        enc.finish()?;
    }
    // SHA-256 of the gzipped bytes -- both the integrity hash AND
    // the filename stem. Cross-edition stability: same bytes ->
    // same name -> same URL.
    let sha256 = {
        let mut h = Sha256::new();
        h.update(&gz);
        format!("{:x}", h.finalize())
    };
    let path = blocks_dir.join(format!("{sha256}.bin"));
    std::fs::write(&path, &gz)?;
    index.push(BlockEntry {
        first_state,
        sha256,
        size: gz.len() as u32,
    });
    Ok(())
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
        // Tiny target so the small test set produces multiple blocks.
        rd.write_blocks(&dir, 128).expect("write_blocks");
        let manifest_path = dir.join("block_index.json");
        assert!(manifest_path.exists(), "manifest must be created");
        let manifest_text = fs::read_to_string(&manifest_path).unwrap();
        let manifest: BlockManifest = serde_json::from_str(&manifest_text).unwrap();
        assert!(!manifest.blocks.is_empty());
        // Non-trivial input should produce more than one block when
        // the target is tiny.
        assert!(
            manifest.blocks.len() >= 2,
            "expected >=2 blocks for {} states, got {}",
            manifest.n_states,
            manifest.blocks.len()
        );

        // Open through a filesystem fetch callback.
        let dir_for_cb = dir.clone();
        let manifest_clone = manifest.clone();
        let fetch = move |block_index: u32| -> io::Result<Vec<u8>> {
            let entry = &manifest_clone.blocks[block_index as usize];
            let name = manifest_clone.block_filename(entry);
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

    /// With `target_block_bytes = 1` every non-root state lives in
    /// its own block, so any non-trivial walk forces fetches across
    /// multiple block boundaries. Guards the cross-block edge-
    /// traversal path.
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
        assert_eq!(manifest.target_block_bytes, 1);
        // Root (state 0) lives in the manifest; one block per
        // remaining state.
        assert_eq!(manifest.blocks.len() as u32, manifest.n_states - 1);

        let dir_for_cb = dir.clone();
        let manifest_clone = manifest.clone();
        let fetch = move |i: u32| {
            let entry = &manifest_clone.blocks[i as usize];
            fs::read(dir_for_cb.join(manifest_clone.block_filename(entry)))
        };
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
    /// for a non-trivial DAFSA with a tiny `target_block_bytes`.
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
        // ~2 states per block at this size given header + record sizes.
        rd.write_blocks(&dir, 48).expect("write_blocks");
        let manifest_text = fs::read_to_string(dir.join("block_index.json")).unwrap();
        let manifest: BlockManifest = serde_json::from_str(&manifest_text).unwrap();
        let total_blocks = manifest.blocks.len();
        assert!(total_blocks >= 2, "test needs at least 2 blocks");

        let dir_for_cb = dir.clone();
        let manifest_clone = manifest.clone();
        let counter = std::sync::Arc::new(std::sync::atomic::AtomicUsize::new(0));
        let counter_for_cb = counter.clone();
        let fetch = move |i: u32| {
            counter_for_cb.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            let entry = &manifest_clone.blocks[i as usize];
            fs::read(dir_for_cb.join(manifest_clone.block_filename(entry)))
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
        // Tiny target so the small set is sliced into multiple blocks.
        original.write_blocks(&dir, 64).expect("write_blocks");
        let manifest_text = fs::read_to_string(dir.join("block_index.json")).unwrap();
        let manifest: BlockManifest = serde_json::from_str(&manifest_text).unwrap();
        // Non-trivial blocking: more than one block, smaller-than-set.
        assert!(manifest.blocks.len() >= 2);

        let dir_for_cb = dir.clone();
        let manifest_clone = manifest.clone();
        let fetch = move |i: u32| {
            let entry = &manifest_clone.blocks[i as usize];
            fs::read(dir_for_cb.join(manifest_clone.block_filename(entry)))
        };
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

    /// Manifests with wrong format / version / scalar / block index
    /// are rejected. Only the manifest is parsed here; no block file
    /// is fetched, so the sha256 strings are placeholders.
    #[test]
    fn manifest_validation() {
        let valid = BlockManifest {
            format: MANIFEST_FORMAT.to_string(),
            version: MANIFEST_VERSION,
            scalar: SCALAR_TAG.to_string(),
            block_format: BLOCK_FORMAT_TAG.to_string(),
            block_version: BLOCK_FORMAT_VERSION,
            target_block_bytes: 1024,
            n_states: 3,
            n_edges: 2,
            n_sequences: 1,
            max_indexed_length: 1,
            root: RootState {
                count: 1,
                is_accept: false,
                edges: vec![RootEdge {
                    label: 1,
                    target: 1,
                }],
            },
            blocks: vec![BlockEntry {
                first_state: 1,
                sha256: "0".repeat(64),
                size: 0,
            }],
            block_base_url: None,
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

        // First block does not start at state 1 (root is 0).
        let mut bad = valid.clone();
        bad.blocks[0].first_state = 2;
        let json = serde_json::to_string(&bad).unwrap();
        assert!(LazyRatDafsa::open(&json, unused).is_err());

        // Non-strictly-ascending block index.
        let mut bad = valid.clone();
        bad.blocks.push(BlockEntry {
            first_state: 1,
            sha256: "1".repeat(64),
            size: 0,
        });
        let json = serde_json::to_string(&bad).unwrap();
        assert!(LazyRatDafsa::open(&json, unused).is_err());

        // Empty index but n_states > 1: orphaned non-root states.
        let mut bad = valid.clone();
        bad.blocks.clear();
        let json = serde_json::to_string(&bad).unwrap();
        assert!(LazyRatDafsa::open(&json, unused).is_err());

        // The valid one must open.
        let json = serde_json::to_string(&valid).unwrap();
        assert!(LazyRatDafsa::open(&json, unused).is_ok());
    }

    /// Upgrade-path immutability: extending an asset with strictly
    /// longer rats (= new length classes added on top) must leave
    /// every full inner block of the smaller asset byte-identical
    /// in the larger asset. Only the manifest's `root` and the
    /// trailing partial block are allowed to change. This is the
    /// concrete property that lets CDN caches stay warm across
    /// asset upgrades: any block file URL that existed before the
    /// upgrade keeps pointing at the same bytes.
    #[test]
    fn upgrade_path_preserves_inner_block_hashes() {
        // Small asset: a handful of short rats with mixed shapes.
        let small_input: Vec<Vec<i8>> = vec![
            vec![-2],
            vec![-1],
            vec![0],
            vec![1],
            vec![2],
            vec![-2, 1, 1],
            vec![-1, 0, 1],
            vec![0, 1, 2],
            vec![1, 2, 3],
            vec![1, 2, 4],
            vec![1, 3, 5],
            vec![2, 2, 2],
        ];
        // Big asset: every small rat + new strictly-longer rats.
        // Length classes 4..6 are entirely new; their sub-DAFSAs
        // are built after the small set's sub-DAFSAs are frozen,
        // so the small set's state numbering survives.
        let mut big_input = small_input.clone();
        big_input.extend(vec![
            vec![1, 2, 3, 4],
            vec![1, 2, 4, 5],
            vec![0, 1, 2, 3],
            vec![1, 2, 3, 4, 5],
            vec![0, 1, 2, 3, 4],
            vec![1, 2, 3, 4, 5, 6],
        ]);

        let small_rd = RatDafsa::from_rats(small_input.iter().map(|v| v.as_slice()));
        let big_rd = RatDafsa::from_rats(big_input.iter().map(|v| v.as_slice()));

        let small_dir = tempdir();
        let big_dir = tempdir();
        // Pick a target small enough that the small asset spans
        // multiple blocks (so we can witness >=1 full block surviving
        // the upgrade), but large enough that not every state is its
        // own block. The exact byte values are not material -- the
        // assertion only relies on small + big using the same target.
        let target = 96u32;
        small_rd.write_blocks(&small_dir, target).expect("small");
        big_rd.write_blocks(&big_dir, target).expect("big");

        let small_manifest: BlockManifest =
            serde_json::from_str(&fs::read_to_string(small_dir.join("block_index.json")).unwrap())
                .unwrap();
        let big_manifest: BlockManifest =
            serde_json::from_str(&fs::read_to_string(big_dir.join("block_index.json")).unwrap())
                .unwrap();

        // Sanity: the small asset has multiple blocks; otherwise the
        // test wouldn't exercise the immutability claim at all.
        assert!(
            small_manifest.blocks.len() >= 2,
            "test needs small asset to have >=2 blocks; got {}",
            small_manifest.blocks.len()
        );
        // Sanity: the big asset is genuinely bigger.
        assert!(big_manifest.n_states > small_manifest.n_states);
        assert!(big_manifest.n_sequences > small_manifest.n_sequences);

        // Every full inner block of the small asset (all but its
        // last block, which may have been only partially full) must
        // appear with the SAME first_state and SAME sha256 in the
        // big asset. That is the cache-immutability guarantee.
        let cutoff = small_manifest.blocks.len() - 1;
        for (i, small_entry) in small_manifest.blocks.iter().take(cutoff).enumerate() {
            let big_entry = &big_manifest.blocks[i];
            assert_eq!(
                big_entry.first_state, small_entry.first_state,
                "block {i}: big.first_state {} != small.first_state {}",
                big_entry.first_state, small_entry.first_state
            );
            assert_eq!(
                big_entry.sha256, small_entry.sha256,
                "block {i} (first_state={}): sha256 diverged after upgrade",
                small_entry.first_state
            );
            // The block file the small asset published is still
            // findable in the big asset's blocks/ dir, byte-identical.
            let small_path = small_dir.join(small_manifest.block_filename(small_entry));
            let big_path = big_dir.join(big_manifest.block_filename(big_entry));
            let small_bytes = fs::read(&small_path).unwrap();
            let big_bytes = fs::read(&big_path).unwrap();
            assert_eq!(
                small_bytes, big_bytes,
                "block {i} (first_state={}): on-disk bytes diverged",
                small_entry.first_state
            );
        }
    }

    /// `block_url` honours `block_base_url`: with the field absent
    /// it returns the relative `blocks/<sha>.bin`; with it set it
    /// returns the joined absolute URL. The integrity check (sha256
    /// in the block index) is unchanged either way -- the field is
    /// purely about resolving the fetch target.
    #[test]
    fn block_url_resolves_with_and_without_base() {
        let entry = BlockEntry {
            first_state: 1,
            sha256: "abcdef0123456789".repeat(4),
            size: 4096,
        };
        // Same fixed length as a real sha256 hex string (64 hex
        // chars) so the URLs are realistic.
        assert_eq!(entry.sha256.len(), 64);

        let mut manifest = BlockManifest {
            format: MANIFEST_FORMAT.to_string(),
            version: MANIFEST_VERSION,
            scalar: SCALAR_TAG.to_string(),
            block_format: BLOCK_FORMAT_TAG.to_string(),
            block_version: BLOCK_FORMAT_VERSION,
            target_block_bytes: 1024,
            n_states: 2,
            n_edges: 0,
            n_sequences: 0,
            max_indexed_length: 0,
            root: RootState {
                count: 0,
                is_accept: false,
                edges: vec![],
            },
            blocks: vec![entry.clone()],
            block_base_url: None,
        };

        // Without block_base_url: relative path matching block_filename.
        assert_eq!(manifest.block_url(&entry), manifest.block_filename(&entry));
        assert!(manifest.block_url(&entry).starts_with("blocks/"));
        assert!(manifest.block_url(&entry).ends_with(".bin"));

        // With block_base_url: joined absolute URL, no `blocks/` prefix.
        manifest.block_base_url = Some(
            "https://github.com/apirogov/tilezz/releases/download/data-zz12-n16-free-v1/"
                .to_string(),
        );
        let url = manifest.block_url(&entry);
        assert!(
            url.starts_with("https://github.com/apirogov/tilezz/releases/download/"),
            "got {url}"
        );
        assert!(url.ends_with(&format!("{}.bin", entry.sha256)), "got {url}");
        assert!(
            !url.contains("/blocks/"),
            "url must skip the `blocks/` prefix when base set; got {url}"
        );

        // The relative on-disk filename is unaffected by block_base_url.
        assert_eq!(
            manifest.block_filename(&entry),
            format!("blocks/{}.bin", entry.sha256)
        );
    }

    /// Round-trip serialisation: a manifest written with
    /// `block_base_url = None` round-trips byte-identically through
    /// serde (no key emitted at all), and a manifest with the field
    /// set survives + parses back to the same value.
    #[test]
    fn block_base_url_serde_roundtrip() {
        let make = |base: Option<&str>| BlockManifest {
            format: MANIFEST_FORMAT.to_string(),
            version: MANIFEST_VERSION,
            scalar: SCALAR_TAG.to_string(),
            block_format: BLOCK_FORMAT_TAG.to_string(),
            block_version: BLOCK_FORMAT_VERSION,
            target_block_bytes: 1024,
            n_states: 2,
            n_edges: 0,
            n_sequences: 0,
            max_indexed_length: 0,
            root: RootState {
                count: 0,
                is_accept: false,
                edges: vec![],
            },
            blocks: vec![BlockEntry {
                first_state: 1,
                sha256: "0".repeat(64),
                size: 0,
            }],
            block_base_url: base.map(str::to_string),
        };

        // None: field omitted from JSON entirely.
        let m = make(None);
        let json = serde_json::to_string(&m).unwrap();
        assert!(
            !json.contains("block_base_url"),
            "field should not serialise when None; got {json}"
        );
        let parsed: BlockManifest = serde_json::from_str(&json).unwrap();
        assert!(parsed.block_base_url.is_none());

        // Some: field round-trips.
        let m = make(Some("https://cdn.example.com/data/"));
        let json = serde_json::to_string(&m).unwrap();
        assert!(json.contains("\"block_base_url\":\"https://cdn.example.com/data/\""));
        let parsed: BlockManifest = serde_json::from_str(&json).unwrap();
        assert_eq!(
            parsed.block_base_url.as_deref(),
            Some("https://cdn.example.com/data/")
        );
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
