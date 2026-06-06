//! DAFSA storage stack: compact, queryable, on-disk representation
//! of large sequence sets.
//!
//! Three layers, each adding a concern on top of the previous one:
//!
//! - [`core::Dafsa`] -- the minimum-state DAFSA itself. Builds from a
//!   strictly-sorted sequence iterator via [`core::Dafsa::from_seqs`]
//!   (streaming-friendly: only the active path plus the register of
//!   finalized states are held in memory during construction). Query
//!   surface: `contains`, `get(i)`, `iter`, `walk_prefix`. JSON and
//!   gzipped-JSON round-trip via [`core::Dafsa::write_json`] /
//!   [`core::Dafsa::read_json`].
//!
//! - [`rat::RatDafsa`] -- a `Dafsa` wrapper specialised for "rats"
//!   (canonical cyclotomic angle sequences from `rat_enum`). Stores
//!   each rat with a leading length-prefix byte so the natural lex
//!   traversal yields rats in `(length asc, lex asc)` order, which is
//!   the assigned-index ordering [`rat::RatDafsa::get`] and
//!   [`rat::RatDafsa::index_of`] expose. Wraps the public-facing JSON
//!   read/write in the format discriminator `tilezz-rat-dafsa`.
//!
//! - [`lazy`] -- block-based, lazy-loadable variant of [`rat::RatDafsa`]
//!   for environments that can't (or shouldn't) materialise the whole
//!   automaton up front. The asset is a manifest plus N gzipped block
//!   files; the reader fetches one block on each cache miss.
//!   [`lazy::LazyRatDafsa`] is the synchronous reader (callback-based
//!   fetcher, used by tests and any sync consumer). [`lazy::LazyRatDafsaAsync`]
//!   is the async cousin that mirrors the query surface for WASM /
//!   async-Rust callers; both share the same on-disk format.
//!
//! # Format discriminator strings
//!
//! - `tilezz-dafsa`        -- single-file [`core::Dafsa`] JSON
//! - `tilezz-rat-dafsa`    -- single-file [`rat::RatDafsa`] JSON
//! - `tilezz-rat-dafsa-blocks` -- block manifest for the lazy readers

pub mod core;
pub mod lazy;
pub mod rat;
pub mod rocrate;

use sha2::{Digest, Sha256};

/// Lowercase hex encoding of a byte slice. One canonical spelling for
/// the whole DAFSA stack so the digest-to-string detail lives in a
/// single place (needed since the `sha2 0.11`/`digest 0.11` bump, where
/// the finalized `Output` stopped implementing `LowerHex`).
pub(crate) fn hex_lower(bytes: &[u8]) -> String {
    data_encoding::HEXLOWER.encode(bytes)
}

/// SHA-256 of `bytes` as lowercase hex. Used for block / manifest
/// digests across the DAFSA stack; streaming file hashing lives in
/// `rocrate::sha256_hex` (which shares [`hex_lower`]).
pub(crate) fn sha256_hex_bytes(bytes: &[u8]) -> String {
    let mut h = Sha256::new();
    h.update(bytes);
    hex_lower(&h.finalize())
}

pub use core::Dafsa;
pub use lazy::{DEFAULT_TARGET_BLOCK_BYTES, LazyRatDafsa, LazyRatDafsaAsync};
pub use rat::RatDafsa;
pub use rocrate::{
    AssetParams, ProducedVia, SequenceCounts, rehost_ro_crate, write_archival_extras,
    write_collection_ro_crate, write_ro_crate,
};
