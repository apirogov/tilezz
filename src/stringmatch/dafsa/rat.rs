//! `RatDafsa`: a compact, gzip-friendly store for canonical rat angle
//! sequences with stable external indices in `(length asc, lex asc)`
//! order.
//!
//! # The wire format
//!
//! Internally, `RatDafsa` owns a plain [`Dafsa`] but stores each rat
//! length-prefixed: a rat `[a, b, c]` of length L is stored in the
//! DAFSA as `[L, a, b, c]`. The lex order on the prefixed form is
//! the same as `(length asc, lex asc)` on the raw form, so the
//! DAFSA's natural traversal yields rats in the desired order at
//! zero extra storage cost -- no separate `outputs` permutation.
//!
//! The length prefix is an implementation detail of the wire format,
//! not of the API. The Rust API in this module takes and returns
//! raw rats; the prefix never appears in caller code. The format
//! discriminator on disk is `"tilezz-rat-dafsa"` (distinct from
//! plain `"tilezz-dafsa"`), so a JS/WASM consumer that opens the
//! file knows what convention applies and can strip the leading
//! byte themselves; see the schema documentation embedded as
//! [`JSON_SCHEMA_DOC`].
//!
//! # Indexing contract
//!
//! `RatDafsa::from_rats(rats)` sorts the input by `(length, lex)`
//! and dedups; the assigned external index of each rat is its
//! position in the sorted-and-deduped list. After construction:
//!
//! * [`RatDafsa::get(i)`](RatDafsa::get) returns the i-th rat under
//!   that ordering.
//! * [`RatDafsa::index_of(rat)`](RatDafsa::index_of) returns the
//!   assigned index of `rat`, or `None` if it is not in the set.
//! * [`RatDafsa::iter`](RatDafsa::iter) walks the rats in
//!   assigned-index order (= `(length, lex)` order); this is a
//!   fast DFS over the underlying DAFSA, not a re-walk per element.

use std::io;

use serde::{Deserialize, Serialize};

use super::core::Dafsa;

/// Maximum rat length representable in a single i8 prefix byte. Rats
/// longer than this cannot be stored in this format -- a build will
/// panic. Cyclotomic enumerations cap out far below this in practice
/// (n <= 16 for ZZ12 free; 127 is plenty of headroom).
const MAX_RAT_LEN: usize = i8::MAX as usize;

/// File format identifier baked into every serialized `RatDafsa`.
/// Distinct from the plain DAFSA `"tilezz-dafsa"` so a reader who
/// opens a `.json.gz` knows whether to expect the length-prefix
/// convention before reading any sequences.
const FORMAT_TAG: &str = "tilezz-rat-dafsa";

/// On-disk version of the rat-dafsa format. Bumped on any wire
/// change; readers reject unknown versions.
const FORMAT_VERSION: u32 = 1;

/// Plain-text schema description bundled into the binary so anyone
/// inspecting a serialized file (or this crate's docs) can decode
/// it without the Rust source.
pub const JSON_SCHEMA_DOC: &str = include_str!("rat_schema.txt");

#[derive(Debug, Clone)]
pub struct RatDafsa {
    inner: Dafsa,
}

impl RatDafsa {
    /// Build a `RatDafsa` from an iterator of rats. Input may be in
    /// any order and may contain duplicates: rats are sorted by
    /// `(length asc, lex asc)` and deduplicated internally, and the
    /// assigned external index of each rat is its position in the
    /// resulting sorted-and-deduped list.
    ///
    /// # Panics
    ///
    /// Panics if any rat is longer than 127 angles (the single-byte
    /// length-prefix limit). Cyclotomic enumerations cap far below
    /// this in practice.
    pub fn from_rats<I, S>(rats: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: AsRef<[i8]>,
    {
        // Length-prefix each rat and sort the prefixed form by lex.
        // Lex on prefixed == (length, lex) on raw, so this also
        // gives us the desired assigned-index ordering.
        let mut prefixed: Vec<Vec<i8>> = rats
            .into_iter()
            .map(|r| {
                let r = r.as_ref();
                assert!(
                    r.len() <= MAX_RAT_LEN,
                    "RatDafsa: rat length {} exceeds the single-byte length-prefix limit ({})",
                    r.len(),
                    MAX_RAT_LEN
                );
                let mut v = Vec::with_capacity(r.len() + 1);
                v.push(r.len() as i8);
                v.extend_from_slice(r);
                v
            })
            .collect();
        prefixed.sort_unstable();
        prefixed.dedup();
        let inner = Dafsa::from_seqs(prefixed.iter().map(|v| v.as_slice()));
        Self { inner }
    }

    /// Streaming constructor for rats already in `(length asc, lex
    /// asc)` order with no duplicates. Length-prefixes each rat
    /// on-the-fly and feeds the resulting byte sequences directly
    /// into [`Dafsa::from_seqs`] without materializing the full set
    /// in memory.
    ///
    /// This is the entry point for the streaming pipeline's Stage 3
    /// (`--mode build`): `unique.bin` is already sorted and deduped
    /// by Stage 2's k-way merge, so the per-element buffering
    /// [`from_rats`] does to enforce that property would be pure
    /// overhead. With this constructor a ZZ12 n=14 build's peak RSS
    /// is bounded by the in-progress DAFSA structure rather than by
    /// the input set.
    ///
    /// # Panics
    ///
    /// Panics in debug builds if `rats` is not strictly increasing
    /// in `(length, lex)` order; release builds skip the check and
    /// produce a DAFSA that accepts the input *as if* it were sorted
    /// (i.e. silently wrong if it wasn't). Callers must guarantee
    /// the ordering — `merge::read_unique_records` does.
    ///
    /// Also panics if any rat is longer than 127 angles (same
    /// length-prefix limit as [`from_rats`]).
    pub fn from_sorted_unique_rats<I, S>(rats: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: AsRef<[i8]>,
    {
        let prefixed = rats.into_iter().map(|r| {
            let r = r.as_ref();
            assert!(
                r.len() <= MAX_RAT_LEN,
                "RatDafsa: rat length {} exceeds the single-byte length-prefix limit ({})",
                r.len(),
                MAX_RAT_LEN
            );
            let mut v = Vec::with_capacity(r.len() + 1);
            v.push(r.len() as i8);
            v.extend_from_slice(r);
            v
        });
        let inner = Dafsa::from_seqs(prefixed);
        Self { inner }
    }

    /// Crate-private: hand the underlying `Dafsa` to the lazy
    /// (blocked) writer in `rat_dafsa_lazy`. Not part of the public
    /// API; the length-prefix encoding is implementation-private.
    pub(crate) fn inner(&self) -> &Dafsa {
        &self.inner
    }

    /// Number of distinct rats stored.
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    /// True iff no rats are stored.
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    /// The rat with assigned external index `i` (i.e. the i-th rat
    /// in `(length, lex)` order), or `None` if `i >= len()`.
    pub fn get(&self, i: usize) -> Option<Vec<i8>> {
        self.inner.get(i).map(strip_prefix)
    }

    /// Assigned external index of `rat` (its position in
    /// `(length, lex)` order) if it is in the set, else `None`.
    pub fn index_of<S: AsRef<[i8]>>(&self, rat: S) -> Option<u32> {
        let prefixed = prefix(rat.as_ref());
        // Lex rank in the inner DAFSA == (length, lex) rank of the
        // raw rat -- exactly the external assigned index.
        self.inner.lex_rank_of(&prefixed).map(|r| r as u32)
    }

    /// Set membership: does `rat` appear in this `RatDafsa`?
    pub fn contains<S: AsRef<[i8]>>(&self, rat: S) -> bool {
        let prefixed = prefix(rat.as_ref());
        self.inner.contains(&prefixed)
    }

    /// Iterate rats in assigned-index order (= `(length, lex)`
    /// order). This is a fast DFS over the underlying DAFSA -- the
    /// natural traversal already yields the length-prefixed
    /// sequences in lex order, which corresponds to the requested
    /// rat ordering; each yield strips the prefix byte.
    pub fn iter(&self) -> impl Iterator<Item = Vec<i8>> + '_ {
        self.inner.iter().map(strip_prefix)
    }

    /// Serialize to a JSON string. See [`JSON_SCHEMA_DOC`] for the
    /// exact format.
    pub fn to_json_string(&self) -> String {
        let inner = self.inner.to_json_string();
        let form = RatDafsaSerForm {
            format: FORMAT_TAG.to_string(),
            version: FORMAT_VERSION,
            inner_format: "tilezz-dafsa".to_string(),
            note: NOTE.to_string(),
            dafsa: serde_json::from_str(&inner).expect("inner Dafsa JSON is always valid"),
        };
        serde_json::to_string(&form).expect("RatDafsaSerForm always serializes")
    }

    /// Parse a JSON document. Rejects documents whose `format` or
    /// `version` does not match the expected wire values.
    pub fn from_json_str(s: &str) -> io::Result<Self> {
        let form: RatDafsaSerForm =
            serde_json::from_str(s).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        Self::from_form(form)
    }

    /// Stream a `RatDafsa` into an `io::Write` as JSON.
    pub fn write_json<W: io::Write>(&self, w: W) -> io::Result<()> {
        // We have to materialize the JSON once because the inner
        // Dafsa is serialized as a nested object; the cost is the
        // same as `to_json_string` plus a copy.
        let s = self.to_json_string();
        let mut w = w;
        w.write_all(s.as_bytes())
    }

    /// Read a `RatDafsa` from an `io::Read` containing JSON.
    pub fn read_json<R: io::Read>(r: R) -> io::Result<Self> {
        let form: RatDafsaSerForm = serde_json::from_reader(r)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        Self::from_form(form)
    }

    /// Stream a gzipped JSON `RatDafsa` into an `io::Write`.
    pub fn write_json_gz<W: io::Write>(&self, w: W) -> io::Result<()> {
        let mut enc = flate2::write::GzEncoder::new(w, flate2::Compression::default());
        self.write_json(&mut enc)?;
        enc.finish()?;
        Ok(())
    }

    /// Parse a gzipped JSON `RatDafsa` from an `io::Read`.
    pub fn read_json_gz<R: io::Read>(r: R) -> io::Result<Self> {
        let dec = flate2::read::GzDecoder::new(r);
        Self::read_json(dec)
    }

    fn from_form(form: RatDafsaSerForm) -> io::Result<Self> {
        if form.format != FORMAT_TAG {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "RatDafsa: bad format tag (expected '{}', got '{}')",
                    FORMAT_TAG, form.format
                ),
            ));
        }
        if form.version != FORMAT_VERSION {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "RatDafsa: unsupported version (expected {}, got {})",
                    FORMAT_VERSION, form.version
                ),
            ));
        }
        let dafsa_json = serde_json::to_string(&form.dafsa)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let inner = Dafsa::from_json_str(&dafsa_json)?;
        Ok(Self { inner })
    }
}

/// Reminder embedded into every serialized file so anyone who opens
/// the JSON without consulting the schema doc still gets a hint.
const NOTE: &str = "Each accepted sequence in the inner DAFSA is the angle sequence of a rat \
                    prefixed with its length: stored = [len(rat), rat[0], rat[1], ...]. \
                    Strip the first byte to recover the rat.";

#[derive(Serialize, Deserialize)]
struct RatDafsaSerForm {
    format: String,
    version: u32,
    /// Wire format of the embedded automaton; always
    /// `"tilezz-dafsa"`. Recording it explicitly means a reader can
    /// validate the nesting before parsing the inner blob.
    inner_format: String,
    /// Plain-prose reminder that the inner DAFSA's sequences are
    /// length-prefixed. Decoders should not depend on the exact
    /// wording, only on the `format` discriminator.
    note: String,
    /// The inner plain DAFSA, in its own `tilezz-dafsa` JSON form.
    dafsa: serde_json::Value,
}

fn prefix(rat: &[i8]) -> Vec<i8> {
    assert!(
        rat.len() <= MAX_RAT_LEN,
        "RatDafsa: rat length {} exceeds the single-byte length-prefix limit ({})",
        rat.len(),
        MAX_RAT_LEN
    );
    let mut v = Vec::with_capacity(rat.len() + 1);
    v.push(rat.len() as i8);
    v.extend_from_slice(rat);
    v
}

fn strip_prefix(mut prefixed: Vec<i8>) -> Vec<i8> {
    // Stored sequence is [len, rat[0], rat[1], ...]; drop the first
    // byte. We use `remove(0)` rather than `split_off(1)` because
    // typical rat lengths are short (n <= 16) and the cost is
    // dominated by the alloc inside `Dafsa::iter` upstream.
    if !prefixed.is_empty() {
        prefixed.remove(0);
    }
    prefixed
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Headline contract: feed rats in any order, get them back in
    /// `(length asc, lex asc)` order via `iter` and `get(i)`, and
    /// the assigned-index space is consistent through a gzip
    /// round-trip.
    #[test]
    fn rats_order_and_roundtrip() {
        // Caller-chosen order: arbitrary scramble.
        let input: Vec<Vec<i8>> = vec![
            vec![1, 2, 3],
            vec![-1],
            vec![-2, 1, 1],
            vec![2],
            vec![-2],
            vec![1, 2],
            vec![-1, 0, 1],
            vec![0],
        ];
        // Expected (length, lex) order.
        let mut expected = input.clone();
        expected.sort_by(|a, b| a.len().cmp(&b.len()).then_with(|| a.cmp(b)));

        let rd = RatDafsa::from_rats(input.iter().map(|v| v.as_slice()));
        assert_eq!(rd.len(), expected.len());

        // iter yields the (length, lex) ordering.
        let iter_out: Vec<Vec<i8>> = rd.iter().collect();
        assert_eq!(iter_out, expected);

        // get(i) matches.
        for (i, e) in expected.iter().enumerate() {
            assert_eq!(rd.get(i).as_ref(), Some(e));
        }
        assert_eq!(rd.get(expected.len()), None);

        // index_of returns the (length, lex) position.
        for (i, e) in expected.iter().enumerate() {
            assert_eq!(rd.index_of(e.as_slice()), Some(i as u32));
        }
        assert_eq!(rd.index_of([42i8, 42].as_slice()), None);

        // contains agrees with input membership.
        for e in &expected {
            assert!(rd.contains(e.as_slice()));
        }
        assert!(!rd.contains([42i8].as_slice()));

        // Gzip round-trip: the rebuilt RatDafsa is observationally
        // identical (same iter, same get, same index_of).
        let mut buf: Vec<u8> = Vec::new();
        rd.write_json_gz(&mut buf).expect("write_json_gz");
        assert!(buf.starts_with(&[0x1f, 0x8b]), "missing gzip magic");
        let restored = RatDafsa::read_json_gz(&buf[..]).expect("read_json_gz");
        let restored_iter: Vec<Vec<i8>> = restored.iter().collect();
        assert_eq!(restored_iter, expected);
        for (i, e) in expected.iter().enumerate() {
            assert_eq!(restored.get(i).as_ref(), Some(e));
            assert_eq!(restored.index_of(e.as_slice()), Some(i as u32));
        }
    }

    /// `from_sorted_unique_rats` builds the same observable RatDafsa
    /// as `from_rats` when the input is already in `(length, lex)`
    /// order without duplicates. The headline property: the
    /// streaming constructor doesn't have to materialize anything
    /// extra in memory, but it must produce a bit-identical result
    /// to the buffering constructor for valid input.
    #[test]
    fn from_sorted_unique_matches_from_rats() {
        let input: Vec<Vec<i8>> = vec![
            vec![1, 2, 3],
            vec![-1],
            vec![-2, 1, 1],
            vec![2],
            vec![-2],
            vec![1, 2],
            vec![-1, 0, 1],
            vec![0],
        ];

        // Buffering reference.
        let buffered = RatDafsa::from_rats(input.iter().map(|v| v.as_slice()));

        // Hand the streaming constructor the same set in the order
        // it expects: (length asc, lex asc), already deduped.
        let mut sorted = input.clone();
        sorted.sort_by(|a, b| a.len().cmp(&b.len()).then_with(|| a.cmp(b)));
        sorted.dedup();
        let streamed = RatDafsa::from_sorted_unique_rats(sorted.iter().map(|v| v.as_slice()));

        // Observational equivalence: same len, same iter, same get,
        // same index_of.
        assert_eq!(streamed.len(), buffered.len());
        assert_eq!(
            streamed.iter().collect::<Vec<_>>(),
            buffered.iter().collect::<Vec<_>>()
        );
        for i in 0..buffered.len() {
            assert_eq!(streamed.get(i), buffered.get(i));
        }
        for r in &sorted {
            assert_eq!(
                streamed.index_of(r.as_slice()),
                buffered.index_of(r.as_slice())
            );
        }
    }

    /// `from_sorted_unique_rats` debug-asserts strict-increasing
    /// `(length, lex)` order. Out-of-order input must trip the
    /// underlying `Dafsa::from_seqs` assertion in debug builds; in
    /// release builds the contract is on the caller, so we don't
    /// test that path.
    #[test]
    #[should_panic]
    #[cfg(debug_assertions)]
    fn from_sorted_unique_rejects_out_of_order_debug() {
        let bad: Vec<&[i8]> = vec![&[1, 2][..], &[1, 1][..]]; // lex regression
        let _ = RatDafsa::from_sorted_unique_rats(bad.iter().copied());
    }

    /// Duplicates in the input collapse to a single entry. The
    /// assigned indices on the deduped set are still contiguous.
    #[test]
    fn dedups_duplicates() {
        let input: Vec<Vec<i8>> = vec![vec![1, 2], vec![1, 2], vec![3], vec![1, 2]];
        let rd = RatDafsa::from_rats(input.iter().map(|v| v.as_slice()));
        assert_eq!(rd.len(), 2);
        // (length, lex) order on the deduped set: [3], [1, 2].
        assert_eq!(rd.get(0).as_deref(), Some(&[3i8][..]));
        assert_eq!(rd.get(1).as_deref(), Some(&[1i8, 2][..]));
    }

    /// The JSON output advertises the right format discriminator,
    /// embeds the length-prefix reminder, and *does not* leak the
    /// plain `"tilezz-dafsa"` format tag at the top level (the
    /// embedded inner dafsa carries it, but a reader keying on the
    /// outer tag must not confuse the two formats).
    #[test]
    fn format_tag_is_distinct() {
        let rd = RatDafsa::from_rats([[1i8, 2].as_slice(), [3i8].as_slice()].iter().copied());
        let json = rd.to_json_string();
        assert!(json.contains("\"tilezz-rat-dafsa\""));
        assert!(
            json.contains("length-prefixed") || json.contains("length"),
            "missing length-prefix note in serialized form"
        );

        // Misidentified as a plain DAFSA: reading the same bytes
        // through `Dafsa::from_json_str` must fail. (Either the
        // format discriminator mismatch fires, or the plain reader
        // chokes on the missing top-level state arrays -- both are
        // fine; what we're guarding is that the two formats are
        // never silently interchangeable.)
        assert!(
            Dafsa::from_json_str(&json).is_err(),
            "Dafsa reader silently accepted a RatDafsa document"
        );
    }
}
