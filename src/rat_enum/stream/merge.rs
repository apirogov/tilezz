//! Stage 2 of the streaming pipeline: k-way merge over the per-thread
//! run files, dedup, and produce a single sorted `unique.bin` plus a
//! `certificate.json` recording the BLAKE3 hash + headline counts.
//!
//! Why a separate stage: Stage 1 writes locally-sorted, locally-deduped
//! runs but the same rat can appear in different runs (different
//! workers reach it via different seeds). Stage 2 does the global
//! dedup with bounded memory by streaming a k-way merge.

use std::collections::BinaryHeap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::rat_enum::stream::runs::list_run_files;

/// Filename written under `out_dir`. Holds the deduplicated, sorted
/// concatenation of every record from `out_dir/runs/*.bin`.
pub const UNIQUE_FILENAME: &str = "unique.bin";

/// Filename written under `out_dir`. Records the BLAKE3 hash of
/// `unique.bin` plus headline counts + the run-file inventory --
/// enough to verify a Stage 1 + Stage 2 result is reproducible
/// before kicking off the (expensive) Stage 3 DAFSA build.
pub const CERTIFICATE_FILENAME: &str = "certificate.json";

/// I/O buffer for each run-file reader. Records are small (typically
/// <= 32 bytes) so a few-KB buffer amortizes syscalls without
/// holding hundreds of MB across many runs.
const READER_BUFFER_BYTES: usize = 64 * 1024;

/// Schema version for `certificate.json`. Bump when the on-disk
/// shape changes incompatibly.
pub const CERTIFICATE_SCHEMA_VERSION: u32 = 1;

/// On-disk layout of `certificate.json`. Serialized via serde_json
/// with pretty formatting so users can `cat` it for sanity checks.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Certificate {
    /// Bumped on incompatible schema changes.
    pub schema_version: u32,
    /// Cyclotomic ring (e.g. 12 for ZZ12).
    pub ring: u8,
    /// Maximum perimeter walked by the DFS.
    pub max_steps: usize,
    /// Direction step (typically 1; > 1 restricts to a subset).
    pub step: i8,
    /// Whether the source enumeration was dihedral-canonical.
    pub dihedral: bool,
    /// Number of run files consumed by the merge.
    pub run_files: usize,
    /// Total records observed across all run files (duplicates
    /// included).
    pub total_input_records: u64,
    /// Number of unique records written to `unique.bin`.
    pub unique_records: u64,
    /// Byte length of `unique.bin`.
    pub unique_bytes: u64,
    /// BLAKE3 hash of `unique.bin` (hex, 64 chars).
    pub unique_blake3: String,
}

/// Min-heap entry for the k-way merge. Compared by `key` (the
/// current record bytes), then by `source_index` to keep the
/// ordering total and stable.
struct HeapEntry {
    key: Vec<u8>,
    source_index: usize,
}

// BinaryHeap is a *max* heap; flip ordering so the smallest record
// (byte-lex) and smallest source index come out first.
impl Ord for HeapEntry {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        other
            .key
            .cmp(&self.key)
            .then_with(|| other.source_index.cmp(&self.source_index))
    }
}
impl PartialOrd for HeapEntry {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl PartialEq for HeapEntry {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == std::cmp::Ordering::Equal
    }
}
impl Eq for HeapEntry {}

/// Streaming reader for one run file. Buffered; reads one record at
/// a time via the length-prefix.
struct RunReader {
    inner: BufReader<File>,
    buf: Vec<u8>,
}

impl RunReader {
    fn open(path: &Path) -> std::io::Result<Self> {
        let f = File::open(path)?;
        Ok(RunReader {
            inner: BufReader::with_capacity(READER_BUFFER_BYTES, f),
            buf: Vec::with_capacity(64),
        })
    }

    /// Reads the next record into the internal buffer and returns a
    /// fresh `Vec<u8>` view of it, or `None` at EOF. Returns the full
    /// `[length, angles...]` byte sequence so byte-lex comparison
    /// works against other reader outputs.
    fn next_record(&mut self) -> std::io::Result<Option<Vec<u8>>> {
        let mut len_buf = [0u8; 1];
        match self.inner.read_exact(&mut len_buf) {
            Ok(()) => {}
            Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => return Ok(None),
            Err(e) => return Err(e),
        }
        let len = len_buf[0] as usize;
        self.buf.clear();
        self.buf.push(len_buf[0]);
        self.buf.resize(1 + len, 0);
        self.inner.read_exact(&mut self.buf[1..1 + len])?;
        Ok(Some(self.buf.clone()))
    }
}

/// K-way merge `out_dir/runs/*.bin` into `out_dir/unique.bin`,
/// dedup as we go. Writes `out_dir/certificate.json` with the
/// BLAKE3 hash of unique.bin and the headline counts.
///
/// `ring`, `max_steps`, `step`, `dihedral` are recorded in the
/// certificate verbatim -- they are not re-derived from the runs
/// (the runs themselves carry no metadata).
#[allow(clippy::too_many_arguments)]
pub fn merge_runs(
    out_dir: &Path,
    ring: u8,
    max_steps: usize,
    step: i8,
    dihedral: bool,
) -> std::io::Result<Certificate> {
    let runs_dir = out_dir.join(crate::rat_enum::stream::enumerate::RUNS_SUBDIR);
    let run_files = list_run_files(&runs_dir)?;
    if run_files.is_empty() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::NotFound,
            format!("no run files under {}", runs_dir.display()),
        ));
    }
    println!(
        "merge: consuming {} run file(s) from {}",
        run_files.len(),
        runs_dir.display()
    );

    let mut readers: Vec<RunReader> = run_files
        .iter()
        .map(|p| RunReader::open(p))
        .collect::<std::io::Result<_>>()?;

    let mut heap = BinaryHeap::with_capacity(readers.len());
    for (i, r) in readers.iter_mut().enumerate() {
        if let Some(rec) = r.next_record()? {
            heap.push(HeapEntry {
                key: rec,
                source_index: i,
            });
        }
    }

    let unique_path = out_dir.join(UNIQUE_FILENAME);
    let unique_file = File::create(&unique_path)?;
    let mut writer = BufWriter::with_capacity(1 << 20, unique_file);
    let mut hasher = blake3::Hasher::new();

    let mut total_in: u64 = 0;
    let mut unique_out: u64 = 0;
    let mut unique_bytes: u64 = 0;
    let mut last: Option<Vec<u8>> = None;

    while let Some(top) = heap.pop() {
        total_in += 1;
        let key = top.key;
        let src = top.source_index;

        let is_new = match last.as_ref() {
            Some(prev) => prev != &key,
            None => true,
        };
        if is_new {
            writer.write_all(&key)?;
            hasher.update(&key);
            unique_out += 1;
            unique_bytes += key.len() as u64;
            last = Some(key);
        }

        if let Some(next_rec) = readers[src].next_record()? {
            heap.push(HeapEntry {
                key: next_rec,
                source_index: src,
            });
        }
    }
    writer.flush()?;
    drop(writer);

    let unique_blake3 = hasher.finalize().to_hex().to_string();

    let cert = Certificate {
        schema_version: CERTIFICATE_SCHEMA_VERSION,
        ring,
        max_steps,
        step,
        dihedral,
        run_files: run_files.len(),
        total_input_records: total_in,
        unique_records: unique_out,
        unique_bytes,
        unique_blake3,
    };
    let cert_path = out_dir.join(CERTIFICATE_FILENAME);
    let cert_file = File::create(&cert_path)?;
    serde_json::to_writer_pretty(BufWriter::new(cert_file), &cert)?;

    println!(
        "merge: {} input records -> {} unique ({} bytes); blake3={}",
        cert.total_input_records, cert.unique_records, cert.unique_bytes, cert.unique_blake3
    );
    println!("merge: wrote {} + {}", unique_path.display(), cert_path.display());

    Ok(cert)
}

/// Iterate the records in a `unique.bin` file. Yields each canonical
/// sequence as a `Vec<i8>` in (length, lex) order. Useful for Stage
/// 3 (DAFSA build) and for cross-validation tests.
pub fn read_unique_records(path: &Path) -> std::io::Result<UniqueRecordIter> {
    let f = File::open(path)?;
    Ok(UniqueRecordIter {
        inner: BufReader::with_capacity(READER_BUFFER_BYTES, f),
    })
}

/// Streaming iterator over a `unique.bin` file. See [`read_unique_records`].
#[derive(Debug)]
pub struct UniqueRecordIter {
    inner: BufReader<File>,
}

impl Iterator for UniqueRecordIter {
    type Item = std::io::Result<Vec<i8>>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut len_buf = [0u8; 1];
        match self.inner.read_exact(&mut len_buf) {
            Ok(()) => {}
            Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => return None,
            Err(e) => return Some(Err(e)),
        }
        let len = len_buf[0] as usize;
        let mut bytes = vec![0u8; len];
        if let Err(e) = self.inner.read_exact(&mut bytes) {
            return Some(Err(e));
        }
        Some(Ok(bytes
            .into_iter()
            .map(|b| (b as i16 - crate::rat_enum::stream::records::ANGLE_BIAS) as i8)
            .collect()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rat_enum::stream::runs::RunWriter;
    use std::path::PathBuf;
    use std::sync::atomic::{AtomicUsize, Ordering};

    fn tempdir() -> PathBuf {
        static C: AtomicUsize = AtomicUsize::new(0);
        let n = C.fetch_add(1, Ordering::Relaxed);
        let pid = std::process::id();
        let path = std::env::temp_dir().join(format!("rat_enum_merge_test_{pid}_{n}"));
        std::fs::create_dir_all(&path).unwrap();
        path
    }

    #[test]
    fn merge_dedupes_across_runs_and_preserves_order() {
        let dir = tempdir();
        let runs_dir = dir.join(super::super::enumerate::RUNS_SUBDIR);
        std::fs::create_dir_all(&runs_dir).unwrap();

        // Two workers that both find the rat [1,2,3] -- merge must
        // emit it exactly once.
        {
            let mut w = RunWriter::with_threshold(&runs_dir, 1, 1000);
            w.record(&[]);
            w.record(&[1, 2, 3]);
            w.record(&[2]);
        }
        {
            let mut w = RunWriter::with_threshold(&runs_dir, 2, 1000);
            w.record(&[1, 2, 3]); // duplicate of worker 1's rat
            w.record(&[-1, 0]);
            w.record(&[1, 1]);
        }

        let cert = merge_runs(&dir, 12, 9, 1, false).unwrap();
        assert_eq!(cert.run_files, 2);
        assert_eq!(cert.total_input_records, 6);
        assert_eq!(cert.unique_records, 5);

        let recs: Vec<Vec<i8>> = read_unique_records(&dir.join(UNIQUE_FILENAME))
            .unwrap()
            .map(|r| r.unwrap())
            .collect();
        assert_eq!(
            recs,
            vec![vec![], vec![2], vec![-1, 0], vec![1, 1], vec![1, 2, 3],]
        );

        // Certificate hash is the BLAKE3 of unique.bin's bytes -- compute
        // it again and check.
        let bytes = std::fs::read(dir.join(UNIQUE_FILENAME)).unwrap();
        assert_eq!(blake3::hash(&bytes).to_hex().to_string(), cert.unique_blake3);
    }

    #[test]
    fn merge_errors_on_empty_runs_dir() {
        let dir = tempdir();
        let runs_dir = dir.join(super::super::enumerate::RUNS_SUBDIR);
        std::fs::create_dir_all(&runs_dir).unwrap();
        let err = merge_runs(&dir, 12, 9, 1, false).unwrap_err();
        assert_eq!(err.kind(), std::io::ErrorKind::NotFound);
    }

    /// Calling `merge_runs` twice on the same `runs/` directory must
    /// be a no-op for outputs: `unique.bin` and `certificate.json`
    /// get cleanly overwritten, BLAKE3 stays the same, no stale
    /// bytes leak through. Catches a class of "merge appends rather
    /// than truncates" bugs.
    #[test]
    fn merge_idempotent_on_rerun() {
        let dir = tempdir();
        let runs_dir = dir.join(super::super::enumerate::RUNS_SUBDIR);
        std::fs::create_dir_all(&runs_dir).unwrap();
        {
            let mut w = RunWriter::with_threshold(&runs_dir, 1, 1000);
            w.record(&[1, 2, 3]);
            w.record(&[2]);
            w.record(&[]);
        }
        {
            let mut w = RunWriter::with_threshold(&runs_dir, 2, 1000);
            w.record(&[1, 2, 3]); // cross-run duplicate
            w.record(&[-1, 0]);
        }

        let cert1 = merge_runs(&dir, 12, 9, 1, false).unwrap();
        let unique_1 = std::fs::read(dir.join(UNIQUE_FILENAME)).unwrap();
        let cert_json_1 = std::fs::read(dir.join(CERTIFICATE_FILENAME)).unwrap();

        let cert2 = merge_runs(&dir, 12, 9, 1, false).unwrap();
        let unique_2 = std::fs::read(dir.join(UNIQUE_FILENAME)).unwrap();
        let cert_json_2 = std::fs::read(dir.join(CERTIFICATE_FILENAME)).unwrap();

        assert_eq!(cert1.unique_blake3, cert2.unique_blake3);
        assert_eq!(cert1.unique_records, cert2.unique_records);
        assert_eq!(cert1.total_input_records, cert2.total_input_records);
        assert_eq!(unique_1, unique_2, "unique.bin not byte-identical across merges");
        assert_eq!(cert_json_1, cert_json_2, "certificate.json not byte-identical across merges");
    }
}
