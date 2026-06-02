//! Per-thread sort-buffer + run-file writer for Stage 1.
//!
//! Each [`RunWriter`] holds a small in-memory buffer of records.
//! When the buffer fills or the writer is dropped, the buffer is
//! sorted (byte-lex == (length, lex)), local duplicates are folded
//! out, and the bytes are appended to a fresh
//! `out_dir/run_tNN_rMM.bin` file. Stage 2's k-way merge eats those
//! files and produces a fully-sorted unique stream.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::rat_enum::stream::records::encode_record;

/// Default sort-buffer threshold: ~1M records per flush. With ~16
/// bytes per record (ZZ12 n=14: 15 bytes), that's ~16 MB per
/// worker. Tune via [`RunWriter::with_threshold`].
pub const DEFAULT_BUFFER_RECORDS: usize = 1 << 20;

/// One worker's sort-buffer + run-file writer.
pub struct RunWriter {
    /// In-memory sort buffer. Each record is a `Vec<u8>` (length
    /// byte + bias-encoded angles).
    buffer: Vec<Vec<u8>>,
    /// Flush threshold in records.
    threshold: usize,
    /// Directory the run files go in. Created on first flush if
    /// missing.
    out_dir: PathBuf,
    /// Worker identifier (used in the filename so files from
    /// different threads don't collide).
    thread_id: usize,
    /// Monotonic per-writer counter for run filenames.
    run_counter: usize,
    /// Total records emitted across all runs (incl. local
    /// duplicates collapsed during flush).
    total_records_emitted: u64,
}

impl RunWriter {
    /// Create a writer that will flush to `out_dir` with the default
    /// buffer threshold.
    pub fn new(out_dir: impl Into<PathBuf>, thread_id: usize) -> Self {
        Self::with_threshold(out_dir, thread_id, DEFAULT_BUFFER_RECORDS)
    }

    /// Create a writer with a custom buffer-record threshold. Mainly
    /// for tests that want to force multiple flushes on small input.
    pub fn with_threshold(
        out_dir: impl Into<PathBuf>,
        thread_id: usize,
        threshold: usize,
    ) -> Self {
        RunWriter {
            buffer: Vec::with_capacity(threshold.min(1 << 14)),
            threshold,
            out_dir: out_dir.into(),
            thread_id,
            run_counter: 0,
            total_records_emitted: 0,
        }
    }

    /// Record a canonical sequence. Flushes when the buffer is full.
    /// Bias-encodes via [`encode_record`].
    pub fn record(&mut self, canonical: &[i8]) {
        let mut rec = Vec::with_capacity(1 + canonical.len());
        encode_record(canonical, &mut rec);
        self.buffer.push(rec);
        if self.buffer.len() >= self.threshold {
            self.flush().expect("RunWriter flush");
        }
    }

    /// Sort + local-dedup the buffer and append it to a fresh run
    /// file. Returns the path written (for telemetry); also
    /// auto-called on Drop so callers don't have to remember.
    pub fn flush(&mut self) -> std::io::Result<Option<PathBuf>> {
        if self.buffer.is_empty() {
            return Ok(None);
        }
        std::fs::create_dir_all(&self.out_dir)?;
        self.buffer.sort();
        self.buffer.dedup();
        let path = self.out_dir.join(format!(
            "run_t{:02}_r{:06}.bin",
            self.thread_id, self.run_counter
        ));
        let mut writer = BufWriter::new(File::create(&path)?);
        for rec in &self.buffer {
            writer.write_all(rec)?;
        }
        writer.flush()?;
        self.total_records_emitted += self.buffer.len() as u64;
        self.run_counter += 1;
        self.buffer.clear();
        Ok(Some(path))
    }

    /// Total records emitted to disk across all runs (after local
    /// dedup but before the cross-run merge). Useful for telemetry.
    pub fn total_records_emitted(&self) -> u64 {
        self.total_records_emitted
    }

    /// Number of run files this writer has produced so far.
    pub fn run_count(&self) -> usize {
        self.run_counter
    }
}

impl Drop for RunWriter {
    fn drop(&mut self) {
        // Flush on drop so callers don't have to remember.
        // I/O errors during drop are reported via panic only in
        // debug builds; production silently drops the records to
        // avoid panicking in destructors.
        if let Err(e) = self.flush() {
            debug_assert!(false, "RunWriter::drop flush failed: {e}");
            let _ = e;
        }
    }
}

/// Enumerate the run files in `runs_dir` in deterministic order.
/// Used by Stage 2's k-way merge.
pub fn list_run_files(runs_dir: &Path) -> std::io::Result<Vec<PathBuf>> {
    let mut out: Vec<PathBuf> = std::fs::read_dir(runs_dir)?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| {
            p.file_name()
                .and_then(|n| n.to_str())
                .map(|s| s.starts_with("run_t") && s.ends_with(".bin"))
                .unwrap_or(false)
        })
        .collect();
    out.sort();
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rat_enum::stream::records::decode_record;
    use std::io::Read;

    #[test]
    fn writes_sorted_run_with_local_dedup() {
        let dir = tempdir();
        let mut w = RunWriter::with_threshold(&dir, 0, 1000);
        // Duplicates and out-of-order input within one buffer.
        w.record(&[1, 2, 3]);
        w.record(&[1, 2, 3]); // dup
        w.record(&[1]);
        w.record(&[]);
        w.record(&[-1, 0]);
        w.record(&[1, 1]);
        w.flush().expect("flush");

        let files = list_run_files(&dir).unwrap();
        assert_eq!(files.len(), 1);

        let mut buf = Vec::new();
        File::open(&files[0]).unwrap().read_to_end(&mut buf).unwrap();

        // Decode the file's records.
        let mut decoded: Vec<Vec<i8>> = Vec::new();
        let mut rest: &[u8] = &buf;
        while let Some((tail, rec)) = decode_record(rest) {
            decoded.push(rec);
            rest = tail;
        }
        assert!(rest.is_empty(), "trailing bytes in run file");

        // Expected: (length, lex) sorted, no dups.
        assert_eq!(
            decoded,
            vec![vec![], vec![1], vec![-1, 0], vec![1, 1], vec![1, 2, 3],]
        );
    }

    #[test]
    fn auto_flushes_when_buffer_full() {
        let dir = tempdir();
        let mut w = RunWriter::with_threshold(&dir, 7, 3);
        w.record(&[1]);
        w.record(&[2]);
        w.record(&[3]);  // hits threshold -> auto-flush
        w.record(&[4]);
        // Drop will flush the second run.
        drop(w);
        let files = list_run_files(&dir).unwrap();
        assert_eq!(files.len(), 2, "expected 2 run files (auto-flush + drop)");
        assert!(files[0].file_name().unwrap().to_str().unwrap().contains("t07"));
    }

    fn tempdir() -> PathBuf {
        use std::sync::atomic::{AtomicUsize, Ordering};
        static C: AtomicUsize = AtomicUsize::new(0);
        let n = C.fetch_add(1, Ordering::Relaxed);
        let pid = std::process::id();
        let path = std::env::temp_dir().join(format!("rat_enum_runs_test_{pid}_{n}"));
        std::fs::create_dir_all(&path).unwrap();
        path
    }
}
