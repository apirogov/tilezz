//! Streaming pipeline for very large enumerations.
//!
//! Three independently invocable stages; all share `-o <output_dir>`:
//!
//! ```text
//!  Stage 1  --mode stream  -> out_dir/runs/run_tNN_rMM.bin (sorted runs)
//!  Stage 2  --mode merge   -> out_dir/unique.bin + certificate.json
//!  Stage 3  --mode build   -> out_dir/dafsa/      (blocked RatDafsa asset)
//! ```
//!
//! What this buys us: peak memory bounded by the per-thread sort
//! buffer (Stage 1) and then by the in-progress DAFSA structure
//! (Stage 3), never by the final set size. ZZ12 n=14+ becomes
//! feasible on a commodity workstation without `--mode list-seeds`
//! orchestration. Stage 3 is the streaming equivalent of
//! `--mode dafsa-blocks` -- both write the same on-disk format, but
//! Stage 3 doesn't materialize the rat set in a `Vec<Vec<i8>>` first.
//!
//! Records are **length-prefixed and bias-encoded**:
//!   * byte 0: length L (0..=255)
//!   * bytes 1..=L: each angle byte = (i8 angle + 128) as u8
//!
//! Bias-encoding makes byte-lex order equal `(length asc, integer-lex
//! asc)`, which is the same order the DAFSA expects. So sorted runs
//! and the merged `unique.bin` are already in the right order for
//! Stage 3 without any in-memory unsort pass.
//!
//! See [`records`] for encode/decode helpers, [`runs`] for the
//! per-thread sort-buffer + flush writer, [`enumerate`] for the
//! `stream_enum_parallel` entry point that wires it all together.

pub mod enumerate;
pub mod merge;
pub mod records;
pub mod runs;

pub use enumerate::{stream_enum_dispatch, stream_enum_parallel};
pub use merge::{
    merge_runs, read_unique_records, Certificate, CERTIFICATE_FILENAME, UNIQUE_FILENAME,
};
pub use records::{decode_record, encode_record, ANGLE_BIAS};
pub use runs::RunWriter;
