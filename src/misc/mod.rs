//! Miscellaneous utilities that don't fit one of the major layers.
//!
//! Currently the lone occupant is [`profile`]: an RAII helper for
//! optional `pprof` flamegraph recording from the CLI binaries.

pub mod profile;
