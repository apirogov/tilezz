//! Cross-cutting utilities used by binaries and various library
//! modules.
//!
//! Currently the lone occupant is [`profile`]: an RAII helper for
//! optional `pprof` flamegraph recording from the CLI binaries.

pub mod profile;
