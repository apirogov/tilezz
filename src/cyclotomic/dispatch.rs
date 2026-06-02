//! Runtime ring dispatch: pick a typed cyclotomic ring from a `u8`.
//!
//! Many call sites need to turn a runtime-known `ring: u8` (from the
//! CLI, the WASM API, a JSON manifest, etc.) into an invocation of a
//! generic function instantiated with the corresponding `ZZn` type.
//! The natural Rust expression for this is a long `match` on the ring
//! number, but with 10 supported rings and 8+ dispatch sites in the
//! crate, the matches were ~80 lines of near-identical boilerplate
//! that all needed editing every time a new ring landed.
//!
//! [`dispatch_ring!`] collapses that to one line per call site. Inside
//! each match arm, the macro shadows two names:
//!
//! - `ZZ` -- a `type` alias pointing at the ring-specific
//!   `ZZ4`/`ZZ6`/.../`ZZ60`.
//! - `PHI` -- the ring's integer-basis storage width (`phi(n)`), as a
//!   `const usize`. Useful for generic functions that take `PHI` as a
//!   const generic parameter (e.g. `extract::<ZZ, PHI, _>(...)`).
//!
//! # Usage
//!
//! ```ignore
//! use tilezz::dispatch_ring;
//!
//! // Simple call with a single type-generic argument:
//! let stats = dispatch_ring!(ring, enumerate_dispatch::<ZZ>(
//!     max_steps, step, n_threads, dihedral, paranoid,
//! ));
//!
//! // Both ZZ (type) and PHI (const) are in scope in the body:
//! let unit_count = dispatch_ring!(ring,
//!     extract::<ZZ, PHI, _>(|x| x.int_coeffs())
//! );
//!
//! // Storing a fn pointer rather than calling immediately:
//! let f: fn(usize) -> Vec<i8> = dispatch_ring!(ring, collect_closure_keys::<ZZ>);
//! ```
//!
//! # Adding a new ring
//!
//! Add one line to the macro body and one `define_integral_zz!`
//! invocation in `rings.rs`. Every existing call site picks it up
//! automatically -- no per-site edits required.

/// Dispatch a runtime `ring: u8` value to a generic body that wants
/// a concrete `ZZn` type. The body is evaluated in a context where:
///
/// - `ZZ` is a `type` alias for the matching `ZZn` ring,
/// - `PHI` is a `const usize` equal to that ring's `phi(n)`
///   (storage width of the integer-basis representation).
///
/// Two forms:
///
/// ```ignore
/// // Panic on unknown ring (typical for CLI / internal call sites):
/// dispatch_ring!(ring, body)
///
/// // Fallback expression on unknown ring (typical for WASM /
/// // user-facing entry points that prefer a structured error):
/// dispatch_ring!(ring, body, else fallback)
/// ```
///
/// The bare form expands to the `else` form with a panic fallback,
/// so both share the same per-ring match table; adding a ring touches
/// exactly one place.
#[macro_export]
macro_rules! dispatch_ring {
    ($ring:expr, $body:expr) => {
        $crate::dispatch_ring!(
            $ring,
            $body,
            else panic!(
                "dispatch_ring: invalid ring {} (expected one of \
                 4, 6, 8, 10, 12, 16, 20, 24, 32, 60)",
                $ring
            )
        )
    };
    ($ring:expr, $body:expr, else $fallback:expr) => {{
        match $ring {
            4 => {
                #[allow(dead_code)] type ZZ = $crate::cyclotomic::ZZ4;
                #[allow(dead_code)] const PHI: usize = 2;
                $body
            }
            6 => {
                #[allow(dead_code)] type ZZ = $crate::cyclotomic::ZZ6;
                #[allow(dead_code)] const PHI: usize = 2;
                $body
            }
            8 => {
                #[allow(dead_code)] type ZZ = $crate::cyclotomic::ZZ8;
                #[allow(dead_code)] const PHI: usize = 4;
                $body
            }
            10 => {
                #[allow(dead_code)] type ZZ = $crate::cyclotomic::ZZ10;
                #[allow(dead_code)] const PHI: usize = 4;
                $body
            }
            12 => {
                #[allow(dead_code)] type ZZ = $crate::cyclotomic::ZZ12;
                #[allow(dead_code)] const PHI: usize = 4;
                $body
            }
            14 => {
                #[allow(dead_code)] type ZZ = $crate::cyclotomic::ZZ14;
                #[allow(dead_code)] const PHI: usize = 6;
                $body
            }
            16 => {
                #[allow(dead_code)] type ZZ = $crate::cyclotomic::ZZ16;
                #[allow(dead_code)] const PHI: usize = 8;
                $body
            }
            18 => {
                #[allow(dead_code)] type ZZ = $crate::cyclotomic::ZZ18;
                #[allow(dead_code)] const PHI: usize = 6;
                $body
            }
            20 => {
                #[allow(dead_code)] type ZZ = $crate::cyclotomic::ZZ20;
                #[allow(dead_code)] const PHI: usize = 8;
                $body
            }
            24 => {
                #[allow(dead_code)] type ZZ = $crate::cyclotomic::ZZ24;
                #[allow(dead_code)] const PHI: usize = 8;
                $body
            }
            32 => {
                #[allow(dead_code)] type ZZ = $crate::cyclotomic::ZZ32;
                #[allow(dead_code)] const PHI: usize = 16;
                $body
            }
            60 => {
                #[allow(dead_code)] type ZZ = $crate::cyclotomic::ZZ60;
                #[allow(dead_code)] const PHI: usize = 16;
                $body
            }
            _ => $fallback,
        }
    }};
}

/// The set of ring numbers `dispatch_ring!` knows about. Useful for
/// CLI argument validation, doc generation, and tests that iterate
/// over every supported ring.
pub const SUPPORTED_RINGS: [u8; 12] = [4, 6, 8, 10, 12, 14, 16, 18, 20, 24, 32, 60];

#[cfg(test)]
mod tests {
    use crate::cyclotomic::SymNum;

    /// Dispatching to a `ZZ::PHI` lookup must agree with the
    /// `SUPPORTED_RINGS` table -- this protects against
    /// dispatch_ring's hard-coded PHI values drifting from the
    /// per-ring `PHI` consts in `rings.rs`.
    #[test]
    fn phi_matches_per_ring_constant() {
        for &ring in &super::SUPPORTED_RINGS {
            let dispatched_phi: usize = crate::dispatch_ring!(ring, ZZ::PHI);
            assert_eq!(
                dispatched_phi,
                expected_phi(ring),
                "dispatch_ring PHI for ZZ{ring} disagrees with the table",
            );
        }
    }

    /// One canonical-element round-trip per supported ring, dispatched
    /// through the macro. Smoke test that the macro can dispatch a
    /// generic operation across every ring without bizarre lifetime
    /// or type-inference issues.
    #[test]
    fn one_is_one_for_every_ring() {
        use num_traits::One;
        for &ring in &super::SUPPORTED_RINGS {
            let c64 = crate::dispatch_ring!(ring, <ZZ as One>::one().complex64());
            assert!(
                (c64.re - 1.0).abs() < 1e-12 && c64.im.abs() < 1e-12,
                "ZZ{ring}::one().complex64() = {c64:?}",
            );
        }
    }

    fn expected_phi(ring: u8) -> usize {
        match ring {
            4 | 6 => 2,
            8 | 10 | 12 => 4,
            14 | 18 => 6,
            16 | 20 | 24 => 8,
            32 | 60 => 16,
            _ => panic!("unreachable"),
        }
    }
}
