#![allow(clippy::doc_overindented_list_items)]
//! Cyclotomic rings.
//!
//! Every `ZZn` ring lives on the **integer-basis** representation in
//! `rings.rs`: storage `[i64; phi(n)]` over the power basis
//! `{1, zeta, ..., zeta^(phi(n)-1)}`, with closed-form exact sign
//! extraction (no f64 fallback). The shared engine lives in
//! `integral_basis.rs` (free functions + `define_integral_zz!` macro).
//!
//! File layout:
//!
//! * `numtraits.rs`       -- foundational numeric traits (`IntRing`,
//!                           `IntField`, `ZSigned`, `InnerIntType`).
//! * `traits.rs`          -- cyclotomic-ring trait surface (`SymNum`,
//!                           `Ccw`, `Conj`, `OneImag`, `ReImSign`,
//!                           `IntersectUnitSegments`, `WithinRadius`,
//!                           `Units`, the `IsRing` hierarchy, the
//!                           `HasZZk` subring markers).
//! * `sign.rs`            -- exact-sign helpers `signum_sum_sqrt_expr_*`.
//! * `integral_basis.rs`  -- the generic engine: free helpers and the
//!                           `define_integral_zz!` macro.
//! * `rings.rs`           -- one `define_integral_zz!` invocation per
//!                           ring, sqrt constants, ring-specific tests.
//! * `constants.rs`       -- `sqrt2/3/5/6` and `zz_units_sum`.
//! * `geometry.rs`        -- `intersect` segment predicate.
//! * `linalg.rs`          -- small linear-algebra utilities.

pub(crate) mod numtraits;

mod sign;

pub(crate) mod traits;

mod rings;

pub mod dispatch;

pub mod integral_basis;

// --------

pub mod constants;

pub mod linalg;

pub mod geometry;

pub use num_traits::{One, Pow, Zero};

pub use numtraits::{IntField, IntRing};

pub use traits::{
    Ccw, CellFloor, Conj, HasZZ4, HasZZ6, HasZZ8, HasZZ10, HasZZ12, IntersectUnitSegments, IsRing,
    IsZZ4, OneImag, ReImSign, SymNum, Units, WithinRadius, ZZComplex,
};

pub use rings::{ZZ4, ZZ6, ZZ8, ZZ10, ZZ12, ZZ14, ZZ16, ZZ18, ZZ20, ZZ24, ZZ32, ZZ60};
