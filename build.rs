//! Build script: injects `TILEZZ_GIT_COMMIT` into the compile-time
//! environment so RO-Crate metadata emitted by `rat_enum --mode
//! dafsa-blocks` can record the source commit that built each asset.
//!
//! Falls back to "unknown" outside a git checkout (e.g. source
//! tarballs, vendored builds). Cargo crate-metadata (name, version,
//! repository, etc.) is already exposed via the `CARGO_PKG_*` env
//! vars at compile time -- no build script needed for those.

use std::process::Command;

fn main() {
    let commit = Command::new("git")
        .args(["rev-parse", "HEAD"])
        .output()
        .ok()
        .filter(|o| o.status.success())
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .map(|s| s.trim().to_string())
        .unwrap_or_else(|| "unknown".to_string());
    println!("cargo:rustc-env=TILEZZ_GIT_COMMIT={commit}");
    // Re-run only when HEAD or its target ref moves; avoids spurious
    // rebuilds on every cargo invocation.
    println!("cargo:rerun-if-changed=.git/HEAD");
    println!("cargo:rerun-if-changed=.git/refs/heads");
}
