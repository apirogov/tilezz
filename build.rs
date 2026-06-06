//! Build script: injects `TILEZZ_GIT_COMMIT` into the compile-time
//! environment so the CLI `--version` output and the RO-Crate metadata
//! emitted by `rat_enum` can record the exact source commit that built
//! each binary / asset.
//!
//! Paranoia: if the working tree has ANY uncommitted change (modified
//! tracked file or untracked non-ignored file), the commit is suffixed
//! `-dirty`. A pristine build therefore self-reports a bare 40-hex
//! commit; anything else is a loud signal that the binary was NOT built
//! from a clean checkout of that commit -- exactly the provenance
//! mix-up this guards against. Build from a `git worktree` pinned at the
//! release tag to get a clean stamp.
//!
//! Falls back to "unknown" outside a git checkout (source tarballs,
//! vendored builds). Cargo crate-metadata (name, version, repository,
//! etc.) is already exposed via the `CARGO_PKG_*` env vars.

use std::process::Command;

fn git(args: &[&str]) -> Option<String> {
    Command::new("git")
        .args(args)
        .output()
        .ok()
        .filter(|o| o.status.success())
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .map(|s| s.trim().to_string())
}

fn main() {
    let commit = git(&["rev-parse", "HEAD"]).unwrap_or_else(|| "unknown".to_string());
    // `--porcelain` lists modified-tracked + untracked-non-ignored
    // entries; empty output means a pristine tree. (.gitignored paths
    // such as target/ never appear, so they don't trip the flag.)
    let dirty = git(&["status", "--porcelain"])
        .map(|s| !s.is_empty())
        .unwrap_or(false);
    let stamp = if commit != "unknown" && dirty {
        format!("{commit}-dirty")
    } else {
        commit
    };
    println!("cargo:rustc-env=TILEZZ_GIT_COMMIT={stamp}");
    // Re-run when HEAD/refs move or the index changes, so the dirty
    // flag and commit stay current without rebuilding on every cargo
    // invocation.
    println!("cargo:rerun-if-changed=.git/HEAD");
    println!("cargo:rerun-if-changed=.git/refs/heads");
    println!("cargo:rerun-if-changed=.git/index");
}
