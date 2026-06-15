#!/usr/bin/env bash
# First-start provisioner for dev-sandbox. The launcher pipes this into the
# container on every start (docker exec -i ... bash -s < provision.sh). Every
# block is idempotent: a guard short-circuits once the tool/toolchain is
# present, so the first start does the work and later starts are a no-op.
# Editing this file needs no image rebuild -- it is read host-side at start.
set -euo pipefail
trap 'echo "[provision] FAILED near line $LINENO" >&2' ERR

echo "[provision] checking agents + toolchains..."

# Claude Code -- per-user install under ~/.local/share with a launcher in
# ~/.local/bin (both on the persisted ~/.local mount).
command -v claude >/dev/null \
  || curl -fsSL https://claude.ai/install.sh | bash

# opencode -- second agent. Its installer hardcodes ~/.opencode/bin (persisted
# via its own mount); auth + data land in ~/.local/share/opencode (rides the
# ~/.local mount), so one `opencode auth login` inside the sandbox persists.
command -v opencode >/dev/null \
  || curl -fsSL https://opencode.ai/install | bash

# elan + Lean -- stable toolchain into the persisted ~/.elan; a project's own
# lean-toolchain auto-installs there on use.
command -v elan >/dev/null \
  || curl -fsSL https://elan.lean-lang.org/elan-init.sh \
     | sh -s -- -y --default-toolchain stable --no-modify-path

# Rust -- rustup + the stable toolchain into the persisted ~/.rustup / ~/.cargo
# (rustup's defaults inside HOME, both bind-mounted).
command -v rustup >/dev/null \
  || curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \
     | sh -s -- -y --no-modify-path --default-toolchain stable
# Components + wasm target (idempotent: no-op when already present).
rustup component add clippy rustfmt rust-analyzer rust-src
rustup target add wasm32-unknown-unknown

# cargo-binstall -- fetches prebuilt release binaries so the cargo tools below
# install in seconds instead of compiling from source.
command -v cargo-binstall >/dev/null \
  || curl -L --proto '=https' --tlsv1.2 -sSf \
       https://raw.githubusercontent.com/cargo-bins/cargo-binstall/main/install-from-binstall-release.sh \
     | bash

# Cargo dev tools -- gated by a sentinel so the set installs once. binstall
# pulls prebuilt binaries (fast). Add tools to this line as needed.
if [ ! -f "$HOME/.cargo/.devsandbox-tools" ]; then
  cargo binstall -y wasm-pack cargo-nextest cargo-watch cargo-edit samply
  touch "$HOME/.cargo/.devsandbox-tools"
fi

# mold linker config -- the same -fuse-ld=mold rustflags the image used to
# bake, relocated into the cargo mount. gcc 12.1+ supports -fuse-ld=mold
# directly. The target triple comes from rustc's own host so it is not
# arch-pinned (wasm32 uses LLVM's lld and is unaffected).
if [ ! -f "$HOME/.cargo/config.toml" ]; then
  triple="$(rustc -Vv | sed -n 's/^host: //p')"
  printf '%s\n' \
    "[target.${triple}]" \
    'rustflags = ["-C", "link-arg=-fuse-ld=mold"]' \
    > "$HOME/.cargo/config.toml"
fi

# Python -- uv-managed 3.13 set as default (uv skips if already installed);
# the pythons + the ruff tool venv land under ~/.local (persisted). uv places
# python / python3 / python3.13 symlinks on PATH so the container behaves like
# any other Python system.
uv python install --default 3.13
command -v ruff >/dev/null || uv tool install ruff

echo "[provision] done."
