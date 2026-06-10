FROM ubuntu:24.04

# Avoid interactive prompts during package install
ENV DEBIAN_FRONTEND=noninteractive

# Basic dev tools + a curated set of CLI utilities. All in one
# apt-install so we land in a single image layer (separate RUNs
# cost disk + break cache more often). Categories:
#   build + version control + editors + file manager + multiplexer:
#       curl git build-essential pkg-config libssl-dev
#       vim nnn tmux rsync sudo libfontconfig1-dev ca-certificates
#   linker + WASM optimiser (`mold` cuts Rust link time ~5x;
#   `binaryen` provides `wasm-opt` so wasm-pack doesn't have to
#   fetch it on first run):
#       mold binaryen
#   task runner used by the project's justfile:
#       just
#   timing + structured-data wrangling:
#       time jq
#   modern grep / find / cat (faster than the GNU classics):
#       ripgrep fd-find bat
#   directory + process + open-file inspection:
#       tree htop lsof
#   man pages + shell completion (the slim image ships without
#   either; bash-completion picks up `cargo`, `git`, `gh`, ...
#   if /etc/bash_completion is sourced at shell startup):
#       man-db manpages bash-completion
#   file inspection + compressor coverage:
#       less file unzip xz-utils zstd wget
#   network + low-level debugging:
#       iputils-ping dnsutils netcat-openbsd strace gdb
#   remote access -- `git` uses this transparently for `git@github.com:...`
#   URLs; without it, fetches/pushes to ssh remotes silently fail
#   with "cannot run ssh: No such file or directory":
#       openssh-client
#
# Not installed but considered: `perf` is awkward in a container --
# Ubuntu 24.04 has no standalone `linux-perf` package, and the
# kernel-coupled `linux-tools-*-generic` chain only matches the
# build host's kernel (not the container's runtime host). Use
# `samply record <bin>` (cargo-installed below) instead -- it
# doesn't need perf, doesn't need root, doesn't care about the
# host kernel, and emits Firefox-profiler-compatible traces.
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl git build-essential pkg-config libssl-dev \
    vim nnn tmux rsync sudo libfontconfig1-dev ca-certificates \
    mold binaryen \
    just \
    time jq \
    ripgrep fd-find bat \
    tree htop lsof \
    man-db manpages bash-completion \
    less file unzip xz-utils zstd wget \
    iputils-ping dnsutils netcat-openbsd strace gdb \
    openssh-client \
    && rm -rf /var/lib/apt/lists/*

# GitHub CLI -- not in the default ubuntu repos, so we add its
# official apt source. Useful for `gh pr create`, `gh release ...`,
# auth-aware `gh api ...` calls, etc.
RUN curl -fsSL https://cli.github.com/packages/githubcli-archive-keyring.gpg \
        -o /usr/share/keyrings/githubcli-archive-keyring.gpg \
    && chmod a+r /usr/share/keyrings/githubcli-archive-keyring.gpg \
    && echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/githubcli-archive-keyring.gpg] https://cli.github.com/packages stable main" \
        > /etc/apt/sources.list.d/github-cli.list \
    && apt-get update \
    && apt-get install -y --no-install-recommends gh \
    && rm -rf /var/lib/apt/lists/*

# Install Rust via rustup (system-wide layout so it survives uid remapping)
ENV RUSTUP_HOME=/usr/local/rustup \
    CARGO_HOME=/usr/local/cargo \
    PATH=/usr/local/cargo/bin:$PATH
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \
    | sh -s -- -y --no-modify-path --default-toolchain stable \
    && rustup component add clippy rustfmt rust-analyzer rust-src \
    && rustup target add wasm32-unknown-unknown \
    && chmod -R a+w /usr/local/rustup /usr/local/cargo

# Use mold for native linking. gcc 12.1+ supports `-fuse-ld=mold`
# directly so we don't need to switch the C compiler to clang. The
# target triple is derived from rustc's own host (works on x86_64 and
# aarch64) so this isn't pinned to one arch -- wasm32 uses LLVM's own
# linker (lld) and isn't affected. Written under CARGO_HOME so every
# user (including the cargo install below + the non-root runtime user)
# picks it up automatically.
RUN triple="$(rustc -Vv | sed -n 's/^host: //p')" \
    && printf '%s\n' \
        "[target.${triple}]" \
        'rustflags = ["-C", "link-arg=-fuse-ld=mold"]' \
        > /usr/local/cargo/config.toml

# Cargo subcommands that round out the Rust dev loop:
#   wasm-pack       - WASM frontend build (the explorer)
#   cargo-nextest   - parallel test runner with better UX
#   cargo-watch     - referenced by the project's `just watch`
#   cargo-edit      - `cargo add` / `cargo rm` / `cargo upgrade`
#   samply          - Firefox-profiler-compatible sampler; no perf,
#                     no root required, container-friendly
# This step runs as root and fills /usr/local/cargo/registry with
# root-owned files, so the chmod must run again afterwards -- the
# one in the rustup step above predates these artifacts. Without it,
# the non-root runtime user gets "Permission denied" the first time
# a build needs a crate that wasn't cached at image-build time.
RUN cargo install --locked \
        wasm-pack cargo-nextest cargo-watch cargo-edit \
        samply \
    && chmod -R a+w /usr/local/cargo

# Install zellij. The release asset is arch-specific, so map the Debian
# arch to zellij's triple instead of hardcoding x86_64 (lets the image
# build on arm64 hosts too).
RUN arch="$(dpkg --print-architecture)" \
    && case "$arch" in \
        amd64) ztriple=x86_64-unknown-linux-musl ;; \
        arm64) ztriple=aarch64-unknown-linux-musl ;; \
        *) echo "unsupported arch for zellij: $arch" >&2; exit 1 ;; \
    esac \
    && curl -fsSL "https://github.com/zellij-org/zellij/releases/latest/download/zellij-${ztriple}.tar.gz" \
    | tar -xz -C /usr/local/bin

# Install uv to /usr/local/bin
COPY --from=ghcr.io/astral-sh/uv:0.11.13 /uv /uvx /bin/

# Install Python 3.13 and set it as the default. uv places
# `python` / `python3` / `python3.13` symlinks on PATH so the
# container behaves like any other Python system; you only need
# to type `uv` for uv-specific operations (`uv pip install ...`,
# `uvx ...`, `uv tool install ...`).
RUN uv python install --default --force 3.13

# Python CLI tools that should be available globally. `uv tool install`
# parks each in its own isolated venv under `~/.local/share/uv/tools/`
# and drops a launcher symlink in `~/.local/bin/` -- the
# `chmod -R a+rX /root/.local` step below makes both readable by the
# non-root runtime user. Add new tools here as needed (mypy, pre-commit,
# ipython, ...) -- one per line keeps diffs surgical.
RUN uv tool install ruff

# Open up /root/.local so the non-root runtime user can run the
# uv-installed tools (python, ruff) whose launchers live in
# /root/.local/bin.
RUN chmod 755 /root && chmod -R a+rX /root/.local

# Pre-create ~/.config and ~/.cache owned by ubuntu. rust-sandbox.sh
# bind-mounts host dirs *under* these (the zellij config under .config,
# the Mathlib olean cache under .cache/mathlib); if a mount target's
# parent is missing, the docker daemon creates it as root:root (the -u
# flag only applies to the container process, not mount setup), which
# then blocks every other tool (nnn, lake, ...) from creating its own
# subdir there with "Permission denied".
RUN mkdir -p /home/ubuntu/.config /home/ubuntu/.cache \
    && chown ubuntu:ubuntu /home/ubuntu/.config /home/ubuntu/.cache

# Neither Claude Code nor the Lean toolchain is baked into the image;
# both are installed on first `start` into host-bind-mounted home dirs
# (see rust-sandbox.sh), so a single install plus its auto-updates
# survive container recreation:
#   ~/.local   Claude Code -- per-user install under .local/share with a
#              launcher symlinked into .local/bin. Baking a copy under
#              /root would just get duplicated (233 MB) into the runtime
#              user's home on first run anyway.
#   ~/.elan    elan + Lean toolchains -- installed via elan-init; each
#              project's pinned toolchain auto-installs here on first use.
# .elan/bin and .local/bin go first so lean/lake/claude resolve to the
# mounts; /root/.local/bin stays on PATH for the uv tools above.
ENV PATH="/home/ubuntu/.elan/bin:/home/ubuntu/.local/bin:/root/.local/bin:$PATH"

ENV SHELL=/bin/bash

CMD ["bash"]
