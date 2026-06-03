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
# directly so we don't need to switch the C compiler to clang. Only
# the x86_64 target is configured -- wasm32 uses LLVM's own linker
# (lld) and isn't affected. Written under CARGO_HOME so every user
# (including the cargo install below + the non-root runtime user)
# picks it up automatically.
RUN printf '%s\n' \
        '[target.x86_64-unknown-linux-gnu]' \
        'rustflags = ["-C", "link-arg=-fuse-ld=mold"]' \
        > /usr/local/cargo/config.toml

# Cargo subcommands that round out the Rust dev loop:
#   wasm-pack       - WASM frontend build (the explorer)
#   cargo-nextest   - parallel test runner with better UX
#   cargo-watch     - referenced by the project's `just watch`
#   cargo-edit      - `cargo add` / `cargo rm` / `cargo upgrade`
#   samply          - Firefox-profiler-compatible sampler; no perf,
#                     no root required, container-friendly
# Build artifacts land under CARGO_HOME which we already chmod a+w
# above, so the runtime non-root user picks them up automatically.
RUN cargo install --locked \
        wasm-pack cargo-nextest cargo-watch cargo-edit \
        samply

# Install zellij
RUN curl -fsSL https://github.com/zellij-org/zellij/releases/latest/download/zellij-x86_64-unknown-linux-musl.tar.gz \
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

# Pre-create ~/.config owned by ubuntu. rust-sandbox.sh bind-mounts
# the host's zellij config to /home/ubuntu/.config/zellij; if the
# mount target's parent is missing, the docker daemon creates it as
# root:root (the -u flag only applies to the container process, not
# mount setup), which then blocks every other tool (nnn, ...) from
# creating its own config dir with "Permission denied".
RUN mkdir -p /home/ubuntu/.config && chown ubuntu:ubuntu /home/ubuntu/.config

# Claude Code is NOT baked into the image: the native launcher manages
# per-user installs under $HOME/.local/share/claude and symlinks
# $HOME/.local/bin/claude, so an image copy under /root just gets
# duplicated (233 MB) into the runtime user's home on first run anyway.
# Instead, /home/ubuntu/.local is bind-mounted from the host (see
# rust-sandbox.sh) and the install happens once on first `start`;
# auto-updates land in the mount and survive container recreation.
# /home/ubuntu/.local/bin comes first so claude resolves there;
# /root/.local/bin stays on PATH for the uv tools above.
ENV PATH="/home/ubuntu/.local/bin:/root/.local/bin:$PATH"

ENV SHELL=/bin/bash

CMD ["bash"]
