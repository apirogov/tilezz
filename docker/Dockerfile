FROM ubuntu:26.04

# Avoid interactive prompts during package install
ENV DEBIAN_FRONTEND=noninteractive

# Basic dev tools + a curated set of CLI utilities, in one apt-install so we
# land in a single image layer. These stay baked: they are the "minimum in
# the image" -- system tooling, plus the compilers/linkers/native libs that
# Rust and Python builds need at runtime. Language toolchains themselves are
# NOT baked; they install on first start into bind mounts (see provision.sh).
#   build + native deps Rust crates link against:
#       curl git build-essential pkg-config libssl-dev libfontconfig1-dev
#   version control + editors + file manager + multiplexer + misc:
#       vim nnn tmux rsync sudo ca-certificates
#   linker + WASM optimiser (mold cuts Rust link time ~5x; binaryen provides
#   wasm-opt so wasm-pack does not fetch it on first run):
#       mold binaryen
#   task runner used by project justfiles:
#       just
#   timing + structured-data wrangling:
#       time jq
#   modern grep / find / cat:
#       ripgrep fd-find bat
#   directory + process + open-file inspection:
#       tree htop lsof
#   man pages + shell completion:
#       man-db manpages bash-completion
#   file inspection + compressor coverage:
#       less file unzip xz-utils zstd wget
#   network + low-level debugging:
#       iputils-ping dnsutils netcat-openbsd strace gdb
#   remote access (git uses this for git@github.com: URLs):
#       openssh-client
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
    && rm -rf /var/lib/apt/lists/* \
    # Ubuntu renames these two binaries to dodge name clashes with other
    # packages (fd-find -> fdfind, bat -> batcat). Symlink the conventional
    # names so `fd` and `bat` resolve as users expect.
    && ln -s /usr/bin/fdfind /usr/local/bin/fd \
    && ln -s /usr/bin/batcat /usr/local/bin/bat

# GitHub CLI -- not in the default ubuntu repos, so add its official apt
# source. Useful for gh pr create, gh release ..., auth-aware gh api ....
RUN curl -fsSL https://cli.github.com/packages/githubcli-archive-keyring.gpg \
        -o /usr/share/keyrings/githubcli-archive-keyring.gpg \
    && chmod a+r /usr/share/keyrings/githubcli-archive-keyring.gpg \
    && echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/githubcli-archive-keyring.gpg] https://cli.github.com/packages stable main" \
        > /etc/apt/sources.list.d/github-cli.list \
    && apt-get update \
    && apt-get install -y --no-install-recommends gh \
    && rm -rf /var/lib/apt/lists/*

# Install zellij. The release asset is arch-specific, so map the Debian arch
# to zellij's triple instead of hardcoding x86_64 (lets the image build on
# arm64 hosts too). A prebuilt binary, not a toolchain, so it stays baked.
RUN arch="$(dpkg --print-architecture)" \
    && case "$arch" in \
        amd64) ztriple=x86_64-unknown-linux-musl ;; \
        arm64) ztriple=aarch64-unknown-linux-musl ;; \
        *) echo "unsupported arch for zellij: $arch" >&2; exit 1 ;; \
    esac \
    && curl -fsSL "https://github.com/zellij-org/zellij/releases/latest/download/zellij-${ztriple}.tar.gz" \
    | tar -xz -C /usr/local/bin

# Install uv to /bin. uv is an installer (like rustup-via-curl), so it stays
# baked; the Python runtime + CLI tools it manages are provisioned on first
# start into the bind-mounted ~/.local (see provision.sh).
COPY --from=ghcr.io/astral-sh/uv:0.11.13 /uv /uvx /bin/

# Pre-create ~/.config and ~/.cache owned by ubuntu. The launcher bind-mounts
# host dirs *under* these (gh/git/opencode config under .config, the Mathlib
# olean cache under .cache); if a mount target's parent is missing, the docker
# daemon creates it root:root (the -u flag only applies to the container
# process, not mount setup), which then blocks other tools from creating
# subdirs there with "Permission denied".
RUN mkdir -p /home/ubuntu/.config /home/ubuntu/.cache \
    && chown ubuntu:ubuntu /home/ubuntu/.config /home/ubuntu/.cache

# Nothing language-specific is baked. On first `start` the launcher pipes
# provision.sh into the container, which installs into host-bind-mounted home
# dirs so installs + auto-updates survive container recreation:
#   ~/.local    Claude Code + uv-managed Python + uv tools (ruff);
#               also opencode auth + data (~/.local/share/opencode)
#   ~/.opencode opencode agent binary
#   ~/.elan     elan + Lean toolchains
#   ~/.rustup   rustup + Rust toolchains
#   ~/.cargo    cargo bin (rust tools) + registry cache + mold config
# These bin dirs go first so the mounted toolchains win over anything else.
ENV PATH="/home/ubuntu/.elan/bin:/home/ubuntu/.cargo/bin:/home/ubuntu/.local/bin:/home/ubuntu/.opencode/bin:$PATH"

ENV SHELL=/bin/bash

CMD ["bash"]
