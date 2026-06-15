#!/usr/bin/env bash
set -euo pipefail

# Optional container-name suffix so the same image can serve several
# workspaces concurrently: `dev-sandbox -n NAME <cmd>` targets the container
# `dev-sandbox-NAME`. With no -n, the single shared `dev-sandbox` container is
# used. All persisted state (toolchains, caches, agent installs, auth) is
# shared across containers regardless of name; only the container name and the
# mounted workspace differ.
SANDBOX_NAME=""
while [ $# -gt 0 ]; do
  case "${1:-}" in
    -n|--name)
      if [ $# -lt 2 ]; then echo "error: $1 needs a NAME" >&2; exit 1; fi
      SANDBOX_NAME="$2"; shift 2 ;;
    --name=*)
      SANDBOX_NAME="${1#--name=}"
      if [ -z "$SANDBOX_NAME" ]; then echo "error: --name= needs a value" >&2; exit 1; fi
      shift ;;
    *) break ;;
  esac
done

CONTAINER_NAME="dev-sandbox"
[ -n "$SANDBOX_NAME" ] && CONTAINER_NAME="dev-sandbox-$SANDBOX_NAME"
IMAGE="dev-sandbox:0.1.0"

WORKSPACE="$(pwd)"
# Resolve the script's own directory so `build` finds its Dockerfile and the
# first-start provisioner regardless of where the script is invoked from.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

cmd_start() {
  if docker ps -q -f name="^${CONTAINER_NAME}$" | grep -q .; then
    echo "Already running. Use 'bash'|'zellij'|'run ...' to connect."
    exit 0
  fi

  # Container exists but stopped -- remove it so docker run starts fresh
  if docker ps -aq -f name="^${CONTAINER_NAME}$" | grep -q .; then
    echo "Removing stopped container..."
    docker rm "$CONTAINER_NAME"
  fi

  # Forward the host's git identity so commits made inside default to the
  # user, not a generic identity. Falls back to the agent name if unset.
  HOST_GIT_NAME="$(git -C "$WORKSPACE" config user.name 2>/dev/null || true)"
  HOST_GIT_EMAIL="$(git -C "$WORKSPACE" config user.email 2>/dev/null || true)"
  GIT_NAME="${HOST_GIT_NAME:-Claude}"
  GIT_EMAIL="${HOST_GIT_EMAIL:-noreply@anthropic.com}"

  # All persisted sandbox state lives under one host-owned base dir. A bind
  # dir (not a named volume) so ownership follows the host uid we run as.
  # Shared across every dev-sandbox container regardless of -n.
  STATE_BASE="$HOME/.dev-sandbox"

  # Agent + tool installs and auth, persisted so one install/login survives
  # container recreation:
  #   local           Claude Code + uv-managed Python + uv tools (ruff);
  #                   also opencode auth + data (~/.local/share/opencode)
  #   opencode        opencode agent binary (~/.opencode/bin)
  #   opencode-config opencode global config (~/.config/opencode)
  #   gh / git        gh CLI auth token + container-side global git config
  #   elan            elan + Lean toolchains
  #   mathlib-cache   Mathlib prebuilt olean download cache
  #   rustup          rustup + Rust toolchains
  #   cargo           cargo bin (rust tools) + registry cache + mold config
  SANDBOX_STATE="$STATE_BASE/local";             mkdir -p "$SANDBOX_STATE"
  OPENCODE_STATE="$STATE_BASE/opencode";         mkdir -p "$OPENCODE_STATE"
  OPENCODE_CONFIG="$STATE_BASE/opencode-config"; mkdir -p "$OPENCODE_CONFIG"
  GH_STATE="$STATE_BASE/gh";                     mkdir -p "$GH_STATE"
  GIT_STATE="$STATE_BASE/git";                   mkdir -p "$GIT_STATE"
  ELAN_STATE="$STATE_BASE/elan";                 mkdir -p "$ELAN_STATE"
  MATHLIB_CACHE="$STATE_BASE/mathlib-cache";     mkdir -p "$MATHLIB_CACHE"
  RUSTUP_STATE="$STATE_BASE/rustup";             mkdir -p "$RUSTUP_STATE"
  CARGO_STATE="$STATE_BASE/cargo";               mkdir -p "$CARGO_STATE"

  # Ensure the host paths bind-mounted as Claude's config exist first: docker
  # materialises a missing bind *source* as a root-owned dir, so an absent
  # ~/.claude.json would mount as a directory and Claude could not read it.
  mkdir -p "$HOME/.claude"
  touch "$HOME/.claude.json"

  # Optional read-only mount of the host's zellij config. Guarded so a host
  # without one does not abort start.
  ZELLIJ_MOUNT=()
  if [ -d "$HOME/.config/zellij" ]; then
    ZELLIJ_MOUNT=(-v "$(realpath "$HOME/.config/zellij")":/home/ubuntu/.config/zellij:ro)
  fi

  # Run as the host uid/gid (files land host-owned). Pin HOME explicitly so the
  # installers and ~/... bind mounts resolve to /home/ubuntu regardless of
  # which uid we run as.
  docker run -d -it \
    --name "$CONTAINER_NAME" \
    -u "$(id -u):$(id -g)" \
    -w /home/ubuntu/workspace \
    -e HOME=/home/ubuntu \
    -e GIT_AUTHOR_NAME="$GIT_NAME" \
    -e GIT_AUTHOR_EMAIL="$GIT_EMAIL" \
    -e GIT_COMMITTER_NAME="$GIT_NAME" \
    -e GIT_COMMITTER_EMAIL="$GIT_EMAIL" \
    "${ZELLIJ_MOUNT[@]+"${ZELLIJ_MOUNT[@]}"}" \
    -v "$GH_STATE":/home/ubuntu/.config/gh \
    -v "$GIT_STATE":/home/ubuntu/.config/git \
    -v "$OPENCODE_CONFIG":/home/ubuntu/.config/opencode \
    -v "$HOME/.claude":/home/ubuntu/.claude \
    -v "$HOME/.claude.json":/home/ubuntu/.claude.json \
    -v "$SANDBOX_STATE":/home/ubuntu/.local \
    -v "$OPENCODE_STATE":/home/ubuntu/.opencode \
    -v "$ELAN_STATE":/home/ubuntu/.elan \
    -v "$MATHLIB_CACHE":/home/ubuntu/.cache/mathlib \
    -v "$RUSTUP_STATE":/home/ubuntu/.rustup \
    -v "$CARGO_STATE":/home/ubuntu/.cargo \
    -v "$WORKSPACE":/home/ubuntu/workspace \
    "$IMAGE"

  # First start with empty state: provision agents + toolchains into the
  # mounts. Piped in from the host-side script so editing provisioning needs
  # no image rebuild. Idempotent -- a no-op on every later start. On failure
  # (e.g. a network blip mid-download) the container is left running but
  # incomplete; toolchain progress persists in the bind mounts, so a
  # `restart` cleanly re-runs the idempotent provisioner and resumes.
  if ! docker exec -i "$CONTAINER_NAME" bash -s < "$SCRIPT_DIR/provision.sh"; then
    echo "Provisioning failed; state is partial. Re-run with 'restart' (add -n NAME if you used it) to retry." >&2
    exit 1
  fi

  echo "Started as ${GIT_NAME} <${GIT_EMAIL}> [${CONTAINER_NAME}]. Use 'bash'|'zellij'|'run ...' to connect."
}

cmd_build() {
  # Build (or rebuild) the image from the Dockerfile next to this script.
  # Layer cache speeds repeat builds; pass --no-cache after `build` to force a
  # clean rebuild.
  docker build "$@" -t "$IMAGE" "$SCRIPT_DIR"
}

cmd_stop() {
  if ! docker ps -q -f name="^${CONTAINER_NAME}$" | grep -q .; then
    echo "Not running."
    # return, not exit: cmd_restart calls this and exit in a function kills
    # the whole script, so a restart of a stopped container would never reach
    # cmd_start.
    return 0
  fi
  docker stop "$CONTAINER_NAME"
  docker rm "$CONTAINER_NAME"
  echo "Stopped."
}

cmd_restart() {
  cmd_stop || true
  cmd_start
}

ensure_running() {
  docker ps -q -f name="^${CONTAINER_NAME}$" | grep -q . && return
  echo "Container not running -- starting it..."
  cmd_start
}

cmd_run() {
  ensure_running
  case "${1:-}" in
    bash)   docker exec -it "$CONTAINER_NAME" bash ;;
    zellij) docker exec -it "$CONTAINER_NAME" zellij attach --create "$CONTAINER_NAME" ;;
    *)      docker exec -it "$CONTAINER_NAME" "$@" ;;
  esac
}

case "${1:-}" in
  build)   shift; cmd_build "$@" ;;
  start)   cmd_start ;;
  stop)    cmd_stop ;;
  restart) cmd_restart ;;
  bash)    cmd_run bash ;;
  zellij)  cmd_run zellij ;;
  run)     shift; cmd_run "$@" ;;
  *)
    echo "Usage: $(basename "$0") [-n NAME] build [docker-build-flags...]|start|stop|restart|bash|zellij|run [CMD]"
    exit 1
    ;;
esac
