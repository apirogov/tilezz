#!/usr/bin/env bash
set -euo pipefail

CONTAINER_NAME="rust-sandbox"
IMAGE="rust-sandbox:0.0.1"

WORKSPACE="$(pwd)"
# Resolve the script's own directory so `build` finds its Dockerfile
# regardless of where the script is invoked from.
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

  # Forward the host's git identity into the container so commits made
  # inside default to the user, not to a generic "Claude" identity.
  # Falls back to the agent name if the host has nothing configured.
  HOST_GIT_NAME="$(git -C "$WORKSPACE" config user.name 2>/dev/null || true)"
  HOST_GIT_EMAIL="$(git -C "$WORKSPACE" config user.email 2>/dev/null || true)"
  GIT_NAME="${HOST_GIT_NAME:-Claude}"
  GIT_EMAIL="${HOST_GIT_EMAIL:-noreply@anthropic.com}"

  # Persist the container's ~/.local on the host. Claude Code installs
  # itself under ~/.local/share/claude and symlinks ~/.local/bin/claude;
  # a durable mount means one install + auto-updates that survive
  # container recreation. Kept separate from the host's own ~/.local so
  # sandbox state never mixes with host tooling. A bind dir (not a named
  # volume) so ownership follows the host uid we run the container as.
  SANDBOX_STATE="$HOME/.rust-sandbox/local"
  mkdir -p "$SANDBOX_STATE"

  docker run -d -it \
    --name "$CONTAINER_NAME" \
    -u $(id -u):$(id -g) \
    -w /home/ubuntu/workspace \
    -e GIT_AUTHOR_NAME="$GIT_NAME" \
    -e GIT_AUTHOR_EMAIL="$GIT_EMAIL" \
    -e GIT_COMMITTER_NAME="$GIT_NAME" \
    -e GIT_COMMITTER_EMAIL="$GIT_EMAIL" \
    -v $(realpath $HOME/.config/zellij):/home/ubuntu/.config/zellij:ro \
    -v $HOME/.claude:/home/ubuntu/.claude \
    -v $HOME/.claude.json:/home/ubuntu/.claude.json \
    -v "$SANDBOX_STATE":/home/ubuntu/.local \
    -v "$WORKSPACE":/home/ubuntu/workspace \
    "$IMAGE"

  # First start with an empty state dir: install Claude Code into the
  # mounted ~/.local. No-op on every later start.
  docker exec "$CONTAINER_NAME" bash -c \
    'command -v claude >/dev/null || curl -fsSL https://claude.ai/install.sh | bash'

  echo "Started as ${GIT_NAME} <${GIT_EMAIL}>. Use 'bash'|'zellij'|'run ...' to connect."
}

cmd_build() {
  # Build (or rebuild) the image from the Dockerfile next to this
  # script. Layer cache speeds repeat builds; pass --no-cache after
  # `build` to force a clean rebuild.
  docker build "$@" -t "$IMAGE" "$SCRIPT_DIR"
}

cmd_stop() {
  if ! docker ps -q -f name="^${CONTAINER_NAME}$" | grep -q .; then
    echo "Not running."
    exit 0
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
    echo "Usage: $(basename $0) build [docker-build-flags...]|start|stop|restart|bash|zellij|run [CMD]"
    exit 1
    ;;
esac
