#! /usr/bin/env bash

# Register the current main of NeuroAnalyzer.jl in the fork of the Julia General registry
# and push the registry branch (called by create_release_STABLE.sh, can be re-run alone)
#
# usage: update_registry.sh

# fail on errors
set -euo pipefail

# the whole script is read before it runs: it lives in the repository, and the checkout
# below may replace this file while bash is still reading it
main() {

# preflight: required commands
require_cmd() {
  if ! command -v "$1" &>/dev/null; then
    echo "Error: required command '$1' not found in PATH." >&2
    exit 1
  fi
}

for cmd in git julia gh mktemp sleep base64 tr; do
  require_cmd "$cmd"
done

# preflight: required environment variables
: "${NA_GH_USER:?NA_GH_USER not set}"
: "${NA_GITHUB_TOKEN:?NA_GITHUB_TOKEN not set}"
REGISTRY_BRANCH="${REGISTRY_BRANCH:-my-update}"

# preflight: required paths (this script is in NeuroAnalyzer.jl/scripts)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NA_DIR="$(dirname "$SCRIPT_DIR")"
UPDATE_REGISTRY_SCRIPT="$SCRIPT_DIR/update_registry.jl"

if [[ ! -e "$UPDATE_REGISTRY_SCRIPT" ]]; then
  echo "Error: expected path not found: $UPDATE_REGISTRY_SCRIPT" >&2
  exit 1
fi

# preflight: gh authentication
if ! gh auth status &>/dev/null; then
  echo "Error: gh is not authenticated. Run 'gh auth login' first." >&2
  exit 1
fi

# the registry is updated from the released (main) code
cd "$NA_DIR"
if [[ -n "$(git status --porcelain)" ]]; then
  echo "Error: working tree in $NA_DIR is not clean. Commit or stash changes first." >&2
  exit 1
fi
git checkout main

# the registry clone is removed on exit
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT
REGISTRY_DIR="$TMP_DIR/JuliaRegistry"

# reuse the fork (synchronized with upstream) instead of deleting and re-forking it
if gh repo view "$NA_GH_USER/General" &>/dev/null; then
  gh repo sync "$NA_GH_USER/General" --source JuliaRegistries/General --force
else
  gh repo fork JuliaRegistries/General --clone=false
  # GitHub provisions forks asynchronously
  for _ in {1..30}; do
    gh repo view "$NA_GH_USER/General" &>/dev/null && break
    sleep 2
  done
fi
git clone --depth 1 "https://github.com/$NA_GH_USER/General" "$REGISTRY_DIR"
git -C "$REGISTRY_DIR" checkout -b "$REGISTRY_BRANCH"

# only necessary if there are new packages added to dependencies
# julia scripts/update_compat.jl

julia "$UPDATE_REGISTRY_SCRIPT" "$REGISTRY_DIR"

# the token is passed in the environment, not on the command line (visible in the process list)
GITHUB_AUTH_HEADER="Authorization: Basic $(printf 'x-access-token:%s' "$NA_GITHUB_TOKEN" | base64 | tr -d '\n')"
GIT_CONFIG_COUNT=1 GIT_CONFIG_KEY_0=http.extraHeader GIT_CONFIG_VALUE_0="$GITHUB_AUTH_HEADER" \
  git -C "$REGISTRY_DIR" push --force "https://github.com/$NA_GH_USER/General.git" "$REGISTRY_BRANCH"
unset GITHUB_AUTH_HEADER
}

main "$@"
exit 0
