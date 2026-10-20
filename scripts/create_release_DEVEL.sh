#! /usr/bin/env bash

# Initiate a DEVEL release of NeuroAnalyzer.jl:
# - merge main into devel
# - set the version to "<version>-dev"
# - switch the artifacts to the "-dev" archives
# - push devel to codeberg
#
# usage: create_release_DEVEL.sh <version> (e.g. 0.25.11; "-dev" is added)

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

for cmd in git grep sed sort; do
  require_cmd "$cmd"
done

# preflight: arguments
VERSION="${1:?Usage: $0 <version> (e.g. 0.25.11)}"

if [[ "$VERSION" == *"-dev"* ]]; then
  echo "Error: VERSION must not contain -dev, just the version (e.g. 0.25.11)." >&2
  exit 1
fi

# MAJOR.MINOR.PATCH, no zero padding (0.25.11, not 0.25.011)
if [[ ! "$VERSION" =~ ^(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)$ ]]; then
  echo "Error: VERSION '$VERSION' is not a semantic version without zero padding (e.g. 0.25.11)." >&2
  exit 1
fi

DEV_VERSION="${VERSION}-dev"

echo "Initiating DEVEL release: $DEV_VERSION"

# preflight: required environment variables
: "${NA_GIT_NAME:?NA_GIT_NAME not set}"
: "${NA_GIT_EMAIL:?NA_GIT_EMAIL not set}"

# commit identity for this run only (the global git configuration is not modified)
export GIT_AUTHOR_NAME="$NA_GIT_NAME" GIT_COMMITTER_NAME="$NA_GIT_NAME"
export GIT_AUTHOR_EMAIL="$NA_GIT_EMAIL" GIT_COMMITTER_EMAIL="$NA_GIT_EMAIL"

# preflight: required paths (this script is in NeuroAnalyzer.jl/scripts)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NA_DIR="$(dirname "$SCRIPT_DIR")"

# version from a "version = ..." line of Project.toml (empty if there is none)
project_version() {
  grep -m1 -E '^version[[:space:]]*=' Project.toml | sed -E 's/^version[[:space:]]*=[[:space:]]*"([^"]*)"/\1/' || true
}

cd "$NA_DIR"

# codeberg remote must exist before we try to push to it
if ! git remote get-url codeberg &>/dev/null; then
  echo "Error: git remote 'codeberg' is not configured in $NA_DIR." >&2
  exit 1
fi

# working tree should be clean before we start merging and rewriting files
if [[ -n "$(git status --porcelain)" ]]; then
  echo "Error: working tree in $NA_DIR is not clean. Commit or stash changes first." >&2
  exit 1
fi

for f in Artifacts.toml Project.toml src/NeuroAnalyzer.jl; do
  if [[ ! -e "$f" ]]; then
    echo "Error: expected file not found: $NA_DIR/$f" >&2
    exit 1
  fi
done

# local branches must contain everything that is on codeberg (pushes are never forced)
git fetch codeberg
for branch in devel main; do
  if ! git merge-base --is-ancestor "codeberg/$branch" "$branch"; then
    echo "Error: local '$branch' does not contain codeberg/$branch. Pull first." >&2
    exit 1
  fi
done

git checkout devel
git merge --no-edit main

# sanity check: after merging main, Project.toml should declare a stable (non "-dev") version
# since main only ever holds stable releases
PROJECT_VERSION="$(project_version)"

if [[ -z "$PROJECT_VERSION" ]]; then
  echo "Error: could not find a 'version = \"...\"' line in Project.toml." >&2
  exit 1
fi

if [[ "$PROJECT_VERSION" == *"-dev" ]]; then
  echo "Error: Project.toml version '$PROJECT_VERSION' is a -dev version, expected a stable version after merging main. Refusing to proceed." >&2
  exit 1
fi

# the new development version must be newer than the last stable release
if [[ "$VERSION" == "$PROJECT_VERSION" ]] ||
   [[ "$(printf '%s\n%s\n' "$PROJECT_VERSION" "$VERSION" | sort -V | tail -n 1)" != "$VERSION" ]]; then
  echo "Error: VERSION '$VERSION' must be newer than the current stable version '$PROJECT_VERSION'." >&2
  exit 1
fi

# "-dev" version strings and artifact archives (idempotent: an existing "-dev" is not doubled)
sed -i -E "s/^(version[[:space:]]*=[[:space:]]*)\"[^\"]*\"/\1\"$DEV_VERSION\"/" Project.toml
sed -i -E "s/^(const VER = v)\"[^\"]*\"/\1\"$DEV_VERSION\"/" src/NeuroAnalyzer.jl
sed -i -E 's/(-dev)?\.tar\.gz/-dev.tar.gz/g' Artifacts.toml

if [[ "$(project_version)" != "$DEV_VERSION" ]] || ! grep -q "^const VER = v\"$DEV_VERSION\"" src/NeuroAnalyzer.jl; then
  echo "Error: could not set version $DEV_VERSION in Project.toml / src/NeuroAnalyzer.jl." >&2
  exit 1
fi

git add --all
git commit -m "$DEV_VERSION initiated"
git push codeberg devel
}

main "$@"
exit 0
