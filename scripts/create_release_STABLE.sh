#! /usr/bin/env bash

# Create a STABLE release of NeuroAnalyzer.jl:
# - format and analyze the code on devel
# - merge devel into main
# - drop "-dev" from the version
# - tag and push to codeberg
# - update the fork of the Julia General registry
# - rebuild and upload the documentation
# - create the ZIP archive
#
# usage: na_create_release_STABLE.sh <version> (e.g. 0.25.11)

# fail on errors
set -euo pipefail

# preflight: required commands
require_cmd() {
  if ! command -v "$1" &>/dev/null; then
    echo "Error: required command '$1' not found in PATH." >&2
    exit 1
  fi
}

JLFMT="$HOME/.julia/bin/jlfmt"

for cmd in git julia gh find grep sed mktemp ncftpput xdg-open base64 tr; do
  require_cmd "$cmd"
done

if [[ ! -x "$JLFMT" ]]; then
  echo "Error: required command '$JLFMT' not found or not executable." >&2
  exit 1
fi

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

echo "Creating STABLE release: $VERSION"

# preflight: required environment variables — fail fast rather than partway through
: "${NA_GIT_NAME:?NA_GIT_NAME not set}"
: "${NA_GIT_EMAIL:?NA_GIT_EMAIL not set}"
: "${NA_GITHUB_TOKEN:?NA_GITHUB_TOKEN not set}"
: "${NA_FTP_USER:?NA_FTP_USER not set}"
: "${NA_FTP_PASS:?NA_FTP_PASS not set}"
THREADS="${THREADS:-auto}"

# commit identity for this run only (the global git configuration is not modified)
export GIT_AUTHOR_NAME="$NA_GIT_NAME" GIT_COMMITTER_NAME="$NA_GIT_NAME"
export GIT_AUTHOR_EMAIL="$NA_GIT_EMAIL" GIT_COMMITTER_EMAIL="$NA_GIT_EMAIL"

# preflight: required paths
CODE_ANALYSIS_SCRIPT="$NA_DIR/scripts/code_analysis.jl"
UPDATE_REGISTRY_SCRIPT="$NA_DIR/scripts/update_registry.jl"
REGISTRY_BRANCH="my-update"

for path in "$NA_DIR" "$CODE_ANALYSIS_SCRIPT" "$UPDATE_REGISTRY_SCRIPT"; do
  if [[ ! -e "$path" ]]; then
    echo "Error: expected path not found: $path" >&2
    exit 1
  fi
done

# preflight: gh authentication
if ! gh auth status &>/dev/null; then
  echo "Error: gh is not authenticated. Run 'gh auth login' first." >&2
  exit 1
fi

# temporary files (registry clone, docs build, FTP login file) are removed on exit
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

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

# working tree should be clean before we start rewriting files and merging
if [[ -n "$(git status --porcelain)" ]]; then
  echo "Error: working tree in $NA_DIR is not clean. Commit or stash changes first." >&2
  exit 1
fi

# local branches must contain everything that is on codeberg (pushes are never forced)
git fetch codeberg
for branch in devel main; do
  if ! git merge-base --is-ancestor "codeberg/$branch" "$branch"; then
    echo "Error: local '$branch' does not contain codeberg/$branch. Pull first." >&2
    exit 1
  fi
done

# the release tag must not exist yet
if git rev-parse -q --verify "refs/tags/$VERSION" &>/dev/null ||
   [[ -n "$(git ls-remote --tags codeberg "refs/tags/$VERSION")" ]]; then
  echo "Error: tag '$VERSION' already exists." >&2
  exit 1
fi

# format and analyze the code on devel; formatting changes are committed there
git checkout devel

PROJECT_VERSION="$(project_version)"
if [[ "$PROJECT_VERSION" != "${VERSION}-dev" ]]; then
  echo "Error: devel Project.toml version '$PROJECT_VERSION' does not match expected '${VERSION}-dev'." >&2
  exit 1
fi

echo "Formatting the code"
find . -name "*.jl" -not -path "./.git/*" -exec "$JLFMT" --prioritize-config-file --inplace --verbose {} +
if [[ -n "$(git status --porcelain)" ]]; then
  git add --all
  git commit -m "update: code formatting"
fi

echo "Running code analysis"
julia -q --color=yes -O2 -g0 --threads "$THREADS" --cpu-target=native "$CODE_ANALYSIS_SCRIPT"

# merge devel into main and drop "-dev" from the version strings and artifact file names
git checkout main
git merge --no-edit devel

PROJECT_VERSION="$(project_version)"
if [[ "$PROJECT_VERSION" != "${VERSION}-dev" ]]; then
  echo "Error: Project.toml version '$PROJECT_VERSION' after the merge does not match expected '${VERSION}-dev'." >&2
  exit 1
fi

sed -i -E "s/^(version[[:space:]]*=[[:space:]]*)\"[^\"]*\"/\1\"$VERSION\"/" Project.toml
sed -i -E "s/^(const VER = v)\"[^\"]*\"/\1\"$VERSION\"/" src/NeuroAnalyzer.jl
sed -i -E 's/-dev\.tar\.gz/.tar.gz/g' Artifacts.toml

if [[ "$(project_version)" != "$VERSION" ]] || ! grep -q "^const VER = v\"$VERSION\"" src/NeuroAnalyzer.jl; then
  echo "Error: could not set version $VERSION in Project.toml / src/NeuroAnalyzer.jl." >&2
  exit 1
fi

git add --all
git commit -m "update: version $VERSION"
git tag -a "$VERSION" -m "stable v$VERSION release"
git push codeberg main devel
git push codeberg "refs/tags/$VERSION"

# update the fork of the Julia General registry and push the registry branch
echo "Update registry"

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
# julia scripts/na_update_compat.jl

julia "$UPDATE_REGISTRY_SCRIPT" "$REGISTRY_DIR"

# the token is passed in the environment, not on the command line (visible in the process list)
GITHUB_AUTH_HEADER="Authorization: Basic $(printf 'x-access-token:%s' "$NA_GITHUB_TOKEN" | base64 | tr -d '\n')"
GIT_CONFIG_COUNT=1 GIT_CONFIG_KEY_0=http.extraHeader GIT_CONFIG_VALUE_0="$GITHUB_AUTH_HEADER" \
  git -C "$REGISTRY_DIR" push --force "https://github.com/$NA_GH_USER/General.git" "$REGISTRY_BRANCH"
unset GITHUB_AUTH_HEADER

# re-generate the documentation from the released code and upload it
echo "Update docs"

DOCS_DIR="$TMP_DIR/update_docs"
git clone --depth 1 --branch "$VERSION" https://codeberg.org/AdamWysokinski/NeuroAnalyzer.jl "$DOCS_DIR/NeuroAnalyzer.jl"
cd "$DOCS_DIR/NeuroAnalyzer.jl/docs"

for f in header.md make_index.sh template.jl make_html.jl; do
  if [[ ! -e "$f" ]]; then
    echo "Error: expected file not found in docs checkout: $f" >&2
    exit 1
  fi
done
if [[ ! -x make_index.sh ]]; then
  echo "Error: make_index.sh is present but not executable." >&2
  exit 1
fi

# make_index.sh removes build/, instantiates the package environment of the checkout and
# writes src/index.md (header.md + @docs blocks generated by template.jl from the docstrings
# of the loaded package)
./make_index.sh
julia --project make_html.jl
mv build docs
sed -i 's/main\/docs\/src\/index.md/main\/docs\//g' docs/index.html

# FTP credentials are read from a private login file, not passed on the command line
FTP_LOGIN="$TMP_DIR/ftp_login.cfg"
( umask 077; printf 'host %s\nuser %s\npass %s\n' "$NA_FTP_HOST" "$NA_FTP_USER" "$NA_FTP_PASS" > "$FTP_LOGIN" )
ncftpput -R -f "$FTP_LOGIN" / docs

# create ZIP archive and upload to codeberg.org and zenodo.org (these have to be done manually)
cd "$NA_DIR"
git archive --format zip --output "../NeuroAnalyzer-$VERSION.zip" "$VERSION"

echo "Opening https://codeberg.org/AdamWysokinski/NeuroAnalyzer.jl/releases"
echo "Stable v$VERSION release. See [changelog](https://neuroanalyzer.org/changelog.html) for the detailed list of changes."
xdg-open https://codeberg.org/AdamWysokinski/NeuroAnalyzer.jl/releases &>/dev/null

echo "Opening https://zenodo.org/records/16964565"
xdg-open https://zenodo.org/records/16964565 &>/dev/null

echo "Opening https://github.com/$NA_GH_USER/General/tree/$REGISTRY_BRANCH"
xdg-open "https://github.com/$NA_GH_USER/General/tree/$REGISTRY_BRANCH" &>/dev/null

exit 0
