#! /usr/bin/env bash

# Create a STABLE release of NeuroAnalyzer.jl:
# - format and analyze the code on devel
# - merge devel into main
# - drop "-dev" from the version
# - publish the "-dev" artifact archives as stable ones and recalculate the checksums in
#   Artifacts.toml (see publish_artifacts.sh)
# - tag and push to codeberg
# - optionally (--registry) update the fork of the Julia General registry
# - rebuild and upload the documentation
# - create the codeberg.org release
# - create the ZIP archive and publish it as a new version on zenodo.org
#
# usage: create_release_STABLE.sh [--registry] [--zenodo-draft] <version> (e.g. 0.25.11)
#   --registry  also update the fork of the Julia General registry (disabled by default:
#               NeuroAnalyzer is no longer published there); needs gh, NA_GITHUB_TOKEN and
#               NA_GH_USER
#   --zenodo-draft  leave the new zenodo.org version as a draft to be reviewed and published
#                   manually

# fail on errors
set -euo pipefail

# the whole script is read before it runs: it lives in the repository, and the checkouts
# below may replace this file while bash is still reading it
main() {

# preflight: required commands
require_cmd() {
  if ! command -v "$1" &>/dev/null; then
    echo "Error: required command '$1' not found in PATH." >&2
    exit 1
  fi
}

JLFMT="$HOME/.julia/bin/jlfmt"

# preflight: arguments
USAGE="Usage: $0 [--registry] [--zenodo-draft] <version> (e.g. 0.25.11)"
REGISTRY=0
ZENODO_OPTS=()
VERSION=""
for arg in "$@"; do
  case "$arg" in
    --registry) REGISTRY=1 ;;
    --zenodo-draft) ZENODO_OPTS=(--draft) ;;
    -*)
      echo "Error: unknown option '$arg'. $USAGE" >&2
      exit 1
      ;;
    *)
      if [[ -n "$VERSION" ]]; then
        echo "Error: more than one version given. $USAGE" >&2
        exit 1
      fi
      VERSION="$arg"
      ;;
  esac
done
if [[ -z "$VERSION" ]]; then
  echo "$USAGE" >&2
  exit 1
fi

for cmd in git julia find grep sed mktemp ncftpput lftp curl tar gzip cmp stat python3; do
  require_cmd "$cmd"
done
[[ "$REGISTRY" -eq 1 ]] && require_cmd gh && require_cmd xdg-open

if [[ ! -x "$JLFMT" ]]; then
  echo "Error: required command '$JLFMT' not found or not executable." >&2
  exit 1
fi

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
: "${NA_CODEBERG_TOKEN:?NA_CODEBERG_TOKEN not set}"
: "${NA_ZENODO_TOKEN:?NA_ZENODO_TOKEN not set}"
: "${NA_FTP_HOST:?NA_FTP_HOST not set}"
: "${NA_FTP_USER:?NA_FTP_USER not set}"
: "${NA_FTP_PASS:?NA_FTP_PASS not set}"
THREADS="${THREADS:-auto}"
if [[ "$REGISTRY" -eq 1 ]]; then
  : "${NA_GITHUB_TOKEN:?NA_GITHUB_TOKEN not set}"
  : "${NA_GH_USER:?NA_GH_USER not set}"
fi

# commit identity for this run only (the global git configuration is not modified)
export GIT_AUTHOR_NAME="$NA_GIT_NAME" GIT_COMMITTER_NAME="$NA_GIT_NAME"
export GIT_AUTHOR_EMAIL="$NA_GIT_EMAIL" GIT_COMMITTER_EMAIL="$NA_GIT_EMAIL"

# preflight: required paths (this script is in NeuroAnalyzer.jl/scripts)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NA_DIR="$(dirname "$SCRIPT_DIR")"
CODE_ANALYSIS_SCRIPT="$SCRIPT_DIR/code_analysis.jl"
UPDATE_REGISTRY_SCRIPT="$SCRIPT_DIR/update_registry.sh"
PUBLISH_ARTIFACTS_SCRIPT="$SCRIPT_DIR/publish_artifacts.sh"
UPDATE_ARTIFACTS_SCRIPT="$SCRIPT_DIR/update_artifacts.jl"
PUBLISH_ZENODO_SCRIPT="$SCRIPT_DIR/publish_zenodo.py"
export NA_ARTIFACTS_DIR="${NA_ARTIFACTS_DIR:-$(dirname "$NA_DIR")/NeuroAnalyzer-artifacts}"
REGISTRY_BRANCH="my-update"

REQUIRED_PATHS=("$CODE_ANALYSIS_SCRIPT")
[[ "$REGISTRY" -eq 1 ]] && REQUIRED_PATHS+=("$UPDATE_REGISTRY_SCRIPT" "$SCRIPT_DIR/update_registry.jl")
for path in "${REQUIRED_PATHS[@]}" "$PUBLISH_ARTIFACTS_SCRIPT" "$UPDATE_ARTIFACTS_SCRIPT" "$PUBLISH_ZENODO_SCRIPT" \
            "$NA_ARTIFACTS_DIR/neuroanalyzer-resources-dev.tar.gz" \
            "$NA_ARTIFACTS_DIR/neuroanalyzer-test-files-dev.tar.gz"; do
  if [[ ! -e "$path" ]]; then
    echo "Error: expected path not found: $path" >&2
    exit 1
  fi
done

# preflight: gh authentication
# gh auth status exits 0 for an invalid token (gh ≤ 2.46); make an authenticated request instead
if [[ "$REGISTRY" -eq 1 ]] && ! gh api user --silent &>/dev/null; then
  echo "Error: gh is not authenticated. Run 'gh auth login' first." >&2
  exit 1
fi

# temporary files (docs build, FTP login file) are removed on exit
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

# the artifact and zenodo scripts are copied: the checkouts below may replace or remove them
mkdir "$TMP_DIR/scripts"
cp "$PUBLISH_ARTIFACTS_SCRIPT" "$UPDATE_ARTIFACTS_SCRIPT" "$PUBLISH_ZENODO_SCRIPT" "$TMP_DIR/scripts/"
PUBLISH_ARTIFACTS_SCRIPT="$TMP_DIR/scripts/publish_artifacts.sh"
UPDATE_ARTIFACTS_SCRIPT="$TMP_DIR/scripts/update_artifacts.jl"
PUBLISH_ZENODO_SCRIPT="$TMP_DIR/scripts/publish_zenodo.py"

# preflight: the zenodo.org token must be valid
python3 "$PUBLISH_ZENODO_SCRIPT" --check

# codeberg.org API: the token is passed in a private header file, not on the command line
CODEBERG_API="https://codeberg.org/api/v1/repos/AdamWysokinski/NeuroAnalyzer.jl"
CODEBERG_AUTH="$TMP_DIR/codeberg_auth.txt"
( umask 077; printf 'Authorization: token %s\n' "$NA_CODEBERG_TOKEN" > "$CODEBERG_AUTH" )

# preflight: the token must be valid and allow creating releases (write access to the repository)
if [[ "$(curl -sf -H @"$CODEBERG_AUTH" "$CODEBERG_API" |
         python3 -c 'import json, sys; print(json.load(sys.stdin)["permissions"]["push"])' 2>/dev/null)" != "True" ]]; then
  echo "Error: NA_CODEBERG_TOKEN is invalid or has no write access to the repository on codeberg.org." >&2
  exit 1
fi

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

# the "-dev" archives must be on the server (devel was tested with them)
echo "Checking the devel artifacts"
NA_DIR="$NA_DIR" bash "$PUBLISH_ARTIFACTS_SCRIPT" --verify-only

# merge devel into main and drop "-dev" from the version strings
git checkout main
git merge --no-edit devel

PROJECT_VERSION="$(project_version)"
if [[ "$PROJECT_VERSION" != "${VERSION}-dev" ]]; then
  echo "Error: Project.toml version '$PROJECT_VERSION' after the merge does not match expected '${VERSION}-dev'." >&2
  exit 1
fi

sed -i -E "s/^(version[[:space:]]*=[[:space:]]*)\"[^\"]*\"/\1\"$VERSION\"/" Project.toml
sed -i -E "s/^(const VER = v)\"[^\"]*\"/\1\"$VERSION\"/" src/NeuroAnalyzer.jl

if [[ "$(project_version)" != "$VERSION" ]] || ! grep -q "^const VER = v\"$VERSION\"" src/NeuroAnalyzer.jl; then
  echo "Error: could not set version $VERSION in Project.toml / src/NeuroAnalyzer.jl." >&2
  exit 1
fi

# stable artifacts: copies of the "-dev" archives; Artifacts.toml gets the stable URLs and
# checksums (--keep-tree: fails if the archives differ from those referenced on devel);
# the archives are uploaded only after Artifacts.toml has been updated successfully
echo "Publishing the stable artifacts"
NA_DIR="$NA_DIR" bash "$PUBLISH_ARTIFACTS_SCRIPT" --stable --no-upload
NA_DIR="$NA_DIR" julia --startup-file=no "$UPDATE_ARTIFACTS_SCRIPT" stable --keep-tree
NA_DIR="$NA_DIR" bash "$PUBLISH_ARTIFACTS_SCRIPT" --stable

git add --all
git commit -m "update: version $VERSION"
git tag -a "$VERSION" -m "stable v$VERSION release"
git push codeberg main devel
git push codeberg "refs/tags/$VERSION"

# create the codeberg.org release for the tag (source archives are added by codeberg);
# skipped if the release already exists (e.g. when the script is re-run)
echo "Creating the codeberg.org release"
if curl -sf -o /dev/null -H @"$CODEBERG_AUTH" "$CODEBERG_API/releases/tags/$VERSION"; then
  echo "Release $VERSION already exists."
else
  RELEASE_JSON="$(python3 -c 'import json, sys; v = sys.argv[1]; print(json.dumps({
      "tag_name": v, "name": v, "target_commitish": "main", "draft": False, "prerelease": False,
      "body": f"Stable v{v} release. See [changelog](https://neuroanalyzer.org/changelog.html) for the detailed list of changes."}))' "$VERSION")"
  if ! curl -sf -o /dev/null -X POST -H @"$CODEBERG_AUTH" -H "Content-Type: application/json" \
       -d "$RELEASE_JSON" "$CODEBERG_API/releases"; then
    echo "Error: could not create the codeberg.org release $VERSION; create it manually." >&2
    exit 1
  fi
fi
echo "Release: https://codeberg.org/AdamWysokinski/NeuroAnalyzer.jl/releases/tag/$VERSION"

# update the fork of the Julia General registry and push the registry branch (main is checked out)
if [[ "$REGISTRY" -eq 1 ]]; then
  echo "Update registry"
  REGISTRY_BRANCH="$REGISTRY_BRANCH" "$UPDATE_REGISTRY_SCRIPT"
fi

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

# create ZIP archive and publish it as a new version on zenodo.org
cd "$NA_DIR"
git archive --format zip --output "../NeuroAnalyzer-$VERSION.zip" "$VERSION"

echo "Update zenodo.org"
python3 "$PUBLISH_ZENODO_SCRIPT" "${ZENODO_OPTS[@]}" "$VERSION" "../NeuroAnalyzer-$VERSION.zip"

if [[ "$REGISTRY" -eq 1 ]]; then
  echo "Opening https://github.com/$NA_GH_USER/General/tree/$REGISTRY_BRANCH"
  xdg-open "https://github.com/$NA_GH_USER/General/tree/$REGISTRY_BRANCH" &>/dev/null || true
fi
}

main "$@"
exit 0
