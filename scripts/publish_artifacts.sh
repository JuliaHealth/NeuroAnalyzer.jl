#! /usr/bin/env bash

# Build the NeuroAnalyzer artifact archives and upload them to neuroanalyzer.org/artifacts/.
#
# devel (default):
#   - pack the directory neuroanalyzer-<name> into neuroanalyzer-<name>-dev.tar.gz
#     (reproducible archive: unchanged contents give an identical file)
#   - upload the archive unless the server already has it
#   - update Artifacts.toml of NeuroAnalyzer.jl if the devel branch is checked out
#     (the change is not committed)
# stable (--stable, used by create_release_STABLE.sh):
#   - copy neuroanalyzer-<name>-dev.tar.gz to neuroanalyzer-<name>.tar.gz, so that the stable
#     release ships exactly the archives that devel was tested with
#   - upload the archive unless the server already has it
#
# usage: publish_artifacts.sh [--stable] [--no-upload] [--verify-only] [--force] [resources] [test-files]
#   --stable       publish the stable archives (see above)
#   --no-upload    only build (or copy) the archives locally
#   --verify-only  only check that the server has the local archives (no build, no upload)
#   --force        upload even if the server already has the archive
#   names          artifacts to publish (default: all)
#
# environment:
#   NA_ARTIFACTS_DIR                     directory with neuroanalyzer-resources/ and
#                                        neuroanalyzer-test-files/ (default:
#                                        ../NeuroAnalyzer-artifacts next to NeuroAnalyzer.jl)
#   NA_DIR                               NeuroAnalyzer.jl checkout (default: parent of the
#                                        script directory)
#   NA_FTP_HOST, NA_FTP_USER, NA_FTP_PASS  FTP login (only needed for uploading)

# fail on errors
set -euo pipefail

BASE_URL="https://neuroanalyzer.org/artifacts"
REMOTE_DIR="/artifacts"
ALL_NAMES=(resources test-files)

# preflight: required commands
require_cmd() {
  if ! command -v "$1" &>/dev/null; then
    echo "Error: required command '$1' not found in PATH." >&2
    exit 1
  fi
}

for cmd in tar gzip curl cmp stat git julia; do
  require_cmd "$cmd"
done

# preflight: arguments
CHANNEL="dev"
UPLOAD=1
VERIFY_ONLY=0
FORCE=0
NAMES=()
for arg in "$@"; do
  case "$arg" in
    --stable)      CHANNEL="stable" ;;
    --no-upload)   UPLOAD=0 ;;
    --verify-only) VERIFY_ONLY=1 ;;
    --force)       FORCE=1 ;;
    resources|test-files) NAMES+=("$arg") ;;
    *)
      echo "Error: unknown argument '$arg'." >&2
      echo "Usage: $0 [--stable] [--no-upload] [--verify-only] [--force] [resources] [test-files]" >&2
      exit 1
      ;;
  esac
done
[[ ${#NAMES[@]} -eq 0 ]] && NAMES=("${ALL_NAMES[@]}")

if [[ "$VERIFY_ONLY" -eq 0 && "$UPLOAD" -eq 1 ]]; then
  require_cmd lftp
  : "${NA_FTP_HOST:?NA_FTP_HOST not set}"
  : "${NA_FTP_USER:?NA_FTP_USER not set}"
  : "${NA_FTP_PASS:?NA_FTP_PASS not set}"
fi

# preflight: required paths (this script is in NeuroAnalyzer.jl/scripts)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export NA_DIR="${NA_DIR:-$(dirname "$SCRIPT_DIR")}"
ARTIFACTS_DIR="${NA_ARTIFACTS_DIR:-$(dirname "$NA_DIR")/NeuroAnalyzer-artifacts}"
export NA_ARTIFACTS_DIR="$ARTIFACTS_DIR"

if [[ ! -d "$ARTIFACTS_DIR" ]]; then
  echo "Error: artifacts directory not found: $ARTIFACTS_DIR (set NA_ARTIFACTS_DIR)." >&2
  exit 1
fi

SUFFIX=""
[[ "$CHANNEL" == "dev" ]] && SUFFIX="-dev"

# size of the file on the server (empty if it is not there)
remote_size() {
  curl -sfIL "$1" | tr -d '\r' | awk 'tolower($1) == "content-length:" { n = $2 } END { print n }' || true
}

# upload under a temporary name, then rename: a download never gets a partial file
# (the login is passed inside the script, not on the command line)
upload() {
  local file="$1"
  local name
  name="$(basename "$file")"
  echo "Uploading $name"
  lftp <<EOF
set cmd:fail-exit yes
set ssl:verify-certificate no
set net:max-retries 3
set net:reconnect-interval-base 5
open -u "$NA_FTP_USER","$NA_FTP_PASS" "$NA_FTP_HOST"
mkdir -p -f "$REMOTE_DIR"
put "$file" -o "$REMOTE_DIR/$name.part"
rm -f "$REMOTE_DIR/$name"
mv "$REMOTE_DIR/$name.part" "$REMOTE_DIR/$name"
bye
EOF
}

# remove an unfinished archive on errors
TMP_FILE=""
trap 'if [[ -n "$TMP_FILE" ]]; then rm -f "$TMP_FILE"; fi' EXIT

FAILED=0
for name in "${NAMES[@]}"; do
  dir="neuroanalyzer-$name"
  archive="$ARTIFACTS_DIR/$dir$SUFFIX.tar.gz"
  url="$BASE_URL/$dir$SUFFIX.tar.gz"

  if [[ "$VERIFY_ONLY" -eq 0 ]]; then
    if [[ "$CHANNEL" == "dev" ]]; then
      if [[ ! -d "$ARTIFACTS_DIR/$dir" ]]; then
        echo "Error: directory not found: $ARTIFACTS_DIR/$dir" >&2
        exit 1
      fi
      echo "Packing $dir"
      # reproducible archive: sorted entries, fixed time stamps and owners, no gzip time stamp
      TMP_FILE="$(mktemp -p "$ARTIFACTS_DIR" ".$dir.XXXXXX")"
      tar --sort=name --mtime='2000-01-01 00:00Z' --owner=0 --group=0 --numeric-owner \
          --exclude-vcs --exclude='*~' --exclude='.DS_Store' \
          -C "$ARTIFACTS_DIR" -cf - "$dir" | gzip -9n > "$TMP_FILE"
      chmod 644 "$TMP_FILE"
      source_file="$TMP_FILE"
    else
      source_file="$ARTIFACTS_DIR/$dir-dev.tar.gz"
      if [[ ! -f "$source_file" ]]; then
        echo "Error: devel archive not found: $source_file" >&2
        exit 1
      fi
    fi

    if [[ -f "$archive" ]] && cmp -s "$source_file" "$archive"; then
      echo "$(basename "$archive"): unchanged"
      [[ -n "$TMP_FILE" ]] && rm -f "$TMP_FILE"
    else
      if [[ "$CHANNEL" == "stable" && -f "$archive" ]]; then
        echo "Note: replacing $(basename "$archive"); previous stable releases referencing the old archive will fail to download it."
      fi
      if [[ -n "$TMP_FILE" ]]; then
        mv -f "$TMP_FILE" "$archive"
      else
        cp -f "$source_file" "$archive"
      fi
      echo "$(basename "$archive"): updated"
    fi
    TMP_FILE=""

    if [[ "$UPLOAD" -eq 1 ]]; then
      if [[ "$FORCE" -eq 1 || "$(remote_size "$url")" != "$(stat -c %s "$archive")" ]]; then
        upload "$archive"
      else
        echo "$(basename "$archive"): already on the server"
      fi
    fi
  fi

  # verify the server copy (by size: downloading the archives takes long)
  if [[ "$VERIFY_ONLY" -eq 1 || "$UPLOAD" -eq 1 ]]; then
    if [[ ! -f "$archive" ]]; then
      echo "Error: local archive not found: $archive" >&2
      exit 1
    fi
    if [[ "$(remote_size "$url")" == "$(stat -c %s "$archive")" ]]; then
      echo "$url: OK"
    else
      echo "Error: $url is missing or differs from $archive." >&2
      FAILED=1
    fi
  fi
done

[[ "$FAILED" -eq 1 ]] && exit 1

# devel: point Artifacts.toml to the new archives
if [[ "$CHANNEL" == "dev" && "$VERIFY_ONLY" -eq 0 ]]; then
  if [[ "$(git -C "$NA_DIR" rev-parse --abbrev-ref HEAD)" == "devel" ]]; then
    julia --startup-file=no "$SCRIPT_DIR/update_artifacts.jl" dev
    echo "Artifacts.toml updated; review and commit it: git -C \"$NA_DIR\" diff Artifacts.toml"
  else
    echo "Note: devel is not checked out in $NA_DIR; Artifacts.toml was not updated."
  fi
fi

exit 0
