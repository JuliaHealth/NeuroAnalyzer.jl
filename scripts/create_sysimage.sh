#! /usr/bin/env bash

# fail on errors
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
rm -rf "$SCRIPT_DIR/sysimage"
mkdir -p "$SCRIPT_DIR/sysimage"
cd "$SCRIPT_DIR/sysimage"
julia "$SCRIPT_DIR/create_sysimage.jl"
chmod 644 NeuroAnalyzer.so
