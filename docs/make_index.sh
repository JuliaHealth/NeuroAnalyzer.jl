#! /usr/bin/env bash

# fail on errors
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
cd "$SCRIPT_DIR"
 
rm -rf build
cat header.md > src/index.md
chmod +x template.sh
./template.sh >> src/index.md
 