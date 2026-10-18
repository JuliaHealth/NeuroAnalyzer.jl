#! /usr/bin/env bash

# fail on errors
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
cd "$SCRIPT_DIR"
 
rm -rf build
cat header.md > src/index.md
# the index is generated from the docstrings of the loaded package (see template.jl);
# the package environment must be instantiated (e.g. in a fresh clone)
julia --project=.. -e 'using Pkg; Pkg.instantiate()'
julia --project=.. template.jl src/index.md
 