#! /usr/bin/env bash

# fail on errors
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

cd "$SCRIPT_DIR/NeuroAnalyzer.jl"
docker build -t neuroanalyzer .
docker login
docker tag neuroanalyzer:latest neuroanalyzer/neuroanalyzer:latest
docker push neuroanalyzer/neuroanalyzer:latest
cd "$SCRIPT_DIR"
