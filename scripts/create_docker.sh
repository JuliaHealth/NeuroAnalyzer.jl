#! /usr/bin/env bash

# fail on errors
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# the Dockerfile is in the NeuroAnalyzer.jl folder (parent of this script)
cd "$(dirname "$SCRIPT_DIR")"
docker build -t neuroanalyzer .
docker login
docker tag neuroanalyzer:latest neuroanalyzer/neuroanalyzer:latest
docker push neuroanalyzer/neuroanalyzer:latest
