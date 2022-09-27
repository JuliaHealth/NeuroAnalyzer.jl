#!/bin/sh

echo
echo "Generating NeuroAnalyzer HTML documentation.."
echo

echo "Delete old data"
rm -rf build/

echo "Generate src/index.md"
cat header.md > src/index.md
./template.sh >> src/index.md

echo "Generate HTML documentation"
julia make_html.jl

mv build docs
mv docs/index.html docs/docs.html