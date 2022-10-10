#!/bin/sh

echo
echo "Generating NeuroAnalyzer Markdown documentation.."
echo

echo "Delete old data"
rm -rf build/
rm -rf src/index.md

echo "Generate src/index.md"
cat header.md > src/index.md
./template.sh >> src/index.md

echo

julia make_md.jl

cp build/index.md ../Documentation.md
rm -rf build/