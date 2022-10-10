#!/bin/sh

echo
echo "Generating NeuroAnalyzer HTML documentation.."
echo

echo "Delete old data"
rm -rf build/

echo "Generate src/index.md"
rm -f src/index.md
cat header.md > src/index.md
./template.sh >> src/index.md

echo

julia make_html.jl

mv build docs
sed -i 's/Edit on GitHub/Edit on Codeberg/g' docs/index.html
sed -i 's///g' docs/index.html
