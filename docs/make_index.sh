#! /bin/sh

rm -rf build
rm -f src/index.md
cat header.md > src/index.md
./template.sh >> src/index.md
