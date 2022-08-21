#!/bin/sh

rm -rf build/
rm -rf src/index.md
./template.sh > src/index.md
julia make.jl
cp build/index.md ../Documentation.md
rm -rf build/