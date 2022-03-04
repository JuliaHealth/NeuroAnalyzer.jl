#!/bin/sh

rm -rf build/
julia make.jl
cp build/index.md ../Documentation.md