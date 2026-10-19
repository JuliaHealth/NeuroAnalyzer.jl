#! /usr/bin/env bash

# fail on errors
set -euo pipefail

cd ~/Documents/Code

rm -rf ~/Documents/Code/JuliaRegistry > /dev/null
# clean up a leftover fork from a previous failed run, if any (fine if it doesn't exist)
gh repo delete AdamWysokinski/General --yes || true

gh repo fork JuliaRegistries/General --clone=false
# GitHub provisions forks asynchronously; give it a moment before cloning
sleep 5
git clone --depth 1 https://github.com/AdamWysokinski/General JuliaRegistry

cd JuliaRegistry
git checkout -b my-update
julia ~/Documents/Code/na_create_localregistry.jl "~/Documents/Code/JuliaRegistry"

cd ~/Documents/Code/NeuroAnalyzer.jl
git checkout main

# only necessary if there are new packages added to dependencies
# julia ~/Documents/Code/na_update_compat.jl
# julia ~/Documents/Code/na_compat_helper.jl

julia ~/Documents/Code/na_update_registry.jl "~/Documents/Code/JuliaRegistry"

cd ~/Documents/Code/JuliaRegistry
git push "https://$GITHUB_TOKEN@github.com/AdamWysokinski/General.git" my-update

rm -rf ~/.julia/registries/LocalRegistry
rm -rf ~/Documents/Code/JuliaRegistry
