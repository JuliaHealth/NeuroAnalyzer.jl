#! /usr/bin/env bash

# fail on errors
set -euo pipefail

# grep exits 1 when a category legitimately has zero matches (e.g. no Base. extensions in a given file)
# that is not a real error, so this wrapper tolerates a clean "no match" result while still letting genuine
# grep errors (bad file, bad option, etc.) propagate under set -e
sgrep() { grep "$@" || [ "$?" = 1 ]; }

echo "# NeuroAnalyzer.jl documentation"
echo ""
echo "This documentation has been generated using [Documenter.jl](https://juliadocs.github.io/Documenter.jl/stable/)."
echo ""
echo "## NeuroAnalyzer"
echo ""
echo "\`\`\`@docs"
cat ../src/na/setup.jl | sgrep ^function | sgrep -v ^"function _" | sgrep -v Base. | sed s/"function "/"NeuroAnalyzer."/g | sed s/"(.*)"//g | sed s/" where {.*}"//g | sed s/"::.*$"//g | sed s/"(;"//g | sed s/"("//g | sort -u
cat ../src/na/plugins.jl | sgrep ^function | sgrep -v ^"function _" | sgrep -v Base. | sed s/"function "/"NeuroAnalyzer."/g | sed s/"(.*)"//g | sed s/" where {.*}"//g | sed s/"::.*$"//g | sed s/"(;"//g | sed s/"("//g | sort -u
echo "\`\`\`"
echo ""
echo "---"
echo ""
echo "## Utils"
echo ""
echo "\`\`\`@docs"
cat *.jl | sgrep ^function | sgrep -v ^"function _" | sgrep Base. | sed s/"function Base."/"NeuroAnalyzer."/g | sed s/"(.*)"//g | sed s/" where {.*}"//g | sed s/"::.*$"//g | sed s/"(;"//g | sed s/"("//g | sort -u
cat ../src/utils/*.jl | sgrep ^function | sgrep -v ^"function _" | sgrep Base. | sed s/"function Base."/"NeuroAnalyzer."/g | sed s/"(.*)"//g | sed s/" where {.*}"//g | sed s/"::.*$"//g | sed s/"(;"//g | sed s/"("//g | sort -u
cat ../src/utils/*.jl | sgrep ^function | sgrep -v ^"function _" | sgrep -v Base. | sed s/"function "/"NeuroAnalyzer."/g | sed s/"(.*)"//g | sed s/" where {.*}"//g | sed s/"::.*$"//g | sed s/"(;"//g | sed s/"("//g | sort -u
echo "\`\`\`"
echo ""
echo "---"
echo ""
echo "## Stats"
echo ""
echo "\`\`\`@docs"
cat ../src/stats/*.jl | sgrep ^function | sgrep -v ^"function _" | sed s/"function "/"NeuroAnalyzer."/g | sed s/"(.*)"//g | sed s/" where {.*}"//g | sed s/"::.*$"//g | sed s/"(;"//g | sed s/"("//g | sort -u
echo "\`\`\`"
echo ""
echo "---"
echo ""
echo "## IO"
echo ""
echo "\`\`\`@docs"
cat ../src/io/*.jl | sgrep ^function | sgrep -v ^"function _" | sed s/"function "/"NeuroAnalyzer."/g | sed s/"(.*)"//g | sed s/" where {.*}"//g | sed s/"::.*$"//g | sed s/"(;"//g | sed s/"("//g | sort -u
echo "\`\`\`"
echo ""
echo "---"
echo ""
echo "## Edit"
echo ""
echo "\`\`\`@docs"
cat ../src/edit/*.jl | sgrep ^function | sgrep -v ^"function _" | sed s/"function "/"NeuroAnalyzer."/g | sed s/"(.*)"//g | sed s/" where {.*}"//g | sed s/"::.*$"//g | sed s/"(;"//g | sed s/"("//g | sort -u
echo "\`\`\`"
echo ""
echo "---"
echo ""
echo "## Process"
echo ""
echo "\`\`\`@docs"
cat ../src/process/*.jl | sgrep ^function | sgrep -v ^"function _" | sed s/"function "/"NeuroAnalyzer."/g | sed s/"(.*)"//g | sed s/" where {.*}"//g | sed s/"::.*$"//g | sed s/"(;"//g | sed s/"("//g | sort -u
echo "\`\`\`"
echo ""
echo "---"
echo ""
echo "## Locs"
echo ""
echo "\`\`\`@docs"
cat ../src/locs/*.jl | sgrep ^function | sed s/"function "/"NeuroAnalyzer."/g | sed s/"(.*)"//g | sed s/" where {.*}"//g | sed s/"::.*$"//g | sed s/"(;"//g | sed s/"("//g | sort -u
echo "\`\`\`"
echo ""
echo "---"
echo ""
echo "## Analyze"
echo ""
echo "\`\`\`@docs"
cat ../src/analyze/*.jl | sgrep ^function | sgrep -v ^"function _" | sgrep Statistics. | sed s/"function Statistics."/"NeuroAnalyzer."/g | sed s/"(.*)"//g | sed s/" where {.*}"//g | sed s/"::.*$"//g | sed s/"(;"//g | sed s/"("//g | sort -u
cat ../src/analyze/*.jl | sgrep ^function | sgrep -v ^"function _" | sgrep -v Statistics. | sed s/"function "/"NeuroAnalyzer."/g | sed s/"(.*)"//g | sed s/" where {.*}"//g | sed s/"::.*$"//g | sed s/"(;"//g | sed s/"("//g | sort -u
echo "\`\`\`"
echo ""
echo "---"
echo ""
echo "## Model"
echo ""
echo "\`\`\`@docs"
cat ../src/model/*.jl | sgrep ^function | sed s/"function "/"NeuroAnalyzer."/g | sed s/"(.*)"//g | sed s/" where {.*}"//g | sed s/"::.*$"//g | sed s/"(;"//g | sed s/"("//g | sort -u
echo "\`\`\`"
echo ""
echo "---"
echo ""
echo "## Plot"
echo ""
echo "\`\`\`@docs"
cat ../src/plots/*.jl | sgrep ^function | sgrep -v ^"function _" | sed s/"function "/"NeuroAnalyzer."/g | sed s/"(.*)"//g | sed s/" where {.*}"//g | sed s/"::.*$"//g | sed s/"(;"//g | sed s/"("//g | sort -u
echo "\`\`\`"
echo ""
echo "---"
echo ""
# echo "## GUI"
# echo ""
# echo "\`\`\`@docs"
# cat ../src/gui/*.jl | sgrep ^function | sgrep -v ^"function _" | sed s/"function "/"NeuroAnalyzer."/g | sed s/"(.*)"//g | sed s/" where {.*}"//g | sed s/"::.*$"//g | sed s/"(;"//g | sed s/"("//g | sort -u
# echo "\`\`\`"
# echo ""
# echo "---"
# echo ""
echo "## NeuroRecorder"
echo ""
echo "\`\`\`@docs"
cat ../src/recorder/*.jl | sgrep ^function | sgrep -v ^"function _" | sed s/"function "/"NeuroAnalyzer."/g | sed s/"(.*)"//g | sed s/" where {.*}"//g | sed s/"::.*$"//g | sed s/"(;"//g | sed s/"("//g | sort -u
echo "\`\`\`"
echo ""
echo "---"
echo ""
echo "## NeuroStim"
echo ""
echo "\`\`\`@docs"
cat ../src/stim/*.jl | sgrep ^function | sgrep -v ^"function _" | sed s/"function "/"NeuroAnalyzer."/g | sed s/"(.*)"//g | sed s/" where {.*}"//g | sed s/"::.*$"//g | sed s/"(;"//g | sed s/"("//g | sort -u
echo "\`\`\`"
echo ""
echo "---"
echo ""
echo "## NeuroTester"
echo ""
echo "\`\`\`@docs"
cat ../src/tester/*.jl | sgrep ^function | sgrep -v ^"function _" | sed s/"function "/"NeuroAnalyzer."/g | sed s/"(.*)"//g | sed s/" where {.*}"//g | sed s/"::.*$"//g | sed s/"(;"//g | sed s/"("//g | sort -u
echo "\`\`\`"
