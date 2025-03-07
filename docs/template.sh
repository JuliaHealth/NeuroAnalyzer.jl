echo "# NeuroAnalyzer.jl documentation"
echo ""
echo "This documentation has been generated using [Documenter.jl](https://juliadocs.github.io/Documenter.jl/stable/)."
echo ""
echo "## NeuroAnalyzer"
echo ""
echo "\`\`\`@docs"
cat ../src/na/setup.jl | grep ^function | grep -v Base. | sed s/"function "/"NeuroAnalyzer."/g | sed s/"(.*)"//g | sed s/" where {.*}"//g | sed s/"::.*$"//g | sort -u
cat ../src/na/plugins.jl | grep ^function | grep -v Base. | sed s/"function "/"NeuroAnalyzer."/g | sed s/"(.*)"//g | sed s/" where {.*}"//g | sed s/"::.*$"//g | sort -u
echo "\`\`\`"
echo ""
echo "## Utils"
echo ""
echo "\`\`\`@docs"
cat ../src/utils/*.jl | grep ^function | grep Base. | sed s/"function Base."/"NeuroAnalyzer."/g | sed s/"(.*)"//g | sed s/" where {.*}"//g | sed s/"::.*$"//g | sort -u
cat ../src/utils/*.jl | grep ^function | grep -v Base. | sed s/"function "/"NeuroAnalyzer."/g | sed s/"(.*)"//g | sed s/" where {.*}"//g | sed s/"::.*$"//g | sort -u
echo "\`\`\`"
echo ""
echo "## IO"
echo ""
echo "\`\`\`@docs"
cat ../src/io/*.jl | grep ^function | sed s/"function "/"NeuroAnalyzer."/g | sed s/"(.*)"//g | sed s/" where {.*}"//g | sed s/"::.*$"//g | sort -u
echo "\`\`\`"
echo ""
echo "## Edit"
echo ""
echo "\`\`\`@docs"
cat ../src/edit/*.jl | grep ^function | sed s/"function "/"NeuroAnalyzer."/g | sed s/"(.*)"//g | sed s/" where {.*}"//g | sed s/"::.*$"//g | sort -u
echo "\`\`\`"
echo ""
echo "## Process"
echo ""
echo "\`\`\`@docs"
cat ../src/process/*.jl | grep ^function | sed s/"function "/"NeuroAnalyzer."/g | sed s/"(.*)"//g | sed s/" where {.*}"//g | sed s/"::.*$"//g | sort -u
echo "\`\`\`"
echo ""
echo "## Locs"
echo ""
echo "\`\`\`@docs"
cat ../src/locs/*.jl | grep ^function | sed s/"function "/"NeuroAnalyzer."/g | sed s/"(.*)"//g | sed s/" where {.*}"//g | sed s/"::.*$"//g | sort -u
echo "\`\`\`"
echo ""
echo "## Analyze"
echo ""
echo "\`\`\`@docs"
cat ../src/analyze/*.jl | grep ^function | grep Statistics. | sed s/"function Statistics."/"NeuroAnalyzer."/g | sed s/"(.*)"//g | sed s/" where {.*}"//g | sed s/"::.*$"//g | sort -u
cat ../src/analyze/*.jl | grep ^function | grep -v Statistics. | sed s/"function "/"NeuroAnalyzer."/g | sed s/"(.*)"//g | sed s/" where {.*}"//g | sed s/"::.*$"//g | sort -u
echo "\`\`\`"
echo ""
echo "## Plot"
echo ""
echo "\`\`\`@docs"
cat ../src/plots/*.jl | grep ^function | sed s/"function "/"NeuroAnalyzer."/g | sed s/"(.*)"//g | sed s/" where {.*}"//g | sed s/"::.*$"//g | sort -u
echo "\`\`\`"
echo ""
echo "## GUI"
echo ""
echo "\`\`\`@docs"
cat ../src/gui/*.jl | grep ^function | sed s/"function "/"NeuroAnalyzer."/g | sed s/"(.*)"//g | sed s/" where {.*}"//g | sed s/"::.*$"//g | sort -u
echo "\`\`\`"
echo ""
echo "## Study"
echo ""
echo "\`\`\`@docs"
cat ../src/study/*.jl | grep ^function | sed s/"function "/"NeuroAnalyzer."/g | sed s/"(.*)"//g | sed s/" where {.*}"//g | sed s/"::.*$"//g | sort -u
echo "\`\`\`"
echo ""
echo "## NeuroRecorder"
echo ""
echo "\`\`\`@docs"
cat ../src/recorder/*.jl | grep ^function | sed s/"function "/"NeuroAnalyzer."/g | sed s/"(.*)"//g | sed s/" where {.*}"//g | sed s/"::.*$"//g | sort -u
echo "\`\`\`"
echo ""
echo "## NeuroStim"
echo ""
echo "\`\`\`@docs"
cat ../src/stim/*.jl | grep ^function | sed s/"function "/"NeuroAnalyzer."/g | sed s/"(.*)"//g | sed s/" where {.*}"//g | sed s/"::.*$"//g | sort -u
echo "\`\`\`"
echo ""
echo "## NeuroTester"
echo ""
echo "\`\`\`@docs"
cat ../src/tester/*.jl | grep ^function | sed s/"function "/"NeuroAnalyzer."/g | sed s/"(.*)"//g | sed s/" where {.*}"//g | sed s/"::.*$"//g | sort -u
echo "\`\`\`"
