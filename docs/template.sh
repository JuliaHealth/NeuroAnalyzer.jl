echo "# NeuroAnalyzer.jl documentation"
echo ""
echo "This documentation has been generated using [Documenter.jl](https://juliadocs.github.io/Documenter.jl/stable/)."
echo ""
echo "## NeuroAnalyzer"
echo ""
echo "\`\`\`@docs"
cat ../src/na.jl | grep ^function | sed s/"function "/"NeuroAnalyzer."/g
echo "\`\`\`"
echo ""
echo "## Utils"
echo ""
echo "\`\`\`@autodocs"
echo "Modules = [NeuroAnalyzer]"
echo "Order   = [:function, :type]"
echo "Pages   = readdir(\"../src/utils\")"
echo "\`\`\`"
echo ""
echo "## IO"
echo ""
echo "\`\`\`@autodocs"
echo "Modules = [NeuroAnalyzer]"
echo "Order   = [:function, :type]"
echo "Pages   = readdir(\"../src/io\")"
echo "\`\`\`"
echo ""
echo "## Edit"
echo ""
echo "\`\`\`@autodocs"
echo "Modules = [NeuroAnalyzer]"
echo "Order   = [:function, :type]"
echo "Pages   = readdir(\"../src/edit\")"
echo "\`\`\`"
echo ""
echo "## Process"
echo ""
echo "\`\`\`@autodocs"
echo "Modules = [NeuroAnalyzer]"
echo "Order   = [:function, :type]"
echo "Pages   = readdir(\"../src/process\")"
echo "\`\`\`"
echo ""
echo "## Locs"
echo ""
echo "\`\`\`@autodocs"
echo "Modules = [NeuroAnalyzer]"
echo "Order   = [:function, :type]"
echo "Pages   = readdir(\"../src/locs\")"
echo "\`\`\`"
echo ""
echo "## Analyze"
echo ""
echo "\`\`\`@autodocs"
echo "Modules = [NeuroAnalyzer]"
echo "Order   = [:function, :type]"
echo "Pages   = readdir(\"../src/analyze\")"
echo "\`\`\`"
echo ""
echo "## Plot"
echo ""
echo "\`\`\`@autodocs"
echo "Modules = [NeuroAnalyzer]"
echo "Order   = [:function, :type]"
echo "Pages   = readdir(\"../src/plot\")"
echo "\`\`\`"
echo ""
echo "## GUI"
echo ""
echo "\`\`\`@autodocs"
echo "Modules = [NeuroAnalyzer]"
echo "Order   = [:function, :type]"
echo "Pages   = readdir(\"../src/gui\")"
echo "\`\`\`"
echo ""
echo "## Statistics"
echo ""
echo "\`\`\`@autodocs"
echo "Modules = [NeuroAnalyzer]"
echo "Order   = [:function, :type]"
echo "Pages   = readdir(\"../src/statistics\")"
echo "\`\`\`"
echo ""
echo "## Study"
echo ""
echo "\`\`\`@autodocs"
echo "Modules = [NeuroAnalyzer]"
echo "Order   = [:function, :type]"
echo "Pages   = readdir(\"../src/study\")"
echo "\`\`\`"
echo ""
echo "## NeuroRecorder"
echo ""
echo "\`\`\`@autodocs"
echo "Modules = [NeuroAnalyzer]"
echo "Order   = [:function, :type]"
echo "Pages   = readdir(\"../src/recorder\")"
echo "\`\`\`"
echo ""
echo "## NeuroStim"
echo ""
echo "\`\`\`@autodocs"
echo "Modules = [NeuroAnalyzer]"
echo "Order   = [:function, :type]"
echo "Pages   = readdir(\"../src/stim\")"
echo "\`\`\`"
echo ""
echo "## NeuroTester"
echo ""
echo "\`\`\`@autodocs"
echo "Modules = [NeuroAnalyzer]"
echo "Order   = [:function, :type]"
echo "Pages   = readdir(\"../src/tester\")"
echo "\`\`\`"
