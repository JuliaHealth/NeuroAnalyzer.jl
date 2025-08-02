@info "Generating HTML documentation"

using Documenter
using NeuroAnalyzer

makedocs(sitename="NeuroAnalyzer.jl",
         modules=[NeuroAnalyzer],
         authors="Adam Wysokiński",
         linkcheck=false,
         remotes=nothing,
         warnonly=true,
         clean=true,
         format=Documenter.HTML(size_threshold=268435456))
