@info "Generating HTML documentation"

using Documenter

using Pkg
Pkg.activate("..")
using NeuroAnalyzer

makedocs(sitename="NeuroAnalyzer.jl",
         modules=[NeuroAnalyzer],
         authors="Adam Wysokiński",
         linkcheck=false,
         remotes=nothing,
         warnonly=true,
         clean=true,
         format=Documenter.HTML(size_threshold=268435456))
