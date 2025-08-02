@info "Generating HTML documentation"

using Pkg
Pkg.add(url="https://codeberg.org/AdamWysokinski/FIRLSFilterDesign.jl")
Pkg.add(url="https://codeberg.org/AdamWysokinski/NeuroAnalyzer.jl.git")
Pkg.instantiate()
using Documenter

makedocs(sitename="NeuroAnalyzer.jl",
         modules=[NeuroAnalyzer],
         authors="Adam Wysokiński",
         linkcheck=false,
         remotes=nothing,
         warnonly=true,
         clean=true,
         format=Documenter.HTML(size_threshold=268435456))
