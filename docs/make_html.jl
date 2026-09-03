@info "Generating HTML documentation"

using Pkg
@info "Activating packages..."
Pkg.add(; url = "https://codeberg.org/AdamWysokinski/NeuroAnalyzer.jl")
Pkg.add("Documenter")
Pkg.instantiate()
using Documenter
using NeuroAnalyzer

makedocs(; sitename = "NeuroAnalyzer.jl",
    modules = [NeuroAnalyzer],
    authors = "Adam Wysokiński",
    linkcheck = false,
    remotes = nothing,
    warnonly = true,
    clean = true,
    format = Documenter.HTML(; size_threshold = 268435456))