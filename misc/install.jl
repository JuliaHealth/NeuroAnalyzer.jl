@info "Installing NeuroAnalyzer"
using Pkg
Pkg.update()
Pkg.add("Revise")
Pkg.add(url="https://codeberg.org/AdamWysokinski/NeuroAnalyzer.jl")
Pkg.instantiate()
Pkg.resolve()
Pkg.update()
@info "Starting NeuroAnalyzer"
using NeuroAnalyzer
