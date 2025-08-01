@info "Generating HTML documentation"

using Documenter
using NeuroAnalyzer

#=
using Pkg
Pkg.add(["Documenter", "Plots", "DataFrames", "Wavelets", "ContinuousWavelets", "StatsModels", "MultivariateStats"])
Pkg.activate("..")
Pkg.add(url="https://codeberg.org/AdamWysokinski/FIRLSFilterDesign.jl")
Pkg.instantiate()

@info "Loading package: Documenter"
using Documenter
@info "Loading package: Plots"
using Plots
@info "Loading package: Cairo"
using Cairo
@info "Loading package: DataFrames"
using DataFrames
@info "Loading package: Wavelets"
using Wavelets
@info "Loading package: ContinuousWavelets"
using ContinuousWavelets
@info "Loading package: StatsModels"
using StatsModels
@info "Loading package: MultivariateStats"
using MultivariateStats
@info "Loading package: NeuroAnalyzer"
using NeuroAnalyzer
println()
=#

makedocs(sitename="NeuroAnalyzer.jl",
         modules=[NeuroAnalyzer],
         authors="Adam Wysokiński",
         linkcheck=false,
         remotes=nothing,
         warnonly=true,
         clean=true,
         format=Documenter.HTML(size_threshold=268435456))
