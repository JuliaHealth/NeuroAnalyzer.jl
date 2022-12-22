using Pkg
Pkg.add("Documenter")
Pkg.add("DocumenterMarkdown")
Pkg.add("Plots")
Pkg.add("GLMakie")
Pkg.add("DataFrames")
Pkg.add("Wavelets")
Pkg.add("ContinuousWavelets")
Pkg.add("StatsModels")
Pkg.add("MultivariateStats")
Pkg.activate(@__DIR__)
Pkg.instantiate()

@info "Loading package: Documenter"
using Documenter
@info "Loading package: DocumenterMarkdown"
using DocumenterMarkdown
@info "Loading package: Plots"
using Plots
@info "Loading package: GLMakie"
using GLMakie
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

makedocs(sitename="NeuroAnalyzer.jl", format=Markdown())