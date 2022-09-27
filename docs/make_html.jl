@info "Loading package: Documenter"
using Documenter
@info "Loading package: Plots"
using Plots
@info "Loading package: GLMakie"
using GLMakie
@info "Loading package: DataFrames"
using DataFrames
@info "Loading package: NeuroAnalyzer"
using NeuroAnalyzer
println()

makedocs(sitename="NeuroAnalyzer.jl")