using Pkg
Pkg.add(["Documenter", "DocumenterMarkdown", "Plots", "GLMakie", "DataFrames", "Wavelets", "ContinuousWavelets", "StatsModels", "MultivariateStats"])
Pkg.add(url="https://codeberg.org/AdamWysokinski/NeuroAnalyzer.jl#devel")
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

## STABLE

## MD
@info "Generate Markdown documentation: STABLE branch"
cd("NA-stable/docs")
makedocs(sitename="NeuroAnalyzer.jl", format=Markdown())
println()
run(`pwd`)
run(`ls -la build/`)
@info "Create Documentation-stable.md file"
run(`cp build/index.md ../../NA-docs/Documentation-stable.md`);
@info "Delete build/ folder"
run(`rm -rf build`);
## HTML
@info "Generate HTML documentation: STABLE branch"
makedocs(sitename="NeuroAnalyzer.jl")
cd("../..")

## DEVEL

## MD
@info "Generate Markdown documentation: DEVEL branch"
cd("NA-devel/docs")
run(`rm -rf src/index.md`)
makedocs(sitename="NeuroAnalyzer.jl", format=Markdown())
@info "Create Documentation-devel.md file"
run(`cp build/index.md ../../NA-docs/Documentation-devel.md`);
@info "Delete build/ folder"
run(`rm -rf build`);
## HTML
@info "Generate HTML documentation: DEVEL branch"
makedocs(sitename="NeuroAnalyzer.jl")
cd("../..")