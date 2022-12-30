@info "Import packages"
using Pkg
Pkg.add(["Documenter", "DocumenterMarkdown", "Plots", "GLMakie", "DataFrames", "Wavelets", "ContinuousWavelets", "StatsModels", "MultivariateStats"])
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

@info "Clone NA-docs"
run(`git clone https://codeberg.org/AdamWysokinski/NeuroAnalyzer-docs NA-docs`);

## STABLE

@info "Clone NA-stable"
run(`git clone -b stable https://codeberg.org/AdamWysokinski/NeuroAnalyzer.jl NA-stable`);
@info "Go to NA-stable/docs/ folder"
cd("NA-stable/docs");
@info "Delete build/ folder"
run(`rm -rf build`);
@info "Delete index.md file"
run(`rm -f src/index.md`);
@info "Create new index.md file (1)"
run(`cat header-stable.md > src/index.md`);
@info "Create new index.md file (2)"
run(`./template.sh >> src/index.md`);
## MD
@info "Generate Markdown documentation"
makedocs(sitename="NeuroAnalyzer.jl", format=Markdown())
@info "Create Documentation-stable.md file"
run(`cp build/index.md ../../NA-docs/Documentation-stable.md`);
@info "Delete build/ folder"
run(`rm -rf build`);
## HTML
@info "Generate HTML documentation"
makedocs(sitename="NeuroAnalyzer.jl")
@info "Rename build => docs-stable"
run(`mv build docs-stable`);
@info "Replace GitHub => Codeberg"
run(`sed -i 's/Edit on GitHub/Edit on Codeberg/g' docs-stable/index.html`);
@info "Remove GitHub logo"
run(`sed -i 's///g' docs-stable/index.html`);
@info "Move docs-stable to outside NA-stable folder"
run(`mv docs-stable ../../`);
@info "Leave NA-stable/docs/ folder"
cd("../..");

## DEVEL

@info "Clone NA-devel"
run(`git clone -b devel https://codeberg.org/AdamWysokinski/NeuroAnalyzer.jl NA-devel`);
@info "Go to NA-devel/docs/ folder"
cd("NA-devel/docs");
@info "Delete build/ folder"
run(`rm -rf build`);
@info "Delete index.md file"
run(`rm -f src/index.md`);
@info "Create new index.md file (1)"
run(`cat header-devel.md > src/index.md`);
@info "Create new index.md file (2)"
run(`./template.sh >> src/index.md`);
## MD
@info "Generate Markdown documentation"
makedocs(sitename="NeuroAnalyzer.jl", format=Markdown())
@info "Create Documentation-devel.md file"
run(`cp build/index.md ../../NA-docs/Documentation-devel.md`);
@info "Delete build/ folder"
run(`rm -rf build`);
## HTML
@info "Generate HTML documentation"
makedocs(sitename="NeuroAnalyzer.jl")
@info "Rename build => docs-devel"
run(`mv build docs-devel`);
@info "Replace GitHub => Codeberg"
run(`sed -i 's/Edit on GitHub/Edit on Codeberg/g' docs-devel/index.html`);
@info "Remove GitHub logo"
run(`sed -i 's///g' docs-devel/index.html`);
@info "Move docs-devel to outside NA-devel folder"
run(`mv docs-devel ../../`);
@info "Leave NA-devel/docs/ folder"
cd("../..");