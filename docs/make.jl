println()
println("Generating NeuroJ Markdown documentation..")
println()
@info "Loading package: Documenter"
using Documenter
@info "Loading package: DocumenterMarkdown"
using DocumenterMarkdown
@info "Loading package: NeuroJ"
using NeuroJ
@info "Loading package: Plots"
using Plots
@info "Loading package: GLMakie"
using GLMakie
println()

makedocs(sitename="NeuroJ.jl", format=Markdown())