println()
println("Generating NeuroJ Markdown documentation..")
println()
println("Loading package: Documenter")
using Documenter
println("Loading package: DocumenterMarkdown")
using DocumenterMarkdown
println("Loading package: NeuroJ")
using NeuroJ
println("Loading package: Plots")
using Plots
println("Loading package: GLMakie")
using GLMakie
println()

makedocs(sitename="NeuroJ.jl", format=Markdown())