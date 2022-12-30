using Pkg
pwd()
Pkg.add(["Documenter", "DocumenterMarkdown", "Plots", "GLMakie", "DataFrames", "Wavelets", "ContinuousWavelets", "StatsModels", "MultivariateStats"])
Pkg.activate(@__DIR__)
Pkg.instantiate()