#!/usr/bin/env julia

using Pkg

Pkg.activate(@__DIR__)
Pkg.instantiate()
Pkg.resolve()
Pkg.update()

using Aqua
using NeuroAnalyzer
using JET

@info "Running NeuroAnalyzer Aqua code analysis.."

Aqua.test_all(NeuroAnalyzer)

@info "Running NeuroAnalyzer JET code analysis.."

JET.report_file("src/NeuroAnalyzer.jl")

@info "Running NeuroAnalyzer tests.."

Pkg.test()

@info "Done."
