#!/usr/bin/env julia

using Pkg

Pkg.activate("..")

using Aqua
using NeuroAnalyzer
using JET

@info "Running NeuroAnalyzer Aqua code analysis.."

Aqua.test_all(NeuroAnalyzer)

@info "Running NeuroAnalyzer JET code analysis.."

JET.report_file("../src/NeuroAnalyzer.jl")

@info "Done."
