#!/usr/bin/env julia

using JET

@info "Running NeuroAnalyzer code analysis.."

report_file("../src/NeuroAnalyzer.jl"; annotate_types = true)

@info "Done."