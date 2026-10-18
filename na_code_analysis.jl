#!/usr/bin/env julia

# Code analysis run before every stable release (na_create_release_STABLE.sh):
# - Aqua quality checks: a failure stops the release
# - JET static analysis: report only (many reports are not errors), saved to jet_report.txt
# - the test suite
#
# Aqua and JET must be installed in the default (global) environment.

using Pkg

Pkg.activate(@__DIR__)
Pkg.instantiate()
# test the release with the newest compatible dependencies (Manifest.toml is not committed,
# so these are the versions users will get)
Pkg.update()

using Aqua
using JET
using NeuroAnalyzer

@info "Running NeuroAnalyzer Aqua code analysis.."

# methods that deliberately extend functions of Base / StatsBase for Base types:
# cis(x::AbstractVector) (confidence interval of SD), summary(x1, x2, xs...) (multi-group
# summary), zscore(x::Real, m, sd) (z-score of a single value)
Aqua.test_all(
    NeuroAnalyzer;
    piracies = (treat_as_own = [Base.cis, Base.summary, NeuroAnalyzer.zscore],),
)

@info "Running NeuroAnalyzer JET code analysis.."

# the report of a whole package is long: it is saved to a file (ignored by git) for review
jet_report = JET.report_package(
    NeuroAnalyzer;
    target_modules = (NeuroAnalyzer,),
    toplevel_logger = nothing,
)
jet_file = joinpath(@__DIR__, "jet_report.txt")
open(io -> show(io, MIME"text/plain"(), jet_report), jet_file, "w")
@info "JET: $(length(JET.get_reports(jet_report))) possible error(s), see $jet_file"

@info "Running NeuroAnalyzer tests.."

Pkg.test()

@info "Done."
