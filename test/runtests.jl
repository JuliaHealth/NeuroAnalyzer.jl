using NeuroAnalyzer
using Test
using Artifacts

global testfiles_path = joinpath(artifact"NeuroAnalyzer_test-files", "neuroanalyzer-test-files")

@testset "NeuroAnalyzer.jl" failfast=true begin

    @info "Running internal.jl tests"
    @test include("internal.jl")

    @info "Running io.jl tests"
    @test include("io.jl")

    @info "Running utils.jl tests"
    @test include("utils.jl")

    @info "Running stats.jl tests"
    @test include("stats.jl")

    @info "Running locs.jl tests"
    @test include("locs.jl")

    @info "Running edit.jl tests"
    @test include("edit.jl")

    @info "Running process.jl tests"
    @test include("process.jl")

    @info "Running analyze.jl tests"
    @test include("analyze.jl")

    @info "Running plots.jl tests"
    @test include("plots.jl")

    @info "Running study.jl tests"
    @test include("study.jl")

    @info "Running stim.jl tests"
    @test include("stim.jl")

end
