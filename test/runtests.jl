using NeuroAnalyzer
using Test

@testset "runtests.jl" begin
    @info "Running statistics.jl tests"
    @test include("statistics.jl")

    @info "Running io.jl tests"
    @test include("io.jl")

    @info "Running plots.jl tests"
    @test include("plots.jl")

    @info "Running stim.jl tests"
    @test include("stim.jl")

    @info "Running study.jl tests"
    @test include("study.jl")

    @info "Running low_level.jl tests"
    @test include("low_level.jl")

    @info "Running eeg.jl tests"
    @test include("eeg.jl")
end