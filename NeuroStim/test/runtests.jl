using NeuroStim
using Test

@testset "NeuroStim.jl" begin

    @info "Running stim.jl tests"
    @test include("stim.jl")

end
