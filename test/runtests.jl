using NeuroAnalyzer
using Test

@testset "runtests.jl" begin
    println("\tRunning statistics.jl tests.. ")
    @test include("statistics.jl")

    println("\tRunning io.jl tests.. ")
    @test include("io.jl")

    println("\tRunning plots.jl tests.. ")
    @test include("plots.jl")

    println("\tRunning stim.jl tests.. ")
    @test include("stim.jl")

    println("\tRunning study.jl tests.. ")
    @test include("study.jl")

    println("\tRunning low_level.jl tests.. ")
    @test include("low_level.jl")

    println("\tRunning eeg.jl tests.. ")
    @test include("eeg.jl")
end