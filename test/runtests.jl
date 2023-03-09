using NeuroAnalyzer
using Test

@testset "runtests.jl" begin
    println("\tRunning io.jl tests.. ")
    @test include("io.jl")

    println("\tRunning low_level.jl tests.. ")
    @test include("low_level.jl")

    println("\tRunning statistics.jl tests.. ")
    @test include("statistics.jl")

    println("\tRunning eeg.jl tests.. ")
    @test include("eeg.jl")

    println("\tRunning eeg_plots.jl tests.. ")
    @test include("eeg_plots.jl")

    println("\tRunning eeg_study.jl tests.. ")
    @test include("eeg_study.jl")

    println("\tRunning tes.jl tests.. ")
    @test include("tes.jl")
end