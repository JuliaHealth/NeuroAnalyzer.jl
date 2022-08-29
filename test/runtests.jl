using NeuroAnalyzer
using Test

@testset "runtests.jl" begin
    #println("\tRunning low_level.jl tests.. ")
    #@test include("low_level.jl")

    #println("\tRunning statistic.jl tests.. ")
    #@test include("statistic.jl")

    #println("\tRunning eeg_io.jl tests.. ")
    #@test include("eeg_io.jl")

    println("\tRunning eeg.jl tests.. ")
    @test include("eeg.jl")

    println("\tRunning eeg_plots.jl tests.. ")
    @test include("eeg_plots.jl")
end