using NeuroJ
using Test

@testset "runtests.jl" begin
    println("\tRunning misc.jl tests.. ")
    @test include("misc.jl")

    println("\tRunning signal.jl tests.. ")
    @test include("signal.jl")

    println("\tRunning eeg_io.jl tests.. ")
    @test include("eeg_io.jl")

    println("\tRunning eeg.jl tests.. ")
    @test include("eeg.jl")

    println("\tRunning plots.jl tests.. ")
    @test include("plots.jl")
end