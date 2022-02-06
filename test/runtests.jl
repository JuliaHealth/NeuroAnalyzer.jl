using NeuroJ
using Test

print("Running misc.jl tests.. ")
@test include("misc.jl")
println("passed")

print("Running signal.jl tests.. ")
@test include("signal.jl")
println("passed")

print("Running eeg_io.jl tests.. ")
@test include("eeg_io.jl")
println("passed")

print("Running eeg.jl tests.. ")
@test include("eeg.jl")
println("passed")

print("Running plots.jl tests.. ")
@test include("plots.jl")
println("passed")