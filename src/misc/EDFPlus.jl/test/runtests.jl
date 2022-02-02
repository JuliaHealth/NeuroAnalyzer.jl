using EDFPlus
using Test
# Run tests

@test include("readtest.jl")
@test include("writetest.jl")
@test include("cover.jl")


