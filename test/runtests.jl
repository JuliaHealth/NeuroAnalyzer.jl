using NeuroJ
using Test

# Run tests
@test include("io_eeg_import_edf.jl")
@test include("misc.jl")
