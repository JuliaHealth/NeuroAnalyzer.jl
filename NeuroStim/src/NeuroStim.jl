"""
Process and model neurostimulation techniques with Julia.

https://neuroanalyzer.org
"""
module NeuroStim

@assert VERSION >= v"1.10.0" "NeuroStim requires Julia 1.10.0 or above."

# set constants

const VER = v"0.24.4-dev"
const allow_wip = occursin("dev", string(VER))  # false for the stable branch, true for the devel branch

# add dependencies

# initialize

function __init__()


end

# load sub-modules

# internal functions
include("internal/check.jl")
# stim
include("stim/ect.jl")
include("stim/tes.jl")
include("stim/tes_model.jl")

end # NeuroStim
