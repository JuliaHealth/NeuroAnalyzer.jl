"""
Design and run psychological studies with Julia.

https://neuroanalyzer.org
"""
module NeuroTester

@assert VERSION >= v"1.10.0" "NeuroTester requires Julia 1.10.0 or above."

# set constants

const VER = v"0.24.4-dev"
const allow_wip = occursin("dev", string(VER))  # false for the stable branch, true for the devel branch

# add dependencies

# initialize

function __init__()


end

# load sub-modules

include("internal/misc.jl")

end # NeuroTester
