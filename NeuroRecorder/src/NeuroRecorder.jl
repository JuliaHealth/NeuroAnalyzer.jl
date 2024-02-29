"""
Record neurophysiological data with Julia.

https://neuroanalyzer.org
"""
module NeuroRecorder

@assert VERSION >= v"1.10.0" "This version of NeuroAnalyzer requires Julia 1.10.0 or above."

# set constants

const VER = v"0.24.4-dev"
const allow_wip = occursin("dev", string(VER))  # false for the stable branch, true for the devel branch

# add dependencies

using Artifacts
using Cairo
using Dates
using Gtk
using WAV
using PiGPIO
using REPL

# initialize

function __init__()

    global res_path = joinpath(artifact"NeuroRecorder_resources", "neurorecorder-resources")

end

# load sub-modules

include("internal.jl")
include("ftt.jl")

end
