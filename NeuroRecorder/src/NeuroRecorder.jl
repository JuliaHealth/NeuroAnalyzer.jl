"""
Record neurophysiological data with Julia.

https://neuroanalyzer.org
"""
module NeuroRecorder

@assert VERSION >= v"1.10.0" "NeuroRecorder requires Julia 1.10.0 or above."

# set constants

const VER = v"0.24.4"
const allow_wip = occursin("dev", string(VER))  # false for the stable branch, true for the devel branch

# initialize preferences

verbose = true

# add dependencies

using Artifacts
using Cairo
using DataFrames
using Dates
using Gtk
using LibSerialPort
using PiGPIO
using Plots
using REPL
using Statistics
using WAV

# initialize

function __init__()

    global res_path = joinpath(artifact"NeuroRecorder_resources", "neurorecorder-resources")

end

# load sub-modules

# internal functions
include("internal/misc.jl")
include("internal/rt_plotter.jl")
# record
include("recorder/ftt.jl")

end # NeuroRecorder
