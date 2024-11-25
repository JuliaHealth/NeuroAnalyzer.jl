"""
Julia toolbox for analyzing neurophysiological data.

https://neuroanalyzer.org
"""
module NeuroAnalyzer

@assert VERSION >= v"1.10.0" "NeuroAnalyzer requires Julia 1.10.0 or above."

# set constants

const VER = v"0.24.12-dev"
const allow_wip = occursin("dev", string(VER))          # false for the stable branch, true for the devel branch
const io = PipeBuffer()                                 # required for interactive preview
const data_types = ["eeg",
                    "ecog",
                    "seeg",
                    "ieeg",
                    "csd",
                    "meg",
                    "nirs",
                    "sensors",
                    "eda",
                    "mep",
                    "erp",
                    "erf"]
const channel_types = ["all",
                       "eeg", "ecog", "seeg", "ieeg",
                       "csd",
                       "meg", "grad", "mag",
                       "csd",
                       "nirs", "nirs_int", "nirs_od", "nirs_dmean", "nirs_dvar", "nirs_dskew", "nirs_mua", "nirs_musp", "nirs_hbo", "nirs_hbr", "nirs_hbt", "nirs_h2o", "nirs_lipid", "nirs_bfi", "nirs_hrf_dod", "nirs_hrf_dmean", "nirs_hrf_dvar", "nirs_hrf_dskew", "nirs_hrf_hbo", "nirs_hrf_hbr", "nirs_hrf_hbt", "nirs_hrf_bfi", "nirs_aux",
                       "ecg",
                       "emg",
                       "eog",
                       "ref",
                       "mrk",
                       "sensors", "accel", "magfld", "orient", "angvel",
                       "mep",
                       "eda",
                       "other"]
const channel_units = ["μV",
                       "mV",
                       "V",
                       "μV/m²",
                       "fT",
                       "fT/cm",
                       "μM/mm",
                       "m/s²",
                       "µT",
                       "°",
                       "μS",
                       "rad/s",
                       ""]
const fiducial_points = (nasion = (0.0, 0.95, -0.2),
                         inion  = (0.0, -0.96, -0.2),
                         lpa    = (-0.98, 0.0, -0.2),
                         rpa    = (0.98, 0.0, -0.2))
begin
    tmp = pwd()
    cd(joinpath(dirname(pathof(NeuroAnalyzer)), ".."))
    global const PATH = pwd()
    cd(tmp)
end

# initialize preferences

use_cuda = nothing
progress_bar = nothing
verbose = nothing

# add dependencies

using Artifacts
using Cairo
using ColorSchemes
using ContinuousWavelets
using CSV
using CubicSplines
using CUDA
using DataFrames
using Dates
using Deconvolution
using DICOM
using Distances
using DSP
using FFTW
using FileIO
using FindPeaks1D
using FIRLSFilterDesign
using FourierTools
using GeometryBasics
using Git
using GLM
using Gtk
using HypothesisTests
using Images
using ImageBinarization
using ImageFiltering
using ImageMorphology
using InformationMeasures
using Interpolations
using Jacobi
using JLD2
using JSON
using LibSerialPort
using LinearAlgebra
using LinRegOutliers
using Loess
using Logging
using MAT
using MLJ
using MultivariateStats
using NPZ
using PiGPIO
using Pkg
using Plots
using Plots.PlotMeasures
using Polynomials
using Preferences
using PrettyTables
using ProgressMeter
using Random
using REPL
using SavitzkyGolay
using ScatteredInterpolation
using Simpson
using StatsFuns
using StatsKit
using StatsModels
using StatsPlots
using TimeZones
using TOML
using WAV
using Wavelets
using WaveletsExt
using XDF

# define structures

mutable struct HEADER
    subject::Dict
    recording::Dict
    experiment::Dict
end

mutable struct NEURO
    header::NeuroAnalyzer.HEADER
    time_pts::Vector{Float64}
    epoch_time::Vector{Float64}
    data::Array{<:Number, 3}
    components::Dict{Any}
    markers::DataFrame
    locs::DataFrame
    history::Vector{String}
end

mutable struct STUDY
    header::Dict{Symbol, Any}
    objects::Vector{NeuroAnalyzer.NEURO}
    groups::Vector{Symbol}
end

mutable struct DIPOLE
    pos::Tuple{Real, Real, Real}
    mag::Tuple{Real, Real, Real}
end

# set package options

Plots.gr_cbar_width[] = 0.01
if Sys.islinux() && Sys.ARCH === :x86_64
    FFTW.set_provider!("mkl")
else
    FFTW.set_provider!("fftw")
end
FFTW.set_num_threads(Sys.CPU_THREADS)
BLAS.set_num_threads(Sys.CPU_THREADS)

# load NA functions

include("na/internal.jl")
include("na/setup.jl")
include("na/plugins.jl")

# load preferences

global use_cuda = @load_preference("use_cuda", false)
global progress_bar = @load_preference("progress_bar", true)
global verbose = @load_preference("verbose", true)
na_set_prefs(use_cuda=use_cuda, progress_bar=progress_bar, verbose=verbose)

# be verbose

_info("NeuroAnalyzer v$(NeuroAnalyzer.VER)")
_info("NeuroAnalyzer path: $(NeuroAnalyzer.PATH)")
_info("Preferences loaded:")
_info(" Use CUDA: $use_cuda")
_info(" Progress bar: $progress_bar")
_info(" Verbose: $verbose")

# setup resources

_info("Preparing resources")
global res_path = joinpath(artifact"NeuroAnalyzer_resources", "neuroanalyzer-resources")

# load plugins

global plugins_path = joinpath(homedir(), "NeuroAnalyzer", "plugins")
if isdir(plugins_path)
    if length(readdir(plugins_path)) > 0
        _info("Loading plugins:")
        na_plugins_reload()
    end
else
    mkpath(plugins_path)
end

# load sub-modules

_load_functions("internal")
_load_functions("utils")
_load_functions("io")
_load_functions("locs")
_load_functions("edit")
_load_functions("process")
_load_functions("analyze")
_load_functions("plots")
_load_functions("gui")
_load_functions("statistics")
_load_functions("study")
_load_functions("recorder")
_load_functions("stim")

end # NeuroAnalyzer
