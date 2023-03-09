__precompile__()

module NeuroAnalyzer

if VERSION < v"1.7.0"
    @error("This version of NeuroAnalyzer requires Julia 1.7.0 or above.")
end

const na_ver = v"0.23.3"

# initialize preferences
use_cuda = nothing
progress_bar = nothing
plugins_path = nothing
verbose = nothing

using ColorSchemes
using CSV
using CubicSplines
using CUDA
using Dates
using DataFrames
using Deconvolution
using DICOM
using Distances
using DSP
using FFTW
using FileIO
using FindPeaks1D
using FourierTools
using GeometryBasics
using Git
using GLM
using GLMakie
using HypothesisTests
using InformationMeasures
using Interpolations
using Jacobi
using JLD2
using LinearAlgebra
using Loess
using MAT
using MultivariateStats
using Pkg
using Plots
using Plots.PlotMeasures
using Polynomials
using Preferences
using ProgressMeter
using Random
using SavitzkyGolay
using ScatteredInterpolation
using Simpson
using StatsFuns
using StatsKit
using StatsModels
using StatsPlots
using TOML
using Wavelets
using WaveletsExt
using ContinuousWavelets

mutable struct HEADER
    subject::Dict
    recording::Dict
    experiment::Dict
    markers::Bool
    components::Vector{Symbol}
    locations::Bool
    history::Vector{String}
end

mutable struct NEURO
    header::NeuroAnalyzer.HEADER
    time_pts::Vector{Float64}
    epoch_time::Vector{Float64}
    data::Union{Array{<:Number, 1}, Array{<:Number, 2}, Array{<:Number, 3}}
    components::Vector{Any}
    events::DataFrame
    locs::DataFrame
end

mutable struct EEG
    eeg_header::Dict
    eeg_time::Vector{Float64}
    eeg_epoch_time::Vector{Float64}
    eeg_signals::Array{Float64, 3}
    eeg_components::Vector{Any}
    eeg_markers::DataFrame
    eeg_locs::DataFrame
end

mutable struct MEG
    meg_header::Dict
    meg_time::Vector{Float64}
    meg_epoch_time::Vector{Float64}
    meg_signals::Array{Float64, 3}
    meg_components::Vector{Any}
    meg_markers::DataFrame
    meg_locs::DataFrame
end

mutable struct STUDY
    study_header::Dict{Symbol, Any}
    study_eeg::Vector{NeuroAnalyzer.NEURO}
    study_group::Vector{Symbol}
end

mutable struct DIPOLE
    loc::Tuple{Real, Real, Real}
end

FFTW.set_provider!("mkl")
FFTW.set_num_threads(Sys.CPU_THREADS)
BLAS.set_num_threads(Sys.CPU_THREADS)

# NA functions
include("utils/na.jl")

function __init__()

    @info "NeuroAnalyzer v$na_ver"

    # load preferences
    @info "Loading preferences..."
    if Sys.isunix() || Sys.isapple()
        def_plugins_path = "$(homedir())/NeuroAnalyzer/plugins/"
    elseif Sys.iswindows()
        def_plugins_path = "$(homedir())\\NeuroAnalyzer\\plugins\\"
    end
    global use_cuda = @load_preference("use_cuda", false)
    global progress_bar = @load_preference("progress_bar", true)
    global plugins_path = @load_preference("plugins_path", def_plugins_path)
    global verbose = @load_preference("verbose", true)
    na_set_prefs(use_cuda=use_cuda, plugins_path=plugins_path, progress_bar=progress_bar, verbose=verbose)

    # load plugins
    @info "Loading plugins..."
    isdir(plugins_path) || mkdir(plugins_path)
    na_plugins_reload()

end

# load sub-modules
# find . -name "*.jl"
# internal functions are exposed outside NA
include("internal/ch_idx.jl")
include("internal/check.jl")
include("internal/create_header.jl")
include("internal/draw_head.jl")
include("internal/fiff.jl")
include("internal/fir_response.jl")
include("internal/gpu.jl")
include("internal/interpolate.jl")
include("internal/io.jl")
include("internal/labeled_matrix.jl")
include("internal/labels.jl")
include("internal/len.jl")
include("internal/locs.jl")
include("internal/make_epochs.jl")
include("internal/map_channels.jl")
include("internal/markers.jl")
include("internal/misc.jl")
include("internal/ml.jl")
include("internal/plots.jl")
include("internal/reflect_chop.jl")
include("internal/select.jl")
include("internal/tester.jl")
include("internal/time.jl")
# analyze
include("analyze/acov.jl")
include("analyze/band_power.jl")
include("analyze/corm.jl")
include("analyze/covm.jl")
include("analyze/cps.jl")
include("analyze/detect_bad.jl")
include("analyze/difference.jl")
include("analyze/dissimilarity.jl")
include("analyze/entropy.jl")
include("analyze/envelopes.jl")
include("analyze/erp.jl")
include("analyze/fcoherence.jl")
include("analyze/ged.jl")
include("analyze/ica.jl")
include("analyze/ispc.jl")
include("analyze/itpc.jl")
include("analyze/mi.jl")
include("analyze/msci95.jl")
include("analyze/pca.jl")
include("analyze/pli.jl")
include("analyze/psd.jl")
include("analyze/psd_mw.jl")
include("analyze/psd_rel.jl")
include("analyze/psdslope.jl")
include("analyze/snr_rms.jl")
include("analyze/spectrogram.jl")
include("analyze/spectrum.jl")
include("analyze/stationarity.jl")
include("analyze/stats.jl")
include("analyze/tcoherence.jl")
include("analyze/tkeo.jl")
include("analyze/total_power.jl")
include("analyze/vartest.jl")
include("analyze/xcov.jl")
# edit
include("edit/channel.jl")
include("edit/delete_channel.jl")
include("edit/delete_epoch.jl")
include("edit/epoch.jl")
include("edit/extract.jl")
include("edit/marker.jl")
include("edit/reflect_chop.jl")
include("edit/trim.jl")
include("edit/vch.jl")
# io
include("io/export_csv.jl")
include("io/fiff.jl")
include("io/import.jl")
include("io/import_alice4.jl")
include("io/import_bdf.jl")
include("io/import_bv.jl")
include("io/import_csv.jl")
include("io/import_digitrack.jl")
include("io/import_edf.jl")
include("io/import_fiff.jl")
include("io/import_set.jl")
include("io/load_locs.jl")
include("io/locs_export.jl")
include("io/locs_import.jl")
include("io/save_load.jl")
# locs
include("locs/add_locs.jl")
include("locs/convert.jl")
include("locs/details.jl")
include("locs/edit.jl")
include("locs/flip.jl")
include("locs/import.jl")
include("locs/rotate.jl")
include("locs/scale.jl")
include("locs/swap.jl")
#process
include("process/add_noise.jl")
include("process/average.jl")
include("process/cbp.jl")
include("process/cwt.jl")
include("process/demean.jl")
include("process/denoise_fft.jl")
include("process/denoise_wavelet.jl")
include("process/denoise_wien.jl")
include("process/derivative.jl")
include("process/detrend.jl")
include("process/dwt.jl")
include("process/dwtsplit.jl")
include("process/erp.jl")
include("process/fbsplit.jl")
include("process/fconv.jl")
include("process/filter.jl")
include("process/gfilter.jl")
include("process/invert.jl")
include("process/lrinterpolate.jl")
include("process/normalize.jl")
include("process/plinterpolate.jl")
include("process/reference.jl")
include("process/resample.jl")
include("process/scale.jl")
include("process/slaplacian.jl")
include("process/standardize.jl")
include("process/taper.jl")
include("process/tconv.jl")
include("process/wbp.jl")
include("process/zero.jl")
# plot
include("plots/misc.jl")
include("plots/plot_connections.jl")
include("plots/plot_dipole3d.jl")
include("plots/plot_electrodes.jl")
include("plots/plot_erp.jl")
include("plots/plot_filter_response.jl")
include("plots/plot_psd.jl")
include("plots/plot_save.jl")
include("plots/plot_signal.jl")
include("plots/plot_spectrogram.jl")
include("plots/plot_topo.jl")
include("plots/plot_varia.jl")
include("plots/plot_weights.jl")
# statistics
include("statistics/dprime.jl")
include("statistics/effsize.jl")
include("statistics/hildebrand_rule.jl")
include("statistics/jaccard_similarity.jl")
include("statistics/linreg.jl")
include("statistics/means.jl")
include("statistics/misc.jl")
include("statistics/ml.jl")
include("statistics/norminv.jl")
include("statistics/outliers.jl")
include("statistics/pred_int.jl")
include("statistics/ranks.jl")
include("statistics/res_norm.jl")
include("statistics/s2cmp.jl")
include("statistics/s2cor.jl")
include("statistics/segments.jl")
include("statistics/sem_diff.jl")
# utils
include("utils/apply.jl")
include("utils/array.jl")
include("utils/components.jl")
include("utils/fft.jl")
include("utils/findpeaks.jl")
include("utils/frequency.jl")
include("utils/generate.jl")
include("utils/header.jl")
include("utils/info.jl")
include("utils/locs_convert.jl")
include("utils/matrix.jl")
include("utils/misc.jl")
include("utils/note.jl")
include("utils/pad.jl")
include("utils/phase.jl")
include("utils/pick.jl")
include("utils/time.jl")
include("utils/trim.jl")
include("utils/vector.jl")
# study
include("study/create.jl")
include("study/info.jl")
# stim
include("stim/tes.jl")
include("stim/ect.jl")

end # NeuroAnalyzer