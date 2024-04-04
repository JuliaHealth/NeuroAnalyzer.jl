"""
Neurophysiological data analysis with Julia.

https://neuroanalyzer.org
"""
module NeuroAnalyzer

@assert VERSION >= v"1.10.0" "NeuroAnalyzer requires Julia 1.10.0 or above."

# set constants

const VER = v"0.24.4"
const allow_wip = occursin("dev", string(VER))  # false for the stable branch, true for the devel branch
const io = PipeBuffer()                         # required for interactive preview
const data_types = ["eeg",
                    "ecog",
                    "seeg",
                    "ieeg",
                    "meg",
                    "nirs",
                    "sensors",
                    "eda",
                    "mep"]
const channel_types = ["all",
                       "eeg", "ecog", "seeg", "ieeg",
                       "meg", "grad", "mag",
                       "csd",
                       "nirs_int", "nirs_od", "nirs_dmean", "nirs_dvar", "nirs_dskew", "nirs_mua", "nirs_musp", "nirs_hbo", "nirs_hbr", "nirs_hbt", "nirs_h2o", "nirs_lipid", "nirs_bfi", "nirs_hrf_dod", "nirs_hrf_dmean", "nirs_hrf_dvar", "nirs_hrf_dskew", "nirs_hrf_hbo", "nirs_hrf_hbr", "nirs_hrf_hbt", "nirs_hrf_bfi", "nirs_aux",
                       "ecg",
                       "emg",
                       "eog",
                       "ref",
                       "mrk",
                       "accel", "magfld", "orient", "angvel",
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
using LinearAlgebra
using Loess
using MAT
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

include("na.jl")

function __init__()

    global use_cuda = @load_preference("use_cuda", false)
    global progress_bar = @load_preference("progress_bar", true)
    global verbose = @load_preference("verbose", true)
    na_set_prefs(use_cuda=use_cuda, progress_bar=progress_bar, verbose=verbose)

    _info("NeuroAnalyzer v$(NeuroAnalyzer.VER)")
    _info("NeuroAnalyzer path: $(NeuroAnalyzer.PATH)")

    # load preferences)
    _info("Preferences loaded:")
    _info(" Use CUDA: $use_cuda")
    _info(" Progress bar: $progress_bar")
    _info(" Verbose: $verbose")

    # setup resources
    _info("Preparing resources")
    global res_path = joinpath(artifact"NeuroAnalyzer_resources", "neuroanalyzer-resources")

    # load plugins
    _info("Loading plugins:")
    global plugins_path = joinpath(homedir(), "NeuroAnalyzer", "plugins")
    isdir(plugins_path) || mkpath(plugins_path)
    na_plugins_reload()

end

# load sub-modules

# internal functions
include("internal/check.jl")
include("internal/components.jl")
include("internal/create_header.jl")
include("internal/draw_head.jl")
include("internal/fiff.jl")
include("internal/fir_response.jl")
include("internal/gpu.jl")
include("internal/interpolate.jl")
include("internal/labeled_matrix.jl")
include("internal/labels.jl")
include("internal/len.jl")
include("internal/locs.jl")
include("internal/gradient.jl")
include("internal/epochs.jl")
include("internal/channels.jl")
include("internal/markers.jl")
include("internal/misc.jl")
include("internal/ml.jl")
include("internal/plots.jl")
include("internal/reflect_chop.jl")
include("internal/select.jl")
include("internal/time.jl")
include("internal/wl2ext.jl")
include("internal/gdf_etp.jl")
include("internal/statistics.jl")
include("internal/recorder.jl")                      
include("internal/rt_plotter.jl")
include("internal/tester.jl")
# analyze
include("analyze/ampdiff.jl")
include("analyze/amp.jl")
include("analyze/rmse.jl")
include("analyze/snr.jl")
include("analyze/erp_peaks.jl")
include("analyze/band_asymmetry.jl")
include("analyze/band_mpower.jl")
include("analyze/band_power.jl")
include("analyze/coherence.jl")
include("analyze/cpsd.jl")
include("analyze/erop.jl")
include("analyze/eros.jl")
include("analyze/frqinst.jl")
include("analyze/hjorth.jl")
include("analyze/peak_frq.jl")
include("analyze/psd.jl")
include("analyze/psd_rel.jl")
include("analyze/psd_slope.jl")
include("analyze/sef.jl")
include("analyze/spectrogram.jl")
include("analyze/spectrum.jl")
include("analyze/total_power.jl")
include("analyze/mep_peaks.jl")
include("analyze/corm.jl")
include("analyze/covm.jl")
include("analyze/dissimilarity.jl")
include("analyze/entropy.jl")
include("analyze/envelopes.jl")
include("analyze/ged.jl")
include("analyze/hrv.jl")
include("analyze/mi.jl")
include("analyze/stationarity.jl")
include("analyze/tkeo.jl")
include("analyze/cph.jl")
include("analyze/ispc.jl")
include("analyze/itpc.jl")
include("analyze/phdiff.jl")
include("analyze/phsd.jl")
include("analyze/pli.jl")
include("analyze/acor.jl")
include("analyze/acov.jl")
include("analyze/axc2frq.jl")
include("analyze/pacor.jl")
include("analyze/xcor.jl")
include("analyze/xcov.jl")
# edit
include("edit/channel.jl")
include("edit/create.jl")
include("edit/delete_channel.jl")
include("edit/delete_optode.jl")
include("edit/delete_epoch.jl")
include("edit/detect_bad.jl")
include("edit/epoch.jl")
include("edit/extract.jl")
include("edit/join.jl")
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
include("io/import_cnt.jl")
include("io/import_csv.jl")
include("io/import_dat.jl")
include("io/import_digitrack.jl")
include("io/import_edf.jl")
include("io/import_fiff.jl")
include("io/import_set.jl")
include("io/load_locs.jl")
include("io/export_locs.jl")
include("io/import_locs.jl")
include("io/save_load.jl")
include("io/import_snirf.jl")
include("io/import_nirs.jl")
include("io/import_nirx.jl")
include("io/import_edf_annotations.jl")
include("io/export_markers.jl")
include("io/import_gdf.jl")
include("io/import_montage.jl")
include("io/import_nwb.jl")
include("io/import_npy.jl")
include("io/import_thymatron.jl")
include("io/import_duomag.jl")
include("io/import_ncs.jl")
include("io/import_xdf.jl")
# locs
include("locs/add_locs.jl")
include("locs/convert.jl")
include("locs/details.jl")
include("locs/edit.jl")
include("locs/flip.jl")
include("locs/rotate.jl")
include("locs/scale.jl")
include("locs/swap.jl")
include("locs/generate.jl")
include("locs/center.jl")
#process
include("process/add_signal.jl")
include("process/average.jl")
include("process/bpsplit.jl")
include("process/cbp.jl")
include("process/ch_zero.jl")
include("process/csd.jl")
include("process/cwt.jl")
include("process/denoise_fft.jl")
include("process/denoise_wavelet.jl")
include("process/denoise_wien.jl")
include("process/derivative.jl")
include("process/detrend.jl")
include("process/dwt.jl")
include("process/dwtsplit.jl")
include("process/edit_montage.jl")
include("process/erp.jl")
include("process/fconv.jl")
include("process/filter.jl")
include("process/filter_g.jl")
include("process/filter_mavg.jl")
include("process/filter_mmed.jl")
include("process/filter_poly.jl")
include("process/filter_sg.jl")
include("process/ica.jl")
include("process/intensity2od.jl")
include("process/invert.jl")
include("process/lrinterpolate.jl")
include("process/normalize.jl")
include("process/normpower.jl")
include("process/npl.jl")
include("process/od2conc.jl")
include("process/pca.jl")
include("process/plinterpolate.jl")
include("process/reference.jl")
include("process/remove_dc.jl")
include("process/remove_pops.jl")
include("process/remove_powerline.jl")
include("process/resample.jl")
include("process/scale.jl")
include("process/standardize.jl")
include("process/taper.jl")
include("process/tconv.jl")
include("process/wbp.jl")
# plot
include("plots/misc.jl")
include("plots/plot_connections.jl")
include("plots/plot_locs.jl")
include("plots/plot_erp.jl")
include("plots/plot_filter_response.jl")
include("plots/plot_mep.jl")
include("plots/plot_psd.jl")
include("plots/plot_save.jl")
include("plots/plot_signal.jl")
include("plots/plot_spectrogram.jl")
include("plots/plot_topo.jl")
include("plots/plot_varia.jl")
include("plots/plot_weights.jl")
include("plots/plot_dipole2d.jl")
include("plots/plot_dipole3d.jl")
include("plots/plot_locs_nirs.jl")
include("plots/cairo.jl")
include("plots/plot_phsd.jl")
include("plots/plot_coherence.jl")
# gui
include("gui/iview.jl")
include("gui/iedit.jl")
include("gui/iedit_ch.jl")
include("gui/iplot.jl")
include("gui/iplot_locs3d.jl")
include("gui/iselect_ts.jl")
include("gui/itopo.jl")
include("gui/ipsd.jl")
include("gui/ispectrogram.jl")
include("gui/iplot_icatopo.jl")
include("gui/iview_plot.jl")
include("gui/iselect_seg.jl")
# statistics
include("statistics/specseg.jl")
include("statistics/stats.jl")
include("statistics/mdiff.jl")
include("statistics/msci95.jl")
include("statistics/dprime.jl")
include("statistics/effsize.jl")
include("statistics/friedman.jl")
include("statistics/hildebrand_rule.jl")
include("statistics/jaccard_similarity.jl")
include("statistics/linreg.jl")
include("statistics/means.jl")
include("statistics/misc.jl")
include("statistics/ml.jl")
include("statistics/norminv.jl")
include("statistics/outliers.jl")
include("statistics/power.jl")
include("statistics/pred_int.jl")
include("statistics/ranks.jl")
include("statistics/res_norm.jl")
include("statistics/cmp_test.jl")
include("statistics/cor_test.jl")
include("statistics/segments.jl")
include("statistics/sem_diff.jl")
include("statistics/summary.jl")
include("statistics/bootstrap.jl")
include("statistics/vartest.jl")
# utils
include("utils/apply.jl")
include("utils/array.jl")
include("utils/components.jl")
include("utils/fft.jl")
include("utils/findpeaks.jl")
include("utils/frequency.jl")
include("utils/generate.jl")
include("utils/gradient.jl")
include("utils/view_header.jl")
include("utils/info.jl")
include("utils/matrix.jl")
include("utils/misc.jl")
include("utils/note.jl")
include("utils/pad.jl")
include("utils/phase.jl")
include("utils/time.jl")
include("utils/vector.jl")
include("utils/make_table.jl")
include("utils/to_df.jl")
include("utils/labels.jl")
# study
include("study/create.jl")
include("study/info.jl")
# recorder
include("recorder/ftt.jl")                                                                                              
# stim
include("stim/ect.jl")                                                                                              
include("stim/tes.jl")                                                                                              
include("stim/tes_model.jl")

end # NeuroAnalyzer