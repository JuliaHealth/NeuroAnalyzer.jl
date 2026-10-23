"""
Julia toolbox for analyzing neurophysiological data.

https://neuroanalyzer.org
"""

module NeuroAnalyzer

VERSION >= v"1.11.0" ||
    throw(ArgumentError("NeuroAnalyzer requires Julia 1.11.0 or above."))

# ---------------------------------------------------------------------------

# set constants

const VER = v"0.26.11-dev"
# false for the stable branch, true for the devel branch
const allow_wip = occursin("dev", string(VER))
const data_types = [
    "eeg",
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
    "erf",
    "tpt",
]
const channel_types = [
    "all",
    "eeg",
    "ecog",
    "seeg",
    "ieeg",
    "csd",
    "meg",
    "grad",
    "mag",
    "nirs",
    "nirs_int",
    "nirs_od",
    "nirs_dmean",
    "nirs_dvar",
    "nirs_dskew",
    "nirs_mua",
    "nirs_musp",
    "nirs_hbo",
    "nirs_hbr",
    "nirs_hbt",
    "nirs_h2o",
    "nirs_lipid",
    "nirs_bfi",
    "nirs_hrf_dod",
    "nirs_hrf_dmean",
    "nirs_hrf_dvar",
    "nirs_hrf_dskew",
    "nirs_hrf_hbo",
    "nirs_hrf_hbr",
    "nirs_hrf_hbt",
    "nirs_hrf_bfi",
    "nirs_aux",
    "ecg",
    "emg",
    "eog",
    "ref",
    "mrk",
    "sensors",
    "accel",
    "magfld",
    "orient",
    "angvel",
    "mep",
    "eda",
    "other",
]
const channel_units = [
    "μV", "mV", "V", "μV/m²", "fT", "fT/cm", "μM/mm", "m/s²", "μT", "°", "μS", "rad/s",
    "",
]
const fiducial_points = (
    nasion = (0.0, 1.03, -0.2),
    inion = (0.0, -1.03, -0.2),
    lpa = (-1.04, 0.2, -0.2),
    rpa = (1.04, 0.2, -0.2),
)
const PATH = normpath(joinpath(@__DIR__, ".."))

# ---------------------------------------------------------------------------

# preferences, resources and plugins are set in __init__() (at load time, not at precompilation)

progress_bar = true
verbose = true
exclude_bads = false
colors = true
res_path = ""
plugins_path = ""

# ---------------------------------------------------------------------------

# load dependencies

const _NA_REQUIRED_PACKAGES = [
    "Artifacts",
    "Cairo",
    "ColorSchemes",
    "ComplexityMeasures",
    "ContinuousWavelets",
    "Crayons",
    "Dates",
    "Deconvolution",
    "DICOM",
    "Dierckx",
    "DSP",
    "Einsum",
    "FFTW",
    "FileIO",
    "FilePaths",
    "FindPeaks1D",
    "FourierTools",
    "FractalDimensions",
    "GeometryBasics",
    "Git",
    "GLMakie",
    "GR",
    "Gtk4",
    "Hurst",
    "ImageBinarization",
    "ImageFiltering",
    "ImageMorphology",
    "Images",
    "InformationMeasures",
    "Interpolations",
    "Jacobi",
    "JLD2",
    "JSON",
    "KernelDensity",
    "LibSerialPort",
    "LinearAlgebra",
    "LinRegOutliers",
    "Logging",
    "MAT",
    "MLJ",
    "NaturalNeighbours",
    "NPZ",
    "Optim",
    "PiGPIO",
    "Pkg",
    "Plots",
    "Polynomials",
    "Preferences",
    "PrettyTables",
    "ProgressMeter",
    "Random",
    "REPL",
    "SavitzkyGolay",
    "ScatteredInterpolation",
    "StatsFuns",
    "StatsKit",
    "StatsModels",
    "StatsPlots",
    "TimeZones",
    "TOML",
    "Trapz",
    "WAV",
    "Wavelets",
    "WaveletsExt",
    "XDF",
]
for pkg in _NA_REQUIRED_PACKAGES
    @eval using $(Symbol(pkg))
end

# functions of Base and of dependencies that NeuroAnalyzer extends with methods
# (for NEURO objects, etc.) instead of defining functions of the same name, so
# that they do not clash when used together; the non-Base ones are re-exported
# (names exported by two packages extend one of them: plot → Makie, not Plots;
# entropy → StatsBase, not Images; transform → DataFrames, not MultivariateStats;
# trim! → StatsBase, not GridLayoutBase)
import Base: chop, cis, filter, filter!, join, summary
import DataFrames: describe, transform
import DSP: coherence, resample, rms, spectrogram
import GLMakie: plot
import Images: reflect, rmse
import LinearAlgebra: normalize, reflect!
import StatsKit.StatsBase: entropy, sem, standardize, trim, trim!, zscore
import Wavelets: downsample, upsample

# ---------------------------------------------------------------------------

# Gtk4 modifiers

const mask_shift = Gtk4.ModifierType_SHIFT_MASK
const mask_ctrl = Gtk4.ModifierType_CONTROL_MASK
const mask_alt = Gtk4.ModifierType_ALT_MASK

# ---------------------------------------------------------------------------
# internal helpers
# ---------------------------------------------------------------------------

"""Planar interpolation methods (`imethod`): ScatteredInterpolation.jl, NaturalNeighbours.jl and spherical spline."""
const _IMETHODS = [:sh, :mq, :imq, :tp, :nn, :ga, :sib, :sib1, :lap, :far, :hiy, :tri, :sph]

""" Columns that hold numeric coordinate data in a locations DataFrame."""
const _LOCS_COORD_COLS = [
    "loc_x", "loc_y", "loc_z",
    "loc_radius", "loc_theta",
    "loc_radius_sph", "loc_theta_sph", "loc_phi_sph",
]

"""Keyword arguments applied to every locked/non-interactive Axis."""
const _AXIS_LOCK_KWARGS = (
    xzoomlock = true,
    yzoomlock = true,
    xpanlock  = true,
    ypanlock  = true,
    xrectzoom = false,
    yrectzoom = false,
)

"""Keyword arguments applied to every time x-axis."""
const _TIME_AXIS_KWARGS = (
    xticks             = LinearTicks(10),
    xminorticksvisible = true,
    xminorticks        = IntervalsBetween(10),
)

# logarithmic frequency axis: major ticks at 1-2-5 steps, minor ticks at 1-9 × 10^k
const _LOG_FRQ_TICKS = [m * 10.0^k for k in -2:3 for m in (1, 2, 5)]
const _LOG_FRQ_MINORTICKS = [m * 10.0^k for k in -2:3 for m in 1:9]

"""Keyword arguments applied to every frequency axis (`axis`: `:x` or `:y`; `frq`: `:lin` or `:log`)."""
function _frq_axis_kwargs(frq::Symbol; axis::Symbol = :x)::NamedTuple
    if frq === :lin
        ticks = LinearTicks(15)
        minorticks = IntervalsBetween(10)
        scale = identity
    else
        ticks = (_LOG_FRQ_TICKS, [isinteger(v) ? string(Int(v)) : string(v) for v in _LOG_FRQ_TICKS])
        minorticks = _LOG_FRQ_MINORTICKS
        scale = log10
    end
    names = Symbol.(string(axis) .* ["ticks", "minorticksvisible", "minorticks", "scale"])
    return NamedTuple{Tuple(names)}((ticks, true, minorticks, scale))
end

"""Apply standard font sizes to an Axis."""
function _style_axis!(ax)
    ax.titlesize      = 18
    ax.xlabelsize     = 18
    ax.ylabelsize     = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12
    return ax
end

# ---------------------------------------------------------------------------

# define structures

mutable struct HEADER
    subject::Dict
    recording::Dict
    experiment::Dict
end

mutable struct NEURO
    header::NeuroAnalyzer.HEADER
    history::Vector{String}
    markers::DataFrame
    locs::DataFrame
    time_pts::Vector{Float64}
    epoch_time::Vector{Float64}
    data::Array{Float64, 3}
end

# equivalent current dipole, as fitted by dipole_fit()
mutable struct DIPOLE
    pos::Matrix{Float64}            # positions: 1 × 3 (fixed dipole) or time points × 3 (moving dipole)
    q::Matrix{Float64}              # moments: 3 × time points
    gof::Vector{Float64}            # goodness of fit (0–1), one value per fit
    sphere_center::Vector{Float64}  # center of the fitted head sphere
    sphere_radius::Float64          # radius of the fitted head sphere
end

# single dipole with position `pos` and moment `q` in a unit sphere (no fit)
DIPOLE(pos::NTuple{3, Real}, q::NTuple{3, Real}) = DIPOLE(
    reshape(Float64[pos...], 1, 3), reshape(Float64[q...], 3, 1), Float64[], zeros(3), 1.0,
)

# ---------------------------------------------------------------------------

# load NA functions

include("na/internal.jl")
include("na/setup.jl")
include("na/plugins.jl")

# load sub-modules

# dependencies
include("deps/PhaseSlopeIndex.jl/PhaseSlopeIndex.jl")
include("deps/FIRLSFilterDesign.jl/FIRLSFilterDesign.jl")
include("deps/Simpson.jl/Simpson.jl")

# internal
include("internal/channels.jl")
include("internal/check.jl")
include("internal/convolution.jl")
include("internal/create_header.jl")
include("internal/edf.jl")
include("internal/epochs.jl")
include("internal/fiff.jl")
include("internal/fir_response.jl")
include("internal/gdf_etp.jl")
include("internal/gradient.jl")
include("internal/gui.jl")
include("internal/interpolate.jl")
include("internal/labeled_matrix.jl")
include("internal/labels.jl")
include("internal/len.jl")
include("internal/locs.jl")
include("internal/map_channels.jl")
include("internal/markers.jl")
include("internal/mesh.jl")
include("internal/misc.jl")
include("internal/ml.jl")
include("internal/tpt.jl")
include("internal/plots.jl")
include("internal/recorder.jl")
include("internal/reflect_chop.jl")
include("internal/select.jl")
include("internal/cm.jl")
include("internal/tester.jl")
include("internal/time.jl")
include("internal/wl2ext.jl")
include("internal/vec.jl")

# utils
include("utils/apply.jl")
include("utils/array.jl")
include("utils/channel.jl")
include("utils/fft.jl")
include("utils/filter.jl")
include("utils/findpeaks.jl")
include("utils/frequency.jl")
include("utils/generate.jl")
include("utils/gradient.jl")
include("utils/info.jl")
include("utils/labels.jl")
include("utils/matrix.jl")
include("utils/misc.jl")
include("utils/note.jl")
include("utils/pad.jl")
include("utils/phase.jl")
include("utils/play.jl")
include("utils/time.jl")
include("utils/to_df.jl")
include("utils/vector.jl")
include("utils/view_header.jl")
include("utils/wavelets.jl")
include("utils/fwhm.jl")
include("utils/mri_coordinates.jl")

# io
include("io/export_csv.jl")
include("io/export_locs.jl")
include("io/export_markers.jl")
include("io/import_alice4.jl")
include("io/import_bdf.jl")
include("io/import_bv.jl")
include("io/import_cnt.jl")
include("io/import_csv.jl")
include("io/import_dat.jl")
include("io/import_digitrack.jl")
include("io/import_duomag.jl")
include("io/import_edf_annotations.jl")
include("io/import_edf.jl")
include("io/import_fiff.jl")
include("io/import_ft.jl")
include("io/import_gdf.jl")
include("io/import.jl")
include("io/import_locs.jl")
include("io/import_montage.jl")
include("io/import_ncs.jl")
include("io/import_nirs.jl")
include("io/import_nirx.jl")
include("io/import_npy.jl")
include("io/import_nwb.jl")
include("io/import_set.jl")
include("io/import_snirf.jl")
include("io/import_thymatron.jl")
include("io/import_xdf.jl")
include("io/load_locs.jl")
include("io/save_load.jl")

# locs
include("locs/add_locs.jl")
include("locs/center.jl")
include("locs/convert.jl")
include("locs/details.jl")
include("locs/edit.jl")
include("locs/flip.jl")
include("locs/generate.jl")
include("locs/origin.jl")
include("locs/rotate.jl")
include("locs/scale.jl")
include("locs/swap.jl")

# edit
include("edit/channel.jl")
include("edit/create.jl")
include("edit/delete_channel.jl")
include("edit/delete_epoch.jl")
include("edit/delete_optode.jl")
include("edit/reject.jl")
include("edit/epoch.jl")
include("edit/extract.jl")
include("edit/join.jl")
include("edit/lrinterpolate.jl")
include("edit/marker.jl")
include("edit/mlinterpolate.jl")
include("edit/plinterpolate.jl")
include("edit/reflect_chop.jl")
include("edit/remove_pops.jl")
include("edit/sphinterpolate.jl")
include("edit/trim.jl")
include("edit/vch.jl")

# process
include("process/add_signal.jl")
include("process/artrem_cwd.jl")
include("process/bpsplit.jl")
include("process/cbp.jl")
include("process/chpi.jl")
include("process/csd.jl")
include("process/cwd.jl")
include("process/denoise_fft.jl")
include("process/denoise_wavelet.jl")
include("process/denoise_wien.jl")
include("process/derivative.jl")
include("process/detrend.jl")
include("process/dwd.jl")
include("process/emd.jl")
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
include("process/maxwell.jl")
include("process/normalize.jl")
include("process/normpower.jl")
include("process/npl.jl")
include("process/od2conc.jl")
include("process/pca.jl")
include("process/reference.jl")
include("process/remove_dc.jl")
include("process/remove_powerline.jl")
include("process/resample.jl")
include("process/ssp.jl")
include("process/standardize.jl")
include("process/taper.jl")
include("process/tconv.jl")
include("process/wbp.jl")

# stats
include("stats/bootstrap.jl")
include("stats/ci.jl")
include("stats/correlation.jl")
include("stats/crit.jl")
include("stats/descriptive.jl")
include("stats/linreg.jl")
include("stats/make_table.jl")
include("stats/misc.jl")
include("stats/outliers.jl")
include("stats/p.jl")
include("stats/pca.jl")
include("stats/se.jl")
include("stats/summary.jl")
include("stats/zscore.jl")

# analyze
include("analyze/acor.jl")
include("analyze/acov.jl")
include("analyze/aecor.jl")
include("analyze/amp.jl")
include("analyze/ampdiff.jl")
include("analyze/axc2frq.jl")
include("analyze/asy_idx.jl")
include("analyze/band_mpower.jl")
include("analyze/band_power.jl")
include("analyze/coherence.jl")
include("analyze/corm.jl")
include("analyze/covm.jl")
include("analyze/cosim.jl")
include("analyze/cph.jl")
include("analyze/corr.jl")
include("analyze/cpsd.jl")
include("analyze/dirinrg.jl")
include("analyze/dissimilarity.jl")
include("analyze/entropy.jl")
include("analyze/envelopes.jl")
include("analyze/erop.jl")
include("analyze/eros.jl")
include("analyze/erp.jl")
include("analyze/fooof.jl")
include("analyze/frqinst.jl")
include("analyze/ftt.jl")
include("analyze/ged.jl")
include("analyze/gfp.jl")
include("analyze/granger.jl")
include("analyze/ghexp.jl")
include("analyze/hfd.jl")
include("analyze/hjorth.jl")
include("analyze/hmspectrum.jl")
include("analyze/hrv.jl")
include("analyze/ispc.jl")
include("analyze/itpc.jl")
include("analyze/lat_idx.jl")
include("analyze/mdiff.jl")
include("analyze/mep_peaks.jl")
include("analyze/mi.jl")
include("analyze/msci95.jl")
include("analyze/pacor.jl")
include("analyze/peak.jl")
include("analyze/phdiff.jl")
include("analyze/phs.jl")
include("analyze/pli.jl")
include("analyze/dpli.jl")
include("analyze/wpli.jl")
include("analyze/plv.jl")
include("analyze/iplv.jl")
include("analyze/psa.jl")
include("analyze/psi.jl")
include("analyze/psd.jl")
include("analyze/psd_rel.jl")
include("analyze/psd_slope.jl")
include("analyze/rms.jl")
include("analyze/sef.jl")
include("analyze/segments.jl")
include("analyze/snr.jl")
include("analyze/specseg.jl")
include("analyze/spectrogram.jl")
include("analyze/stationarity.jl")
include("analyze/stats.jl")
include("analyze/std.jl")
include("analyze/sumsim.jl")
include("analyze/sym_idx.jl")
include("analyze/tkeo.jl")
include("analyze/total_power.jl")
include("analyze/tpt.jl")
include("analyze/transform.jl")
include("analyze/vartest.jl")
include("analyze/xcor.jl")
include("analyze/xcov.jl")
include("analyze/zip.jl")

# plots
include("plots/cairo.jl")
include("plots/misc.jl")
include("plots/plot_filter.jl")
include("plots/plot_locs.jl")
include("plots/plot_locs3d.jl")
include("plots/plot_locs_nirs.jl")
include("plots/plot_gridlocs.jl")
include("plots/plot_save.jl")
include("plots/plot_topo.jl")
include("plots/plot_connectivity_circle.jl")
include("plots/plot_coherence.jl")
include("plots/plot_efield.jl")
include("plots/plot_stats.jl")
include("plots/plot_matrix.jl")
include("plots/plot_phase.jl")
include("plots/plot_tf.jl")
include("plots/plot_dipole2d.jl")
include("plots/plot_dipole3d.jl")
include("plots/plot_erp.jl")
include("plots/plot_mep.jl")
include("plots/plot_phs.jl")
include("plots/plot_psd.jl")
include("plots/plot_spectrogram.jl")
include("plots/plot.jl")
include("plots/plot_signal_cont1.jl")
include("plots/plot_signal_cont2.jl")
include("plots/plot_signal_ep.jl")

# gui
# include("gui/iedit.jl")
# include("gui/ispectrogram.jl")
# include("gui/itopo.jl")
# include("gui/iview_ica.jl")

# recorder
include("recorder/edar.jl")
include("recorder/rt_plotter.jl")

# tester
include("tester/ftt.jl")
include("tester/iavh.jl")
include("tester/tpt.jl")
include("tester/test_run.jl")

# stim
include("stim/ect.jl")
include("stim/tes.jl")

# modeling
include("model/dipole_fit.jl")
include("model/efield.jl")
include("model/tes_model.jl")

# ---------------------------------------------------------------------------

# code in __init__() runs every time the package is loaded
# top-level code runs only once, when the package is precompiled
function __init__()
    # load preferences
    global progress_bar = @load_preference("progress_bar", true)
    global verbose = @load_preference("verbose", true)
    global exclude_bads = @load_preference("exclude_bads", false)
    global colors = @load_preference("colors", true)

    # show major parameters
    @info "NeuroAnalyzer v$(NeuroAnalyzer.VER)"
    @info "NeuroAnalyzer path: $(NeuroAnalyzer.PATH)"
    @info "Preferences:"
    @info " Progress bar: $progress_bar"
    @info "      Verbose: $verbose"
    @info " Exclude bads: $exclude_bads"
    @info "       Colors: $colors"
    @info "System info:"
    @info "  Free memory: $(round(Sys.free_memory() / 2^20, digits = 1)) MB"

    # set package options
    # FFTW: one thread - most transforms are small (Welch / STFT segments, wavelets) and many run
    # inside Julia-threaded loops; multi-threaded FFTW makes them ~100× slower (thread overhead and
    # oversubscription), and large single transforms gain nothing
    FFTW.set_num_threads(1)
    BLAS.set_num_threads(max(1, Sys.CPU_THREADS - 1))

    # setup resources
    @info "Preparing resources"
    global res_path = joinpath(artifact"NeuroAnalyzer_resources", "neuroanalyzer-resources")

    # load plugins (not while another package that depends on NeuroAnalyzer is
    # being precompiled: plugins are included into this module)
    global plugins_path = joinpath(homedir(), "NeuroAnalyzer", "plugins")
    if ccall(:jl_generating_output, Cint, ()) == 0
        @info "Loading plugins"
        try
            isdir(plugins_path) || mkpath(plugins_path)
            isempty(readdir(plugins_path)) || plugins_reload()
        catch err
            @warn "Plugins could not be loaded: $err"
        end
    end

    return nothing
end

end # NeuroAnalyzer
