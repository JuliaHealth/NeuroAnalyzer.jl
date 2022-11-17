__precompile__()

module NeuroAnalyzer

const na_ver = v"0.22.11"

if VERSION < v"1.7.0"
    @error("This version of NeuroAnalyzer requires Julia 1.7.0 or above.")
end

using ColorSchemes
using CSV
using CubicSplines
using CUDA
using DataFrames
using Deconvolution
using Distances
using DSP
using FFTW
using FileIO
using FindPeaks1D
using Git
using GLM
using GLMakie
using HypothesisTests
using InformationMeasures
using Interpolations
using JLD2
using LinearAlgebra
using Loess
using MultivariateStats
using Pkg
using Plots
using Plots.PlotMeasures
using Polynomials
using Preferences
using ProgressMeter
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

mutable struct EEG
    eeg_header::Dict
    eeg_time::Vector{Float64}
    eeg_epochs_time::Vector{Float64}
    eeg_signals::Array{Float64, 3}
    eeg_components::Vector{Any}
    eeg_markers::DataFrame
    eeg_locs::DataFrame
end

mutable struct STUDY
    study_header::Dict{Symbol, Any}
    study_eeg::Vector{NeuroAnalyzer.EEG}
    study_group::Vector{Symbol}
end

FFTW.set_num_threads(Sys.CPU_THREADS)
BLAS.set_num_threads(Sys.CPU_THREADS)

include("na.jl")
export na_info
export na_plugins_reload
export na_plugins_list
export na_plugins_remove
export na_plugins_install
export na_plugins_update
export na_set_use_cuda
export na_set_progress_bar
export na_set_plugins_path
export na_set_prefs
export na_set_verbose
export na_version

# preferences
if Sys.isunix() || Sys.isapple()
    def_plugins_path = "$(homedir())/NeuroAnalyzer/plugins/"
elseif Sys.iswindows()
    def_plugins_path = "$(homedir())\\NeuroAnalyzer\\plugins\\"
end
const use_cuda = @load_preference("use_cuda", false)
const progress_bar = @load_preference("progress_bar", true)
const plugins_path = @load_preference("plugins_path", def_plugins_path)
const verbose = @load_preference("verbose", true)
isdir(plugins_path) || mkdir(plugins_path)
na_set_prefs(use_cuda=use_cuda, plugins_path=plugins_path, progress_bar=progress_bar, verbose=verbose)

# reload plugins
na_plugins_reload()

include("internal.jl")

include("low_level.jl")
export linspace
export logspace
export m_pad0
export cmax
export cmin
export vsearch
export cart2pol
export pol2cart
export sph2cart
export cart2sph
export generate_window
export fft0
export ifft0
export nextpow2
export vsplit
export s_rms
export generate_sine
export s_freqs
export m_sortperm
export m_sort
export pad0
export pad2
export hz2rads
export rads2hz
export generate_sinc
export generate_morlet
export generate_gaussian
export tuple_order
export s2_rmse
export m_norm
export s_cov
export s2_cov
export s_dft
export s_msci95
export s2_mean
export s2_difference
export s_acov
export s2_xcov
export s_spectrum
export s_total_power
export s_band_power
export s_taper
export s_detrend
export s_demean
export s_normalize_zscore
export s_normalize_minmax
export s_normalize_max
export s_normalize_log
export s_add_noise
export s_resample
export s_invert_polarity
export s_derivative
export s_tconv
export s_filter
export s_psd
export s_stationarity_hilbert
export s_stationarity_mean
export s_stationarity_var
export s_trim
export s2_mi
export s_entropy
export s_negentropy
export s_average
export s2_average
export s2_tcoherence
export s_pca
export s_pca_reconstruct
export s_fconv
export s_ica
export s_ica_reconstruct
export s_spectrogram
export s_detect_epoch_flat
export s_detect_epoch_rmse
export s_detect_epoch_rmsd
export s_detect_epoch_euclid
export s_detect_epoch_p2p
export s_snr
export s_findpeaks
export s_wdenoise
export s2_ispc
export s_itpc
export s2_pli
export s2_ged
export s_frqinst
export s_hspectrum
export t2f
export f2t
export s_wspectrogram
export s_fftdenoise
export s_gfilter
export s_ghspectrogram
export s_tkeo
export s_wspectrum
export a2_cmp
export s_fcoherence
export s2_fcoherence
export a2_l1
export a2_l2
export s_cums
export s_gfp
export s_gfp_norm
export s2_diss
export generate_morlet_fwhm
export f_nearest
export s_band_mpower
export s_rel_psd
export s_wbp
export s_normalize_gauss
export s_cbp
export s_specseg
export s_denoise_wien
export s2_cps
export s2_phdiff
export s_normalize_log10
export s_normalize_neglog
export s_normalize_neglog10
export s_normalize_neg
export s_normalize_pos
export s_normalize_perc
export s_normalize
export s_phases
export s_cwtspectrogram
export s_dwt
export s_idwt
export s_normalize_invroot
export s_cwt
export s_icwt
export t2s
export s2t

include("statistics.jl")
export hildebrand_rule
export jaccard_similarity
export z_score
export k_categories
export effsize
export infcrit
export grubbs
export outlier_detect
export seg_mean
export seg2_mean
export binom_prob
export binom_stat
export cvar_mean
export cvar_median
export cvar
export meang
export meanh
export meanw
export moe
export rng
export se
export pred_int
export sem_diff
export prank
export linreg
export s2_cmp
export s2_cor

include("eeg_io.jl")
export eeg_export_csv
export eeg_import
export eeg_import_edf
export eeg_import_bdf
export eeg_import_ced
export eeg_import_locs
export eeg_import_elc
export eeg_import_tsv
export eeg_import_sfp
export eeg_import_csd
export eeg_load
export eeg_load_electrodes
export eeg_load_electrodes!
export eeg_save
export eeg_save_electrodes
export eeg_add_electrodes
export eeg_add_electrodes!
export eeg_import_digitrack
export eeg_import_bv

include("eeg_edit.jl")
export eeg_copy
export eeg_add_component
export eeg_add_component!
export eeg_list_components
export eeg_extract_component
export eeg_delete_component
export eeg_delete_component!
export eeg_reset_components
export eeg_reset_components!
export eeg_component_idx
export eeg_component_type
export eeg_rename_component
export eeg_rename_component!
export eeg_delete_channel
export eeg_delete_channel!
export eeg_keep_channel
export eeg_keep_channel!
export eeg_get_channel
export eeg_rename_channel
export eeg_rename_channel!
export eeg_extract_channel
export eeg_history
export eeg_labels
export eeg_sr
export eeg_channel_n
export eeg_epoch_n
export eeg_signal_len
export eeg_epoch_len
export eeg_info
export eeg_epochs
export eeg_epochs!
export eeg_extract_epoch
export eeg_trim
export eeg_trim!
export eeg_edit_header
export eeg_edit_header!
export eeg_show_header
export eeg_delete_epoch
export eeg_delete_epoch!
export eeg_keep_epoch
export eeg_keep_epoch!
export eeg_detect_bad_epochs
export eeg_add_labels
export eeg_add_labels!
export eeg_edit_channel
export eeg_edit_channel!
export eeg_keep_channel_type
export eeg_keep_channel_type!
export eeg_view_note
export eeg_epochs_time
export eeg_epochs_time!
export eeg_add_note
export eeg_add_note!
export eeg_delete_note
export eeg_delete_note!
export eeg_replace_channel
export eeg_replace_channel!
export eeg_interpolate_channel
export eeg_interpolate_channel!
export eeg_loc_flipy
export eeg_loc_flipy!
export eeg_loc_flipx
export eeg_loc_flipx!
export eeg_loc_flipz
export eeg_loc_flipz!
export eeg_channel_type
export eeg_channel_type!
export eeg_edit_electrode
export eeg_edit_electrode!
export eeg_electrode_loc
export eeg_loc_swapxy
export eeg_loc_swapxy!
export eeg_loc_sph2cart
export eeg_loc_sph2cart!
export eeg_loc_cart2sph
export eeg_loc_cart2sph!
export eeg_view_markers
export eeg_delete_marker
export eeg_delete_marker!
export eeg_add_marker
export eeg_add_marker!
export eeg_channel_idx
export eeg_vch
export eeg_edit_marker
export eeg_edit_marker!

include("eeg_process.jl")
export eeg_reference_ch
export eeg_reference_ch!
export eeg_reference_car
export eeg_reference_car!
export eeg_reference_a
export eeg_reference_a!
export eeg_reference_m
export eeg_reference_m!
export eeg_derivative
export eeg_derivative!
export eeg_detrend
export eeg_detrend!
export eeg_taper
export eeg_taper!
export eeg_demean
export eeg_demean!
export eeg_normalize
export eeg_normalize!
export eeg_add_noise
export eeg_add_noise!
export eeg_filter
export eeg_filter!
export eeg_pca
export eeg_pca_reconstruct
export eeg_pca_reconstruct!
export eeg_ica
export eeg_ica_reconstruct
export eeg_ica_reconstruct!
export eeg_average
export eeg_average!
export eeg_average
export eeg_invert_polarity
export eeg_invert_polarity!
export eeg_resample
export eeg_resample!
export eeg_upsample
export eeg_upsample!
export eeg_downsample
export eeg_downsample!
export eeg_wdenoise
export eeg_wdenoise!
export eeg_fftdenoise
export eeg_fftdenoise!
export eeg_reference_plap
export eeg_reference_plap!
export eeg_zero
export eeg_zero!
export eeg_wbp
export eeg_wbp!
export eeg_cbp
export eeg_cbp!
export eeg_denoise_wien
export eeg_denoise_wien!
export eeg_scale
export eeg_scale!

include("eeg_analyze.jl")
export eeg_total_power
export eeg_band_power
export eeg_cov
export eeg_cor
export eeg_xcov
export eeg_psd
export eeg_stationarity
export eeg_mi
export eeg_mi
export eeg_entropy
export eeg_negentropy
export eeg_band
export eeg_tcoherence
export eeg_freqs
export eeg_difference
export eeg_pick
export eeg_epochs_stats
export eeg_spectrogram
export eeg_spectrum
export eeg_s2t
export eeg_t2s
export eeg_channels_stats
export eeg_snr
export eeg_standardize
export eeg_standardize!
export eeg_fconv
export eeg_tconv
export eeg_dft
export eeg_mean
export eeg_difference
export eeg_acov
export eeg_tenv
export eeg_tenv_mean
export eeg_tenv_median
export eeg_penv
export eeg_penv_mean
export eeg_penv_median
export eeg_senv
export eeg_senv_mean
export eeg_senv_median
export eeg_ispc
export eeg_itpc
export eeg_pli
export eeg_ispc_m
export eeg_aec
export eeg_ged
export eeg_frqinst
export eeg_itpc_s
export eeg_wspectrogram
export eeg_tkeo
export eeg_wspectrum
export eeg_fcoherence
export eeg_vartest
export eeg_band_mpower
export eeg_rel_psd
export eeg_fbsplit
export eeg_chdiff
export eeg_cps
export eeg_phdiff
export eeg_ampdiff
export eeg_dwt
export eeg_cwt
export eeg_psdslope
export eeg_henv
export eeg_henv_mean
export eeg_henv_median

include("eeg_study.jl")
export eeg_study_create

include("eeg_plots.jl")
export eeg_plot_save
export plot_signal
export plot_signal_avg
export plot_signal_butterfly
export eeg_plot
export plot_psd
export plot_psd_avg
export plot_psd_butterfly
export plot_psd_3d
export plot_psd_topo
export eeg_plot_psd
export plot_spectrogram
export eeg_plot_spectrogram
export plot_electrodes
export plot_electrodes3d
export eeg_plot_electrodes
export plot_matrix
export plot_covmatrix
export plot_filter_response
export plot_weights
export eeg_plot_weights
export plot_histogram
export plot_connections
export eeg_plot_connections
export eeg_plot_stats
export plot_topo
export eeg_plot_topo
export eeg_plot_compose
export plot_empty

include("nstim.jl")
export tes_dose

end # NeuroAnalyzer