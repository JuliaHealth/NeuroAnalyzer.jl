module NeuroJ

using CSV
using DataFrames
using Distances
using DSP
using FFTW
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
using ScatteredInterpolation
using Simpson
using StatsKit
using StatsPlots

mutable struct EEG
    eeg_header::Dict
    eeg_time::Vector{Float64}
    eeg_signals::Array{Float64, 3}
    eeg_components::Vector{Any}
end

if VERSION < v"1.0.0"
    @warn("This version of NeuroJ requires Julia 1.0 or above.")
end

include("neuroj.jl")
export neuroj_reload_plugins
export neuroj_version
neuroj_reload_plugins()

include("misc.jl")
export cart2pol
export cmax
export cmin
export demean
export fft0
export freqs
export generate_gaussian
export generate_hanning
export generate_morlet
export generate_sinc
export generate_sine
export hildebrand_rule
export hz2rads
export ifft0
export jaccard_similarity
export k_categories
export linspace
export logspace
export matrix_sort
export matrix_sortperm
export nextpow2
export pad0
export pol2cart
export rads2hz
export rms
export rmse
export sph2cart
export tuple_order
export vsearch
export vsplit
export z_score
export zero_pad

include("signal.jl")
export signal_add_noise
export signal_autocov
export signal_average
export signal_band
export signal_band_power
export signal_ci95
export signal_coherence
export signal_cor
export signal_cov
export signal_crosscov
export signal_delete_channel
export signal_demean
export signal_derivative
export signal_detect_epoch_euclid
export signal_detect_epoch_flat
export signal_detect_epoch_p2p
export signal_detect_epoch_rmsd
export signal_detect_epoch_rmse
export signal_detrend
export signal_difference
export signal_downsample
export signal_entropy
export signal_epochs
export signal_epochs_stats
export signal_fconv
export signal_filter
export signal_ica
export signal_ica_reconstruct
export signal_invert_polarity
export signal_make_spectrum
export signal_mean
export signal_mi
export signal_normalize_minmax
export signal_normalize_zscore
export signal_pca
export signal_psd
export signal_reference_car
export signal_reference_channel
export signal_spectrogram
export signal_spectrum
export signal_stationarity
export signal_stationarity_hilbert
export signal_stationarity_mean
export signal_stationarity_var
export signal_taper
export signal_tconv
export signal_total_power
export signal_trim
export signal_upsample

include("eeg_io.jl")
export eeg_export_csv
export eeg_import_ced
export eeg_import_edf
export eeg_import_elc
export eeg_import_locs
export eeg_load
export eeg_load_electrodes
export eeg_load_electrodes!
export eeg_save

include("eeg_edit.jl")
export eeg_reset_components
export eeg_reset_components!
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
export eeg_list_components
export eeg_extract_component
export eeg_delete_component
export eeg_delete_component!
export eeg_add_component
export eeg_add_component!
export eeg_detect_bad_epochs
export eeg_delete_bad_epochs
export eeg_delete_bad_epochs!
export eeg_add_labels
export eeg_add_labels!
export eeg_edit_channel
export eeg_edit_channel!
export eeg_keep_eeg_channels
export eeg_keep_eeg_channels!

include("eeg_process.jl")
export eeg_reference_channel
export eeg_reference_channel!
export eeg_reference_car
export eeg_reference_car!
export eeg_derivative
export eeg_derivative!
export eeg_detrend
export eeg_detrend!
export eeg_taper
export eeg_taper!
export eeg_demean
export eeg_demean!
export eeg_normalize_zscore
export eeg_normalize_zscore!
export eeg_normalize_minmax
export eeg_normalize_minmax!
export eeg_upsample
export eeg_upsample!
export eeg_tconv
export eeg_tconv!
export eeg_filter
export eeg_filter!
export eeg_downsample
export eeg_downsample!
export eeg_pca
export eeg_pca!
export eeg_fconv
export eeg_fconv!
export eeg_ica
export eeg_ica!
export eeg_average
export eeg_average!
export eeg_ica_reconstruct
export eeg_ica_reconstruct!
export eeg_resample
export eeg_resample!
export eeg_invert_polarity
export eeg_invert_polarity!

include("eeg_analyze.jl")
export eeg_total_power
export eeg_total_power!
export eeg_band_power
export eeg_cov
export eeg_cov!
export eeg_cor
export eeg_cor!
export eeg_autocov
export eeg_autocov!
export eeg_crosscov
export eeg_crosscov!
export eeg_crosscov
export eeg_psd
export eeg_psd!
export eeg_stationarity
export eeg_stationarity!
export eeg_mi
export eeg_mi!
export eeg_mi
export eeg_entropy
export eeg_entropy!
export eeg_band
export eeg_coherence
export eeg_coherence
export eeg_freqs
export eeg_freqs!
export eeg_difference
export eeg_pick
export eeg_epochs_stats
export eeg_epochs_stats!
export eeg_spectrogram
export eeg_spectrogram!
export eeg_spectrum
export eeg_spectrum!
export eeg_s2t
export eeg_t2s

include("eeg_plots.jl")
export eeg_draw_head
export eeg_plot
export eeg_plot_avg
export eeg_plot_butterfly
export eeg_plot_covmatrix
export eeg_plot_electrodes
export eeg_plot_filter_response
export eeg_plot_histogram
export eeg_plot_ica
export eeg_plot_matrix
export eeg_plot_psd
export eeg_plot_spectrogram
export eeg_plot_topo
export eeg_plot_band
export eeg_plot_save
export signal_plot
export signal_plot_avg
export signal_plot_band
export signal_plot_butterfly
export signal_plot_histogram
export signal_plot_ica
export signal_plot_psd
export signal_plot_spectrogram

include("nstim.jl")
export tes_dose

end