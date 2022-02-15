module NeuroJ

using CSV
using DataFrames
using Distances
using DSP
using FFTW
using InformationMeasures
using Interpolations
using JLD2
using LinearAlgebra
using Pkg
using Plots
using Plots.PlotMeasures
using Polynomials
using Simpson
using StatsKit

mutable struct EEG
    eeg_header::Dict
    eeg_time::Vector{Float64}
    eeg_signals::Array{Float64, 3}
end

if VERSION < v"1.0.0"
    @warn("This version of EDFPlus requires Julia 1.0 or above.")
end

function neuroj_version()
    m = Pkg.Operations.Context().env.manifest
    println("NeuroJ version: $(m[findfirst(v->v.name=="NeuroJ", m)].version)")
    println("Imported packages:")
    required_packages = ["CSV",
                         "DataFrames",
                         "Distances",
                         "DSP",
                         "FFTW",
                         "InformationMeasures",
                         "Interpolations",
                         "JLD2",
                         "LinearAlgebra",
                         "Plots",
                         "Polynomials",
                         "Simpson",
                         "StatsKit"]
    for idx in 1:length(required_packages)
        pkg = lpad(required_packages[idx], 20 - length(idx), " ")
        pkg_ver = m[findfirst(v->v.name==required_packages[idx], m)].version
        println("$pkg $pkg_ver")
    end
end

export neuroj_version

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
export k
export linspace
export logspace
export matrix_sort
export matrix_sortperm
export nextpow2
export pad0
export pol2cart
export rads2hz
export rms
export vsearch
export vsplit
export z_score
export zero_pad

include("signal.jl")
export signal_add_noise
export signal_autocov
export signal_average
export signal_band_power
export signal_ci95
export signal_coherence
export signal_cor
export signal_cov
export signal_crosscov
export signal_delete_channel
export signal_demean
export signal_derivative
export signal_detrend
export signal_difference
export signal_downsample
export signal_entropy
export signal_epochs
export signal_fconv
export signal_filter
export signal_make_spectrum
export signal_mean
export signal_mi
export signal_normalize_minmax
export signal_normalize_zscore
export signal_pca
export signal_psd
export signal_reference_car
export signal_reference_channel
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

include("eeg.jl")
export eeg_autocov
export eeg_band
export eeg_band_power
export eeg_coherence
export eeg_cor
export eeg_cov
export eeg_crosscov
export eeg_delete_channel
export eeg_demean
export eeg_derivative
export eeg_detrend
export eeg_difference
export eeg_downsample
export eeg_entropy
export eeg_epochs
export eeg_extract_channel
export eeg_extract_epoch
export eeg_fconv
export eeg_filter
export eeg_freqs
export eeg_get_channel
export eeg_history
export eeg_info
export eeg_keep_channel
export eeg_labels
export eeg_mi
export eeg_normalize_minmax
export eeg_normalize_zscore
export eeg_pca
export eeg_pca
export eeg_psd
export eeg_reference_car
export eeg_reference_channel
export eeg_rename_channel
export eeg_sr
export eeg_stationarity
export eeg_taper
export eeg_tconv
export eeg_total_power
export eeg_trim
export eeg_upsample
export eeg_edit
export eeg_show_header
export eeg_delete_epoch
export eeg_keep_epoch
export eeg_pick

include("eeg_io.jl")
export eeg_import_edf
export eeg_load
export eeg_load_electrode_positions
export eeg_save
export eeg_export_csv

include("plots.jl")
export eeg_draw_head
export eeg_plot
export eeg_plot_avg
export eeg_plot_butterfly
export eeg_plot_covmatrix
export eeg_plot_electrodes
export eeg_plot_matrix
export eeg_plot_psd
export eeg_plot_spectrogram
export filter_response
export signal_plot
export signal_plot_avg
export signal_plot_butterfly
export signal_plot_psd
export signal_plot_spectrogram

include("nstim.jl")
export tes_dose

end