module NeuroJ

using DataFrames
using DSP
using FFTW
using JLD2
using LinearAlgebra
using Plots
using Simpson
using StatsKit

struct EEG
    eeg_object_header::Dict
    eeg_signal_header::Dict
    eeg_signals::Matrix
end

include("eeg.jl")
export eeg_derivative
export eeg_detrend
export eeg_drop_channel
export eeg_filter_butter
export eeg_plot
export eeg_total_power
export eeg_band_power
export eeg_make_spectrum
export eeg_reference_channel
export eeg_reference_car
export eeg_load
export eeg_save
export eeg_get_channel_idx
export eeg_get_channel_name
export eeg_rename_channel
export eeg_taper
export eeg_demean
export eeg_normalize_mean
export eeg_normalize_minmax
export eeg_get_channel
export eeg_cov
export eeg_cor

include("eeg_load_edf.jl")
export eeg_load_edf

include("signal.jl")
export signal_autocov
export signal_total_power
export signal_band_power
export signal_ci95
export signal_crosscov
export signal_derivative
export signal_detrend
export signal_difference
export signal_drop_channel
export signal_epoch
export signal_filter_butter
export signal_make_spectrum
export signal_mean
export signal_plot
export signal_spectrum
export signal_reference_channel
export signal_reference_car
export signal_taper
export signal_demean
export signal_normalize_mean
export signal_normalize_minmax
export signal_cov
export signal_cor
export signal_add_noise

include("misc.jl")
export cart2pol
export cvangle
export db
export demean
export fft0
export frequencies
export hann
export hildebrand_rule
export hz2rads
export ifft0
export jaccard_similarity
export k
export linspace
export logspace
export matrix_sort
export matrix_sortperm
export nexpow2
export pad0
export pol2cart
export rads2hz
export rms
export sine
export vsearch
export vsplit
export z_score
export zero_pad
export cmin
export cmax
export sinc

include("nstim.jl")
export tes_dose

end
