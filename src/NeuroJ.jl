module NeuroJ

using DataFrames
using DSP
using FFTW
using LinearAlgebra
using Plots
using StatsKit
using Simpson

struct EEG
    eeg_file_header::Dict
    eeg_signal_header::Dict
    eeg_signals::Matrix
end

include("eeg.jl")
export eeg_drop_channel

include("eeg_load_edf.jl")
export eeg_load_edf

include("signal.jl")
export signal_autocov
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
export normalize_mean
export normalize_minmax
export pad0
export pol2cart
export rads2hz
export rms
export sine
export vsearch
export vsplit
export z_score
export zero_pad

include("nstim.jl")
export tes_dose

end
