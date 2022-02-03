module NeuroJ

using DataFrames
using DSP
using FFTW
using Interpolations
using JLD2
using LinearAlgebra
using Plots
using Simpson
using StatsKit

struct EEG
    eeg_header::Dict
    eeg_time::Vector{Float64}
    eeg_signals::Array{Float64, 3}
end

include("eeg.jl")
export eeg_derivative
export eeg_detrend
export eeg_delete_channel
export eeg_keep_channel
export eeg_filter
export eeg_total_power
export eeg_band_power
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
export eeg_upsample
export eeg_history
export eeg_info
export eeg_epochs
export eeg_get_epoch
export eeg_labels
export eeg_samplingrate
export eeg_tconv
export eeg_downsample
export eeg_autocov
export eeg_crosscov
export eeg_psd

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
export signal_epochs
export signal_filter
export signal_make_spectrum
export signal_mean
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
export signal_upsample
export signal_tconv
export signal_downsample
export signal_psd

include("misc.jl")
export generate_time
export cart2pol
export cvangle
export db
export demean
export fft0
export freqs
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
export morlet

include("nstim.jl")
export tes_dose

include("plots.jl")
export signal_plot
export eeg_plot
export eeg_draw_head
export filter_response
export signal_plot_avg
export eeg_plot_avg
export signal_plot_butterfly
export eeg_plot_butterfly
export signal_plot_psd
export eeg_plot_psd

end