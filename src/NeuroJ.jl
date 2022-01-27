module NeuroJ

using DSP
using FFTW
using LinearAlgebra
using StatsKit
using Simpson

struct EEG
    eeg_file_header::Array
    eeg_signal_header::DataFrame
    eeg_signals::Matrix
end

include("eeg.jl")
export signal_derivative
export signal_band_power
export signal_make_spectrum
export signal_detrend
export signals_ci95
export signals_mean
export signals_difference
export signal_autocov
export signals_crosscov
export signal_spectrum
export eeg_load

include("misc.jl")
export linspace
export logspace
export zero_pad
export vsearch
export cart2pol
export pol2cart
export minmax_scaler
export cvangle
export hann
export hildebrand_rule
export jaccard_similarity
export fft0
export ifft0
export nexpow2
export vsplit
export rms
export db
export sine
export frequencies
export matrix_sortperm
export matrix_sort
export pad0
export hz2rads

include("nstim.jl")
export tes_dose

end
