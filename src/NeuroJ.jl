module NeuroJ

import FFTW
import DSP
import Statistics
import LinearAlgebra
import Simpson

include("eeg.jl")
export signal_derivative
export signal_band_power
export signal_make_spectrum
export signal_detrend
export signals_ci95
export signals_mean
export signals_difference
export signal_autocov

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
export nexpow2
export vsplit
export rms
export db

include("nstim.jl")
export tes_dose

end
