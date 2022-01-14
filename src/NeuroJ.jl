module NeuroJ

import FFTW
import DSP
import Statistics

include("eeg.jl")
export signal_derivative
export band_power
export make_spectrum
export signal_detrend

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
export rms
export nexpow2
export vsplit


include("nstim.jl")
export tes_dose

end
