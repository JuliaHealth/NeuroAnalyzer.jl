module NeuroJ

import FFTW
import DSP
import Statistics

include("eeg.jl")
export signal_derivative
export band_power
export make_spectrum

include("misc.jl")
export linspace
export logspace
export zero_pad
export vsearch(x::Vector, y::Number)
export vsearch(x::Vector, y::Vector)
export cart2pol
export pol2cart
export minmax_scaler
export cvangle
export hann
export hildebrand_rule
export jaccard_similarity
export fft0
export next_power_of_2

include("nstim.jl")
export tes_dose

end