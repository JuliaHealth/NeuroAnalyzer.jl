module NeuroJ

include("eeg.jl")
export signal_derivative
export band_power
export make_spectrum

include("misc.jl")
export zero_pad
export vsearch
export cart2pol
export pol2cart
export minmax_scaler
export cvangle
export hann
export hildebrand_rule
export jaccard_similarity
export manhattan_distance
export fft0
export simpson

include("nstim.jl")
export tes_dose

end