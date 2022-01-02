module NeuroJ

include("eeg.jl")
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
export signal_derivative
export simpson
export band_power

include("nstim.jl")
export tes_dose

end