export hmspectrum

"""
    hmspectrum(obj; <keyword arguments>)

Calculate Hilbert marginal spectrum for a NEURO object.

The Hilbert marginal spectrum is computed from the Hilbert-Huang Transform (HHT) spectrogram by integrating (summing) the instantaneous power over time for each frequency, independently per epoch. The result is a frequency × epochs power matrix analogous to a classical power spectrum.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel name; must resolve to exactly one channel

# Returns

Named tuple:

- `p::Matrix{Float64}`: Hilbert marginal spectra, shape (frequency, epochs)
- `f::Vector{Float64}`: frequencies

# References

Huang et al. (1998), "The empirical mode decomposition and the Hilbert spectrum for nonlinear and non-stationary time series analysis."
"""
function hmspectrum(
    obj::NeuroAnalyzer.NEURO;
    ch::String
)::@NamedTuple{
    p::Matrix{Float64},
    f::Vector{Float64}
}

    # validate
    length(get_channel(obj; ch=ch)) == 1 || throw(ArgumentError("ch must resolve to exactly one channel."))

    # compute HHT time-frequency spectrogram with dB normalization
    spec = NeuroAnalyzer.spectrogram(obj; ch = ch, method = :hht, db = false)

    p = dropdims(spec.p, dims = 3)
    p = dropdims(sum(p, dims = 2), dims = 2)
    f = spec.f

    return (; p, f)

end
