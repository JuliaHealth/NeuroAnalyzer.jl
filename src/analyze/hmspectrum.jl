export hmspectrum

"""
    hmspectrum(obj; ch)

Calculate Hilbert marginal spectrum.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::String`: channel name

# Returns

Named tuple containing:
- `p::Array{Float64, 3}`: Hilbert marginal spectra for each epoch
- `t::Vector{Float64}`: time points
"""
function hmspectrum(obj; ch::String)::@NamedTuple{p::Array{Float64, 3}, t::Vector{Float64}}

    p, _, t = NeuroAnalyzer.spectrogram(obj, ch=ch, method=:hht, db=false)
    p = dropdims(sum(p, dims=1), dims=3)

    return (p=p, t=t)

end
