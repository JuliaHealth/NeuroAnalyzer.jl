export amp

"""
    amp(s)

Computes amplitude descriptors.

# Arguments

- `s::AbstractVector`: signal vector

# Returns

Named tuple:

- `p::Float64`: peak amplitude (`max(|s|)`)
- `r::Float64`: RMS amplitude (`p / √2`; exact only for a pure sinusoid)
- `p2p::Float64`: peak-to-peak amplitude (`max(s) - min(s)`)
- `semi_p2p::Float64`: half of the peak-to-peak amplitude
- `msa::Float64`: mean square amplitude (`mean(s²)`)
- `rmsa::Float64`: root mean square amplitude (`p2p / √2`; exact only for a pure sinusoid)
- `es::Float64`: total signal energy (`Σ s²`)
- `rmsq::Float64`: root mean square (`√mean(s²)`)
"""
function amp(
    s::AbstractVector
)::@NamedTuple{
    p::Float64,
    r::Float64,
    p2p::Float64,
    semi_p2p::Float64,
    msa::Float64,
    rmsa::Float64,
    es::Float64,
    rmsq::Float64
}

    p = maximum(abs, s)
    r = p / sqrt(2)
    s_min, s_max = extrema(s)
    p2p = s_max - s_min
    semi_p2p = p2p / 2
    msa = sum(abs2, s) / length(s)
    rmsa = p2p / sqrt(2)
    es = sum(abs2, s)
    rmsq = rms(s)

    return (; p, r, p2p, semi_p2p, msa, rmsa, es, rmsq)

end

"""
    amp(s)

Computes amplitude descriptors.

# Arguments

- `s::AbstractArray`: signal array (channels, samples, epochs)

# Returns

Named tuple:

- `p::Matrix{Float64}`: peak amplitude (`max(|s|)`)
- `r::Matrix{Float64}`: RMS amplitude (`p / √2`; exact only for a pure sinusoid)
- `p2p::Matrix{Float64}`: peak-to-peak amplitude (`max(s) - min(s)`)
- `semi_p2p::Matrix{Float64}`: half of the peak-to-peak amplitude
- `msa::Matrix{Float64}`: mean square amplitude (`mean(s²)`)
- `rmsa::Matrix{Float64}`: root mean square amplitude (`p2p / √2`; exact only for a pure sinusoid)
- `es::Matrix{Float64}`: total signal energy (`Σ s²`)
- `rmsq::Matrix{Float64}`: root mean square (`√mean(s²)`)
"""
function amp(
    s::AbstractArray
)::@NamedTuple{
    p::Matrix{Float64},
    r::Matrix{Float64},
    p2p::Matrix{Float64},
    semi_p2p::Matrix{Float64},
    msa::Matrix{Float64},
    rmsa::Matrix{Float64},
    es::Matrix{Float64},
    rmsq::Matrix{Float64},
}

    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate output
    p = zeros(ch_n, ep_n)
    r = zeros(ch_n, ep_n)
    p2p = zeros(ch_n, ep_n)
    semi_p2p = zeros(ch_n, ep_n)
    msa = zeros(ch_n, ep_n)
    rmsa = zeros(ch_n, ep_n)
    es = zeros(ch_n, ep_n)
    rmsq = zeros(ch_n, ep_n)

    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        amp_data = amp(@view(s[ch_idx, :, ep_idx]))
        p[ch_idx, ep_idx] = amp_data.p
        r[ch_idx, ep_idx] = amp_data.r
        p2p[ch_idx, ep_idx] = amp_data.p2p
        semi_p2p[ch_idx, ep_idx] = amp_data.semi_p2p
        msa[ch_idx, ep_idx] = amp_data.msa
        rmsa[ch_idx, ep_idx] = amp_data.rmsa
        nrg[ch_idx, ep_idx] = amp_data.es
        rmsq[ch_idx, ep_idx] = amp_data.rmsq
    end

    return (; p, r, p2p, semi_p2p, msa, rmsa, es, rmsq)

end

"""
    amp(obj; <keyword arguments>)

Calculate amplitudes.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)

# Returns

Named tuple:

- `p::Matrix{Float64}`: peak amplitude
- `r::Matrix{Float64}`: RMS amplitude
- `p2p::Matrix{Float64}`: peak-to-peak amplitude
- `semi_p2p::Matrix{Float64}`: half of the peak-to-peak amplitude
- `msa::Matrix{Float64}`: mean square amplitude
- `rmsa::Matrix{Float64}`: root mean square amplitude
- `es::Matrix{Float64}`: total signal energy
- `rmsq::Matrix{Float64}`: root mean square
"""
function amp(
    obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}
)::@NamedTuple{
    p::Matrix{Float64},
    r::Matrix{Float64},
    p2p::Matrix{Float64},
    semi_p2p::Matrix{Float64},
    msa::Matrix{Float64},
    rmsa::Matrix{Float64},
    es::Matrix{Float64},
    rmsq::Matrix{Float64},
}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    return amp(@view(obj.data[ch, :, :]))

end
