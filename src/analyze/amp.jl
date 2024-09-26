export amp

"""
    amp(s)

Calculate amplitudes.

# Arguments

- `s::AbstractVector`

# Returns

Named tuple containing:
- `p::Float64`: peak amplitude
- `r::Float64`: RMS amplitude
- `p2p::Float64`: peak-to-peak amplitude
- `semi_p2p::Float64`: half of the peak-to-peak amplitude
- `msa::Float64`: mean square amplitude
- `rmsa::Float64`: root mean square amplitude
- `energy::Float64`: total signal energy
- `rms::Float64`: root mean square
"""
function amp(s::AbstractVector)::NamedTuple{(:p, :r, :p2p, :semi_p2p, :msa, :rmsa, :energy, :rms), Tuple{Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64}}

    p = maximum(abs.(s))
    r = p / sqrt(2)
    p2p = abs(maximum(s)) + abs(minimum(s))
    semi_p2p = p2p / 2
    msa = 1/length(s) * sum(s.^2)
    rmsa = p2p / sqrt(2)
    nrg = sum(s.^2)
    # rms = sqrt(1/length(s) * sum(s.^2))
    rms = norm(s) / sqrt(length(s))

    return (p=p, r=r, p2p=p2p, semi_p2p=semi_p2p, msa=msa, rmsa=rmsa, energy=nrg, rms=rms)

end

"""
    amp(s)

Calculate amplitudes.

# Arguments

- `s::AbstractArray`

# Returns

Named tuple containing:
- `p::Matrix{Float64}`: peak amplitude
- `r::Matrix{Float64}`: RMS amplitude
- `p2p::Matrix{Float64}`: peak-to-peak amplitude
- `semi_p2p::Matrix{Float64}`: half of the peak-to-peak amplitude
- `msa::Matrix{Float64}`: mean square amplitude
- `rmsa::Matrix{Float64}`: root mean square amplitude
- `energy::Matrix{Float64}`: total signal energy
"""
function amp(s::AbstractArray)::NamedTuple{(:p, :r, :p2p, :semi_p2p, :msa, :rmsa, :energy, :rms), Tuple{Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}}}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)

    p = zeros(ch_n, ep_n)
    r = zeros(ch_n, ep_n)
    p2p = zeros(ch_n, ep_n)
    semi_p2p = zeros(ch_n, ep_n)
    msa = zeros(ch_n, ep_n)
    rmsa = zeros(ch_n, ep_n)
    nrg = zeros(ch_n, ep_n)
    rms = zeros(ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        for ch_idx in 1:ch_n
            p[ch_idx, ep_idx], r[ch_idx, ep_idx], p2p[ch_idx, ep_idx], semi_p2p[ch_idx, ep_idx], msa[ch_idx, ep_idx], rmsa[ch_idx, ep_idx], nrg[ch_idx, ep_idx], rms[ch_idx, ep_idx] = @views amp(s[ch_idx, :, ep_idx])
        end
    end

    return (p=p, r=r, p2p=p2p, semi_p2p=semi_p2p, msa=msa, rmsa=rmsa, energy=nrg, rms=rms)

end

"""
    amp(obj; <keyword arguments>)

Calculate amplitudes.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}}`: channel name or list of channel names

# Returns

Named tuple containing:
- `p::Matrix{Float64}`: peak amplitude
- `r::Matrix{Float64}`: RMS amplitude
- `p2p::Matrix{Float64}`: peak-to-peak amplitude
- `semi_p2p::Matrix{Float64}`: half of the peak-to-peak amplitude
- `msa::Matrix{Float64}`: mean square amplitude
- `rmsa::Matrix{Float64}`: root mean square amplitude
- `energy::Matrix{Float64}`: total signal energy
"""
function amp(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}})::NamedTuple{(:p, :r, :p2p, :semi_p2p, :msa, :rmsa, :energy, :rms), Tuple{Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}}}

    ch = get_channel(obj, ch=ch)

    p, r, p2p, semi_p2p, msa, rmsa, nrg, rms = @views amp(obj.data[ch, :, :])

    return (p=p, r=r, p2p=p2p, semi_p2p=semi_p2p, msa=msa, rmsa=rmsa, energy=nrg, rms=rms)

end
