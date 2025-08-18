export phsd

"""
    phsd(s; <keyword arguments>)

Calculate phase spectral density.

# Arguments
- `s::Vector{Float64}`
- `fs::Int64`: sampling rate

# Returns

Named tuple containing:
- `ph::Vector{Float64}`: phases (in radians)
- `f::Vector{Float64}`: frequencies
"""
function phsd(s::AbstractVector; fs::Int64)::@NamedTuple{ph::Vector{Float64}, f::Vector{Float64}}

    @assert fs >= 1 "fs must be ≥ 1."

    _, _, _, sph = NeuroAnalyzer.ftransform(s)
    f, _ = freqs(s, fs)

    return (ph=sph, f=f)

end

"""
    phsd(s; <keyword arguments>)

Calculate phase spectral density.

# Arguments

- `s::AbstractMatrix`
- `fs::Int64`: sampling rate

# Returns

Named tuple containing:
- `ph::Matrix{Float64}`: phases (in radians)
- `f::Vector{Float64}`: frequencies
"""
function phsd(s::AbstractMatrix; fs::Int64)::@NamedTuple{ph::Matrix{Float64}, f::Vector{Float64}}

    ch_n = size(s, 1)
    _, f = phsd(s[1, :], fs=fs)

    ph = zeros(ch_n, length(f))

    @inbounds for ch_idx in 1:ch_n
        ph[ch_idx, :], _ = phsd(s[ch_idx, :], fs=fs)
    end

    return (ph=ph, f=f)

end

"""
    phsd(s; <keyword arguments>)

Calculate phase spectral density.

# Arguments
- `s::AbstractArray`
- `fs::Int64`: sampling rate

# Returns

Named tuple containing:
- `ph::Array{Float64, 3}`: phases (in radians)
- `f::Vector{Float64}`: frequencies
"""
function phsd(s::AbstractArray; fs::Int64)::@NamedTuple{ph::Array{Float64, 3}, f::Vector{Float64}}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)

    _, f = phsd(s[1, :, 1], fs=fs)

    ph = zeros(ch_n, length(f), ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads :greedy for ch_idx in 1:ch_n
            ph[ch_idx, :, ep_idx], _ = phsd(s[ch_idx, :, ep_idx], fs=fs)
        end
    end

    return (ph=ph, f=f)

end

"""
    phsd(obj; <keyword arguments>)

Calculate phase spectral density.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names

# Returns

Named tuple containing:
- `ph::Array{Float64, 3}`: phases (in radians)
- `f::Vector{Float64}`: frequencies
"""
function phsd(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::@NamedTuple{ph::Array{Float64, 3}, f::Vector{Float64}}

    ch = exclude_bads ? get_channel(obj, ch=ch, exclude="bad") : get_channel(obj, ch=ch, exclude="")
    ph, f = phsd(obj.data[ch, :, :], fs=sr(obj))

    return (ph=ph, f=f)

end
