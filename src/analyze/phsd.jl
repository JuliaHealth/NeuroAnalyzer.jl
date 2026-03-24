export phsd

"""
    phsd(s; <keyword arguments>)

Calculate phase spectral density.

# Arguments

- `s::Vector{Float64}`
- `fs::Int64`: sampling rate in Hz; must be ≥ 1

# Returns

Named tuple:

- `ph::Vector{Float64}`: phases (in radians)
- `f::Vector{Float64}`: frequencies
"""
function phsd(
    s::AbstractVector;
    fs::Int64
)::@NamedTuple{
    ph::Vector{Float64},
    f::Vector{Float64}
}

    # validate
    fs >= 1 || throw(ArgumentError("fs must be ≥ 1."))

    ft_data = NeuroAnalyzer.ftransform(s)
    ph = ft_data.ph
    f = freqs(s, fs)[1]

    return (; ph, f)

end

"""
    phsd(s; <keyword arguments>)

Calculate phase spectral density.

# Arguments

- `s::AbstractMatrix`
- `fs::Int64`: sampling rate in Hz; must be ≥ 1

# Returns

Named tuple:

- `ph::Matrix{Float64}`: phases (in radians)
- `f::Vector{Float64}`: frequencies
"""
function phsd(
    s::AbstractMatrix;
    fs::Int64
)::@NamedTuple{
    ph::Matrix{Float64},
    f::Vector{Float64}
}

    ch_n = size(s, 1)
    phsd_data = phsd(s[1, :], fs = fs)
    f = phsd_data.f

    ph = zeros(ch_n, length(f))

    @inbounds for ch_idx in 1:ch_n
        phsd_data = phsd(@view(s[ch_idx, :]), fs = fs)
        ph[ch_idx, :] = phsd_data.ph
    end

    return (; ph, f)

end

"""
    phsd(s; <keyword arguments>)

Calculate phase spectral density.

# Arguments

- `s::AbstractArray`
- `fs::Int64`: sampling rate in Hz; must be ≥ 1

# Returns

Named tuple:

- `ph::Array{Float64, 3}`: phases (in radians)
- `f::Vector{Float64}`: frequencies
"""
function phsd(
    s::AbstractArray;
    fs::Int64
)::@NamedTuple{
    ph::Array{Float64, 3},
    f::Vector{Float64}
}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)

    phsd_data = phsd(s[1, :, 1], fs = fs)
    f = phsd_data.f

    ph = zeros(ch_n, length(f), ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads :dynamic for ch_idx in 1:ch_n
            phsd_data = phsd(@view(s[ch_idx, :, ep_idx]), fs = fs)
            ph[ch_idx, :, ep_idx] = phsd_data.ph
        end
    end

    return (; ph, f)

end

"""
    phsd(obj; <keyword arguments>)

Calculate phase spectral density.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)

# Returns

Named tuple:

- `ph::Array{Float64, 3}`: phases (in radians)
- `f::Vector{Float64}`: frequencies
"""
function phsd(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex}
)::@NamedTuple{
    ph::Array{Float64, 3},
    f::Vector{Float64}
}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    return phsd(@view(obj.data[ch, :, :]), fs = sr(obj))

end
