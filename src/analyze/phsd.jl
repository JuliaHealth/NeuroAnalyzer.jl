export phsd

"""
    phsd(s; <keyword arguments>)

Calculate phase spectral density for a 1-D signal vector.

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
    fs::Int64,
)::@NamedTuple{
    ph::Vector{Float64},
    f::Vector{Float64},
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
    fs::Int64,
)::@NamedTuple{
    ph::Matrix{Float64},
    f::Vector{Float64},
}
    ch_n = size(s, 1)
    phsd_data = phsd(s[1, :]; fs = fs)
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

Calculate phase spectral density for a 3-D signal array.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `fs::Int64`: sampling rate in Hz; must be ≥ 1

# Returns

Named tuple:

- `ph::Array{Float64, 3}`: phases (in radians)
- `f::Vector{Float64}`: frequencies
"""
function phsd(
    s::AbstractArray;
    fs::Int64,
)::@NamedTuple{
    ph::Array{Float64, 3},
    f::Vector{Float64},
}

    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    f = phsd(s[1, :, 1]; fs = fs).f

    ph = zeros(ch_n, length(f), ep_n)

    # calculate over channel and epochs
    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        phsd_data = phsd(@view(s[ch_idx, :, ep_idx]), fs = fs)
        ph[ch_idx, :, ep_idx] = phsd_data.ph
    end

    return (; ph, f)
end

"""
    phsd(obj; <keyword arguments>)

Calculate phase spectral density for a NEURO object.

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
    ch::Union{String, Vector{String}, Regex},
)::@NamedTuple{
    ph::Array{Float64, 3},
    f::Vector{Float64},
}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ?
        get_channel(obj; ch = ch, exclude = "bad") :
        get_channel(obj; ch = ch, exclude = "")

    return phsd(@view(obj.data[ch, :, :]); fs = sr(obj))
end
