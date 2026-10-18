export phs

"""
    _phs(s; <keyword arguments>)

Calculate the phase spectrum (one-sided FFT phase) of a 1-D signal vector.

# Arguments

- `s::AbstractVector`: signal vector
- `fs::Int64`: sampling rate in Hz; must be ≥ 1

# Returns

Named tuple:

- `ph::Vector{Float64}`: phases in radians, in (-π, π]
- `f::Vector{Float64}`: frequencies (0 Hz to Nyquist)

# Notes

Phases of bins with near-zero amplitude are numerically meaningless.
"""
function _phs(
    s::AbstractVector;
    fs::Int64,
)::@NamedTuple{
    ph::Vector{Float64},
    f::Vector{Float64},
}
    # validate
    fs >= 1 || throw(ArgumentError("fs must be ≥ 1."))
    length(s) >= 2 || throw(ArgumentError("s must contain at least 2 samples."))

    ph = NeuroAnalyzer.ftransform(s).ph
    f = _freqs(s, fs)[1]

    return (; ph, f)
end

"""
    _phs(s; <keyword arguments>)

Calculate the phase spectrum for each channel of a matrix.

# Arguments

- `s::AbstractMatrix`: signal matrix (channels, samples)
- `fs::Int64`: sampling rate in Hz; must be ≥ 1

# Returns

Named tuple:

- `ph::Matrix{Float64}`: phases in radians, shape (channels, frequencies)
- `f::Vector{Float64}`: frequencies
"""
function _phs(
    s::AbstractMatrix;
    fs::Int64,
)::@NamedTuple{
    ph::Matrix{Float64},
    f::Vector{Float64},
}
    ch_n = size(s, 1)
    ch_n >= 1 || throw(ArgumentError("s must contain at least one channel."))

    # first channel also determines the frequency vector
    phs_data = _phs(@view(s[1, :]); fs = fs)
    f = phs_data.f
    ph = zeros(ch_n, length(f))
    ph[1, :] = phs_data.ph

    Threads.@threads :static for ch_idx = 2:ch_n
        ph[ch_idx, :] = _phs(@view(s[ch_idx, :]); fs = fs).ph
    end

    return (; ph, f)
end

"""
    _phs(s; <keyword arguments>)

Calculate the phase spectrum for a 3-D signal array.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `fs::Int64`: sampling rate in Hz; must be ≥ 1

# Returns

Named tuple:

- `ph::Array{Float64, 3}`: phases in radians, shape (channels, frequencies, epochs)
- `f::Vector{Float64}`: frequencies
"""
function _phs(
    s::AbstractArray;
    fs::Int64,
)::@NamedTuple{
    ph::Array{Float64, 3},
    f::Vector{Float64},
}
    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    # first channel/epoch also determines the frequency vector
    phs_data = _phs(@view(s[1, :, 1]); fs = fs)
    f = phs_data.f
    ph = zeros(ch_n, length(f), ep_n)
    ph[1, :, 1] = phs_data.ph

    # calculate over channels and epochs (skipping the one already done)
    idxs = vec(CartesianIndices((ch_n, ep_n)))[2:end]
    Threads.@threads :static for idx in idxs
        ch_idx, ep_idx = idx[1], idx[2]
        ph[ch_idx, :, ep_idx] = _phs(@view(s[ch_idx, :, ep_idx]); fs = fs).ph
    end

    return (; ph, f)
end

"""
    phs(obj; <keyword arguments>)

Calculate the phase spectrum for a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)

# Returns

Named tuple:

- `ph::Array{Float64, 3}`: phases in radians, shape (channels, frequencies, epochs)
- `f::Vector{Float64}`: frequencies
"""
function phs(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
)::@NamedTuple{
    ph::Array{Float64, 3},
    f::Vector{Float64},
}
    # resolve channel names to integer indices; skip bad channels if the global
    # `exclude_bads` preference is set
    ch =
        exclude_bads ?
        get_channel(obj; ch = ch, exclude = "bad") :
        get_channel(obj; ch = ch, exclude = "")
    isempty(ch) && throw(ArgumentError("No channels selected."))

    return _phs(@view(obj.data[ch, :, :]); fs = sr(obj))
end