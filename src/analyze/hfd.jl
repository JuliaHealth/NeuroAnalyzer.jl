export hfd

"""
    hfd(s)

Calculate the Higuchi fractal dimension for a 1-D signal vector.

# Arguments

- `s::AbstractVector`: signal vector

# Returns

- `Float64`: the Higuchi fractal dimension

# Notes

Wraps FractalDimensions.higuchi_dim() to compute the Higuchi fractal dimension (Higuchi, 1988) for neural signals.

The Higuchi FD estimates the fractal dimension directly from the time series by measuring how the length of the signal scales with the lag parameter k. Values range from 1 (smooth/linear) to 2 (maximally irregular/random).
"""
function hfd(s::AbstractVector)::Float64
    return higuchi_dim(s)
end

"""
    hfd(s)

Calculate the Higuchi fractal dimension for a 3-D signal array.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)

# Returns

- `Matrix{Float64}`: the Higuchi fractal dimension, shape (channels, epochs)

# Notes

Wraps FractalDimensions.higuchi_dim() to compute the Higuchi fractal dimension (Higuchi, 1988) for neural signals.

The Higuchi FD estimates the fractal dimension directly from the time series by measuring how the length of the signal scales with the lag parameter k. Values range from 1 (smooth/linear) to 2 (maximally irregular/random).
"""
function hfd(s::AbstractArray)::Matrix{Float64}

    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate output
    hd = zeros(ch_n, ep_n)

    # calculate over channel and epochs
    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        hd[ch_idx, ep_idx] = hfd(@view(s[ch_idx, :, ep_idx]))
    end

    return hd
end

"""
    hfd(obj; <keyword arguments>)

Calculate the Higuchi fractal dimension for a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)

# Returns

- `Matrix{Float64}`

# Notes

Wraps FractalDimensions.higuchi_dim() to compute the Higuchi fractal dimension (Higuchi, 1988) for neural signals.

The Higuchi FD estimates the fractal dimension directly from the time series by measuring how the length of the signal scales with the lag parameter k. Values range from 1 (smooth/linear) to 2 (maximally irregular/random).
"""
function hfd(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
)::Matrix{Float64}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch =
        exclude_bads ?
        get_channel(obj; ch = ch, exclude = "bad") :
        get_channel(obj; ch = ch, exclude = "")
    isempty(ch) && throw(ArgumentError("No channels selected."))

    return hfd(@view(obj.data[ch, :, :]))
end
