export hfd

"""
    _hfd(s)

Calculate the Higuchi fractal dimension for a 1-D signal vector.

# Arguments

- `s::AbstractVector`: signal vector

# Returns

- `Float64`: the Higuchi fractal dimension

# Notes

Wraps FractalDimensions.higuchi_dim() to compute the Higuchi fractal dimension (Higuchi, 1988) for neural signals.

The Higuchi FD estimates the fractal dimension directly from the time series by measuring how the length of the signal scales with the lag parameter k. Values range from 1 (smooth/linear) to 2 (maximally irregular/random).
"""
function _hfd(s::AbstractVector)::Float64
    return higuchi_dim(s)
end

"""
    _hfd(s)

Calculate the Higuchi fractal dimension for a 3-D signal array.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)

# Returns

- `Matrix{Float64}`: the Higuchi fractal dimension, shape (channels, epochs)

# Notes

Wraps FractalDimensions.higuchi_dim() to compute the Higuchi fractal dimension (Higuchi, 1988) for neural signals.

The Higuchi FD estimates the fractal dimension directly from the time series by measuring how the length of the signal scales with the lag parameter k. Values range from 1 (smooth/linear) to 2 (maximally irregular/random).
"""
function _hfd(s::AbstractArray)::Matrix{Float64}
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
        hd[ch_idx, ep_idx] = _hfd(@view(s[ch_idx, :, ep_idx]))
    end

    return hd
end

"""
    hfd(obj; <keyword arguments>)

Calculate the Higuchi fractal dimension for a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `exclude_bad_segments::Bool=NeuroAnalyzer.exclude_bad_segments`: if `true`, bad segments (markers of type `"BAD"` of all channels or of the selected channels) of a continuous signal are not used: their samples are removed and the remaining parts are joined (default: the `exclude_bad_segments` preference, see `set_exclude_bad_segments()`)
- `exclude_bad_channels::Bool=NeuroAnalyzer.exclude_bad_channels`: if `true`, channels marked as bad are not used (default: the `exclude_bad_channels` preference, see `set_exclude_bad_channels()`)

# Returns

- `Matrix{Float64}`

# Notes

Wraps FractalDimensions.higuchi_dim() to compute the Higuchi fractal dimension (Higuchi, 1988) for neural signals.

The Higuchi FD estimates the fractal dimension directly from the time series by measuring how the length of the signal scales with the lag parameter k. Values range from 1 (smooth/linear) to 2 (maximally irregular/random).
"""
function hfd(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    exclude_bad_channels::Bool = NeuroAnalyzer.exclude_bad_channels,
    exclude_bad_segments::Bool = NeuroAnalyzer.exclude_bad_segments,
)::Matrix{Float64}
    # resolve channel names to integer indices, optionally skipping bad channels
    ch = _get_ch(obj, ch, exclude_bad_channels)
    isempty(ch) && throw(ArgumentError("No channels selected."))

    return _hfd(_seg_data(obj, ch, exclude_bad_segments))
end
