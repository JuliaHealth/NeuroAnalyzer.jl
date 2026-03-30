export sym_idx

"""
    sym_idx(s)

Calculate signal symmetry index (ratio of positive to negative amplitudes) for a 1-D signal vector.

Perfectly symmetrical signal has symmetry of 1.0. Symmetry above 1.0 indicates there are more positive amplitudes.

# Arguments

- `s::AbstractVector`: signal vector

# Returns

- `Float64`: symmetry index
"""
function sym_idx(s::AbstractVector)::Float64
    sym = sum(s .< 0) == 0 ? sum(s .>= 0) : sum(s .>= 0) / sum(s .< 0)

    return sym
end

"""
    sym_idx(s)

Calculate signal symmetry index (ratio of positive to negative amplitudes) for a 3-D signal array.

Perfectly symmetrical signal has symmetry of 1.0. Symmetry above 1.0 indicates there are more positive amplitudes.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)

# Returns

- `Matrix{Float64}`: symmetry index
"""
function sym_idx(s::AbstractArray)::Matrix{Float64}

    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate output
    sym = zeros(ch_n, ep_n)

    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        sym[ch_idx, ep_idx] = sym_idx(@view(s[ch_idx, :, ep_idx]))
    end

    return sym
end

"""
    sym_idx(obj; <keyword arguments>)

Calculate signal symmetry index (ratio of positive to negative amplitudes). Perfectly symmetrical signal has symmetry of 1.0. Symmetry above 1.0 indicates there are more positive amplitudes.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)

# Returns

- `Matrix{Float64}`: symmetry index
"""
function sym_idx(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
)::Matrix{Float64}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch =
        exclude_bads ?
        get_channel(obj; ch = ch, exclude = "bad") :
        get_channel(obj; ch = ch, exclude = "")

    return sym_idx(@view(obj.data[ch, :, :]))
end
