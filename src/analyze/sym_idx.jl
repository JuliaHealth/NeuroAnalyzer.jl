export sym_idx

"""
    sym_idx(s)

Calculate signal symmetry index (ratio of positive to negative amplitudes). Perfectly symmetrical signal has symmetry of 1.0. Symmetry above 1.0 indicates there are more positive amplitudes.

# Arguments

- `s::AbstractVector`: signal vector

# Returns

- `sym::Float64`: symmetry index
"""
function sym_idx(s::AbstractVector)::Float64

    sym = sum(s .< 0) == 0 ? sum(s .>= 0) : sum(s .>= 0) / sum(s .< 0)

    return sym

end

"""
    sym_idx(s)

Calculate signal symmetry index (ratio of positive to negative amplitudes). Perfectly symmetrical signal has symmetry of 1.0. Symmetry above 1.0 indicates there are more positive amplitudes.

# Arguments

- `s::AbstractArray`

# Returns

- `sym::Matrix{Float64}`: symmetry index
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

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads :static for ch_idx in 1:ch_n
            sym[ch_idx, ep_idx] = @views sym_idx(s[ch_idx, :, ep_idx])
        end
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

- `sym::Matrix{Float64}`: symmetry index
"""
function sym_idx(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::Matrix{Float64}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    return @views sym_idx(obj.data[ch, :, :])

end
