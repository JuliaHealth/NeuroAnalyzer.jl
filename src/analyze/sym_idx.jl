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

    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)

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

    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")
    sym = @views sym_idx(obj.data[ch, :, :])

    return sym

end
