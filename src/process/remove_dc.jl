export remove_dc

"""
    remove_dc(signal)

Remove mean value (DC offset).

# Arguments

- `signal::AbstractVector`

# Returns

- `s_demeaned::Vector{Float64}`
"""
function remove_dc(signal::AbstractVector)
    return signal .- mean(signal)
end

"""
    remove_dc(obj; channel)

Remove mean value (DC offset).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function remove_dc(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj)))

    _check_channels(obj, channel)

    ep_n = epoch_n(obj)

    obj_new = deepcopy(obj)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in eachindex(channel)
            @views obj_new.data[channel[ch_idx], :, ep_idx] = remove_dc(obj_new.data[channel[ch_idx], :, ep_idx])
        end
    end

    reset_components!(obj_new)
    push!(obj_new.header.history, "remove_dc(OBJ, channel=$channel)")

    return obj_new
end

"""
    remove_dc!(obj; channel)

Remove mean value (DC offset).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
"""
function remove_dc!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj)))

    obj_tmp = remove_dc(obj, channel=channel)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    obj.components = obj_tmp.components

    return nothing
end
