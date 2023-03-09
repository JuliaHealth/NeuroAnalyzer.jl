export derivative
export derivative!

"""
    derivative(signal)

Return derivative of the same length.

# Arguments

- `signal::AbstractVector`
"""
function derivative(signal::AbstractVector)
    s_der = diff(signal)
    return vcat(s_der, s_der[end])
end

"""
    derivative(obj; channel)

Return the derivative of channel(s) with length same as the signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function derivative(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj)))

    ep_n = epoch_n(obj)
    
    obj_new = deepcopy(obj)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in eachindex(channel)
            @views obj_new.data[channel[ch_idx], :, ep_idx] = derivative(obj_new.data[channel[ch_idx], :, ep_idx])
        end
    end

    reset_components!(obj_new)
    push!(obj_new.header.history, "derivative(OBJ, channel=$channel)")

    return obj_new
end

"""
    derivative!(obj; channel)

Return the derivative of channel(s) with length same as the signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
"""
function derivative!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj)))

    obj_tmp = derivative(obj, channel=channel)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    reset_components!(obj)

    return nothing
end
