export derivative
export derivative!

"""
    derivative(s)

Return derivative of the same length.

# Arguments

- `s::AbstractVector`

# Returns

- `s_new::AbstractVector`
"""
function derivative(s::AbstractVector)

    s_der = diff(s)
    
    return vcat(s_der, s_der[end])

end

"""
    derivative(s)

Return derivative of the same length.

# Arguments

- `s::AbstractArray`

# Returns

- `s_new::Array{Float64, 3}`
"""
function derivative(s::AbstractArray)

    ch_n = size(s, 1)
    ep_n = size(s, 3)
    
    s_new = similar(s)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s_new[ch_idx, :, ep_idx] = @views derivative(s[ch_idx, :, ep_idx])
        end
    end

    return s_new

end

"""
    derivative(obj; ch)

Return derivative of the same length.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function derivative(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)))

    _check_channels(obj, ch)

    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = derivative(obj.data[ch, :, :])
    reset_components!(obj_new)
    push!(obj_new.header.history, "derivative(OBJ, ch=$ch)")

    return obj_new

end

"""
    derivative!(obj; ch)

Return derivative of the same length.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
"""
function derivative!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)))

    obj_tmp = derivative(obj, ch=ch)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    obj.components = obj_tmp.components

    return nothing

end
