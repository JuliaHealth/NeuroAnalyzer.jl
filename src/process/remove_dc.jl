export remove_dc

"""
    remove_dc(s)

Remove mean value (DC offset).

# Arguments

- `s::AbstractVector`

# Returns

- `s_new::Vector{Float64}`
"""
function remove_dc(s::AbstractVector)
    
    s_new = s .- mean(s)

    return s_new

end

"""
    remove_dc(s; ch)

Remove mean value (DC offset).

# Arguments

- `s::AbstractArray`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels

# Returns

- `s::Array{Float64, 3}`
"""
function remove_dc(s::AbstractArray)

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    s_new = similar(s)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            s[ch_idx, :, ep_idx] = remove_dc(s[ch_idx, :, ep_idx])
        end
    end

    return s_new

end

"""
    remove_dc(obj; ch)

Remove mean value (DC offset).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function remove_dc(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)))

    _check_channels(obj, ch)

    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = @views remove_dc(obj.data[ch, :, :])
    reset_components!(obj_new)
    push!(obj_new.header.history, "remove_dc(OBJ, ch=$ch)")

    return obj_new

end

"""
    remove_dc!(obj; ch)

Remove mean value (DC offset).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
"""
function remove_dc!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)))

    obj_tmp = remove_dc(obj, ch=ch)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    obj.components = obj_tmp.components

    return nothing

end
