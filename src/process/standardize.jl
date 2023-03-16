export standardize
export standardize!

"""
    standardize(obj; channel)

Standardize channels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

# Returns

- `new::NeuroAnalyzer.NEURO`: standardized OBJ
- `scaler::Matrix{Float64}`: standardizing matrix
"""
function standardize(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj))
    
    _check_channels(obj, channel)
    ep_n = epoch_n(obj)

    scaler = Vector{Any}()

    obj_new = deepcopy(obj)
    @inbounds @simd for ep_idx in 1:ep_n
        @views push!(scaler, StatsBase.fit(ZScoreTransform, obj.data[channel, :, ep_idx], dims=2)) 
        @views obj_new.data[channel,:, ep_idx] = StatsBase.transform(scaler[ep_idx], obj.data[channel, :, ep_idx])
    end

    reset_components!(obj_new)
    push!(obj_new.history, "standardize(OBJ)")

    return obj_new, scaler

end

"""
    standardize!(obj; channel)

Standardize channels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

# Returns

- `scaler::Matrix{Float64}`: standardizing matrix
"""
function standardize!(obj::NeuroAnalyzer.NEURO)

    obj_new, scaler = standardize(obj, channel=channel)
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components

    return scaler

end
