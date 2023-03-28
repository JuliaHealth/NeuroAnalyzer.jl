export standardize
export standardize!

"""
    standardize(s)

Standardize channels.

# Arguments

- `s::AbstractArray`

# Returns

- `s_new::NeuroAnalyzer.NEURO`:
- `scaler::Vector{Any}`: standardizing matrix
"""
function standardize(s::AbstractArray)
    
    ep_n = size(s, 3)

    scaler = Vector{ZScoreTransform{Float64, Vector{Float64}}}()

    s_new = similar(s)
    @inbounds @simd for ep_idx in 1:ep_n
        @views push!(scaler, StatsBase.fit(ZScoreTransform, s[:, :, ep_idx], dims=2)) 
        @views s_new[:, :, ep_idx] = StatsBase.transform(scaler[ep_idx], s[:, :, ep_idx])
    end

    return s_new, scaler

end

"""
    standardize(obj; ch)

Standardize channels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
- `scaler::Vector{Any}`: standardizing matrix
"""
function standardize(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj))
    
    _check_channels(obj, ch)

    obj_new = deepcopy(obj)
    obj_new.data, scaler = standardize(obj.data[ch, :, :])
    reset_components!(obj_new)
    push!(obj_new.history, "standardize(OBJ)")

    return obj_new, scaler

end

"""
    standardize!(obj; ch)

Standardize channels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

# Returns

- `scaler::Matrix{Float64}`: standardizing matrix
"""
function standardize!(obj::NeuroAnalyzer.NEURO)

    obj_new, scaler = standardize(obj, ch=ch)
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components

    return scaler

end
