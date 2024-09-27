export standardize
export standardize!

"""
    standardize(s)

Standardize channels.

# Arguments

- `s::AbstractArray`

# Returns

- `s_new::Array{Float64, 3}`:
- `scaler::Vector{ZScoreTransform{Float64, Vector{Float64}}}`
"""
function standardize(s::AbstractArray)::Tuple{Array{Float64, 3}, Vector{ZScoreTransform{Float64, Vector{Float64}}}}

    _chk3d(s)
    ep_n = size(s, 3)

    scaler = Vector{ZScoreTransform{Float64, Vector{Float64}}}()

    s_new = similar(s)
    @inbounds for ep_idx in 1:ep_n
        @views push!(scaler, StatsBase.fit(ZScoreTransform, s[:, :, ep_idx], dims=2))
        @views s_new[:, :, ep_idx] = StatsBase.transform(scaler[ep_idx], s[:, :, ep_idx])
    end

    return s_new, scaler

end

"""
    standardize(obj; <keyword arguments>)

Standardize channels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}}`: channel name or list of channel names

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
- `scaler::Vector{ZScoreTransform{Float64, Vector{Float64}}}`
"""
function standardize(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}})::Tuple{NeuroAnalyzer.NEURO, Vector{ZScoreTransform{Float64, Vector{Float64}}}}

    ch = get_channel(obj, ch=ch)
    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :], scaler = standardize(obj.data[ch, :, :])
    reset_components!(obj_new)
    push!(obj_new.history, "standardize(OBJ)")

    return obj_new, scaler

end

"""
    standardize!(obj; <keyword arguments>)

Standardize channels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}}`: channel name or list of channel names

# Returns

- `scaler::Vector{ZScoreTransform{Float64, Vector{Float64}}}`
"""
function standardize!(obj::NeuroAnalyzer.NEURO)::Vector{ZScoreTransform{Float64, Vector{Float64}}}

    obj_new, scaler = standardize(obj, ch=ch)
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components

    return scaler

end
