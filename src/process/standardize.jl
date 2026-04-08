export standardize
export standardize!

"""
    standardize(s)

Standardize channels of a 3-D signal array.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)

# Returns

- `Array{Float64, 3}`
- `Vector{ZScoreTransform{Float64, Vector{Float64}}}`
"""
function standardize(
    s::AbstractArray,
)::Tuple{Array{Float64, 3}, Vector{ZScoreTransform{Float64, Vector{Float64}}}}
    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of epochs
    ep_n = size(s, 3)

    scaler = Vector{ZScoreTransform{Float64, Vector{Float64}}}()

    s_new = similar(s, Float64)
    @inbounds for ep_idx in 1:ep_n
        push!(scaler, StatsBase.fit(ZScoreTransform, @view(s[:, :, ep_idx]), dims = 2))
        s_new[:, :, ep_idx] = StatsBase.transform(scaler[ep_idx], @view(s[:, :, ep_idx]))
    end

    return s_new, scaler
end

"""
    standardize(obj; <keyword arguments>)

Standardize channels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
- `Vector{ZScoreTransform{Float64, Vector{Float64}}}`
"""
function standardize(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
)::Tuple{NeuroAnalyzer.NEURO, Vector{ZScoreTransform{Float64, Vector{Float64}}}}
    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))

    # create new dataset
    obj_new = deepcopy(obj)

    obj_new.data[ch, :, :], scaler = standardize(obj.data[ch, :, :])
    push!(obj_new.history, "standardize(obj)")

    return obj_new, scaler
end

"""
    standardize!(obj; <keyword arguments>)

Standardize channels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)

# Returns

- `Vector{ZScoreTransform{Float64, Vector{Float64}}}`
"""
function standardize!(
    obj::NeuroAnalyzer.NEURO,
)::Vector{ZScoreTransform{Float64, Vector{Float64}}}
    obj_new, scaler = standardize(obj; ch = ch)
    obj.data = obj_new.data
    obj.history = obj_new.history

    return scaler
end
