export invert_polarity
export invert_polarity!

"""
    invert_polarity(obj; <keyword arguments>)

Invert polarity.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function invert_polarity(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
)::NeuroAnalyzer.NEURO

    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))

    # create new dataset
    obj_new = deepcopy(obj)

    obj_new.data[ch, :, :] = .- obj_new.data[ch, :, :]
    push!(obj_new.history, "invert_polarity(obj; ch=$ch)")

    return obj_new
end

"""
    invert_polarity!(obj; <keyword arguments>)

Invert polarity.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel(s) to invert

# Returns

- `Nothing`
"""
function invert_polarity!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
)::Nothing
    obj_new = invert_polarity(obj; ch = ch)
    obj.data = obj_new.data
    obj.history = obj_new.history

    return nothing
end
