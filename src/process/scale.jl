export scale
export scale!

"""
    scale(obj; <keyword arguments>)

Multiply channel(s) by `factor`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `factor::Real`: signal is multiplied by `factor`

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function scale(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    factor::Real,
)::NeuroAnalyzer.NEURO

    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)

    # create new dataset
    obj_new = deepcopy(obj)

    obj_new.data[ch, :, :] = @view(obj_new.data[ch, :, :]) .* factor
    push!(obj_new.history, "scale(obj; ch=$ch, factor=$factor)")

    return obj_new
end

"""
    scale!(obj; <keyword arguments>)

Multiply channel(s) by `factor`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `factor::Real`: signal is multiplied by `factor`

# Returns

- `Nothing`
"""
function scale!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    factor::Real,
)::Nothing
    obj_new = scale(obj; ch = ch, factor = factor)
    obj.data = obj_new.data
    obj.history = obj_new.history

    return nothing
end
