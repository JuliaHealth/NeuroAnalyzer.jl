export add_signal
export add_signal!

"""
    add_signal(s1, s2)

Add two equal-length signal vectors element-wise.

# Arguments

- `s1::AbstractVector`: target signal
- `s2::AbstractVector`: signal to add to `s1`; must have the same length as `s1`

# Returns

- `AbstractVector`: element-wise sum `s1 .+ s2`

# Throws

- `ArgumentError`: if `s1` and `s2` have different lengths

# See also

[`add_signal(::NeuroAnalyzer.NEURO)`](@ref)
"""
function add_signal(s1::AbstractVector, s2::AbstractVector)::AbstractVector

    length(s1) == length(s2) || throw(ArgumentError("s1 and s2 must have the same length."))

    return s1 .+ s2

end

"""
    add_signal(obj; <keyword arguments>)

Add a signal vector to selected channels of a NEURO object.

The same signal `s` is added to every selected channel in every epoch. `s` must have the same length as one epoch of the signal (`epoch_len(obj)`). Channels are processed in parallel using `Threads.@threads`.


# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `s::AbstractVector`: signal to add; length must equal `epoch_len(obj)`

# Returns

- `NeuroAnalyzer.NEURO`: new object with `s` added to the selected channels

# Throws

- `ArgumentError`: if `length(s) ≠ epoch_len(obj)`

# See also

[`add_signal!`](@ref), [`add_signal(::AbstractVector, ::AbstractVector)`](@ref)
"""
function add_signal(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    s::AbstractVector
)::NeuroAnalyzer.NEURO

    # validate s length against epoch length before any allocation
    length(s) == epoch_len(obj) || throw(ArgumentError("Length of s must equal epoch_len(obj) ($(epoch_len(obj)))."))

    # resolve channel names to integer indices
    ch = get_channel(obj, ch = ch)
    # number of channels
    ch_n = length(ch)
    # number of epochs
    ep_n = nepochs(obj)

    obj_new = deepcopy(obj)

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        obj_new.data[ch[ch_idx], :, ep_idx] =
            add_signal(@view(obj.data[ch[ch_idx], :, ep_idx]), s)
    end

    push!(obj_new.history, "add_signal(OBJ, ch=$ch)")

    return obj_new

end

"""
    add_signal!(obj; <keyword arguments>)

Add a signal vector to selected channels of a NEURO object in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; modified in-place
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `s::AbstractVector`: signal to add; length must equal `epoch_len(obj)`

# Returns

- `Nothing`

# Throws
- `ArgumentError`: if `length(s) ≠ epoch_len(obj)`

# See also

[`add_signal`](@ref)
"""
function add_signal!(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, s::AbstractVector)::Nothing

    obj_new = add_signal(obj, ch = ch, s = s)
    obj.data = obj_new.data
    obj.history = obj_new.history

    return nothing

end
