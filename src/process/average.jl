export average

"""
    average(s)

Average all channels.

# Arguments

- `s::AbstractArray`

# Returns

- `average::AbstractArray`
"""
function average(s::AbstractArray)::AbstractArray

    return mean(s, dims=1)

end

"""
    average(s1, s2)

Averages two signals.

# Arguments

- `s1::AbstractArray`
- `s2::AbstractArray`

# Returns

- `average::AbstractArray`
"""
function average(s1::AbstractArray, s2::AbstractArray)::AbstractArray

    return mean(hcat(s1, s2), dims=2)

end

"""
    average(obj; <keyword arguments>)

Return the average signal of channels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function average(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::NeuroAnalyzer.NEURO

    ch = get_channel(obj, ch=ch)
    obj_new = deepcopy(obj)
    keep_channel!(obj_new, ch=labels(obj)[1])
    obj_new.data = @views average(obj.data[ch, :, :])
    obj_new.header.recording[:label]=["averaged ch"]
    reset_components!(obj_new)
    push!(obj_new.history, "average(OBJ, ch=$ch)")

    return obj_new

end

"""
    average!(obj; <keyword arguments>)

Return the average signal of channels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names

# Returns

Nothing
"""
function average!(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::Nothing

    obj_new = average(obj, ch=ch)
    obj.header = obj_new.header
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    return nothing

end

"""
    average(obj1, obj2)

Return the average signal of two objects.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function average(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO)::NeuroAnalyzer.NEURO

    @assert size(obj1.data) == size(obj2.data) "Both signals must have the same size."

    ch_n = nchannels(obj1)
    ep_n = nepochs(obj1)

    obj_new = deepcopy(obj1)
    @inbounds for ep_idx in 1:ep_n
        Threads.@threads :static for ch_idx in 1:ch_n
            obj_new.data[ch_idx, :, ep_idx] = @views average(obj1.data[ch_idx, :, ep_idx], obj2.data[ch_idx, :, ep_idx])
        end
    end

    reset_components!(obj_new)
    push!(obj_new.history, "average(OBJ1, OBJ2)")

    return obj_new

end
