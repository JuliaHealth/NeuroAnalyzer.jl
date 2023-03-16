export delete_channel
export delete_channel!
export keep_channel
export keep_channel!
export keep_channel_type
export keep_channel_type!

"""
    delete_channel(obj; ch)

Delete channel(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}`: channel number(s) to be removed

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function delete_channel(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange})

    typeof(ch) <: AbstractRange && (ch = collect(ch))
    ch_n = channel_n(obj)
    length(ch) > 1 && (ch = sort!(ch, rev=true))
    length(ch) == ch_n && throw(ArgumentError("You cannot delete all channels."))

    _check_channels(obj, ch)

    obj_new = deepcopy(obj)

    # update headers
    for idx in ch
        loc = findfirst(isequal(lowercase(obj_new.header.recording[:labels][idx])), lowercase.(string.(obj_new.locs[!, :labels])))
        loc !== nothing && deleteat!(obj_new.locs, loc)
        deleteat!(obj_new.header.recording[:labels], idx)
        deleteat!(obj_new.header.recording[:channel_type], idx)
        deleteat!(obj_new.header.recording[:units], idx)
        deleteat!(obj_new.header.recording[:prefiltering], idx)
        if obj_new.header.recording[:data_type] === "eeg"
            deleteat!(obj_new.header.recording[:transducers], idx)
            deleteat!(obj_new.header.recording[:gain], idx)
        elseif obj_new.header.recording[:data_type] === "meg"
            deleteat!(obj_new.header.recording[:coils], idx)
            deleteat!(obj_new.header.recording[:magnetometers], idx)
            deleteat!(obj_new.header.recording[:gradiometers], idx)
            deleteat!(obj_new.header.recording[:gradiometers_axial], idx)
            deleteat!(obj_new.header.recording[:gradiometers_planar], idx)
        end
    end

    # remove channel
    obj_new.data =obj_new.data[setdiff(1:end, (ch)), :, :]

    reset_components!(obj_new)
    push!(obj_new.history, "delete_channel(OBJ, ch=$ch)")

    return obj_new

end

"""
    delete_channel!(obj; ch)

Delete channel(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}`: channel number(s) to be removed
"""
function delete_channel!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange})

    obj_new = delete_channel(obj, ch=ch)
    obj.data = obj_new.data
    obj.header = obj_new.header
    obj.history = obj_new.history

    return nothing

end

"""
    keep_channel(obj; ch)

Keep channel(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}`: channel number(s) to keep

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function keep_channel(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange})

    typeof(ch) <: AbstractRange && (ch = collect(ch))
    _check_channels(obj, ch)

    ch_n = channel_n(obj)
    chs_to_remove = setdiff(collect(1:ch_n), ch)
    length(chs_to_remove) == ch_n && throw(ArgumentError("You cannot delete all channels."))

    obj_new = delete_channel(obj, ch=chs_to_remove)

    return obj_new

end

"""
    keep_channel!(obj; ch)

Keep channel(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}`: channel number(s) to keep
"""
function keep_channel!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange})

    obj_new = keep_channel(obj, ch=ch)
    obj.data = obj_new.data
    obj.header = obj_new.header
    obj.history = obj_new.history

    return nothing

end

"""
    keep_channel_type(obj; type)

Keep channel(s) of `type` type.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `type::Symbol=:eeg`: type of channels to keep

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function keep_channel_type(obj::NeuroAnalyzer.NEURO; type::Symbol=:eeg)

    _check_var(type, [:all, :eeg, :meg, :ecg, :eog, :emg, :ref, :mrk], "type")

    chs_idx = Vector{Int64}()
    for idx in 1:channel_n(obj, type=:all)
        obj.header.recording[:channel_type][idx] == string(type) && push!(chs_idx, idx)
    end

    obj_new = keep_channel(obj, ch=chs_idx)

    return obj_new

end

"""
    keep_channel_type!(obj; type)

Keep OBJ channels of `type` type.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `type::Symbol=:eeg`: type of channels to keep
"""
function keep_channel_type!(obj::NeuroAnalyzer.NEURO; type::Symbol=:eeg)

    obj_new = keep_channel_type(obj, type=type)
    obj.data = obj_new.data
    obj.header = obj_new.header
    obj.history = obj_new.history

    return nothing

end
