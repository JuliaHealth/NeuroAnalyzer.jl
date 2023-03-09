export delete_channel
export delete_channel!
export keep_channel
export keep_channel!
export keep_channel_type
export keep_channel_type!

"""
    delete_channel(obj; channel)

Delete channel(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel number(s) to be removed

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function delete_channel(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange})

    typeof(channel) <: AbstractRange && (channel = collect(channel))
    ch_n = channel_n(obj)
    length(channel) > 1 && (channel = sort!(channel, rev=true))
    length(channel) == ch_n && throw(ArgumentError("You cannot delete all channels."))

    _check_channels(obj, channel)

    obj_new = deepcopy(obj)

    # update headers
    for idx in channel
        loc = findfirst(isequal(lowercase(obj_new.header.recording[:labels][idx])), lowercase.(string.obj_new.locs[!, :labels]))
        loc !== nothing && deleteat!( obj_new.locs, loc)
        if obj_new.header.recording[:data_type] === "eeg"
            deleteat!(obj_new.header.recording[:labels], idx)
            deleteat!(obj_new.header.recording[:channel_type], idx)
            deleteat!(obj_new.header.recording[:transducers], idx)
            deleteat!(obj_new.header.recording[:units], idx)
            deleteat!(obj_new.header.recording[:prefiltering], idx)
            deleteat!(obj_new.header.recording[:gain], idx)
        elseif obj_new.header.recording[:data_type] === "meg"
            deleteat!(obj_new.header.recording[:labels], idx)
            deleteat!(obj_new.header.recording[:channel_type], idx)
            deleteat!(obj_new.header.recording[:units], idx)
            deleteat!(obj_new.header.recording[:prefiltering], idx)
            deleteat!(obj_new.header.recording[:coils], idx)
            deleteat!(obj_new.header.recording[:magnetometers], idx)
            deleteat!(obj_new.header.recording[:gradiometers], idx)
            deleteat!(obj_new.header.recording[:gradiometers_axial], idx)
            deleteat!(obj_new.header.recording[:gradiometers_planar], idx)
        end
    end
    obj_new.header.recording[:ch_n] -= length(channel)

    # remove channel
    obj_new.data =obj_new.data[setdiff(1:end, (channel)), :, :]

    reset_components!(obj_new)
    push!(obj_new.header.history, "delete_channel(OBJ, $channel)")

    returnobj_new
end

"""
    delete_channel!(obj; channel)

Delete channel(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel number(s) to be removed
"""
function delete_channel!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange})

    typeof(channel) <: AbstractRange && (channel = collect(channel))
    ch_n = channel_n(obj)
    length(channel) > 1 && (channel = sort!(channel, rev=true))
    length(channel) == ch_n && throw(ArgumentError("You cannot delete all channels."))

    _check_channels(obj, channel)

    # update headers
    for idx in channel
        loc = findfirst(isequal(lowercase(obj.header.recording[:labels][idx])), lowercase.(string.(obj.locs[!, :labels])))
        loc !== nothing && deleteat!(obj.locs, loc)
        if obj.header.recording[:data_type] === "eeg"
            deleteat!(obj.header.recording[:labels], idx)
            deleteat!(obj.header.recording[:channel_type], idx)
            deleteat!(obj.header.recording[:transducers], idx)
            deleteat!(obj.header.recording[:units], idx)
            deleteat!(obj.header.recording[:prefiltering], idx)
            deleteat!(obj.header.recording[:gain], idx)
        elseif obj.header.recording[:data_type] === "meg"
            deleteat!(obj.header.recording[:labels], idx)
            deleteat!(obj.header.recording[:channel_type], idx)
            deleteat!(obj.header.recording[:units], idx)
            deleteat!(obj.header.recording[:prefiltering], idx)
            deleteat!(obj.header.recording[:coils], idx)
            deleteat!(obj.header.recording[:magnetometers], idx)
            deleteat!(obj.header.recording[:gradiometers], idx)
            deleteat!(obj.header.recording[:gradiometers_axial], idx)
            deleteat!(obj.header.recording[:gradiometers_planar], idx)
        end
    end
    obj.header[:ch_n] -= length(channel)

    # remove channel
    obj.data = obj.data[setdiff(1:end, (channel)), :, :]

    reset_components!(obj)
    push!(obj.header[:history], "delete_channel!(OBJ, $channel)")

    return nothing
end

"""
    keep_channel(obj; channel)

Keep channel(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel number(s) to keep

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function keep_channel(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange})

    typeof(channel) <: AbstractRange && (channel = collect(channel))
    _check_channels(obj, channel)

    ch_n = channel_n(obj)
    channels_to_remove = setdiff(collect(1:ch_n), channel)
    length(channels_to_remove) == ch_n && throw(ArgumentError("You cannot delete all channels."))

    return delete_channel(obj, channel=channels_to_remove)
end

"""
    keep_channel!(obj; channel)

Keep channel(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel number(s) to keep
"""
function keep_channel!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange})

    typeof(channel) <: AbstractRange && (channel = collect(channel))
    _check_channels(obj, channel)

    ch_n = channel_n(obj)
    channels_to_remove = setdiff(collect(1:ch_n), channel)
    length(channels_to_remove) == ch_n && throw(ArgumentError("You cannot delete all channels."))

    delete_channel!(obj, channel=channels_to_remove)
end

"""
    keep_channel_type(eeg; type)

Keep OBJ channels of `type` type.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `type::Symbol=:eeg`: type of channels to keep

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function keep_channel_type(obj::NeuroAnalyzer.NEURO; type::Symbol=:eeg)

    _check_var(type, [:all, :eeg, :meg, :ecg, :eog, :emg, :ref, :mrk], "type")

    channels_idx = Vector{Int64}()
    for idx in 1:channel_n(eeg, type=:all)
        obj.header[:channel_type][idx] == string(type) && push!(channels_idx, idx)
    end
    obj_new = keep_channel(eeg, channel=channels_idx)
    reset_components!(obj_new)
    push!(obj_new.header.recording.history, "keep_channel_type(OBJ, type=$type")

    return obj_new
end

"""
    keep_channel_type!(eeg; type)

Keep OBJ channels of `type` type.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `type::Symbol=:eeg`: type of channels to keep
"""
function keep_channel_type!(obj::NeuroAnalyzer.NEURO; type::Symbol=:eeg)

    obj_tmp = keep_channel_type(eeg, type=type)
    obj.header = obj_tmp.header
    obj.data = obj_tmp.data
    reset_components!(eeg)

    return nothing
end
