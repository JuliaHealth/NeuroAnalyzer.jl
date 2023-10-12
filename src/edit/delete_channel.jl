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
    ch_n = nchannels(obj)
    length(ch) > 1 && (ch = sort!(ch, rev=true))
    @assert length(ch) < ch_n "Number of channels to delete ($(length(ch))) must be smaller than number of all channels ($ch_n)."

    _check_channels(obj, ch)

    obj_new = deepcopy(obj)

    # remove channel locations
    for ch_idx in ch
        if labels(obj_new)[ch_idx] in obj_new.locs[!, :labels]
            if length(NeuroAnalyzer._find_bylabel(obj_new.locs, labels(obj_new)[ch_idx])) == 1
                deleteat!(obj_new.locs, NeuroAnalyzer._find_bylabel(obj_new.locs, labels(obj_new)[ch_idx]))
            else
                deleteat!(obj_new.locs, sort(NeuroAnalyzer._find_bylabel(obj_new.locs, labels(obj_new)[ch_idx])))
            end
        end
    end
    
    # update headers
    for idx in ch
        loc = findfirst(isequal(lowercase(obj_new.header.recording[:labels][idx])), lowercase.(string.(obj_new.locs[!, :labels])))
        deleteat!(obj_new.header.recording[:labels], idx)
        deleteat!(obj_new.header.recording[:channel_type], idx)
        deleteat!(obj_new.header.recording[:units], idx)
        if obj_new.header.recording[:data_type] == "eeg"
            idx in obj_new.header.recording[:prefiltering] && deleteat!(obj_new.header.recording[:prefiltering], idx)
            idx in obj_new.header.recording[:transducers] && deleteat!(obj_new.header.recording[:transducers], idx)
            idx in obj_new.header.recording[:gain] && deleteat!(obj_new.header.recording[:gain], idx)
        elseif obj_new.header.recording[:data_type] == "meg"
            idx in obj_new.header.recording[:prefiltering] && deleteat!(obj_new.header.recording[:prefiltering], idx)
            idx in obj_new.header.recording[:coils] && deleteat!(obj_new.header.recording[:coils], idx)
            idx in obj_new.header.recording[:magnetometers] && deleteat!(obj_new.header.recording[:magnetometers], idx)
            idx in obj_new.header.recording[:gradiometers] && deleteat!(obj_new.header.recording[:gradiometers], idx)
            idx in obj_new.header.recording[:gradiometers_axial] && deleteat!(obj_new.header.recording[:gradiometers_axial], idx)
            idx in obj_new.header.recording[:gradiometers_planar] && deleteat!(obj_new.header.recording[:gradiometers_planar], idx)
        elseif obj_new.header.recording[:data_type] == "nirs"
            idx in 1:length(obj_new.header.recording[:wavelength_index]) && (deleteat!(obj_new.header.recording[:wavelength_index], idx))
            if idx in obj_new.header.recording[:channel_pairs]
                chp1 = obj_new.header.recording[:channel_pairs][:, 1]
                chp2 = obj_new.header.recording[:channel_pairs][:, 2]
                for chp_idx in size(chp1, 1):-1:1
                    if idx in chp1[chp_idx] == idx || idx in chp2[chp_idx] == idx
                        deleteat!(chp1, chp_idx)
                        deleteat!(chp2, chp_idx)
                    end
                end
                obj_new.header.recording[:channel_pairs] = hcat(chp1, chp2)
            end
            _warn("TO DO: remove optode_labels if contains removed channel")
            # deleteat!(obj_new.header.recording[:optode_labels], idx)
        end
    end

    # remove channel
    obj_new.data = obj_new.data[setdiff(_c(ch_n), ch), :, :]

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
    obj.header = obj_new.header
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components
    obj.locs = obj_new.locs

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

    ch_n = nchannels(obj)
    chs_to_remove = setdiff(_c(ch_n), ch)
    @assert length(chs_to_remove) < ch_n "Number of channels to delete ($(length(chs_to_remove))) must be smaller than number of all channels ($ch_n)."

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
    obj.header = obj_new.header
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components
    obj.locs = obj_new.locs

    return nothing

end

"""
    keep_channel_type(obj; type)

Keep channel(s) of `type` type.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `type::String="eeg"`: type of channels to keep

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function keep_channel_type(obj::NeuroAnalyzer.NEURO; type::String="eeg")

    _check_var(type, channel_types, "type")

    chs_idx = Vector{Int64}()
    for idx in 1:nchannels(obj, type="all")
        obj.header.recording[:channel_type][idx] == type && push!(chs_idx, idx)
    end

    obj_new = keep_channel(obj, ch=chs_idx)

    return obj_new

end

"""
    keep_channel_type!(obj; type)

Keep OBJ channels of `type` type.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `type::String="eeg"`: type of channels to keep
"""
function keep_channel_type!(obj::NeuroAnalyzer.NEURO; type::String="eeg")

    obj_new = keep_channel_type(obj, type=type)
    obj.header = obj_new.header
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components
    obj.locs = obj_new.locs

    return nothing

end
