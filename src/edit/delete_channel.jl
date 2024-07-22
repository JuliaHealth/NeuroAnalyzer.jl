export delete_channel
export delete_channel!
export keep_channel
export keep_channel!

"""
    delete_channel(obj; <keyword arguments>)

Delete channel(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}}`: channels to be removed
- `del_opt::Bool=false`: for NIRS data is set as `true` if called from `remove_optode()`

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function delete_channel(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}}, del_opt::Bool=false)

    ch_n = nchannels(obj)
    ch = _ch_idx(obj, ch)
    length(ch) > 1 && (ch = sort!(ch, rev=true))
    @assert length(ch) < ch_n "Number of channels to delete ($(length(ch))) must be smaller than number of all channels ($ch_n)."
    obj_new = deepcopy(obj)

    # update headers
    for idx in ch
        # remove channel locations
        loc_idx = _find_bylabel(obj_new.locs, labels(obj)[idx])
        !isnothing(loc_idx) && deleteat!(obj_new.locs, loc_idx)
        deleteat!(obj_new.header.recording[:label], idx)
        deleteat!(obj_new.header.recording[:channel_type], idx)
        obj_new.header.recording[:bad_channel] = obj_new.header.recording[:bad_channel][1:end .!= idx, :]
        deleteat!(obj_new.header.recording[:unit], idx)
        if obj_new.header.recording[:data_type] == "eeg"
            deleteat!(obj_new.header.recording[:prefiltering], idx)
            deleteat!(obj_new.header.recording[:transducers], idx)
            deleteat!(obj_new.header.recording[:gain], idx)
        elseif obj_new.header.recording[:data_type] == "meg"
            deleteat!(obj_new.header.recording[:prefiltering], idx)
            deleteat!(obj_new.header.recording[:coils], idx)
            idx_tmp = findfirst(isequal(idx), obj_new.header.recording[:gradiometers])
            !isnothing(idx_tmp) && deleteat!(obj_new.header.recording[:gradiometers], idx_tmp)
            idx_tmp = findfirst(isequal(idx), obj_new.header.recording[:magnetometers])
            !isnothing(idx_tmp) && deleteat!(obj_new.header.recording[:magnetometers], idx_tmp)
            deleteat!(obj_new.header.recording[:coil_type], idx)
        elseif obj_new.header.recording[:data_type] == "nirs"
            if !del_opt && idx in 1:length(obj_new.header.recording[:optode_labels])
                @warn "NIRS signal channels must be deleted using delete_optode()."
                return nothing
            end
            idx in 1:length(obj_new.header.recording[:wavelength_index]) && deleteat!(obj_new.header.recording[:wavelength_index], idx)
            chp1 = obj_new.header.recording[:optode_pairs][:, 1]
            chp2 = obj_new.header.recording[:optode_pairs][:, 2]
            if idx in 1:size(obj_new.header.recording[:optode_pairs], 1)
                deleteat!(chp1, idx)
                deleteat!(chp2, idx)
                obj_new.header.recording[:optode_pairs] = hcat(chp1, chp2)
            end
        end
    end

    obj_new.header.recording[:channel_order] = _sort_channels(obj_new.header.recording[:channel_type])

    # remove channel
    obj_new.data = obj_new.data[setdiff(_c(ch_n), ch), :, :]

    reset_components!(obj_new)
    push!(obj_new.history, "delete_channel(OBJ, ch=$ch)")

    return obj_new

end

"""
    delete_channel!(obj; <keyword arguments>)

Delete channels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}}`: channels to be removed
- `del_opt::Bool=false`: for NIRS data is set as `true` if called from `remove_optode()`
"""
function delete_channel!(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}}, del_opt::Bool=false)

    obj_new = delete_channel(obj, ch=ch, del_opt=del_opt)
    obj.header = obj_new.header
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components
    obj.locs = obj_new.locs

    return nothing

end

"""
    keep_channel(obj; <keyword arguments>)

Keep channels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}}`: channels to keep

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function keep_channel(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}})

    ch_n = nchannels(obj)
    chs_to_remove = labels(obj)[setdiff(_c(ch_n), _ch_idx(obj, ch))]
    @assert length(chs_to_remove) < ch_n "Number of channels to delete ($(length(chs_to_remove))) must be smaller than number of all channels ($ch_n)."

    obj_new = delete_channel(obj, ch=chs_to_remove)

    return obj_new

end

"""
    keep_channel!(obj; <keyword arguments>)

Keep channels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}}`: channels to keep
"""
function keep_channel!(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}})

    obj_new = keep_channel(obj, ch=ch)
    obj.header = obj_new.header
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components
    obj.locs = obj_new.locs

    return nothing

end
