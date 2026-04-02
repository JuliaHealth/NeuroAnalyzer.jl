export delete_channel
export delete_channel!
export keep_channel
export keep_channel!

"""
    delete_channel(obj; <keyword arguments>)

Return a copy of `obj` with the specified channel(s) removed.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s) to remove
- `del_opt::Bool=false`: set `true` only when called from `delete_optode()`; bypasses the NIRS guard that prevents direct signal-channel deletion

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function delete_channel(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    del_opt::Bool = false,
)::NeuroAnalyzer.NEURO
    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))

    # number of channels
    ch_n = nchannels(obj)

    # validate
    length(ch) > 1 && (ch = sort(ch; rev = true))
    length(ch) < ch_n ||
        throw(
            ArgumentError(
                "Number of channels to delete ($(length(ch))) must be smaller than number of all channels ($ch_n).",
            ),
        )

    # sort descending so deleteat! on earlier header vectors doesn't shift indices used by later iterations
    ch_idx = sort(ch; rev = true)

    # create new dataset
    obj_new = deepcopy(obj)

    (datatype(obj) == "meg" && size(obj.header.recording[:ssp_data]) != (0,)) && _warn(
        "OBJ contains SSP projections data, you should apply them before modifying OBJ data.",
    )

    # update headers
    for idx in ch_idx
        # remove matching location entry, if any
        loc_idx = _find_bylabel(obj_new.locs, labels(obj)[idx])
        if loc_idx isa Int64
            deleteat!(obj_new.locs, loc_idx)
        elseif !isempty(loc_idx)
            deleteat!(obj_new.locs, loc_idx[1])
        end
 
        # remove from universal header vectors
        deleteat!(obj_new.header.recording[:label],        idx)
        deleteat!(obj_new.header.recording[:channel_type], idx)
        deleteat!(obj_new.header.recording[:bad_channel],  idx)
        deleteat!(obj_new.header.recording[:unit],         idx)
 
        # remove from type-specific header vectors
        dt = obj_new.header.recording[:data_type]
        if dt == "eeg" || dt == "seeg" || dt == "ecog"
            deleteat!(obj_new.header.recording[:prefiltering], idx)
            deleteat!(obj_new.header.recording[:transducers],  idx)
            deleteat!(obj_new.header.recording[:gain],         idx)
        elseif dt == "meg"
            deleteat!(obj_new.header.recording[:prefiltering], idx)
            deleteat!(obj_new.header.recording[:coil_type],    idx)
            for field in (:gradiometers, :magnetometers)
                tmp = findfirst(isequal(idx), obj_new.header.recording[field])
                isnothing(tmp) || deleteat!(obj_new.header.recording[field], tmp)
            end
        elseif dt == "nirs"
            if !del_opt && idx in eachindex(obj_new.header.recording[:optode_labels])
                throw(ArgumentError(
                    "NIRS signal channels must be deleted using delete_optode().",
                ))
            end
            if idx in eachindex(obj_new.header.recording[:wavelength_index])
                deleteat!(obj_new.header.recording[:wavelength_index], idx)
            end
            chp1 = obj_new.header.recording[:optode_pairs][:, 1]
            chp2 = obj_new.header.recording[:optode_pairs][:, 2]
            if idx in axes(obj_new.header.recording[:optode_pairs], 1)
                deleteat!(chp1, idx)
                deleteat!(chp2, idx)
                obj_new.header.recording[:optode_pairs] = hcat(chp1, chp2)
            end
        end
    end
 
    obj_new.header.recording[:channel_order] = _sort_channels(
        obj_new.header.recording[:channel_type],
    )
    obj_new.data = obj_new.data[setdiff(_c(ch_n), ch_idx), :, :]
 
    push!(obj_new.history, "delete_channel(obj; ch=$(labels(obj)[sort(ch_idx)]))")
 
    return obj_new
end

"""
    delete_channel!(obj; <keyword arguments>)

Delete channel(s) in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel nanme(s) to remove
- `del_opt::Bool=false`: set `true` only when called from `delete_optode()`; bypasses the NIRS guard that prevents direct signal-channel deletion

# Returns

- `Nothing`
"""
function delete_channel!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    del_opt::Bool = false,
)::Nothing
    obj_new = delete_channel(obj; ch = ch, del_opt = del_opt)
    obj.header  = obj_new.header
    obj.data    = obj_new.data
    obj.history = obj_new.history
    obj.locs    = obj_new.locs

    return nothing
end

"""
    keep_channel(obj; <keyword arguments>)

Return a copy of `obj` retaining only the specified channel(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s) to keep

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function keep_channel(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
)::NeuroAnalyzer.NEURO
    # number of channels
    ch_n = nchannels(obj)

    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))
    length(ch) == ch_n && (return obj)

    chs_to_remove = labels(obj)[setdiff(_c(ch_n), ch)]
    length(chs_to_remove) < ch_n || throw(
        ArgumentError(
            "Number of channels to delete ($(length(chs_to_remove))) must be less than " *
            "the total number of channels ($ch_n).",
        ),
    )
 
    return delete_channel(obj; ch = chs_to_remove)
end

"""
    keep_channel!(obj; <keyword arguments>)

Keep only the specified channel(s) in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s) to keep

# Returns

- `Nothing`
"""
function keep_channel!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
)::Nothing
    obj_new = keep_channel(obj; ch = ch)
    obj.header  = obj_new.header
    obj.data    = obj_new.data
    obj.history = obj_new.history
    obj.locs    = obj_new.locs

    return nothing
end
