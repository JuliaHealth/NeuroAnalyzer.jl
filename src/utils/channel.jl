export get_channel

"""
    get_channel(obj; <keyword arguments>)

Return list of channel names of specified type or their numbers if names are specified.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}=""`: channel name or list of channel names
- `type::Union{String, Vector{String}}="all"`: channels types
- `wl::Real`: return NIRS channels for wavelength (in nm)
- `exclude::Union{String, Vector{String}, Regex}=""`: channel name or list of channel names to exclude from the list

# Returns

- `ch::Union{Vector{String}, Vector{Int64}}`
"""
function get_channel(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex} = "",
    type::Union{String, Vector{String}} = "all",
    wl::Real = 0,
    exclude::Union{String, Vector{String}, Regex} = "",
)::Union{Vector{String}, Vector{Int64}}
    if ch != ""
        exclude = _ch_idx(obj, exclude)
        ch = _ch_idx(obj, ch)
        if isnothing(exclude)
            return sort(ch)
        else
            chs = setdiff(ch, exclude)
            return sort(chs)
        end
    end

    # return channel names
    ch = String[]
    isa(type, String) && (type = [type])
    [_check_var(idx, channel_types, "type") for idx in type]

    l = labels(obj)
    if wl == 0
        if type == ["all"]
            ch = l
        else
            for type_idx in eachindex(type)
                for ch_idx in eachindex(obj.header.recording[:channel_type])
                    obj.header.recording[:channel_type][ch_idx] == type[type_idx] &&
                        push!(ch, l[ch_idx])
                end
            end
        end
    else
        _check_datatype(obj, ["nirs"])
        wl in obj.header.recording[:wavelengths] ||
            throw(
                ArgumentError(
                    "OBJ does not contain data for $wl wavelength. Available wavelengths: $(obj.header.recording[:wavelengths]).",
                ),
            )
        wl_idx = findfirst(isequal(wl), obj.header.recording[:wavelengths])
        for ch_idx in eachindex(obj.header.recording[:wavelength_index])
            obj.header.recording[:wavelength_index][ch_idx] == wl_idx &&
                push!(ch, l[ch_idx])
        end
    end

    exclude = exclude == "" ? Int64[] : _ch_idx(obj, exclude)
    ch = exclude == [] ? ch : setdiff(ch, labels(obj)[exclude])

    chs = unique(ch)
    return sort(chs)
end
