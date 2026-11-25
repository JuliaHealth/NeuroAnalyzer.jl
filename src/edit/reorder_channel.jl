export reorder_channel
export reorder_channel!

"""
    reorder_channel(obj; <keyword arguments>)

Return a copy of `obj` with the channels in a given order: the data and all per-channel header fields are reordered.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `order::Union{Vector{String}, Symbol}`: new order:
    - channel labels, e.g. the order of a montage; labels that are not channels of `obj` are ignored (with a warning); labels are matched exactly or, if there is no exact match, ignoring case
    - `:type`: by channel type (the order used by the plots, see `channel_order()`)
    - `:locs`: the order of the channel locations (`obj.locs`, e.g. loaded from a montage file with `load_locs()`)
- `keep_rest::Bool=true`: channels not in `order` are placed after the ordered channels (in their current order); if `false` they are deleted

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object

# Notes

Channel locations and markers refer to the channels by label and are not changed. For NIRS objects, the channels described by the optode pairs (intensity or OD) must stay the first channels.
"""
function reorder_channel(
    obj::NeuroAnalyzer.NEURO;
    order::Union{Vector{String}, Symbol},
    keep_rest::Bool = true,
)::NeuroAnalyzer.NEURO
    obj_new = deepcopy(obj)
    reorder_channel!(obj_new; order = order, keep_rest = keep_rest)

    return obj_new
end

"""
    reorder_channel!(obj; <keyword arguments>)

Put the channels in a given order in-place: the data and all per-channel header fields are reordered.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `order::Union{Vector{String}, Symbol}`: new order:
    - channel labels, e.g. the order of a montage; labels that are not channels of `obj` are ignored (with a warning); labels are matched exactly or, if there is no exact match, ignoring case
    - `:type`: by channel type (the order used by the plots, see `channel_order()`)
    - `:locs`: the order of the channel locations (`obj.locs`, e.g. loaded from a montage file with `load_locs()`)
- `keep_rest::Bool=true`: channels not in `order` are placed after the ordered channels (in their current order); if `false` they are deleted

# Returns

- `Nothing`

# Notes

Channel locations and markers refer to the channels by label and are not changed. For NIRS objects, the channels described by the optode pairs (intensity or OD) must stay the first channels.
"""
function reorder_channel!(
    obj::NeuroAnalyzer.NEURO;
    order::Union{Vector{String}, Symbol},
    keep_rest::Bool = true,
)::Nothing
    idx = _reorder_idx(obj, order)
    isempty(idx) && throw(ArgumentError("None of the channels in order exists."))
    ch_n = nchannels(obj)
    rest = setdiff(1:ch_n, idx)

    n0 = length(obj.history)
    if keep_rest
        idx = vcat(idx, rest)
    elseif !isempty(rest)
        # deleted as by delete_channel() (locations, markers); indices of the others shift
        lbl = labels(obj)
        _delete_channel!(obj, rest)
        new = cumsum(.!in(Set(rest)).(1:ch_n))
        idx = new[idx]
        _info("Channels not in order deleted: $(lbl[rest])")
    end
    idx == 1:nchannels(obj) || _select_channels!(obj, idx)

    _history!(obj, "reorder_channel", (; order, keep_rest); from = n0)

    return nothing
end

# channel indices in the requested order
function _reorder_idx(obj::NeuroAnalyzer.NEURO, order::Union{Vector{String}, Symbol})::Vector{Int64}
    clabels = labels(obj)
    if order isa Symbol
        _check_var(order, [:type, :locs], "order")
        order === :type && return copy(channel_order(obj))
        _has_locs(obj)
        order = String.(obj.locs[!, :label])
    end
    allunique(order) || throw(ArgumentError("order contains duplicated labels."))

    idx = Int64[]
    missing_lbl = String[]
    for l in order
        i = findall(==(l), clabels)
        isempty(i) && (i = findall(c -> lowercase(c) == lowercase(l), clabels))
        if isempty(i)
            push!(missing_lbl, l)
        else
            length(i) == 1 ||
                throw(ArgumentError("Label $l matches several channels: $(clabels[i])."))
            i[1] in idx && throw(ArgumentError("Channel $(clabels[i[1]]) is given twice."))
            push!(idx, i[1])
        end
    end
    isempty(missing_lbl) || _warn("Labels not found in the channels (ignored): $missing_lbl")

    return idx
end
