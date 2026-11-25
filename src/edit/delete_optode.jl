export delete_optode
export delete_optode!
export delete_pair
export delete_pair!

"""
    delete_optode(obj; <keyword arguments>)

Return a copy of `obj` with the specified NIRS optode(s) and all their associated channels removed.

Source or detector labels, optode-pair mappings, and channel location entries are all updated consistently. All channels using the optode(s) are removed (intensity, OD and concentration channels).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `opt::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}}`: optode number(s) to remove

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function delete_optode(
    obj::NeuroAnalyzer.NEURO;
    opt::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}},
)::NeuroAnalyzer.NEURO
    obj_new = deepcopy(obj)
    delete_optode!(obj_new; opt = opt)

    return obj_new
end

"""
    delete_optode!(obj; <keyword arguments>)

Delete optode(s) in-place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `opt::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}}`: optode number(s) to be removed

# Returns

- `Nothing`
"""
function delete_optode!(
    obj::NeuroAnalyzer.NEURO;
    opt::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}},
)::Nothing
    # validate
    _check_datatype(obj, "nirs")
    # highest index first, so that lower indices stay valid while deleting
    opt_vec = sort(unique(collect(opt isa Int64 ? [opt] : opt)); rev = true)
    opt_n   = length(obj.header.recording[:optode_labels])
    for idx in opt_vec
        idx in 1:opt_n ||
            throw(ArgumentError("opt index $idx is out of range [1, $opt_n]."))
    end
    length(opt_vec) < opt_n || throw(
        ArgumentError(
            "Number of optodes to delete ($(length(opt_vec))) must be less than " *
            "the total number of optodes ($opt_n).",
        ),
    )

    rec = obj.header.recording

    # remove channel locations
    for idx in opt_vec
        lbl = optode_labels(obj)[idx] # optode labels (header) are not changed in this loop
        loc_result = _find_bylabel(obj.locs, lbl)
        if loc_result isa Int64
            deleteat!(obj.locs, loc_result)
        elseif !isempty(loc_result)
            deleteat!(obj.locs, sort(loc_result))
        end
    end

    # update headers and build list of channels to delete
    # optode_pairs holds the *source* index (column 1) and the *detector* index (column 2),
    # i.e. positions in src_labels / det_labels, not positions in optode_labels
    chs_to_delete = Int64[]
    for idx in opt_vec
        ol = rec[:optode_labels][idx]
        deleteat!(rec[:optode_labels], idx)
        col, field = if ol in rec[:src_labels]
            1, :src_labels
        elseif ol in rec[:det_labels]
            2, :det_labels
        else
            continue
        end
        k = findfirst(isequal(ol), rec[field])
        chp = rec[:optode_pairs][:, col]
        append!(chs_to_delete, findall(isequal(k), chp))
        deleteat!(rec[field], k)
        chp[chp .== k] .= 0
        chp[chp .> k] .-= 1
        rec[:optode_pairs][:, col] = chp
    end

    _history!(obj, "delete_optode", (; opt))
    chs_to_delete = sort(unique(chs_to_delete))
    if !isempty(chs_to_delete)
        _info("Deleting NIRS channels: $(labels(obj)[chs_to_delete])")
        # by index: the channels are identified by their optode pairs, not by labels
        h = copy(obj.history)
        _delete_channel!(obj, chs_to_delete)
        obj.history = h
    end

    return nothing
end

"""
    delete_pair(obj; <keyword arguments>)

Return a copy of `obj` with all channels of the specified NIRS source-detector pair(s) removed: intensity, OD and concentration channels at all wavelengths. The optodes (labels and locations) are kept, as they may be used by other pairs; use `delete_optode()` to remove optodes.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `pair::Union{String, Vector{String}}`: pair(s) given as `"source_detector"` labels, e.g. `"S1_D2"`

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function delete_pair(
    obj::NeuroAnalyzer.NEURO;
    pair::Union{String, Vector{String}},
)::NeuroAnalyzer.NEURO
    obj_new = deepcopy(obj)
    delete_pair!(obj_new; pair = pair)

    return obj_new
end

"""
    delete_pair!(obj; <keyword arguments>)

Delete all channels of the specified NIRS source-detector pair(s) in-place. See `delete_pair()`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `pair::Union{String, Vector{String}}`: pair(s) given as `"source_detector"` labels, e.g. `"S1_D2"`

# Returns

- `Nothing`
"""
function delete_pair!(
    obj::NeuroAnalyzer.NEURO;
    pair::Union{String, Vector{String}},
)::Nothing
    # validate
    _check_datatype(obj, "nirs")
    pairs = pair isa String ? [pair] : pair
    isempty(pairs) && throw(ArgumentError("No pairs selected."))

    # channels described by :optode_pairs (source index, detector index) come first
    rec = obj.header.recording
    n_opt = size(rec[:optode_pairs], 1)
    ch_idx = Int64[]
    for p in pairs
        sd = split(p, "_")
        length(sd) == 2 || throw(
            ArgumentError(
                "Pair must be given as \"source_detector\", e.g. \"S1_D1\" (got \"$p\").",
            ),
        )
        s = findfirst(==(sd[1]), rec[:src_labels])
        d = findfirst(==(sd[2]), rec[:det_labels])
        isnothing(s) && throw(ArgumentError("Source $(sd[1]) not found."))
        isnothing(d) && throw(ArgumentError("Detector $(sd[2]) not found."))
        idx = [j for j = 1:n_opt if rec[:optode_pairs][j, :] == [s, d]]
        isempty(idx) && throw(ArgumentError("No channels of pair $p found."))
        append!(ch_idx, idx)
    end

    n0 = length(obj.history)
    _delete_channel!(obj, ch_idx)
    _history!(obj, "delete_pair", (; pair); from = n0)

    return nothing
end
