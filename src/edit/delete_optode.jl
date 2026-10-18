export delete_optode
export delete_optode!

"""
    delete_optode(obj; <keyword arguments>)

Return a copy of `obj` with the specified NIRS optode(s) and all their associated channels removed.

Source or detector labels, optode-pair mappings, and channel location entries are all updated consistently. Associated NIRS signal channels are removed via `delete_channel!`.

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

    # create new dataset
    obj_tmp = deepcopy(obj)
    rec = obj_tmp.header.recording

    # remove channel locations
    for idx in opt_vec
        lbl = optode_labels(obj)[idx] # use original obj - obj_tmp changes each iteration
        loc_result = _find_bylabel(obj_tmp.locs, lbl)
        if loc_result isa Int64
            deleteat!(obj_tmp.locs, loc_result)
        elseif !isempty(loc_result)
            deleteat!(obj_tmp.locs, sort(loc_result))
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

    push!(obj_tmp.history, "delete_optode(obj; opt=$opt)")
    chs_to_delete = labels(obj_tmp)[sort(unique(chs_to_delete))]
    isempty(chs_to_delete) || _info("Deleting NIRS channels: $chs_to_delete")
    isempty(chs_to_delete) || delete_channel!(obj_tmp; ch = chs_to_delete, del_opt = true)

    return obj_tmp
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
    obj_tmp = delete_optode(obj; opt = opt)
    obj.header = obj_tmp.header
    obj.data = obj_tmp.data
    obj.history = obj_tmp.history
    obj.locs = obj_tmp.locs
    obj_tmp = nothing

    return nothing
end
