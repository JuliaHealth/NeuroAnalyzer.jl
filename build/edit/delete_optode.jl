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

    opt_vec = sort(collect(opt isa Int64 ? [opt] : opt); rev = true)
    opt_n   = length(obj.header.recording[:optode_labels])

    for idx in opt_vec
        idx in 1:opt_n || throw(ArgumentError("opt index $idx is out of range [1, $opt_n]."))
    end
    length(opt_vec) < opt_n || throw(
        ArgumentError(
            "Number of optodes to delete ($(length(opt_vec))) must be less than " *
            "the total number of optodes ($opt_n).",
        ),
    )

    # create new dataset
    obj_new = deepcopy(obj)

    # remove channel locations
    for idx in opt_vec
        lbl    = optode_labels(obj)[idx] # use original obj - obj_new changes each iteration
        loc_result = _find_bylabel(obj_new.locs, lbl)
        if loc_result isa Int64
            deleteat!(obj_new.locs, loc_result)
        elseif !isempty(loc_result)
            deleteat!(obj_new.locs, sort(loc_result))
        end
    end

    # update headers and build list of channels to delete
    chs_to_delete = Int64[]
    for idx in opt_vec
        ol = optode_labels(obj_new)[idx]
        deleteat!(obj_new.header.recording[:optode_labels], idx)

        if ol in source_labels(obj_new)
            chp = obj_new.header.recording[:optode_pairs][:, 1]
            append!(chs_to_delete, findall(isequal(idx), chp))
            deleteat!(
                obj_new.header.recording[:src_labels],
                obj_new.header.recording[:src_labels] .== ol,
            )
            chp[chp .== idx] .= 0
            chp[chp .>  idx] .-= 1
            obj_new.header.recording[:optode_pairs][:, 1] = chp

        elseif ol in detector_labels(obj_new)
            chp = obj_new.header.recording[:optode_pairs][:, 2]
            append!(chs_to_delete, findall(isequal(idx), chp))
            deleteat!(
                obj_new.header.recording[:det_labels],
                obj_new.header.recording[:det_labels] .== ol,
            )
            chp[chp .== idx] .= 0
            chp[chp .>  idx] .-= 1
            obj_new.header.recording[:optode_pairs][:, 2] = chp
        end
    end

    push!(obj_new.history, "delete_optode(obj; opt=$opt)")

    chs_to_delete = labels(obj_new)[sort(unique(chs_to_delete))]
    isempty(chs_to_delete) || _info("Deleting NIRS channels: $chs_to_delete")
    isempty(chs_to_delete) || delete_channel!(obj_new; ch = chs_to_delete, del_opt = true)

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
    obj_new = delete_optode(obj; opt = opt)
    obj.header  = obj_new.header
    obj.data    = obj_new.data
    obj.history = obj_new.history
    obj.locs    = obj_new.locs

    return nothing
end
