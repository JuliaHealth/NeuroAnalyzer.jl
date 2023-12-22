export delete_optode
export delete_optode!

"""
    delete_optode(obj; opt, delete_channels)

Delete optode(s) and channels associated with removed optodes.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `opt::Union{Int64, Vector{Int64}, <:AbstractRange}`: optode number(s) to be removed

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function delete_optode(obj::NeuroAnalyzer.NEURO; opt::Union{Int64, Vector{Int64}, <:AbstractRange})

    _check_datatype(obj, "nirs")

    typeof(opt) <: AbstractRange && (opt = collect(opt))
    opt_n = length(obj.header.recording[:optode_labels])
    length(opt) > 1 && (opt = sort!(opt, rev=true))
    @assert length(opt) < opt_n "Number of optodes to delete ($(length(opt))) must be smaller than number of all optodes ($opt_n)."
    @assert opt in 1:opt_n "Opt must be in [1, $opt_n]."

    obj_new = deepcopy(obj)

    # remove channel locations
    for idx in opt
        if optode_labels(obj_new)[idx] in obj_new.locs[!, :labels]
            if length(_find_bylabel(obj_new.locs, optode_labels(obj)[idx])) == 1
                deleteat!(obj_new.locs, _find_bylabel(obj_new.locs, optode_labels(obj)[idx]))
            else
                deleteat!(obj_new.locs, sort(_find_bylabel(obj_new.locs, optode_labels(obj)[idx])))
            end
        end
    end

    # update headers
    chs_to_delete = Int64[]
    for idx in opt
        ol = optode_labels(obj_new)[idx]
        deleteat!(obj_new.header.recording[:optode_labels], idx)
        if ol in source_labels(obj_new)
            chp = obj_new.header.recording[:optode_pairs][:, 1]
            chs_to_delete = vcat(chs_to_delete, findall(isequal(idx), chp))
            deleteat!(obj_new.header.recording[:src_labels], obj_new.header.recording[:src_labels] .== ol)
            chp[chp .== idx] .= 0
            chp[chp .> idx] .-= 1
            obj_new.header.recording[:optode_pairs][:, 1] = chp
        elseif ol in detector_labels(obj_new)
            chp = obj_new.header.recording[:optode_pairs][:, 2]
            chs_to_delete = vcat(chs_to_delete, findall(isequal(idx), chp))
            deleteat!(obj_new.header.recording[:det_labels], obj_new.header.recording[:det_labels] .== ol)
            chp[chp .== idx] .= 0
            chp[chp .> idx] .-= 1
            obj_new.header.recording[:optode_pairs][:, 2] = chp
        end
    end

    reset_components!(obj_new)
    push!(obj_new.history, "delete_optode(OBJ, opt=$opt)")

    chs_to_delete = sort(unique(chs_to_delete))
    _info("Deleting the following NIRS channels: $(_v2s(chs_to_delete))")
    delete_channel!(obj_new, ch=chs_to_delete, del_opt=true)
    
    return obj_new

end

"""
    delete_optode!(obj; opt)

Delete optopode(s).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `opt::Union{Int64, Vector{Int64}, <:AbstractRange}`: optopode number(s) to be removed
"""
function delete_optode!(obj::NeuroAnalyzer.NEURO; opt::Union{Int64, Vector{Int64}, <:AbstractRange})

    obj_new = delete_optode(obj, opt=opt)
    obj.header = obj_new.header
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components
    obj.locs = obj_new.locs

    return nothing

end