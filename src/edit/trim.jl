export trim
export trim!

"""
    trim(s; seg, inverse)

Remove segment from the signal.

# Arguments

- `v::AbstractVector`
- `seg::Tuple{Int64, Int64}`: segment (from, to) in samples
- `inverse::Bool=false`: if true, keep the segment

# Returns

- `trim::Vector{Float64}`
"""
function trim(v::AbstractVector; seg::Tuple{Int64, Int64}, inverse::Bool=false)

    _check_segment(v, seg[1], seg[2])
    
    if inverse == false
        return vcat(v[1:seg[1] - 1], v[(seg[2] + 1):end])
    else
        return v[seg[1]:seg[2]]
    end

end

"""
    trim(m; seg, inverse)

Remove segment from the signal.

# Arguments

- `m::AbstractMatrix`
- `seg::Tuple{Int64, Int64}`: segment (from, to) in samples
- `inverse::Bool=false`: if true, keep the segment

# Returns

- `trim::Array{Float64}`
"""
function trim(m::AbstractMatrix; seg::Tuple{Int64, Int64}, inverse::Bool=false)
    
    _check_segment(m[1, :], seg[1], seg[2])

    if inverse == false
        return hcat(m[:, 1:(seg[1] - 1)], m[:, (seg[2] + 1):end])
    else
        return m[:, seg[1]:seg[2]]
    end
end

"""
    trim(a; seg, inverse)

Remove segment from the signal.

# Arguments

- `a::AbstractArray`
- `seg::Tuple{Int64, Int64}`: segment (from, to) in samples
- `inverse::Bool=false`: if true, keep the segment

# Returns

- `trim::Array{Float64}`
"""
function trim(a::AbstractArray; seg::Tuple{Int64, Int64}, inverse::Bool=false)

    _check_segment(a[1, :, 1], seg[1], seg[2])

    if inverse == false
        return hcat(a[:, 1:(seg[1] - 1), :], a[:, (seg[2] + 1):end, :])
    else
        return a[:, seg[1]:seg[2], :]
    end

end

"""
    trim(obj; seg, inverse, remove_epochs)

Trim signal by removing parts of the signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `seg::Tuple{Real, Real}`: segment to be removed (from, to) in seconds
- `inverse::Bool=false`: if true, keep the segment
- `remove_epochs::Bool=true`: if true, remove epochs containing signal to trim or remove signal and re-epoch trimmed signal

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function trim(obj::NeuroAnalyzer.NEURO; seg::Tuple{Real, Real}, inverse::Bool=false, remove_epochs::Bool=true)

    _check_segment(obj, seg)
    seg = (vsearch(seg[1], obj.time_pts), vsearch(seg[2], obj.time_pts))

    if remove_epochs == true
        epoch_n(obj) == 1 && throw(ArgumentError("OBJ has only one epoch, cannot use remove_epochs=true."))
        # seg = (vsearch(seg[1], obj.time_pts), vsearch(seg[2], obj.time_pts))
        eps = NeuroAnalyzer._s2epoch(obj, seg[1], seg[2])
        if inverse == false
            _info("Removing epochs: $eps.")
            obj_new = delete_epoch(obj, ep=eps)
        else
            _info("Keeping epochs: $eps.")
            obj_new = keep_epoch(obj, ep=eps)
        end
    else
        obj_new = deepcopy(obj)
        epoch_n(obj) > 1 && (epoch!(obj_new, ep_n=1))

        obj_new.data = trim(obj_new.data, seg=seg, inverse=inverse)
        obj_new.time_pts, obj_new.epoch_time = _get_t(obj_new)

        if epoch_n(obj) > 1
            if epoch_len(obj) <= signal_len(obj_new)
                epoch!(obj_new, ep_len=epoch_len(obj))
            else
                _info("Cannot apply original epoch length, returning single-epoch OBJ.")
            end
        end
        
        obj_new.markers = _delete_markers(obj_new.markers, seg)
        obj_new.markers = _shift_markers(obj_new.markers, seg[1], length(seg[1]:seg[2]))

    end

    reset_components!(obj_new)
    push!(obj_new.history, "trim(OBJ, seg=$seg, remove_epochs=$remove_epochs)")

    return obj_new
end

"""
    trim!(obj; seg, inverse, remove_epochs)

Trim signal by removing parts of the signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `seg::Tuple{Real, Real}`: segment to be removed (from, to) in seconds
- `inverse::Bool=false`: if true, keep the segment
- `remove_epochs::Bool=true`: if true, remove epochs containing signal to trim or remove signal and re-epoch trimmed signal
"""
function trim!(obj::NeuroAnalyzer.NEURO; seg::Tuple{Real, Real}, inverse::Bool=false, remove_epochs::Bool=true)

    obj_new = trim(obj, seg=seg, inverse=inverse, remove_epochs=remove_epochs)
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components
    obj.time_pts = obj_new.time_pts
    obj.epoch_time = obj_new.epoch_time

    return nothing

end
