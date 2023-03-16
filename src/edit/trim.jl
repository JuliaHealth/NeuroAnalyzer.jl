export trim
export trim!

"""
    trim(s; seg)

Remove segment from the signal.

# Arguments

- `v::AbstractVector`
- `seg::Tuple{Int64, Int64}`: segment (from, to) in samples

# Returns

- `trim::Vector{Float64}`
"""
function trim(v::AbstractVector; seg::Tuple{Int64, Int64})

    _check_segment(v, seg[1], seg[2])
    
    return vcat(v[1:seg[1] - 1], v[(seg[2] + 1):end])

end

"""
    trim(m; seg)

Remove segment from the signal.

# Arguments

- `m::AbstractMatrix`
- `seg::Tuple{Int64, Int64}`: segment (from, to) in samples

# Returns

- `trim::Array{Float64}`
"""
function trim(m::AbstractMatrix; seg::Tuple{Int64, Int64})
    
    _check_segment(m[1, :], seg[1], seg[2])

    return hcat(m[:, 1:(seg[1] - 1)], m[:, (seg[2] + 1):end])

end

"""
    trim(a; seg)

Remove segment from the signal.

# Arguments

- `a::AbstractArray`
- `seg::Tuple{Int64, Int64}`: segment (from, to) in samples

# Returns

- `trim::Array{Float64}`
"""
function trim(a::AbstractArray; seg::Tuple{Int64, Int64})

    _check_segment(a[1, :, 1], seg[1], seg[2])

    return hcat(a[:, 1:(seg[1] - 1), :], a[:, (seg[2] + 1):end, :])

end

"""
    trim(obj; seg, remove_epochs)

Trim signal by removing parts of the signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `seg::Tuple{Int64, Int64}`: segment to be removed (from, to) in samples
- `remove_epochs::Bool=true`: if true, remove epochs containing signal to trim or remove signal and re-epoch trimmed signal

# Returns

- OBJ
"""
function trim(obj::NeuroAnalyzer.NEURO; seg::Tuple{Int64, Int64}, remove_epochs::Bool=true)

    if remove_epochs == true
        epoch_n(obj) == 1 && throw(ArgumentError("OBJ has only one epoch, cannot use remove_epochs=true."))
        eps = _s2epoch(obj, seg[1], seg[2])
        _info("Removing epochs: $eps.")
        obj_new = delete_epoch(obj, ep=eps)
    else
        obj_new = deepcopy(obj)
        epoch_n(obj) > 1 && (epoch!(obj_new, ep_n=1))
        _check_segment(obj_new, seg[1], seg[2])
        obj_new.data = trim(obj_new.data, seg=seg)
        obj_new.time_pts = round.(collect(obj_new.time_pts[1]:(1 / sr(obj_new)):((size(obj_new.data, 2) * size(obj_new.data, 3)) / sr(obj_new))) .- (1 / sr(obj_new)), digits=3)
        obj_new.epoch_time = round.(collect(obj_new.time_pts[1]:(1 / sr(obj_new)):(size(obj_new.data, 2) / sr(obj_new))) .- (1 / sr(obj_new)), digits=3)
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
    trim!(obj; seg, remove_epochs)

Trim signal by removing parts of the signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `seg::Tuple{Int64, Int64}`: segment to be removed (from, to) in samples
- `remove_epochs::Bool=true`: remove epochs containing signal to trim (remove_epochs=true) or remove signal and remove epoching
"""
function trim!(obj::NeuroAnalyzer.NEURO; seg::Tuple{Int64, Int64}, remove_epochs::Bool=true)

    obj_new = trim(obj, seg=seg, remove_epochs=remove_epochs)
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components
    obj.time_pts = obj_new.time_pts
    obj.epoch_time = obj_new.epoch_time

    return nothing

end
