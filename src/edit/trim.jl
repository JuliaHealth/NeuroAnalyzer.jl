export trim
export trim!

"""
    trim(s; <keyword arguments>)

Remove segment from the signal.

# Arguments

- `v::AbstractVector`
- `seg::Tuple{Int64, Int64}`: segment (from, to) in samples
- `inverse::Bool=false`: if true, keep the segment

# Returns

- `trim::Vector{Float64}`
"""
function trim(v::AbstractVector; seg::Tuple{Int64, Int64}, inverse::Bool=false)::Vector{Float64}

    _check_segment(v, seg[1], seg[2])

    if inverse
        return v[seg[1]:seg[2]]
    else
        return vcat(v[1:seg[1] - 1], v[(seg[2] + 1):end])
    end

end

"""
    trim(m; <keyword arguments>)

Remove segment from the signal.

# Arguments

- `m::AbstractMatrix`
- `seg::Tuple{Int64, Int64}`: segment (from, to) in samples
- `inverse::Bool=false`: if true, keep the segment

# Returns

- `trim::Matrix{Float64}`
"""
function trim(m::AbstractMatrix; seg::Tuple{Int64, Int64}, inverse::Bool=false)::Matrix{Float64}

    _check_segment(m[1, :], seg[1], seg[2])

    if inverse
        return m[:, seg[1]:seg[2]]
    else
        return hcat(m[:, 1:(seg[1] - 1)], m[:, (seg[2] + 1):end])
    end

end

"""
    trim(a; <keyword arguments>)

Remove segment from the signal.

# Arguments

- `a::AbstractArray`
- `seg::Tuple{Int64, Int64}`: segment (from, to) in samples
- `inverse::Bool=false`: if true, keep the segment

# Returns

- `trim::Array{Float64, 3}`
"""
function trim(a::AbstractArray; seg::Tuple{Int64, Int64}, inverse::Bool=false)::Array{Float64, 3}

    _chk3d(a)
    _check_segment(a[1, :, 1], seg[1], seg[2])

    if inverse
        return a[:, seg[1]:seg[2], :]
    else
        return hcat(a[:, 1:(seg[1] - 1), :], a[:, (seg[2] + 1):end, :])
    end

end

"""
    trim(obj; <keyword arguments>)

Trim signal by removing parts of the signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `seg::Tuple{Real, Real}`: segment to be removed (from, to) in seconds
- `inverse::Bool=false`: if true, keep the segment
- `remove_epochs::Bool=true`: if true, remove epochs containing signal to trim or remove signal and re-epoch trimmed signal

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function trim(obj::NeuroAnalyzer.NEURO; seg::Tuple{Real, Real}, inverse::Bool=false, remove_epochs::Bool=true)::NeuroAnalyzer.NEURO

    _check_segment(obj, seg)

    seg = (vsearch(seg[1], obj.time_pts), vsearch(seg[2], obj.time_pts))

    (datatype(obj) == "meg" && size(obj.header.recording[:ssp_data]) != (0,)) && _warn("OBJ contains SSP projections data, you should apply them before modifying OBJ data.")

    if remove_epochs
        @assert nepochs(obj) > 1 "OBJ has only one epoch, cannot use remove_epochs=true."
        # seg = (vsearch(seg[1], obj.time_pts), vsearch(seg[2], obj.time_pts))
        eps = _s2epoch(obj, seg[1], seg[2])
        if !inverse
            _info("Removing epochs: $eps")
            obj_new = delete_epoch(obj, ep=eps)
        else
            _info("Keeping epochs: $eps")
            obj_new = keep_epoch(obj, ep=eps)
        end
    else
        obj_new = deepcopy(obj)
        nepochs(obj) > 1 && (epoch!(obj_new, ep_n=1))

        obj_new.data = trim(obj_new.data, seg=seg, inverse=inverse)
        obj_new.time_pts, obj_new.epoch_time = _get_t(obj_new)

        if nepochs(obj) > 1
            if epoch_len(obj) <= signal_len(obj_new)
                epoch!(obj_new, ep_len=epoch_len(obj) / sr(obj))
            else
                _warn("Cannot apply original epoch length, returning single-epoch OBJ.")
            end
        end

        obj_new.markers = _delete_markers(obj_new.markers, seg, sr(obj))
        obj_new.markers = _shift_markers(obj_new.markers, seg[1], length(seg[1]:seg[2]), sr(obj))

    end

    reset_components!(obj_new)
    push!(obj_new.history, "trim(OBJ, seg=$seg, remove_epochs=$remove_epochs)")

    return obj_new

end

"""
    trim!(obj; <keyword arguments>)

Trim signal by removing parts of the signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `seg::Tuple{Real, Real}`: segment to be removed (from, to) in seconds
- `inverse::Bool=false`: if true, keep the segment
- `remove_epochs::Bool=true`: if true, remove epochs containing signal to trim or remove signal and re-epoch trimmed signal

# Returns

Nothing
"""
function trim!(obj::NeuroAnalyzer.NEURO; seg::Tuple{Real, Real}, inverse::Bool=false, remove_epochs::Bool=true)::Nothing

    obj_new = trim(obj, seg=seg, inverse=inverse, remove_epochs=remove_epochs)
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components
    obj.time_pts = obj_new.time_pts
    obj.epoch_time = obj_new.epoch_time
    obj.markers = obj_new.markers

    return nothing

end
