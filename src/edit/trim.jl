export trim
export trim!

"""
    trim(s; <keyword arguments>)

Remove segment from the signal.

# Arguments

- `v::AbstractVector`
- `seg::Tuple{Int64, Int64}`: segment (from, to) in samples
- `keep::Bool=false`: if true, keep the segment

# Returns

- `trim::Vector{Float64}`
"""
function trim(v::AbstractVector; seg::Tuple{Int64, Int64}, keep::Bool=false)::Vector{Float64}

    _check_segment(v, seg[1], seg[2])

    if keep
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
- `keep::Bool=false`: if true, keep the segment

# Returns

- `trim::Matrix{Float64}`
"""
function trim(m::AbstractMatrix; seg::Tuple{Int64, Int64}, keep::Bool=false)::Matrix{Float64}

    _check_segment(m[1, :], seg[1], seg[2])

    if keep
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
- `keep::Bool=false`: if true, keep the segment

# Returns

- `trim::Array{Float64, 3}`
"""
function trim(a::AbstractArray; seg::Tuple{Int64, Int64}, keep::Bool=false)::Array{Float64, 3}

    _chk3d(a)
    _check_segment(a[1, :, 1], seg[1], seg[2])

    if keep
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
- `keep::Bool=false`: if true, keep the segment
- `remove_epochs::Bool=false`: if true, remove epochs containing signal to trim or remove signal and re-epoch trimmed signal

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function trim(obj::NeuroAnalyzer.NEURO; seg::Tuple{Real, Real}, keep::Bool=false)::NeuroAnalyzer.NEURO

    @assert nepochs(obj) == 1 "trim() must be applied to a continuous object."
    NeuroAnalyzer._check_segment(obj, seg)

    s_idx = findfirst(x -> x == seg[1], obj.time_pts) - 1
    seg = (vsearch(seg[1], obj.time_pts), vsearch(seg[2], obj.time_pts))

    (datatype(obj) == "meg" && size(obj.header.recording[:ssp_data], 1) != 0) && _warn("OBJ contains SSP projections data, you should apply them before modifying OBJ data.")

    obj_new = deepcopy(obj)
    obj_new.data = trim(obj_new.data, seg=seg, keep=keep)
    obj_new.time_pts, obj_new.epoch_time = _get_t(obj_new)

    obj_new.markers = _delete_markers(obj_new.markers, seg, sr(obj))
    obj_new.markers = _shift_markers(obj_new.markers, seg[1], length(seg[1]:seg[2]), sr(obj))

    add_marker!(obj_new, id="NA", start=obj_new.time_pts[s_idx - 1], value="DELETED")

    reset_components!(obj_new)
    push!(obj_new.history, "trim(OBJ, seg=$seg, keep=$keep")

    return obj_new

end

"""
    trim!(obj; <keyword arguments>)

Trim signal by removing parts of the signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `seg::Tuple{Real, Real}`: segment to be removed (from, to) in seconds
- `keep::Bool=false`: if true, keep the segment

# Returns

Nothing
"""
function trim!(obj::NeuroAnalyzer.NEURO; seg::Tuple{Real, Real}, keep::Bool=false)::Nothing

    @assert nepochs(obj) == 1 "trim!() must be applied to a continuous object."

    obj_new = trim(obj, seg=seg, keep=keep)
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components
    obj.time_pts = obj_new.time_pts
    obj.epoch_time = obj_new.epoch_time
    obj.markers = obj_new.markers

    return nothing

end
