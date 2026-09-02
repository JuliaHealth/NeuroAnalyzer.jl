export trim
export trim!
export crop
export crop!

"""
    trim(s; <keyword arguments>)

Remove segment from the signal.

# Arguments

- `s::AbstractVector`: signal vector
- `seg::Tuple{Int64, Int64}`: segment (from, to) in samples
- `keep::Bool=false`: if `true`, keep the segment

# Returns

- `Vector{Float64}`
"""
function trim(
    s::AbstractVector;
    seg::Tuple{Int64, Int64},
    keep::Bool = false,
)::Vector{Float64}
    # validate
    _check_segment(s, seg[1], seg[2])

    if keep
        return s[seg[1]:seg[2]]
    else
        return vcat(s[1:(seg[1] - 1)], s[(seg[2] + 1):end])
    end
end

"""
    trim(s; <keyword arguments>)

Remove segment from the signal.

# Arguments

- `s::AbstractMatrix`: signal matrix, shape (channel, samples)
- `seg::Tuple{Int64, Int64}`: segment (from, to) in samples
- `keep::Bool=false`: if `true`, keep the segment

# Returns

- `Matrix{Float64}`
"""
function trim(
    s::AbstractMatrix;
    seg::Tuple{Int64, Int64},
    keep::Bool = false,
)::Matrix{Float64}
    # validate
    _check_segment(s[1, :], seg[1], seg[2])

    if keep
        return s[:, seg[1]:seg[2]]
    else
        return hcat(s[:, 1:(seg[1] - 1)], s[:, (seg[2] + 1):end])
    end
end

"""
    trim(s; <keyword arguments>)

Remove segment from a 3-D signal array.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `seg::Tuple{Int64, Int64}`: segment (from, to) in samples
- `keep::Bool=false`: if `true`, keep the segment

# Returns

- `Array{Float64, 3}`
"""
function trim(
    s::AbstractArray;
    seg::Tuple{Int64, Int64},
    keep::Bool = false,
)::Array{Float64, 3}
    # validate that the input is s proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # validate
    _check_segment(s[1, :, 1], seg[1], seg[2])

    if keep
        return s[:, seg[1]:seg[2], :]
    else
        return hcat(s[:, 1:(seg[1] - 1), :], s[:, (seg[2] + 1):end, :])
    end
end

"""
    trim(obj; <keyword arguments>)

Trim signal by removing parts of the signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `seg::Tuple{Real, Real}`: segment to be removed (from, to) in seconds
- `keep::Bool=false`: if `true`, keep the segment
- `remove_epochs::Bool=false`: if `true`, remove epochs containing signal to trim or remove signal and re-epoch trimmed signal

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function trim(
    obj::NeuroAnalyzer.NEURO;
    seg::Tuple{Real, Real},
    keep::Bool = false,
)::NeuroAnalyzer.NEURO
    # validate
    nepochs(obj) == 1 ||
        throw(ArgumentError("trim() must be applied to continuous object."))
    NeuroAnalyzer._check_segment(obj, seg)

    s_idx = vsearch(seg[1], obj.time_pts)
    seg_tpos = (vsearch(seg[1], obj.time_pts), vsearch(seg[2], obj.time_pts))

    (datatype(obj) == "meg" && size(obj.header.recording[:ssp_data], 1) != 0) && _warn(
        "OBJ contains SSP projections data, you should apply them before modifying OBJ data.",
    )

    # create new dataset
    obj_tmp = deepcopy(obj)

    obj_tmp.data = trim(obj_tmp.data; seg = seg_tpos, keep = keep)

    if keep
        obj_tmp.time_pts = obj.time_pts[seg_tpos[1]:seg_tpos[2]]
        obj_tmp.epoch_time = obj.time_pts[seg_tpos[1]:seg_tpos[2]]
    else
        obj_tmp.time_pts, obj_tmp.epoch_time = _get_t(obj_tmp)
    end

    obj_tmp.markers = _delete_markers(obj_tmp.markers, seg)
    obj_tmp.markers = _shift_markers(obj_tmp.markers, seg)

    if keep
        obj_tmp.time_pts = obj.time_pts[1:size(obj_tmp.data, 2)]
        obj_tmp.epoch_time = obj.time_pts[1:size(obj_tmp.data, 2)]
    end

    if !keep
        if s_idx <= length(obj_tmp.time_pts)
            add_marker!(
                obj_tmp; id = "NA", start = obj_tmp.time_pts[s_idx], len = 0.0,
                value = "DELETED",
            )
            obj_tmp.markers = unique(obj_tmp.markers)
        else
            # if the terminal part is removed the marker is placed on the time point
            add_marker!(
                obj_tmp; id = "NA", start = obj_tmp.time_pts[s_idx - 1], len = 0.0,
                value = "DELETED",
            )
            obj_tmp.markers = unique(obj_tmp.markers)
        end
    end

    push!(obj_tmp.history, "trim(obj, seg=$seg, keep=$keep")

    return obj_tmp
end

"""
    trim!(obj; <keyword arguments>)

Trim signal by removing parts of the signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `seg::Tuple{Real, Real}`: segment to be removed (from, to) in seconds
- `keep::Bool=false`: if `true`, keep the segment

# Returns

- `Nothing`
"""
function trim!(
    obj::NeuroAnalyzer.NEURO;
    seg::Tuple{Real, Real},
    keep::Bool = false,
)::Nothing
    # validate
    nepochs(obj) == 1 ||
        throw(ArgumentError("trim!() must be applied to continuous object."))

    obj_tmp = trim(obj; seg = seg, keep = keep)
    obj.data = obj_tmp.data
    obj.history = obj_tmp.history
    obj.time_pts = obj_tmp.time_pts
    obj.epoch_time = obj_tmp.epoch_time
    obj.markers = obj_tmp.markers
    obj_tmp = nothing
    obj_tmp = nothing

    return nothing
end

"""
    crop(obj; <keyword arguments>)

Crop signal by removing parts of the signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `seg::Tuple{Real, Real}`: segment to be cropped (from, to) in seconds

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function crop(obj::NeuroAnalyzer.NEURO; seg::Tuple{Real, Real})::NeuroAnalyzer.NEURO
    # validate
    nepochs(obj) == 1 ||
        throw(ArgumentError("crop() must be applied to continuous object."))

    obj_tmp = trim(obj; seg = seg, keep = true)

    return obj_tmp
end

"""
    crop!(obj; <keyword arguments>)

Crop signal by removing parts of the signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `seg::Tuple{Real, Real}`: segment to be cropped (from, to) in seconds

# Returns

- `Nothing`
"""
function crop!(obj::NeuroAnalyzer.NEURO; seg::Tuple{Real, Real})::Nothing
    nepochs(obj) == 1 ||
        throw(ArgumentError("crop!() must be applied to continuous object."))

    obj_tmp = trim(obj; seg = seg, keep = true)
    obj.data = obj_tmp.data
    obj.history = obj_tmp.history
    obj.time_pts = obj_tmp.time_pts
    obj.epoch_time = obj_tmp.epoch_time
    obj.markers = obj_tmp.markers
    obj_tmp = nothing
    obj_tmp = nothing

    return nothing
end
