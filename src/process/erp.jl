export erp
export erp!
export sort_epochs
export sort_epochs!

"""
    erp(obj; bl, blfirst)

Average epochs. Non-signal channels are removed. `OBJ.header.recording[:data_type]` becomes `erp`. First epoch is the ERP.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `bl::Tuple{Real, Real}=(0, 0)`: baseline is from `b[1]` to `b[2]` seconds; if `bl` is greater than (0, 0), DC value is calculated as mean of the `b[1]` to `b[2]` seconds and subtracted from the signal
- `blfirst::Bool=false`: if true, subtract the baseline segment prior to averaging

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function erp(obj::NeuroAnalyzer.NEURO; bl::Tuple{Real, Real}=(0, 0), blfirst::Bool=false)

    nchannels(obj) > length(signal_channels(obj)) && _warn("Non-signal channels will be removed.")

    obj_new = keep_channel(obj, ch=signal_channels(obj))

    # remove DC
    if blfirst
        if bl != (0, 0)
            _check_tuple(bl, "bl", (obj.epoch_time[1], obj.epoch_time[end]))
            bl = (vsearch(bl[1], obj.epoch_time), vsearch(bl[2], obj.epoch_time))
            obj_new.data[:, :, :] = remove_dc(obj_new.data[:, :, :], bl)
        end
    end

    obj_new.data = cat(mean(obj_new.data, dims=3)[:, :, :], obj_new.data, dims=3)
    obj_new.header.recording[:data_type] = "erp"
    obj_new.time_pts, obj_new.epoch_time = _get_t(obj_new)

    # remove DC
    if !blfirst
        if bl != (0, 0)
            _check_tuple(bl, "bl", (obj.epoch_time[1], obj.epoch_time[end]))
            bl = (vsearch(bl[1], obj.epoch_time), vsearch(bl[2], obj.epoch_time))
            obj_new.data[:, :, 1] = remove_dc(obj_new.data[:, :, 1], bl)
        end
    end

    # remove markers of deleted epochs
    for marker_idx in nrow(obj_new.markers):-1:1
        obj_new.markers[marker_idx, :start] > size(obj_new.data, 2) && deleteat!(obj_new.markers, marker_idx)
    end
    obj_new.markers[!, :start] .+= (obj_new.epoch_time[1] * sr(obj_new))

    reset_components!(obj_new)
    push!(obj_new.history, "erp(OBJ, bl=$bl)")

    return obj_new

end

"""
    erp!(obj; bl, blfirst)

Average epochs. Non-signal channels are removed. `OBJ.header.recording[:data_type]` becomes `erp`. First epoch is the ERP.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `bl::Tuple{Real, Real}=(0, 0)`: baseline is the first `bl` seconds; if `bl` is greater than 0, DC value is calculated as mean of the first `n` samples and subtracted from the signal
- `blfirst::Bool=false`: if true, subtract the baseline segment prior to averaging
"""
function erp!(obj::NeuroAnalyzer.NEURO; bl::Tuple{Real, Real}=(0, 0), blfirst::Bool=false)

    obj_new = erp(obj, bl=bl, blfirst=blfirst)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history
    obj.header = obj_new.header
    obj.time_pts = obj_new.time_pts
    obj.epoch_time = obj_new.epoch_time
    obj.markers = obj_new.markers

    return nothing

end

"""
    sort_epochs(obj; s)

Sort epochs.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `s::Vector{Int64}`: vector of sorting indices

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function sort_epochs(obj::NeuroAnalyzer.NEURO; s::Vector{Int64})

    _check_datatype(obj, "erp")
    @assert length(s) == nepochs(obj) - 1 "Length of the sorting vector must be equal to the number of epochs."

    obj_new = deepcopy(obj)
    obj_new.data[:, :, 2:end] = obj.data[:, :, s]

    _warn("Markers are not sorted.")

    reset_components!(obj_new)
    push!(obj_new.history, "sort_epochs(OBJ, s=$s)")

    return obj_new

end

"""
    sort_epochs(obj; s)

Sort epochs.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `s::Vector{Int64}`: vector of sorting indices

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function sort_epochs!(obj::NeuroAnalyzer.NEURO; s::Vector{Int64})

    obj_new = sort_epochs(obj, s=s)
    obj.data = obj_new.data
    obj.components = obj_new.components
    obj.history = obj_new.history

    # to do: markers should be sorted
    _warn("Markers are not sorted.")
    obj.markers = obj_new.markers

    reset_components!(obj_new)
    push!(obj_new.history, "sort_epochs(OBJ, s=$s)")

    return obj_new

end