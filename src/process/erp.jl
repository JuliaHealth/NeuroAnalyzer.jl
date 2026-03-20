export average_epochs
export average_epochs!
export sort_epochs
export sort_epochs!

"""
    average_epochs(obj; <keyword arguments>)

Average EEG/MEG epochs and prepend the average as epoch 1. Non-signal channels are removed. obj.header.recording[:data_type]` is updated to `"erp"` (EEG) or `"erf"` (MEG).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `bl::Tuple{Real, Real}=(0, 0)`: baseline window in seconds; when not `(0, 0)` the mean of this window is subtracted from the signal
- `blfirst::Bool=false`: when `true`, subtract the baseline from each epoch before averaging; when `false`, subtract the baseline from the averaged epoch only

# Returns

- `NeuroAnalyzer.NEURO`
"""
function average_epochs(
    obj::NeuroAnalyzer.NEURO;
    bl::Tuple{Real, Real} = (0, 0),
    blfirst::Bool = false
)::NeuroAnalyzer.NEURO

    _check_datatype(obj, ["eeg", "meg"])

    nchannels(obj) > length(get_channel(obj, type = datatype(obj))) &&
        _warn("Non-signal channels will be removed.")

    obj_new = if datatype(obj) == "eeg"
        keep_channel(obj, ch = get_channel(obj, type = datatype(obj)))
    else
        keep_channel(obj, ch = ["meg", "mag", "grad"])
    end

    # remove baseline prior to averaging
    bl_samples = if bl != (0, 0)
        _check_tuple(bl, (obj.epoch_time[1], obj.epoch_time[end]), "bl")
        (vsearch(bl[1], obj.epoch_time), vsearch(bl[2], obj.epoch_time))
    else
        nothing
    end

    # subtract baseline from each epoch before averaging (pre-average correction)
    if blfirst && !isnothing(bl_samples)
        obj_new.data[:, :, :] = remove_dc(obj_new.data[:, :, :], bl_samples)
    end

    # prepend the trial average as epoch 1; original epochs follow
    obj_new.data = cat(mean(obj_new.data; dims = 3), obj_new.data; dims = 3)

    obj_new.header.recording[:data_type] = datatype(obj) == "eeg" ? "erp" : "erf"
    obj_new.time_pts, obj_new.epoch_time = _get_t(obj_new)

    # subtract baseline from the averaged epoch only (post-average correction)
    if !blfirst && !isnothing(bl_samples)
        obj_new.data[:, :, 1] = remove_dc(obj_new.data[:, :, [1]], bl_samples)
    end

    # ------------------------------------------------------------------ #
    # update markers                                                     #
    # remove markers that fall outside the remaining data and shift the  #
    # remaining ones to account for the prepended average epoch.         #
    # ------------------------------------------------------------------ #
    ep_len = size(obj_new.data, 2)
    for idx in DataFrames.nrow(obj_new.markers):-1:1
        obj_new.markers[idx, :start] > ep_len && deleteat!(obj_new.markers, idx)
    end
    obj_new.markers[!, :start] .+= ep_len

    push!(obj_new.history, "average_epochs(OBJ, bl=$bl, blfirst=$blfirst)")

    return obj_new

end

"""
    average_epochs!(obj; <keyword arguments>)

Average EEG/MEG epochs in-place and prepend the average as epoch 1. Non-signal channels are removed. obj.header.recording[:data_type]` is updated to `"erp"` (EEG) or `"erf"` (MEG).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `bl::Tuple{Real, Real}=(0, 0)`: baseline window in seconds; when not `(0, 0)` the mean of this window is subtracted from the signal
- `blfirst::Bool=false`: when `true`, subtract the baseline from each epoch before averaging; when `false`, subtract the baseline from the averaged epoch only

# Returns

- `Nothing`
"""
function average_epochs!(
    obj::NeuroAnalyzer.NEURO;
    bl::Tuple{Real, Real} = (0, 0),
    blfirst::Bool = false
)::Nothing

    obj_new = average_epochs(obj, bl = bl, blfirst = blfirst)
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.header = obj_new.header
    obj.time_pts = obj_new.time_pts
    obj.epoch_time = obj_new.epoch_time
    obj.markers = obj_new.markers

    return nothing

end

"""
    sort_epochs(obj; <keyword arguments>)

Sort epochs 2:end of an ERP/ERF object according to a permutation vector. Epoch 1 (the average) is always kept in place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object (must be ERP or ERF)
- `s::Vector{Int64}`: permutation vector of length `nepochs(obj) - 1`; values must be in `2:nepochs(obj)` to avoid overwriting the average epoch

# Returns

- `NeuroAnalyzer.NEURO`
"""
function sort_epochs(obj::NeuroAnalyzer.NEURO; s::Vector{Int64})::NeuroAnalyzer.NEURO

    _check_datatype(obj, ["erp", "erf"])
    length(s) == nepochs(obj) - 1 ||
        throw(ArgumentError(
            "Length of s must be $(nepochs(obj) - 1) (number of non-average epochs)."))

    all(i -> 2 <= i <= nepochs(obj), s) ||
        throw(ArgumentError(
            "All values in s must be in 2:$(nepochs(obj)); epoch 1 is the average and cannot be reordered."))

    obj_new = deepcopy(obj)
    obj_new.data[:, :, 2:end] = obj.data[:, :, s]

    _warn("Markers are not sorted when epochs are reordered.")
    push!(obj_new.history, "sort_epochs(OBJ, s=$s)")

    return obj_new

end

"""
    sort_epochs(obj; <keyword arguments>)

Sort epochs 2:end in-place of an ERP/ERF object according to a permutation vector. Epoch 1 (the average) is always kept in place.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object (must be ERP or ERF)
- `s::Vector{Int64}`: permutation vector of length `nepochs(obj) - 1`; values must be in `2:nepochs(obj)` to avoid overwriting the average epoch

# Returns

- `Nothing`
"""
function sort_epochs!(obj::NeuroAnalyzer.NEURO; s::Vector{Int64})::Nothing

    obj_new = sort_epochs(obj, s = s)
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.markers = obj_new.markers

    return nothing

end
