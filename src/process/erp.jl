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
    blfirst::Bool = false,
)::NeuroAnalyzer.NEURO
    # validate
    _check_datatype(obj, ["eeg", "meg"])

    nchannels(obj) > length(get_channel(obj; type = datatype(obj))) &&
        _warn("Non-signal channels will be removed.")

    # signal channels (get_channel(; type) returns names)
    obj_tmp = if datatype(obj) == "eeg"
        keep_channel(obj; ch = get_channel(obj; type = "eeg"))
    else
        keep_channel(obj; ch = get_channel(obj; type = ["meg", "mag", "grad"]))
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
        obj_tmp.data[:, :, :] = _remove_dc(obj_tmp.data[:, :, :], bl_samples)
    end

    # prepend the trial average as epoch 1; original epochs follow
    obj_tmp.data = cat(
        mean(
            obj_tmp.data; dims = 3,
        ), obj_tmp.data; dims = 3,
    )

    obj_tmp.header.recording[:data_type] = datatype(obj) == "eeg" ? "erp" : "erf"
    obj_tmp.time_pts, obj_tmp.epoch_time = _get_t(obj_tmp)

    # subtract baseline from the averaged epoch only (post-average correction)
    if !blfirst && !isnothing(bl_samples)
        obj_tmp.data[:, :, 1] = _remove_dc(obj_tmp.data[:, :, [1]], bl_samples)
    end

    # markers (in seconds) follow their epochs, which are shifted by the prepended average epoch
    obj_tmp.markers[!, :start] .+= n2t(size(obj_tmp.data, 2), sr(obj))

    push!(obj_tmp.history, "average_epochs(obj, bl=$bl, blfirst=$blfirst)")

    return obj_tmp
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
    blfirst::Bool = false,
)::Nothing
    obj_tmp = average_epochs(obj; bl = bl, blfirst = blfirst)
    obj.data = obj_tmp.data
    obj.history = obj_tmp.history
    obj.header = obj_tmp.header
    obj.time_pts = obj_tmp.time_pts
    obj.epoch_time = obj_tmp.epoch_time
    obj.markers = obj_tmp.markers
    obj.locs = obj_tmp.locs
    obj_tmp = nothing

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
    # validate
    _check_datatype(obj, ["erp", "erf"])
    length(s) == nepochs(obj) - 1 ||
        throw(
            ArgumentError(
                "Length of s must be $(nepochs(obj) - 1) (number of non-average epochs).",
            ),
        )

    all(i -> 2 <= i <= nepochs(obj), s) ||
        throw(
            ArgumentError(
                "All values in s must be in 2:$(nepochs(obj)); epoch 1 is the average and cannot be reordered.",
            ),
        )
    allunique(s) || throw(ArgumentError("s must be a permutation (no repeated epochs)."))

    # create new dataset
    obj_tmp = deepcopy(obj)

    obj_tmp.data[:, :, 2:end] = obj.data[:, :, s]

    # markers move with their epochs: epoch s[i] becomes epoch i + 1
    ep_len = epoch_len(obj)
    new_pos = Dict(e => i + 1 for (i, e) in enumerate(s))
    obj_tmp.markers = _remap_markers(
        obj.markers,
        sr(obj),
        q -> begin
            e = (q - 1) ÷ ep_len + 1
            e == 1 ? q : (new_pos[e] - 1) * ep_len + (q - 1) % ep_len + 1
        end,
    )
    push!(obj_tmp.history, "sort_epochs(obj, s=$s)")

    return obj_tmp
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
    obj_tmp = sort_epochs(obj; s = s)
    obj.data = obj_tmp.data
    obj.history = obj_tmp.history
    obj.markers = obj_tmp.markers
    obj_tmp = nothing

    return nothing
end
