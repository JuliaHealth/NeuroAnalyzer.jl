export average_epochs
export average_epochs!
export sort_epochs
export sort_epochs!

"""
    average_epochs(obj; <keyword arguments>)

Average EEG/MEG/NIRS epochs and prepend the average as epoch 1. Non-signal channels are removed. `obj.header.recording[:data_type]` is updated to `"erp"` (EEG) or `"erf"` (MEG); for NIRS it stays `"nirs"` (epoch 1 is the averaged hemodynamic response, the auxiliary channels are removed).

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
    obj_new = deepcopy(obj)
    average_epochs!(obj_new; bl = bl, blfirst = blfirst)

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
    blfirst::Bool = false,
)::Nothing
    # validate
    _check_datatype(obj, ["eeg", "meg", "nirs"])

    # signal channels (get_channel(; type) returns names)
    sig_ch = if datatype(obj) == "eeg"
        get_channel(obj; type = "eeg")
    elseif datatype(obj) == "meg"
        get_channel(obj; type = ["meg", "mag", "grad"])
    else
        # NIRS: all optical and derived channels (auxiliary channels are removed)
        labels(obj)[obj.header.recording[:channel_type] .!= "nirs_aux"]
    end
    isempty(sig_ch) && throw(ArgumentError("OBJ does not contain signal channels."))
    nchannels(obj) > length(sig_ch) && _warn("Non-signal channels will be removed.")
    # data type before it is changed (EEG → ERP, MEG → ERF)
    dt = datatype(obj)
    # non-signal channels are removed
    nchannels(obj) > length(sig_ch) && keep_channel!(obj; ch = sig_ch)

    # remove baseline prior to averaging
    bl_samples = if bl != (0, 0)
        _check_tuple(bl, (obj.epoch_time[1], obj.epoch_time[end]), "bl")
        (vsearch(bl[1], obj.epoch_time), vsearch(bl[2], obj.epoch_time))
    else
        nothing
    end

    # subtract baseline from each epoch before averaging (pre-average correction)
    if blfirst && !isnothing(bl_samples)
        obj.data[:, :, :] = _remove_dc(obj.data[:, :, :], bl_samples)
    end

    # epoch metadata: the average epoch has the condition and source of the averaged
    # epochs (if they are the same), no onset (NaN)
    ep_tbl = _epochs_table(obj)
    _uq(x) = length(unique(x)) == 1 ? x[1] : join(unique(x), ";")
    obj.epochs = vcat(
        DataFrame(
            :condition => _uq(ep_tbl.condition),
            :source    => _uq(ep_tbl.source),
            :reject    => false,
            :onset     => NaN,
        ),
        ep_tbl,
    )

    # prepend the trial average as epoch 1; original epochs follow
    obj.data = cat(
        mean(
            obj.data; dims = 3,
        ), obj.data; dims = 3,
    )

    # EEG → ERP, MEG → ERF; NIRS stays NIRS (epoch 1 is the averaged hemodynamic response)
    dt == "eeg" && (obj.header.recording[:data_type] = "erp")
    dt == "meg" && (obj.header.recording[:data_type] = "erf")
    obj.time_pts, obj.epoch_time = _get_t(obj)

    # subtract baseline from the averaged epoch only (post-average correction)
    if !blfirst && !isnothing(bl_samples)
        obj.data[:, :, 1] = _remove_dc(obj.data[:, :, [1]], bl_samples)
    end

    # markers (in seconds) follow their epochs, which are shifted by the prepended average epoch
    obj.markers[!, :start] .+= n2t(size(obj.data, 2), sr(obj))

    push!(obj.history, "average_epochs(obj, bl=$bl, blfirst=$blfirst)")

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
    obj_new = deepcopy(obj)
    sort_epochs!(obj_new; s = s)

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

    obj.data[:, :, 2:end] = obj.data[:, :, s]
    obj.epochs = _epochs_table(obj)[vcat(1, s), :]

    # markers move with their epochs: epoch s[i] becomes epoch i + 1
    ep_len = epoch_len(obj)
    new_pos = Dict(e => i + 1 for (i, e) in enumerate(s))
    obj.markers = _remap_markers(
        obj.markers,
        sr(obj),
        q -> begin
            e = (q - 1) ÷ ep_len + 1
            e == 1 ? q : (new_pos[e] - 1) * ep_len + (q - 1) % ep_len + 1
        end,
    )
    push!(obj.history, "sort_epochs(obj, s=$s)")

    return nothing
end
