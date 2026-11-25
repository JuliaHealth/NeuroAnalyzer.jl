"""
    _epochs_default(ep_n, ep_dur; <keyword arguments>)

Default epoch metadata table: `ep_n` epochs of `ep_dur` seconds, without conditions, not rejected, onsets of consecutive epochs starting at `onset`.

# Returns

- `DataFrame`: columns `condition::String`, `source::String`, `reject::Bool`, `onset::Float64` (seconds)
"""
function _epochs_default(
    ep_n::Int64,
    ep_dur::Real;
    onset::Real = 0.0,
    source::String = "",
)::DataFrame
    return DataFrame(
        :condition => fill("", ep_n),
        :source    => fill(source, ep_n),
        :reject    => falses(ep_n),
        :onset     => Float64[onset + (idx - 1) * ep_dur for idx = 1:ep_n],
    )
end

# source recording of the epochs: the name of the imported file (without the folder)
_epochs_source(hdr::NeuroAnalyzer.HEADER)::String =
    basename(string(get(hdr.recording, :file_name, "")))

"""
    _epochs_upgrade!(obj)

Convert an object saved before v0.26.12: the `:epoch_id` header field is removed; the ID of marker-based epochs (the marker value) becomes the epoch condition.

# Returns

- `Nothing`
"""
function _epochs_upgrade!(obj::NeuroAnalyzer.NEURO)::Nothing
    haskey(obj.header.recording, :epoch_id) || return nothing
    ep_id = string(obj.header.recording[:epoch_id])
    delete!(obj.header.recording, :epoch_id)
    # fixed-length IDs: "length_<s>s", "full"
    (ep_id == "" || ep_id == "full" || startswith(ep_id, "length_")) && return nothing
    obj.epochs = copy(_epochs_table(obj))
    obj.epochs[!, :condition] .= ep_id
    return nothing
end

"""
    _epochs_table(obj)

Epoch metadata of `obj`; if the table does not match the number of epochs (e.g. after the data were replaced directly), the default table is returned.

# Returns

- `DataFrame`
"""
function _epochs_table(obj::NeuroAnalyzer.NEURO)::DataFrame
    DataFrames.nrow(obj.epochs) == nepochs(obj) && return obj.epochs
    return _epochs_default(
        nepochs(obj), epoch_duration(obj);
        source = _epochs_source(obj.header),
    )
end

function _make_epochs(
    s::AbstractMatrix;
    ep_len::Int64,
)::Array{Float64, 3}
    # validate
    ep_len >= 1 || throw(ArgumentError("ep_len must be ≥ 1."))
    ep_len <= size(s, 2) || throw(ArgumentError("ep_len must be ≤ $(size(s, 2))."))

    ch_n = size(s, 1)
    ep_n = size(s, 2) ÷ ep_len

    # trim trailing samples that don't fill a complete epoch, then reshape
    return reshape(s[:, 1:(ep_len * ep_n)], ch_n, ep_len, ep_n)
end

function _make_epochs(
    s::AbstractArray;
    ep_len::Int64,
)::Array{Float64, 3}
    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # validate
    ep_len >= 1 || throw(ArgumentError("ep_len must be ≥ 1."))
    ep_len <= size(s, 2) || throw(ArgumentError("ep_len must be ≤ $(size(s, 2))."))

    n_samp = size(s, 2) * size(s, 3)
    ep_n = n_samp ÷ ep_len

    # flatten epochs dimension, trim, then re-epoch
    s_flat = reshape(s, ch_n, n_samp)
    return reshape(s_flat[:, 1:(ep_len * ep_n)], ch_n, ep_len, ep_n)
end

function _make_epochs_bymarkers(
    s::AbstractArray;
    marker::String,
    markers::DataFrame,
    marker_start::Vector{Int64},
    offset::Int64,
    ep_len::Int64,
    fs::Int64,
)::Tuple{Array{Float64, 3}, DataFrame}
    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # validate
    # offset = 0 is valid: the epoch starts at the marker
    offset >= 0 || throw(ArgumentError("offset must be ≥ 0."))
    ep_len >= 1 || throw(ArgumentError("ep_len must be ≥ 1."))
    fs >= 1 || throw(ArgumentError("fs must be ≥ 1."))

    sig_len = size(s, 2)
    mrk_n = length(marker_start)

    # compute epoch sample boundaries
    ep_start = marker_start .- offset
    ep_end = ep_start .+ ep_len .- 1

    # remove epochs that fall (even partially) outside the signal
    # iterate backwards to allow safe deleteat!
    for idx = mrk_n:-1:1
        if ep_start[idx] < 1 || ep_end[idx] > sig_len
            deleteat!(ep_start, idx)
            deleteat!(ep_end, idx)
        end
    end

    mrk_n = length(ep_start)
    epochs = zeros(size(s, 1), ep_len, mrk_n)

    @inbounds for mrk_idx = 1:mrk_n
        # flatten the epoch dimension of s before slicing (s is 3-D with 1 epoch)
        epochs[:, :, mrk_idx] = reshape(
            s[:, ep_start[mrk_idx]:ep_end[mrk_idx], :],
            size(s, 1), ep_len,
        )
    end

    # markers (of any type) inside each extracted epoch, mapped to the epoched timeline:
    # sample m of epoch k (window ep_start[k]:ep_end[k]) becomes (k - 1) * ep_len + m - ep_start[k] + 1
    # (marker :start is in seconds; ep_start / ep_end are 1-based sample numbers); markers in
    # overlapping windows appear in each epoch
    mrk_s = [t2s(markers[idx, :start], fs) for idx = 1:DataFrames.nrow(markers)]
    rows = Int64[]
    new_start = Float64[]
    @inbounds for ep_idx = 1:mrk_n, idx in eachindex(mrk_s)
        if ep_start[ep_idx] <= mrk_s[idx] <= ep_end[ep_idx]
            push!(rows, idx)
            push!(
                new_start,
                s2t((ep_idx - 1) * ep_len + mrk_s[idx] - ep_start[ep_idx] + 1, fs),
            )
        end
    end
    mrk_tmp = markers[rows, :]
    mrk_tmp[!, :start] = new_start
    return epochs, mrk_tmp
end

function _epochs_tps(obj::NeuroAnalyzer.NEURO)::Matrix{Float64}
    ep_l = epoch_len(obj)
    ep_n = nepochs(obj)
    tps = zeros(2, ep_n)

    # epoch start times: sample indices 1, ep_l+1, 2ep_l+1, …
    tps[1, :] = obj.time_pts[1:ep_l:(ep_n * ep_l - ep_l + 1)]

    # epoch end times: sample indices ep_l, 2ep_l, 3ep_l, …
    tps[2, :] = obj.time_pts[ep_l:ep_l:(ep_n * ep_l)]

    return tps
end

function _markers_epochs(obj::NeuroAnalyzer.NEURO)::Vector{Int64}
    mrk_start = obj.markers[!, :start]
    mrk_epoch = zeros(Int64, length(mrk_start))
    ep_tps = _epochs_tps(obj)
    # cache: avoids repeated header lookups in the inner loop
    ep_n = nepochs(obj)

    for mrk_idx in eachindex(mrk_start)
        # cache: avoids repeated indexing in the inner loop
        t = mrk_start[mrk_idx]
        for ep_idx = 1:ep_n
            if t >= ep_tps[1, ep_idx] && t <= ep_tps[2, ep_idx]
                mrk_epoch[mrk_idx] = ep_idx
                break # a marker belongs to exactly one epoch; no need to continue
            end
        end
    end

    return mrk_epoch
end

# thresholds per channel type: a number applies to the channels of the data type of the object
function _ptp_thresholds(
    obj::NeuroAnalyzer.NEURO,
    t::Union{Nothing, Real, AbstractDict},
    name::String,
)::Dict{String, Float64}
    isnothing(t) && return Dict{String, Float64}()
    d = t isa Real ? Dict(datatype(obj) => Float64(t)) :
        Dict{String, Float64}(String(k) => Float64(v) for (k, v) in t)
    ct = obj.header.recording[:channel_type]
    for (k, v) in d
        k in ct || throw(
            ArgumentError(
                "$name: no channels of type $k" *
                (t isa Real ? " (give the thresholds per channel type, e.g. Dict(\"mag\" => …, \"grad\" => …))." : "."),
            ),
        )
        v > 0 || throw(ArgumentError("$name thresholds must be > 0."))
    end
    return d
end

"""
    _ptp_bad_epochs(obj, reject, flat)

Epochs with a peak-to-peak amplitude above `reject` or below `flat` in at least one channel of the given type(s) (as `reject` and `flat` of MNE-Python); channels marked as bad are not checked.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: epoched object
- `reject::Union{Nothing, Real, AbstractDict}`: maximum peak-to-peak amplitude, a number (channels of the data type of the object) or a dictionary of channel types and thresholds (in the units of the channels)
- `flat::Union{Nothing, Real, AbstractDict}`: minimum peak-to-peak amplitude (as `reject`)

# Returns

Named tuple:

- `ep::Vector{Int64}`: epochs to drop
- `n_reject::Int64`: number of epochs above `reject`
- `n_flat::Int64`: number of epochs below `flat` (and not above `reject`)
"""
function _ptp_bad_epochs(
    obj::NeuroAnalyzer.NEURO,
    reject::Union{Nothing, Real, AbstractDict},
    flat::Union{Nothing, Real, AbstractDict},
)::@NamedTuple{ep::Vector{Int64}, n_reject::Int64, n_flat::Int64}
    rj = _ptp_thresholds(obj, reject, "reject")
    fl = _ptp_thresholds(obj, flat, "flat")
    ct = obj.header.recording[:channel_type]
    bad = obj.header.recording[:bad_channel]
    # peak-to-peak amplitude of the good channels of type k, shape (channels, epochs)
    function p2p(k)
        ch = findall(i -> ct[i] == k && !bad[i], eachindex(ct))
        isempty(ch) && return zeros(0, nepochs(obj))
        d = @view obj.data[ch, :, :]
        return dropdims(maximum(d; dims = 2) .- minimum(d; dims = 2); dims = 2)
    end
    is_rj = falses(nepochs(obj))
    for (k, t) in rj
        is_rj .|= vec(any(p2p(k) .> t; dims = 1))
    end
    is_fl = falses(nepochs(obj))
    for (k, t) in fl
        is_fl .|= vec(any(p2p(k) .< t; dims = 1))
    end
    return (
        ep = findall(is_rj .| is_fl),
        n_reject = count(is_rj),
        n_flat = count(is_fl .& .!is_rj),
    )
end
