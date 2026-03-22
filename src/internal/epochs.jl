function _make_epochs(
    s::AbstractMatrix;
    ep_len::Int64
)::Array{Float64, 3}
    !(ep_len >= 1) && throw(ArgumentError("ep_len must be ≥ 1."))
    !(ep_len <= size(s, 2)) && throw(ArgumentError("ep_len must be ≤ $(size(s, 2))."))

    ch_n  = size(s, 1)
    ep_n  = size(s, 2) ÷ ep_len

    # trim trailing samples that don't fill a complete epoch, then reshape
    return reshape(s[:, 1:(ep_len * ep_n)], ch_n, ep_len, ep_n)
end

function _make_epochs(
    s::AbstractArray;
    ep_len::Int64
)::Array{Float64, 3}
    _chk3d(s)
    !(ep_len >= 1           ) && throw(ArgumentError("ep_len must be ≥ 1."))
    !(ep_len <= size(s, 2)  ) && throw(ArgumentError("ep_len must be ≤ $(size(s, 2))."))

    ch_n    = size(s, 1)
    n_samp  = size(s, 2) * size(s, 3)
    ep_n    = n_samp ÷ ep_len

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
    fs::Int64
)::Tuple{Array{Float64, 3}, DataFrame}

    _chk3d(s)
    !(offset >= 0) && throw(ArgumentError("offset must be ≥ 0."))  # was: ≥ 1 — offset=0 is valid (epoch starts at marker)
    !(ep_len >= 1) && throw(ArgumentError("ep_len must be ≥ 1."))
    !(fs >= 1) && throw(ArgumentError("fs must be ≥ 1."))

    sig_len = size(s, 2)
    mrk_n = length(marker_start)

    # compute epoch sample boundaries
    ep_start = marker_start .- offset
    ep_end = ep_start .+ ep_len .- 1

    # remove epochs that fall (even partially) outside the signal 
    # iterate backwards to allow safe deleteat!
    for idx in mrk_n:-1:1
        if ep_start[idx] < 1 || ep_end[idx] > sig_len
            deleteat!(ep_start, idx)
            deleteat!(ep_end,   idx)
        end
    end

    mrk_n  = length(ep_start)
    epochs = zeros(size(s, 1), ep_len, mrk_n)

    @inbounds for mrk_idx in 1:mrk_n
        # flatten the epoch dimension of s before slicing (s is 3-D with 1 epoch)
        epochs[:, :, mrk_idx] = reshape(
            s[:, ep_start[mrk_idx]:ep_end[mrk_idx], :],
            size(s, 1), ep_len,
        )
    end

    # remove markers that lie outside all extracted epoch windows
    @inbounds for mrk_idx in DataFrames.nrow(markers):-1:1
        within = any(
            _in(markers[mrk_idx, :start] * fs, (ep_start[ep_idx], ep_end[ep_idx]))
            for ep_idx in 1:mrk_n
        )
        !within && deleteat!(markers, mrk_idx)
    end

    # keep only markers of the requested type and recalculate
    # their time offsets relative to the start of the first epoch
    mrk_tmp = markers[markers[!, :value] .== marker, :]

    if DataFrames.nrow(mrk_tmp) > 0
        mrk_tmp[1, :start] = offset / fs
        for idx in 2:DataFrames.nrow(mrk_tmp)
            mrk_tmp[idx, :start] = mrk_tmp[idx - 1, :start] + ep_len / fs
        end
    end

    # remove stale commented-out block (marker offset recalculation was replaced above)
    return epochs, mrk_tmp
end

function _epochs_tps(obj::NeuroAnalyzer.NEURO)::Matrix{Float64}
    ep_l = epoch_len(obj)
    ep_n = nepochs(obj)
    tps = zeros(2, n)

    # epoch start times: sample indices 1, ep_l+1, 2ep_l+1, …
    tps[1, :] = obj.time_pts[1:ep_l:(ep_n * ep_l - ep_l + 1)]

    # epoch end times: sample indices ep_l, 2ep_l, 3ep_l, …
    tps[2, :] = obj.time_pts[ep_l:ep_l:(ep_n * ep_l)]

    return tps
end

function _markers_epochs(obj::NeuroAnalyzer.NEURO)::Vector{Int64}
    mrk_start  = obj.markers[!, :start]
    mrk_epoch  = zeros(Int64, length(mrk_start))
    ep_tps = _ep_tps(obj)
    # cache: avoids repeated header lookups in the inner loop
    ep_n = nepochs(obj)

    for mrk_idx in eachindex(mrk_start)
        # cache: avoids repeated indexing in the inner loop
        t = mrk_start[mrk_idx]
        for ep_idx in 1:ep_n
            if t >= ep_tps[1, ep_idx] && t <= ep_tps[2, ep_idx]
                mrk_epoch[mrk_idx] = ep_idx
                break # a marker belongs to exactly one epoch; no need to continue
            end
        end
    end

    return mrk_epoch
end

