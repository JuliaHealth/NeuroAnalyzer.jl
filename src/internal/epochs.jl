function _make_epochs(s::Matrix{<:Real}; ep_n::Union{Int64, Nothing}=nothing, ep_len::Union{Int64, Nothing}=nothing)

    ep_len === nothing && @assert ep_n !== nothing "Either ep_n or ep_len must be specified."
    ep_len !== nothing && @assert ep_n === nothing "Either ep_n or ep_len must be specified."
    ep_len === nothing && @assert ep_n !== nothing "Both ep_n and ep_len cannot be specified."
    ep_n === nothing && @assert ep_len !== nothing "Both ep_n and ep_len cannot be specified."
    @assert !(ep_len !== nothing && ep_len < 1) "ep_len must be ≥ 1."
    @assert !(ep_n !== nothing && ep_n < 1) "ep_n must be ≥ 1."

    ch_n = size(s, 1)
    if ep_n === nothing
        ep_n = size(s, 2) ÷ ep_len
    else
        ep_len = size(s, 2) ÷ ep_n
    end

    @assert ep_len <= size(s, 2) "ep_len must be ≤ $(size(s, 2))."
    @assert ep_len >= 1 "ep_len must be ≥ 1."
    @assert ep_n >= 1 "ep_n must be ≥ 1."

    epochs = reshape(s[:, 1:(ep_len * ep_n)], ch_n, ep_len, ep_n)

    return epochs
end

function _make_epochs(s::Array{<:Real, 3}; ep_n::Union{Int64, Nothing}=nothing, ep_len::Union{Int64, Nothing}=nothing)

    ep_len === nothing && @assert ep_n !== nothing "Either ep_n or ep_len must be specified."
    ep_len !== nothing && @assert ep_n === nothing "Either ep_n or ep_len must be specified."
    ep_len === nothing && @assert ep_n !== nothing "Both ep_n and ep_len cannot be specified."
    ep_n === nothing && @assert ep_len !== nothing "Both ep_n and ep_len cannot be specified."
    @assert !(ep_len !== nothing && ep_len < 1) "ep_len must be ≥ 1."
    @assert !(ep_n !== nothing && ep_n < 1) "ep_n must be ≥ 1."

    ch_n = size(s, 1)
    if ep_n === nothing
        ep_n = size(s, 2) * size(s, 3) ÷ ep_len
    else
        ep_len = size(s, 2) * size(s, 3) ÷ ep_n
    end
    s = reshape(s, ch_n, (size(s, 2) * size(s, 3)), 1)
    epochs = reshape(s[:, 1:(ep_len * ep_n), 1], ch_n, ep_len, ep_n)

    return epochs
end

function _make_epochs_bymarkers(s::Array{<:Real, 3}; marker::String, markers::DataFrame, marker_start::Vector{Int64}, offset::Int64, ep_len::Int64, fs::Int64)

    if size(s, 3) > 1
        _warn("Signal has already been epoched, parts of the signal might have been removed.")
        s = reshape(s, ch_n, (size(s, 2) * size(s, 3)), 1)
    end

    @assert offset >= 1 "offset must be ≥ 1."
    @assert ep_len >= 1 "ep_len must be ≥ 1."
    @assert offset + ep_len <= size(s, 2) "offset + ep_len must be ≤ signal length $(size(s, 2))."

    mrk_n = length(marker_start)
    ep_start = @. marker_start - offset
    ep_end = @. ep_start + ep_len - 1

    # delete epochs outside signal limits
    for mrk_idx in mrk_n:-1:1
        if ep_start[mrk_idx] < 1
            deleteat!(ep_start, mrk_idx)
            deleteat!(ep_end, mrk_idx)
        elseif ep_end[mrk_idx] > size(s, 2)
            deleteat!(ep_start, mrk_idx)
            deleteat!(ep_end, mrk_idx)
        end
    end
    mrk_n = length(ep_start)

    epochs = zeros(size(s, 1), ep_len, mrk_n)
    @inbounds for mrk_idx in 1:mrk_n
        epochs[:, :, mrk_idx] = @views reshape(s[:, ep_start[mrk_idx]:ep_end[mrk_idx], :],
                                               size(s, 1),
                                               ep_len)
    end

    # remove markers outside epoch limits
    @inbounds for mrk_idx in nrow(markers):-1:1
        within_epoch = false
        for ep_idx in 1:mrk_n
            markers[mrk_idx, :start] * fs in ep_start[ep_idx]:ep_end[ep_idx] && (within_epoch = true)
        end
        !within_epoch && deleteat!(markers, mrk_idx)
    end

    # keep only markers of the given type and shift their offsets
    mrk_tmp = markers[markers[!, :description] .== marker, :]
    mrk_tmp[1, :start] = offset / fs
    if nrow(markers) > 1
        @inbounds for idx in 2:nrow(mrk_tmp)
            mrk_tmp[idx, :start] = mrk_tmp[(idx - 1), :start] + (ep_len / fs)
        end
    end

    # calculate marker offsets
    # @inbounds for ep_idx in 1:mrk_n
    #     for mrk_idx in 1:nrow(markers)
    #         if markers[mrk_idx, :start] in ep_start[ep_idx]:ep_end[ep_idx]
    #             markers[mrk_idx, :start] = (ep_idx - 1) * ep_len + markers[mrk_idx, :start] .- ep_start[ep_idx]
    #         end
    #     end
    # end

    return epochs, mrk_tmp

end

function _epochs_tps(obj::NeuroAnalyzer.NEURO)
    # return time borders of epochs
    tps = zeros(2, nepochs(obj))
    # epoch start
    tps[1, :] = obj.time_pts[1:epoch_len(obj):(length(obj.time_pts) - epoch_len(obj) + 1)]
    # epoch end
    tps[2, :] = obj.time_pts[(epoch_len(obj)):epoch_len(obj):length(obj.time_pts)]
    return tps
end

function _markers_epochs(obj::NeuroAnalyzer.NEURO)
    # return epoch numbers of markers
    mrk_start = obj.markers[!, :start]
    mrk_epoch = zeros(Int64, length(mrk_start))
    epochs_tps = NeuroAnalyzer._epochs_tps(obj)
    for mrk_idx in eachindex(mrk_start)
        for ep_idx in 1:nepochs(obj)
            mrk_start[mrk_idx] >= epochs_tps[1, ep_idx] && mrk_start[mrk_idx] <= epochs_tps[2, ep_idx] && (mrk_epoch[mrk_idx] = ep_idx)
        end
    end
    return mrk_epoch
end