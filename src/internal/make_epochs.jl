function _make_epochs(s::Matrix{<:Real}; ep_n::Union{Int64, Nothing}=nothing, ep_len::Union{Int64, Nothing}=nothing)

    (ep_len === nothing && ep_n === nothing) && throw(ArgumentError("Either ep_n or ep_len must be specified."))
    (ep_len !== nothing && ep_n !== nothing) && throw(ArgumentError("Both ep_n and ep_len cannot be specified."))
    (ep_len !== nothing && ep_len < 1) && throw(ArgumentError("ep_len must be ≥ 1."))
    (ep_n !== nothing && ep_n < 1) && throw(ArgumentError("ep_n must be ≥ 1."))

    ch_n = size(s, 1)
    if ep_n === nothing
        ep_n = size(s, 2) ÷ ep_len
    else
        ep_len = size(s, 2) ÷ ep_n
    end

    ep_len > size(s, 2) && throw(ArgumentError("ep_len must be ≤ $(size(s, 2))."))
    ep_len < 1 && throw(ArgumentError("ep_len must be ≥ 1."))
    ep_n < 1 && throw(ArgumentError("ep_n must be ≥ 1."))

    epochs = reshape(s[:, 1:(ep_len * ep_n)], ch_n, ep_len, ep_n)

    return epochs
end

function _make_epochs(s::Array{<:Real, 3}; ep_n::Union{Int64, Nothing}=nothing, ep_len::Union{Int64, Nothing}=nothing)

    (ep_len === nothing && ep_n === nothing) && throw(ArgumentError("Either ep_n or ep_len must be specified."))
    (ep_len !== nothing && ep_n !== nothing) && throw(ArgumentError("Both ep_n and ep_len cannot be specified."))
    (ep_len !== nothing && ep_len < 1) && throw(ArgumentError("ep_len must be ≥ 1."))
    (ep_n !== nothing && ep_n < 1) && throw(ArgumentError("ep_n must be ≥ 1."))

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

    offset < 1 && throw(ArgumentError("offset must be ≥ 1."))
    ep_len < 1 && throw(ArgumentError("ep_len must be ≥ 1."))
    (offset + ep_len > size(s, 2)) && throw(ArgumentError("offset + ep_len must be ≤ signal length $(size(s, 2))."))

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
        within_epoch == false && deleteat!(markers, mrk_idx)
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
