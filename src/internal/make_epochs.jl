function _make_epochs(signal::Matrix{<:Real}; ep_n::Union{Int64, Nothing}=nothing, ep_len::Union{Int64, Nothing}=nothing)

    (ep_len === nothing && ep_n === nothing) && throw(ArgumentError("Either ep_n or ep_len must be specified."))
    (ep_len !== nothing && ep_n !== nothing) && throw(ArgumentError("Both ep_n and ep_len cannot be specified."))
    (ep_len !== nothing && ep_len < 1) && throw(ArgumentError("ep_len must be ≥ 1."))
    (ep_n !== nothing && ep_n < 1) && throw(ArgumentError("ep_n must be ≥ 1."))

    ch_n = size(signal, 1)
    if ep_n === nothing
        ep_n = size(signal, 2) ÷ ep_len
    else
        ep_len = size(signal, 2) ÷ ep_n
    end

    ep_len > size(signal, 2) && throw(ArgumentError("ep_len must be ≤ $(size(signal, 2))."))
    ep_len < 1 && throw(ArgumentError("ep_len must be ≥ 1."))
    ep_n < 1 && throw(ArgumentError("ep_n must be ≥ 1."))

    epochs = reshape(signal[:, 1:(ep_len * ep_n)], ch_n, ep_len, ep_n)

    return epochs
end

function _make_epochs(signal::Array{<:Real, 3}; ep_n::Union{Int64, Nothing}=nothing, ep_len::Union{Int64, Nothing}=nothing)

    (ep_len === nothing && ep_n === nothing) && throw(ArgumentError("Either ep_n or ep_len must be specified."))
    (ep_len !== nothing && ep_n !== nothing) && throw(ArgumentError("Both ep_n and ep_len cannot be specified."))
    (ep_len !== nothing && ep_len < 1) && throw(ArgumentError("ep_len must be ≥ 1."))
    (ep_n !== nothing && ep_n < 1) && throw(ArgumentError("ep_n must be ≥ 1."))

    ch_n = size(signal, 1)
    if ep_n === nothing
        ep_n = size(signal, 2) * size(signal, 3) ÷ ep_len
    else
        ep_len = size(signal, 2) * size(signal, 3) ÷ ep_n
    end
    signal = reshape(signal, ch_n, (size(signal, 2) * size(signal, 3)), 1)
    epochs = reshape(signal[:, 1:(ep_len * ep_n), 1], ch_n, ep_len, ep_n)

    return epochs
end

function _make_epochs_bymarkers(signal::Array{<:Real, 3}; markers::DataFrame, marker_start::Vector{Int64}, epoch_offset::Int64, ep_len::Int64)

    if size(signal, 3) > 1
        _info("Signal has already been epoched, parts of the signal might have been removed.")
        signal = reshape(signal, ch_n, (size(signal, 2) * size(signal, 3)), 1)
    end

    marker_n = length(marker_start)
    epoch_start = marker_start .- epoch_offset
    epoch_end = epoch_start .+ ep_len .- 1

    # delete epochs outside signal limits
    for marker_idx in marker_n:-1:1
        if epoch_start[marker_idx] < 1
            deleteat!(epoch_start, marker_idx)
            deleteat!(epoch_end, marker_idx)
            marker_n -= 1
        elseif epoch_end[marker_idx] > size(signal, 2)
            deleteat!(epoch_start, marker_idx)
            deleteat!(epoch_end, marker_idx)
            marker_n -= 1
        end
    end

    (epoch_offset < 1) && throw(ArgumentError("epoch_offset must be ≥ 1."))
    (ep_len !== nothing && ep_len < 1) && throw(ArgumentError("ep_len must be ≥ 1."))

    (epoch_offset + ep_len > size(signal, 2)) && throw(ArgumentError("epoch_offset + ep_len must be ≤ signal length $(size(signal, 2))."))

    epochs = zeros(size(signal, 1), ep_len, marker_n)
    for marker_idx in 1:marker_n
        epochs[:, :, marker_idx] = reshape(signal[:, epoch_start[marker_idx]:epoch_end[marker_idx], :],
                                           size(signal, 1), 
                                           ep_len)
    end

    # remove markers outside epoch limits
    for marker_idx in nrow(markers):-1:1
        within_epoch = false
        for epoch_idx in 1:marker_n
            markers[!, :start][marker_idx] in epoch_start[epoch_idx]:epoch_end[epoch_idx] && (within_epoch = true)
        end
        within_epoch == false && deleteat!(markers, marker_idx)
    end

    # calculate marker offsets
    for epoch_idx in 1:marker_n
        for marker_idx in 1:nrow(markers)
            if markers[!, :start][marker_idx] in epoch_start[epoch_idx]:epoch_end[epoch_idx]
                markers[!, :start][marker_idx] = (epoch_idx - 1) * ep_len + markers[!, :start][marker_idx] .- epoch_start[epoch_idx]
            end
        end
    end

    return epochs, markers
end
