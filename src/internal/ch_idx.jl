function _get_channel_idx(clabels::Vector{String}, channel::Union{String, Int64})
    if typeof(channel) == String
        channel_found = nothing
        for idx in eachindex(clabels)
            if channel == clabels[idx]
                channel_found = idx
            end
        end
        if channel_found === nothing
            throw(ArgumentError("channel name does not match signal labels."))
        end
    else
        channel < 1 || channel > length(clabels) && throw(ArgumentError("channel index does not match signal channels."))
        channel_found = channel
    end

    return channel_found
end
