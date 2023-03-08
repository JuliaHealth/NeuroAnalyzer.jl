function _get_channel_idx(labels::Vector{String}, channel::Union{String, Int64})
    if typeof(channel) == String
        channel_found = nothing
        for idx in eachindex(labels)
            if channel == labels[idx]
                channel_found = idx
            end
        end
        if channel_found === nothing
            throw(ArgumentError("channel name does not match signal labels."))
        end
    else
        channel < 1 || channel > length(labels) && throw(ArgumentError("channel index does not match signal channels."))
        channel_found = channel
    end

    return channel_found
end
