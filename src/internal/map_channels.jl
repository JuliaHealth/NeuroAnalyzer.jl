function _map_channels(channel::Union{Int64, Vector{Int64}, <:AbstractRange}, channels=Vector{Int64})
    channel_orig = channel
    if typeof(channel) == Int64
        channel = vsearch(channel, channels)
    else
        for idx in eachindex(channel)
            channel[idx] = vsearch(channel[idx], channels)
        end
    end
    return channel, channel_orig
end
