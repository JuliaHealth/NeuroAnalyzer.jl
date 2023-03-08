function _gen_clabels(eeg::NeuroAnalyzer.EEG, c::Symbol)
    c, _ = _get_component(eeg, c)
    clabels = Vector{String}()
    for idx in 1:size(c, 1)
        push!(clabels, lpad(string(idx), length(string(size(c, 1))), "0"))
    end
    return clabels
end

function _gen_clabels(c::Array{Float64, 3})
    clabels = Vector{String}()
    for idx in 1:size(c, 1)
        push!(clabels, lpad(string(idx), length(string(size(c, 1))), "0"))
    end
    return clabels
end

function _channel2channel_name(channel::Union{Int64, Vector{Int64}, AbstractRange})
    if typeof(channel) == Int64
        return channel
    else
        if collect(channel[1]:channel[end]) == channel
            channel_name = string(channel[1]) * ":" * string(channel[end])
        else
            channel_name = ""
            for idx in 1:(length(channel) - 1)
                channel_name *= string(channel[idx])
                channel_name *= ", "
            end
            channel_name *= string(channel[end])
        end
    end
    return channel_name
end
