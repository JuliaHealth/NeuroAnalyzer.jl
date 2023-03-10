export extract_channel
export extract_epoch
export extract_data
export extract_time
export extract_etime

"""
    extract_channel(obj; channel)

Extract channel data.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, String}`: channel number or name

# Returns

- `extract_channel::Vector{Float64}`
"""
function extract_channel(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, String})

    clabels = labels(obj)
    if typeof(channel) == String
        # get channel by name
        ch_idx = nothing
        for idx in eachindex(clabels)
            if channel == clabels[idx]
                ch_idx = idx
            end
        end
        if ch_idx === nothing
            throw(ArgumentError("Channel name ($channel )does not match channel labels."))
        end
        return reshape(obj.data[ch_idx, :, :], 1, epoch_len(obj), epoch_n(obj))
    else
        # get channel by number
        _check_channels(obj, channel)
        return reshape(obj.data[channel, :, :], 1, epoch_len(obj), epoch_n(obj))
    end    
end

"""
    extract_epoch(obj; epoch)

Extract epoch.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `epoch::Int64`: epoch index

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function extract_epoch(obj::NeuroAnalyzer.NEURO; epoch::Int64)

    _check_epochs(obj, epoch)

    s_new = reshape(obj.data[:, :, epoch], channel_n(obj), signal_len(obj), 1)
    obj_new = deepcopy(obj)
    obj_new.data = s_new
    obj_new.epoch_time = obj.epoch_time
    obj_new.header.recording[:epoch_n] = 1
    obj_new.header.recording[:duration_samples] = obj_new.header.recording[:epoch_duration_samples]
    obj_new.header.recording[:duration_seconds] = obj_new.header.recording[:epoch_duration_seconds]

    reset_components!(obj_new)
    push!(obj_new.header.history, "extract_epoch(OBJ, epoch=$epoch)")

    return obj_new
end

"""
    extract_data(obj; channel)

Extract data.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}=epoch_n(obj)`: index of epochs, default is all epochs
- `time::Bool=false`: return time vector
- `etime::Bool=false`: return epoch time vector

# Returns

- `signal::Array{Float64, 3}`
- `time::Vector{Float64}`
- `etime::Vector{Float64}`
"""
function extract_data(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), epoch::Union{Int64, Vector{Int64}, AbstractRange}=epoch_n(obj), time::Bool=false, etime::Bool=false)
    _check_channels(obj, channel)
    _check_epochs(obj, epoch)
    if time == false && etime == false
        return obj.data[channel, :, epoch][:, :, :]
    elseif time == true && etime == false
        return obj.data[channel, :, epoch][:, :, :], obj.time_pts
    elseif time == false && etime == true
        return obj.data[channel, :, epoch][:, :, :], obj.epoch_time
    else
        return obj.data[channel, :, epoch][:, :, :], obj.time_pts, obj.epoch_time
    end
end

"""
    extract_time(obj)

Extract time.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `time::Array{Float64, 3}`
"""
function extract_time(obj::NeuroAnalyzer.NEURO)
    return obj.time_pts
end

"""
    extract_etime(obj)

Extract epochs time.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `time::Array{Float64, 3}`
"""
function extract_etime(obj::NeuroAnalyzer.NEURO)
    return obj.epoch_time
end
