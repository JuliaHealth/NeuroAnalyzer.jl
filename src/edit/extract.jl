export extract_data
export extract_time
export extract_etime

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
