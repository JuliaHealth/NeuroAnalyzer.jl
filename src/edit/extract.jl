export extract_channel
export extract_epoch
export extract_epoch!
export extract_data
export extract_time
export extract_eptime

"""
    extract_channel(obj; ch)

Extract channel data.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, String}`: channel number or name

# Returns

- `extract_channel::Vector{Float64}`
"""
function extract_channel(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, String})

    clabels = labels(obj)
    if typeof(ch) == String
        # get ch by name
        ch_idx = nothing
        for idx in eachindex(clabels)
            if ch == clabels[idx]
                ch_idx = idx
            end
        end
        if ch_idx === nothing
            throw(ArgumentError("Channel name ($ch )does not match channel labels."))
        end
        return reshape(obj.data[ch_idx, :, :], 1, epoch_len(obj), epoch_n(obj))
    else
        # get channel by number
        _check_channels(obj, ch)
        return reshape(obj.data[ch, :, :], 1, epoch_len(obj), epoch_n(obj))
    end

end

"""
    extract_epoch(obj; ep)

Extract epoch.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ep::Int64`: epoch index

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function extract_epoch(obj::NeuroAnalyzer.NEURO; ep::Int64)

    _check_epochs(obj, ep)

    obj_new = deepcopy(obj)

    obj_new.data = reshape(obj.data[:, :, ep], channel_n(obj), epoch_len(obj), 1)
    obj_new.time_pts = obj.epoch_time
    obj_new.epoch_time = obj.epoch_time

    reset_components!(obj_new)
    push!(obj_new.history, "extract_epoch(OBJ, ep=$ep)")

    return obj_new

end

"""
    extract_epoch!(obj; ep)

Extract epoch.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ep::Int64`: epoch index
"""
function extract_epoch!(obj::NeuroAnalyzer.NEURO; ep::Int64)

    obj_new = extract_epoch(obj, ep=ep)
    obj.header = obj_new.header
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components
    obj.time_pts = obj_new.time_pts

    return nothing

end

"""
    extract_data(obj; ch)

Extract data.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `ep::Union{Int64, Vector{Int64}, <:AbstractRange}=1:epoch_n(obj)`: index of epochs, default is all epochs
- `time::Bool=false`: return time vector
- `etime::Bool=false`: return epoch time vector

# Returns

- `signal::Array{Float64, 3}`
- `time::Vector{Float64}`
- `etime::Vector{Float64}`
"""
function extract_data(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), ep::Union{Int64, Vector{Int64}, <:AbstractRange}=1:epoch_n(obj), time::Bool=false, etime::Bool=false)

    _check_channels(obj, ch)
    _check_epochs(obj, ep)

    if time == false && etime == false
        return obj.data[ch, :, ep][:, :, :]
    elseif time == true && etime == false
        return obj.data[ch, :, ep][:, :, :], obj.time_pts
    elseif time == false && etime == true
        return obj.data[ch, :, ep][:, :, :], obj.epoch_time
    else
        return obj.data[ch, :, ep][:, :, :], obj.time_pts, obj.epoch_time
    end

end

"""
    extract_time(obj)

Extract time.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `tpts::Array{Float64, 3}`
"""
function extract_time(obj::NeuroAnalyzer.NEURO)

    tpts = obj.time_pts
    
    return tpts

end

"""
    extract_eptime(obj)

Extract epochs time.

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `et::Array{Float64, 3}`
"""
function extract_eptime(obj::NeuroAnalyzer.NEURO)

    et = obj.epoch_time

    return et

end
