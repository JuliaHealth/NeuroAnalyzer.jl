export extract_channel
export extract_epoch
export extract_epoch!
export extract_data
export extract_time
export extract_eptime

"""
    extract_channel(obj; <keyword arguments>)

Extract channel data.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::String`: channel name

# Returns

- `d::Vector{Float64}`
"""
function extract_channel(obj::NeuroAnalyzer.NEURO; ch::String)::Vector{Float64}

    ch = get_channel(obj, ch=ch)
    d = reshape(obj.data[ch, :, :], 1, epoch_len(obj), nepochs(obj))

    return d

end

"""
    extract_epoch(obj; <keyword arguments>)

Extract epoch.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ep::Int64`: epoch index

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function extract_epoch(obj::NeuroAnalyzer.NEURO; ep::Int64)::NeuroAnalyzer.NEURO

    _check_epochs(obj, ep)

    obj_new = deepcopy(obj)

    obj_new.data = reshape(obj.data[:, :, ep], nchannels(obj), epoch_len(obj), 1)
    obj_new.time_pts = obj.epoch_time
    obj_new.epoch_time = obj.epoch_time

    reset_components!(obj_new)
    push!(obj_new.history, "extract_epoch(OBJ, ep=$ep)")

    return obj_new

end

"""
    extract_epoch!(obj; <keyword arguments>)

Extract epoch.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ep::Int64`: epoch index

# Returns

Nothing
"""
function extract_epoch!(obj::NeuroAnalyzer.NEURO; ep::Int64)::Nothing

    obj_new = extract_epoch(obj, ep=ep)
    obj.header = obj_new.header
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components
    obj.time_pts = obj_new.time_pts

    return nothing

end

"""
    extract_data(obj; <keyword arguments>)

Extract data.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}}`: channel name or list of channel names
- `ep::Union{Int64, Vector{Int64}, <:AbstractRange}=1:nepochs(obj)`: index of epochs, default is all epochs
- `time::Bool=false`: return time vector
- `etime::Bool=false`: return epoch time vector

# Returns

- `signal::Array{Float64, 3}`
- `time::Vector{Float64}`
- `etime::Vector{Float64}`
"""
function extract_data(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}}, ep::Union{Int64, Vector{Int64}, <:AbstractRange}=1:nepochs(obj), time::Bool=false, etime::Bool=false)::Union{Array{Float64, 3}, Tuple{Array{Float64, 3}, Vector{Float64}}, Tuple{Array{Float64, 3}, Vector{Float64}, Vector{Float64}}}

    ch = get_channel(obj, ch=ch)
    _check_epochs(obj, ep)
    isa(ep, Int64) && (ep = [ep])

    if !time && !etime
        return obj.data[ch, :, ep][:, :, :]
    elseif time && !etime
        return obj.data[ch, :, ep][:, :, :], obj.time_pts
    elseif !time && etime
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
function extract_time(obj::NeuroAnalyzer.NEURO)::Array{Float64, 3}

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
function extract_eptime(obj::NeuroAnalyzer.NEURO)::Array{Float64, 3}

    et = obj.epoch_time

    return et

end
