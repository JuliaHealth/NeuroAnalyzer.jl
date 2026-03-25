export extract_channel
export extract_epoch
export extract_epoch!
export extract_data

"""
    extract_channel(obj; <keyword arguments>)

Extract channel data.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::String`: channel name; must resolve to exactly one channel

# Returns

- `Array{Float64, 3}`
"""
function extract_channel(obj::NeuroAnalyzer.NEURO; ch::String)::Array{Float64, 3}

    # resolve channel names to integer indices
    ch = get_channel(obj, ch = ch)
    length(ch) == 1 || throw(ArgumentError("ch must resolve to exactly one channel."))
    ch = ch[1]

    return reshape(obj.data[ch, :, :], 1, epoch_len(obj), nepochs(obj))

end

"""
    extract_epoch(obj; <keyword arguments>)

Extract epoch.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ep::Int64`: epoch index

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function extract_epoch(obj::NeuroAnalyzer.NEURO; ep::Int64)::NeuroAnalyzer.NEURO

    # validate
    _check_epochs(obj, ep)

    # create new dataset
    obj_new = deepcopy(obj)

    obj_new.data = reshape(obj.data[:, :, ep], nchannels(obj), epoch_len(obj), 1)
    obj_new.time_pts = obj.epoch_time
    obj_new.epoch_time = obj.epoch_time

    push!(obj_new.history, "extract_epoch(OBJ, ep=$ep)")

    return obj_new

end

"""
    extract_epoch!(obj; <keyword arguments>)

Extract epoch.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ep::Int64`: epoch index

# Returns

- `Nothing`
"""
function extract_epoch!(obj::NeuroAnalyzer.NEURO; ep::Int64)::Nothing

    obj_new = extract_epoch(obj, ep = ep)
    obj.header = obj_new.header
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.time_pts = obj_new.time_pts

    return nothing

end

"""
    extract_data(obj; <keyword arguments>)

Extract data.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `ep::Union{Int64, Vector{Int64}, AbstractRange}=1:nepochs(obj)`: index of epochs, default is all epochs
- `time::Bool=false`: return time vector
- `etime::Bool=false`: return epoch time vector

# Returns

- `Array{Float64, 3}`
- `Vector{Float64}`
- `Vector{Float64}`
"""
function extract_data(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    ep::Union{Int64, Vector{Int64}, AbstractRange} = 1:nepochs(obj),
    time::Bool = false,
    etime::Bool = false
)::Union{
    Array{Float64, 3},
    Tuple{Array{Float64, 3}, Vector{Float64}},
    Tuple{Array{Float64, 3}, Vector{Float64}, Vector{Float64}},
}

    # resolve channel names to integer indices
    ch = get_channel(obj, ch = ch)

    # validate
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
