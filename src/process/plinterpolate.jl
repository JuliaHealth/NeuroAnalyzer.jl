export plinterpolate_channel
export plinterpolate_channel!

"""
    plinterpolate_channel(obj; channel, epoch, m, q)

Interpolate OBJ channel(s) using planar interpolation.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel number(s) to interpolate
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}`: epoch number(s) within to interpolate
- `imethod::Symbol=:sh`: interpolation method:
    - `:sh`: Shepard
    - `:mq`: Multiquadratic
    - `:imq`: InverseMultiquadratic
    - `:tp`: ThinPlate
    - `:nn`: NearestNeighbour
    - `:ga`: Gaussian
- `interpolation_factor::Int64=100`: interpolation quality

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function plinterpolate_channel(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}, epoch::Union{Int64, Vector{Int64}, AbstractRange}, imethod::Symbol=:sh, interpolation_factor::Int64=100)

    for idx in channel
        idx in get_channel_bytype(obj, type=Symbol(obj.header.recording[:data_type])) || throw(ArgumentError("channel must be OBJ/MEG signal channel(s); cannot interpolate non-OBJ/MEG channels."))
    end

    _check_var(imethod, [:sh, :mq, :imq, :tp, :nn, :ga], "imethod")
    obj.header.recording[:channel_locations] == false && throw(ArgumentError("Electrode locations not available, use load_electrodes() or add_electrodes() first."))

    typeof(channel) == Vector{Int64} && sort!(channel, rev=true)

    obj_new = deepcopy(obj)
    obj_tmp = deepcopy(obj)
    _check_channels(obj, channel)
    _check_epochs(obj, epoch)

    locs_x1 = obj.locs[!, :loc_x]
    locs_y1 = obj.locs[!, :loc_y]
    
    obj_tmp = delete_channel(obj, channel=channel)
    locs_x2 = obj_tmp.locs[!, :loc_x]
    locs_y2 = obj_tmp.locs[!, :loc_y]
    channels = get_channel_bytype(obj_tmp, type=Symbol(obj.header.recording[:data_type]))

    ep_n = length(epoch)
    ep_len = epoch_len(obj_tmp)

    s_interpolated = zeros(Float64, length(channel), ep_len, ep_n)

    # initialize progress bar
    progress_bar == true && (p = Progress(ep_n * ep_len, 1))

    @inbounds @simd for epoch_idx in eachindex(epoch)
        Threads.@threads for length_idx in 1:ep_len
            s_tmp, x, y = @views _interpolate(obj_tmp.data[channels, length_idx, epoch[epoch_idx]], locs_x2, locs_y2, interpolation_factor, imethod, :none)
            for channel_idx in eachindex(channel)
                x_idx = vsearch(locs_x1[channel[channel_idx]], x)
                y_idx = vsearch(locs_y1[channel[channel_idx]], y)
                s_interpolated[channel_idx, length_idx, epoch_idx] = s_tmp[x_idx, y_idx]
            end

            # update progress bar
            progress_bar == true && next!(p)
        end
    end

    obj_new.data[channel, :, epoch] = s_interpolated

    reset_components!(obj_new)
    push!(obj_new.header.history, "plinterpolate_channel(OBJ, channel=$channel, epoch=$epoch, imethod=$imethod, interpolation_factor=$interpolation_factor)")

    return obj_new
end

"""
    plinterpolate_channel!(obj; channel, epoch, imethod, interpolation_factor)

Interpolate OBJ channel(s) using planar interpolation.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel number(s) to interpolate
- `epoch::Union{Int64, Vector{Int64}, AbstractRange}`: epoch number(s) within to interpolate
- `imethod::Symbol=:sh`: interpolation method Shepard (`:sh`), Multiquadratic (`:mq`), InverseMultiquadratic (`:imq`), ThinPlate (`:tp`), NearestNeighbour (`:nn`), Gaussian (`:ga`)
- `interpolation_factor::Int64=100`: interpolation quality
"""
function plinterpolate_channel!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}}, epoch::Union{Int64, Vector{Int64}, AbstractRange}, imethod::Symbol=:shepard, interpolation_factor::Int64=100)

    obj_tmp = plinterpolate_channel(obj, channel=channel, epoch=epoch, imethod=imethod, interpolation_factor=interpolation_factor)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    reset_components!(obj)

    return nothing
end
