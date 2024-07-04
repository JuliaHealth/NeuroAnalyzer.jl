export plinterpolate_channel
export plinterpolate_channel!
export plinterpolate

"""
    plinterpolate_channel(obj; ch, ep, m, q)

Interpolate channel using planar interpolation.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Int64`: channel number to interpolate
- `ep::Union{Int64, Vector{Int64}, <:AbstractRange}`: epoch number(s) within to interpolate
- `imethod::Symbol=:sh`: interpolation method:
    - `:sh`: Shepard
    - `:mq`: Multiquadratic
    - `:imq`: InverseMultiquadratic
    - `:tp`: ThinPlate
    - `:nn`: NearestNeighbour
    - `:ga`: Gaussian
- `ifactor::Int64=100`: interpolation quality

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function plinterpolate_channel(obj::NeuroAnalyzer.NEURO; ch::Int64, ep::Union{Int64, Vector{Int64}, <:AbstractRange}, imethod::Symbol=:sh, ifactor::Int64=100)

    channels = get_channel_bytype(obj, type=obj.header.recording[:data_type])
    @assert length(channels) > 1 "OBJ must contain > 1 signal channel."
    @assert ch in channels "ch must be a signal channel; cannot interpolate non-signal channels."

    _check_var(imethod, [:sh, :mq, :imq, :tp, :nn, :ga], "imethod")
    @assert _has_locs(obj) "Electrode locations not available, use load_locs() or add_locs() first."

    _check_channels(obj, ch)
    _check_epochs(obj, ep)
    isa(ep, Int64) && (ep = [ep])

    obj_new = deepcopy(obj)
    obj_tmp = deepcopy(obj)
    delete_channel!(obj_tmp, ch=get_channel_bytype(obj_tmp, type="ref"))
    delete_channel!(obj_tmp, ch=get_channel_bytype(obj_tmp, type="eog"))

    locs_x1 = obj_tmp.locs[!, :loc_x]
    locs_y1 = obj_tmp.locs[!, :loc_y]

    delete_channel!(obj_tmp, ch=ch)
    locs_x2 = obj_tmp.locs[!, :loc_x]
    locs_y2 = obj_tmp.locs[!, :loc_y]
    chs = get_channel_bytype(obj_tmp, type=obj.header.recording[:data_type])

    ep_n = length(ep)
    ep_len = epoch_len(obj_tmp)

    s_interpolated = zeros(Float64, length(ch), ep_len, ep_n)

    # initialize progress bar
    progress_bar && (progbar = Progress(ep_n * ep_len, dt=1, barlen=20, color=:white))

    @inbounds for ep_idx in eachindex(ep)
        Threads.@threads for length_idx in 1:ep_len
            s_tmp, x, y = @views _interpolate2d(obj_tmp.data[chs, length_idx, ep[ep_idx]], locs_x2, locs_y2, ifactor, imethod, :none)
            for ch_idx in eachindex(ch)
                x_idx = vsearch(locs_x1[ch[ch_idx]], x)
                y_idx = vsearch(locs_y1[ch[ch_idx]], y)
                s_interpolated[ch_idx, length_idx, ep_idx] = s_tmp[x_idx, y_idx]
            end

            # update progress bar
            progress_bar && next!(progbar)
        end
    end

    obj_new.data[ch, :, ep] = s_interpolated

    reset_components!(obj_new)
    push!(obj_new.history, "plinterpolate_channel(OBJ, ch=$ch, ep=$ep, imethod=$imethod, ifactor=$ifactor)")

    return obj_new

end

"""
    plinterpolate_channel!(obj; ch, ep, imethod, ifactor)

Interpolate channel using planar interpolation.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Int64`: channel number to interpolate
- `ep::Union{Int64, Vector{Int64}, <:AbstractRange}`: epoch number(s) within to interpolate
- `imethod::Symbol=:sh`: interpolation method Shepard (`:sh`), Multiquadratic (`:mq`), InverseMultiquadratic (`:imq`), ThinPlate (`:tp`), NearestNeighbour (`:nn`), Gaussian (`:ga`)
- `ifactor::Int64=100`: interpolation quality
"""
function plinterpolate_channel!(obj::NeuroAnalyzer.NEURO; ch::Int64, ep::Union{Int64, Vector{Int64}, <:AbstractRange}, imethod::Symbol=:shepard, ifactor::Int64=100)

    obj_new = plinterpolate_channel(obj, ch=ch, ep=ep, imethod=imethod, ifactor=ifactor)
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components

    return nothing

end

"""
    plinterpolate(s; locs, ch, imethod, nmethod, cart)

Interpolate channel using planar interpolation.

# Arguments

- `s::Matrix{Float64}`: values to plot (one value per channel)
- `locs::DataFrame`: columns: channel, labels, loc_radius, loc_theta, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
- `ch::Int64`: channel number to interpolate
- `imethod::Symbol=:sh`: interpolation method:
    - `:sh`: Shepard
    - `:mq`: Multiquadratic
    - `:imq`: InverseMultiquadratic
    - `:tp`: ThinPlate
    - `:nn`: NearestNeighbour
    - `:ga`: Gaussian
- `nmethod::Symbol=:minmax`: method for normalization, see `normalize()`
- `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates for XY plane and spherical coordinates for XZ and YZ planes
- `ifactor::Int64=100`: interpolation quality

# Returns

- `int_s::Matrix{Float64}`: interpolated signal
- `int_x::Vector{Float64}`: X-axis coordinates
- `int_y::Vector{Float64}`: Y-axis coordinates
"""
function plinterpolate(s::Matrix{Float64}; locs::DataFrame, ch::Int64, imethod::Symbol=:sh, nmethod::Symbol=:minmax, cart::Bool=false, ifactor::Int64=100)

    @assert ch in 1:size(s, 1) "ch must be in [1, $(size(s, 1))"
    _check_var(imethod, [:sh, :mq, :imq, :tp, :nn, :ga], "imethod")

    locs = locs[ch, :]

    if !cart
        loc_x = zeros(length(ch))
        loc_y = zeros(length(ch))
        for idx in eachindex(ch)
            loc_x[idx], loc_y[idx] = pol2cart(locs[!, :loc_radius][idx], locs[!, :loc_theta][idx])
        end
    else
        loc_x = locs[ch, :loc_x]
        loc_y = locs[ch, :loc_y]
    end

    loc_x = _s2v(loc_x)
    loc_y = _s2v(loc_y)

    s_interpolated, interpolated_x, interpolated_y = _interpolate2d(s, loc_x, loc_y, ifactor, imethod, nmethod)

    return (int_s=s_interpolated, int_x = interpolated_x, int_y = interpolated_y)

end