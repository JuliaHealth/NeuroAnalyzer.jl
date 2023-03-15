export plinterpolate_channel
export plinterpolate_channel!

"""
    plinterpolate_channel(obj; ch, ep, m, q)

Interpolate channel(s) using planar interpolation.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}`: channel number(s) to interpolate
- `ep::Union{Int64, Vector{Int64}, <:AbstractRange}`: epoch number(s) within to interpolate
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
function plinterpolate_channel(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}, ep::Union{Int64, Vector{Int64}, <:AbstractRange}, imethod::Symbol=:sh, interpolation_factor::Int64=100)

    for idx in ch
        idx in get_channel_bytype(obj, type=Symbol(obj.header.recording[:data_type])) || throw(ArgumentError("channel must be signal channel(s); cannot interpolate non-signal channels."))
    end

    _check_var(imethod, [:sh, :mq, :imq, :tp, :nn, :ga], "imethod")
    obj.header.has_locs == false && throw(ArgumentError("Electrode locations not available, use load_locs() or add_locs() first."))

    typeof(ch) == Vector{Int64} && sort!(ch, rev=true)

    obj_new = deepcopy(obj)
    obj_tmp = deepcopy(obj)
    _check_channels(obj, ch)
    _check_epochs(obj, ep)

    locs_x1 = obj.locs[!, :loc_x]
    locs_y1 = obj.locs[!, :loc_y]
    
    obj_tmp = delete_channel(obj, ch=ch)
    locs_x2 = obj_tmp.locs[!, :loc_x]
    locs_y2 = obj_tmp.locs[!, :loc_y]
    chs = get_channel_bytype(obj_tmp, type=Symbol(obj.header.recording[:data_type]))

    ep_n = length(ep)
    ep_len = epoch_len(obj_tmp)

    s_interpolated = zeros(Float64, length(ch), ep_len, ep_n)

    # initialize progress bar
    progress_bar == true && (p = Progress(ep_n * ep_len, 1))

    @inbounds @simd for ep_idx in eachindex(ep)
        Threads.@threads for length_idx in 1:ep_len
            s_tmp, x, y = @views _interpolate(obj_tmp.data[chs, length_idx, ep[ep_idx]], locs_x2, locs_y2, interpolation_factor, imethod, :none)
            for ch_idx in eachindex(ch)
                x_idx = vsearch(locs_x1[ch[ch_idx]], x)
                y_idx = vsearch(locs_y1[ch[ch_idx]], y)
                s_interpolated[ch_idx, length_idx, ep_idx] = s_tmp[x_idx, y_idx]
            end

            # update progress bar
            progress_bar == true && next!(p)
        end
    end

    obj_new.data[ch, :, ep] = s_interpolated

    reset_components!(obj_new)
    push!(obj_new.header.history, "plinterpolate_channel(OBJ, ch=$ch, ep=$ep, imethod=$imethod, interpolation_factor=$interpolation_factor)")

    return obj_new

end

"""
    plinterpolate_channel!(obj; ch, ep, imethod, interpolation_factor)

Interpolate channel(s) using planar interpolation.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}`: channel number(s) to interpolate
- `ep::Union{Int64, Vector{Int64}, <:AbstractRange}`: epoch number(s) within to interpolate
- `imethod::Symbol=:sh`: interpolation method Shepard (`:sh`), Multiquadratic (`:mq`), InverseMultiquadratic (`:imq`), ThinPlate (`:tp`), NearestNeighbour (`:nn`), Gaussian (`:ga`)
- `interpolation_factor::Int64=100`: interpolation quality
"""
function plinterpolate_channel!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}}, ep::Union{Int64, Vector{Int64}, <:AbstractRange}, imethod::Symbol=:shepard, interpolation_factor::Int64=100)

    obj_new = plinterpolate_channel(obj, ch=ch, ep=ep, imethod=imethod, interpolation_factor=interpolation_factor)
    obj.data = obj_new.data
    obj.header = obj_new.header
    obj.components = obj_new.components

    return nothing
    
end
