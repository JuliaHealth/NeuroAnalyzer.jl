export plot_topo

"""
    plot_topo(c; <keyword arguments>)

Plot topographical view.

# Arguments

- `signal::Vector{<:Real}`: values to plot (one value per channel)
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}`: channel(s) to plot
- `locs::DataFrame`: columns: channel, labels, loc_theta, loc_radius, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
- `cb::Bool=true`: plot color bar
- `cb_label::String="[A.U.]"`: color bar label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `imethod::Symbol=:sh`: interpolation method:
    - `:sh`: Shepard
    - `:mq`: Multiquadratic
    - `:imq`: InverseMultiquadratic
    - `:tp`: ThinPlate
    - `:nn`: NearestNeighbour
    - `:ga`: Gaussian
- `nmethod::Symbol=:minmax`: method for normalization, see `normalize()`
- `plot_size::Int64=800`: plot dimensions in pixels (size × size)
- `plot_contours::Bools=true`: plot contours over topo plot
- `plot_electrodes::Bools=true`: plot electrodes over topo plot
- `head_labels::Bool=false`: plot head labels
- `head_details::Bool=true`: draw nose and ears
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_topo(signal::Vector{<:Real}; ch::Union{Int64, Vector{Int64}, <:AbstractRange}, locs::DataFrame, cb::Bool=true, cb_label::String="[A.U.]", title::String="default", mono::Bool=false, imethod::Symbol=:sh, nmethod::Symbol=:minmax, plot_contours::Bool=true, plot_electrodes::Bool=true, plot_size::Int64=800, head_labels::Bool=false, head_details::Bool=true, kwargs...)
    
    pal = mono == true ? :grays : :darktest
    _check_var(imethod, [:sh, :mq, :imq, :tp, :nn, :ga], "imethod")

    loc_x = zeros(size(locs, 1))
    loc_y = zeros(size(locs, 1))
    for idx in 1:size(locs, 1)
        loc_x[idx], loc_y[idx] = pol2cart(locs[!, :loc_radius][idx], locs[!, :loc_theta][idx])
    end
    # loc_x, loc_y = _locnorm(loc_x, loc_y)
    loc_x = loc_x[ch]
    loc_y = loc_y[ch]
    loc_x = _s2v(loc_x)
    loc_y = _s2v(loc_y)

    s_interpolated, interpolated_x, interpolated_y = _interpolate(signal, loc_x, loc_y, 100, imethod, nmethod)

    p = Plots.plot(grid=true,
                   framestyle=:none,
                   palette=pal,
                   size=(plot_size, plot_size),
                   border=:none,
                   aspect_ratio=1,
                   left_margin=-20 * Plots.px,
                   titlefontsize=8,
                   xlabelfontsize=6,
                   ylabelfontsize=6,
                   xtickfontsize=4,
                   ytickfontsize=4,
                   title=title;
                   kwargs...)

    p = Plots.plot!(interpolated_x,
                    interpolated_y,
                    s_interpolated,
                    fill=:darktest,
                    seriestype=:heatmap,
                    seriescolor=pal,
                    colorbar=cb,
                    colorbar_title=cb_label,
                    levels=10,
                    linewidth=0)
    if plot_contours
        p = Plots.plot!(interpolated_x,
                        interpolated_y,
                        s_interpolated,
                        fill=:darktest,
                        seriestype=:contour,
                        seriescolor=pal,
                        colorbar=cb,
                        colorbar_title=cb_label,
                        levels=5,
                        linecolor=:black,
                        linewidth=0.2)
    end
    if plot_electrodes
        p = Plots.plot!((loc_x, loc_y),
                        color=:black,
                        seriestype=:scatter,
                        grid=true,
                        label="",
                        markersize=2,
                        markeralpha=0.5,
                        markerstrokewidth=0,
                        markerstrokealpha=0)
    end

    # draw head
    hd = _draw_head(p, head_labels=head_labels, head_details=head_details, topo=true)
    p = Plots.plot!(hd)

    return p

end

"""
    plot_topo(obj; <keyword arguments>)

Topographical plot.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ep::Union{Int64, AbstractRange}=0`: epoch to display
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `seg::Tuple{Int64, Int64}=(1, 10*sr(obj))`: segment (from, to) in samples to display, default is 10 seconds or less if single epoch is shorter
- `title::String="default"`: plot title, default is Amplitude topographical plot [channels: 1:19, epoch: 1, time window: 0 ms:20 s]
- `mono::Bool=false`: use color or grey palette
- `cb::Bool=true`: plot color bar
- `cb_label::String="[A.U.]"`: color bar label
- `amethod::Symbol=:mean`: averaging method:
    - `:mean`
    - `:median`
- `imethod::Symbol=:sh`: interpolation method:
    - `:sh`: Shepard
    - `:mq`: Multiquadratic
    - `:imq`: InverseMultiquadratic
    - `:tp`: ThinPlate
    - `:nn`: NearestNeighbour
    - `:ga`: Gaussian
- `nmethod::Symbol=:minmax`: method for normalization, see `normalize()`
- `plot_size::Int64=800`: plot dimensions in pixels (size × size)
- `plot_contours::Bools=true`: plot contours over topo plot
- `plot_electrodes::Bools=true`: plot electrodes over topo plot
- `head_labels::Bool=false`: plot head labels
- `head_details::Bool=true`: draw nose and ears
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_topo(obj::NeuroAnalyzer.NEURO; ep::Union{Int64, AbstractRange}=0, ch::Union{Vector{Int64}, AbstractRange}=signal_channels(obj), seg::Tuple{Int64, Int64}=(1, 10*sr(obj)), title::String="default", mono::Bool=false, cb::Bool=true, cb_label::String="default", amethod::Symbol=:mean, imethod::Symbol=:sh, nmethod::Symbol=:minmax, plot_contours::Bool=true, plot_electrodes::Bool=true, plot_size::Int64=800, head_labels::Bool=false, head_details::Bool=true, kwargs...)

    signal_len(obj) < 10 * sr(obj) && seg == (1, 10*sr(obj)) && (seg = (1, signal_len(obj)))

    _has_locs(obj) == false && throw(ArgumentError("Electrode locations not available, use load_locs() or add_locs() first."))
    _check_var(imethod, [:sh, :mq, :imq, :tp, :nn, :ga], "imethod")
    _check_var(amethod, [:mean, :median], "amethod")
    _check_segment(obj, seg[1], seg[2])

    if ep != 0
        _check_epochs(obj, ep)
        if epoch_n(obj) == 1
            ep = 0
        else
            seg = (((ep[1] - 1) * epoch_len(obj) + 1), seg[2])
            if typeof(ep) == Int64
                seg = (seg[1], (seg[1] + epoch_len(obj) - 1))
            else
                seg = (seg[1], (ep[end] * epoch_len(obj)))
            end
            ep = 0
        end
    end

    # remove non-signal channels
    obj_tmp = deepcopy(obj)
    keep_channel_type!(obj_tmp, type=Symbol(obj_tmp.header.recording[:data_type]))

    length(ch) < 2 && throw(ArgumentError("plot_topo() requires ≥ 2 channels."))
    _check_channels(obj_tmp, ch)

    length(ch) > nrow(obj_tmp.locs) && throw(ArgumentError("Some channels do not have locations."))

    # get time vector
    if seg[2] <= epoch_len(obj_tmp)
        signal = obj_tmp.data[ch, seg[1]:seg[2], 1]
    else
        signal = ep(obj_tmp, ep_n=1).data[ch, seg[1]:seg[2], 1]
    end
    # t = _get_t(seg[1], seg[2], sr(obj_tmp))
    t = obj.time_pts[seg[1]:seg[2]]
    _, t_s1, _, t_s2 = _convert_t(t[1], t[end])
    ep = _s2epoch(obj_tmp, seg[1], seg[2])
    
    # average signal and convert to vector
    if size(signal, 2) > 1
        if amethod === :mean
            signal = vec(mean(signal, dims=2))
        elseif amethod === :median
            signal = vec(median(signal, dims=2))
        end
    else
        signal = vec(signal)
    end

    if seg[2] != seg[1] + 1
        title == "default" && (title = "Amplitude topographical plot\n[channel$(_pl(length(ch))): $(_channel2channel_name(ch)), epoch$(_pl(length(ep))): $ep, averaged ($(string(amethod))) over time window: $t_s1:$t_s2]")
    else
        title == "default" && (title = "Amplitude topographical plot\n[channel$(_pl(length(ch))): $(_channel2channel_name(ch)), epoch$(_pl(length(ep))): $ep, time point: $t_s1]")
    end
    cb_label == "default" && (cb_label = "[A.U.]")

    p = plot_topo(signal, ch=ch, locs=obj_tmp.locs, cb=cb, cb_label=cb_label, title=title, mono=mono, imethod=imethod, nmethod=nmethod, plot_contours=plot_contours, plot_electrodes=plot_electrodes, plot_size=plot_size, head_labels=head_labels, head_details=head_details, kwargs=kwargs)

    Plots.plot(p)

    return p

end

"""
    plot_topo(obj; <keyword arguments>)

Topographical plot of embedded or external component.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `c::Union{Symbol, AbstractArray}`: component to plot
- `ep::Union{Int64, AbstractRange}=0`: epoch to display
- `c_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: component channel to display, default is all component channels
- `seg::Tuple{Int64, Int64}=(1, 10*sr(obj))`: segment (from, to) in samples to display, default is 10 seconds or less if single epoch is shorter
- `title::String="default"`: plot title, default is Amplitude topographical plot [channels: 1:19, epoch: 1, time window: 0 ms:20 s]
- `mono::Bool=false`: use color or grey palette
- `cb::Bool=true`: plot color bar
- `cb_label::String="[A.U.]"`: color bar label
- `amethod::Symbol=:mean`: averaging method:
    - `:mean`
    - `:median`
- `imethod::Symbol=:sh`: interpolation method:
    - `:sh`: Shepard
    - `:mq`: Multiquadratic
    - `:imq`: InverseMultiquadratic
    - `:tp`: ThinPlate
    - `:nn`: NearestNeighbour
    - `:ga`: Gaussian
- `nmethod::Symbol=:minmax`: method for normalization, see `normalize()`
- `plot_size::Int64=800`: plot dimensions in pixels (size × size)
- `plot_contours::Bools=true`: plot contours over topo plot
- `plot_electrodes::Bools=true`: plot electrodes over topo plot
- `head_labels::Bool=false`: plot head labels
- `head_details::Bool=true`: draw nose and ears
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_topo(obj::NeuroAnalyzer.NEURO, c::Union{Symbol, AbstractArray}; ep::Union{Int64, AbstractRange}=0, c_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0, seg::Tuple{Int64, Int64}=(1, 10*sr(obj)), title::String="default", mono::Bool=false, cb::Bool=true, cb_label::String="default", amethod::Symbol=:mean, imethod::Symbol=:sh, nmethod::Symbol=:minmax, plot_contours::Bool=true, plot_electrodes::Bool=true, plot_size::Int64=800, head_labels::Bool=false, head_details::Bool=true, kwargs...)

    signal_len(obj) < 10 * sr(obj) && seg == (1, 10*sr(obj)) && (seg = (1, signal_len(obj)))

    _has_locs(obj) == false && throw(ArgumentError("Electrode locations not available, use load_locs() or add_locs() first."))
    _check_var(imethod, [:sh, :mq, :imq, :tp, :nn, :ga], "imethod")
    _check_var(amethod, [:mean, :median], "amethod")

    no_timepoint = false
    if typeof(c) == Matrix{Float64}
        c = reshape(c, size(c, 1), size(c, 2), 1)
        seg = (1, size(c, 2))
        no_timepoint = true
    elseif typeof(c) == Vector{Float64}
        c = reshape(c, length(c), 1, 1)
        seg = (1, 1)
        no_timepoint = true
    elseif size(c, 2) == 1
        no_timepoint = true
        seg = (1, 1)
    end

    _check_segment(obj, seg[1], seg[2])

    if ep != 0
        _check_epochs(obj, ep)
        if epoch_n(obj) == 1
            ep = 0
        else
            seg = (((ep[1] - 1) * epoch_len(obj) + 1), seg[2])
            if typeof(ep) == Int64
                seg = (seg[1], (seg[1] + epoch_len(obj) - 1))
            else
                seg = (seg[1], (ep[end] * epoch_len(obj)))
            end
            ep = 0
        end
    end

    # remove non-signal channels
    obj_tmp = deepcopy(obj)
    keep_channel_type!(obj_tmp, type=Symbol(obj_tmp.header.recording[:data_type]))

    # select component channels, default is all channels
    typeof(c) == Symbol && (c = _get_component(obj_tmp, c).c)
    c_idx == 0 && (c_idx = _select_cidx(c, c_idx))
    _check_cidx(c, c_idx)
    clabels = _gen_clabels(c)[c_idx]
    length(c_idx) == 1 && (clabels = [clabels])

    length(c_idx) < 2 && throw(ArgumentError("plot_topo() requires ≥ 2 channels."))
    length(ch) > nrow(obj_tmp.locs) && throw(ArgumentError("Some channels do not have locations."))

    # get time vector
    if seg[2] <= epoch_len(obj_tmp)
        signal = c[c_idx, seg[1]:seg[2], 1]
    else
        signal = _make_epochs(c, ep_n=1)[c_idx, seg[1]:seg[2], 1]
    end
    if seg[1] != seg[2]
        t = _get_t(seg[1], seg[2], sr(obj_tmp))
    else
        t = _get_t(seg[1], seg[2] + 1, sr(obj_tmp))
    end
    _, t_s1, _, t_s2 = _convert_t(t[1], t[end])
    ep = _s2epoch(obj_tmp, seg[1], seg[2])
    
    # average signal and convert to vector
    if size(signal, 2) > 1
        if amethod === :mean
            signal = vec(mean(signal, dims=2))
        elseif amethod === :median
            signal = vec(median(signal, dims=2))
        end
    else
        signal = vec(signal)
    end

    if seg[2] != seg[1]
        if no_timepoint != true
            title == "default" && (title = "Amplitude topographical plot\n[component$(_pl(length(c_idx))): $(_channel2channel_name(c_idx)), epoch$(_pl(length(ep))): $ep, averaged ($(string(amethod))) over time window: $t_s1:$t_s2]")
        else
            title == "default" && (title = "Amplitude topographical plot\n[component$(_pl(length(c_idx))): $(_channel2channel_name(c_idx)), averaged ($(string(amethod))) over $(length(seg[1]:seg[2])) time point$(_pl(length(seg)))]")
        end
    else
        if no_timepoint != true
            title == "default" && (title = "Amplitude topographical plot\n[component$(_pl(length(c_idx))): $(_channel2channel_name(c_idx)), epoch$(_pl(length(ep))): $ep, time point: $t_s1]")
        else
            title == "default" && (title = "Amplitude topographical plot\n[component$(_pl(length(c_idx))): $(_channel2channel_name(c_idx)), $(size(c, 2)) time point$(_pl(size(c, 2)))]")
        end
    end
    cb_label == "default" && (cb_label = "[A.U.]")

    p = plot_topo(signal, ch=c_idx, locs=obj_tmp.locs, cb=cb, cb_label=cb_label, title=title, mono=mono, imethod=imethod, nmethod=nmethod, plot_contours=plot_contours, plot_electrodes=plot_electrodes, plot_size=plot_size, head_labels=head_labels, head_details=head_details, kwargs=kwargs)

    Plots.plot(p)

    return p
    
end
