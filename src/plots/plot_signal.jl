export plot_signal
export plot_signal_avg
export plot_signal_butterfly
export plot

"""
    plot_signal(t, s; <keyword arguments>)

Plot amplitude of single- or multi-channel `s`.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
- `s::Union{AbstractVector, AbstractArray}`: data to plot
- `clabels::Vector{String}=[""]`: signal channel labels vector
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `scale::Bool=true`: draw scale
- `units::String="μV"`: units of the scale
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_signal(t::Union{AbstractVector, AbstractRange}, s::Union{AbstractVector, AbstractArray}; clabels::Vector{String}=[""], xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, scale::Bool=true, units::String="μV", kwargs...)

    # convert single-channel signal to single-row matrix
    ndims(s) == 1 && (s = reshape(s, 1, length(s)))
    ch_n = size(s, 1)

    # reverse so 1st channel is on top
    s = @views reverse(s[:, 1:length(t)], dims = 1)
    # also, reverse colors if palette is not mono
    if mono == true
        pal = :grays
        channel_color = Vector{Symbol}()
        for idx in 1:ch_n
            push!(channel_color, :black)
        end
    else
        pal = :darktest
        channel_color = ch_n:-1:1
    end

    # get range of the original signal for the scale
    range = _get_range(s)

    # normalize and shift so all channels are visible
    # each channel is between -1.0 and +1.0
    for idx in 1:ch_n
        # scale by 0.5 so maxima do not overlap
        s[idx, :] = @views normalize(s[idx, :], method=:minmax) .* 0.5 .+ (idx - 1)
    end

    # prepare plot
    ch_n == 1 && (plot_size = (1200, 500))
    ch_n > 1 && (plot_size = (1200, 800))
    p = Plots.plot(xlabel=xlabel,
                   ylabel=ylabel,
                   xlims=_xlims(t),
                   xticks=_ticks(t),
                   ylims=(-1, ch_n),
                   title=title,
                   palette=pal,
                   size=plot_size,
                   margins=20Plots.px,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=6,
                   ytickfontsize=6;
                   kwargs...)

    # plot zero line
    p = Plots.hline!(collect((ch_n - 1):-1:0),
                     color=:grey,
                     lw=0.5,
                     labels="")

    # plot channels
    for idx in 1:ch_n
        p = @views Plots.plot!(t,
                               s[idx, :],
                               linewidth=1,
                               label="",
                               color=channel_color[idx])
    end

    # plot labels
    p = Plots.plot!(yticks=((ch_n - 1):-1:0, clabels))

    # draw scale
    if scale == true
        p = Plots.plot!([t[1], t[1]], [(ch_n - 1.5), (ch_n - 0.5)], color=:red, linewidth=5, label="")
        p = Plots.plot!(annotation=(t[1], (ch_n - 1), Plots.text("$range$units", pointsize=6, halign=:center, valign=:bottom, rotation=90)), label=false)
    end

    return p

end

"""
    plot_signal(t, s, bad; <keyword arguments>)

Plot amplitude of single- or multi-channel `s`.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
- `s::Union{AbstractVector, AbstractArray}`: data to plot
- `norm::Bool=false`: normalize signal for butterfly and averaged plots
- `bad::Vector{Bool}}`: list of bad channels
- `clabels::Vector{String}=[""]`: signal channel labels vector
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `scale::Bool=true`: draw scale
- `units::String="μV"`: units of the scale
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_signal(t::Union{AbstractVector, AbstractRange}, s::Union{AbstractVector, AbstractArray}, bad::Vector{Bool}; clabels::Vector{String}=[""], xlabel::String="", ylabel::String="", title::String="", scale::Bool=true, units::String="μV", kwargs...)

    length(bad) == size(s, 1) || throw(ArgumentError("Length of bad channels vector and number of channels must be equal."))

    # convert single-channel signal to single-row matrix
    ndims(s) == 1 && (s = reshape(s, 1, length(s)))
    ch_n = size(s, 1)

    # reverse so 1st channel is on top
    s = @views reverse(s[:, 1:length(t)], dims = 1)
    bad = reverse(bad)

    pal = mono == true ? :grays : :darktest

    # get range of the original signal for the scale
    range = _get_range(s)

    # normalize and shift so all channels are visible
    # each channel is between -1.0 and +1.0
    for idx in 1:ch_n
        # scale by 0.5 so maxima do not overlap
        s[idx, :] = @views normalize(s[idx, :], method=:minmax) .* 0.5 .+ (idx - 1)
    end

    # prepare plot
    ch_n == 1 && (plot_size = (1200, 500))
    ch_n > 1 && (plot_size = (1200, 800))
    p = Plots.plot(xlabel=xlabel,
                   ylabel=ylabel,
                   xlims=_xlims(t),
                   xticks=_ticks(t),
                   ylims=(-1, ch_n),
                   title=title,
                   palette=pal,
                   size=plot_size,
                   margins=20Plots.px,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=6,
                   ytickfontsize=6;
                   kwargs...)
    
    # plot zero line
    p = Plots.hline!(collect((ch_n - 1):-1:0),
                     color=:grey,
                     lw=0.5,
                     labels="")

    # plot channels
    for idx in 1:ch_n
        if bad[idx] == true
            p = @views Plots.plot!(t,
                                   s[idx, :],
                                   linewidth=1,
                                   label="",
                                   color=:red)
        else
            p = @views Plots.plot!(t,
                                   s[idx, :],
                                   linewidth=1,
                                   label="",
                                   color=:black)
        end
    end

    # plot labels
    p = Plots.plot!(yticks=((ch_n - 1):-1:0, clabels))

    # draw scale
    if scale == true
        p = Plots.plot!([t[1], t[1]], [(ch_n - 1.5), (ch_n - 0.5)], color=:red, linewidth=5, label="")
        p = Plots.plot!(annotation=(t[1], (ch_n - 1), Plots.text("$range$units", pointsize=6, halign=:center, valign=:bottom, rotation=90)), label=false)
    end

    return p

end

"""
    plot_signal_avg(t, signal; <keyword arguments>)

Plot amplitude mean and ±95% CI of averaged `signal` channels.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
- `s::AbstractArray`: data to plot
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `scale::Bool=true`: draw scale
- `units::String="μV"`: units of the scale
- `norm::Bool=false`: normalize to -1 .. +1
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_signal_avg(t::Union{AbstractVector, AbstractRange}, s::AbstractArray; xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, scale::Bool=true, units::String="μV", norm::Bool=false, kwargs...)

    pal = mono == true ? :grays : :darktest

    # get range of the original signal for the scale
    range = _get_range(s)

    # get mean and 95%CI
    s_m, _, s_u, s_l = msci95(s)

    # get limits
    if norm != true
        ylim = (floor(minimum(s_l), digits=0), ceil(maximum(s_u), digits=0))
        ylim = _tuple_max(ylim)
        yticks = [ylim[1], 0, ylim[2]]
    else
        s_m = normalize(s_m, method=:minmax)
        s_u = normalize(s_u, method=:minmax)
        s_l = normalize(s_l, method=:minmax)
        ylim = (-1.0, 1.0)
        yticks = [0]
    end

    # prepare plot
    p = Plots.plot(xlabel=xlabel,
                   ylabel=ylabel,
                   xlims=_xlims(t),
                   xticks=_ticks(t),
                   ylims=ylim,
                   yticks=yticks,
                   title=title,
                   palette=pal,
                   size=(1200, 500),
                   margins=20Plots.px,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=6,
                   ytickfontsize=6;
                   kwargs...)
    # plot upper 95% CI
    p = Plots.plot!(t,
                    s_u,
                    fillrange=s_l,
                    fillalpha=0.35, 
                    label=false,
                    t=:line,
                    c=:grey,
                    lw=0.5)
    # plot lower 95% CI
    p = Plots.plot!(t,
                    s_l,
                    label=false,
                    t=:line,
                    c=:grey,
                    lw=0.5)
    # plot mean
    p = Plots.plot!(t,
                    s_m,
                    label=false,
                    t=:line,
                    c=:black,
                    lw=0.5)

    # draw scale
    if norm == true && scale == true
        p = Plots.plot!([t[1], t[1]], [-1, 1], color=:red, linewidth=5, label=false)
        p = Plots.plot!(annotation=(t[1], 0, Plots.text("$range$units", pointsize=6, halign=:center, valign=:bottom, rotation=90)), label=false)
    end

    return p

end

"""
    plot_signal_butterfly(t, s; <keyword arguments>)

Butterfly plot of `s` channels.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
- `s::AbstractArray`: data to plot
- `clabels::Vector{String}=[""]`: signal channel labels vector
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `scale::Bool=true`: draw scale
- `units::String="μV"`: units of the scale
- `mono::Bool=false`: use color or grey palette
- `norm::Bool=false`: normalize to -1 .. +1
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_signal_butterfly(t::Union{AbstractVector, AbstractRange}, s::AbstractArray; clabels::Vector{String}=[""], xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, scale::Bool=true, units::String="μV", norm::Bool=false, kwargs...)

    pal = mono == true ? :grays : :darktest

    # get range of the original signal for the scale
    range = _get_range(s)

    ch_n = size(s, 1)

    # get limits
    if norm != true
        ylim = (floor(minimum(s), digits=0), ceil(maximum(s), digits=0))
        ylim = _tuple_max(ylim)
        yticks = [ylim[1], 0, ylim[2]]
    else
        s = normalize(s, method=:minmax)
        ylim = (-1.0, 1.0)
        yticks = [0]
    end

    # channel labels
    clabels == [""] && (clabels = repeat([""], ch_n))
    
    # plot channels
    p = Plots.plot(xlabel=xlabel,
                   ylabel=ylabel,
                   xlims=_xlims(t),
                   xticks=_ticks(t),
                   ylims=ylim,
                   yticks=yticks,
                   title=title,
                   palette=pal,
                   size=(1200, 500),
                   margins=20Plots.px,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=6,
                   ytickfontsize=6;
                   kwargs...)
    for idx in 1:ch_n
        p = Plots.plot!(t,
                        s[idx, :],
                        t=:line,
                        linecolor=idx,
                        linewidth=0.5,
                        label=clabels[idx],
                        legend=true)
    end
    
    # draw scale
    if norm == true && scale == true
        p = Plots.plot!([t[1], t[1]], [-1, 1], color=:red, linewidth=5, label=false)
        p = Plots.plot!(annotation=(t[1], 0, Plots.text("$range$units", pointsize=6, halign=:center, valign=:bottom, rotation=90)), label=false)
    end

    return p

end

"""
    plot(obj; <keyword arguments>)

Plot signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ep::Union{Int64, AbstractRange}=0`: epoch to display
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj))`: channel(s) to plot, default is all channels
- `seg::Tuple{Int64, Int64}=(1, 10*sr(obj))`: segment (from, to) in samples to display, default is 10 seconds or less if single epoch is shorter
- `xlabel::String="default"`: x-axis label, default is Time [s]
- `ylabel::String="default"`: y-axis label, default is no label
- `title::String="default"`: plot title, default is Amplitude [channels: 1:2, epochs: 1:2, time window: 0 ms:20 s]
- `mono::Bool=false`: use color or grey palette
- `emarkers::Bool`: draw epoch markers if available
- `markers::Bool`: draw markers if available
- `scale::Bool=true`: draw scale
- `units::String="μV"`: units of the scale
- `type::Symbol=:normal`: plot type: `:normal`, mean ± 95%CI (`:mean`), butterfly plot (`:butterfly`)
- `norm::Bool=false`: normalize signal for butterfly and averaged plots
- `bad::Union{Bool, Matrix{Bool}}=false`: list of bad channels; if not empty - plot bad channels using this list
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot(obj::NeuroAnalyzer.NEURO; ep::Union{Int64, AbstractRange}=0, ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)), seg::Tuple{Int64, Int64}=(1, 10*sr(obj)), xlabel::String="default", ylabel::String="default", title::String="default", mono::Bool=false, emarkers::Bool=true, markers::Bool=true, scale::Bool=true, units::String="μV", type::Symbol=:normal, norm::Bool=false, bad::Union{Bool, Matrix{Bool}}=false, kwargs...)

    signal_len(obj) < 10 * sr(obj) && seg == (1, 10*sr(obj)) && (seg = (1, signal_len(obj)))

    _check_var(type, [:normal, :butterfly, :mean], "type")
    _check_segment(obj, seg[1], seg[2])

    if ep != 0
        _check_epochs(obj, ep)
        seg = (((ep[1] - 1) * epoch_len(obj) + 1), seg[2])
        if typeof(ep) == Int64
            seg = (seg[1], (seg[1] + epoch_len(obj) - 1))
        else
            seg = (seg[1], (ep[end] * epoch_len(obj)))
        end
    end

    # do not show epoch markers if there are no epochs
    epoch_n(obj) == 1 && (emarkers = false)
    if emarkers == true
        epoch_markers = _get_epoch_markers(obj)
    end

    # check channels
    _check_channels(obj, ch)
    clabels = labels(obj)[ch]
    length(ch) == 1 && (clabels = [clabels])

    # get time vector
    if seg[2] <= epoch_len(obj)
        s = obj.data[ch, seg[1]:seg[2], 1]
    else
        s = epoch(obj, ep_n=1).data[ch, seg[1]:seg[2], 1]
    end
    #t = _get_t(seg[1], seg[2], sr(obj))
    t = obj.time_pts[seg[1]:seg[2]]

    _, t_s1, _, t_s2 = _convert_t(t[1], t[end])
    ep = _s2epoch(obj, seg[1], seg[2])

    if type === :normal
        if bad == false
            xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "Channel$(_pl(length(ch))) $(_channel2channel_name(ch)) amplitude\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
            p = plot_signal(t,
                            s,
                            clabels=clabels,
                            xlabel=xlabel,
                            ylabel=ylabel,
                            title=title,
                            scale=scale,
                            units=units,
                            mono=mono;
                            kwargs...)
        else
            xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "Bad channel$(_pl(length(ch))) $(_channel2channel_name(ch))\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
            length(ch) > size(bad, 1) && throw(ArgumentError("Number of channels cannot be larger than number of bad channels rows."))
            ep > size(bad, 2) && throw(ArgumentError("Epoch number cannot be larger than number of bad channels columns."))
            p = plot_signal(t,
                            s,
                            bad[ch, ep],
                            clabels=clabels,
                            xlabel=xlabel,
                            ylabel=ylabel,
                            title=title,
                            scale=scale,
                            units=units;
                            kwargs...)
        end
    elseif type === :butterfly
        size(s, 1) == 1 && throw(ArgumentError("For type=:butterfly plot the signal must contain ≥ 2 channels."))
        xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [s]", "Amplitude [μV]", "Channels $(_channel2channel_name(ch)) amplitude\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        p = plot_signal_butterfly(t,
                                  s,
                                  clabels=clabels,
                                  xlabel=xlabel,
                                  ylabel=ylabel,
                                  title=title,
                                  scale=scale,
                                  units=units,
                                  norm=norm,
                                  mono=mono;
                                  kwargs...)
    elseif type === :mean
        size(s, 1) == 1 && throw(ArgumentError("For type=:mean plot the signal must contain ≥ 2 channels."))
        xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [s]", "Amplitude [μV]", "Averaged channels $(_channel2channel_name(ch)) amplitude [mean ± 95%CI]\n [epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        p = plot_signal_avg(t,
                            s,
                            xlabel=xlabel,
                            ylabel=ylabel,
                            title=title,
                            scale=scale,
                            units=units,
                            norm=norm,
                            mono=mono;
                            kwargs...)
    end

    # add epochs markers
    # TODO: draw epoch numbers
    if emarkers == true
        p = Plots.vline!(epoch_markers,
                         linestyle=:dash,
                         linewidth=0.5,
                         linecolor=:blue,
                         label="")
    end

    # plot markers if available
    # TODO: draw markers length
    if markers == true && _has_markers(obj) == true
        markers_pos = obj.markers[!, :start] ./ sr(obj)
        markers_desc = obj.markers[!, :description]
        p = Plots.vline!(markers_pos,
                         linestyle=:dash,
                         linewidth=0.5,
                         linecolor=:black,
                         label=false)
        for idx in eachindex(markers_desc)
            p = Plots.plot!(annotation=(markers_pos[idx], -0.92, Plots.text("$(markers_desc[idx])", pointsize=5, halign=:left, valign=:top, rotation=90)), label=false)
        end
    end

    Plots.plot(p)

    return p

end

"""
    plot(obj, c; <keyword arguments>)

Plot embedded or external component.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `c::Union{Symbol, AbstractArray}`: component to plot
- `ep::Union{Int64, AbstractRange}=0`: epoch to display
- `c_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0`: component channel to display, default is all component channels
- `seg::Tuple{Int64, Int64}=(1, 10*sr(obj))`: segment (from, to) in samples to display, default is 10 seconds or less if single epoch is shorter
- `xlabel::String="default"`: x-axis label, default is Time [s]
- `ylabel::String="default"`: y-axis label, default is no label
- `title::String="default"`: plot title, default is Amplitude [channels: 1:2, epochs: 1:2, time window: 0 ms:20 s]
- `mono::Bool=false`: use color or grey palette
- `emarkers::Bool`: draw epoch markers if available
- `markers::Bool`: draw markers if available
- `scale::Bool=true`: draw scale
- `units::String=""`: units of the scale
- `type::Symbol=:normal`: plot type: `:normal`, mean ± 95%CI (`:mean`), butterfly plot (`:butterfly`)
- `norm::Bool=false`: normalize signal for butterfly and averaged plots
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot(obj::NeuroAnalyzer.NEURO, c::Union{Symbol, AbstractArray}; ep::Union{Int64, AbstractRange}=0, c_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0, seg::Tuple{Int64, Int64}=(1, 10*sr(obj)), xlabel::String="default", ylabel::String="default", title::String="default", mono::Bool=false, emarkers::Bool=true, markers::Bool=true, scale::Bool=true, units::String="", type::Symbol=:normal, norm::Bool=false, kwargs...)

    signal_len(obj) < 10 * sr(obj) && seg == (1, 10*sr(obj)) && (seg = (1, signal_len(obj)))

    _check_var(type, [:normal, :butterfly, :mean], "type")
    _check_segment(obj, seg[1], seg[2])

    if ep != 0
        _check_epochs(obj, ep)
        seg = (((ep[1] - 1) * epoch_len(obj) + 1), seg[2])
        if typeof(ep) == Int64
            seg = (seg[1], (seg[1] + epoch_len(obj) - 1))
        else
            seg = (seg[1], (ep[end] * epoch_len(obj)))
        end
    end

    # do not show epoch markers if there are no epochs
    epoch_n(obj) == 1 && (emarkers = false)
    if emarkers == true
        epoch_markers = _get_epoch_markers(obj)
    end

    # select component channels, default is all channels
    typeof(c) == Symbol && (c = _get_component(obj, c).c)
    c_idx == 0 && (c_idx = _select_cidx(c, c_idx))
    _check_cidx(c, c_idx)
    clabels = _gen_clabels(c)[c_idx]
    length(c_idx) == 1 && (clabels = [clabels])

    # get time vector
    if seg[2] <= epoch_len(obj)
        s = c[c_idx, seg[1]:seg[2], 1]
    else
        s = _make_epochs(c, ep_n=1)[c_idx, seg[1]:seg[2], 1]
    end
    t = _get_t(seg[1], seg[2], sr(obj))

    _, t_s1, _, t_s2 = _convert_t(t[1], t[end])
    ep = _s2epoch(obj, seg[1], seg[2])

    if type === :normal
        xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "Component$(_pl(length(c_idx))) $(_channel2channel_name(c_idx)) amplitude\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        p = plot_signal(t,
                        s,
                        clabels=clabels,
                        xlabel=xlabel,
                        ylabel=ylabel,
                        title=title,
                        scale=scale,
                        units=units,
                        mono=mono;
                        kwargs...)
    elseif type === :butterfly
        size(s, 1) == 1 && throw(ArgumentError("For type=:butterfly plot the signal must contain ≥ 2 channels."))
        xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [s]", "Amplitude [μV]", "Components $(_channel2channel_name(c_idx)) amplitude\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        p = plot_signal_butterfly(t,
                                  s,
                                  clabels=clabels,
                                  xlabel=xlabel,
                                  ylabel=ylabel,
                                  title=title,
                                  scale=scale,
                                  units=units,
                                  norm=norm,
                                  mono=mono;
                                  kwargs...)
    elseif type === :mean
        size(s, 1) == 1 && throw(ArgumentError("For type=:mean plot the signal must contain ≥ 2 channels."))
        xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [s]", "Amplitude [μV]", "Averaged components $(_channel2channel_name(c_idx)) amplitude [mean ± 95%CI]\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        p = plot_signal_avg(t,
                            s,
                            clabels=clabels,
                            xlabel=xlabel,
                            ylabel=ylabel,
                            title=title,
                            scale=scale,
                            units=units,
                            norm=norm,
                            mono=mono;
                            kwargs...)
    end

    # add epochs markers
    # TODO: draw epoch numbers
    if emarkers == true
        p = Plots.vline!(epoch_markers,
                         linestyle=:dash,
                         linewidth=0.5,
                         linecolor=:blue,
                         label="")
    end

    # plot markers if available
    # TODO: draw markers length
    if markers == true && _has_markers(obj) == true
        markers_pos = obj.markers[!, :start] ./ sr(obj)
        markers_desc = obj.markers[!, :description]
        p = Plots.vline!(markers_pos,
                         linestyle=:dash,
                         linewidth=0.5,
                         linecolor=:black,
                         label=false)
        for idx in eachindex(markers_desc)
            p = Plots.plot!(annotation=(markers_pos[idx], -0.92, Plots.text("$(markers_desc[idx])", pointsize=5, halign=:left, valign=:top, rotation=90)), label=false)
        end
    end

    Plots.plot(p)

    return p
    
end
