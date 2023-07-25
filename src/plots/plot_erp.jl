export plot_erp
export plot_erp_butterfly
export plot_erp_avg
export plot_erp_topo
export plot_erp_stack

"""
    plot_erp(t, s, bad; <keyword arguments>)

Plot ERP.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
- `s::AbstractVector`: data to plot
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `yrev::Bool=false`: reverse Y axis
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_erp(t::Union{AbstractVector, AbstractRange}, s::AbstractVector; xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, yrev::Bool=false, kwargs...)

    pal = mono == true ? :grays : :darktest

    # get limits
    ylim = (floor(minimum(s) * 1.1, digits=0), ceil(maximum(s) * 1.1, digits=0))
    ylim = _tuple_max(ylim)
    yticks = [ylim[1], 0, ylim[2]]

    # prepare plot
    p = Plots.plot(xlabel=xlabel,
                   ylabel=ylabel,
                   xlims=_xlims(t),
                   xticks=(_erpticks(t), string.(_erpticks(t) .* 1000)),
                   ylims=ylim,
                   yticks=yticks,
                   title=title,
                   palette=pal,
                   size=(1200, 400),
                   margins=20Plots.px,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=6,
                   ytickfontsize=6;
                   kwargs...)
    # reverse Y axis
    yrev == true && yflip!(true)

    # plot 0 h-line
    p = Plots.hline!([0],
                     color=:grey,
                     lw=0.5,
                     labels="")

    # plot ERP
    p = Plots.plot!(t,
                    s,
                    linewidth=1,
                    label="",
                    color=:black)

    # plot 0 v-line
    p = Plots.vline!([0],
                     linestyle=:dash,
                     linewidth=0.5,
                     linecolor=:black,
                     label=false)

    return p

end

"""
    plot_erp_butterfly(t, s; <keyword arguments>)

Butterfly plot of ERP.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
- `s::AbstractArray`: data to plot
- `clabels::Vector{String}=[""]`: signal channel labels vector
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `avg::Bool=false`: plot average ERP
- `yrev::Bool=false`: reverse Y axis
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_erp_butterfly(t::Union{AbstractVector, AbstractRange}, s::AbstractArray; clabels::Vector{String}=[""], xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, avg::Bool=true, yrev::Bool=false, kwargs...)

    pal = mono == true ? :grays : :darktest

    ch_n = size(s, 1)

    # get limits
    ylim = (floor(minimum(s) * 1.1, digits=0), ceil(maximum(s) * 1.1, digits=0))
    ylim = _tuple_max(ylim)
    yticks = [ylim[1], 0, ylim[2]]

    # plot channels
    p = Plots.plot(xlabel=xlabel,
                   ylabel=ylabel,
                   xlims=_xlims(t),
                   xticks=(_erpticks(t), string.(_erpticks(t) .* 1000)),
                   ylims=ylim,
                   yticks=yticks,
                   title=title,
                   palette=pal,
                   size=(1200, 400),
                   margins=20Plots.px,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=6,
                   ytickfontsize=6;
                   kwargs...)

    # reverse Y axis
    yrev == true && yflip!(true)

    # plot 0 h-line
    p = Plots.hline!([0],
                     color=:grey,
                     lw=0.5,
                     labels="")

    # plot signals
    for idx in 1:ch_n
        if clabels == [""]
            p = Plots.plot!(t,
                            s[idx, :],
                            t=:line,
                            linecolor=idx,
                            linewidth=0.2,
                            alpha=0.2,
                            legend=false)
        else
            if clabels == repeat([""], ch_n)
                p = Plots.plot!(t,
                                s[idx, :],
                                t=:line,
                                legend=false,
                                linecolor=idx,
                                linewidth=0.5,
                                alpha=0.5)
            else
                p = Plots.plot!(t,
                                s[idx, :],
                                t=:line,
                                label=clabels[idx],
                                linecolor=idx,
                                linewidth=0.5,
                                alpha=0.5)
            end
        end
    end

    # plot averaged ERP
    if avg == true
        if ch_n == 1
            s = mean(s, dims=2)[:]
        else
            s = mean(s, dims=1)[:]
        end
        p = Plots.plot!(t,
                        s,
                        linewidth=1,
                        linecolor=:black,
                        label=false)
    end

    # plot 0 v-line
    p = Plots.vline!([0],
                     linestyle=:dash,
                     linewidth=0.5,
                     linecolor=:black,
                     label=false)

    return p

end

"""
    plot_erp_avg(t, s; <keyword arguments>)

Plot ERP amplitude mean and ±95% CI.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
- `s::AbstractArray`: data to plot
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `yrev::Bool=false`: reverse Y axis
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_erp_avg(t::Union{AbstractVector, AbstractRange}, s::AbstractArray; xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, yrev::Bool=false, kwargs...)

    pal = mono == true ? :grays : :darktest

    # get mean and 95%CI
    s_m, _, s_u, s_l = msci95(s)

    # get limits
    ylim = (floor(minimum(s_l) * 1.1, digits=0), ceil(maximum(s_u) * 1.1, digits=0))
    ylim = _tuple_max(ylim)
    yticks = [ylim[1], 0, ylim[2]]

    # prepare plot
    p = Plots.plot(xlabel=xlabel,
                   ylabel=ylabel,
                   xlims=_xlims(t),
                   xticks=(_erpticks(t), string.(_erpticks(t) .* 1000)),
                   ylims=ylim,
                   yticks=yticks,
                   title=title,
                   palette=pal,
                   size=(1200, 400),
                   margins=20Plots.px,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=6,
                   ytickfontsize=6;
                   kwargs...)

    # reverse Y axis
    yrev == true && yflip!(true)

    # plot 0 h-line
    p = Plots.hline!([0],
                     color=:grey,
                     lw=0.5,
                     labels="")

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

    # plot 0 v-line
    p = Plots.vline!([0],
                     linestyle=:dash,
                     linewidth=0.5,
                     linecolor=:black,
                     label=false)

    return p

end

"""
    plot_erp_topo(locs, t, s; <keyword arguments>)

Plot topographical map ERPs.

# Arguments

- `locs::DataFrame`: columns: channel, labels, loc_theta, loc_radius, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
- `t::Vector{Float64}`: time vector
- `s::Array{Float64, 2}`: ERPs
- `ch::Union{Vector{Int64}, AbstractRange}`: which channels to plot
- `clabels::Vector{String}=[""]`: signal channel labels vector
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `yrev::Bool=false`: reverse Y axis
- `polar::Bool=true`: if true, use polar coordinates, otherwise use Cartesian spherical x and y coordinates
- `mono::Bool=false`: use color or grey palette

- `kwargs`: optional arguments for plot() function

# Returns

- `fig::GLMakie.Figure`
"""
function plot_erp_topo(locs::DataFrame, t::Vector{Float64}, s::Array{Float64, 2}; ch=Union{Vector{Int64}, AbstractRange}, clabels::Vector{String}=[""], xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, yrev::Bool=false, polar::Bool=true, kwargs...)

    size(s, 2) == length(t) || throw(ArgumentError("Signal length and time length must be equal."))
    length(ch) > nrow(locs) && throw(ArgumentError("Some channels do not have locations."))

    pal = mono == true ? :grays : :darktest
    
    # channel labels
    clabels == [""] && (clabels = repeat([""], size(s, 1)))

    # get limits
    ylim = (floor(minimum(s) * 1.1, digits=0), ceil(maximum(s) * 1.1, digits=0))
    ylim = _tuple_max(ylim)

    # plot parameters
    plot_size = 1200
    marker_size = (150, 75)
    
    # get locations
    if polar == true
        loc_x = zeros(nrow(locs))
        loc_y = zeros(nrow(locs))
        for idx in 1:nrow(locs)
            loc_x[idx], loc_y[idx] = pol2cart(locs[!, :loc_radius][idx], locs[!, :loc_theta][idx])
        end
        loc_x = loc_x[ch]
        loc_y = loc_y[ch]
    else
        loc_x = locs[ch, :loc_x]
        loc_y = locs[ch, :loc_y]
    end
    loc_x = _s2v(loc_x)
    loc_y = _s2v(loc_y)
    # get marker centers
    loc_x .*= ((plot_size / 2) - marker_size[1] / 2)
    loc_y .*= ((plot_size / 2) - marker_size[2] / 2)

    fig = Figure(; resolution=(plot_size, plot_size))
    fig_axis = Axis(fig[1, 1])
    fig_axis.aspect = AxisAspect(1)
    fig_axis.title = title
    GLMakie.xlims!(fig_axis, [-plot_size / 1.75, plot_size / 1.75])
    GLMakie.ylims!(fig_axis, [-plot_size / 1.75, plot_size / 1.75])
    hidedecorations!(fig_axis, grid=true, ticks=true)

    for idx in 1:size(s, 1)
        p = Plots.plot(xlabel=xlabel,
                       ylabel=ylabel,
                       legend=false,
                       xticks=false,
                       yticks=false,
                       grid=false,
                       border=:none, 
                       xlims=(t[1], t[end]),
                       ylims=ylim,
                       title=clabels[idx],
                       palette=pal,
                       size=marker_size,
                       titlefontsize=8,
                       xlabelfontsize=8,
                       ylabelfontsize=8,
                       xtickfontsize=6,
                       ytickfontsize=6;
                       kwargs...)
        # reverse Y axis
        yrev == true && yflip!(true)

        # plot 0 h-line
        p = Plots.hline!([0],
                         color=:grey,
                         linewidth=0.5,
                         labels="")

        # plot ERP
        p = Plots.plot!(t,
                        s[idx, :],
                        t=:line,
                        color=:black,
                        linewidth=0.5,
                        alpha=0.75)

        # plot 0 v-line
        p = Plots.vline!([0],
                         linestyle=:dash,
                         linewidth=0.5,
                         linecolor=:black,
                         label=false)

        marker_img = tempname() * ".png"
        savefig(p, marker_img)
        marker = FileIO.load(marker_img)
        GLMakie.scatter!(fig_axis, (loc_x[idx], loc_y[idx]), marker=marker, markersize=marker_size)
        rm(marker_img)
    end

    return fig

end

"""
    plot_erp_stack(s; <keyword arguments>)

Plot EPRs stacked by channels or by epochs.

# Arguments

- `t::AbstractVector`: x-axis values
- `s::AbstractArray`
- `clabels::Vector{String}=[""]`: signal channel labels vector
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `cb::Bool=true`: plot color bar
- `cb_title::String=""`: color bar title
- `mono::Bool=false`: use color or grey palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_erp_stack(t::AbstractVector, s::AbstractArray; clabels::Vector{String}=[""], xlabel::String="", ylabel::String="", title::String="", cb::Bool=true, cb_title::String="", mono::Bool=false, kwargs...)

    ndims(s) == 2 || throw(ArgumentError("signal must have 2 dimensions."))
    length(t) == size(s, 2) || throw(ArgumentError("Number of signal columns ($(size(s, 2))) must be equal to length of x-axis values ($(length(t)))."))

    pal = mono == true ? :grays : :darktest

    if clabels == [""]
        yticks = round.(Int64, range(1, size(s, 1), length=10))
    else
        yticks = (1:size(s, 1), clabels)
    end
    p = Plots.heatmap(t,
                      1:size(s, 1),
                      s,
                      size=(1200, 500),
                      margins=20Plots.px,
                      legend=false,
                      xticks=(_erpticks(t), string.(_erpticks(t) .* 1000)),
                      yticks=yticks,
                      xlabel=xlabel,
                      ylabel=ylabel,
                      cb=cb,
                      cbtitle=cb_title,
                      title=title,
                      seriescolor=pal,
                      linewidth=0.5,
                      titlefontsize=8,
                      xlabelfontsize=8,
                      ylabelfontsize=8,
                      xtickfontsize=8,
                      ytickfontsize=8;
                      kwargs...)

    # plot 0 v-line
    p = Plots.vline!([0],
                     linestyle=:dash,
                     linewidth=0.5,
                     linecolor=:black,
                     label=false)

    Plots.plot(p)

    return p

end

"""
    plot_erp(obj; <keyword arguments>)

Plot ERP.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}`: channel(s) to plot
- `tm::Union{Int64, Vector{Int64}}=0`: time markers (in miliseconds) to plot as vertical lines, useful for adding topoplots at these time points 
- `xlabel::String="default"`: x-axis label, default is Time [ms]
- `ylabel::String="default"`: y-axis label, default is Amplitude [units] 
- `title::String="default"`: plot title, default is ERP amplitude [channel: 1, epochs: 1:2, time window: -0.5 s:1.5 s]
- `cb::Bool=true`: plot color bar
- `cb_title::String="default"`: color bar title, default is Amplitude [units] 
- `mono::Bool=false`: use color or grey palette
- `peaks::Bool=true`: draw peaks
- `channel_labels::Bool=true`: draw labels legend (using channel labels) for multi-channel `:butterfly` plot
- `type::Symbol=:normal`: plot type: `:normal`, butterfly plot (`:butterfly`), topographical plot of ERPs (`:topo`) or stacked epochs/channels (`:stack`)
- `yrev::Bool=false`: reverse Y axis
- `avg::Bool=false`: plot average ERP for `:butterfly` plot
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_erp(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}, tm::Union{Int64, Vector{Int64}}=0, xlabel::String="default", ylabel::String="default", title::String="default", cb::Bool=true, cb_title::String="default", mono::Bool=false, peaks::Bool=true, channel_labels::Bool=true, type::Symbol=:normal, yrev::Bool=false, avg::Bool=true, kwargs...)

    _check_datatype(obj, :erp)

    # check channels
    _check_channels(obj, ch)

    # set units
    units = _set_units(obj, ch[1])

    _check_var(type, [:normal, :butterfly, :mean, :topo, :stack], "type")
    (length(ch) > 1 && length(unique(obj.header.recording[:channel_type][ch])) > 1) && throw(ArgumentError("All channels must be of the same type."))

    type in [:normal] && length(ch) > 1 && throw(ArgumentError("For :normal plot type, only one channel must be specified."))
    # type in [:butterfly, :stack] && length(ch) < 2 && throw(ArgumentError("For :butterfly and :stack plot type ≥ 2 channels must be specified."))

    # get data
    # ch=1
    # ch=1:10

    ep_n = epoch_n(obj) - 1

    if type in [:normal, :topo]
        s = obj.data[ch, :, 1]
    else
        if ch isa Int64
            s = obj.data[ch, :, 2:end]'
            channel_labels = false
            # s = reshape(obj.data[ch, :, 2:end], 1, :, ep_n)
        else
            s = obj.data[ch, :, 1]
        end
    end

    # get time vector
    t = obj.epoch_time
    _, t_s1, _, t_s2 = _convert_t(t[1], t[end])

    if tm != 0
        for tm_idx in eachindex(tm)
            tm[tm_idx] / 1000 < t[1] && throw(ArgumentError("tm value ($(tm[tm_idx])) is out of epoch time segment ($(t[1]):$(t[end]))."))
            tm[tm_idx] / 1000 > t[end] && throw(ArgumentError("tm value ($(tm[tm_idx])) is out of epoch time segment ($(t[1]):$(t[end]))."))
            tm[tm_idx] = vsearch(tm[tm_idx] / 1000, t)
        end
    end

    if type === :normal
        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [ms]", "Amplitude [$units]", "ERP amplitude channel $(_channel2channel_name(ch))\n[averaged epochs: $ep_n, time window: $t_s1:$t_s2]")
        p = plot_erp(t,
                     s,
                     xlabel=xl,
                     ylabel=yl,
                     title=tt,
                     mono=mono,
                     yrev=yrev;
                     kwargs...)
    elseif type === :butterfly
        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [ms]", "Amplitude [$units]", "ERP amplitude channel$(_pl(length(ch))) $(_channel2channel_name(ch))\n[averaged epochs: $ep_n, time window: $t_s1:$t_s2]")
        if channel_labels == true
            clabels = labels(obj)[ch]
        else
            clabels = repeat([""], length(ch))
        end
        p = plot_erp_butterfly(t,
                               s,
                               xlabel=xl,
                               ylabel=yl,
                               title=tt,
                               clabels=clabels,
                               mono=mono,
                               avg=avg,
                               yrev=yrev;
                               kwargs...)
    elseif type === :mean
        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [ms]", "Amplitude [$units]", "ERP amplitude [mean ± 95%CI] channel$(_pl(length(ch))) $(_channel2channel_name(ch))\n[averaged epochs: $ep_n, time window: $t_s1:$t_s2]")
        p = plot_erp_avg(t,
                         s,
                         xlabel=xl,
                         ylabel=yl,
                         title=tt,
                         mono=mono,
                         yrev=yrev;
                         kwargs...)
    elseif type === :topo
        _has_locs(obj) == false && throw(ArgumentError("Electrode locations not available."))
        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "", "", "ERP amplitude channel$(_pl(length(ch))) $(_channel2channel_name(ch))\n[averaged epochs: $ep_n, time window: $t_s1:$t_s2]")
        peaks = false
        ndims(s) == 1 && (s = reshape(s, 1, length(s)))
        clabels = labels(obj)[ch]
        clabels isa String && (clabels = [clabels])
        p = plot_erp_topo(obj.locs,
                          t,
                          s,
                          ch=ch,
                          clabels=clabels,
                          xlabel=xl,
                          ylabel=yl,
                          title=tt,
                          mono=mono,
                          yrev=yrev;
                          kwargs...)
    elseif type === :stack
        peaks = false
        cb_title == "default" && (cb_title = "Amplitude [$units]")
        if ch isa Int64
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [ms]", "Epochs", "ERP amplitude channel$(_pl(length(ch))) $(_channel2channel_name(ch))\n[averaged epochs: $ep_n, time window: $t_s1:$t_s2]")
        else
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [ms]", "Channels", "ERP amplitude channel$(_pl(length(ch))) $(_channel2channel_name(ch))\n[averaged epochs: $ep_n, time window: $t_s1:$t_s2]")
        end
        if channel_labels == true
            clabels = labels(obj)[ch]
        else
            clabels = repeat([""], length(ch))
        end
        p = plot_erp_stack(t,
                           s,
                           xlabel=xl,
                           ylabel=yl,
                           title=tt,
                           clabels=clabels,
                           cb=cb,
                           cb_title=cb_title,
                           mono=mono;
                           kwargs...)
    end

    # draw time markers
    if tm != 0
        for tm_idx in eachindex(tm)
            p = Plots.vline!([t[tm[tm_idx]]],
                             linewidth=1,
                             linecolor=:black,
                             label=false)
        end
    end

    # draw peaks
    if peaks == true && avg == true
        if ch isa Int64
            pp = erp_peaks(obj)
            if mono == false
                Plots.scatter!((t[pp[ch, 1]], obj.data[ch, pp[ch, 1]]), marker=:xcross, markercolor=:red, markersize=3, label=false)
                Plots.scatter!((t[pp[ch, 2]], obj.data[ch, pp[ch, 2]]), marker=:xcross, markercolor=:blue, markersize=3, label=false)
            else
                Plots.scatter!((t[pp[ch, 1]], obj.data[ch, pp[ch, 1]]), marker=:xcross, markercolor=:black, markersize=3, label=false)
                Plots.scatter!((t[pp[ch, 2]], obj.data[ch, pp[ch, 2]]), marker=:xcross, markercolor=:black, markersize=3, label=false)
            end
            _info("Positive peak time: $(round(t[pp[ch, 1]] * 1000, digits=0)) ms")
            _info("Positive peak amplitude: $(round(obj.data[ch, pp[ch, 1]], digits=2)) $units")
            _info("Negative peak time: $(round(t[pp[ch, 2]] * 1000, digits=0)) ms")
            _info("Negative peak amplitude: $(round(obj.data[ch, pp[ch, 2]], digits=2)) $units")
        else
            erp_tmp = mean(mean(obj.data[ch, :, 2:end], dims=1), dims=3)
            obj_tmp = keep_channel(obj, ch=1)
            obj_tmp.data = erp_tmp
            pp = erp_peaks(obj_tmp)
            if mono == false
                Plots.scatter!((t[pp[1, 1]], erp_tmp[pp[1, 1]]), marker=:xcross, markercolor=:red, markersize=3, label=false)
                Plots.scatter!((t[pp[1, 2]], erp_tmp[pp[1, 2]]), marker=:xcross, markercolor=:blue, markersize=3, label=false)
            else
                Plots.scatter!((t[pp[1, 1]], erp_tmp[pp[1, 1]]), marker=:xcross, markercolor=:black, markersize=3, label=false)
                Plots.scatter!((t[pp[1, 2]], erp_tmp[pp[1, 2]]), marker=:xcross, markercolor=:black, markersize=3, label=false)
            end
            _info("Positive peak time: $(round(t[pp[1, 1]] * 1000, digits=0)) ms")
            _info("Positive peak amplitude: $(round(erp_tmp[pp[1, 1]], digits=2)) $units")
            _info("Negative peak time: $(round(t[pp[1, 2]] * 1000, digits=0)) ms")
            _info("Negative peak amplitude: $(round(erp_tmp[pp[1, 2]], digits=2)) $units")
        end
    end

    if type !== :topo
        Plots.plot(p)
    else
        GLMakie.show(p)
    end

    return p

end
