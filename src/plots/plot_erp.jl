export plot_erp
export plot_erp_avg
export plot_erp_butterfly
export plot_erp_topo
export plot_erp_stack

"""
    plot_erp(t, signal, bad; <keyword arguments>)

Plot ERP.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
- `signal::AbstractVector`: data to plot
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `yrev::Bool=false`: reverse Y axis
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_erp(t::Union{AbstractVector, AbstractRange}, signal::AbstractVector; xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, yrev::Bool=false, kwargs...)

    pal = mono == true ? :grays : :darktest

    # get limits
    ylim = (floor(minimum(signal) * 1.1, digits=0), ceil(maximum(signal) * 1.1, digits=0))
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
                    signal,
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
    plot_erp_avg(t, signal; <keyword arguments>)

Plot ERP amplitude mean and ±95% CI.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
- `signal::AbstractArray`: data to plot
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or grey palette
- `yrev::Bool=false`: reverse Y axis
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_erp_avg(t::Union{AbstractVector, AbstractRange}, signal::AbstractArray; xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, yrev::Bool=false, kwargs...)

    pal = mono == true ? :grays : :darktest

    # get mean and 95%CI
    s_m, _, s_u, s_l = s_msci95(signal')

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
    plot_erp_butterfly(t, signal; <keyword arguments>)

Butterfly plot of ERP.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
- `signal::AbstractArray`: data to plot
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
function plot_erp_butterfly(t::Union{AbstractVector, AbstractRange}, signal::AbstractArray; clabels::Vector{String}=[""], xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, avg::Bool=true, yrev::Bool=false, kwargs...)

    pal = mono == true ? :grays : :darktest

    ch_n = size(signal, 1)

    # get limits
    ylim = (floor(minimum(signal) * 1.1, digits=0), ceil(maximum(signal) * 1.1, digits=0))
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
                            signal[idx, :],
                            t=:line,
                            linecolor=idx,
                            linewidth=0.2,
                            alpha=0.2,
                            legend=false)
        else
            if clabels == repeat([""], ch_n)
                p = Plots.plot!(t,
                                signal[idx, :],
                                t=:line,
                                legend=false,
                                linecolor=idx,
                                linewidth=0.5,
                                alpha=0.5)
            else
                p = Plots.plot!(t,
                                signal[idx, :],
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
            signal = mean(signal, dims=2)[:]
        else
            signal = mean(signal, dims=1)[:]
        end
        p = Plots.plot!(t,
                        signal,
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
    plot_erp_topo(locs, t, erp; <keyword arguments>)

Plot topographical map ERPs. It uses polar :loc_radius and :loc_theta locations, which are translated into Cartesian x and y positions.

# Arguments

- `locs::DataFrame`: columns: channel, labels, loc_theta, loc_radius, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
- `t::Vector{Float64}`: time vector
- `signal::Array{Float64, 2}`: ERPs
- `channels::Union{Vector{Int64}, AbstractRange}`: which channels to plot
- `clabels::Vector{String}=[""]`: signal channel labels vector
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `yrev::Bool=false`: reverse Y axis
- `mono::Bool=false`: use color or grey palette

- `kwargs`: optional arguments for plot() function

# Returns

- `fig::GLMakie.Figure`
"""
function plot_erp_topo(locs::DataFrame, t::Vector{Float64}, signal::Array{Float64, 2}; channel=Union{Vector{Int64}, AbstractRange}, clabels::Vector{String}=[""], xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, yrev::Bool=false, kwargs...)

    size(signal, 2) == length(t) || throw(ArgumentError("Length of powers vector must equal length of frequencies vector."))
    length(channel) > nrow(locs) && throw(ArgumentError("Some channels do not have locations."))

    pal = mono == true ? :grays : :darktest
    
    # channel labels
    clabels == [""] && (clabels = repeat([""], size(s_pow, 1)))

    # get limits
    ylim = (floor(minimum(signal) * 1.1, digits=0), ceil(maximum(signal) * 1.1, digits=0))
    ylim = _tuple_max(ylim)

    # plot parameters
    plot_size = 1200
    marker_size = (150, 75)
    
    # get locations
    loc_x = zeros(size(locs, 1))
    loc_y = zeros(size(locs, 1))
    for idx in 1:size(locs, 1)
        loc_x[idx], loc_y[idx] = pol2cart(locs[!, :loc_radius][idx], locs[!, :loc_theta][idx])
    end
    # loc_x, loc_y = _locnorm(loc_x, loc_y)
    loc_x = loc_x[channel]
    loc_y = loc_y[channel]
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

    for idx in 1:size(signal, 1)
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
                        signal[idx, :],
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
        marker = load(marker_img)
        GLMakie.scatter!(fig_axis, (loc_x[idx], loc_y[idx]), marker=marker, markersize=marker_size)
        rm(marker_img)
    end

    return fig
end

"""
    plot_erp(obj; <keyword arguments>)

Plot ERP.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `channel::Union{Int64, Vector{Int64}, AbstractRange}`: channel(s) to plot
- `tm::Union{Int64, Vector{Int64}}=0`: time markers (in miliseconds) to plot as vertical lines, useful for adding topoplots at these time points 
- `xlabel::String="default"`: x-axis label, default is Time [ms]
- `ylabel::String="default"`: y-axis label, default is Amplitude [μV] 
- `title::String="default"`: plot title, default is ERP amplitude [channel: 1, epochs: 1:2, time window: -0.5 s:1.5 s]
- `cb::Bool=true`: plot color bar
- `cb_title::String="default"`: color bar title, default is Amplitude [μV] 
- `mono::Bool=false`: use color or grey palette
- `peaks::Bool=true`: draw peaks
- `channel_labels::Bool=true`: draw labels legend (using channel labels) for multi-channel `:butterfly` plot
- `type::Symbol=:normal`: plot type: `:normal`, mean ± 95%CI (`:mean`), butterfly plot (`:butterfly`), topographical plot of ERPs (`:topo`) or stacked epochs/channels (`:stack`)
- `yrev::Bool=false`: reverse Y axis
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_erp(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}, tm::Union{Int64, Vector{Int64}}=0, xlabel::String="default", ylabel::String="default", title::String="default", cb::Bool=true, cb_title::String="default", mono::Bool=false, peaks::Bool=true, channel_labels::Bool=true, type::Symbol=:normal, yrev::Bool=false, kwargs...)

    _check_var(type, [:normal, :butterfly, :mean, :topo, :stack], "type")

    type in [:normal, :mean] && length(channel) > 1 && throw(ArgumentError("For :normal and :mean plot types, only one channel must be specified."))

    # check channels
    _check_channels(obj, channel)

    # average all epochs
    epoch = 1:epoch_n(obj)

    signal = obj.data[channel, :, epoch]

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
        xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [ms]", "Amplitude [μV]", "ERP amplitude channel $(_channel2channel_name(channel))\n[averaged epochs: $epoch, time window: $t_s1:$t_s2]")
        signal = mean(signal, dims=2)[:]
        p = plot_erp(t,
                     signal,
                     xlabel=xlabel,
                     ylabel=ylabel,
                     title=title,
                     mono=mono,
                     yrev=yrev;
                     kwargs...)
    elseif type === :butterfly
        if length(channel) > 1
            xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [ms]", "Amplitude [μV]", "ERP amplitude channel$(_pl(length(channel))) $(_channel2channel_name(channel))\n[averaged epochs: $epoch, time window: $t_s1:$t_s2]")
            signal = mean(signal, dims=3)[:, :]
            if channel_labels == true
                clabels = labels(obj)[channel]
            else
                clabels = repeat([""], length(channel))
            end
        else
            xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [ms]", "Amplitude [μV]", "ERP amplitude channel$(_pl(length(channel))) $(_channel2channel_name(channel))\n[epochs: $epoch, time window: $t_s1:$t_s2]")
            signal = signal'
            clabels = [""]
        end
        p = plot_erp_butterfly(t,
                               signal,
                               xlabel=xlabel,
                               ylabel=ylabel,
                               title=title,
                               clabels=clabels,
                               mono=mono,
                               yrev=yrev;
                               kwargs...)
    elseif type === :mean
        xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [ms]", "Amplitude [μV]", "ERP amplitude [mean ± 95%CI] channel $(_channel2channel_name(channel))\n[averaged epoch$(_pl(length(epoch))): $epoch, time window: $t_s1:$t_s2]")
        p = plot_erp_avg(t,
                         signal,
                         xlabel=xlabel,
                         ylabel=ylabel,
                         title=title,
                         mono=mono,
                         yrev=yrev;
                         kwargs...)
    elseif type === :topo
        obj.header[:channel_locations] == false && throw(ArgumentError("Electrode locations not available."))
        xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "", "", "ERP amplitude channel$(_pl(length(channel))) $(_channel2channel_name(channel))\n[averaged epochs: $epoch, time window: $t_s1:$t_s2]")
        peaks = false
        signal = mean(signal, dims=3)[:, :]
        ndims(signal) == 1 && (signal = reshape(signal, 1, length(signal)))
        clabels = labels(obj)[channel]
        typeof(clabels) == String && (clabels = [clabels])
        p = plot_erp_topo(obj.locs,
                          t,
                          signal,
                          channel=channel,
                          clabels=clabels,
                          xlabel=xlabel,
                          ylabel=ylabel,
                          title=title,
                          mono=mono,
                          yrev=yrev;
                          kwargs...)
    elseif type === :stack
        peaks = false
        cb_title == "default" && (cb_title = "Amplitude [μV]")
        if length(channel) > 1
            xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [ms]", "", "ERP amplitude channel$(_pl(length(channel))) $(_channel2channel_name(channel))\n[averaged epochs: $epoch, time window: $t_s1:$t_s2]")
            signal = mean(signal, dims=3)[:, :]
            if channel_labels == true
                clabels = labels(obj)[channel]
            else
                clabels = repeat([""], length(channel))
            end
            ylabel = "Channel"
        else
            xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [ms]", "", "ERP amplitude channel$(_pl(length(channel))) $(_channel2channel_name(channel))\n[epochs: $epoch, time window: $t_s1:$t_s2]")
            signal = signal'
            clabels = [""]
            ylabel = "Epoch"
        end
        p = plot_erp_stack(t,
                           signal,
                           xlabel=xlabel,
                           ylabel=ylabel,
                           title=title,
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
    if peaks == true
        if length(channel) == 1
            erp = erp(obj).data
            pp = erp_peaks(obj)
            if mono == false
                Plots.scatter!((t[pp[channel, 1]], erp[channel, pp[channel, 1]]), marker=:xcross, markercolor=:red, markersize=3, label=false)
                Plots.scatter!((t[pp[channel, 2]], erp[channel, pp[channel, 2]]), marker=:xcross, markercolor=:blue, markersize=3, label=false)
            else
                Plots.scatter!((t[pp[channel, 1]], erp[channel, pp[channel, 1]]), marker=:xcross, markercolor=:black, markersize=3, label=false)
                Plots.scatter!((t[pp[channel, 2]], erp[channel, pp[channel, 2]]), marker=:xcross, markercolor=:black, markersize=3, label=false)
            end
            _info("Positive peak time: $(round(t[pp[channel, 1]] * 1000, digits=0)) ms")
            _info("Positive peak amplitude: $(round(erp[channel, pp[channel, 1]], digits=2)) μV")
            _info("Negative peak time: $(round(t[pp[channel, 2]] * 1000, digits=0)) ms")
            _info("Negative peak amplitude: $(round(erp[channel, pp[channel, 2]], digits=2)) μV")
        else
            erp = mean(erp(obj).data[channel, :], dims=1)[:]
            obj_tmp = keep_channel(obj, channel=1)
            obj_tmp.data = reshape(erp, 1, length(erp), 1)
            pp = erp_peaks(obj_tmp)
            if mono == false
                Plots.scatter!((t[pp[1, 1]], erp[pp[1, 1]]), marker=:xcross, markercolor=:red, markersize=3, label=false)
                Plots.scatter!((t[pp[1, 2]], erp[pp[1, 2]]), marker=:xcross, markercolor=:blue, markersize=3, label=false)
            else
                Plots.scatter!((t[pp[1, 1]], erp[pp[1, 1]]), marker=:xcross, markercolor=:black, markersize=3, label=false)
                Plots.scatter!((t[pp[1, 2]], erp[pp[1, 2]]), marker=:xcross, markercolor=:black, markersize=3, label=false)
            end
            _info("Positive peak time: $(round(t[pp[1, 1]] * 1000, digits=0)) ms")
            _info("Positive peak amplitude: $(round(erp[pp[1, 1]], digits=2)) μV")
            _info("Negative peak time: $(round(t[pp[1, 2]] * 1000, digits=0)) ms")
            _info("Negative peak amplitude: $(round(erp[pp[1, 2]], digits=2)) μV")
        end
    end

    if type !== :topo
        Plots.plot(p)
    else
        GLMakie.show(p)
    end

    return p
end

"""
    plot_erp_stack(signal; <keyword arguments>)

Plot EPRs stacked by channels or by epochs.

# Arguments

- `t::AbstractVector`: x-axis values
- `signal::AbstractArray`
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
function plot_erp_stack(t::AbstractVector, signal::AbstractArray; clabels::Vector{String}=[""], xlabel::String="", ylabel::String="", title::String="", cb::Bool=true, cb_title::String="", mono::Bool=false, kwargs...)

    ndims(signal) == 2 || throw(ArgumentError("signal must have 2 dimensions."))
    length(t) == size(signal, 2) || throw(ArgumentError("Number of signal columns ($(size(signal, 2))) must be equal to length of x-axis values ($(length(t)))."))

    pal = mono == true ? :grays : :darktest

    if clabels == [""]
        yticks = round.(Int64, range(1, size(signal, 1), length=10))
    else
        yticks = (1:size(signal, 1), clabels)
    end
    p = Plots.heatmap(t,
                      1:size(signal, 1),
                      signal,
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
    plot_erp(obj, c; <keyword arguments>)

Plot ERP.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `c::Union{Symbol, AbstractArray}`: component to plot
- `c_idx::Union{Int64, Vector{Int64}, AbstractRange}=0`: component channel to display, default is all component channels
- `tm::Union{Int64, Vector{Int64}}=0`: time markers (in miliseconds) to plot as vertical lines, useful for adding topoplots at these time points 
- `xlabel::String="default"`: x-axis label, default is Time [ms]
- `ylabel::String="default"`: y-axis label, default is Amplitude [μV] 
- `title::String="default"`: plot title, default is ERP amplitude [component: 1, epochs: 1:2, time window: -0.5 s:1.5 s]
- `cb::Bool=true`: plot color bar
- `cb_title::String="default"`: color bar title
- `mono::Bool=false`: use color or grey palette
- `peaks::Bool=true`: draw peaks
- `channel_labels::Bool=true`: draw labels legend (using component labels) for multi-channel `:butterfly` plot
- `type::Symbol=:normal`: plot type: `:normal`, mean ± 95%CI (`:mean`), butterfly plot (`:butterfly`), topographical plot of ERPs (`:topo`) or stacked epochs/channels (`:stack`)
- `yrev::Bool=false`: reverse Y axis
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_erp(obj::NeuroAnalyzer.NEURO, c::Union{Symbol, AbstractArray}; c_idx::Union{Int64, Vector{Int64}, AbstractRange}=0, tm::Union{Int64, Vector{Int64}}=0, xlabel::String="default", ylabel::String="default", title::String="default", cb::Bool=true, cb_title::String="default", mono::Bool=false, peaks::Bool=true, channel_labels::Bool=true, type::Symbol=:normal, yrev::Bool=false, kwargs...)

    _check_var(type, [:normal, :butterfly, :mean, :topo, :stack], "type")

    # select component channels, default is all channels
    typeof(c) == Symbol && (c = _get_component(obj, c).c)
    c_idx == 0 && (c_idx = _select_cidx(c, c_idx))
    _check_cidx(c, c_idx)
    if clabels == true
        clabels = _gen_clabels(c)[c_idx]
    else
        clabels = repeat([""], length(c_idx))
    end
    length(c_idx) == 1 && (clabels = [clabels])

    type in [:normal, :mean] && length(c_idx) > 1 && throw(ArgumentError("For :normal and :mean plot types, only one component channel must be specified."))

    # average all epochs
    epoch = 1:epoch_n(obj)

    signal = c[c_idx, :, epoch]

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
        xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [ms]", "Amplitude [μV]", "ERP amplitude component $(_channel2channel_name(c_idx))\n[averaged epochs: $epoch, time window: $t_s1:$t_s2]")
        signal = mean(signal, dims=2)[:]
        p = plot_erp(t,
                     signal,
                     xlabel=xlabel,
                     ylabel=ylabel,
                     title=title,
                     mono=mono,
                     yrev=yrev;
                     kwargs...)
    elseif type === :butterfly
        if length(c_idx) > 1
            xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [ms]", "Amplitude [μV]", "ERP amplitude component$(_pl(length(c_idx))) $(_channel2channel_name(c_idx))\n[averaged epochs: $epoch, time window: $t_s1:$t_s2]")
            signal = mean(signal, dims=3)[:, :]
            if channel_labels == true
                clabels = _gen_clabels(c)[c_idx]
            else
                clabels = repeat([""], length(c_idx))
            end
        else
            xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [ms]", "Amplitude [μV]", "ERP amplitude component$(_pl(length(c_idx))) $(_channel2channel_name(c_idx))\n[epochs: $epoch, time window: $t_s1:$t_s2]")
            signal = signal'
            clabels = [""]
        end
        p = plot_erp_butterfly(t,
                               signal,
                               xlabel=xlabel,
                               ylabel=ylabel,
                               title=title,
                               clabels=clabels,
                               mono=mono,
                               yrev=yrev;
                               kwargs...)
    elseif type === :mean
        xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [ms]", "Amplitude [μV]", "ERP amplitude [mean ± 95%CI] component $(_channel2channel_name(c_idx))\n[averaged epoch$(_pl(length(epoch))): $epoch, time window: $t_s1:$t_s2]")
        p = plot_erp_avg(t,
                         signal,
                         xlabel=xlabel,
                         ylabel=ylabel,
                         title=title,
                         mono=mono,
                         yrev=yrev;
                         kwargs...)
    elseif type === :topo
        obj.header[:channel_locations] == false && throw(ArgumentError("Electrode locations not available."))
        xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "", "", "ERP amplitude component$(_pl(length(c_idx))) $(_channel2channel_name(c_idx))\n[averaged epochs: $epoch, time window: $t_s1:$t_s2]")
        peaks = false
        signal = mean(signal, dims=3)[:, :]
        ndims(signal) == 1 && (signal = reshape(signal, 1, length(signal)))
        typeof(clabels) == String && (clabels = [clabels])
        p = plot_erp_topo(obj.locs,
                          t,
                          signal,
                          channel=c_idx,
                          clabels=clabels,
                          xlabel=xlabel,
                          ylabel=ylabel,
                          title=title,
                          mono=mono;
                          kwargs...)
    elseif type === :stack
        cb_title == "default" && (cb_title = "Amplitude [μV]")
        peaks = false
        if length(c_idx) > 1
            xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [ms]", "", "ERP amplitude component$(_pl(length(c_idx))) $(_channel2channel_name(c_idx))\n[averaged epochs: $epoch, time window: $t_s1:$t_s2]")
            signal = mean(signal, dims=3)[:, :]
            ylabel = "Component channel"
        else
            xlabel, ylabel, title = _set_defaults(xlabel, ylabel, title, "Time [ms]", "", "ERP amplitude component$(_pl(length(c_idx))) $(_channel2channel_name(c_idx))\n[epochs: $epoch, time window: $t_s1:$t_s2]")
            signal = signal'
            clabels = [""]
            ylabel = "Epoch"
        end
        p = plot_erp_stack(t,
                           signal,
                           xlabel=xlabel,
                           ylabel=ylabel,
                           title=title,
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
    if peaks == true
        signal = c[c_idx, :, epoch]
        if length(c_idx) == 1
            erp = mean(signal, dims=2)[:]
            obj_tmp = keep_channel(obj, channel=1)
            obj_tmp.data = reshape(erp, 1, length(erp), 1)
            pp = erp_peaks(obj_tmp)
            if mono == false
                Plots.scatter!((t[pp[1, 1]], erp[pp[1, 1]]), marker=:xcross, markercolor=:red, markersize=3, label=false)
                Plots.scatter!((t[pp[1, 2]], erp[pp[1, 2]]), marker=:xcross, markercolor=:blue, markersize=3, label=false)
            else
                Plots.scatter!((t[pp[1, 1]], erp[pp[1, 1]]), marker=:xcross, markercolor=:black, markersize=3, label=false)
                Plots.scatter!((t[pp[1, 2]], erp[pp[1, 2]]), marker=:xcross, markercolor=:black, markersize=3, label=false)
            end
            _info("Positive peak time: $(round(t[pp[1, 1]] * 1000, digits=0)) ms")
            _info("Positive peak amplitude: $(round(erp[pp[1, 1]], digits=2)) μV")
            _info("Negative peak time: $(round(t[pp[1, 2]] * 1000, digits=0)) ms")
            _info("Negative peak amplitude: $(round(erp[pp[1, 2]], digits=2)) μV")
        else         
            erp = mean(mean(signal, dims=3), dims=1)[:]
            obj_tmp = keep_channel(obj, channel=1)
            obj_tmp.data = reshape(erp, 1, length(erp), 1)
            pp = erp_peaks(obj_tmp)
            if mono == false
                Plots.scatter!((t[pp[1, 1]], erp[pp[1, 1]]), marker=:xcross, markercolor=:red, markersize=3, label=false)
                Plots.scatter!((t[pp[1, 2]], erp[pp[1, 2]]), marker=:xcross, markercolor=:blue, markersize=3, label=false)
            else
                Plots.scatter!((t[pp[1, 1]], erp[pp[1, 1]]), marker=:xcross, markercolor=:black, markersize=3, label=false)
                Plots.scatter!((t[pp[1, 2]], erp[pp[1, 2]]), marker=:xcross, markercolor=:black, markersize=3, label=false)
            end
            _info("Positive peak time: $(round(t[pp[1, 1]] * 1000, digits=0)) ms")
            _info("Positive peak amplitude: $(round(erp[pp[1, 1]], digits=2)) μV")
            _info("Negative peak time: $(round(t[pp[1, 2]] * 1000, digits=0)) ms")
            _info("Negative peak amplitude: $(round(erp[pp[1, 2]], digits=2)) μV")
        end
    end

    if type !== :topo
        Plots.plot(p)
    else
        GLMakie.show(p)
    end

    return p
end
