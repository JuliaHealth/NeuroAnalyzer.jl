export plot_mep
export plot_mep_butterfly
export plot_mep_avg
export plot_mep_stack

"""
    plot_mep(t, s, bad; <keyword arguments>)

Plot MEP.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
- `s::AbstractVector`: data to plot
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `yrev::Bool=false`: reverse Y axis
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_mep(t::Union{AbstractVector, AbstractRange}, s::AbstractVector; xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, yrev::Bool=false, kwargs...)

    pal = mono ? :grays : :darktest

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
    yrev && yflip!(true)

    # plot 0 h-line
    p = Plots.hline!([0],
                     color=:grey,
                     lw=0.5,
                     labels="")

    # plot MEP
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
    plot_mep_butterfly(t, s; <keyword arguments>)

Butterfly plot of MEP.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
- `s::AbstractArray`: data to plot
- `clabels::Vector{String}=[""]`: signal channel labels vector
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `avg::Bool=false`: plot average MEP
- `yrev::Bool=false`: reverse Y axis
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_mep_butterfly(t::Union{AbstractVector, AbstractRange}, s::AbstractArray; clabels::Vector{String}=[""], xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, avg::Bool=true, yrev::Bool=false, kwargs...)

    pal = mono ? :grays : :darktest

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
    yrev && yflip!(true)

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

    # plot averaged MEP
    if avg
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
    plot_mep_avg(t, s; <keyword arguments>)

Plot MEP amplitude mean and ±95% CI.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
- `s::AbstractArray`: data to plot
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `yrev::Bool=false`: reverse Y axis
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_mep_avg(t::Union{AbstractVector, AbstractRange}, s::AbstractArray; xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, yrev::Bool=false, kwargs...)

    pal = mono ? :grays : :darktest

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
    yrev && yflip!(true)

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
    plot_mep_stack(s; <keyword arguments>)

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
- `mono::Bool=false`: use color or gray palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_mep_stack(t::AbstractVector, s::AbstractArray; clabels::Vector{String}=[""], xlabel::String="", ylabel::String="", title::String="", cb::Bool=true, cb_title::String="", mono::Bool=false, kwargs...)

    @assert ndims(s) == 2 "signal must have 2 dimensions."
    @assert length(t) == size(s, 2) "Number of signal columns ($(size(s, 2))) must be equal to length of x-axis values ($(length(t)))."

    pal = mono ? :grays : :darktest

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
    plot_mep(obj; <keyword arguments>)

Plot MEP.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{String, Vector{String}}`: list of channels
- `tm::Union{Int64, Vector{Int64}}=0`: time markers (in miliseconds) to plot as vertical lines, useful for adding topoplots at these time points
- `xlabel::String="default"`: x-axis label, default is Time [ms]
- `ylabel::String="default"`: y-axis label, default is Amplitude [units]
- `title::String="default"`: plot title, default is MEP amplitude [channel: 1, epochs: 1:2, time window: -0.5 s:1.5 s]
- `cb::Bool=true`: plot color bar
- `cb_title::String="default"`: color bar title, default is Amplitude [units]
- `mono::Bool=false`: use color or gray palette
- `peaks::Bool=true`: draw peaks
- `peaks_detect::Bool=true`: if true, detect MEP peaks, otherwise use embedded
- `channel_labels::Bool=true`: draw labels legend (using channel labels) for multi-channel `:butterfly` plot
- `type::Symbol=:normal`: plot type: `:normal`, butterfly plot (`:butterfly`), topographical plot of ERPs (`:topo`) or stacked epochs/channels (`:stack`)
- `yrev::Bool=false`: reverse Y axis
- `avg::Bool=false`: plot average MEP for `:butterfly` plot
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_mep(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}}, tm::Union{Int64, Vector{Int64}}=0, xlabel::String="default", ylabel::String="default", title::String="default", cb::Bool=true, cb_title::String="default", mono::Bool=false, peaks::Bool=true, peaks_detect::Bool=true, channel_labels::Bool=true, type::Symbol=:normal, yrev::Bool=false, avg::Bool=true, kwargs...)

    _check_datatype(obj, "mep")

    # check channels
    ch = get_channel(obj, ch=ch)

    # set units
    units = _ch_units(obj, ch[1])

    _check_var(type, [:normal, :butterfly, :mean, :stack], "type")
    @assert !(length(ch) > 1 && length(unique(obj.header.recording[:channel_type][ch])) > 1) "All channels must be of the same type."

    # get data
    if type === :normal
        s = obj.data[ch, :, 1]
    else
        if ch isa Int64
            s = obj.data[ch, :, 1]'
            channel_labels = false
            # s = reshape(obj.data[ch, :, 1], 1, :, ep_n)
        else
            s = obj.data[ch, :, 1]
        end
    end

    # get time vector
    t = obj.time_pts
    _, t_s1, _, t_s2 = _convert_t(t[1], t[end])

    if tm != 0
        for tm_idx in eachindex(tm)
            @assert tm[tm_idx] / 1000 >= t[1] "tm value ($(tm[tm_idx])) is out of time segment ($(t[1]):$(t[end]))."
            @assert tm[tm_idx] / 1000 <= t[end] "tm value ($(tm[tm_idx])) is out of time segment ($(t[1]):$(t[end]))."
            tm[tm_idx] = vsearch(tm[tm_idx] / 1000, t)
        end
    end

    if type === :normal
        @assert ch isa Int64 "For :normal plot type, only one channel must be specified."
        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [ms]", "Amplitude [$units]", "MEP amplitude $(_channel2channel_name(ch))\n[time window: $t_s1:$t_s2]")
        p = plot_mep(t,
                     s,
                     xlabel=xl,
                     ylabel=yl,
                     title=tt,
                     mono=mono,
                     yrev=yrev;
                     kwargs...)
    elseif type === :butterfly
        @assert !(ch isa Int64) "For :butterfly plot type, more than one channel must be specified."
        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [ms]", "Amplitude [$units]", "MEP amplitude$(_pl(length(ch))) $(_channel2channel_name(ch))\n[time window: $t_s1:$t_s2]")
        if channel_labels
            clabels = labels(obj)[ch]
        else
            clabels = repeat([""], length(ch))
        end
        p = plot_mep_butterfly(t,
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
        @assert !(ch isa Int64) "For :mean plot type, more than one channel must be specified."
        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [ms]", "Amplitude [$units]", "MEP amplitude [mean ± 95%CI] signal$(_pl(length(ch))) $(_channel2channel_name(ch))\n[time window: $t_s1:$t_s2]")
        p = plot_mep_avg(t,
                         s,
                         xlabel=xl,
                         ylabel=yl,
                         title=tt,
                         mono=mono,
                         yrev=yrev;
                         kwargs...)
    elseif type === :stack
        @assert !(ch isa Int64) "For :stack plot type, more than one channel must be specified."
        peaks = false
        cb_title == "default" && (cb_title = "Amplitude [$units]")
        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [ms]", "", "MEP amplitude$(_pl(length(ch))) $(_channel2channel_name(ch))\n[time window: $t_s1:$t_s2]")
        if channel_labels
            clabels = labels(obj)[ch]
        else
            clabels = repeat([""], length(ch))
        end
        p = plot_mep_stack(t,
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
    if peaks
        if peaks_detect
            if ch isa Int64
                pp = mep_peaks(obj)
                if !mono
                    Plots.scatter!((t[pp[ch, 1]], obj.data[ch, pp[ch, 1], 1]), marker=:xcross, markercolor=:red, markersize=3, label=false)
                    Plots.scatter!((t[pp[ch, 2]], obj.data[ch, pp[ch, 2], 1]), marker=:xcross, markercolor=:blue, markersize=3, label=false)
                else
                    Plots.scatter!((t[pp[ch, 1]], obj.data[ch, pp[ch, 1], 1]), marker=:xcross, markercolor=:black, markersize=3, label=false)
                    Plots.scatter!((t[pp[ch, 2]], obj.data[ch, pp[ch, 2], 1]), marker=:xcross, markercolor=:black, markersize=3, label=false)
                end
                _info("Positive peak time: $(round(t[pp[ch, 1]] * 1000, digits=0)) ms")
                _info("Positive peak amplitude: $(round(obj.data[ch, pp[ch, 1]], digits=2)) $units")
                _info("Negative peak time: $(round(t[pp[ch, 2]] * 1000, digits=0)) ms")
                _info("Negative peak amplitude: $(round(obj.data[ch, pp[ch, 2]], digits=2)) $units")
            elseif (type === :butterfly && avg) || type === :mean
                mep_tmp = mean(mean(obj.data[ch, :, :], dims=1), dims=3)
                obj_tmp = keep_channel(obj, ch=1)
                obj_tmp.data = mep_tmp
                pp = mep_peaks(obj_tmp)
                if !mono
                    Plots.scatter!((t[pp[1, 1]], mep_tmp[pp[1, 1]]), marker=:xcross, markercolor=:red, markersize=3, label=false)
                    Plots.scatter!((t[pp[1, 2]], mep_tmp[pp[1, 2]]), marker=:xcross, markercolor=:blue, markersize=3, label=false)
                else
                    Plots.scatter!((t[pp[1, 1]], mep_tmp[pp[1, 1]]), marker=:xcross, markercolor=:black, markersize=3, label=false)
                    Plots.scatter!((t[pp[1, 2]], mep_tmp[pp[1, 2]]), marker=:xcross, markercolor=:black, markersize=3, label=false)
                end
                _info("Positive peak time: $(round(t[pp[1, 1]] * 1000, digits=0)) ms")
                _info("Positive peak amplitude: $(round(mep_tmp[pp[1, 1]], digits=2)) $units")
                _info("Negative peak time: $(round(t[pp[1, 2]] * 1000, digits=0)) ms")
                _info("Negative peak amplitude: $(round(mep_tmp[pp[1, 2]], digits=2)) $units")
            end
        else
            if ch isa Int64
                pp = hcat(obj.header.recording[:markers_pos], obj.header.recording[:markers_neg])
                if pp[ch, 1] != 0 && pp[ch, 2] != 0
                    if !mono
                        Plots.scatter!((t[pp[ch, 1]], obj.data[ch, pp[ch, 1], 1]), marker=:xcross, markercolor=:red, markersize=3, label=false)
                        Plots.scatter!((t[pp[ch, 2]], obj.data[ch, pp[ch, 2], 1]), marker=:xcross, markercolor=:blue, markersize=3, label=false)
                    else
                        Plots.scatter!((t[pp[ch, 1]], obj.data[ch, pp[ch, 1], 1]), marker=:xcross, markercolor=:black, markersize=3, label=false)
                        Plots.scatter!((t[pp[ch, 2]], obj.data[ch, pp[ch, 2], 1]), marker=:xcross, markercolor=:black, markersize=3, label=false)
                    end
                    _info("Positive peak time: $(round(t[pp[ch, 1]] * 1000, digits=0)) ms")
                    _info("Positive peak amplitude: $(round(obj.data[ch, pp[ch, 1]], digits=2)) $units")
                    _info("Negative peak time: $(round(t[pp[ch, 2]] * 1000, digits=0)) ms")
                    _info("Negative peak amplitude: $(round(obj.data[ch, pp[ch, 2]], digits=2)) $units")
                else
                    _warn("Peaks information not available, use `peaks_detect=true` parameter.")
                end
            else
                _warn("Cannot use embedded peaks for multichannel plot when `avg=false`.")
            end
        end
    end

    Plots.plot(p)

    return p

end
