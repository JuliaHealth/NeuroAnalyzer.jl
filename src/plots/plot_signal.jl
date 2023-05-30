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
- `units::String=""`: units of the scale
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_signal(t::Union{AbstractVector, AbstractRange}, s::Union{AbstractVector, AbstractArray}; clabels::Vector{String}=[""], xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, scale::Bool=true, units::String="", kwargs...)

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
    ch_n <= 2 && (plot_size = (1200, 400))
    ch_n > 2 && (plot_size = (1200, 800))
    p = Plots.plot(xlabel=xlabel,
                   ylabel=ylabel,
                   xlims=_xlims(t),
                   xticks=_ticks(t),
                   ylims=(-1, ch_n),
                   title=title,
                   palette=pal,
                   size=plot_size,
                   top_margin=0Plots.px,
                   bottom_margin=15Plots.px,
                   right_margin=10Plots.px,
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
        p = Plots.plot!([t[1], t[1]], [(ch_n - 1.5), (ch_n - 0.5)], color=:red, linewidth=2, label="")
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
- `units::String=""`: units of the scale
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_signal(t::Union{AbstractVector, AbstractRange}, s::Union{AbstractVector, AbstractArray}, bad::Vector{Bool}; clabels::Vector{String}=[""], xlabel::String="", ylabel::String="", title::String="", scale::Bool=true, units::String="", kwargs...)

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
    ch_n <= 2 && (plot_size = (1200, 400))
    ch_n > 2 && (plot_size = (1200, 800))
    p = Plots.plot(xlabel=xlabel,
                   ylabel=ylabel,
                   xlims=_xlims(t),
                   xticks=_ticks(t),
                   ylims=(-1, ch_n),
                   title=title,
                   palette=pal,
                   size=plot_size,
                   margins=10Plots.px,
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
- `units::String=""`: units of the scale
- `norm::Bool=false`: normalize to -1 .. +1
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_signal_avg(t::Union{AbstractVector, AbstractRange}, s::AbstractArray; xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, scale::Bool=true, units::String="", norm::Bool=false, kwargs...)

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
- `units::String=""`: units of the scale
- `mono::Bool=false`: use color or grey palette
- `norm::Bool=false`: normalize to -1 .. +1
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_signal_butterfly(t::Union{AbstractVector, AbstractRange}, s::AbstractArray; clabels::Vector{String}=[""], xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, scale::Bool=true, units::String="", norm::Bool=false, kwargs...)

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
- `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
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
- `bad::Union{Bool, Matrix{Bool}}=false`: list of bad channels; if not empty - plot bad channels using this list
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot(obj::NeuroAnalyzer.NEURO; ep::Union{Int64, AbstractRange}=0, ch::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(channel_n(obj)), seg::Tuple{Real, Real}=(0, 10), xlabel::String="default", ylabel::String="default", title::String="default", mono::Bool=false, emarkers::Bool=true, markers::Bool=true, scale::Bool=true, units::String="", type::Symbol=:normal, norm::Bool=false, bad::Union{Bool, Matrix{Bool}}=false, kwargs...)

    if signal_len(obj) < 10 * sr(obj) && seg == (0, 10)
        seg = (0, obj.time_pts[end])
    else
        _check_segment(obj, seg)
    end
    seg = (vsearch(seg[1], obj.time_pts), vsearch(seg[2], obj.time_pts))

    _check_var(type, [:normal, :butterfly, :mean, :stack], "type")

    if ep != 0
        _check_epochs(obj, ep)
        if epoch_n(obj) == 1
            ep = 0
        else
            seg = (((ep[1] - 1) * epoch_len(obj) + 1), seg[2])
            if ep isa Int64
                seg = (seg[1], (seg[1] + epoch_len(obj) - 1))
            else
                seg = (seg[1], (ep[end] * epoch_len(obj)))
            end
            ep = 0
        end
    end

    # do not show epoch markers if there are no epochs
    epoch_n(obj) == 1 && (emarkers = false)
    if emarkers == true
        epoch_markers = _get_epoch_markers(obj)
    end

    # check channels
    _check_channels(obj, ch)
    clabels = labels(obj)

    # set units
    units = _set_units(obj, ch[1])

    # get time vector
    if seg[2] <= epoch_len(obj)
        s = obj.data[:, seg[1]:seg[2], 1]
    else
        s = epoch(obj, ep_n=1).data[:, seg[1]:seg[2], 1]
    end
    #t = _get_t(seg[1], seg[2], sr(obj))
    t = obj.time_pts[seg[1]:seg[2]]

    _, t_s1, _, t_s2 = _convert_t(t[1], t[end])
    ep = _s2epoch(obj, seg[1], seg[2])

    ch_t = obj.header.recording[:channel_type]
    ch_tmp = Vector{Vector{Int64}}()
    ch_t_uni = nothing
    if length(ch) > 1
        ch_t_uni = unique(ch_t[ch])
        for cht_idx in 1:length(ch_t_uni)
            ch_tmp2 = Vector{Int64}()
            for ch_idx in 1:length(ch)
                ch_t[ch[ch_idx]] == ch_t_uni[cht_idx] && push!(ch_tmp2, ch[ch_idx])
            end
            push!(ch_tmp, ch_tmp2)
        end
    else
        ch_t_uni = ch_t[ch]
        ch_tmp = [[ch]]
    end

    xl, yl, tt = "", "", ""

    p = Plots.Plot[]

    if type === :normal
        if bad == false
            if length(ch_tmp) > 1
                for cht_idx in 1:length(ch_t_uni)
                    units = _set_units(obj, ch_tmp[cht_idx][1])
                    if ch_t[ch_tmp[cht_idx][1]] == "eeg"
                        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "EEG channel$(_pl(length(ch_tmp[cht_idx]))) ($(_channel2channel_name(ch_tmp[cht_idx])))")
                    end
                    if ch_t[ch_tmp[cht_idx][1]] == "grad"
                        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "Gradiometer$(_pl(length(ch_tmp[cht_idx]))) ($(_channel2channel_name(ch_tmp[cht_idx])))")
                    end
                    if ch_t[ch_tmp[cht_idx][1]] == "mag"
                        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "Magnetometer$(_pl(length(ch_tmp[cht_idx]))) ($(_channel2channel_name(ch_tmp[cht_idx])))")
                    end
                    if ch_t[ch_tmp[cht_idx][1]] == "nirs_int"
                        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "NIRS intensity channel$(_pl(length(ch_tmp[cht_idx]))) ($(_channel2channel_name(ch_tmp[cht_idx])))")
                    end
                    if ch_t[ch_tmp[cht_idx][1]] == "nirs_od"
                        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "NIRS optical density channel$(_pl(length(ch_tmp[cht_idx]))) ($(_channel2channel_name(ch_tmp[cht_idx])))")
                    end
                    if ch_t[ch_tmp[cht_idx][1]] == "nirs_hbo"
                        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "NIRS HbO concentration channel$(_pl(length(ch_tmp[cht_idx]))) ($(_channel2channel_name(ch_tmp[cht_idx])))")
                    end
                    if ch_t[ch_tmp[cht_idx][1]] == "nirs_hbr"
                        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "NIRS HbR concentration channel$(_pl(length(ch_tmp[cht_idx]))) ($(_channel2channel_name(ch_tmp[cht_idx])))")
                    end
                    if ch_t[ch_tmp[cht_idx][1]] == "nirs_hbt"
                        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "NIRS HbT concentration channel$(_pl(length(ch_tmp[cht_idx]))) ($(_channel2channel_name(ch_tmp[cht_idx])))")
                    end
                    if ch_t[ch_tmp[cht_idx][1]] == "nirs_aux"
                        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "NIRS AUX channel$(_pl(length(ch_tmp[cht_idx]))) ($(_channel2channel_name(ch_tmp[cht_idx])))")
                    end
                    if ch_t[ch_tmp[cht_idx][1]] == "ref"
                        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "Reference channel$(_pl(length(ch_tmp[cht_idx]))) ($(_channel2channel_name(ch_tmp[cht_idx])))")
                    end
                    if ch_t[ch_tmp[cht_idx][1]] == "eog"
                        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "EOG channel$(_pl(length(ch_tmp[cht_idx]))) ($(_channel2channel_name(ch_tmp[cht_idx])))")
                    end
                    if ch_t[ch_tmp[cht_idx][1]] == "emg"
                        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "EMG channel$(_pl(length(ch_tmp[cht_idx]))) ($(_channel2channel_name(ch_tmp[cht_idx])))")
                    end
                    if ch_t[ch_tmp[cht_idx][1]] == "ecg"
                        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "ECG channel ($(_channel2channel_name(ch_tmp[cht_idx])))")
                    end
                    if ch_t[ch_tmp[cht_idx][1]] == "other"
                        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "Other channel$(_pl(length(ch_tmp[cht_idx]))) ($(_channel2channel_name(ch_tmp[cht_idx])))")
                    end
                    if ch_t[ch_tmp[cht_idx][1]] == "mrk"
                        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "Marker$(_pl(length(ch_tmp[cht_idx]))) ($(_channel2channel_name(ch_tmp[cht_idx])))")
                    end

                    cht_idx < length(ch_t_uni) && (xl = "")

                    p_tmp = plot_signal(t,
                                        s[ch_tmp[cht_idx], :],
                                        clabels=clabels[ch_tmp[cht_idx]],
                                        xlabel=xl,
                                        ylabel=yl,
                                        title=tt,
                                        scale=scale,
                                        units=units,
                                        mono=mono;
                                        kwargs...)
                    push!(p, p_tmp)

                end
            else

                if ch_t[ch_tmp[1][1]] == "eeg"
                    xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "EEG channel$(_pl(length(ch_tmp[1]))) ($(_channel2channel_name(ch_tmp[1])))\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
                end
                if ch_t[ch_tmp[1][1]] == "grad"
                    xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "Gradiometer$(_pl(length(ch_tmp[1]))) ($(_channel2channel_name(ch_tmp[1])))\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
                end
                if ch_t[ch_tmp[1][1]] == "mag"
                    xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "Magnetometer$(_pl(length(ch_tmp[1]))) ($(_channel2channel_name(ch_tmp[1])))\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
                end
                if ch_t[ch_tmp[1][1]] == "nirs_int"
                    xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "NIRS intensity channel$(_pl(length(ch_tmp[1]))) ($(_channel2channel_name(ch_tmp[1])))\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
                end
                if ch_t[ch_tmp[1][1]] == "nirs_od"
                    xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "NIRS optical density channel$(_pl(length(ch_tmp[1]))) ($(_channel2channel_name(ch_tmp[1])))\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
                end
                if ch_t[ch_tmp[1][1]] == "nirs_hbo"
                    xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "NIRS HbO concentration channel$(_pl(length(ch_tmp[1]))) ($(_channel2channel_name(ch_tmp[1])))\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
                end
                if ch_t[ch_tmp[1][1]] == "nirs_hbr"
                    xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "NIRS HbR concentration channel$(_pl(length(ch_tmp[1]))) ($(_channel2channel_name(ch_tmp[1])))\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
                end
                if ch_t[ch_tmp[1][1]] == "nirs_hbt"
                    xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "NIRS HbT concentration channel$(_pl(length(ch_tmp[1]))) ($(_channel2channel_name(ch_tmp[1])))\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
                end
                if ch_t[ch_tmp[1][1]] == "nirs_aux"
                    xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "NIRS AUX channel$(_pl(length(ch_tmp[1]))) ($(_channel2channel_name(ch_tmp[1])))\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
                end
                if ch_t[ch_tmp[1][1]] == "ref"
                    xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "Reference channel$(_pl(length(ch_tmp[1]))) ($(_channel2channel_name(ch_tmp[1])))\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
                end
                if ch_t[ch_tmp[1][1]] == "eog"
                    xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "EOG channel$(_pl(length(ch_tmp[1]))) ($(_channel2channel_name(ch_tmp[1])))\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
                end
                if ch_t[ch_tmp[1][1]] == "emg"
                    xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "EMG channel$(_pl(length(ch_tmp[1]))) ($(_channel2channel_name(ch_tmp[1])))\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
                end
                if ch_t[ch_tmp[1][1]] == "ecg"
                    xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "ECG channel ($(_channel2channel_name(ch_tmp[1])))")
                end
                if ch_t[ch_tmp[1][1]] == "other"
                    xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "Other channel$(_pl(length(ch_tmp[1]))) ($(_channel2channel_name(ch_tmp[1])))\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
                end
                if ch_t[ch_tmp[1][1]] == "mrk"
                    xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "Marker$(_pl(length(ch_tmp[1]))) ($(_channel2channel_name(ch_tmp[1])))\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
                end
                if length(ch) == 1
                    p = plot_signal(t,
                                    s[ch, :],
                                    clabels=[clabels[ch_tmp[1][1]]],
                                    xlabel=xl,
                                    ylabel=yl,
                                    title=tt,
                                    scale=scale,
                                    units=units,
                                    mono=mono;
                                    kwargs...)
                else
                    p = plot_signal(t,
                                    s[ch, :],
                                    clabels=clabels[ch_tmp[1][:]],
                                    xlabel=xl,
                                    ylabel=yl,
                                    title=tt,
                                    scale=scale,
                                    units=units,
                                    mono=mono;
                                    kwargs...)
                end
            end
        else
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "Bad channel$(_pl(length(ch))) $(_channel2channel_name(ch))\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
            length(ch) > size(bad, 1) && throw(ArgumentError("Number of channels cannot be larger than number of bad channels rows."))
            ep > size(bad, 2) && throw(ArgumentError("Epoch number cannot be larger than number of bad channels columns."))
            p = plot_signal(t,
                            s[ch, :],
                            bad[ch, ep],
                            clabels=clabels[ch],
                            xlabel=xl,
                            ylabel=yl,
                            title=tt,
                            scale=scale,
                            units=units;
                            kwargs...)
        end
    end
    
    if type === :butterfly
        length(ch_t_uni) > 1 && throw(ArgumentError("For type=:butterfly plot all channels should be of the same type."))
        size(s, 1) == 1 && throw(ArgumentError("For type=:butterfly plot the signal must contain ≥ 2 channels."))
        units = _set_units(obj, ch[1])
        if ch_t[ch[1]] == "eeg"
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "Amplitude [$units]", "EEG channel$(_pl(length(ch))) $(_channel2channel_name(ch))\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        end
        if ch_t[ch[1]] == "grad"
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "Amplitude [$units]", "Gradiometer$(_pl(length(ch))) $(_channel2channel_name(ch))\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        end
        if ch_t[ch[1]] == "mag"
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "Amplitude [$units]", "Magnetometer$(_pl(length(ch))) $(_channel2channel_name(ch))\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        end
        if ch_t[ch[1]] == "nirs_int"
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "Intensity [$units]", "NIRS intensity channel$(_pl(length(ch))) $(_channel2channel_name(ch))\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        end
        if ch_t[ch[1]] == "nirs_od"
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "OD [$units]", "NIRS optical density channel$(_pl(length(ch))) $(_channel2channel_name(ch))\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        end
        if ch_t[ch[1]] == "nirs_hbo"
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "HbO concentration [$units]", "NIRS HbO concentration channel$(_pl(length(ch))) $(_channel2channel_name(ch))\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        end
        if ch_t[ch[1]] == "nirs_hbr"
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "HbR concentration [$units]", "NIRS HbR concentration channel$(_pl(length(ch))) $(_channel2channel_name(ch))\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        end
        if ch_t[ch[1]] == "nirs_hbt"
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "HbT concentration [$units]", "NIRS HbT concentration channel$(_pl(length(ch))) $(_channel2channel_name(ch))\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        end
        if ch_t[ch[1]] == "nirs_aux"
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "Amplitude [$units]", "NIRS AUX channel$(_pl(length(ch))) $(_channel2channel_name(ch))\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        end
        if ch_t[ch[1]] == "ref"
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "Amplitude [$units]", "Reference channel$(_pl(length(ch))) $(_channel2channel_name(ch))\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        end
        if ch_t[ch[1]] == "eog"
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "Amplitude [$units]", "EOG channel$(_pl(length(ch))) $(_channel2channel_name(ch))\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        end
        if ch_t[ch[1]] == "emg"
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "Amplitude [$units]", "EMG channel$(_pl(length(ch))) $(_channel2channel_name(ch))\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        end
        if ch_t[ch[1]] == "ecg"
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "Amplitude [$units]", "ECG channel\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        end
        if ch_t[ch[1]] == "other"
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "Amplitude [$units]", "Other channel$(_pl(length(ch))) $(_channel2channel_name(ch))\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        end
        if ch_t[ch[1]] == "mrk"
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "Amplitude [$units]", "Marker$(_pl(length(ch))) $(_channel2channel_name(ch))\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        end
        p = plot_signal_butterfly(t,
                                  s[ch, :],
                                  clabels=clabels[ch],
                                  xlabel=xl,
                                  ylabel=yl,
                                  title=tt,
                                  scale=scale,
                                  units=units,
                                  norm=norm,
                                  mono=mono;
                                  kwargs...)
    end
    
    if type === :mean
        length(ch_t_uni) > 1 && throw(ArgumentError("For type=:mean plot all channels should be of the same type."))
        size(s, 1) == 1 && throw(ArgumentError("For type=:mean plot the signal must contain ≥ 2 channels."))

        units = _set_units(obj, ch[1])
        if ch_t[ch[1]] == "eeg"
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "Amplitude [$units]", "Averaged EEG channel$(_pl(length(ch))) $(_channel2channel_name(ch)) [mean ± 95%CI]\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        end
        if ch_t[ch[1]] == "grad"
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "Amplitude [$units]", "Averaged gradiometer$(_pl(length(ch))) $(_channel2channel_name(ch)) [mean ± 95%CI]\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        end
        if ch_t[ch[1]] == "mag"
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "Amplitude [$units]", "Averaged magnetometer$(_pl(length(ch))) $(_channel2channel_name(ch)) [mean ± 95%CI]\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        end
        if ch_t[ch[1]] == "nirs_int"
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "Intensity [$units]", "Averaged NIRS intensity channel$(_pl(length(ch))) $(_channel2channel_name(ch)) [mean ± 95%CI]\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        end
        if ch_t[ch[1]] == "nirs_od"
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "OD [$units]", "Averaged NIRS optical density channel$(_pl(length(ch))) $(_channel2channel_name(ch)) [mean ± 95%CI]\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        end
        if ch_t[ch[1]] == "nirs_hbo"
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "HbO concentration [$units]", "Averaged NIRS HbO concentration channel$(_pl(length(ch))) $(_channel2channel_name(ch)) [mean ± 95%CI]\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        end
        if ch_t[ch[1]] == "nirs_hbr"
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "HbR concentration [$units]", "Averaged NIRS HbR concentration channel$(_pl(length(ch))) $(_channel2channel_name(ch)) [mean ± 95%CI]\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        end
        if ch_t[ch[1]] == "nirs_hbt"
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "HbT concentration [$units]", "Averaged NIRS HbT concentration channel$(_pl(length(ch))) $(_channel2channel_name(ch)) [mean ± 95%CI]\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        end
        if ch_t[ch[1]] == "nirs_aux"
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "Amplitude [$units]", "Averaged NIRS AUX channel$(_pl(length(ch))) $(_channel2channel_name(ch)) [mean ± 95%CI]\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        end
        if ch_t[ch[1]] == "ref"
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "Amplitude [$units]", "Averaged Reference channel$(_pl(length(ch))) $(_channel2channel_name(ch)) [mean ± 95%CI]\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        end
        if ch_t[ch[1]] == "eog"
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "Amplitude [$units]", "Averaged EOG channel$(_pl(length(ch))) $(_channel2channel_name(ch)) [mean ± 95%CI]\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        end
        if ch_t[ch[1]] == "emg"
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "Amplitude [$units]", "Averaged EMG channel$(_pl(length(ch))) $(_channel2channel_name(ch)) [mean ± 95%CI]\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        end
        if ch_t[ch[1]] == "other"
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "Amplitude [$units]", "Averaged  other channel$(_pl(length(ch))) $(_channel2channel_name(ch)) [mean ± 95%CI]\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        end
        if ch_t[ch[1]] == "mrk"
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "Amplitude [$units]", "Averaged  marker$(_pl(length(ch))) $(_channel2channel_name(ch)) [mean ± 95%CI]\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        end
        p = plot_signal_avg(t,
                            s[ch, :],
                            xlabel=xl,
                            ylabel=yl,
                            title=tt,
                            scale=scale,
                            units=units,
                            norm=norm,
                            mono=mono;
                            kwargs...)
    end

    # add epochs markers
    # TODO: draw epoch numbers
    if emarkers == true
        if length(p) > 1
            for p_idx in 1:length(p)
                p[p_idx] = Plots.vline!(p[p_idx],
                                        epoch_markers,
                                        linestyle=:dash,
                                        linewidth=0.5,
                                        linecolor=:blue,
                                        label="")
            end
        else
            p = Plots.vline!(epoch_markers,
                             linestyle=:dash,
                             linewidth=0.5,
                             linecolor=:blue,
                             label="")
        end
    end

    # plot markers if available
    # TODO: draw markers length
    if markers == true && _has_markers(obj) == true
        markers_pos = obj.markers[!, :start]
        markers_id = obj.markers[!, :id]
        markers_desc = obj.markers[!, :description]
        if length(p) > 1
            for p_idx in 1:length(p)
                p[p_idx] = Plots.vline!(p[p_idx],
                                        markers_pos,
                                        linestyle=:dash,
                                        linewidth=0.5,
                                        linecolor=:black,
                                        label=false)
                for idx in eachindex(markers_desc)
                    p[p_idx] = Plots.plot!(p[p_idx], annotation=(markers_pos[idx], -0.92, Plots.text("$(markers_id[idx])/$(markers_desc[idx])", pointsize=5, halign=:left, valign=:top, rotation=90)), label=false)
                end
            end
        else
            p = Plots.vline!(markers_pos,
                             linestyle=:dash,
                             linewidth=0.5,
                             linecolor=:black,
                             label=false)
            for idx in eachindex(markers_desc)
                p = Plots.plot!(annotation=(markers_pos[idx], -0.92, Plots.text("$(markers_id[idx])/$(markers_desc[idx])", pointsize=5, halign=:left, valign=:top, rotation=90)), label=false)
            end
        end
    end

    if length(p) > 1
        if obj.header.recording[:data_type] == "eeg"
            h_primary = (0.75 + 0.05 * (5 - length(p)))
            h_secondary = 1.0 - h_primary
            h = vcat(h_primary, repeat([h_secondary / (length(ch_t_uni) - 1)], (length(ch_t_uni) - 1)))
            p = Plots.plot!(p..., layout=grid(length(ch_t_uni), 1, heights=h); kwargs...)
        end
        if obj.header.recording[:data_type] == "meg"
            h_primary = (0.75 + 0.05 * (5 - length(p)))
            h_secondary = 1.0 - h_primary
            h = vcat(h_primary, repeat([h_secondary / (length(ch_t_uni) - 1)], (length(ch_t_uni) - 1)))
            p = Plots.plot!(p..., layout=grid(length(ch_t_uni), 1, heights=h); kwargs...)
        end
        if obj.header.recording[:data_type] == "nirs"
            p = Plots.plot!(p..., layout=(length(ch_t_uni), 1); kwargs...)
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
- `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
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
function plot(obj::NeuroAnalyzer.NEURO, c::Union{Symbol, AbstractArray}; ep::Union{Int64, AbstractRange}=0, c_idx::Union{Int64, Vector{Int64}, <:AbstractRange}=0, seg::Tuple{Int64, Int64}=(0, 10), xlabel::String="default", ylabel::String="default", title::String="default", mono::Bool=false, emarkers::Bool=true, markers::Bool=true, scale::Bool=true, units::String="a.u.", type::Symbol=:normal, norm::Bool=false, kwargs...)

    if signal_len(obj) < 10 * sr(obj) && seg == (0, 10)
        seg = (0, obj.time_pts[end])
    else
        _check_segment(obj, seg)
    end
    seg = (vsearch(seg[1], obj.time_pts), vsearch(seg[2], obj.time_pts))

    _check_var(type, [:normal, :butterfly, :mean, :stack], "type")

    if ep != 0
        _check_epochs(obj, ep)
        seg = (((ep[1] - 1) * epoch_len(obj) + 1), seg[2])
        if ep isa Int64
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
    c_name = ""
    if typeof(c) == Symbol
        c_name = string(c)
        c = _get_component(obj, c)
    end
    c_idx == 0 && (c_idx = _select_cidx(c, c_idx))
    _check_cidx(c, c_idx)
    if size(c, 1) == 1
        clabels = c_name
    else
        clabels = _gen_clabels(c)[c_idx]
        clabels = c_name .* clabels
    end
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
        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "", "Component$(_pl(length(c_idx))) $(_channel2channel_name(c_idx)) amplitude\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        p = plot_signal(t,
                        s,
                        clabels=clabels,
                        xlabel=xl,
                        ylabel=yl,
                        title=tt,
                        scale=scale,
                        units=units,
                        mono=mono;
                        kwargs...)
    elseif type === :butterfly
        size(s, 1) == 1 && throw(ArgumentError("For type=:butterfly plot the signal must contain ≥ 2 channels."))
        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "Amplitude [$units]", "Components $(_channel2channel_name(c_idx)) amplitude\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        p = plot_signal_butterfly(t,
                                  s,
                                  clabels=clabels,
                                  xlabel=xl,
                                  ylabel=yl,
                                  title=tt,
                                  scale=scale,
                                  units=units,
                                  norm=norm,
                                  mono=mono;
                                  kwargs...)
    elseif type === :mean
        size(s, 1) == 1 && throw(ArgumentError("For type=:mean plot the signal must contain ≥ 2 channels."))
        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "Amplitude [$units]", "Averaged components $(_channel2channel_name(c_idx)) amplitude [mean ± 95%CI]\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        p = plot_signal_avg(t,
                            s,
                            clabels=clabels,
                            xlabel=xl,
                            ylabel=yl,
                            title=tt,
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
        markers_pos = obj.markers[!, :start]
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