export plot_signal
export plot_signal_avg
export plot_signal_butterfly
export plot_2signals
export plot

"""
    plot_signal(t, s; <keyword arguments>)

Plot amplitude of single-channel signal.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
- `s::AbstractVector`: data to plot
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `bad::Bool=false`: is this a bad channel
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_signal(t::Union{AbstractVector, AbstractRange}, s::AbstractVector; xlabel::String="", ylabel::String="", title::String="", bad::Bool=false, kwargs...)::Plots.Plot{Plots.GRBackend}

    # prepare plot
    plot_size = (1200, 400)
    p = Plots.plot(xlabel=xlabel,
                   ylabel=ylabel,
                   xlims=_xlims(t),
                   xticks=_ticks(t),
                   ytick_direction=:out,
                   xtick_direction=:out,
                   ylims=_ylims(s),
                   title=title,
                   size=plot_size,
                   margins=20Plots.px,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=6,
                   ytickfontsize=6;
                   kwargs...)

    # plot signal
    if bad
        p = Plots.plot!(t,
                        s,
                        linewidth=1,
                        label="",
                        alpha=0.2,
                        color=:black)
    else
        p = Plots.plot!(t,
                        s,
                        linewidth=1,
                        label="",
                        color=:black)
    end

    return p

end

"""
    plot_signal(t, s; <keyword arguments>)

Plot amplitude of multi-channel signal.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
- `s::AbstractArray`: data to plot
- `clabels::Vector{String}=repeat([""], size(s, 1))`: channel labels
- `ctypes:::Vector{String}=repeat([""], size(s, 1))`: channel types
- `cunits::Vector{String}=repeat([""], size(s, 1))`: channel units
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `scale::Bool=true`: draw scale
- `bad::Vector{Bool}=zeros(Bool, size(s, 1))`: list of bad channels
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_signal(t::Union{AbstractVector, AbstractRange}, s::AbstractArray; clabels::Vector{String}=repeat([""], size(s, 1)), ctypes::Vector{String}=repeat([""], size(s, 1)), cunits::Vector{String}=repeat([""], size(s, 1)), xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, scale::Bool=true, bad::Vector{Bool}=zeros(Bool, size(s, 1)), kwargs...)::Plots.Plot{Plots.GRBackend}

    ch_n = size(s, 1)

    ctypes_uni = unique(ctypes)
    t_pos = zeros(Int64, length(ctypes_uni))
    for idx in eachindex(ctypes_uni)
         t_pos[idx] = findfirst(isequal(ctypes_uni[idx]), ctypes)
    end
    ctypes_uni_pos = zeros(Int64, length(ctypes))
    ctypes_uni_pos[t_pos] .= 1

    # set colors
    if mono
        pal = :grays
        channel_color = repeat([:black], ch_n)
    else
        pal = :darktest
        channel_color = zeros(Int64, ch_n)
        for idx in eachindex(ctypes_uni)
            channel_color[ctypes .== ctypes_uni[idx]] .= idx
        end
    end

    # get ranges of the original signal for the scales
    # normalize and shift so all channels are visible
    r = Float64[]
    @inbounds for ch_idx in eachindex(ctypes_uni)
        push!(r, round(_get_range(s[ctypes .== ctypes_uni[ch_idx], :])))
        s[ctypes .== ctypes_uni[ch_idx], :] = @views normalize(s[ctypes .== ctypes_uni[ch_idx], :], method=:minmax)
    end
    # reverse so 1st channel is on top
    s = @views reverse(s[:, eachindex(t)], dims = 1)
    channel_color = reverse(channel_color)
    bad = reverse(bad)
    # each channel is between -1.0 and +1.0
    # scale by 0.5 so maxima do not overlap
    @inbounds for idx in 1:ch_n
        s[idx, :] = @views (s[idx, :] .* 0.5) .+ (idx - 1)
    end

    # prepare plot
    plot_size = (1200, 100 + 40 * ch_n)
    p = Plots.plot(ylabel=ylabel,
                   xlims=_xlims(t),
                   xticks=(_ticks(t), []),
                   ylims=(-1, ch_n),
                   yticks=((ch_n - 1):-1:0, []),
                   ytick_direction=:none,
                   xtick_direction=:out,
                   title=title,
                   palette=pal,
                   size=plot_size,
                   top_margin=10Plots.px,
                   bottom_margin=40Plots.px,
                   right_margin=20Plots.px,
                   left_margin=60Plots.px,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=6,
                   ytickfontsize=6;
                   kwargs...)

    # plot channels
    for idx in 1:ch_n
        if !bad[idx]
            p = @views Plots.plot!(t,
                                   s[idx, :],
                                   linewidth=0.75,
                                   label="",
                                   color=channel_color[idx])
        else
            p = @views Plots.plot!(t,
                                   s[idx, :],
                                   linewidth=0.75,
                                   alpha=0.2,
                                   label="",
                                   color=channel_color[idx])
        end
    end

    # draw labels
    for idx in 1:ch_n
        s_pos = ch_n - idx
        p = Plots.plot!(annotations=(_xlims(t)[1], (s_pos), Plots.text("$(clabels[idx])   ", pointsize=8, halign=:right, valign=:center)), label=false)
    end

    # draw ticks
    xt = collect(_ticks(t))
    xt_s = string.(xt)
    for idx in eachindex(xt)
        p = Plots.plot!(annotations=(xt[idx], (-1.5), Plots.text("$(xt_s[idx])", pointsize=6, halign=:center, valign=:center)), label=false)
    end
    p = Plots.plot!(annotations=(mean(xt), (-2), Plots.text("$xlabel", pointsize=8, halign=:center, valign=:center)), label=false)

    # draw scales
    if scale
        idx2 = 1
        for idx1 in 1:ch_n
            if ctypes_uni_pos[idx1] == 1
                s_pos = ch_n - idx1 + 1
                p = Plots.plot!([_xlims(t)[1], _xlims(t)[1]], [(s_pos - 1.5), (s_pos - 0.5)], color=:red, linewidth=5, label="")
                p = Plots.plot!(annotations=(_xlims(t)[1], (s_pos - 0.5), Plots.text("$(r[idx2]) $(cunits[idx1])  ", pointsize=5, halign=:right, valign=:bottom, rotation=90)), label=false)
                idx2 += 1
            end
        end
    end

    return p

end

"""
    plot_signal(t, s1, s2; <keyword arguments>)

Plot amplitude of single-channel signal.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
- `s1::AbstractVector`: data to plot
- `s2::AbstractVector`: data to plot
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_signal(t::Union{AbstractVector, AbstractRange}, s1::AbstractVector, s2::AbstractVector; xlabel::String="", ylabel::String="", title::String="", kwargs...)::Plots.Plot{Plots.GRBackend}

    @assert length(s1) == length(s2) "s1 and s2 must have the same length."

    # prepare plot
    plot_size = (1200, 400)
    p = Plots.plot(xlabel=xlabel,
                   ylabel=ylabel,
                   xlims=_xlims(t),
                   xticks=_ticks(t),
                   ytick_direction=:out,
                   xtick_direction=:out,
                   ylims=_ylims(s1)[1] < _ylims(s2)[1] ? _ylims(s1) : _ylims(s2),
                   title=title,
                   size=plot_size,
                   margins=20Plots.px,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=6,
                   ytickfontsize=6;
                   kwargs...)

    # plot signals
    p = Plots.plot!(t,
                    s1,
                    linewidth=1,
                    label="",
                    alpha=0.5,
                    color=:black)
    p = Plots.plot!(t,
                    s2,
                    linewidth=1,
                    label="",
                    alpha=0.5,
                    color=:blue)

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
- `mono::Bool=false`: use color or gray palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_signal_avg(t::Union{AbstractVector, AbstractRange}, s::AbstractArray; xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)::Plots.Plot{Plots.GRBackend}

    pal = mono ? :grays : :darktest

    # get mean and 95%CI
    s_m, _, s_u, s_l = msci95(s)

    # get limits
    ylim = (round(minimum(s_l) * 1.5, digits=0), round(maximum(s_u) * 1.5, digits=0))
    ylim = _tuple_max(ylim)

    # prepare plot
    p = Plots.plot(xlabel=xlabel,
                   ylabel=ylabel,
                   xlims=_xlims(t),
                   xticks=_ticks(t),
                   ylims=ylim,
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
- `avg::Bool=false`: plot average channels
- `mono::Bool=false`: use color or gray palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_signal_butterfly(t::Union{AbstractVector, AbstractRange}, s::AbstractArray; clabels::Vector{String}=[""], xlabel::String="", ylabel::String="", title::String="", avg::Bool=true, mono::Bool=false, kwargs...)::Plots.Plot{Plots.GRBackend}

    pal = mono ? :grays : :darktest

    ch_n = size(s, 1)

    # get limits
    ylim = (round(minimum(s) * 1.5, digits=0), round(maximum(s) * 1.5, digits=0))
    ylim = _tuple_max(ylim)

    # channel labels
    clabels == [""] && (clabels = repeat([""], ch_n))

    # plot channels
    plot_size = (1200, 500)
    p = Plots.plot(xlabel=xlabel,
                   ylabel=ylabel,
                   xlims=_xlims(t),
                   xticks=_ticks(t),
                   ylims=ylim,
                   title=title,
                   palette=pal,
                   size=plot_size,
                   legend=false,
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
                        label=clabels[idx])
    end

    # plot averaged channels
    if avg
        s = mean(s, dims=1)[:]
        p = Plots.plot!(t,
                        s,
                        linewidth=2,
                        linecolor=:black,
                        label=false)
    end

    return p

end

"""
    plot_signal(t, s1, s2; <keyword arguments>)

Plot amplitude of multi-channel signals.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
- `s1::AbstractArray`: data to plot
- `s2::AbstractArray`: data to plot
- `clabels::Vector{String}=repeat([""], size(s1, 1))`: channel labels
- `ctypes:::Vector{String}=repeat([""], size(s1, 1))`: channel types
- `cunits::Vector{String}=repeat([""], size(s1, 1))`: channel units
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `scale::Bool=true`: draw scale
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_signal(t::Union{AbstractVector, AbstractRange}, s1::AbstractArray, s2::AbstractArray; clabels::Vector{String}=repeat([""], size(s1, 1)), ctypes::Vector{String}=repeat([""], size(s1, 1)), cunits::Vector{String}=repeat([""], size(s1, 1)), xlabel::String="", ylabel::String="", title::String="", scale::Bool=true, kwargs...)::Plots.Plot{Plots.GRBackend}

    ch_n = size(s1, 1)

    ctypes_uni = unique(ctypes)
    t_pos = zeros(Int64, length(ctypes_uni))
    for idx in eachindex(ctypes_uni)
         t_pos[idx] = findfirst(isequal(ctypes_uni[idx]), ctypes)
    end
    ctypes_uni_pos = zeros(Int64, length(ctypes))
    ctypes_uni_pos[t_pos] .= 1

    # get ranges of the original signal for the scales
    # normalize and shift so all channels are visible
    r = Float64[]
    @inbounds for ch_idx in eachindex(ctypes_uni)
        push!(r, round(_get_range(s1[ctypes .== ctypes_uni[ch_idx], :])))
        s1[ctypes .== ctypes_uni[ch_idx], :] = @views normalize(s1[ctypes .== ctypes_uni[ch_idx], :], method=:minmax)
        s2[ctypes .== ctypes_uni[ch_idx], :] = @views normalize(s2[ctypes .== ctypes_uni[ch_idx], :], method=:minmax)
    end
    # reverse so 1st channel is on top
    s1 = @views reverse(s1[:, eachindex(t)], dims = 1)
    s2 = @views reverse(s2[:, eachindex(t)], dims = 1)
    # each channel is between -1.0 and +1.0
    # scale by 0.5 so maxima do not overlap
    @inbounds for idx in 1:ch_n
        s1[idx, :] = @views (s1[idx, :] .* 0.5) .+ (idx - 1)
        s2[idx, :] = @views (s2[idx, :] .* 0.5) .+ (idx - 1)
    end

    # prepare plot
    plot_size = (1200, 100 + 40 * ch_n)
    p = Plots.plot(ylabel=ylabel,
                   xlims=_xlims(t),
                   xticks=(_ticks(t), []),
                   ylims=(-1, ch_n),
                   yticks=((ch_n - 1):-1:0, []),
                   ytick_direction=:none,
                   xtick_direction=:out,
                   title=title,
                   palette=:darktest,
                   size=plot_size,
                   top_margin=10Plots.px,
                   bottom_margin=40Plots.px,
                   right_margin=20Plots.px,
                   left_margin=60Plots.px,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=6,
                   ytickfontsize=6;
                   kwargs...)

    # plot channels
    for idx in 1:ch_n
        p = @views Plots.plot!(t,
                               s1[idx, :],
                               linewidth=1,
                               label="",
                               color=:black,
                               alpha=0.5)
        p = @views Plots.plot!(t,
                               s2[idx, :],
                               linewidth=1,
                               label="",
                               color=:blue,
                               alpha=0.5)
    end

    # draw labels
    for idx in 1:ch_n
        s_pos = ch_n - idx
        p = Plots.plot!(annotations=(_xlims(t)[1], (s_pos), Plots.text("$(clabels[idx])  ", pointsize=8, halign=:right, valign=:center)), label=false)
    end

    # draw ticks
    xt = collect(_ticks(t))
    xt_s = string.(xt)
    for idx in eachindex(xt)
        p = Plots.plot!(annotations=(xt[idx], (-1.5), Plots.text("$(xt_s[idx])", pointsize=6, halign=:center, valign=:center)), label=false)
    end
    p = Plots.plot!(annotations=(mean(xt), (-2), Plots.text("$xlabel", pointsize=8, halign=:center, valign=:center)), label=false)

    # draw scales
    if scale
        idx2 = 1
        for idx1 in 1:ch_n
            if ctypes_uni_pos[idx1] == 1
                s_pos = ch_n - idx1 + 1
                p = Plots.plot!([_xlims(t)[1], _xlims(t)[1]], [(s_pos - 1.5), (s_pos - 0.5)], color=:red, linewidth=5, label="")
                p = Plots.plot!(annotations=(_xlims(t)[1], (s_pos - 1.5), Plots.text("$(r[idx2]) $(cunits[idx1])  ", pointsize=5, halign=:right, valign=:bottom)), label=false)
                idx2 += 1
            end
        end
    end

    return p

end

"""
    plot(obj; <keyword arguments>)

Plot signal.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ep::Union{Int64, AbstractRange}=0`: epoch to display
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
- `xlabel::String="default"`: x-axis label, default is Time [s]
- `ylabel::String="default"`: y-axis label, default is no label
- `title::String="default"`: plot title
- `mono::Bool=false`: use color or gray palette
- `emarkers::Bool`: draw epoch markers if available
- `markers::Bool`: draw markers if available
- `scale::Bool=true`: draw scale
- `type::Symbol=:normal`: plot type:
    - `:normal`
    - `:mean`: mean ± 95%CI
    - `:butterfly`: butterfly plot
- `avg::Bool=false`: plot average EDA
- `bad::Bool=false`: plot bad channels
- `s_pos::Tuple{Real, Real}=(0, 0)`: draw segment borders if different than (0, 0), used by `iedit()`
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot(obj::NeuroAnalyzer.NEURO; ep::Union{Int64, AbstractRange}=0, ch::Union{String, Vector{String}, Regex}, seg::Tuple{Real, Real}=(0, 10), xlabel::String="default", ylabel::String="default", title::String="default", mono::Bool=false, emarkers::Bool=true, markers::Bool=true, scale::Bool=true, type::Symbol=:normal, avg::Bool=true, bad::Bool=true, s_pos::Tuple{Real, Real}=(0, 0), kwargs...)::Plots.Plot{Plots.GRBackend}

    datatype(obj) == "erp" && _warn("For ERP objects, use plot_erp()")
    datatype(obj) == "erf" && _warn("For ERF objects, use plot_erp()")
    datatype(obj) == "mep" && _warn("For MEP objects, use plot_mep()")

    if signal_len(obj) <= 10 * sr(obj) && seg == (0, 10)
        seg = (obj.time_pts[1], obj.time_pts[end])
    else
        _check_segment(obj, seg)
    end
    seg = (vsearch(seg[1], obj.time_pts), vsearch(seg[2], obj.time_pts))

    _check_var(type, [:normal, :butterfly, :mean], "type")

    if ep != 0
        _check_epochs(obj, ep)
        if nepochs(obj) == 1
            ep = 1
        else
            if ep isa Int64
                seg = (((ep - 1) * epoch_len(obj) + 1), (ep * epoch_len(obj)))
            else
                seg = (((ep[1] - 1) * epoch_len(obj) + 1), (ep[end] * epoch_len(obj)))
            end
        end
    end

    # do not show epoch markers if there are no epochs
    nepochs(obj) == 1 && (emarkers = false)
    if emarkers
        epoch_markers = _get_epoch_markers(obj)
    end

    # check channels
    ch = get_channel(obj, ch=ch, exclude="")

    # get time vector
    if seg[2] <= epoch_len(obj)
        s = obj.data[:, seg[1]:seg[2], 1]
    else
        s = epoch(obj, ep_n=1).data[:, seg[1]:seg[2], 1]
    end
    t = obj.time_pts[seg[1]:seg[2]]

    _, t_s1, _, t_s2 = _convert_t(t[1], t[end])
    ep = _s2epoch(obj, seg[1], seg[2])

    length(ch) == 1 && (ch = ch[1])

    xl, yl, tt = "", "", ""

    ch_init = ch
    bm = obj.header.recording[:bad_channel]
    ctypes = obj.header.recording[:channel_type]
    clabels = labels(obj)
    # sort channels by their type
    if !isa(ch, Int64)
        ctypes = ctypes[ch]
        clabels = clabels[ch]
        cunits = obj.header.recording[:unit][ch]
        if bad
            bm = bm[ch, ep]
            if isa(bm, Matrix{Bool})
                bm_tmp = zeros(Bool, size(bm, 1))
                [sum(bm[idx, :]) > 0 && (bm_tmp[idx] = true) for idx in 1:size(bm, 1)]
                bm = bm_tmp
            end
        else
            bm = zeros(Bool, length(ch))
        end
    else
        if bad
            bm = bm[ch]
        else
            bm = false
        end
    end

    if type === :normal
        if isa(ch, Int64)
            ch_name = _ch_rename(ctypes[ch])
            xl, yl, tt = _set_defaults(xlabel,
                                       ylabel,
                                       title,
                                       "Time [s]",
                                       "",
                                       "Channel $ch: $(clabels[ch]) ($ch_name)")
            if datatype(obj) == "eda"
                ylabel == "default" && (yl = "Impedance [μS]")
                p = plot_eda(t,
                             s[ch, :],
                             xlabel=xl,
                             ylabel=yl,
                             title=tt,
                             mono=mono;
                             kwargs...)
            else
                ylabel == "default" && (yl = "Amplitude [$(_ch_units(obj, clabels[ch]))]")
                p = plot_signal(t,
                                s[ch, :],
                                xlabel=xl,
                                ylabel=yl,
                                title=tt,
                                bad=bm,
                                mono=mono;
                                kwargs...)
            end
        else
            xl, yl, tt = _set_defaults(xlabel,
                                       ylabel,
                                       title,
                                       "Time [s]",
                                       "",
                                       "")
            p = plot_signal(t,
                            s[ch, :],
                            ctypes=ctypes,
                            clabels=clabels,
                            cunits=cunits,
                            xlabel=xl,
                            ylabel=yl,
                            title=tt,
                            bad=bm,
                            scale=scale,
                            mono=mono;
                            kwargs...)
        end
    end

    if type === :butterfly
        @assert length(unique(ctypes)) == 1 "For plot type=:butterfly all channels must be of the same type."
        @assert size(s, 1) >= 2 "For plot type=:butterfly the signal must contain ≥ 2 channels."
        xl, yl, tt = _set_defaults(xlabel,
                                   ylabel,
                                   title,
                                   "Time [s]",
                                   "Amplitude [$(obj.header.recording[:unit][ch[1]])]",
                                   "$(size(s[ch, :], 1)) $(uppercase(unique(ctypes)[1])) channels")
        if datatype(obj) == "eda"
            (datatype(obj) == "eda" && ylabel == "default") && (yl = "Impedance [μS]")
            p = plot_eda_butterfly(t,
                                   s[ch, :],
                                   clabels=clabels,
                                   xlabel=xl,
                                   ylabel=yl,
                                   title=tt,
                                   avg=avg,
                                   mono=mono;
                                   kwargs...)
        else
            p = plot_signal_butterfly(t,
                                      s[ch, :],
                                      clabels=clabels,
                                      xlabel=xl,
                                      ylabel=yl,
                                      title=tt,
                                      avg=avg,
                                      mono=mono;
                                      kwargs...)
        end
    end

    if type === :mean
        @assert length(unique(ctypes)) == 1 "For plot type=:mean all channels must be of the same type."
        @assert size(s, 1) >= 2 "For plot type=:mean the signal must contain ≥ 2 channels."
        xl, yl, tt = _set_defaults(xlabel,
                                   ylabel,
                                   title,
                                   "Time [s]",
                                   "Amplitude [$(obj.header.recording[:unit][ch[1]])]",
                                   "$(size(s[ch, :], 1)) $(uppercase(unique(ctypes)[1])) channels")
        if datatype(obj) == "eda"
            (datatype(obj) == "eda" && ylabel == "default") && (yl = "Impedance [μS]")
            p = plot_eda_avg(t,
                             s[ch, :],
                             xlabel=xl,
                             ylabel=yl,
                             title=tt,
                             mono=mono;
                             kwargs...)
        else
            p = plot_signal_avg(t,
                                s[ch, :],
                                xlabel=xl,
                                ylabel=yl,
                                title=tt,
                                mono=mono;
                                kwargs...)
        end
    end

    # add epochs markers
    # TODO: draw epoch numbers
    if emarkers
        p = Plots.vline!(epoch_markers,
                         linestyle=:dot,
                         linewidth=0.5,
                         linecolor=:blue,
                         label="")
    end

    # plot markers if available
    # TODO: draw markers length
    if markers && _has_markers(obj)
        markers_pos = obj.markers[!, :start]
        markers_id = obj.markers[!, :id]
        markers_desc = obj.markers[!, :value]
        p = Plots.vline!(p,
                         markers_pos,
                         linestyle=:dash,
                         linewidth=1,
                         linecolor=:black,
                         label=false)
        for idx in eachindex(markers_desc)
            p = Plots.plot!(p, annotations=(markers_pos[idx] + 0.1, -0.92, Plots.text("$(markers_id[idx]) / $(markers_desc[idx])", pointsize=5, halign=:left, valign=:top, rotation=90)), label=false)
        end
    end

    # draw segment borders
    p = Plots.vline!(p,
                     [s_pos[1]],
                     color=:black,
                     lw=1,
                     labels="")
    p = Plots.vline!(p,
                     [s_pos[2]],
                     color=:black,
                     lw=1,
                     labels="")

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
- `c_idx::Union{Int64, Vector{Int64}, AbstractRange}=0`: component channel to display, default is all component channels
- `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
- `xlabel::String="default"`: x-axis label, default is Time [s]
- `ylabel::String="default"`: y-axis label, default is no label
- `title::String="default"`: plot title
- `mono::Bool=false`: use color or gray palette
- `emarkers::Bool`: draw epoch markers if available
- `scale::Bool=true`: draw scale
- `type::Symbol=:normal`: plot type:
    - `:normal`
    - `:mean`: mean ± 95%CI
    - `:butterfly`: butterfly plot
- `avg::Bool=false`: plot average EDA
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot(obj::NeuroAnalyzer.NEURO, c::Union{Symbol, AbstractArray}; ep::Union{Int64, AbstractRange}=0, c_idx::Union{Int64, Vector{Int64}, AbstractRange}=0, seg::Tuple{Real, Real}=(0, 10), xlabel::String="default", ylabel::String="default", title::String="default", mono::Bool=false, emarkers::Bool=true, scale::Bool=true, type::Symbol=:normal, avg::Bool=true, kwargs...)::Plots.Plot{Plots.GRBackend}

    if signal_len(obj) < 10 * sr(obj) && seg == (0, 10)
        seg = (0, obj.time_pts[end])
    else
        _check_segment(obj, seg)
    end
    seg = (vsearch(seg[1], obj.time_pts), vsearch(seg[2], obj.time_pts))

    _check_var(type, [:normal, :butterfly, :mean], "type")

    if ep != 0
        _check_epochs(obj, ep)
        if nepochs(obj) == 1
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
    nepochs(obj) == 1 && (emarkers = false)
    if emarkers
        epoch_markers = _get_epoch_markers(obj)
    end

    # select component channels, default is all channels
    c_name = ""
    if c isa Symbol
        c_name = string(c)
        c = _get_component(obj, c)
    end
    c_idx == 0 && (c_idx = _select_cidx(c, c_idx))
    (c_idx isa(Vector{Int64}) && length(c_idx) == 1) && (c_idx = c_idx[1])
    _check_cidx(c, c_idx)
    if size(c, 1) == 1
        clabels = c_name
    else
        clabels = _gen_clabels(c)[c_idx]
        clabels = c_name .* clabels
    end
    c_idx isa Int64 && (clabels = [clabels])

    # get time vector
    if seg[2] <= epoch_len(obj)
        s = c[c_idx, seg[1]:seg[2], 1]
    else
        s = _make_epochs(c, ep_n=1)[c_idx, seg[1]:seg[2], 1]
    end
    t = obj.time_pts[seg[1]:seg[2]]

    _, t_s1, _, t_s2 = _convert_t(t[1], t[end])
    ep = _s2epoch(obj, seg[1], seg[2])

    if type === :normal
        if isa(c_idx, Int64)
            xl, yl, tt = _set_defaults(xlabel,
                                       ylabel,
                                       title,
                                       "Time [s]",
                                       "Amplitude",
                                       "Component$(_pl(length(c_idx))) $(c_idx) amplitude\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
                                       #"Component$(_pl(length(c_idx))) $(_channel2channel_name(c_idx)) amplitude\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
            p = plot_signal(t,
                            s,
                            xlabel=xl,
                            ylabel=yl,
                            title=tt,
                            mono=mono;
                            kwargs...)
        else
            xl, yl, tt = _set_defaults(xlabel,
                                       ylabel,
                                       title,
                                       "Time [s]",
                                       "",
                                       "")
            p = plot_signal(t,
                            s,
                            clabels=clabels,
                            xlabel=xl,
                            ylabel=yl,
                            title=tt,
                            scale=scale,
                            mono=mono;
                            kwargs...)
        end
    end

    if type === :butterfly
        @assert size(s, 1) >= 2 "For type=:butterfly plot the signal must contain ≥ 2 channels."
        xl, yl, tt = _set_defaults(xlabel,
                                   ylabel,
                                   title,
                                   "Time [s]",
                                   "Amplitude [$units]",
                                   "Components $(_channel2channel_name(c_idx)) amplitude\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        p = plot_signal_butterfly(t,
                                  s,
                                  clabels=clabels,
                                  xlabel=xl,
                                  ylabel=yl,
                                  title=tt,
                                  avg=avg,
                                  mono=mono;
                                  kwargs...)
    end

    if type === :mean
        @assert size(s, 1) >= 2 "For type=:mean plot the signal must contain ≥ 2 channels."
        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "Amplitude [$units]", "Averaged components $(_channel2channel_name(c_idx)) amplitude [mean ± 95%CI]\n[epoch$(_pl(length(ep))): $ep, time window: $t_s1:$t_s2]")
        p = plot_signal_avg(t,
                            s,
                            xlabel=xl,
                            ylabel=yl,
                            title=tt,
                            mono=mono;
                            kwargs...)
    end

    # add epochs markers
    # TODO: draw epoch numbers
    if emarkers
        p = Plots.vline!(epoch_markers,
                         linestyle=:dash,
                         linewidth=0.5,
                         linecolor=:blue,
                         label="")
    end

    Plots.plot(p)

    return p

end

"""
    plot(obj1, obj2; <keyword arguments>)

Plot signal.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `obj2::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ep::Union{Int64, AbstractRange}=0`: epoch to display
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
- `xlabel::String="default"`: x-axis label, default is Time [s]
- `ylabel::String="default"`: y-axis label, default is no label
- `title::String="default"`: plot title
- `scale::Bool=true`: draw scale
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ep::Union{Int64, AbstractRange}=0, ch::Union{String, Vector{String}, Regex}, seg::Tuple{Real, Real}=(0, 10), xlabel::String="default", ylabel::String="default", title::String="default", scale::Bool=true, kwargs...)::Plots.Plot{Plots.GRBackend}

    @assert sr(obj1) == sr(obj2) "OBJ1 and OBJ2 must have the same sampling rate."
    @assert size(obj1.data) == size(obj2.data) "Signals of OBJ1 and OBJ2 must have the same size."
    @assert datatype(obj1) == obj2.header.recording[:data_type] "OBJ1 and OBJ2 must have the same data type."

    if signal_len(obj1) <= 10 * sr(obj1) && seg == (0, 10)
        seg = (obj1.time_pts[1], obj1.time_pts[end])
    else
        _check_segment(obj1, seg)
    end
    seg = (vsearch(seg[1], obj1.time_pts), vsearch(seg[2], obj1.time_pts))

    if ep != 0
        _check_epochs(obj1, ep)
        if nepochs(obj1) == 1
            ep = 0
        else
            seg = (((ep[1] - 1) * epoch_len(obj1) + 1), seg[2])
            if ep isa Int64
                seg = (seg[1], (seg[1] + epoch_len(obj1) - 1))
            else
                seg = (seg[1], (ep[end] * epoch_len(obj1)))
            end
            ep = 0
        end
    end

    # check channels
    _ = get_channel(obj2, ch=ch, exclude="")
    ch = get_channel(obj1, ch=ch, exclude="")
    clabels = labels(obj1)

    # get time vector
    if seg[2] <= epoch_len(obj1)
        s1 = obj1.data[:, seg[1]:seg[2], 1]
        s2 = obj2.data[:, seg[1]:seg[2], 1]
    else
        s1 = epoch(obj1, ep_n=1).data[:, seg[1]:seg[2], 1]
        s2 = epoch(obj2, ep_n=1).data[:, seg[1]:seg[2], 1]
    end
    t = obj1.time_pts[seg[1]:seg[2]]

    _, t_s1, _, t_s2 = _convert_t(t[1], t[end])
    ep = _s2epoch(obj1, seg[1], seg[2])

    (ch isa(Vector{Int64}) && length(ch) == 1) && (ch = ch[1])

    xl, yl, tt = "", "", ""

    # sort channels by their type
    ctypes = obj1.header.recording[:channel_type]
    if !isa(ch, Int64)
        s1 = @views s1[ch, :]
        s2 = @views s2[ch, :]
        ctypes = ctypes[ch]
        clabels = clabels[ch]
        cunits = obj1.header.recording[:unit][ch]
    end

    if isa(ch, Int64)
        xl, yl, tt = _set_defaults(xlabel,
                                   ylabel,
                                   title,
                                   "Time [s]",
                                   "",
                                   "")
        ylabel == "default" && (yl = "Amplitude [$(_ch_units(obj1, labels(obj1)[ch]))]")
        p = plot_signal(t,
                        vec(s1[ch, :]),
                        vec(s2[ch, :]),
                        xlabel=xl,
                        ylabel=yl,
                        title=tt,
                        kwargs...)
    else
        xl, yl, tt = _set_defaults(xlabel,
                                   ylabel,
                                   title,
                                   "Time [s]",
                                   "",
                                   "")
        p = plot_signal(t,
                        s1,
                        s2,
                        ctypes=ctypes,
                        clabels=clabels,
                        cunits=cunits,
                        xlabel=xl,
                        ylabel=yl,
                        title=tt,
                        scale=scale;
                        kwargs...)
    end

    Plots.plot(p)

    return p

end
