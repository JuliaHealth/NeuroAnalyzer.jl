export plot_erp
export plot_erp_butterfly
export plot_erp_avg
export plot_erp_topo
export plot_erp_stack

"""
    plot_erp(t, s, bad; <keyword arguments>)

Plot ERP/ERF.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
- `s::AbstractVector`: data to plot
- `rt::Union{Nothing, Real}=nothing`:: response time value
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `yrev::Bool=false`: reverse Y axis
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_erp(t::Union{AbstractVector, AbstractRange}, s::AbstractVector; rt::Union{Nothing, Real}=nothing, xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, yrev::Bool=false, kwargs...)::Plots.Plot{Plots.GRBackend}

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

    # plot RT v-line
    if !isnothing(rt)
        p = Plots.vline!([rt],
                         linewidth=1.0,
                         linecolor=:red,
                         label=false)
    end

    return p

end

"""
    plot_erp_butterfly(t, s; <keyword arguments>)

Butterfly plot of ERP.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
- `s::AbstractArray`: data to plot
- `rt::Union{Nothing, Real}=nothing`:: response time value
- `clabels::Vector{String}=[""]`: signal channel labels vector
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `avg::Bool=false`: plot average ERP
- `yrev::Bool=false`: reverse Y axis
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_erp_butterfly(t::Union{AbstractVector, AbstractRange}, s::AbstractArray; rt::Union{Nothing, Real}=nothing, clabels::Vector{String}=[""], xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, avg::Bool=true, yrev::Bool=false, kwargs...)::Plots.Plot{Plots.GRBackend}

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
                   legend=ch_n < 20,
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
                            alpha=0.2)
        else
            if clabels == repeat([""], ch_n)
                p = Plots.plot!(t,
                                s[idx, :],
                                t=:line,
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
    if avg
        if ch_n == 1
            s = mean(s, dims=2)[:]
        else
            s = mean(s, dims=1)[:]
        end
        p = Plots.plot!(t,
                        s,
                        linewidth=2,
                        linecolor=:black,
                        label=false)
    end

    # plot 0 v-line
    p = Plots.vline!([0],
                     linestyle=:dash,
                     linewidth=0.5,
                     linecolor=:black,
                     label=false)

    # plot RT v-line
    if !isnothing(rt)
        p = Plots.vline!([rt],
                         linewidth=1.0,
                         linecolor=:red,
                         label=false)
    end

    return p

end

"""
    plot_erp_avg(t, s; <keyword arguments>)

Plot ERP/ERF amplitude mean and ±95% CI.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
- `s::AbstractArray`: data to plot
- `rt::Union{Nothing, Real}=nothing`:: response time value
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `yrev::Bool=false`: reverse Y axis
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_erp_avg(t::Union{AbstractVector, AbstractRange}, s::AbstractArray; rt::Union{Nothing, Real}=nothing, xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, yrev::Bool=false, kwargs...)::Plots.Plot{Plots.GRBackend}

    pal = mono ? :grays : :darktest

    # get mean and 95%CI
    s_m, _, s_u, s_l = NeuroAnalyzer.msci95(s)

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

    # plot RT v-line
    if !isnothing(rt)
        p = Plots.vline!([rt],
                         linewidth=1.0,
                         linecolor=:red,
                         label=false)
    end

    return p

end

"""
    plot_erp_topo(locs, t, s; <keyword arguments>)

Plot topographical map ERPs.

# Arguments

- `locs::DataFrame`: columns: channel, labels, loc_radius, loc_theta, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
- `t::Vector{Float64}`: time vector
- `s::Matrix{Float64}`: ERPs
- `ch::Union{Vector{Int64}, AbstractRange}`: which channels to plot
- `clabels::Vector{String}=[""]`: signal channel labels vector
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `yrev::Bool=false`: reverse Y axis
- `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates for XY plane and spherical coordinates for XZ and YZ planes
- `mono::Bool=false`: use color or gray palette
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_erp_topo(locs::DataFrame, t::Vector{Float64}, s::Matrix{Float64}; ch=Union{Vector{Int64}, AbstractRange}, clabels::Vector{String}=[""], xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, yrev::Bool=false, cart::Bool=false, kwargs...)::Plots.Plot{Plots.GRBackend}

    @assert size(s, 2) == length(t) "Signal length and time length must be equal."

    pal = mono ? :grays : :darktest

    # channel labels
    clabels == [""] && (clabels = repeat([""], size(s, 1)))

    # get limits
    ylim = (floor(minimum(s) * 1.1, digits=0), ceil(maximum(s) * 1.1, digits=0))
    ylim = _tuple_max(ylim)

    # plot parameters
    if size(s, 1) <= 64
        plot_size = 800
        marker_size = (120, 80)
        rx = 1
        ry = 1
        offset = 0
    elseif size(s, 1) <= 100
        plot_size = 1200
        marker_size = (120, 80)
        rx = 0.9
        ry = 0.7
        offset = 150
    else
        plot_size = 2500
        marker_size = (90, 50)
        rx = 0.6
        ry = 0.4
        offset = 250
    end

    # get locations
    if !cart
        loc_x = zeros(size(locs, 1))
        loc_y = zeros(size(locs, 1))
        for idx in axes(locs, 1)
            loc_x[idx], loc_y[idx] = pol2cart(locs[!, :loc_radius][idx], locs[!, :loc_theta][idx])
        end
    else
        loc_x = locs[!, :loc_x]
        loc_y = locs[!, :loc_y]
    end
    loc_x = _s2v(loc_x)
    loc_y = _s2v(loc_y)
    # get marker centers
    loc_x .*= ((plot_size / 2) - marker_size[1] / 2) .* rx
    loc_y .*= ((plot_size / 2) - marker_size[2] / 2) .* ry
    # origin is in the left top corner, convert positions
    loc_x = round.(Int64, loc_x .+ (plot_size / 2) .- marker_size[1] / 2)
    loc_y = (plot_size - marker_size[2]) .- round.(Int64, loc_y .+ (plot_size / 2) .- marker_size[2] / 2)

    c = CairoRGBSurface(plot_size, plot_size - 3 * offset)
    cr = CairoContext(c)
    Cairo.set_source_rgb(cr, 256, 256, 256)
    Cairo.rectangle(cr, 0.0, 0.0, plot_size, plot_size - 3 * offset)
    Cairo.fill(cr)
    for idx in axes(s, 1)
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
        yrev && yflip!(true)

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

        show(io, MIME("image/png"), p)
        img = read_from_png(io)
        Cairo.set_source_surface(cr, img, loc_x[idx], loc_y[idx] - 1.5 * offset)
        Cairo.paint(cr)
    end

    img_png = tempname() * ".png"
    Cairo.write_to_png(c, img_png)
    img = FileIO.load(img_png)
    p = nothing
    p = Plots.plot(img,
                   size=(plot_size, plot_size - 3 * offset),
                   title=title,
                   titlefontsize=12,
                   border=:none)
    rm(img_png)

    return p

end

"""
    plot_erp_stack(s; <keyword arguments>)

Plot EPRs stacked by channels or by epochs.

# Arguments

- `t::AbstractVector`: x-axis values
- `s::AbstractArray`
- `rt::Union{Nothing, AbstractVector}=nothing`: response time for each epoch; if provided, the response time line will be plotted over the `:stack` plot
- `clabels::Vector{String}=[""]`: signal channel labels vector
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `cb::Bool=true`: plot color bar
- `cb_title::String=""`: color bar title
- `mono::Bool=false`: use color or gray palette
- `smooth::Bool=false`: smooth the image using Gaussian blur
- `n::Int64=3`: kernel size of the Gaussian blur (larger kernel means more smoothing)
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_erp_stack(t::AbstractVector, s::AbstractArray, rt::Union{Nothing, AbstractVector}=nothing; clabels::Vector{String}=[""], xlabel::String="", ylabel::String="", title::String="", cb::Bool=true, cb_title::String="", mono::Bool=false, smooth::Bool=false, n::Int64=3, kwargs...)::Plots.Plot{Plots.GRBackend}

    _chk2d(s)
    @assert length(t) == size(s, 2) "Number of s columns ($(size(s, 2))) must be equal to length of t ($(length(t)))."

    if !isnothing(rt)
        @assert length(rt) == size(s, 1) "Length of the rt vector must be the same as the number of ERP epochs ($(size(s, 1)))."
    end

    pal = mono ? :grays : :darktest

    if clabels == [""]
        yticks = round.(Int64, range(1, size(s, 1), length=10))
    else
        yticks = (axes(s, 1), clabels)
    end

    if smooth
        s = imfilter(s, Kernel.gaussian(n))
    end

    p = Plots.heatmap(t,
                      axes(s, 1),
                      s,
                      size=size(s, 1) <= 64 ? (1200, 800) : (1200, 1200),
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
                      ytickfontsize=size(s, 1) <= 64 ? 8 : 5;
                      kwargs...)

    # plot 0 v-line
    p = Plots.vline!([0],
                     linestyle=:dash,
                     linewidth=0.5,
                     linecolor=:black,
                     label=false)

    # plot RT line
    if !isnothing(rt)
        p = Plots.plot!(rt, 1:length(rt),
                        linewidth=1,
                        linecolor=:black,
                        label=false)
    end

    Plots.plot(p)

    return p

end

"""
    plot_erp(obj; <keyword arguments>)

Plot ERP/ERF.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `tm::Union{Int64, Vector{Int64}}=0`: time markers (in miliseconds) to plot as vertical lines, useful for adding topoplots at these time points
- `xlabel::String="default"`: x-axis label, default is Time [ms]
- `ylabel::String="default"`: y-axis label, default is Amplitude [units]
- `title::String="default"`: plot title, default is ERP amplitude [channel: 1, epochs: 1:2, time window: -0.5 s:1.5 s]
- `cb::Bool=true`: plot color bar
- `cb_title::String="default"`: color bar title, default is Amplitude [units]
- `mono::Bool=false`: use color or gray palette
- `peaks::Bool=true`: draw peaks
- `channel_labels::Bool=true`: draw labels legend (using channel labels) for multi-channel `:butterfly` plot
- `type::Symbol=:normal`: plot type:
    - `:normal`
    - `:butterfly`: butterfly plot
    - `:topo`: topographical plot of ERPs
    - `:stack`: stacked epochs/channels
- `yrev::Bool=false`: reverse Y axis
- `avg::Bool=false`: plot average ERP for `:butterfly` plot
- `smooth::Bool=false`: smooth the image using Gaussian blur
- `n::Int64=3`: kernel size of the Gaussian blur (larger kernel means more smoothing)
- `rt::Union{Nothing, Real, AbstractVector}=nothing`: response time for each epoch; if provided, the response time line will be plotted over the `:stack` plot
- `sort_epochs::Bool=false`:: sort epochs by rt vector
- `kwargs`: optional arguments for plot() function

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_erp(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, tm::Union{Int64, Vector{Int64}}=0, xlabel::String="default", ylabel::String="default", title::String="default", cb::Bool=true, cb_title::String="default", mono::Bool=false, peaks::Bool=true, channel_labels::Bool=true, type::Symbol=:normal, yrev::Bool=false, avg::Bool=true, smooth::Bool=false, n::Int64=3, rt::Union{Nothing, Real, AbstractVector}=nothing, sort_epochs::Bool=false, kwargs...)::Plots.Plot{Plots.GRBackend}

    _check_datatype(obj, ["erp", "erf"])

    # check channels
    ch = get_channel(obj, ch=ch)
    length(ch) == 1  && (ch = ch[1])

    # set units
    units = _ch_units(obj, labels(obj)[ch[1]])

    _check_var(type, [:normal, :butterfly, :mean, :topo, :stack], "type")
    @assert !(length(ch) > 1 && length(unique(obj.header.recording[:channel_type][ch])) > 1) "All channels must be of the same type."

    # get data
    ep_n = nepochs(obj) - 1

    if type in [:normal, :topo]
        s = obj.data[ch, :, 1]
    else
        if ch isa Int64
            s = obj.data[ch, :, 2:end]'
            channel_labels = false
        else
            s = obj.data[ch, :, 1]
        end
    end

    # get time vector
    t = obj.epoch_time
    _, t_s1, _, t_s2 = _convert_t(t[1], t[end])

    if tm != 0
        for tm_idx in eachindex(tm)
            @assert tm[tm_idx] / 1000 >= t[1] "tm value ($(tm[tm_idx])) is out of epoch time segment ($(t[1]):$(t[end]))."
            @assert tm[tm_idx] / 1000 <= t[end] "tm value ($(tm[tm_idx])) is out of epoch time segment ($(t[1]):$(t[end]))."
            tm[tm_idx] = vsearch(tm[tm_idx] / 1000, t)
        end
    end

    if type === :normal
        @assert ch isa Int64 "For :normal plot type, only one channel must be specified."
        xl, yl, tt = _set_defaults(xlabel,
                                   ylabel,
                                   title,
                                   "Time [ms]",
                                   "Amplitude [$units]",
                                   "$(uppercase(datatype(obj))) amplitude\n[averaged epochs: $ep_n, time window: $t_s1:$t_s2]")
        p = plot_erp(t,
                     s,
                     xlabel=xl,
                     ylabel=yl,
                     title=tt,
                     mono=mono,
                     rt=rt,
                     yrev=yrev;
                     kwargs...)
    elseif type === :butterfly
        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [ms]", "Amplitude [$units]", "$(uppercase(datatype(obj))) amplitude\n[averaged epochs: $ep_n, time window: $t_s1:$t_s2]")
        if channel_labels
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
                               rt=rt,
                               yrev=yrev;
                               kwargs...)
    elseif type === :mean
        xl, yl, tt = _set_defaults(xlabel,
                                   ylabel,
                                   title,
                                   "Time [ms]",
                                   "Amplitude [$units]",
                                   "$(uppercase(datatype(obj))) amplitude [mean ± 95%CI]\n[averaged epochs: $ep_n, time window: $t_s1:$t_s2]")
        p = plot_erp_avg(t,
                         s,
                         xlabel=xl,
                         ylabel=yl,
                         title=tt,
                         mono=mono,
                         rt=rt,
                         yrev=yrev;
                         kwargs...)
    elseif type === :topo
        _has_locs(obj)
        xl, yl, tt = _set_defaults(xlabel,
                                   ylabel,
                                   title,
                                   "",
                                   "",
                                   "$(uppercase(datatype(obj))) amplitude\n[averaged epochs: $ep_n, time window: $t_s1:$t_s2]")
        chs = intersect(obj.locs[!, :label], labels(obj)[ch])
        locs = Base.filter(:label => in(chs), obj.locs)
        _check_ch_locs(ch, labels(obj), obj.locs[!, :label])
        peaks = false
        ndims(s) == 1 && (s = reshape(s, 1, length(s)))
        clabels = labels(obj)[ch]
        clabels isa String && (clabels = [clabels])
        p = plot_erp_topo(locs,
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
            xl, yl, tt = _set_defaults(xlabel,
                                       ylabel,
                                       title,
                                       "Time [ms]",
                                       "Epochs",
                                       "$(uppercase(datatype(obj))) amplitude\n[averaged epochs: $ep_n, time window: $t_s1:$t_s2]")
        else
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [ms]", "", "$(uppercase(datatype(obj))) amplitude\n[averaged epochs: $ep_n, time window: $t_s1:$t_s2]")
        end
        if channel_labels
            clabels = labels(obj)[ch]
        else
            clabels = repeat([""], length(ch))
        end

        if ch isa Int64
            if sort_epochs
                rt_idx = sortperm(rt)
                rt = rt[rt_idx]
                s = s[rt_idx, :]
            end
        else
            rt = nothing
            sort_epochs = false
        end

        p = plot_erp_stack(t,
                           s,
                           rt,
                           xlabel=xl,
                           ylabel=yl,
                           title=tt,
                           clabels=clabels,
                           cb=cb,
                           cb_title=cb_title,
                           smooth=smooth,
                           n=n,
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
        if ch isa Int64
            pp = erp_peaks(obj)
            if !mono
                Plots.scatter!((t[pp[ch, 1]], obj.data[ch, pp[ch, 1], 1]), marker=:xcross, markercolor=:red, markersize=3, label=false)
                Plots.scatter!((t[pp[ch, 2]], obj.data[ch, pp[ch, 2], 1]), marker=:xcross, markercolor=:blue, markersize=3, label=false)
            else
                Plots.scatter!((t[pp[ch, 1]], obj.data[ch, pp[ch, 1], 1]), marker=:xcross, markercolor=:black, markersize=3, label=false)
                Plots.scatter!((t[pp[ch, 2]], obj.data[ch, pp[ch, 2], 1]), marker=:xcross, markercolor=:black, markersize=3, label=false)
            end
            _info("Positive peak time: $(round(t[pp[ch, 1]] * 1000, digits=0)) ms")
            _info("Positive peak amplitude: $(round(obj.data[ch, pp[ch, 1], 1], digits=2)) $units")
            _info("Negative peak time: $(round(t[pp[ch, 2]] * 1000, digits=0)) ms")
            _info("Negative peak amplitude: $(round(obj.data[ch, pp[ch, 2], 1], digits=2)) $units")
        elseif (type === :butterfly && avg) || type === :mean
            erp_tmp = mean(mean(obj.data[ch, :, 2:end], dims=1), dims=3)
            obj_tmp = keep_channel(obj, ch=labels(obj)[1])
            obj_tmp.data = erp_tmp
            pp = erp_peaks(obj_tmp)
            if !mono
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

    Plots.plot(p)

    return p

end
