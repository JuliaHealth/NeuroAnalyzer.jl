export plot_erp
export plot_erp_topo
export plot_erp_stack
export plot_gfp

"""
    plot_erp(t, s; <keyword arguments>)

Plot ERP/ERF (single channel).

# Arguments

  - `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
  - `s::AbstractVector`: data to plot
  - `rt::Union{Nothing, Real}=nothing`: response time value(s)
  - `xlabel::String=""`: x-axis label
  - `ylabel::String=""`: y-axis label
  - `title::String=""`: plot title
  - `yrev::Bool=false`: reverse y-axis
  - `zl::Bool=true`: draw line at t = 0
  - `mono::Bool=false`: use color or gray palette

# Returns

  - `p::GLMakie.Figure`
"""
function plot_erp(
    t::Union{AbstractVector, AbstractRange},
    s::AbstractVector;
    rt::Union{Nothing, Real} = nothing,
    xlabel::String = "",
    ylabel::String = "",
    title::String = "",
    yrev::Bool = false,
    zl::Bool = true,
    mono::Bool = false,
)::GLMakie.Figure

    # prepare plot
    GLMakie.activate!(title = "plot_erp()")
    plot_size = (900, 450)
    p = GLMakie.Figure(size = plot_size)
    ax = GLMakie.Axis(
        p[1, 1];
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        xticks = LinearTicks(10),
        xminorticksvisible = true,
        xminorticks = IntervalsBetween(10),
        yticks = LinearTicks(10),
        yminorticksvisible = true,
        yminorticks = IntervalsBetween(10),
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        yreversed = yrev,
        xzoomlock = true,
        yzoomlock = true,
        xpanlock = true,
        ypanlock = true,
        xrectzoom = false,
        yrectzoom = false,
    )
    GLMakie.ylims!(ax, yrev ? reverse(_ylims(s) .* 1.5) : (_ylims(s) .* 1.5))
    ax.titlesize = 18
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    # plot 0 v-line
    if zl
        GLMakie.vlines!(ax, 0, color = :gray, linestyle = :dash, linewidth = 2)
    end

    # plot ERP
    GLMakie.lines!(ax, t, s, color = :black, linewidth = 1)

    # plot RT v-line
    if !isnothing(rt)
        if rt >= t[1] && rt <= t[end]
            GLMakie.vlines!(ax, rt, linewidth = 1, color = mono ? :black : :red)
        end
    end

    return p

end

"""
    plot_erp(t, s; <keyword arguments>)

Plot ERP/ERF (multi-channel).

# Arguments

  - `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
  - `s::AbstractMatrix`: data to plot
  - `rt::Union{Nothing, Real}=nothing`: response time value(s)
  - `clabels::Vector{String}=string.(1:size(s, 1))`: signal channel labels vector
  - `xlabel::String=""`: x-axis label
  - `ylabel::String=""`: y-axis label
  - `title::String=""`: plot title
  - `yrev::Bool=false`: reverse y-axis
  - `avg::Bool=true`: if true, plot averaged ERP
  - `ci95::Bool=false`: if true, plot mean and ±95% CI
  - `leg::Bool=true`: if true, add legend with channel labels
  - `zl::Bool=true`: draw line at t = 0
  - `mono::Bool=false`: use color or gray palette

# Returns

  - `p::GLMakie.Figure`
"""
function plot_erp(
    t::Union{AbstractVector, AbstractRange},
    s::AbstractMatrix;
    rt::Union{Nothing, Real} = nothing,
    clabels::Vector{String} = string.(1:size(s, 1)),
    xlabel::String = "",
    ylabel::String = "",
    title::String = "",
    yrev::Bool = false,
    avg::Bool = true,
    ci95::Bool = false,
    leg::Bool = true,
    zl::Bool = true,
    mono::Bool = false,
)::GLMakie.Figure

    pal = mono ? :grays : :darktest

    ch_n = size(s, 1)

    # prepare plot
    GLMakie.activate!(title = "plot_erp()")
    plot_size = (900, 450)
    p = GLMakie.Figure(size = plot_size)
    ax = GLMakie.Axis(
        p[1, 1];
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        xticks = LinearTicks(10),
        xminorticksvisible = true,
        xminorticks = IntervalsBetween(10),
        yticks = LinearTicks(10),
        yminorticksvisible = true,
        yminorticks = IntervalsBetween(10),
        yreversed = yrev,
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        xzoomlock = true,
        yzoomlock = true,
        xpanlock = true,
        ypanlock = true,
        xrectzoom = false,
        yrectzoom = false,
    )
    GLMakie.ylims!(ax, yrev ? reverse(_ylims(s) .* 1.5) : (_ylims(s) .* 1.5))
    ax.titlesize = 18
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    # plot 0 v-line
    if zl
        GLMakie.vlines!(ax, 0, color = :gray, linestyle = :dash, linewidth = 2)
    end

    # plot ERPs
    if ci95
        avg = false
        leg = false
        s_m, _, s_u, s_l = NeuroAnalyzer.msci95(s)
        # draw 95% CI
        Makie.band!(ax, t, s_u, s_l, alpha = 0.25, color = :grey, strokewidth = 0.5)

        # draw mean
        Makie.lines!(ax, t, s_m, color = :black, linewidth = 2)
    else
        cmap = GLMakie.resample_cmap(pal, ch_n)
        for idx in 1:ch_n
            GLMakie.lines!(
                ax,
                t,
                s[idx, :];
                color = cmap[idx],
                colormap = pal,
                colorrange = 1:ch_n,
                linewidth = 1,
                alpha = avg ? 0.25 : 1.0,
                label = clabels[idx],
            )
        end
    end

    # plot averaged ERP
    if avg
        GLMakie.lines!(ax, t, mean(s, dims = 1)[:], color = :black, linewidth = 2)
    end

    # plot RT v-line
    if !isnothing(rt)
        if rt >= t[1] && rt <= t[end]
            GLMakie.vlines!(ax, rt, linewidth = 1.0, color = mono ? :black : :red)
        end
    end

    (leg && ch_n < 30) && axislegend(position = :rt, colormap = pal)

    return p

end

"""
    plot_erp_topo(locs, t, s; <keyword arguments>)

Plot topographical map ERPs.

# Arguments

  - `locs::DataFrame`: columns: channel, labels, loc_radius, loc_theta, loc_x, loc_y, loc_z, loc_radius_sph, loc_theta_sph, loc_phi_sph
  - `t::Vector{Float64}`: time vector
  - `s::Matrix{Float64}`: ERPs
  - `rt::Union{Nothing, Real}=nothing`: response time value(s)
  - `clabels::Vector{String}=string.(1:size(s, 1))`: signal channel labels vector
  - `xlabel::String=""`: x-axis label
  - `ylabel::String=""`: y-axis label
  - `title::String=""`: plot title
  - `yrev::Bool=false`: reverse y-axis
  - `cart::Bool=false`: if true, use Cartesian coordinates, otherwise use polar coordinates for XY plane and spherical coordinates for XZ and YZ planes
  - `head::Bool=true`: plot head shape
  - `zl::Bool=true`: draw line at t = 0
  - `mono::Bool=false`: use color or gray palette

# Returns

  - `p::GLMakie.Figure`
"""
function plot_erp_topo(
    locs::DataFrame,
    t::Vector{Float64},
    s::Matrix{Float64};
    rt::Union{Nothing, Real} = nothing,
    clabels::Vector{String} = string.(1:size(s, 1)),
    title::String = "",
    xlabel::String = "",
    ylabel::String = "",
    yrev::Bool = false,
    cart::Bool = false,
    head::Bool = true,
    zl::Bool = true,
    mono::Bool = false,
)::GLMakie.Figure

    @assert size(s, 2) == length(t) "Signal length must equal time length."

    pal = mono ? :grays : :darktest

    pos = collect(1:DataFrames.nrow(locs))

    # plot parameters
    if size(s, 1) <= 64
        plot_size = (1000, 1000)
        marker_size = (150, 75)
        xl = 1.2
        yl = 1.2
    elseif _in(size(s, 1), (64, 100))
        plot_size = (1200, 1200)
        marker_size = (110, 55)
        xl = 1.5
        yl = 1.5
    else
        plot_size = (1400, 1400)
        marker_size = (90, 45)
        xl = 1.5
        yl = 1.5
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

    # prepare ERP/ERF plots
    pp_vec = GLMakie.Figure[]
    pp_full_vec = GLMakie.Figure[]
    for idx in axes(s, 1)
        pp = GLMakie.Figure(size = marker_size, figure_padding = 0)
        ax = GLMakie.Axis(
            pp[1, 1];
            xlabel = "",
            ylabel = "",
            title = locs[idx, :label],
            yreversed = yrev,
            xautolimitmargin = (0, 0),
            yautolimitmargin = (0.1, 0.1),
        )
        hidedecorations!(ax)
        ax.titlesize = 8
        # plot ERPs
        GLMakie.hlines!(ax, 0, color = :black, linewidth = 1)
        GLMakie.vlines!(ax, 0, color = :gray, linestyle = :dash, linewidth = 1)
        GLMakie.lines!(ax, t, s[idx, :], linewidth = 1, color = :black)
        # plot RT v-line
        if !isnothing(rt)
            if rt >= t[1] && rt <= t[end]
                GLMakie.vlines!(ax, rt, linewidth = 1, color = mono ? :black : :red)
            end
        end
        push!(pp_vec, pp)
        pp_full = plot_erp(t, s[idx, :], xlabel = xlabel, ylabel = ylabel, title = title, rt = rt, yrev = yrev)
        push!(pp_full_vec, pp_full)
    end

    # prepare plot
    GLMakie.activate!(title = "plot_erp()")
    p = GLMakie.Figure(
                    size = plot_size,
                    figure_padding = 0,
                )
    ax = GLMakie.Axis(
        p[1, 1];
        xlabel = "",
        ylabel = "",
        title = title,
        aspect = 1,
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        xzoomlock = true,
        yzoomlock = true,
        xpanlock = true,
        ypanlock = true,
        xrectzoom = false,
        yrectzoom = false,
    )
    GLMakie.xlims!(ax, (-xl, xl))
    GLMakie.ylims!(ax, (-yl, yl))
    hidespines!(ax)
    hidedecorations!(ax)
    ax.titlesize = 18

    if head
        # nose
        GLMakie.lines!(ax, [-0.2, 0], [0.980, 1.08], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [0.2, 0], [0.980, 1.08], linewidth = 3, color = :black)

        # ears
        # left
        GLMakie.lines!(ax, [-0.995, -1.03], [0.1, 0.15], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [-1.03, -1.06], [0.15, 0.16], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [-1.06, -1.1], [0.16, 0.14], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [-1.1, -1.12], [0.14, 0.05], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [-1.12, -1.10], [0.05, -0.1], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [-1.10, -1.13], [-0.1, -0.3], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [-1.13, -1.09], [-0.3, -0.37], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [-1.09, -1.02], [-0.37, -0.39], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [-1.02, -0.98], [-0.39, -0.33], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [-0.98, -0.975], [-0.33, -0.22], linewidth = 3, color = :black)
        # right
        GLMakie.lines!(ax, [0.995, 1.03], [0.1, 0.15], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [1.03, 1.06], [0.15, 0.16], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [1.06, 1.1], [0.16, 0.14], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [1.1, 1.12], [0.14, 0.05], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [1.12, 1.10], [0.05, -0.1], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [1.10, 1.13], [-0.1, -0.3], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [1.13, 1.09], [-0.3, -0.37], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [1.09, 1.02], [-0.37, -0.39], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [1.02, 0.98], [-0.39, -0.33], linewidth = 3, color = :black)
        GLMakie.lines!(ax, [0.98, 0.975], [-0.33, -0.22], linewidth = 3, color = :black)

        # head
        GLMakie.arc!(ax, (0, 0), 1, 0, 2pi, linewidth = 3, color = :black)
    end

    for idx in axes(s, 1)
        io = IOBuffer()
        show(io, MIME"image/png"(), pp_vec[idx])
        pp = FileIO.load(io)
        GLMakie.scatter!(loc_x[idx], loc_y[idx], marker = pp, markersize = marker_size, markerspace = :pixel)
    end

    loc_x_range = Tuple{Float64, Float64}[]
    loc_y_range = Tuple{Float64, Float64}[]
    for idx in eachindex(loc_x)
        push!(loc_x_range, (loc_x[idx] - 0.15, loc_x[idx] + 0.15))
        push!(loc_y_range, (loc_y[idx] - 0.1, loc_y[idx] + 0.1))
    end
    on(events(p).mousebutton) do event
        if event.button == Mouse.left
            if event.action == Mouse.press
                ax_x = mouseposition(ax)[1]
                ax_y = mouseposition(ax)[2]
                for idx in eachindex(loc_x)
                    if ax_x >= loc_x_range[idx][1] &&
                        ax_x <= loc_x_range[idx][2] &&
                        ax_y >= loc_y_range[idx][1] &&
                        ax_y <= loc_y_range[idx][2]
                        display(GLMakie.Screen(), pp_full_vec[idx])
                        break
                    end
                end
            end
        end
    end

    return p

end

"""
    plot_erp_stack(t, s; <keyword arguments>)

Plot EPRs stacked by channels or by epochs.

# Arguments

  - `t::AbstractVector`: x-axis values
  - `s::AbstractMatrix`
  - `rt::Union{Nothing, AbstractVector}=nothing`: response time for each epoch; if provided, the response time line will be plotted over the `:stack` plot
  - `clabels::Vector{String}=string.(1:size(s, 1))`: signal channel labels vector
  - `xlabel::String=""`: x-axis label
  - `ylabel::String=""`: y-axis label
  - `title::String=""`: plot title
  - `cb::Bool=true`: plot color bar
  - `cb_title::String=""`: color bar title
  - `smooth::Bool=false`: smooth the image using Gaussian blur
  - `ks::Int64=3`: kernel size of the Gaussian blur (larger kernel means more smoothing)
  - `zl::Bool=true`: draw line at t = 0
  - `mono::Bool=false`: use color or gray palette

# Returns

  - `p::GLMakie.Figure`
"""
function plot_erp_stack(
    t::AbstractVector,
    s::AbstractMatrix;
    rt::Union{Nothing, AbstractVector} = nothing,
    clabels::Vector{String} = string.(1:size(s, 1)),
    xlabel::String = "",
    ylabel::String = "",
    title::String = "",
    cb::Bool = true,
    cb_title::String = "",
    smooth::Bool = false,
    ks::Int64 = 3,
    zl::Bool = true,
    mono::Bool = false,
)::GLMakie.Figure

    @assert length(t) == size(s, 2) "Number of s columns ($(size(s, 2))) must equal length of t ($(length(t)))."

    if !isnothing(rt)
        @assert length(rt) == size(s, 1) "Length of the rt vector must equal number of ERP epochs ($(size(s, 1)))."
    end

    pal = mono ? :grays : :darktest

    if smooth
        s = imfilter(s, Kernel.gaussian(ks))
    end

    # prepare plot
    GLMakie.activate!(title = "plot_erp()")
    plot_size = size(s, 1) <= 64 ? (900, 600) : (900, 900)
    p = GLMakie.Figure(size = plot_size)
    ax = GLMakie.Axis(
        p[1, 1];
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        xticks = LinearTicks(10),
        xminorticksvisible = true,
        xminorticks = IntervalsBetween(10),
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        yticks = size(s, 1) <= 30 ? (1:length(clabels), clabels) : (5:5:length(clabels), clabels[5:5:end]),
        yticklabelsize = size(s, 1) <= 64 ? 8 : 5,
        xzoomlock = true,
        yzoomlock = true,
        xpanlock = true,
        ypanlock = true,
        xrectzoom = false,
        yrectzoom = false,
    )
    ax.titlesize = 18
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    hm = GLMakie.heatmap!(ax, t, axes(s, 1), rotr90(s), colormap = pal)

    # plot 0 v-line
    if zl
        GLMakie.vlines!(ax, 0, color = :white, linestyle = :dash, linewidth = 2)
    end

    # plot RT v-line
    if !isnothing(rt)
        for idx in eachindex(rt)
            if rt[idx] >= t[1] && rt[idx] <= t[end]
                GLMakie.lines!(ax, rt[idx], idx, linewidth = 1, color = mono ? :black : :red)
            end
        end
    end

    # draw colorbar
    if cb
        Colorbar(p[1, 2], hm, label = cb_title, labelsize = 16)
    end

    return p

end

"""
    plot_gfp(t, s; <keyword arguments>)

Plot Global Field Power.

# Arguments

  - `t::Union{AbstractVector, AbstractRange}`: x-axis values (usually time)
  - `g::AbstractVector`: data to plot
  - `rt::Union{Nothing, Real}=nothing`: response time value(s)
  - `xlabel::String=""`: x-axis label
  - `ylabel::String=""`: y-axis label
  - `title::String=""`: plot title
  - `zl::Bool=true`: draw line at t = 0
  - `mono::Bool=false`: use color or gray palette

# Returns

  - `p::GLMakie.Figure`
"""
function plot_gfp(
    t::Union{AbstractVector, AbstractRange},
    g::AbstractVector;
    rt::Union{Nothing, Real} = nothing,
    xlabel::String = "",
    ylabel::String = "",
    title::String = "",
    yrev::Bool = false,
    zl::Bool = true,
    mono::Bool = false,
)::GLMakie.Figure

    # prepare plot
    GLMakie.activate!(title = "plot_gfp()")
    plot_size = (900, 450)
    p = GLMakie.Figure(size = plot_size)
    ax = GLMakie.Axis(
        p[1, 1];
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        xticks = LinearTicks(10),
        xminorticksvisible = true,
        xminorticks = IntervalsBetween(10),
        yticks = LinearTicks(10),
        yminorticksvisible = true,
        yminorticks = IntervalsBetween(2),
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        xzoomlock = true,
        yzoomlock = true,
        xpanlock = true,
        ypanlock = true,
        xrectzoom = false,
        yrectzoom = false,
    )
    GLMakie.ylims!(ax, 0, maximum(g) * 1.5)
    ax.titlesize = 18
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    # plot 0 v-line
    if zl
        GLMakie.vlines!(ax, 0, color = :gray, linestyle = :dash, linewidth = 2)
    end

    # plot GFP
    GLMakie.band!(ax, t, 0, g, color = :gray)
    GLMakie.lines!(ax, t, g, color = :black, linewidth = 1)

    # plot RT v-line
    if !isnothing(rt)
        if rt >= t[1] && rt <= t[end]
            GLMakie.vlines!(ax, rt, linewidth = 1, color = mono ? :black : :red)
        end
    end

    return p

end


"""
    plot_erp(obj; <keyword arguments>)

Plot ERP/ERF.

# Arguments

  - `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  - `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
  - `tm::Union{Nothing, Int64, Vector{Int64}}=nothing`: time markers (in milliseconds) to plot as vertical lines, useful for adding topoplots at these time points
  - `xlabel::String="default"`: x-axis label
  - `ylabel::String="default"`: y-axis label
  - `title::String="default"`: plot title
  - `cb::Bool=true`: plot color bar
  - `cb_title::String="default"`: color bar title
  - `peaks::Bool=true`: draw peaks
  - `leg::Bool=true`: if true, add legend with channel labels
  - `type::Symbol=:normal`: plot type:
      + `:normal`
      + `:gfp`: plot Global Field Power
      + `:stack`: stacked epochs/channels
      + `:topo`: topographical plot of ERPs
  - `yrev::Bool=false`: reverse y-axis
  - `avg::Bool=true`: if true, plot averaged ERP
  - `ci95::Bool=false`: if true, plot mean and ±95% CI
  - `smooth::Bool=false`: smooth the image using Gaussian blur
  - `ks::Int64=3`: kernel size of the Gaussian blur (larger kernel means more smoothing)
  - `rt::Union{Nothing, Real, AbstractVector}=nothing`: response time for each epoch; if provided, the response time line will be plotted over the `:stack` plot
  - `sort_epochs::Bool=false`:: sort epochs by rt vector
  - `zl::Bool=true`: draw line at t = 0
  - `mono::Bool=false`: use color or gray palette
  - `gui::Bool=false`: ignored

# Returns

  - `p::GLMakie.Figure`
"""
function plot_erp(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    tm::Union{Nothing, Int64, Vector{Int64}} = nothing,
    xlabel::String = "default",
    ylabel::String = "default",
    title::String = "default",
    cb::Bool = true,
    cb_title::String = "default",
    peaks::Bool = true,
    leg::Bool = true,
    type::Symbol = :normal,
    yrev::Bool = false,
    avg::Bool = true,
    ci95::Bool = false,
    smooth::Bool = false,
    ks::Int64 = 3,
    rt::Union{Nothing, Real, AbstractVector} = nothing,
    sort_epochs::Bool = false,
    zl::Bool = true,
    mono::Bool = false,
    gui::Bool = false,
)::GLMakie.Figure

    _check_datatype(obj, ["erp", "erf"])
    _check_var(type, [:normal, :topo, :stack, :gfp], "type")

    # check channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")
    @assert !(length(ch) > 1 && length(unique(obj.header.recording[:channel_type][ch])) > 1) "All channels must be of the same type."
    length(ch) > 1 && (eavg = false)
    type === :gfp && @assert length(ch) > 1 "More than 1 channel must be selected."

    # set units
    units = _ch_units(obj, labels(obj)[ch[1]])

    # get data
    ep_n = nepochs(obj) - 1
    if length(ch) == 1
        if type === :stack
            s = obj.data[ch[1], :, 2:end]'
        else
            s = obj.data[ch, :, 1][:]
        end
    else
        s = obj.data[ch, :, 1]
    end

    # get labels
    clabels = labels(obj)[ch]

    # get time vector
    t = obj.epoch_time

    if length(ch) == 1

        if type === :stack

            xl, yl, tt = _set_defaults(
                xlabel, ylabel, title, "Time [ms]", "Epochs", "ERP amplitude, $(clabels[1]), avgₑ: $ep_n"
            )
            cb_title == "default" && (cb_title = "Amplitude [$units]")

            if sort_epochs
                if !isnothing(rt)
                    rt_idx = sortperm(rt)
                    rt = rt[rt_idx]
                    s = s[rt_idx, :]
                end
            end

            p = plot_erp_stack(
                t,
                s;
                rt = rt,
                xlabel = xl,
                ylabel = yl,
                title = tt,
                cb_title = cb_title,
                smooth = smooth,
                ks = ks,
                zl = zl,
                mono = mono,
            )

        elseif type === :normal

            xl, yl, tt = _set_defaults(
                xlabel, ylabel, title, "Time [ms]", "Amplitude [$units]", "ERP amplitude, $(clabels[1]), avgₑ: $ep_n"
            )
            p = plot_erp(t, s, xlabel = xl, ylabel = yl, title = tt, rt = rt, yrev = yrev, zl = zl, mono = mono)

        end

    elseif type === :normal

        avg == false && (peaks = false)
        xl, yl, tt = _set_defaults(
            xlabel,
            ylabel,
            title,
            "Time [ms]",
            "Amplitude [$units]",
            "ERP amplitude, $(length(ch)) channels, avgₑ: $ep_n",
        )
        p = plot_erp(
            t,
            s;
            xlabel = xl,
            ylabel = yl,
            title = tt,
            clabels = clabels,
            rt = rt,
            yrev = yrev,
            avg = avg,
            ci95 = ci95,
            leg = leg,
            zl = zl,
            mono = mono,
        )

    elseif type === :stack

        cb_title == "default" && (cb_title = "Amplitude [$units]")
        xl, yl, tt = _set_defaults(
            xlabel, ylabel, title, "Time [ms]", "", "ERP amplitude, $(length(ch)) channels, avgₑ: $ep_n"
        )
        p = plot_erp_stack(
            t,
            s;
            rt = rt,
            xlabel = xl,
            ylabel = yl,
            title = tt,
            clabels = clabels,
            cb = cb,
            cb_title = cb_title,
            smooth = smooth,
            ks = ks,
            zl = zl,
            mono = mono,
        )

    elseif type === :gfp

        g = erp_gfp(obj, ch = labels(obj)[ch])
        xl, yl, tt = _set_defaults(
            xlabel,
            ylabel,
            title,
            "Time [ms]",
            "GFP [$units]",
            "Global Field Power, $(length(ch)) channels, avgₑ: $ep_n",
        )
        p = plot_gfp(t, g, xlabel = xl, ylabel = yl, title = tt, rt = rt, zl = zl, mono = mono)

    elseif type === :topo

        _has_locs(obj)
        xl, yl, tt = _set_defaults(
            xlabel, ylabel, title, "Time [ms]", "Amplitude [$units]", "ERP amplitude, avgₑ: $ep_n"
        )
        _check_ch_locs(ch, labels(obj), obj.locs[!, :label])
        @assert length(unique(obj.header.recording[:channel_type][ch])) == 1 "For multi-channel topo plot all channels must be of the same type."
        _has_locs(obj)
        chs = intersect(obj.locs[!, :label], labels(obj)[ch])
        locs = Base.filter(:label => in(chs), obj.locs)
        _check_ch_locs(ch, labels(obj), obj.locs[!, :label])
        p = plot_erp_topo(
            locs,
            t,
            s;
            xlabel = xl,
            ylabel = yl,
            title = tt,
            clabels = clabels,
            rt = rt,
            yrev = yrev,
            mono = mono,
            zl = zl,
        )
    end

    # draw time markers
    if !isnothing(tm)
        for idx in eachindex(tm)
            @assert tm[idx] / 1000 >= t[1] "tm value ($(tm[idx])) is out of epoch time segment ($(t[1]):$(t[end]))."
            @assert tm[idx] / 1000 <= t[end] "tm value ($(tm[idx])) is out of epoch time segment ($(t[1]):$(t[end]))."
            tm[idx] = vsearch(tm[idx] / 1000, t)
            GLMakie.vlines!(p[1, 1], t[tm[idx]], linewidth = 0.5, color = :black)
        end
    end

    # draw peaks
    if peaks
        if length(ch) == 1 && type === :normal
            pp = erp_peaks(obj)
            GLMakie.scatter!(
                p[1, 1],
                t[pp[ch, 1]][1],
                obj.data[ch, pp[ch, 1], 1][1];
                marker = :xcross,
                color = mono ? :black : :red,
                markersize = 15,
            )
            GLMakie.scatter!(
                p[1, 1],
                t[pp[ch, 2]][1],
                obj.data[ch, pp[ch, 2], 1][1];
                marker = :xcross,
                color = mono ? :black : :blue,
                markersize = 15,
            )
            _info("Positive peak time: $(round(t[pp[ch, 1]][1] * 1000, digits=0)) ms")
            _info("Positive peak amplitude: $(round(obj.data[ch, pp[ch, 1], 1][1], digits=2)) $units")
            _info("Negative peak time: $(round(t[pp[ch, 2]][1] * 1000, digits=0)) ms")
            _info("Negative peak amplitude: $(round(obj.data[ch, pp[ch, 2], 1][1], digits=2)) $units")
        elseif length(ch) > 1 && type === :normal
            mep_tmp = mean(obj.data[ch, :, 1], dims = 1)[:, :, :]
            obj_tmp = keep_channel(obj, ch = labels(obj)[1])
            obj_tmp.data = mep_tmp
            pp = erp_peaks(obj_tmp)
            GLMakie.scatter!(
                p[1, 1], t[pp[1, 1]], mep_tmp[pp[1, 1]]; marker = :xcross, color = mono ? :black : :red, markersize = 15
            )
            GLMakie.scatter!(
                p[1, 1],
                t[pp[1, 2]],
                mep_tmp[pp[1, 2]];
                marker = :xcross,
                color = mono ? :black : :blue,
                markersize = 15,
            )
            _info("Positive peak time: $(round(t[pp[1, 1]] * 1000, digits=0)) ms")
            _info("Positive peak amplitude: $(round(mep_tmp[pp[1, 1]], digits=2)) $units")
            _info("Negative peak time: $(round(t[pp[1, 2]] * 1000, digits=0)) ms")
            _info("Negative peak amplitude: $(round(mep_tmp[pp[1, 2]], digits=2)) $units")
        end
    end

    return p

end
