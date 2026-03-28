export plot_erp
export plot_erp_topo
export plot_erp_stack
export plot_gfp

"""
    plot_erp(t, s; <keyword arguments>)

Plot Event-Related Potential/Field (single channel).

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: time values (in seconds)
- `s::AbstractVector`: signal vector (ERP/ERF data)
- `rt::Union{Nothing, Real}=nothing`: response time (in milliseconds)
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `yrev::Bool=false`: if `true`, reverse the y-axis
- `zl::Bool`: if `true`, draw vertical line at t = 0
- `mono::Bool=false`: if `true`, use a monochrome palette

# Returns

- `GLMakie.Figure`: the plotted figure
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

    # validate
    length(t) == length(s) ||
        throw(ArgumentError("Time and signal vectors must have the same length."))
    all(isfinite, t) ||
        throw(ArgumentError("Time values must be finite."))
    all(isfinite, s) ||
        throw(ArgumentError("Signal values must be finite."))

    # prepare plot
    GLMakie.activate!(; title = "plot_erp()")
    plot_size = (900, 450)
    fig = GLMakie.Figure(; size = plot_size)

    # create axis with customizable properties
    ax = GLMakie.Axis(
        fig[1, 1];
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

    # draw zero line if requested
    if zl
        GLMakie.vlines!(ax, 0; color = :gray, linestyle = :dash, linewidth = 2)
    end

    # plot ERP signal
    GLMakie.lines!(ax, t, s; color = :black, linewidth = 1)

    # plot response time if provided
    if !isnothing(rt) && rt / 1000 ∈ t
        GLMakie.vlines!(
            ax,
            rt / 1000;
            linewidth = 1.5,
            color = mono ? :black : :red,
        )
    end

    return fig
end

"""
    plot_erp(t, s; <keyword arguments>)

Plot multi-channel Event-Related Potential/Field with optional averaging and confidence intervals.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: time values (in seconds).
- `s::AbstractMatrix`: signal data, shape (channels, samples)
- `rt::Union{Nothing, Real}=nothing`: response time (in milliseconds)
- `clabels::Vector{String}=string.(1:size(s, 1))`: channel labels (default: auto-generated)
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `yrev::Bool=false`: if `true`, reverse the y-axis
- `avg::Bool=true`: if `true`, plot averaged ERP
- `ci95::Bool=false`: if `true`, plot mean and ±95% CI
- `leg::Bool=true`: if `true`, add legend with channel labels
- `zl::Bool`: if `true`, draw vertical line at t = 0
- `mono::Bool=false`: if `true`, use a monochrome palette

# Returns

- `GLMakie.Figure`: the plotted figure
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

    # validate
    size(s, 2) == length(t) ||
        throw(ArgumentError("Signal matrix columns must match time vector length"))
    length(clabels) == size(s, 1) ||
        throw(ArgumentError("Number of channel labels must match number of channels"))
    all(isfinite, s) ||
        throw(ArgumentError("Signal contains non-finite values"))

    # select color palette
    pal = mono ? :grays : :darktest

    # number of channels
    ch_n = size(s, 1)

    # prepare plot
    GLMakie.activate!(; title = "plot_erp()")
    plot_size = (900, 450)
    fig = GLMakie.Figure(; size = plot_size)

    # create axis with customizable properties
    ax = GLMakie.Axis(
        fig[1, 1];
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

    # draw zero line if requested
    if zl
        GLMakie.vlines!(ax, 0; color = :gray, linestyle = :dash, linewidth = 2)
    end

    # plot ERPs
    if ci95
        avg = false
        leg = false

        # calculate mean and 95% CI
        msci95_data = NeuroAnalyzer.msci95(s)
        s_m = msci95_data.sm
        s_u = msci95_data.ul
        s_l = msci95_data.ll

        # plot confidence interval
        GLMakie.Makie.band!(ax, t, s_l, s_u; alpha = 0.25, color = :grey)

        # plot mean signal
        GLMakie.Makie.lines!(ax, t, s_m; color = :black, linewidth = 2)
    else
        # plot individual channels
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

    # plot averaged signal if requested
    if avg
        GLMakie.lines!(ax, t, mean(s; dims = 1)[:]; color = :black, linewidth = 2)
    end

    # plot response time if provided
    if !isnothing(rt) && rt / 1000 ∈ t
        GLMakie.vlines!(ax, rt / 1000; linewidth = 1.5, color = mono ? :black : :red)
    end

    # add legend if requested and not too many channels
    if leg && ch_n < 30
        axislegend(ax; position = :rt, colormap = pal)
    end

    return fig
end

"""
    plot_erp_topo(locs, t, s; <keyword arguments>)

Plot topographical maps of Event-Related Potentials/Fields.

# Arguments

- `locs::DataFrame`: channel locations with columns: `channel`, `labels`, `loc_radius`, `loc_theta`, `loc_x`, `loc_y`, ``loc_z`, `loc_radius_sph`, `loc_theta_sph`, `loc_phi_sph`
- `t::Vector{Float64}`: time points (in seconds)
- `s::Matrix{Float64}`: ERP/ERF data, shape (channels, samples)
- `rt::Union{Nothing, Real}=nothing`: response time (in milliseconds)
- `clabels::Vector{String}=string.(1:size(s, 1))`: channel labels (default: auto-generated)
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `yrev::Bool=false`: if `true`, reverse the y-axis
- `cart::Bool=false`: if `true`, use Cartesian coordinates, otherwise use polar coordinates for XY plane and spherical coordinates for XZ and YZ planes
- `head::Bool=true`: if `true`, draw head outline
- `zl::Bool`: if `true`, draw vertical line at t = 0
- `mono::Bool=false`: if `true`, use a monochrome palette

# Returns

- `GLMakie.Figure`: the plotted figure
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

    # validate
    size(s, 2) == length(t) ||
        throw(ArgumentError("Signal matrix columns must match time vector length"))
    size(s, 1) == nrow(locs) ||
        throw(ArgumentError("Number of channels must match channel locations"))
    all(isfinite, s) ||
        throw(ArgumentError("Signal contains non-finite values"))

    # select color palette
    pal = mono ? :grays : :darktest

    # determine plot size and marker size based on number of channels
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
            loc_x[idx], loc_y[idx] = pol2cart(locs.loc_radius[idx], locs.loc_theta[idx])
        end
    else
        loc_x = locs.loc_x
        loc_y = locs.loc_y
    end

    # create individual ERP plots for each channel
    pp_vec = GLMakie.Figure[]
    pp_full_vec = GLMakie.Figure[]
    for idx in axes(s, 1)
        pp = GLMakie.Figure(; size = marker_size, figure_padding = 0)
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

        # plot ERP with zero line
        GLMakie.hlines!(ax, 0; color = :black, linewidth = 1)
        if zl
            GLMakie.vlines!(ax, 0; color = :gray, linestyle = :dash, linewidth = 1)
        end
        GLMakie.lines!(ax, t, s[idx, :]; linewidth = 1, color = :black)
        # plot response time if provided
        if !isnothing(rt) && rt / 1000 ∈ t
            GLMakie.vlines!(ax, rt / 1000; linewidth = 1.5, color = mono ? :black : :red)
        end
        push!(pp_vec, pp)
        pp_full = plot_erp(
            t,
            s[idx, :];
            xlabel = xlabel,
            ylabel = ylabel,
            title = title,
            rt = rt,
            yrev = yrev,
        )
        push!(pp_full_vec, pp_full)
    end

    # prepare main topographical plot
    GLMakie.activate!(; title = "plot_erp()")
    fig = GLMakie.Figure(;
        size = plot_size,
        figure_padding = 0,
    )

    # create axis with customizable properties
    ax = GLMakie.Axis(
        fig[1, 1];
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

    # Draw head outline if requested
    if head
        draw_head_outline!(ax)
    end

    # Draw channel markers with embedded ERP plots
    for (idx, (x, y)) in enumerate(zip(loc_x, loc_y))
        io = IOBuffer()
        show(io, MIME"image/png"(), pp_vec[idx])
        marker_img = FileIO.load(io)
        GLMakie.scatter!(
            ax,
            x, y;
            marker = marker_img,
            markersize = marker_size,
            markerspace = :pixel,
        )
    end

    # generate clickable areas
    loc_x_range = Tuple{Float64, Float64}[]
    loc_y_range = Tuple{Float64, Float64}[]
    for idx in eachindex(loc_x)
        push!(loc_x_range, (loc_x[idx] - 0.15, loc_x[idx] + 0.15))
        push!(loc_y_range, (loc_y[idx] - 0.1, loc_y[idx] + 0.1))
    end

    # mouse events
    on(events(fig).mousebutton) do event
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

    return fig
end

"""
    plot_erp_stack(t, s; <keyword arguments>)

Plot stacked Event-Related Potentials/Fields.

# Arguments

- `t::AbstractVector`: time points (in seconds)
- `s::AbstractMatrix`: ERP/ERF data, shape (epochs, samples)
- `rt::Union{Nothing, AbstractVector}=nothing`: response time for each epoch; if provided, the response time line will be plotted over the `:stack` plot
- `clabels::Vector{String}=string.(1:size(s, 1))`: channel labels (default: auto-generated)
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `cb::Bool=true`: if `true`, show color bar
- `cb_title::String=""`: color bar title
- `smooth::Bool=false`: if `true`, apply Gaussian blur smoothing
- `ks::Int64=3`: smoothing kernel size; larger kernel means more smoothing
- `zl::Bool`: if `true`, draw vertical line at t = 0
- `mono::Bool=false`: if `true`, use a monochrome palette

# Returns

- `GLMakie.Figure`: the plotted figure
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

    # validate
    length(t) == size(s, 2) ||
        throw(
            ArgumentError(
                "Number of s columns ($(size(s, 2))) must equal length of t ($(length(t))).",
            ),
        )
    !isnothing(rt) && (length(rt) == size(s, 1)) ||
        throw(
            ArgumentError(
                "Length of the rt vector must equal number of ERP epochs ($(size(s, 1))).",
            ),
        )

    # select color palette
    pal = mono ? :grays : :darktest

    # apply smoothing if requested
    if smooth
        s = imfilter(s, Kernel.gaussian(ks))
    end

    # prepare plot
    GLMakie.activate!(; title = "plot_erp()")
    plot_size = size(s, 1) <= 64 ? (900, 600) : (900, 900)
    fig = GLMakie.Figure(; size = plot_size)

    # create axis with customizable properties
    ax = GLMakie.Axis(
        fig[1, 1];
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        xticks = LinearTicks(10),
        xminorticksvisible = true,
        xminorticks = IntervalsBetween(10),
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        yticks = (1:size(s, 1), size(s, 1) <= 30 ? clabels : clabels[1:5:end]),
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

    # create heatmap of ERP data
    hm = GLMakie.heatmap!(ax, t, axes(s, 1), rotr90(s); colormap = pal)

    # draw zero line if requested
    if zl
        GLMakie.vlines!(ax, 0; color = :white, linestyle = :dash, linewidth = 2)
    end

    # plot response times if provided
    if !isnothing(rt)
        for (i, rt_val) in enumerate(rt)
            if rt_val / 1000 ∈ t
                GLMakie.lines!(
                    ax,
                    [rt_val / 1000, rt_val / 1000],
                    [i - 0.5, i + 0.5];
                    linewidth = 1.5,
                    color = mono ? :black : :red,
                )
            end
        end
    end

    # add colorbar if requested
    if cb
        Colorbar(fig[1, 2], hm; label = cb_title, labelsize = 16)
    end

    return fig
end

"""
    plot_gfp(t, s; <keyword arguments>)

Plot Global Field Power (GFP).

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: time points (in seconds)
- `g::AbstractVector`: GFP values
- `rt::Union{Nothing, Real}=nothing`: response time (in milliseconds)
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `zl::Bool`: if `true`, draw vertical line at t = 0
- `mono::Bool=false`: if `true`, use a monochrome palette

# Returns

- `GLMakie.Figure`: the plotted figure
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

    # validate
    length(t) == length(g) ||
        throw(ArgumentError("Time vector and GFP must have equal length."))
    all(isfinite, g) ||
        throw(ArgumentError("GFP contains non-finite values."))
    all(g .>= 0) ||
        throw(ArgumentError("GFP values must be non-negative."))

    # prepare plot
    GLMakie.activate!(; title = "plot_gfp()")
    plot_size = (900, 450)
    fig = GLMakie.Figure(; size = plot_size)

    # create axis with customizable properties
    ax = GLMakie.Axis(
        fig[1, 1];
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
    GLMakie.ylims!(ax, 0, maximum(g) * 1.15)
    ax.titlesize = 18
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    # draw zero line if requested
    if zl
        GLMakie.vlines!(ax, 0; color = :gray, linestyle = :dash, linewidth = 2)
    end

    # plot GFP as filled band with outline
    GLMakie.band!(ax, t, 0, g; color = :gray)
    GLMakie.lines!(ax, t, g; color = :black, linewidth = 1)

    # plot response time if provided
    if !isnothing(rt) && rt / 1000 ∈ t
        GLMakie.vlines!(ax, rt / 1000; linewidth = 1.5, color = mono ? :black : :red)
    end

    return fig
end

"""
    plot_erp(obj; <keyword arguments>)

Plot Event-Related Potential/Field (ERP/ERF) from a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object containing ERP/ERF data
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `tm::Union{Nothing, Int64, Vector{Int64}}=nothing`: time markers (in milliseconds) to plot as vertical lines, useful for adding topoplots at these time points
- `xlabel::String="default"`: x-axis label
- `ylabel::String="default"`: y-axis label
- `title::String="default"`: plot title
- `cb::Bool=true`: if `true`, show color bar
- `cb_title::String="default"`: color bar title
- `peaks::Bool=true`: draw peak markers
- `leg::Bool=true`: if `true`, add legend with channel labels
- `type::Symbol=:normal`: plot type:
    - `:normal`: standard ERP plot
    - `:gfp`: plot Global Field Power
    - `:stack`: stacked epochs/channels
    - `:topo`: topographical plot
- `yrev::Bool=false`: if `true`, reverse the y-axis
- `avg::Bool=true`: if `true`, plot averaged ERP
- `ci95::Bool=false`: if `true`, plot mean and ±95% CI
- `smooth::Bool=false`: if `true`, apply Gaussian blur smoothing
- `ks::Int64=3`: smoothing kernel size; larger kernel means more smoothing
- `rt::Union{Nothing, Real, AbstractVector}=nothing`: response time for each epoch; if provided, the response time line will be plotted over the `:stack` plot
- `sort_epochs::Bool=false`:: sort epochs by rt vector
- `zl::Bool`: if `true`, draw vertical line at t = 0
- `mono::Bool=false`: if `true`, use a monochrome palette
- `gui::Bool=false`: ignored

# Returns

- `GLMakie.Figure`: the plotted figure
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

    # validate
    _check_datatype(obj, ["erp", "erf"])
    _check_var(type, [:normal, :topo, :stack, :gfp], "type")

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ?
        get_channel(obj; ch = ch, exclude = "bad") :
        get_channel(obj; ch = ch, exclude = "")
    length(ch) > 1 && length(unique(obj.header.recording[:channel_type][ch])) > 1 ||
        throw(ArgumentError("All channels must be of the same type."))
    length(ch) > 1 && (eavg = false)

    # for GFP plots, we need at least 2 channels
    if type === :gfp
        length(ch) > 1 || throw(ArgumentError("More than 1 channel must be selected."))
    end

    # get channel labels and units
    units = _ch_units(obj, labels(obj)[ch[1]])
    clabels = labels(obj)[ch]

    # extract data based on number of channels
    ep_n = nepochs(obj) - 1
    if length(ch) == 1
        if type === :stack
            s = obj.data[ch[1], :, 2:end]' # transpose for stacked plot
        else
            s = obj.data[ch, :, 1][:] # single channel data
        end
    else
        s = obj.data[ch, :, 1] # multi-channel data
    end

    # get time vector
    t = obj.epoch_time

    if length(ch) == 1
        if type === :stack

            # prepare default labels/titles
            xl, yl, tt = _set_defaults(
                xlabel, ylabel, title,
                "Time [ms]", "Amplitude [$units]",
                "ERP amplitude, $(length(clabels)) channels, avgₑ: $ep_n",
            )
            cb_title = cb_title == "default" ? "Amplitude [$units]" : cb_title

            # handle sorting by response time if requested
            if sort_epochs && !isnothing(rt)
                rt_idx = sortperm(rt)
                rt = rt[rt_idx]
                s = s[rt_idx, :]
            end

            fig = plot_erp_stack(
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

            # prepare default labels/titles
            xl, yl, tt = _set_defaults(
                xlabel,
                ylabel,
                title,
                "Time [ms]",
                "Amplitude [$units]",
                "ERP amplitude, $(clabels[1]), avgₑ: $ep_n",
            )
            fig = plot_erp(
                t,
                s;
                xlabel = xl,
                ylabel = yl,
                title = tt,
                rt = rt,
                yrev = yrev,
                zl = zl,
                mono = mono,
            )
        end

    elseif type === :normal

        # turn off peaks if not averaging
        avg == false && (peaks = false)
        xl, yl, tt = _set_defaults(
            xlabel,
            ylabel,
            title,
            "Time [ms]",
            "Amplitude [$units]",
            "ERP amplitude, $(length(ch)) channels, avgₑ: $ep_n",
        )
        fig = plot_erp(
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
            xlabel,
            ylabel,
            title,
            "Time [ms]",
            "",
            "ERP amplitude, $(length(ch)) channels, avgₑ: $ep_n",
        )
        fig = plot_erp_stack(
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
        g = erp_gfp(obj; ch = clabels)
        xl, yl, tt = _set_defaults(
            xlabel,
            ylabel,
            title,
            "Time [ms]",
            "GFP [$units]",
            "Global Field Power, $(length(ch)) channels, avgₑ: $ep_n",
        )
        fig = plot_gfp(
            t,
            g;
            xlabel = xl,
            ylabel = yl,
            title = tt,
            rt = rt,
            zl = zl,
            mono = mono,
        )

    elseif type === :topo
        xl, yl, tt = _set_defaults(
            xlabel,
            ylabel,
            title,
            "Time [ms]",
            "Amplitude [$units]",
            "ERP amplitude, avgₑ: $ep_n",
        )
        _check_ch_locs(ch, labels(obj), obj.locs.label)
        length(unique(obj.header.recording[:channel_type][ch])) == 1 ||
            throw(
                ArgumentError(
                    "For multi-channel topo plot all channels must be of the same type.",
                ),
            )
        _has_locs(obj)
        chs = intersect(obj.locs.label, labels(obj)[ch])
        locs = Base.filter(:label => in(chs), obj.locs)
        _check_ch_locs(ch, labels(obj), obj.locs.label)
        fig = plot_erp_topo(
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

    # add time markers if requested
    if !isnothing(tm)
        for (i, marker) in enumerate(tm)
            marker / 1000 >= t[1] && marker / 1000 <= t[end] ||
                throw(ArgumentError("Time marker $marker is out of epoch range"))
            marker_idx = vsearch(marker / 1000, t)
            GLMakie.vlines!(
                fig[1, 1],
                t[marker_idx];
                linewidth = 0.5,
                color = :black,
            )
        end
    end

    # add peak markers if requested
    if peaks
        if length(ch_indices) == 1 && type === :normal
            pp = erp_peaks(obj)
            pos_time = t[pp[ch_indices[1], 1]][1] * 1000
            pos_amp = obj.data[ch_indices[1], pp[ch_indices[1], 1], 1][1]
            neg_time = t[pp[ch_indices[1], 2]][1] * 1000
            neg_amp = obj.data[ch_indices[1], pp[ch_indices[1], 2], 1][1]

            GLMakie.scatter!(
                fig[1, 1],
                pos_time / 1000,
                pos_amp;
                marker = :xcross,
                color = mono ? :black : :red,
                markersize = 15,
            )
            GLMakie.scatter!(
                fig[1, 1],
                neg_time / 1000,
                neg_amp;
                marker = :xcross,
                color = mono ? :black : :blue,
                markersize = 15,
            )

            @info "Positive peak: $(round(pos_time, digits = 0)) ms, $(round(pos_amp, digits = 2)) $units"
            @info "Negative peak: $(round(neg_time, digits = 0)) ms, $(round(neg_amp, digits = 2)) $units"

        elseif length(ch_indices) > 1 && type === :normal

            # for multi-channel average
            mep_tmp = mean(obj.data[ch_indices, :, 1]; dims = 1)[:, :, :]
            obj_tmp = keep_channel(obj; ch = clabels[1])
            obj_tmp.data = mep_tmp
            pp = erp_peaks(obj_tmp)

            pos_time = t[pp[1, 1]] * 1000
            pos_amp = mep_tmp[pp[1, 1]]
            neg_time = t[pp[1, 2]] * 1000
            neg_amp = mep_tmp[pp[1, 2]]

            GLMakie.scatter!(
                fig[1, 1],
                pos_time / 1000,
                pos_amp;
                marker = :xcross,
                color = mono ? :black : :red,
                markersize = 15,
            )
            GLMakie.scatter!(
                fig[1, 1],
                neg_time / 1000,
                neg_amp;
                marker = :xcross,
                color = mono ? :black : :blue,
                markersize = 15,
            )

            @info "Positive peak: $(round(pos_time, digits = 0)) ms, $(round(pos_amp, digits = 2)) $units"
            @info "Negative peak: $(round(neg_time, digits = 0)) ms, $(round(neg_amp, digits = 2)) $units"
        end
    end

    return fig
end
