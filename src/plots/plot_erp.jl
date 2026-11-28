export plot_erp
export plot_erp_topo
export plot_erp_stack
export plot_gfp

"""
    plot_erp(t, s; <keyword arguments>)

Plot Event-Related Potential/Field (single channel).

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: time values (in seconds; the time axis is displayed in milliseconds)
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

    # time is displayed in milliseconds
    tms = t .* 1000

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
        _TIME_AXIS_KWARGS...,
        yticks = LinearTicks(10),
        yminorticksvisible = true,
        yminorticks = IntervalsBetween(10),
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        yreversed = yrev,
        _AXIS_LOCK_KWARGS...,
    )
    GLMakie.ylims!(ax, yrev ? reverse(_ylims(s) .* 1.5) : (_ylims(s) .* 1.5))
    _style_axis!(ax)

    # draw zero line if requested
    zl && GLMakie.vlines!(ax, 0; color = :gray, linestyle = :dash, linewidth = 2)

    # plot ERP signal
    GLMakie.lines!(ax, tms, s; color = :black, linewidth = 1)

    # plot response time if provided
    if !isnothing(rt) && t[1] <= rt / 1000 <= t[end]
        GLMakie.vlines!(
            ax,
            rt;
            linewidth = 1.5,
            color = :black,
        )
    end

    return fig
end

"""
    plot_erp(t, s; <keyword arguments>)

Plot multi-channel Event-Related Potential/Field with optional averaging and confidence intervals.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: time values (in seconds; the time axis is displayed in milliseconds)
- `s::AbstractMatrix`: signal data, shape (channels, samples)
- `rt::Union{Nothing, Real}=nothing`: response time (in milliseconds)
- `clabels::Vector{String}=string.(1:size(s, 1))`: channel labels (default: auto-generated)
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `yrev::Bool=false`: if `true`, reverse the y-axis
- `avg::Bool=true`: if `true`, plot averaged ERP
- `ci95::Bool=false`: if `true`, plot mean and ±95% CI
- `leg::Bool=true`: if `true`, draw a legend with channel labels (outside the plot, on the right; not drawn if it does not fit: too many channels)
- `zl::Bool`: if `true`, draw vertical line at t = 0
- `mono::Bool=false`: if `true`, use a monochrome palette
- `bad::AbstractVector{Bool}=falses(size(s, 1))`: bad channels, drawn in light gray with a dashed line and labeled "(bad)"
- `exclude_bad::Bool=false`: if `true`, the bad channels are not used for the average and the 95% CI

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
    bad::AbstractVector{Bool} = falses(size(s, 1)),
    exclude_bad::Bool = false,
)::GLMakie.Figure
    # validate
    size(s, 2) == length(t) ||
        throw(ArgumentError("Signal matrix columns must match time vector length"))
    length(clabels) == size(s, 1) ||
        throw(ArgumentError("Number of channel labels must match number of channels"))
    length(bad) == size(s, 1) ||
        throw(ArgumentError("Length of bad must match number of channels"))
    # channels used for the average and the 95% CI
    s_avg = exclude_bad && any(bad) ? s[.!bad, :] : s
    size(s_avg, 1) > 0 || throw(ArgumentError("All channels are bad."))
    all(isfinite, s) ||
        throw(ArgumentError("Signal contains non-finite values"))

    # time is displayed in milliseconds
    tms = t .* 1000

    # select color palette
    pal = _palette(:darktest; mono = mono)

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
        _TIME_AXIS_KWARGS...,
        yticks = LinearTicks(10),
        yminorticksvisible = true,
        yminorticks = IntervalsBetween(10),
        yreversed = yrev,
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        _AXIS_LOCK_KWARGS...,
    )
    GLMakie.ylims!(ax, yrev ? reverse(_ylims(s) .* 1.5) : (_ylims(s) .* 1.5))
    _style_axis!(ax)

    # draw zero line if requested
    zl && GLMakie.vlines!(ax, 0; color = :gray, linestyle = :dash, linewidth = 2)

    # plot ERPs
    if ci95
        avg = false
        leg = false

        # calculate mean and 95% CI
        msci95_data = NeuroAnalyzer._msci95(s_avg)
        s_m = msci95_data.sm
        s_u = msci95_data.ul
        s_l = msci95_data.ll

        # plot confidence interval
        GLMakie.Makie.band!(ax, tms, s_l, s_u; alpha = 0.25, color = :grey)

        # plot mean signal
        GLMakie.Makie.lines!(ax, tms, s_m; color = :black, linewidth = 2)
    else
        # plot individual channels
        cmap = _trace_colors(ch_n; mono = mono)
        for idx = 1:ch_n
            GLMakie.lines!(
                ax,
                tms,
                s[idx, :];
                color = bad[idx] ? :lightgray : cmap[idx],
                colormap = pal,
                linewidth = 1,
                linestyle = bad[idx] ? :dash : :solid,
                alpha = avg ? 0.25 : 1.0,
                label = bad[idx] ? clabels[idx] * " (bad)" : clabels[idx],
            )
        end
    end

    # plot averaged signal if requested
    avg && GLMakie.lines!(ax, tms, mean(s_avg; dims = 1)[:]; color = :black, linewidth = 2)

    # plot response time if provided
    if !isnothing(rt) && t[1] <= rt / 1000 <= t[end]
        GLMakie.vlines!(ax, rt; linewidth = 1.5, color = :black)
    end

    # add legend if requested (outside the plot; not drawn if it does not fit)
    leg && _side_legend!(fig, ax)

    return fig
end

"""
    plot_erp_topo(locs, t, s; <keyword arguments>)

Plot topographical maps of Event-Related Potentials/Fields.

# Arguments

- `locs::DataFrame`: channel locations with columns: `channel`, `labels`, `loc_radius`, `loc_theta`, `loc_x`, `loc_y`, ``loc_z`, `loc_radius_sph`, `loc_theta_sph`, `loc_phi_sph`
- `t::Vector{Float64}`: time points (in seconds; the time axis is displayed in milliseconds)
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
- `gui::Bool=true`: if `true`, keep window open and interactive
- `bad::AbstractVector{Bool}=falses(size(s, 1))`: bad channels (one value per row of `locs`), drawn in light gray and labeled "(bad)"

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
    gui::Bool = true,
    bad::AbstractVector{Bool} = falses(size(s, 1)),
)::GLMakie.Figure
    # locations are stored in their original units: scale the montage to the unit sphere
    locs = _locs_unit(locs)
    # validate
    size(s, 2) == length(t) ||
        throw(ArgumentError("Signal matrix columns must match time vector length"))
    size(s, 1) == nrow(locs) ||
        throw(ArgumentError("Number of channels must match channel locations"))
    all(isfinite, s) ||
        throw(ArgumentError("Signal contains non-finite values"))

    # select color palette
    pal = _palette(:darktest; mono = mono)

    # number of channels
    ch_n = size(s, 1)

    # plot parameters
    if ch_n <= 64
        xl = 1.2
        yl = 1.2
    else
        xl = 1.3
        yl = 1.3
    end
    if gui
        if ch_n <= 64
            plot_size = (1000, 1000)
            marker_size = (150, 75)
        elseif ch_n <= 100
            plot_size = (1200, 1200)
            marker_size = (110, 55)
        else
            plot_size = (1400, 1400)
            marker_size = (90, 45)
        end
    else
        if ch_n <= 64
            plot_size = (500, 500)
            marker_size = (70, 35)
        elseif ch_n <= 100
            plot_size = (700, 700)
            marker_size = (60, 30)
        else
            plot_size = (900, 900)
            marker_size = (70, 35)
        end
    end

    # get locations
    loc_x = zeros(size(locs, 1))
    loc_y = zeros(size(locs, 1))
    if cart
        loc_x = locs.loc_x
        loc_y = locs.loc_y
    else
        for idx in axes(locs, 1)
            loc_x[idx], loc_y[idx] = pol2cart(locs.loc_radius[idx], locs.loc_theta[idx])
        end
    end

    # create individual ERP plots for each channel
    pp_vec = GLMakie.Figure[]
    pp_full_vec = GLMakie.Figure[]
    for idx in axes(s, 1)
        pp = GLMakie.Figure(; size = marker_size, figure_padding = 0)
        gl = GridLayout(pp[1, 1]; alignmode = Outside(5))
        ax = GLMakie.Axis(
            gl[1, 1];
            xlabel = "",
            ylabel = "",
            title = bad[idx] ? locs[idx, :label] * " (bad)" : locs[idx, :label],
            yreversed = yrev,
            xautolimitmargin = (0, 0),
            yautolimitmargin = (0.1, 0.1),
        )
        hidedecorations!(ax)
        ax.titlesize = 8
        GLMakie.autolimits!(ax)

        # plot ERP with zero line
        GLMakie.hlines!(ax, 0; color = :black, linewidth = 1)
        zl && GLMakie.vlines!(ax, 0; color = :gray, linestyle = :dash, linewidth = 1)
        GLMakie.lines!(ax, t .* 1000, s[idx, :]; linewidth = 1, color = bad[idx] ? :lightgray : :black)

        # plot response time if provided
        if !isnothing(rt) && t[1] <= rt / 1000 <= t[end]
            GLMakie.vlines!(ax, rt; linewidth = 1.5, color = :black)
        end

        # draw black box around the entire figure
        GLMakie.Box(pp[:, :];
            color = :transparent,
            strokecolor = :black,
            strokewidth = 1.5,
        )

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
        _AXIS_LOCK_KWARGS...,
    )
    GLMakie.xlims!(ax, (-xl, xl))
    GLMakie.ylims!(ax, (-yl, yl))
    hidespines!(ax)
    hidedecorations!(ax)
    ax.titlesize = gui ? 18 : 12

    # draw head outline if requested
    head && _draw_head_outline!(ax; lw = 2)

    # draw channel markers with embedded ERP plots
    for (idx, (x, y)) in enumerate(zip(loc_x, loc_y))
        io = IOBuffer()
        show(io, MIME"image/png"(), pp_vec[idx])
        marker_img = FileIO.load(io)
        GLMakie.scatter!(
            ax,
            x, y;
            marker = marker_img,
            markersize = (0.3, 0.15), # data units; matches your loc_x_range/loc_y_range
            markerspace = :data,
            # markerspace = :pixel,
        )
    end

    if gui
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

        screen = display(fig)
        GLMakie.GLFW.MaximizeWindow(screen.glscreen)
        wait(screen)
    end

    return fig
end

"""
    plot_erp_stack(t, s; <keyword arguments>)

Plot stacked Event-Related Potentials/Fields.

# Arguments

- `t::AbstractVector`: time points (in seconds; the time axis is displayed in milliseconds)
- `s::AbstractMatrix`: ERP/ERF data, shape (epochs, samples)
- `rt::Union{Nothing, AbstractVector}=nothing`: response time (in milliseconds) for each epoch; if provided, the response time line will be plotted over the stacked epochs (outside the epoch time range: not plotted)
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
    if !isnothing(rt)
        length(rt) == size(s, 1) ||
            throw(
                ArgumentError(
                    "Length of the rt vector must equal number of ERP epochs ($(size(s, 1))).",
                ),
            )
    end

    # select color palette
    pal = _palette(:darktest; mono = mono)

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
        _TIME_AXIS_KWARGS...,
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        yticks = size(s, 1) <= 30 ? (1:size(s, 1), clabels) :
                 (1:5:size(s, 1), clabels[1:5:end]),
        yticklabelsize = size(s, 1) <= 64 ? 8 : 5,
        _AXIS_LOCK_KWARGS...,
    )
    _style_axis!(ax)

    # create heatmap of ERP data
    # heatmap!(x, y, z) expects z of shape (length(x), length(y)): row i of s at y = i
    # (time displayed in milliseconds)
    hm = GLMakie.heatmap!(ax, t .* 1000, axes(s, 1), permutedims(s); colormap = pal)

    # draw zero line if requested
    zl && GLMakie.vlines!(ax, 0; color = :white, linestyle = :dash, linewidth = 2)

    # plot response times if provided
    if !isnothing(rt)
        for (i, rt_val) in enumerate(rt)
            if t[1] <= rt_val / 1000 <= t[end]
                GLMakie.lines!(
                    ax,
                    [rt_val, rt_val],
                    [i - 0.5, i + 0.5];
                    linewidth = 3,
                    color = :black,
                )
            end
        end
    end

    # add colorbar if requested
    cb && Colorbar(fig[1, 2], hm; label = cb_title, labelsize = 16, ticks = _cb_tickfun())

    return fig
end

"""
    plot_gfp(t, s; <keyword arguments>)

Plot Global Field Power (GFP).

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: time points (in seconds; the time axis is displayed in milliseconds)
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
        _TIME_AXIS_KWARGS...,
        yticks = LinearTicks(10),
        yminorticksvisible = true,
        yminorticks = IntervalsBetween(2),
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        _AXIS_LOCK_KWARGS...,
    )
    GLMakie.ylims!(ax, 0, maximum(g) * 1.15)
    _style_axis!(ax)

    # draw zero line if requested
    zl && GLMakie.vlines!(ax, 0; color = :gray, linestyle = :dash, linewidth = 2)

    # plot GFP as filled band with outline
    # (time displayed in milliseconds)
    GLMakie.band!(ax, t .* 1000, 0, g; color = :gray)
    GLMakie.lines!(ax, t .* 1000, g; color = :black, linewidth = 1)

    # plot response time if provided
    if !isnothing(rt) && t[1] <= rt / 1000 <= t[end]
        GLMakie.vlines!(ax, rt; linewidth = 1.5, color = :black)
    end

    return fig
end

"""
    plot_erp(obj; <keyword arguments>)

Plot Event-Related Potential/Field (ERP/ERF) from a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object containing ERP/ERF data
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `exclude_bad_channels::Bool=NeuroAnalyzer.exclude_bad_channels`: if `true`, bad channels are not used for the average, the 95% CI (`type = :normal`) and the GFP (`type = :gfp`); all selected channels are shown, bad channels are marked (default: the `exclude_bad_channels` preference, see `set_exclude_bad_channels()`)
- `tm::Union{Nothing, Int64, Vector{Int64}}=nothing`: time markers (in milliseconds) to plot as vertical lines
- `tm_topo::Bool=false`: if `true`, draw small topographical maps at the time markers `tm` above the plot, with lines pointing to these time points (only for `type = :normal` and `type = :gfp`); the maps are evenly spaced, show all channels of the type of the plotted channels that have locations (bad channels are marked; with `exclude_bad_channels` they are not interpolated) and share one symmetric color range (colorbar on the right)
- `xlabel::String="default"`: x-axis label
- `ylabel::String="default"`: y-axis label
- `title::String="default"`: plot title
- `cb::Bool=true`: if `true`, show color bar
- `cb_title::String="default"`: color bar title
- `peaks::Bool=true`: if `true`, draw peak markers
- `leg::Bool=true`: if `true`, draw a legend with channel labels (outside the plot, on the right; not drawn if it does not fit: too many channels)
- `type::Symbol=:normal`: plot type:
    - `:normal`: standard ERP plot
    - `:gfp`: plot Global Field Power
    - `:stack`: stacked epochs/channels
    - `:topo`: topographical plot
- `yrev::Bool=false`: if `true`, reverse the y-axis
- `avg::Bool=true`: if `true`, plot averaged ERP
- `ci95::Bool=false`: if `true`, plot the mean and its 95% CI (over channels or epochs, see `ci_over`)
- `ci_over::Symbol=:channels`: with several channels, the 95% CI is calculated over:
    - `:channels`: the channels of the ERP (epoch 1)
    - `:epochs`: the epochs following the ERP (the single trials of `average_epochs()` or the subjects of `study_erp()`), each the mean of the channels (with `exclude_bad_channels` without the bad channels)

  with one channel, the CI is always calculated over the epochs
- `smooth::Bool=false`: if `true`, apply Gaussian blur smoothing
- `ks::Int64=3`: smoothing kernel size; larger kernel means more smoothing
- `rt::Union{Nothing, Real, AbstractVector}=nothing`: response time (in milliseconds): a single value is drawn as a vertical line (`:normal`, `:gfp`, `:topo`), a vector (one value per epoch) over the stacked epochs (`:stack`); values outside the epoch time range are not plotted
- `sort_epochs::Bool=false`:: sort epochs by rt vector
- `zl::Bool`: if `true`, draw vertical line at t = 0
- `mono::Bool=false`: if `true`, use a monochrome palette
- `gui::Bool=true`: if `true`, keep window open and interactive
- `ch_locs::Bool=false`: if `true`, draw a head outline with the locations of the plotted channels in the colors of their traces (black for a single channel) next to the plot (only for `type = :normal`)
- `ch_locs_pos::Symbol=:outside`: position of the channel locations inset:
    - `:outside`: to the right of the plot
    - `:rt`, `:lt`, `:rb`, `:lb`: inside the plot, in the right-top, left-top, right-bottom or left-bottom corner

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_erp(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    tm::Union{Nothing, Int64, Vector{Int64}} = nothing,
    tm_topo::Bool = false,
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
    ci_over::Symbol = :channels,
    smooth::Bool = false,
    ks::Int64 = 3,
    rt::Union{Nothing, Real, AbstractVector} = nothing,
    sort_epochs::Bool = false,
    zl::Bool = true,
    mono::Bool = false,
    ch_locs::Bool = false,
    ch_locs_pos::Symbol = :outside,
    gui::Bool = true,
    exclude_bad_channels::Bool = NeuroAnalyzer.exclude_bad_channels,
)::GLMakie.Figure
    _check_var(ch_locs_pos, [:outside, :rt, :lt, :rb, :lb], "ch_locs_pos")

    # validate
    _check_datatype(obj, ["erp", "erf"])
    _check_var(type, [:normal, :topo, :stack, :gfp], "type")
    _check_var(ci_over, [:channels, :epochs], "ci_over")
    # the CI over epochs needs ≥ 2 epochs after the ERP
    ci95 && type === :normal && (length(ch) == 1 || ci_over === :epochs) && nepochs(obj) < 3 &&
        throw(ArgumentError("ci95 over epochs requires ≥ 2 epochs after the ERP (epoch 1)."))
    if tm_topo
        isnothing(tm) && throw(ArgumentError("tm must be specified for tm_topo."))
        type in (:normal, :gfp) ||
            throw(ArgumentError("tm_topo can only be used for type = :normal or :gfp."))
    end

    # all selected channels are shown (bad channels are marked); with exclude_bad_channels
    # the bad channels are not used for the average, the 95% CI and the GFP
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))
    bad = obj.header.recording[:bad_channel][ch]
    if length(ch) > 1
        length(unique(obj.header.recording[:channel_type][ch])) == 1 ||
            throw(ArgumentError("All channels must be of the same type."))
    end

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
                "ERP amplitude, $(clabels[1])$(bad[1] ? " (bad)" : ""), avgₑ: $ep_n" *
                (ci95 ? ", 95% CI over epochs" : ""),
            )
            fig = if ci95
                # mean and 95% CI over the epochs following the ERP
                plot_erp(
                    t,
                    Matrix(obj.data[ch[1], :, 2:end]');
                    xlabel = xl,
                    ylabel = yl,
                    title = tt,
                    clabels = string.(1:ep_n),
                    rt = rt,
                    yrev = yrev,
                    avg = false,
                    ci95 = true,
                    leg = false,
                    zl = zl,
                    mono = mono,
                )
            else
                plot_erp(
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
            ch_locs && _ch_locs_inset!(fig, obj, ch, [:black]; pos = ch_locs_pos)
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
            "ERP amplitude, $(length(ch)) channels, avgₑ: $ep_n" *
            (ci95 ? ", 95% CI over $(ci_over)" : ""),
        )
        fig = if ci95 && ci_over === :epochs
            # each epoch following the ERP: the mean of the (good) channels
            ch_m = exclude_bad_channels && !all(bad) ? ch[.!bad] : ch
            S = dropdims(mean(obj.data[ch_m, :, 2:end]; dims = 1); dims = 1)
            plot_erp(
                t,
                Matrix(S');
                xlabel = xl,
                ylabel = yl,
                title = tt,
                clabels = string.(1:ep_n),
                rt = rt,
                yrev = yrev,
                avg = false,
                ci95 = true,
                leg = false,
                zl = zl,
                mono = mono,
            )
        else
            plot_erp(
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
                bad = bad,
                exclude_bad = exclude_bad_channels,
            )
        end
        ch_locs && _ch_locs_inset!(
            fig, obj, ch, _trace_colors(length(ch); mono = mono, black = ci95);
            pos = ch_locs_pos,
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
            clabels = [b ? l * " (bad)" : l for (l, b) in zip(clabels, bad)],
            cb = cb,
            cb_title = cb_title,
            smooth = smooth,
            ks = ks,
            zl = zl,
            mono = mono,
        )

    elseif type === :gfp
        g = erp_gfp(obj; ch = clabels, exclude_bad_channels = exclude_bad_channels)
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
        locs = Base.filter(:label => in(chs), _locs_unit(obj.locs))
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
            gui = gui,
            bad = in(Set(labels(obj)[ch][bad])).(locs[!, :label]),
        )
    end

    # add time markers if requested
    if !isnothing(tm)
        for (i, marker) in enumerate(tm)
            marker / 1000 >= t[1] && marker / 1000 <= t[end] ||
                throw(ArgumentError("Time marker $marker is out of epoch range"))
            marker_idx = vsearch(marker / 1000, t)
            GLMakie.vlines!(
                _main_axis(fig),
                t[marker_idx] * 1000;
                linewidth = 0.5,
                color = :black,
            )
        end
    end

    # add peak markers if requested
    if peaks
        if length(ch) == 1 && type === :normal
            pp = erp_peaks(obj)
            pos_time = t[pp[ch[1], 1]][1] * 1000
            pos_amp = obj.data[ch[1], pp[ch[1], 1], 1][1]
            neg_time = t[pp[ch[1], 2]][1] * 1000
            neg_amp = obj.data[ch[1], pp[ch[1], 2], 1][1]

            GLMakie.scatter!(
                _main_axis(fig),
                pos_time,
                pos_amp;
                marker = :xcross,
                color = mono ? :black : :red,
                markersize = 15,
            )
            GLMakie.scatter!(
                _main_axis(fig),
                neg_time,
                neg_amp;
                marker = :xcross,
                color = mono ? :black : :blue,
                markersize = 15,
            )

            _info("Positive peak: $(round(pos_time, digits = 0)) ms, $(round(pos_amp, digits = 2)) $units")
            _info("Negative peak: $(round(neg_time, digits = 0)) ms, $(round(neg_amp, digits = 2)) $units")

        elseif length(ch) > 1 && type === :normal

            # for multi-channel average
            mep_tmp = mean(obj.data[ch, :, 1]; dims = 1)[:, :, :]
            obj_tmp = keep_channel(obj; ch = clabels[1])
            obj_tmp.data = mep_tmp
            pp = erp_peaks(obj_tmp)

            pos_time = t[pp[1, 1]] * 1000
            pos_amp = mep_tmp[pp[1, 1]]
            neg_time = t[pp[1, 2]] * 1000
            neg_amp = mep_tmp[pp[1, 2]]

            GLMakie.scatter!(
                _main_axis(fig),
                pos_time,
                pos_amp;
                marker = :xcross,
                color = mono ? :black : :red,
                markersize = 15,
            )
            GLMakie.scatter!(
                _main_axis(fig),
                neg_time,
                neg_amp;
                marker = :xcross,
                color = mono ? :black : :blue,
                markersize = 15,
            )

            _info("Positive peak: $(round(pos_time, digits = 0)) ms, $(round(pos_amp, digits = 2)) $units")
            _info("Negative peak: $(round(neg_time, digits = 0)) ms, $(round(neg_amp, digits = 2)) $units")
        end
    end

    # topographical maps at the time markers, above the plot (drawn last: the plot moves
    # to the second layout row)
    tm_topo && _erp_topo_row!(
        fig, obj, ch, t, tm isa Int64 ? [tm] : tm;
        units = units, mono = mono, exclude_bad_channels = exclude_bad_channels,
    )

    return fig
end

# small topographical maps at the time points `tm` (in milliseconds) in a row above the
# plot of `fig`, connected by lines with the time points on the time axis of the plot;
# channels: all channels of the type of the plotted channels `ch` that have locations
function _erp_topo_row!(
    fig::GLMakie.Figure,
    obj::NeuroAnalyzer.NEURO,
    ch::AbstractVector{<:Integer},
    t::AbstractVector,
    tm::Vector{Int64};
    units::String,
    mono::Bool,
    exclude_bad_channels::Bool,
)::Nothing
    _has_locs(obj)
    ch_type = obj.header.recording[:channel_type][ch[1]]
    t_lbl = labels(obj)[obj.header.recording[:channel_type] .== ch_type]
    locs = Base.filter(:label => in(Set(t_lbl)), _locs_unit(obj.locs))
    # data rows of the location rows (matched by label)
    ch_idx = [findfirst(==(l), labels(obj)) for l in locs[!, :label]]
    bad_rows = findall(obj.header.recording[:bad_channel][ch_idx])
    ch_rows = exclude_bad_channels ? setdiff(1:length(ch_idx), bad_rows) :
              collect(1:length(ch_idx))
    length(ch_rows) >= 3 || throw(
        ArgumentError(
            "tm_topo requires ≥ 3 channels of type $ch_type with locations (not bad).",
        ),
    )

    # values at the time points (averaged epoch), one color range for all maps
    t_idx = [vsearch(m / 1000, t) for m in tm]
    v = obj.data[ch_idx, t_idx, 1]
    v_max = maximum(abs, v[ch_rows, :])
    v_max > 0 || (v_max = 1.0)
    crange = (-v_max, v_max)
    pal = _palette(:bluesreds; mono = mono)

    # map size [px]; the figure is made higher by the row of maps
    w = 110
    fw, fh = fig.scene.viewport[].widths
    GLMakie.resize!(fig.scene, round(Int64, fw), round(Int64, fh + w + 30))

    # (bottom margin: room for the lines between the maps and the plot)
    top = GLMakie.GridLayout(fig[0, 1]; alignmode = Outside(0, 0, 20, 0))
    # the plot title is moved above the maps (the lines would cross it)
    ax = _main_axis(fig)
    r = 1
    if !isempty(ax.title[])
        GLMakie.Label(
            top[1, 1:(length(tm) + 1)], ax.title[];
            fontsize = ax.titlesize[], font = ax.titlefont[], tellwidth = false,
        )
        ax.title = ""
        r = 2
    end
    axs = GLMakie.Axis[]
    for (i, m) in enumerate(tm)
        pl = plot_topo(
            v[ch_rows, i];
            locs      = locs,
            ch        = ch_rows,
            bad       = bad_rows,
            cb        = false,
            ps        = :m,
            nmethod   = :none,
            crange    = crange,
            mono      = mono,
            ch_labels = false,
        )
        ax_t = GLMakie.Axis(
            top[r, i];
            height          = w,
            aspect          = DataAspect(),
            title           = "$m ms",
            titlesize       = 11,
            titlegap        = 2,
            titlefont       = :regular,
            backgroundcolor = :transparent,
            _AXIS_LOCK_KWARGS...,
        )
        hidedecorations!(ax_t)
        hidespines!(ax_t)
        GLMakie.image!(ax_t, rotr90(_figure_image(pl)))
        push!(axs, ax_t)
    end
    GLMakie.Colorbar(
        top[r, length(tm) + 1];
        colormap      = pal,
        limits        = crange,
        label         = "[$units]",
        labelsize     = 10,
        ticklabelsize = 9,
        ticks         = _cb_tickfun(; n = 3),
        height        = round(Int64, 0.8 * w),
        width         = 8,
        tellheight    = false,
    )

    # lines from the bottom of each map to its time point at the top of the plot
    for (i, ax_t) in enumerate(axs)
        # (the time axis of the plot is in milliseconds)
        tp = t[t_idx[i]] * 1000
        pts = GLMakie.lift(ax_t.scene.viewport, ax.scene.viewport, ax.finallimits) do vt, vm, lim
            x = vm.origin[1] + (tp - lim.origin[1]) / lim.widths[1] * vm.widths[1]
            return [
                Point2f(vt.origin[1] + vt.widths[1] / 2, vt.origin[2]),
                Point2f(x, vm.origin[2] + vm.widths[2]),
            ]
        end
        GLMakie.lines!(fig.scene, pts; space = :pixel, color = :black, linewidth = 0.75)
    end

    return nothing
end
