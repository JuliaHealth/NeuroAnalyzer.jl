export plot_histogram
export plot_bar
export plot_line
export plot_box
export plot_violin
export plot_dots
export plot_paired
export plot_ci

# ---------------------------------------------------------------------------
# shared helpers
# ---------------------------------------------------------------------------

"""
Compute padded y-limits for a data array.
Lower bound is 0 when all values are positive, otherwise 1.5× the minimum.
"""
function _ylims_padded(s::AbstractArray)
    lo = minimum(s)
    hi = maximum(s)
    lo_lim = lo > 0 ? 0 : floor(Int64, round(lo * 1.5; digits = 1))
    hi_lim = ceil(Int64, round(hi * 1.5; digits = 1))
    return (lo_lim, hi_lim)
end

"""
    plot_histogram(s; <keyword arguments>)

Plot histogram of signal data with optional statistical overlays and comparison.

# Arguments

- `s::AbstractVector`: signal vector
- `x::Union{Nothing, Real}=nothing`: optional value to plot as vertical reference line
- `type::Symbol`: type of histogram to plot:
    - `:hist`: Standard histogram (default)
    - `:kd`: Kernel density estimate
- `bins::Int64=15`: number of histogram bins
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `draw_mean::Bool=true`: if `true`, draw mean value line
- `draw_median::Bool=true`: if `true`, draw median value line
- `mono::Bool=false`: if `true`, use a monochrome palette

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_histogram(
    s::AbstractVector,
    x::Union{Nothing, Real} = nothing;
    type::Symbol = :hist,
    bins::Int64 = 15,
    xlabel::String = "",
    ylabel::String = "",
    title::String = "",
    draw_mean::Bool = true,
    draw_median::Bool = true,
    mono::Bool = false,
)::GLMakie.Figure
    # validate
    _check_var(type, [:hist, :kd], "type")

    type === :kd && (type = :density)

    # set color palette
    pal = mono ? :grays : :darktest

    # prepare plot
    GLMakie.activate!(; title = "plot_histogram()")
    fig = GLMakie.Figure(; size = (800, 500))

    # create axis with customizable properties
    ax = GLMakie.Axis(
        fig[1, 1];
        xlabel             = xlabel,
        ylabel             = ylabel,
        title              = title,
        xminorticksvisible = true,
        xminorticks        = IntervalsBetween(10),
        xautolimitmargin   = (0.1, 0.1),
        yautolimitmargin   = (0.0, 0.1),
        _AXIS_LOCK_KWARGS...,
    )
    _style_axis!(ax)

    # plot histogram
    GLMakie.hist!(
        s;
        bins        = bins,
        colormap    = pal,
        strokecolor = :black,
        color       = :grey,
        alpha       = 0.5,
    )

    # plot vertical line at mean
    draw_mean && GLMakie.vlines!(
        mean(s);
        linestyle = :dot,
        color     = :black,
        label     = "mean",
    )
    # plot vertical line at median
    draw_median && GLMakie.vlines!(
        median(s);
        linestyle = :dash,
        color     = :black,
        label     = "median",
    )

    if !isnothing(x)
        GLMakie.vlines!(
            [x];
            linewidth = 2,
            linestyle = :dash,
            color     = mono ? :black : :red,
            label     = "test value",
        )
        prop = round(cmp_stat(s, x); digits = 3)
        _info("Proportion of values > $x: $prop")
        _info("Proportion of values < $x: $(1 - prop)")
    end

    # plot legend
    (draw_median || draw_mean || !isnothing(x)) && axislegend(; position = :rt)

    return fig
end

"""
    plot_bar(s; <keyword arguments>)

Create a bar plot from signal data with customizable labels and styling.

# Arguments

- `s::AbstractVector`: signal vector
- `glabels::Vector{String}`: group labels for x-axis ticks
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: if `true`, use a monochrome palette

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_bar(
    s::AbstractVector;
    glabels::Vector{String},
    xlabel::String = "",
    ylabel::String = "",
    title::String = "",
    mono::Bool = false,
)::GLMakie.Figure
    # validate
    length(s) == length(glabels) || throw(
        ArgumentError(
            "Lengths of signal ($(length(s))) and glabels ($(length(glabels))) must be equal.",
        ),
    )

    # set color palette
    pal   = mono ? :grays : :darktest
    color = mono ? :lightgrey : :lightblue

    # set y-axis limits
    yl = _ylims_padded(s)

    # prepare plot
    GLMakie.activate!(; title = "plot_bar()")
    fig = GLMakie.Figure(; size = (800, 500))

    # create axis with customizable properties
    ax = GLMakie.Axis(
        fig[1, 1];
        xlabel           = xlabel,
        ylabel           = ylabel,
        title            = title,
        xticks           = (eachindex(glabels), glabels),
        xautolimitmargin = (0.01, 0.01),
        yautolimitmargin = (0, 0),
        _AXIS_LOCK_KWARGS...,
    )
    GLMakie.ylims!(ax, yl)
    _style_axis!(ax)

    GLMakie.barplot!(s; color = color, colormap = pal)

    return fig
end

"""
    plot_line(s; <keyword arguments>)

Create a line plot from signal data with customizable labels and styling.

# Arguments

- `s::AbstractVector`: signal vector containing y-values to plot
- `glabels::Vector{String}`: group labels for x-axis ticks
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_line(
    s::AbstractVector;
    glabels::Vector{String},
    xlabel::String = "",
    ylabel::String = "",
    title::String = "",
)::GLMakie.Figure
    # validate
    length(s) == length(glabels) || throw(
        ArgumentError(
            "Lengths of signal ($(length(s))) and glabels ($(length(glabels))) must be equal.",
        ),
    )

    # set y-axis limits
    yl = _ylims_padded(s)

    # prepare plot
    GLMakie.activate!(; title = "plot_line()")
    fig = GLMakie.Figure(; size = (800, 500))

    # create axis with customizable properties
    ax = GLMakie.Axis(
        fig[1, 1];
        xlabel           = xlabel,
        ylabel           = ylabel,
        title            = title,
        xticks           = (eachindex(glabels), glabels),
        xautolimitmargin = (0.1, 0.1),
        yautolimitmargin = (0.1, 0.1),
        _AXIS_LOCK_KWARGS...,
    )
    GLMakie.ylims!(ax, yl)
    _style_axis!(ax)

    GLMakie.lines!(eachindex(glabels), s; color = :black)

    return fig
end

"""
    plot_line(s; <keyword arguments>)

Create a multi-line plot from multiple signal vectors (matrix rows) with customizable labels and styling.

# Arguments

- `s::AbstractMatrix`: matrix where each row represents a signal to plot
- `rlabels::Vector{String}`: signal row labels for legend
- `glabels::Vector{String}`: group labels for x-axis ticks
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: if `true`, use a monochrome palette

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_line(
    s::AbstractMatrix;
    rlabels::Vector{String},
    glabels::Vector{String},
    xlabel::String = "",
    ylabel::String = "",
    title::String = "",
    mono::Bool = false,
)::GLMakie.Figure
    # validate
    size(s, 1) == length(rlabels) || throw(
        ArgumentError(
            "Number of s rows ($(size(s, 1))) and length of rlabels ($(length(rlabels))) must be equal.",
        ),
    )
    size(s, 2) == length(glabels) || throw(
        ArgumentError(
            "Number of s columns ($(size(s, 2))) and length of glabels ($(length(glabels))) must be equal.",
        ),
    )

    # set color palette
    pal = mono ? :grays : :darktest

    # set y-axis limits
    yl = _ylims_padded(s)

    # prepare plot
    GLMakie.activate!(; title = "plot_line()")
    fig = GLMakie.Figure(; size = (800, 500))

    # create axis with customizable properties
    ax = GLMakie.Axis(
        fig[1, 1];
        xlabel           = xlabel,
        ylabel           = ylabel,
        title            = title,
        xticks           = (eachindex(glabels), glabels),
        xautolimitmargin = (0.1, 0.1),
        yautolimitmargin = (0.1, 0.1),
        _AXIS_LOCK_KWARGS...,
    )
    GLMakie.ylims!(ax, yl)
    _style_axis!(ax)

    cmap = GLMakie.resample_cmap(pal, max(2, size(s, 1)))
    for idx in axes(s, 1)
        GLMakie.lines!(
            eachindex(glabels),
            s[idx, :];
            label      = rlabels[idx],
            color      = cmap[idx],
            colormap   = pal,
            colorrange = eachindex(glabels),
        )
    end

    # plot legend
    axislegend(; position = :rt)

    return fig
end

"""
    plot_box(s; <keyword arguments>)

Create a box plot from matrix data with customizable labels and styling.

# Arguments

- `s::AbstractMatrix`: matrix where each row represents a group for which to plot a box plot
- `glabels::Vector{String}`: group labels for x-axis ticks
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: if `true`, use a monochrome palette

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_box(
    s::AbstractMatrix;
    glabels::Vector{String},
    xlabel::String = "",
    ylabel::String = "",
    title::String = "",
    mono::Bool = false,
)::GLMakie.Figure
    # validate
    size(s, 1) == length(glabels) || throw(
        ArgumentError(
            "Number of signal rows ($(size(s, 1))) and length of glabels ($(length(glabels))) must be equal.",
        ),
    )

    # set color palette
    pal   = mono ? :grays : :darktest
    color = mono ? :lightgrey : :lightblue

    # set y-axis limits
    yl = _ylims_padded(s)

    # prepare plot
    GLMakie.activate!(; title = "plot_box()")
    fig = GLMakie.Figure(; size = (800, 500))

    # create axis with customizable properties
    ax = GLMakie.Axis(
        fig[1, 1];
        xlabel           = xlabel,
        ylabel           = ylabel,
        title            = title,
        xticks           = (eachindex(glabels), glabels),
        xautolimitmargin = (0.01, 0.01),
        yautolimitmargin = (0, 0),
        _AXIS_LOCK_KWARGS...,
    )
    GLMakie.ylims!(ax, yl)
    _style_axis!(ax)

    GLMakie.boxplot!(
        repeat(eachindex(glabels), size(s, 2)),
        s[:];
        color    = color,
        colormap = pal,
    )

    return fig
end

"""
    plot_violin(s; <keyword arguments>)

Create a violin plot from matrix data with customizable labels and styling.

# Arguments

- `s::AbstractMatrix`: matrix where each row represents a group for which to create a violin plot
- `glabels::Vector{String}`: group labels for x-axis ticks
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: if `true`, use a monochrome palette

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_violin(
    s::AbstractMatrix;
    glabels::Vector{String},
    xlabel::String = "",
    ylabel::String = "",
    title::String = "",
    mono::Bool = false,
)::GLMakie.Figure
    # validate
    size(s, 1) == length(glabels) || throw(
        ArgumentError(
            "Number of s columns ($(size(s, 1))) and length of glabels ($(length(glabels))) must be equal.",
        ),
    )

    # set color palette
    pal   = mono ? :grays : :darktest
    color = mono ? :lightgrey : :lightblue

    # set y-axis limits
    yl = _ylims_padded(s)

    # prepare plot
    GLMakie.activate!(; title = "plot_violin()")
    fig = GLMakie.Figure(; size = (800, 500))

    # create axis with customizable properties
    ax = GLMakie.Axis(
        fig[1, 1];
        xlabel           = xlabel,
        ylabel           = ylabel,
        title            = title,
        xticks           = (eachindex(glabels), glabels),
        xautolimitmargin = (0.01, 0.01),
        yautolimitmargin = (0, 0),
        _AXIS_LOCK_KWARGS...,
    )
    GLMakie.ylims!(ax, yl)
    _style_axis!(ax)

    GLMakie.violin!(
        repeat(eachindex(glabels), size(s, 2)),
        s[:];
        strokecolor = :black,
        strokewidth = 0.25,
        color       = color,
    )

    return fig
end

"""
    plot_dots(s; <keyword arguments>)

Create a dots plot from matrix data with customizable labels and styling.

# Arguments

- `s::AbstractArray`: matrix where each row represents a group for which to create a violin plot
- `glabels::Vector{String}`: group labels for x-axis ticks
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: if `true`, use a monochrome palette

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_dots(
    s::AbstractArray;
    glabels::Vector{String},
    xlabel::String = "",
    ylabel::String = "",
    title::String = "",
    mono::Bool = false,
)::GLMakie.Figure
    # validate
    size(s, 1) == length(glabels) || throw(
        ArgumentError(
            "Number of signal rows ($(size(s, 1))) and length of glabels ($(length(glabels))) must be equal.",
        ),
    )

    # set color palette
    pal = mono ? :grays : :darktest

    # set y-axis limits
    yl = _ylims_padded(s)

    # prepare plot
    GLMakie.activate!(; title = "plot_dots()")
    fig = GLMakie.Figure(; size = (800, 500))

    # create axis with customizable properties
    ax = GLMakie.Axis(
        fig[1, 1];
        xlabel           = xlabel,
        ylabel           = ylabel,
        title            = title,
        xticks           = (eachindex(glabels), glabels),
        xautolimitmargin = (0.25, 0.25),
        yautolimitmargin = (0, 0),
        _AXIS_LOCK_KWARGS...,
    )
    GLMakie.ylims!(ax, yl)
    _style_axis!(ax)

    cmap = GLMakie.resample_cmap(pal, max(2, length(glabels)))
    for idx in eachindex(glabels)
        if mono
            GLMakie.scatter!(repeat([idx], size(s, 2)), s[idx, :]; color = :black)
        else
            GLMakie.scatter!(
                repeat([idx], size(s, 2)), s[idx, :];
                color      = cmap[idx],
                colormap   = pal,
                colorrange = eachindex(glabels),
            )
        end
    end

    return fig
end

"""
    plot_paired(signal; <keyword arguments>)

Create a paired data plot showing connections between paired observations (grouped by rows).

# Arguments

- `s::AbstractArray`: paired data to plot
- `glabels::Vector{String}`: group labels for x-axis ticks
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: if `true`, use a monochrome palette

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_paired(
    s::AbstractArray;
    glabels::Vector{String},
    xlabel::String = "",
    ylabel::String = "",
    title::String = "",
    mono::Bool = false,
)::GLMakie.Figure
    # validate
    size(s, 1) == length(glabels) || throw(
        ArgumentError(
            "Number of signal rows ($(size(s, 1))) and length of glabels ($(length(glabels))) must be equal.",
        ),
    )

    # set color palette
    pal = mono ? :grays : :darktest

    # set y-axis limits
    yl = _ylims_padded(s)

    # prepare plot
    GLMakie.activate!(; title = "plot_paired()")
    fig = GLMakie.Figure(; size = (800, 500))

    # create axis with customizable properties
    ax = GLMakie.Axis(
        fig[1, 1];
        xlabel           = xlabel,
        ylabel           = ylabel,
        title            = title,
        xticks           = (eachindex(glabels), glabels),
        xautolimitmargin = (0.25, 0.25),
        yautolimitmargin = (0, 0),
        xzoomlock        = true,
        yzoomlock        = true,
        xpanlock         = true,
        ypanlock         = true,
        xrectzoom        = false,
        yrectzoom        = false,
    )
    GLMakie.ylims!(ax, yl)
    _style_axis!(ax)

    cmap = GLMakie.resample_cmap(pal, max(2, length(glabels)))

    for idx in eachindex(glabels)
        if mono
            GLMakie.scatter!(repeat([idx], size(s, 2)), s[idx, :]; color = :black)
        else
            GLMakie.scatter!(
                repeat([idx], size(s, 2)), s[idx, :];
                color      = cmap[idx],
                colormap   = pal,
                colorrange = eachindex(glabels),
            )
        end
    end

    for idx in axes(s, 2)
        GLMakie.lines!(eachindex(glabels), s[:, idx]; color = :black, linewidth = 0.5)
    end

    return fig
end

"""
    plot_ci(s, s_ci_l, s_ci_h, t; <keyword arguments>)

Plot a signal with its confidence interval (shaded region) over time.

# Arguments

- `s::AbstractVector`: signal vector containing the mean values
- `s_l::AbstractVector`: lower bound of the confidence interval
- `s_u::AbstractVector`: upper bound of the confidence interval
- `t::AbstractVector`: time points corresponding to the signal values
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: if `true`, use a monochrome palette

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_ci(
    s::AbstractVector,
    s_l::AbstractVector,
    s_u::AbstractVector,
    t::AbstractVector;
    xlabel::String = "",
    ylabel::String = "",
    title::String = "",
    mono::Bool = false,
)::GLMakie.Figure
    # validate
    length(s) == length(s_l) == length(s_u) || throw(
        ArgumentError("All input signals must be of the same length."),
    )

    # set y-axis limits
    yl     = (floor(minimum(s_l); digits = 0), ceil(maximum(s_u); digits = 0))
    yl     = _tuple_max(yl)
    yticks = [yl[1], 0, yl[2]]

    # prepare plot
    GLMakie.activate!(; title = "plot_ci()")
    fig = GLMakie.Figure(; size = (800, 500))

    # create axis with customizable properties
    ax = GLMakie.Axis(
        fig[1, 1];
        xlabel             = xlabel,
        ylabel             = ylabel,
        title              = title,
        xticks             = LinearTicks(10),
        yticks             = yticks,
        xminorticksvisible = true,
        xminorticks        = IntervalsBetween(10),
        _AXIS_LOCK_KWARGS...,
    )
    GLMakie.ylims!(ax, yl)
    _style_axis!(ax)

    GLMakie.band!(t, s_u, s_l; alpha = 0.25, color = :grey, strokewidth = 0.5)
    GLMakie.lines!(t, s; color = :black, linewidth = 2)

    return fig
end
