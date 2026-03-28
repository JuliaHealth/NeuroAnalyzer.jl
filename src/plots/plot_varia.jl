export plot_matrix
export plot_xac
export plot_histogram
export plot_bar
export plot_line
export plot_box
export plot_violin
export plot_dots
export plot_paired
export plot_polar
export plot_eros
export plot_erop
export plot_icatopo
export plot_ci
export plot_heatmap
export plot_imf
export plot_fi
export plot_phase
export plot_polezero
export plot_dwc

# ---------------------------------------------------------------------------
# shared helpers
# ---------------------------------------------------------------------------

# Keyword arguments applied to every locked/non-interactive Axis.
const _AXIS_LOCK_KWARGS = (
    xzoomlock  = true,
    yzoomlock  = true,
    xpanlock   = true,
    ypanlock   = true,
    xrectzoom  = false,
    yrectzoom  = false,
)

# Apply standard font sizes to an Axis.
function _style_axis!(ax)
    ax.titlesize      = 18
    ax.xlabelsize     = 18
    ax.ylabelsize     = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12
    return ax
end

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
    plot_matrix(m; <keyword arguments>)

Plot matrix.

# Arguments

- `m::Matrix{<:Real}`: matrix to plot
- `xlabels::Vector{String}`: labels for x-axis ticks
- `ylabels::Vector{String}`: labels for y-axis ticks
- `xlabel::String=""`: x-axis labels
- `ylabel::String=""`: y-axis labels
- `title::String=""`: plot title
- `cb::Bool=true`: if `true`, show colorbar
- `cb_title::String=""`: colorbar title
- `xrot::Int64=90`: x-axis label rotation in degrees
- `mono::Bool=false`: if `true`, use a monochrome palette

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_matrix(
    m::AbstractMatrix;
    xlabels::Vector{String},
    ylabels::Vector{String},
    xlabel::String = "",
    ylabel::String = "",
    title::String = "",
    cb::Bool = true,
    cb_title::String = "",
    xrot::Int64 = 90,
    mono::Bool = false,
)::GLMakie.Figure
    # validate
    size(m, 1) == size(m, 2) || throw(ArgumentError("Matrix must be square."))
    length(xlabels) == length(ylabels) || throw(
        ArgumentError(
            "Lengths of xlabels ($(length(xlabels))) and ylabels ($(length(ylabels))) must be equal.",
        ),
    )
    length(xlabels) == size(m, 1) || throw(
        ArgumentError(
            "Length of xlabels ($(length(xlabels))) and matrix size $(size(m)) must be equal.",
        ),
    )
    length(ylabels) == size(m, 2) || throw(
        ArgumentError(
            "Length of ylabels ($(length(ylabels))) and matrix size $(size(m)) must be equal.",
        ),
    )


    n = size(m, 1)

    # set color palette
    pal = mono ? :grays : :bluesreds


    # prepare plot
    GLMakie.activate!(; title = "plot_matrix()")
    fig = GLMakie.Figure(; size = (800, 800))

    # create axis with customizable properties
    ax = GLMakie.Axis(
        fig[1, 1];
        xlabel               = xlabel,
        ylabel               = ylabel,
        title                = title,
        xticks               = (1:n, xlabels),
        xticklabelrotation   = deg2rad(xrot),
        xticksvisible        = false,
        yticks               = (1:n, ylabels),
        yticksvisible        = false,
        xautolimitmargin     = (0, 0),
        yautolimitmargin     = (0, 0),
        _AXIS_LOCK_KWARGS...,
    )
    _style_axis!(ax)


    hm = GLMakie.heatmap!(m'; colormap = pal)
    cb && GLMakie.Colorbar(fig[1, 2], hm; label = cb_title, labelsize = 16)


    return fig
end

"""
    plot_xac(m, lags; <keyword arguments>)

Plot cross-correlation, auto-correlation, or covariance between signal(s) and lags.

# Arguments

- `m::AbstractVector`: cross/auto-covariance/correlation matrix of
- `lags::AbstractVector`: vector of lag values in seconds
- `xlabel::String="Lag [s]"`: x-axis labels
- `ylabel::String=""`: y-axis labels
- `title::String=""`: plot title

# Returns

- `GLMakie.Figure`: the plotted figure

# Notes

- Positive lags represent future correlation, negative lags represent past correlation
"""
function plot_xac(
    m::AbstractVector,
    lags::AbstractVector;
    xlabel::String = "Lag [s]",
    ylabel::String = "",
    title::String = "",
)::GLMakie.Figure

    # prepare plot
    GLMakie.activate!(; title = "plot_xac()")
    fig = GLMakie.Figure(; size = (800, 300))

    # create axis with customizable properties
    ax  = GLMakie.Axis(
        fig[1, 1];
        xlabel               = xlabel,
        ylabel               = ylabel,
        title                = title,
        xminorticksvisible   = true,
        xminorticks          = IntervalsBetween(10),
        xautolimitmargin     = (0, 0),
        yautolimitmargin     = (0, 0),
        _AXIS_LOCK_KWARGS...,
    )
    _style_axis!(ax)


    GLMakie.lines!(lags, m; linewidth = 1, color = :black)


    return fig
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

    xticks = if !isnothing(x)
        [
            round(minimum(s); digits = 2),
            round(mean(s);    digits = 2),
            round(median(s);  digits = 2),
            round(x;          digits = 2),
            round(maximum(s); digits = 2),
        ]
    else
        [
            round(minimum(s); digits = 2),
            round(mean(s);    digits = 2),
            round(median(s);  digits = 2),
            round(maximum(s); digits = 2),
        ]
    end


    !draw_median && deleteat!(xticks, 3)
    !draw_mean   && deleteat!(xticks, 2)
    sort!(unique(xticks))


    # prepare plot
    GLMakie.activate!(; title = "plot_histogram()")
    fig = GLMakie.Figure(; size = (800, 500))

    # create axis with customizable properties
    ax  = GLMakie.Axis(
        fig[1, 1];
        xlabel             = xlabel,
        ylabel             = ylabel,
        title              = title,
        xticks             = xticks,
        xticklabelrotation = pi / 2,
        xautolimitmargin     = (0, 0),
        yautolimitmargin     = (0, 0),
        _AXIS_LOCK_KWARGS...,
    )
    GLMakie.xlims!(ax, extrema(xticks))
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
        round(mean(s); digits = 2);
        linestyle = :dot,
        color     = :black,
        label     = "mean",
    )
    # plot vertical line at median
    draw_median && GLMakie.vlines!(
        round(median(s); digits = 2);
        linestyle = :dash,
        color     = :grey,
        label     = "median",
    )


    if !isnothing(x)
        GLMakie.vlines!(
            [x];
            linewidth = 2,
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
    yl    = _ylims_padded(s)

    # prepare plot
    GLMakie.activate!(; title = "plot_bar()")
    fig = GLMakie.Figure(; size = (800, 500))

    # create axis with customizable properties
    ax = GLMakie.Axis(
        fig[1, 1];
        xlabel             = xlabel,
        ylabel             = ylabel,
        title              = title,
        xticks             = (eachindex(glabels), glabels),
        xautolimitmargin   = (0.01, 0.01),
        yautolimitmargin   = (0, 0),
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
    yl  = _ylims_padded(s)


    # prepare plot
    GLMakie.activate!(; title = "plot_line()")
    fig = GLMakie.Figure(; size = (800, 500))

    # create axis with customizable properties
    ax  = GLMakie.Axis(
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


    cmap = GLMakie.resample_cmap(pal, size(s, 1))
    for idx in axes(s, 1)
        GLMakie.lines!(
            eachindex(glabels),
            s[idx, :];
            label       = rlabels[idx],
            color       = cmap[idx],
            colormap    = pal,
            colorrange  = eachindex(glabels),
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
    yl    = _ylims_padded(s)


    # prepare plot
    GLMakie.activate!(; title = "plot_box()")
    fig = GLMakie.Figure(; size = (800, 500))

    # create axis with customizable properties
    ax  = GLMakie.Axis(
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
    yl    = _ylims_padded(s)


    # prepare plot
    GLMakie.activate!(; title = "plot_violin()")
    fig = GLMakie.Figure(; size = (800, 500))

    # create axis with customizable properties
    ax  = GLMakie.Axis(
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
    yl  = _ylims_padded(s)

    # prepare plot
    GLMakie.activate!(; title = "plot_dots()")
    fig = GLMakie.Figure(; size = (800, 500))

    # create axis with customizable properties
    ax  = GLMakie.Axis(
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


    cmap = GLMakie.resample_cmap(pal, length(glabels))
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
    yl  = _ylims_padded(s)

    # prepare plot
    GLMakie.activate!(; title = "plot_paired()")
    fig = GLMakie.Figure(; size = (800, 500))

    # create axis with customizable properties
    ax  = GLMakie.Axis(
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


    cmap = GLMakie.resample_cmap(pal, length(glabels))


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
    plot_polar(s; <keyword arguments>)

Create a polar plot from signal data with customizable styling.

# Arguments

- `s::Union{AbstractVector, AbstractMatrix}`: input data to plot; vector must contain 2 values: phases and lengths; matrix must contain 2 columns: phases and lengths
- `m::Tuple{Real, Real}=(0, 0)`: major value to plot
- `title::String=""`: plot title
- `mono::Bool=false`: if `true`, use a monochrome palette
- `ticks::Bool=false`: draw x- and y-axis ticks

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_polar(
    s::Union{AbstractVector, AbstractMatrix};
    m::Tuple{Real, Real} = (0, 0),
    title::String = "",
    mono::Bool = false,
    ticks::Bool = true,
)::GLMakie.Figure
    # validate
    size(s, 1) == 2 && (s = s')
    length(m) == 2 ||
        throw(ArgumentError("m must contain 2 values: phases and lengths."))
    ndims(s) > 1 && size(s, 2) != 2 &&
        throw(ArgumentError("signal must contain 2 columns: phases and lengths."))


    # prepare plot
    GLMakie.activate!(; title = "plot_polar()")
    fig = GLMakie.Figure(; size = (800, 800))

    # create axis with customizable properties
    ax  = GLMakie.PolarAxis(
        fig[1, 1];
        title         = title,
        thetazoomlock = true,
        rzoomlock     = true,
    )
    !ticks && hidespines!(ax)


    if ndims(s) == 1
        for idx in eachindex(s)
            GLMakie.lines!([0, s[idx]], [0, 1]; linewidth = 2, color = :black)
        end
    else
        for idx in axes(s, 1)
            GLMakie.lines!([0, s[idx, 1]], [0, s[idx, 2]]; linewidth = 2, color = :black)
        end
    end


    if m != (0, 0)
        GLMakie.lines!(
            [0, m[1]], [0, m[2]];
            linewidth = 2,
            color     = mono ? :darkgray : :red,
        )
    end


    return fig
end

"""
    plot_eros(s, f, t; <keyword arguments>)

Plot an Event-Related Oscillations (ERO) spectrogram with customizable visualization options.

# Arguments

- `sp::AbstractArray`: ERO spectrogram power values, shape (frequency, time, powers)
- `sf::AbstractVector`: vector of frequency values in Hz
- `st::AbstractVector`: vector of time values in seconds
- `db::Bool=true`: whether to display power values in decibels
- `frq::Symbol=:lin`: frequency scaling (`:lin` for linear, `:log` for logarithmic)
- `flim::Tuple{Real, Real}=(f[1], f[end])`: frequency limits for the plot
- `tm::Union{Nothing, Int64, Vector{Int64}} = nothing`: time markers (in milliseconds) to be plot as vertical lines, useful for adding topoplots at these time points
- `xlabel::String="default"`: x-axis label
- `ylabel::String="default"`: y-axis label
- `title::String="default"`: plot title
- `cb::Bool=true`: if `true`, show colorbar
- `mono::Bool=false`: if `true`, use a monochrome palette
- `units::String="μV"`: power units
- `smooth::Bool=false`: if `true`, apply Gaussian blur smoothing
- `ks::Int64=3`: smoothing kernel size; larger kernel means more smoothing

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_eros(
    sp::AbstractArray,
    sf::AbstractVector,
    st::AbstractVector;
    db::Bool = true,
    frq::Symbol = :lin,
    flim::Tuple{Real, Real} = (sf[1], sf[end]),
    tm::Union{Nothing, Int64, Vector{Int64}} = nothing,
    xlabel::String = "default",
    ylabel::String = "default",
    title::String = "default",
    cb::Bool = true,
    mono::Bool = false,
    units::String = "μV",
    smooth::Bool = false,
    ks::Int64 = 3,
)::GLMakie.Figure
    # validate
    size(sp, 1) == length(sf) || throw(
        ArgumentError(
            "Length of sf ($(length(sf))) and number of spectrogram rows ($(size(sp, 1))) must be equal.",
        ),
    )
    size(sp, 2) == length(st) || throw(
        ArgumentError(
            "Length of st ($(length(st))) and number of spectrogram columns ($(size(sp, 2))) must be equal.",
        ),
    )
    ndims(sp) == 3 || throw(ArgumentError("sp must have 3 dimensions."))
    size(sp, 3) <= 2 || throw(ArgumentError("sp must contain ≤ 2 epochs."))
    ks > 0 || throw(ArgumentError("ks must be ≥ 1."))
     _check_var(frq, [:lin, :log], "frq")
    _check_tuple(flim, extrema(sf), "flim")

    # set color palette
    pal      = mono ? :grays : :darktest

    # colorbar title
    cb_title = db ? "[dB $units^2/Hz]" : "[$units^2/Hz]"

    # set frequency limits
    if frq === :lin
        yt = flim[2] > 100 ? (flim[1]:10:flim[2]) : (flim[1]:5:flim[2])
    else
        if flim[1] == 0
            _warn("Lower frequency bound truncated to $(sf[2]) Hz")
            flim = (sf[2], flim[2])
        end
        yt = round.(logspace(flim[1], flim[2], nfrq); digits = 1)
    end


    # apply Gaussian filter if requested
    if smooth
        for idx in axes(sp, 3)
            sp[:, :, idx] = imfilter(@view(sp[:, :, idx]), Kernel.gaussian(ks))
        end
    end


    # resolve time markers to indices into a new local variable — do not mutate the input
    tm_indices = if !isnothing(tm)
        markers = tm isa Int64 ? [tm] : tm
        for val in markers
            val / 1000 >= st[1] || throw(
                ArgumentError(
                    "tm value ($val) is out of epoch time segment ($(st[1]):$(st[end])).",
                ),
            )
            val / 1000 <= st[end] || throw(
                ArgumentError(
                    "tm value ($val) is out of epoch time segment ($(st[1]):$(st[end])).",
                ),
            )
        end
        [vsearch(val / 1000, st) for val in markers]
    else
        Int64[]
    end


    function _draw_axis(fig, pos, xl, yl, tt, nticks)
        ax = GLMakie.Axis(
            fig[pos...];
            xlabel             = xl,
            ylabel             = yl,
            title              = tt,
            xticks             = LinearTicks(nticks),
            xminorticksvisible = true,
            xminorticks        = IntervalsBetween(10),
            yticks             = yt,
            yscale             = frq === :lin ? identity : log,
            _AXIS_LOCK_KWARGS...,
        )
        GLMakie.ylims!(ax, flim)
        _style_axis!(ax)
        return ax
    end


    function _draw_markers!(fig, pos, indices)
        for i in indices
            GLMakie.vlines!(fig[pos...], [st[i]]; color = :black, linewidth = 1)
        end
    end


    if size(sp, 3) == 1

        # set default values
        xl, yl, tt = _set_defaults(
            xlabel, ylabel, title,
            "Time [ms]", "Frequency [Hz]", "Averaged spectrograms of epochs",
        )


        # prepare plot
        GLMakie.activate!(; title = "plot_eros()")
        fig = GLMakie.Figure(; size = (900, 450))

        # create axis with customizable properties
        ax  = _draw_axis(fig, (1, 1), xl, yl, tt, 15)
        hm  = GLMakie.heatmap!(ax, st, sf, sp[:, :, 1]'; colormap = pal)
        cb && GLMakie.Colorbar(fig[1, 2], hm; label = cb_title, labelsize = 16)
        _draw_markers!(fig, (1, 1), tm_indices)


    else

        # set default values
        xl1, yl1, tt1 = _set_defaults(
            xlabel, ylabel, title,
            "Time [ms]", "Frequency [Hz]", "ERP spectrogram",
        )
        xl2, yl2, tt2 = _set_defaults(
            xlabel, ylabel, title,
            "Time [ms]", "Frequency [Hz]", "Averaged spectrograms of ERP epochs",
        )


        # prepare plot
        GLMakie.activate!(; title = "plot_eros()")
        fig = GLMakie.Figure(; size = (1200, 800))


        # create axis with customizable properties
        ax1 = _draw_axis(fig, (1, 1), xl1, yl1, tt1, 10)
        hm1 = GLMakie.heatmap!(ax1, st, sf, sp[:, :, 1]'; colormap = pal)
        cb && GLMakie.Colorbar(fig[1, 2], hm1; label = cb_title, labelsize = 16)
        GLMakie.vlines!(ax1, [0]; linestyle = :dash, linewidth = 0.5, color = :black)
        _draw_markers!(fig, (1, 1), tm_indices)


        # create axis with customizable properties
        ax2 = _draw_axis(fig, (2, 1), xl2, yl2, tt2, 10)
        hm2 = GLMakie.heatmap!(ax2, st, sf, sp[:, :, 2]'; colormap = pal)
        cb && GLMakie.Colorbar(fig[2, 2], hm2; label = cb_title, labelsize = 16)
        GLMakie.vlines!(ax2, [0]; linestyle = :dash, linewidth = 0.5, color = :black)
        _draw_markers!(fig, (2, 1), tm_indices)
    end


    return fig
end

"""
    plot_erop(fig, f; <keyword arguments>)

Plot the power spectrum of Event-Related Oscillations (ERO) with customizable visualization.

# Arguments

- `sp::AbstractArray`: ERO power values, shape (frequency, powers)
- `sf::AbstractVector`: vector of frequency values in Hz
- `db::Bool=true`: whether to display power values in decibels
- `xlabel::String="default"`: x-axis label
- `ylabel::String="default"`: y-axis label
- `title::String="default"`: plot title
- `frq::Symbol=:lin`: frequency scaling (`:lin` for linear, `:log` for logarithmic)
- `flim::Tuple{Real, Real}=(f[1], f[end])`: frequency limits for the plot
- `cb::Bool=true`: if `true`, show colorbar
- `units::String="μV"`: power units
- `mono::Bool=false`: if `true`, use a monochrome palette

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_erop(
    sp::AbstractArray,
    sf::AbstractVector;
    db::Bool = true,
    xlabel::String = "default",
    ylabel::String = "default",
    title::String = "default",
    flim::Tuple{Real, Real} = (sf[1], sf[end]),
    frq::Symbol = :lin,
    units::String = "μV",
    mono::Bool = false,
)::GLMakie.Figure
    # validate
    _in(flim[2], (sf[1], sf[end]), "flim")
    size(sp, 1) == length(sf) || throw(
        ArgumentError(
            "Length of sf ($(length(sf))) and number of power rows ($(size(sp, 1))) must be equal.",
        ),
    )
    ndims(sp) == 2 || throw(ArgumentError("sp must have 2 dimensions."))
    size(sp, 2) <= 2 || throw(ArgumentError("sp must contain ≤ 2 epochs."))
    _check_var(frq, [:lin, :log], "frq")


    # set frequency limits
    if frq === :log && flim[1] == 0
        _warn("Lower frequency bound truncated to $(sf[2]) Hz")
        flim = (sf[2], flim[2])
    end


    # set y-axis label
    power_ylabel = db ? "Power [dB $units^2/Hz]" : "Power [$units^2/Hz]"


    function _make_power_axis(fig, pos, xl, yl, tt)
        ax = GLMakie.Axis(
            fig[pos...];
            xlabel             = xl,
            ylabel             = yl,
            title              = tt,
            xminorticksvisible = true,
            xminorticks        = IntervalsBetween(10),
            xscale             = frq === :lin ? identity : log,
            _AXIS_LOCK_KWARGS...,
        )
        GLMakie.xlims!(ax, flim)
        _style_axis!(ax)
        return ax
    end


    # prepare plot
    GLMakie.activate!(; title = "plot_erop()")


    if size(sp, 2) == 1

        # set default values
        xl, _, tt = _set_defaults(
            xlabel, ylabel, title,
            "Frequency [Hz]", power_ylabel, "Averaged power-spectra of epochs",
        )


        # prepare plot
        fig = GLMakie.Figure(; size = (900, 450))

        # create axis with customizable properties
        ax  = _make_power_axis(fig, (1, 1), xl, power_ylabel, tt)
        GLMakie.lines!(ax, sf, sp[:, 1]; color = :black)


    else

        # set default values
        xl, _, tt1 = _set_defaults(
            xlabel, ylabel, title,
            "Frequency [Hz]", power_ylabel, "ERP power-spectrum",
        )
        _, _, tt2 = _set_defaults(
            xlabel, ylabel, title,
            "Frequency [Hz]", power_ylabel, "Averaged power-spectra of ERP epochs",
        )


        # prepare plot
        fig = GLMakie.Figure(; size = (1200, 800))

        # create axis with customizable properties
        ax1 = _make_power_axis(fig, (1, 1), xl, power_ylabel, tt1)
        GLMakie.lines!(ax1, sf, sp[:, 1]; color = :black)


        ax2 = _make_power_axis(fig, (2, 1), xl, power_ylabel, tt2)
        GLMakie.lines!(ax2, sf, sp[:, 2]; color = :black)

    end


    return fig
end

"""
    plot_icatopo(obj; <keyword arguments>)

Create a topographical plot of Independent Component Analysis (ICA) components from EEG/MEG data.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object containing channel locations and data
- `ic::Matrix{Float64}`: ICA component matrix IC(1)..IC(n) containing spatial patterns
- `ic_mw::Matrix{Float64}`: weighting matrix for ICA components
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `ic_idx::Union{Int64, Vector{Int64}, AbstractRange}=axes(ic_idx, 1)`: component indices to plot, default is all components
- `tpos::Union{Nothing, Real, AbstractVector}=nothing`: time point(s) in seconds to plot, ignored if `data` is provided
- `imethod::Symbol=:sh`: interpolation method:
    - `:sh`: Shepard
    - `:mq`: Multiquadratic
    - `:imq`: Inverse Multiquadratic
    - `:tp`: ThinPlate
    - `:nn`: Nearest Neighbour
    - `:ga`: Gaussian
- `nmethod::Symbol=:minmax`: method for normalization, see `normalize()`
- `contours::Int64=0`: number of contour levels to plot
- `electrodes::Bools=true`: if `true`, plot electrode locations over topography
- `ps::Symbol`: plot size:
    - `:l`: large (800×800 px)
    - `:m`: medium (300×300 px)
    - `:s`: small (100×100 px)

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_icatopo(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    ic::Matrix{Float64},
    ic_mw::Matrix{Float64},
    ic_idx::Union{Int64, Vector{Int64}, AbstractRange} = axes(ic, 1),
    tpos::Union{Nothing, Real, AbstractVector},
    imethod::Symbol = :sh,
    nmethod::Symbol = :minmax,
    contours::Int64 = 0,
    electrodes::Bool = true,
    ps::Symbol = :l,
)::GLMakie.Figure
    fig_topo = GLMakie.Figure[]
    for idx in eachindex(ic_idx)
        obj_tmp = ica_reconstruct(
            obj; ch = ch, ic = ic, ic_mw = ic_mw, ic_idx = idx, keep = true,
        )
        fig_tmp = plot_topo(
            obj_tmp;
            ch         = ch,
            tpos       = tpos,
            title      = "IC $idx",
            imethod    = imethod,
            nmethod    = nmethod,
            contours   = contours,
            electrodes = electrodes,
            ps         = ps,
            cb         = true,
        )
        push!(fig_topo, fig_tmp)
    end


    return plot_compose(fig_topo; layout = (1, length(ic_idx)))
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
    ax  = GLMakie.Axis(
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

"""
    plot_heatmap(m; <keyword arguments>)

Plot a heatmap with customizable labels, styling, and optional threshold highlighting.

# Arguments

- `m::AbstractMatrix`: matrix containing values to plot as heatmap
- `x::AbstractVector`: x-axis coordinates for each column of the matrix
- `y::AbstractVector`: y-axis coordinates for each row of the matrix
- `xlabel::String=""`: x-axis label
- `ylabel::String=""`: y-axis label
- `title::String=""`: plot title
- `mono::Bool=false`: if `true`, use a monochrome palette
- `cb::Bool=true`: if `true`, show colorbar
- `cb_title::String=""`: colorbar title
- `threshold::Union{Nothing, Real, Tuple{Real, Real}}=nothing`: threshold for marking regions
    - if `Real`, use a single threshold value
    - if `Tuple{Real, Real}`, use a range for `:in` or `:bin` thresholding
- `threshold_type::Symbol=:neq`: rule for threshold-based highlighting:
    - `:eq`: draw region where values are not equal to threshold
    - `:neq`: draw region where values are equal to threshold
    - `:geq`: draw region where values are ≥ threshold
    - `:leq`: draw region where values are ≤ threshold
    - `:g`: draw region where values are > threshold
    - `:l`: draw region where values are < threshold

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_heatmap(
    m::AbstractMatrix;
    x::AbstractVector,
    y::AbstractVector,
    xlabel::String = "",
    ylabel::String = "",
    title::String = "",
    mono::Bool = false,
    cb::Bool = true,
    cb_title::String = "",
    threshold::Union{Nothing, Real, Tuple{Real, Real}} = nothing,
    threshold_type::Symbol = :neq,
)::GLMakie.Figure
    # validate
    size(m, 1) == length(y) || throw(
        ArgumentError(
            "Number of m rows ($(size(m, 1))) and y length ($(length(y))) must be equal.",
        ),
    )
    size(m, 2) == length(x) || throw(
        ArgumentError(
            "Number of m columns ($(size(m, 2))) and x length ($(length(x))) must be equal.",
        ),
    )


    # set color palette
    pal = mono ? :grays : :bluesreds


    # prepare plot
    GLMakie.activate!(; title = "plot_heatmap()")
    fig = GLMakie.Figure(; size = (800, 500))

    # create axis with customizable properties
    ax  = GLMakie.Axis(
        fig[1, 1];
        xlabel = xlabel,
        ylabel = ylabel,
        title  = title,
        xticks = LinearTicks(10),
        yticks = LinearTicks(10),
        _AXIS_LOCK_KWARGS...,
    )
    _style_axis!(ax)

    hm = GLMakie.heatmap!(x, y, m'; colormap = pal)

    # plot colorbar
    cb && GLMakie.Colorbar(fig[1, 2], hm; label = cb_title, labelsize = 16)


    # apply thresholding
    if !isnothing(threshold)
        _, bm  = seg_extract(m; threshold = threshold, threshold_type = threshold_type)
        reg    = ones(size(m)) .* minimum(m)
        reg[bm] .= maximum(m)
        GLMakie.contour!(ax, x, y, reg'; levels = 1, color = :black, linewidth = 2)
    end


    return fig
end

"""
    plot_imf(imf; <keyword arguments>)

Plot intrinsic mode functions (IMF), the residual, and the reconstructed signal from Empirical Mode Decomposition (EMD).

# Arguments

- `imf::Matrix{Float64}`: matrix where each row represents an IMF component from EMD
- `n::Int64=size(imf, 1) - 1`: number of IMF components to plot
- `t::AbstractVector`: time points corresponding to the signal values

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_imf(
    imf::Matrix{Float64};
    n::Int64 = size(imf, 1) - 1,
    t::AbstractVector,
)::GLMakie.Figure
    # validate
    n > 0 || throw(ArgumentError("n must be ≥ 1."))
    n + 1 <= size(imf, 1) || throw(ArgumentError("n must be ≤ $(size(imf, 1) - 1)."))
    size(imf, 2) == length(t) || throw(
        ArgumentError(
            "Length of t ($(length(t))) and number of imf columns ($(size(imf, 2))) must be equal.",
        ),
    )


    # the last row of imf is the residual; the reconstruction is the sum of all rows
    s_restored = sum(imf; dims = 1)[:]
    imf_plot   = vcat(imf, s_restored')


    # set y-axis limits
    ylim   = (floor(minimum(imf_plot); digits = 0), ceil(maximum(imf_plot); digits = 0))
    ylim   = _tuple_max(ylim)
    yticks = [ylim[1], 0, ylim[2]]


    # prepare plot
    GLMakie.activate!(; title = "plot_imf()")
    fig = GLMakie.Figure(; size = (1200, 800))
    nr  = ceil(Int64, (n + 1) / 2)


    idx  = 1
    cidx = 1
    for idx1 in 1:nr
        cidx = 1
        for idx2 in 1:2
            if idx <= n + 1
                label = idx == n + 1 ? "Residual" : "IMF: $idx"
                # create axis with customizable properties
                ax = GLMakie.Axis(
                    fig[idx1, idx2];
                    xlabel             = "Time [s]",
                    ylabel             = "",
                    title              = label,
                    xticks             = LinearTicks(10),
                    xminorticksvisible = true,
                    xminorticks        = IntervalsBetween(10),
                    yticks             = yticks,
                    _AXIS_LOCK_KWARGS...,
                )
                GLMakie.ylims!(ax, ylim)
                _style_axis!(ax)
                GLMakie.lines!(ax, t, imf[idx, :]; color = :black)
                idx  += 1
                cidx += 1
            end
        end
    end


    # place the reconstructed signal spanning both columns
    row = cidx == 1 ? nr : nr + 1
    ax  = GLMakie.Axis(
        fig[row, 1:2];
        xlabel             = "Time [s]",
        ylabel             = "",
        title              = "Reconstructed signal",
        xticks             = LinearTicks(10),
        xminorticksvisible = true,
        xminorticks        = IntervalsBetween(10),
        yticks             = yticks,
        _AXIS_LOCK_KWARGS...,
    )
    GLMakie.ylims!(ax, ylim)
    _style_axis!(ax)
    GLMakie.lines!(ax, t, s_restored; color = :black)


    return fig
end

"""
    plot_fi(fi, st; <keyword arguments>)

Plot instantaneous frequency over time from time-frequency analysis (e.g., Hilbert-Huang Transform).

# Arguments

- `fi::Vector{Float64}`: vector of instantaneous frequency values in Hz
- `st::Vector{Float64}`: vector of time points in seconds corresponding to frequency values
- `xlabel::String="default"`: x-axis label, default is Time [s]
- `ylabel::String="default"`: y-axis label, default is Power [μV^2/Hz]
- `title::String="default"`: plot title

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_fi(
    fi::Vector{Float64},
    st::Vector{Float64};
    xlabel::String = "default",
    ylabel::String = "default",
    title::String = "default",
)::GLMakie.Figure
    # validate
    length(fi) == length(st) || throw(
        ArgumentError(
            "Length of frequencies ($(length(fi))) and time points ($(length(st))) must be equal.",
        ),
    )

    # set default values
    xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "Frequency [Hz]", "")


    # prepare plot
    GLMakie.activate!(; title = "plot_fi()")
    fig = GLMakie.Figure(; size = (900, 450))

    # create axis with customizable properties
    ax  = GLMakie.Axis(
        fig[1, 1];
        xlabel             = xl,
        ylabel             = yl,
        title              = tt,
        xticks             = LinearTicks(10),
        xminorticksvisible = true,
        xminorticks        = IntervalsBetween(10),
        _AXIS_LOCK_KWARGS...,
    )
    GLMakie.xlims!(ax, _xlims(st))
    _style_axis!(ax)


    GLMakie.lines!(st, fi; linewidth = 1, color = :black)


    return fig
end

"""
    plot_phase(ph, sf; <keyword arguments>)

Plot phase values from time-frequency analysis with customizable styling and units.

# Arguments

- `ph::Vector{Float64}`: vector of phase values (radians or degrees)
- `sf::Vector{Float64}`: vector of corresponding frequencies or time points
- `unit::Symbol=:rad`: phase unit specification (`:rad` radians or `:deg` degrees)
- `type::Symbol=:line`: plot type (`:line`: line plot, `:stem`: stem plot with markers)
- `xlabel::String="default"`: x-axis label
- `ylabel::String="default"`: y-axis label
- `title::String="default"`: plot title

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_phase(
    ph::Vector{Float64},
    sf::Vector{Float64};
    unit::Symbol = :rad,
    type::Symbol = :line,
    xlabel::String = "default",
    ylabel::String = "default",
    title::String = "default",
)::GLMakie.Figure
    # validate
    _check_var(unit, [:rad, :deg], "unit")
    _check_var(type, [:line, :stem], "type")
    length(ph) == length(sf) || throw(
        ArgumentError(
            "Length of phases ($(length(ph))) and frequencies ($(length(sf))) must be equal.",
        ),
    )

    # set default values
    xl, yl, tt = _set_defaults(
        xlabel, ylabel, title,
        "Frequency [Hz]",
        unit === :rad ? "Phase [rad]" : "Phase [°]",
        "",
    )


    # prepare plot
    GLMakie.activate!(; title = "plot_phase()")
    fig = GLMakie.Figure(; size = (900, 450))

    # create axis with customizable properties
    ax  = GLMakie.Axis(
        fig[1, 1];
        xlabel             = xl,
        ylabel             = yl,
        title              = tt,
        xminorticksvisible = true,
        xminorticks        = IntervalsBetween(10),
        xautolimitmargin   = (0, 0),
        yautolimitmargin   = (0.05, 0.05),
        xzoomlock          = true,
        yzoomlock          = true,
        xpanlock           = true,
        ypanlock           = true,
        xrectzoom          = false,
        yrectzoom          = false,
    )
    GLMakie.xlims!(ax, _xlims(sf))
    _style_axis!(ax)


    if type === :line
        GLMakie.lines!(sf, ph; linewidth = 1, color = :black)
    else
        GLMakie.stem!(sf, ph; markersize = 10, color = :black)
    end


    return fig
end

"""
    plot_polezero(pol, zer; <keyword arguments>)

Plot a polar pole-zero map for digital filter analysis showing poles and zeros in the complex plane.

# Arguments

- `fig::Vector{Complex{Float64}}`: vector of complex pole locations in the z-plane
- `z::Vector{Complex{Float64}}`: vector of complex zero locations in the z-plane
- `m::Tuple{Real, Real}=(0, 0)`: major value to plot
- `title::String=""`: plot title
- `ticks::Bool=false`: if `true`, draw x- and y-axis ticks
- `ms::Symbol=:circle`: marker shape for drawing complex numbers (`:circle` or `:xcross`)
- `mono::Bool=false`: if `true`, use a monochrome palette

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_polezero(
    pol::Vector{Complex{Float64}},
    zer::Vector{Complex{Float64}};
    title::String = "default",
    mono::Bool = false,
)::GLMakie.Figure

    # prepare plot
    GLMakie.activate!(; title = "plot_polezero()")
    fig = GLMakie.Figure(; size = (600, 600))

    # create axis with customizable properties
    ax  = GLMakie.Axis(
        fig[1, 1];
        xlabel   = "Real",
        ylabel   = "Imag",
        aspect   = 1,
        title    = title == "default" ? "Pole-zero map" : title,
        xzoomlock = true,
        yzoomlock = true,
        xpanlock  = true,
        ypanlock  = true,
        xrectzoom = false,
        yrectzoom = false,
    )

    # plot poles
    GLMakie.scatter!(
        ax, real.(pol), imag.(pol);
        markersize = 15,
        color      = mono ? :black : :blue,
        marker     = :xcross,
    )

    # plot zeros
    GLMakie.scatter!(
        ax, real.(zer), imag.(zer);
        markersize   = 15,
        strokecolor  = mono ? :black : :blue,
        strokewidth  = 2,
        color        = :transparent,
        marker       = :circle,
    )

    # plot unit circle
    GLMakie.arc!(Point2f(0), 1, -pi, pi; linestyle = :dot, linewidth = 0.5, color = :black)


    return fig
end

"""
    plot_dwc(dc; <keyword arguments>)

Plot discrete wavelet decomposition (DWC) coefficients showing signal decomposition across different scales.

# Arguments

- `dc::Matrix{Float64}`: matrix where each row represents wavelet coefficients at a specific decomposition level
- `n::Int64=size(dc, 1) - 1`: number of decomposition levels to plot
- `t::AbstractVector`: time points corresponding to the signal values

# Returns

- `GLMakie.Figure`: the plotted figure
"""
function plot_dwc(
    dc::Matrix{Float64};
    n::Int64 = size(dc, 1) - 1,
    t::AbstractVector,
)::GLMakie.Figure
    # validate
    n > 1 || throw(ArgumentError("n must be > 1."))
    n <= size(dc, 1) - 1 || throw(ArgumentError("n must be ≤ $(size(dc, 1) - 1)."))
    size(dc, 2) == length(t) || throw(
        ArgumentError(
            "Length of t ($(length(t))) and number of dc columns ($(size(dc, 2))) must be equal.",
        ),
    )

    # set y-axis limits
    ylim   = (floor(minimum(dc); digits = 0), ceil(maximum(dc); digits = 1))
    ylim   = _tuple_max(ylim)
    yticks = unique([ylim[1], 0, ylim[2]])


    # prepare plot
    GLMakie.activate!(; title = "plot_dwc()")
    fig = GLMakie.Figure(; size = (1200, 800))
    nr  = ceil(Int64, (n + 1) / 2)


    # dc[1, :] is the original signal; coefficients start at dc[2, :]
    idx  = 2
    cidx = 1
    for idx1 in 1:nr
        cidx = 1
        for idx2 in 1:2
            if idx < n + 2
                ax = GLMakie.Axis(
                    fig[idx1, idx2];
                    xlabel             = "Time [s]",
                    ylabel             = "",
                    title              = "Coefficient #$(idx - 1)",
                    xticks             = LinearTicks(10),
                    xminorticksvisible = true,
                    xminorticks        = IntervalsBetween(10),
                    yticks             = yticks,
                    _AXIS_LOCK_KWARGS...,
                )
                GLMakie.ylims!(ax, ylim)
                _style_axis!(ax)
                GLMakie.lines!(ax, t, dc[idx, :]; color = :black)
                idx  += 1
                cidx += 1
            end
        end
    end


    # place the original signal spanning both columns
    row = cidx == 1 ? nr : nr + 1
    ax  = GLMakie.Axis(
        fig[row, 1:2];
        xlabel             = "Time [s]",
        ylabel             = "",
        title              = "Original signal",
        xticks             = LinearTicks(10),
        xminorticksvisible = true,
        xminorticks        = IntervalsBetween(10),
        yticks             = yticks,
        _AXIS_LOCK_KWARGS...,
    )
    GLMakie.xlims!(ax, _xlims(t))
    GLMakie.ylims!(ax, ylim)
    _style_axis!(ax)
    GLMakie.lines!(ax, t, dc[1, :]; color = :black)


    return fig
end
