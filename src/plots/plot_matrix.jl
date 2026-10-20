export plot_matrix
export plot_heatmap
export plot_xac

# ---------------------------------------------------------------------------

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
        xlabel             = xlabel,
        ylabel             = ylabel,
        title              = title,
        xticks             = (1:n, xlabels),
        xticklabelrotation = deg2rad(xrot),
        xticksvisible      = false,
        yticks             = (1:n, ylabels),
        yticksvisible      = false,
        xautolimitmargin   = (0, 0),
        yautolimitmargin   = (0, 0),
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
    ax = GLMakie.Axis(
        fig[1, 1];
        xlabel             = xlabel,
        ylabel             = ylabel,
        title              = title,
        xminorticksvisible = true,
        xminorticks        = IntervalsBetween(10),
        xautolimitmargin   = (0, 0),
        yautolimitmargin   = (0, 0),
        _AXIS_LOCK_KWARGS...,
    )
    _style_axis!(ax)

    GLMakie.lines!(lags, m; linewidth = 1, color = :black)

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
    ax = GLMakie.Axis(
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
        _, bm = seg_extract(m; threshold = threshold, threshold_type = threshold_type)
        reg = ones(size(m)) .* minimum(m)
        reg[bm] .= maximum(m)
        GLMakie.contour!(ax, x, y, reg'; levels = 1, color = :black, linewidth = 2)
    end

    return fig
end
