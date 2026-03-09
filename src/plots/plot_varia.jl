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
export plot_hs
export plot_fi
export plot_phase
export plot_polezero
export plot_dwc

"""
    plot_matrix(m; <keyword arguments>)

Plot matrix.

# Arguments

  - `m::Matrix{<:Real}`
  - `xlabels::Vector{String}`
  - `ylabels::Vector{String}`
  - `xlabel::String=""`
  - `ylabel::String=""`
  - `title::String=""`
  - `cb::Bool=true`: draw colorbar
  - `cb_title::String=""`: colorbar title
  - `xrot::Int64=90`: rotate xlabels (in degrees)
  - `mono::Bool=false`: use color or gray palette

# Returns

  - `p::GLMakie.Figure`
"""
function plot_matrix(
        m::Matrix{<:Real};
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

    @assert size(m, 1) == size(m, 2) "Matrix must be square."
    @assert length(xlabels) == length(ylabels) "Lengths of xlabels ($(length(xlabels))) and ylabels ($(length(ylabels))) must be equal."
    @assert length(xlabels) == size(m, 1) "Length of xlabels ($(length(xlabels))) and matrix size $(size(m)) must be equal."
    @assert length(ylabels) == size(m, 2) "Length of ylabels ($(length(xlabels))) and matrix size $(size(m)) must be equal."

    n = size(m, 1)
    pal = mono ? :grays : :bluesreds

    # prepare plot
    GLMakie.activate!(title = "plot_matrix()")
    plot_size = (800, 800)
    p = GLMakie.Figure(size = plot_size)
    ax = GLMakie.Axis(
        p[1, 1],
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        xticks = (1:n, xlabels),
        xticklabelrotation = deg2rad(xrot),
        xticksvisible = false,
        yticks = (1:n, ylabels),
        yticksvisible = false,
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
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

    hm = GLMakie.heatmap!(m'; colormap = pal)
    if cb
        Colorbar(p[1, 2], hm; label = cb_title, labelsize = 16)
    end

    return p

end

"""
    plot_xac(m, lags; <keyword arguments>)

Plot cross/auto-covariance/correlation.

# Arguments

  - `m::AbstractVector`
  - `lags::AbstractVector`
  - `xlabel::String="Lag [s]"`
  - `ylabel::String=""`
  - `title::String=""`

# Returns

  - `p::GLMakie.Figure`
"""
function plot_xac(
        m::AbstractVector,
        lags::AbstractVector;
        xlabel::String = "Lag [s]",
        ylabel::String = "",
        title::String = ""
    )::GLMakie.Figure

    # prepare plot
    GLMakie.activate!(title = "plot_xac()")
    plot_size = (800, 300)
    p = GLMakie.Figure(size = plot_size)
    ax = GLMakie.Axis(
        p[1, 1],
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        xminorticksvisible = true,
        xminorticks = IntervalsBetween(10),
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
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

    GLMakie.lines!(lags, m; linewidth = 1, color = :black)

    return p

end

"""
    plot_histogram(s; <keyword arguments>)

Plot histogram.

# Arguments

  - `s::AbstractVector`
  - `x::Union{Nothing, Real}=nothing`: value to plot against the histogram
  - `type::Symbol`: type of histogram: regular (`:hist`) or kernel density (`:kd`)
  - `bins::Int64=15`: histogram bins: number of bins
  - `xlabel::String=""`: x-axis label
  - `ylabel::String=""`: y-axis label
  - `title::String=""`: plot title
  - `draw_mean::Bool=true`
  - `draw_median::Bool=true`
  - `mono::Bool=false`: use color or gray palette

# Returns

  - `p::GLMakie.Figure`
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

    _check_var(type, [:hist, :kd], "type")

    type === :kd && (type = :density)

    pal = mono ? :grays : :darktest

    if !isnothing(x)
        xticks = [
            round(minimum(s), digits = 2),
            round(mean(s), digits = 2),
            round(median(s), digits = 2),
            round(x, digits = 2),
            round(maximum(s), digits = 2),
        ]
    else
        xticks = [
            round(minimum(s), digits = 2),
            round(mean(s), digits = 2),
            round(median(s), digits = 2),
            round(maximum(s), digits = 2),
        ]
    end

    !draw_median && deleteat!(xticks, 3)
    !draw_mean && deleteat!(xticks, 2)
    sort!(unique(xticks))

    # prepare plot
    GLMakie.activate!(title = "plot_histogram()")
    plot_size = (800, 500)
    p = GLMakie.Figure(size = plot_size)
    ax = GLMakie.Axis(
        p[1, 1],
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        xticks = xticks,
        xticklabelrotation = pi / 2,
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        xzoomlock = true,
        yzoomlock = true,
        xpanlock = true,
        ypanlock = true,
        xrectzoom = false,
        yrectzoom = false,
    )
    GLMakie.xlims!(ax, extrema(xticks))
    ax.titlesize = 18
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12
    GLMakie.hist!(s; bins = bins, colormap = pal, strokecolor = :black, color = :grey, alpha = 0.5)

    draw_mean && (GLMakie.vlines!(round(mean(s), digits = 2); linestyle = :dot, color = :black, label = "mean"))
    draw_median && (GLMakie.vlines!(round(median(s), digits = 2); linestyle = :dash, color = :grey, label = "median"))

    if isnothing(x) != true
        if mono
            GLMakie.vlines!(x; linewidth = 2, color = :black, label = "test value")
        else
            GLMakie.vlines!([x], linewidth = 2, color = :red, label = "test value")
        end
        prop = round(cmp_stat(s, x), digits = 3)
        _info("Proportion of values > $x: $prop")
        _info("Proportion of values < $x: $(1 - prop)")
    end

    (draw_median || draw_mean || !isnoting(x)) && axislegend(; position = :rt)

    return p

end

"""
    plot_bar(s; <keyword arguments>)

Bar plot.

# Arguments

  - `s::AbstractVector`
  - `xlabels::Vector{String}`: x-ticks labels
  - `xlabel::String=""`: x-axis label
  - `ylabel::String=""`: y-axis label
  - `title::String=""`: plot title
  - `mono::Bool=false`: use color or gray palette

# Returns

  - `p::GLMakie.Figure`
"""
function plot_bar(
        s::AbstractVector;
        xlabels::Vector{String},
        xlabel::String = "",
        ylabel::String = "",
        title::String = "",
        mono::Bool = false,
    )::GLMakie.Figure

    @assert length(s) == length(xlabels) "Lengths of signal ($(length(s))) and xlabels ($(length(xlabels))) must be equal."

    pal = mono ? :grays : :darktest
    color = mono ? :lightgrey : :lightblue

    yl = if minimum(s) > 0
        (0, ceil(Int64, round(maximum(s) * 1.5, digits = 1)))
    else
        (floor(Int64, round(minimum(s) * 1.5, digits = 1)), ceil(Int64, round(maximum(s) * 1.5, digits = 1)))
    end

    # prepare plot
    GLMakie.activate!(title = "plot_bar()")
    plot_size = (800, 500)
    p = GLMakie.Figure(size = plot_size)
    ax = GLMakie.Axis(
        p[1, 1],
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        xticks = (eachindex(xlabels), xlabels),
        xautolimitmargin = (0.01, 0.01),
        yautolimitmargin = (0, 0),
        xzoomlock = true,
        yzoomlock = true,
        xpanlock = true,
        ypanlock = true,
        xrectzoom = false,
        yrectzoom = false,
    )
    GLMakie.ylims!(ax, yl)
    ax.titlesize = 18
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    GLMakie.barplot!(s; color = color, colormap = pal)

    return p

end

"""
    plot_line(s; <keyword arguments>)

Line plot.

# Arguments

  - `s::AbstractVector`
  - `xlabels::Vector{String}`: x-ticks labels
  - `xlabel::String=""`: x-axis label
  - `ylabel::String=""`: y-axis label
  - `title::String=""`: plot title

# Returns

  - `p::GLMakie.Figure`
"""
function plot_line(
        s::AbstractVector;
        xlabels::Vector{String},
        xlabel::String = "",
        ylabel::String = "",
        title::String = ""
    )::GLMakie.Figure

    @assert length(s) == length(xlabels) "Lengths of signal ($(length(s))) and xlabels ($(length(xlabels))) must be equal."

    yl = if minimum(s) > 0
        (0, ceil(Int64, round(maximum(s) * 1.5, digits = 1)))
    else
        (floor(Int64, round(minimum(s) * 1.5, digits = 1)), ceil(Int64, round(maximum(s) * 1.5, digits = 1)))
    end

    # prepare plot
    GLMakie.activate!(title = "plot_line()")
    plot_size = (800, 500)
    p = GLMakie.Figure(size = plot_size)
    ax = GLMakie.Axis(
        p[1, 1],
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        xticks = (eachindex(xlabels), xlabels),
        xautolimitmargin = (0.1, 0.1),
        yautolimitmargin = (0.1, 0.1),
        xzoomlock = true,
        yzoomlock = true,
        xpanlock = true,
        ypanlock = true,
        xrectzoom = false,
        yrectzoom = false,
    )
    GLMakie.ylims!(ax, yl)
    ax.titlesize = 18
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    GLMakie.lines!(eachindex(xlabels), s; color = :black)

    return p

end

"""
    plot_line(s; <keyword arguments>)

Line plot.

# Arguments

  - `s::AbstractArray`
  - `rlabels::Vector{String}`: signal rows labels
  - `xlabels::Vector{String}`: x-ticks labels
  - `xlabel::String=""`: x-axis label
  - `ylabel::String=""`: y-axis label
  - `title::String=""`: plot title
  - `mono::Bool=false`: use color or gray palette

# Returns

  - `p::GLMakie.Figure`
"""
function plot_line(
        s::AbstractArray;
        rlabels::Vector{String},
        xlabels::Vector{String},
        xlabel::String = "",
        ylabel::String = "",
        title::String = "",
        mono::Bool = false,
    )::GLMakie.Figure

    _chk2d(s)
    @assert size(s, 1) == length(rlabels) "Number of s columns ($(size(s, 1))) and length or rlabels ($(length(rlabels))) must be equal."
    @assert size(s, 2) == length(xlabels) "Number of s columns ($(size(s, 2))) and length of xlabels ($(length(xlabels))) must be equal."

    pal = mono ? :grays : :darktest

    yl = if minimum(s) > 0
        (0, ceil(Int64, round(maximum(s) * 1.5, digits = 1)))
    else
        (floor(Int64, round(minimum(s) * 1.5, digits = 1)), ceil(Int64, round(maximum(s) * 1.5, digits = 1)))
    end

    # prepare plot
    GLMakie.activate!(title = "plot_line()")
    plot_size = (800, 500)
    p = GLMakie.Figure(size = plot_size)
    ax = GLMakie.Axis(
        p[1, 1],
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        xticks = (eachindex(xlabels), xlabels),
        xautolimitmargin = (0.1, 0.1),
        yautolimitmargin = (0.1, 0.1),
        xzoomlock = true,
        yzoomlock = true,
        xpanlock = true,
        ypanlock = true,
        xrectzoom = false,
        yrectzoom = false,
    )
    GLMakie.ylims!(ax, yl)
    ax.titlesize = 18
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    cmap = GLMakie.resample_cmap(pal, length(xlabels))
    for idx in axes(s, 1)
        GLMakie.lines!(
            eachindex(xlabels),
            s[idx, :],
            label = rlabels[idx],
            color = cmap[idx],
            colormap = pal,
            colorrange = eachindex(xlabels),
        )
    end

    axislegend(; position = :rt)

    return p

end

"""
    plot_box(s; <keyword arguments>)

Box plot.

# Arguments

  - `s::AbstractArray`
  - `xlabels::Vector{String}`: group labels (X ticks)
  - `xlabel::String=""`: x-axis label
  - `ylabel::String=""`: y-axis label
  - `title::String=""`: plot title
  - `mono::Bool=false`: use color or gray palette

# Returns

  - `p::GLMakie.Figure`
"""
function plot_box(
        s::AbstractArray;
        xlabels::Vector{String},
        xlabel::String = "",
        ylabel::String = "",
        title::String = "",
        mono::Bool = false,
    )::GLMakie.Figure

    _chk2d(s)
    @assert size(s, 1) == length(xlabels) "Number of signal columns ($(size(s, 1))) and length of xlabels ($(length(xlabels))) must be equal."

    pal = mono ? :grays : :darktest
    color = mono ? :lightgrey : :lightblue

    yl = if minimum(s) > 0
        (0, ceil(Int64, round(maximum(s) * 1.5, digits = 1)))
    else
        (floor(Int64, round(minimum(s) * 1.5, digits = 1)), ceil(Int64, round(maximum(s) * 1.5, digits = 1)))
    end

    # prepare plot
    GLMakie.activate!(title = "plot_box()")
    plot_size = (800, 500)
    p = GLMakie.Figure(size = plot_size)
    ax = GLMakie.Axis(
        p[1, 1],
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        xticks = (eachindex(xlabels), xlabels),
        xautolimitmargin = (0.01, 0.01),
        yautolimitmargin = (0, 0),
        xzoomlock = true,
        yzoomlock = true,
        xpanlock = true,
        ypanlock = true,
        xrectzoom = false,
        yrectzoom = false,
    )
    GLMakie.ylims!(ax, yl)
    ax.titlesize = 18
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    GLMakie.boxplot!(repeat(eachindex(xlabels), size(s, 2)), s[:], color = color, colormap = pal)

    return p

end

"""
    plot_violin(s; <keyword arguments>)

Violin plot.

# Arguments

  - `s::AbstractArray`
  - `glabels::Vector{String}`: group labels (X ticks)
  - `xlabel::String=""`: x-axis label
  - `ylabel::String=""`: y-axis label
  - `title::String=""`: plot title
  - `mono::Bool=false`: use color or gray palette

# Returns

  - `p::GLMakie.Figure`
"""
function plot_violin(
        s::AbstractArray;
        xlabels::Vector{String},
        xlabel::String = "",
        ylabel::String = "",
        title::String = "",
        mono::Bool = false,
    )::GLMakie.Figure

    _chk2d(s)
    @assert size(s, 1) == length(xlabels) "Number of s columns ($(size(s, 1))) and length of xlabels ($(length(xlabels))) must be equal."

    pal = mono ? :grays : :darktest
    color = mono ? :lightgrey : :lightblue

    yl = if minimum(s) > 0
        (0, ceil(Int64, round(maximum(s) * 1.5, digits = 1)))
    else
        (floor(Int64, round(minimum(s) * 1.5, digits = 1)), ceil(Int64, round(maximum(s) * 1.5, digits = 1)))
    end

    # prepare plot
    GLMakie.activate!(title = "plot_violin()")
    plot_size = (800, 500)
    p = GLMakie.Figure(size = plot_size)
    ax = GLMakie.Axis(
        p[1, 1],
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        xticks = (eachindex(xlabels), xlabels),
        xautolimitmargin = (0.01, 0.01),
        yautolimitmargin = (0, 0),
        xzoomlock = true,
        yzoomlock = true,
        xpanlock = true,
        ypanlock = true,
        xrectzoom = false,
        yrectzoom = false,
    )
    GLMakie.ylims!(ax, yl)
    ax.titlesize = 18
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    GLMakie.violin!(
        repeat(eachindex(xlabels), size(s, 2)),
        s[:],
        strokecolor = :black,
        strokewidth = 0.25,
        #colormap=pal,
        color = color,
    )

    return p

end

"""
    plot_dots(s; <keyword arguments>)

Dots plot.

# Arguments

  - `s::AbstractArray`
  - `xlabels::Vector{String}`: group labels (X ticks)
  - `xlabel::String=""`: x-axis label
  - `ylabel::String=""`: y-axis label
  - `title::String=""`: plot title
  - `mono::Bool=false`: use color or gray palette

# Returns

  - `p::GLMakie.Figure`
"""
function plot_dots(
        s::AbstractArray;
        xlabels::Vector{String},
        xlabel::String = "",
        ylabel::String = "",
        title::String = "",
        mono::Bool = false,
    )::GLMakie.Figure

    @assert size(s, 1) == length(xlabels) "Number of signal columns ($(size(s, 1))) and length of xlabels ($(length(xlabels))) must be equal."

    pal = mono ? :grays : :darktest

    yl = if minimum(s) > 0
        (0, ceil(Int64, round(maximum(s) * 1.5, digits = 1)))
    else
        (floor(Int64, round(minimum(s) * 1.5, digits = 1)), ceil(Int64, round(maximum(s) * 1.5, digits = 1)))
    end

    # prepare plot
    GLMakie.activate!(title = "plot_dots()")
    plot_size = (800, 500)
    p = GLMakie.Figure(size = plot_size)
    ax = GLMakie.Axis(
        p[1, 1],
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        xticks = (eachindex(xlabels), xlabels),
        xautolimitmargin = (0.25, 0.25),
        yautolimitmargin = (0, 0),
        xzoomlock = true,
        yzoomlock = true,
        xpanlock = true,
        ypanlock = true,
        xrectzoom = false,
        yrectzoom = false,
    )
    GLMakie.ylims!(ax, yl)
    ax.titlesize = 18
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    cmap = GLMakie.resample_cmap(pal, length(xlabels))
    for idx in eachindex(xlabels)
        if mono
            GLMakie.scatter!(repeat([idx], size(s, 2)), s[idx, :], color = :black)
        else
            GLMakie.scatter!(
                repeat([idx], size(s, 2)), s[idx, :], color = cmap[idx], colormap = pal, colorrange = eachindex(xlabels)
            )
        end
    end

    return p

end

"""
    plot_paired(signal; <keyword arguments>)

Plot paired data.

# Arguments

  - `s::AbstractArray`
  - `xlabels::Vector{String}`: group labels (x-axis ticks)
  - `xlabel::String=""`: x-axis label
  - `ylabel::String=""`: y-axis label
  - `title::String=""`: plot title
  - `mono::Bool=false`: use color or gray palette

# Returns

  - `p::GLMakie.Figure`
"""
function plot_paired(
        s::AbstractArray;
        xlabels::Vector{String},
        xlabel::String = "",
        ylabel::String = "",
        title::String = "",
        mono::Bool = false,
    )::GLMakie.Figure

    @assert size(s, 1) == length(xlabels) "Number of signal columns ($(size(s, 1))) and length of xlabels ($(length(xlabels))) must be equal."

    pal = mono ? :grays : :darktest
    yl = if minimum(s) > 0
        (0, ceil(Int64, round(maximum(s) * 1.5, digits = 1)))
    else
        (floor(Int64, round(minimum(s) * 1.5, digits = 1)), ceil(Int64, round(maximum(s) * 1.5, digits = 1)))
    end

    # prepare plot
    GLMakie.activate!(title = "plot_paired()")
    plot_size = (800, 500)
    p = GLMakie.Figure(size = plot_size)
    ax = GLMakie.Axis(
        p[1, 1],
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        xticks = (eachindex(xlabels), xlabels),
        xautolimitmargin = (0.25, 0.25),
        yautolimitmargin = (0, 0),
        xzoomlock = true,
        yzoomlock = true,
        xpanlock = true,
        ypanlock = true,
        xrectzoom = false,
        yrectzoom = false,
    )
    GLMakie.ylims!(ax, yl)
    ax.titlesize = 18
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    cmap = GLMakie.resample_cmap(pal, length(xlabels))
    for idx in eachindex(xlabels)
        if mono
            GLMakie.scatter!(repeat([idx], size(s, 2)), s[idx, :], color = :black)
        else
            GLMakie.scatter!(
                repeat([idx], size(s, 2)), s[idx, :], color = cmap[idx], colormap = pal, colorrange = eachindex(xlabels)
            )
        end
    end

    cmap = GLMakie.resample_cmap(pal, length(xlabels))
    for idx in eachindex(xlabels)
        if mono
            GLMakie.scatter!(repeat([idx], size(s, 2)), s[idx, :], color = :black)
        else
            GLMakie.scatter!(
                repeat([idx], size(s, 2)), s[idx, :], color = cmap[idx], colormap = pal, colorrange = eachindex(xlabels)
            )
        end
    end
    for idx in axes(s, 2)
        GLMakie.lines!(eachindex(xlabels), s[:, idx], color = :black, linewidth = 0.5)
    end

    return p

end

"""
    plot_polar(s; <keyword arguments>)

Polar plot.

# Arguments

  - `s::Union{AbstractVector, AbstractMatrix}`
  - `m::Tuple{Real, Real}=(0, 0)`: major value to plot
  - `title::String=""`: plot title
  - `mono::Bool=false`: use color or gray palette
  - `ticks::Bool=false`: draw x- and y-axis ticks

# Returns

  - `p::GLMakie.Figure`
"""
function plot_polar(
        s::Union{AbstractVector, AbstractMatrix};
        m::Tuple{Real, Real} = (0, 0),
        title::String = "",
        mono::Bool = false,
        ticks::Bool = true,
    )::GLMakie.Figure

    size(s, 1) == 2 && (s = s')
    @assert length(m) == 2 "m must have exactly 2 values: phases and lengths."
    ndims(s) > 1 && @assert size(s, 2) == 2 "signal must have exactly 2 columns: phases and lengths."

    pal = mono ? :grays : :darktest

    # prepare plot
    GLMakie.activate!(title = "plot_polar()")
    plot_size = (800, 800)
    p = GLMakie.Figure(size = plot_size)
    ax = GLMakie.PolarAxis(
        p[1, 1],
        title = title,
        thetazoomlock = true,
        rzoomlock = true,
    )
    !ticks && hidespines!(ax)

    if ndims(s) == 1
        GLMakie.lines!([0, s[1]], [0, 1], linewidth = 2, color = :black)
        for idx in eachindex(s)[(begin + 1):end]
            GLMakie.lines!([0, s[idx]], [0, 1], linewidth = 2, color = :black)
        end
    else
        GLMakie.lines!([0, s[1, 1]], [0, s[1, 2]], linewidth = 2, color = :black)
        for idx in axes(s, 1)[(begin + 1):end]
            GLMakie.lines!([0, s[idx, 1]], [0, s[idx, 2]], linewidth = 2, color = :black)
        end

    end

    if m != (0, 0)
        GLMakie.lines!([0, m[1]], [0, m[2]], linewidth = 2, color = mono ? :darkgray : :red)
    end

    return p

end

"""
    plot_eros(s, f, t; <keyword arguments>)

Plot ERO (Event-Related Oscillations) spectrogram.

# Arguments

  - `sp::AbstractArray`: ERO spectrogram
  - `sf::AbstractVector`: ERO frequencies
  - `st::AbstractVector`: ERO time
  - `db::Bool=true`: whether ERO powers are normalized to dB
  - `frq::Symbol=:lin`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  - `flim::Tuple{Real, Real}=(f[1], f[end])`: frequency limit
  - `tm::Union{Int64, Vector{Int64}}=0`: time markers (in milliseconds) to be plot as vertical lines, useful for adding topoplots at these time points
  - `xlabel::String="default"`
  - `ylabel::String="default"`
  - `title::String="default"`
  - `cb::Bool=true`: draw colorbar
  - `mono::Bool=false`: use color or gray palette
  - `units::String="μV"`
  - `smooth::Bool=false`: smooth the image using Gaussian blur
  - `n::Int64=3`: kernel size of the Gaussian blur (larger kernel means more smoothing)

# Returns

  - `p::GLMakie.Figure`
"""
function plot_eros(
        sp::AbstractArray,
        sf::AbstractVector,
        st::AbstractVector;
        db::Bool = true,
        frq::Symbol = :lin,
        flim::Tuple{Real, Real} = (sf[1], sf[end]),
        tm::Union{Int64, Vector{Int64}} = 0,
        xlabel::String = "default",
        ylabel::String = "default",
        title::String = "default",
        cb::Bool = true,
        mono::Bool = false,
        units::String = "μV",
        smooth::Bool = false,
        n::Int64 = 3,
    )::GLMakie.Figure

    @assert size(sp, 1) == length(sf) "Length of sf ($(length(sf))) and number of spectrogram rows ($(size(sp, 1))) must be equal."
    @assert size(sp, 2) == length(st) "Length of st ($(length(st))) and number of spectrogram columns ($(size(sp, 2))) must be equal."
    @assert ndims(sp) == 3 "sp must have 3 dimensions."
    @assert size(sp, 3) <= 2 "sp must contain ≤ 2 epochs."
    @assert n > 0 "n must be ≥ 1."

    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(flim, extrema(sf), "flim")

    pal = mono ? :grays : :darktest
    cb_title = db ? "[dB $units^2/Hz]" : "[$units^2/Hz]"

    if frq === :lin
        if flim[2] > 100
            yt = flim[1]:10:flim[2]
        else
            yt = flim[1]:5:flim[2]
        end
    else
        if flim[1] == 0
            _warn("Lower frequency bound truncated to $(sf[2]) Hz")
            flim = (sf[2], flim[2])
        end
        yt = round.(logspace(flim[1], flim[2], nfrq), digits = 1)
    end

    if smooth
        for idx in axes(sp, 3)
            sp[:, :, idx] = @views imfilter(sp[:, :, idx], Kernel.gaussian(n))
        end
    end

    # set time markers
    if tm != 0
        if length(tm) > 1
            for tm_idx in eachindex(tm)
                @assert tm[tm_idx] / 1000 >= st[1] "tm value ($(tm[tm_idx])) is out of epoch time segment ($(st[1]):$(st[end]))."
                @assert tm[tm_idx] / 1000 <= st[end] "tm value ($(tm[tm_idx])) is out of epoch time segment ($(st[1]):$(st[end]))."
                tm[tm_idx] = vsearch(tm[tm_idx] / 1000, st)
            end
        else
            tm = vsearch(tm / 1000, st)
        end
    end

    if size(sp, 3) == 1
        xl, yl, tt = _set_defaults(
            xlabel, ylabel, title, "Time [ms]", "Frequency [Hz]", "Averaged spectrograms of epochs"
        )

        # prepare plot
        GLMakie.activate!(title = "plot_eros()")
        plot_size = (900, 450)
        p = GLMakie.Figure(size = plot_size)
        ax = GLMakie.Axis(
            p[1, 1],
            xlabel = xl,
            ylabel = yl,
            title = tt,
            xticks = LinearTicks(15),
            xminorticksvisible = true,
            xminorticks = IntervalsBetween(10),
            yticks = yt,
            yscale = frq === :lin ? identity : log,
            xautolimitmargin = (0, 0),
            yautolimitmargin = (0, 0),
            xzoomlock = true,
            yzoomlock = true,
            xpanlock = true,
            ypanlock = true,
            xrectzoom = false,
            yrectzoom = false,
        )
        GLMakie.ylims!(ax, flim)
        ax.titlesize = 18
        ax.xlabelsize = 18
        ax.ylabelsize = 18
        ax.xticklabelsize = 12
        ax.yticklabelsize = 12

        hm = GLMakie.heatmap!(ax, st, sf, sp[:, :, 1]'; colormap = pal)
        if cb
            Colorbar(p[1, 2], hm; label = cb_title, labelsize = 16)
        end

        # draw time markers
        if tm != 0
            for tm_idx in eachindex(tm)
                GLMakie.vlines!(p[1, 1], [st[tm[tm_idx]]], color = :black, linewidth = 1)
            end
        end
    else
        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [ms]", "Frequency [Hz]", "ERP spectrogram")

        # prepare plot
        GLMakie.activate!(title = "plot_eros()")
        plot_size = (1200, 800)
        p = GLMakie.Figure(size = plot_size)
        ax1 = GLMakie.Axis(
            p[1, 1],
            xlabel = xl,
            ylabel = yl,
            title = tt,
            xticks = LinearTicks(10),
            xminorticksvisible = true,
            xminorticks = IntervalsBetween(10),
            yticks = yt,
            yscale = frq === :lin ? identity : log,
            xautolimitmargin = (0, 0),
            yautolimitmargin = (0, 0),
            xzoomlock = true,
            yzoomlock = true,
            xpanlock = true,
            ypanlock = true,
            xrectzoom = false,
            yrectzoom = false,
        )
        GLMakie.ylims!(ax1, flim)
        ax1.titlesize = 18
        ax1.xlabelsize = 18
        ax1.ylabelsize = 18
        ax1.xticklabelsize = 12
        ax1.yticklabelsize = 12

        hm1 = GLMakie.heatmap!(ax1, st, sf, sp[:, :, 1]'; colormap = pal)

        if cb
            Colorbar(p[1, 2], hm1; label = cb_title, labelsize = 16)
        end

        # plot 0 v-line
        GLMakie.vlines!(ax1, [0], linestyle = :dash, linewidth = 0.5, color = :black)

        # draw time markers
        if tm != 0
            for tm_idx in eachindex(tm)
                GLMakie.vlines!(p[1, 1], [st[tm[tm_idx]]], color = :black, linewidth = 1)
            end
        end

        xl, yl, tt = _set_defaults(
            xlabel, ylabel, title, "Time [ms]", "Frequency [Hz]", "Averaged spectrograms of ERP epochs"
        )
        ax2 = GLMakie.Axis(
            p[2, 1],
            xlabel = xl,
            ylabel = yl,
            title = tt,
            xticks = LinearTicks(10),
            xminorticksvisible = true,
            xminorticks = IntervalsBetween(10),
            yticks = yt,
            yscale = frq === :lin ? identity : log,
            xautolimitmargin = (0, 0),
            yautolimitmargin = (0, 0),
            xzoomlock = true,
            yzoomlock = true,
            xpanlock = true,
            ypanlock = true,
            xrectzoom = false,
            yrectzoom = false,
        )
        GLMakie.ylims!(ax2, flim)
        ax2.titlesize = 18
        ax2.xlabelsize = 18
        ax2.ylabelsize = 18
        ax2.xticklabelsize = 12
        ax2.yticklabelsize = 12

        hm2 = GLMakie.heatmap!(ax2, st, sf, sp[:, :, 2]'; colormap = pal)

        if cb
            Colorbar(p[2, 2], hm2; label = cb_title, labelsize = 16)
        end

        # plot 0 v-line
        GLMakie.vlines!(ax2, [0], linestyle = :dash, linewidth = 0.5, color = :black)

        # draw time markers
        if tm != 0
            for tm_idx in eachindex(tm)
                GLMakie.vlines!(p[2, 1], [st[tm[tm_idx]]], color = :black, linewidth = 1)
            end
        end

    end

    return p

end

"""
    plot_erop(p, f; <keyword arguments>)

Plot ERO (Event-Related Oscillations) power-spectrum.

# Arguments

  - `p::AbstractArray`: ERO powers
  - `f::AbstractVector`: ERO frequencies
  - `db::Bool=true`: whether ERO powers are normalized to dB
  - `xlabel::String="default"`
  - `ylabel::String="default"`
  - `title::String="default"`
  - `flim::Tuple{Real, Real}=(f[1], f[end])`: frequency limit
  - `frq::Symbol=:lin`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
  - `units::String="μV"`
  - `mono::Bool=false`: use color or gray palette

# Returns

  - `p::GLMakie.Figure`
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

    _in(flim[1], (sf[1], sf[end]), "flim")
    _in(flim[2], (sf[1], sf[end]), "flim")
    @assert size(sp, 1) == length(sf) "Length of sf ($(length(sf))) and number of powers rows ($(size(sp, 1)))) must be equal."
    @assert ndims(sp) == 2 "sp must have 2 dimensions."
    @assert size(sp, 2) <= 2 "sp must contain ≤ 2 epochs."

    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(flim, extrema(sf), "flim")

    pal = mono ? :grays : :darktest

    if frq === :log && flim[1] == 0
        _warn("Lower frequency bound truncated to $(sf[2]) Hz")
        flim = (sf[2], flim[2])
    end

    if size(sp, 2) == 1
        if db
            xl, yl, tt = _set_defaults(
                xlabel, ylabel, title, "Frequency [Hz]", "Power [dB $units^2/Hz]", "Averaged power-spectra of epochs"
            )
        else
            xl, yl, tt = _set_defaults(
                xlabel, ylabel, title, "Frequency [Hz]", "Power [$units^2/Hz]", "Averaged power-spectra of epochs"
            )
        end

        # prepare plot
        GLMakie.activate!(title = "plot_erop()")
        plot_size = (900, 450)
        p = GLMakie.Figure(size = plot_size)
        ax = GLMakie.Axis(
            p[1, 1],
            xlabel = xl,
            ylabel = yl,
            title = tt,
            xminorticksvisible = true,
            xminorticks = IntervalsBetween(10),
            xscale = frq === :lin ? identity : log,
            xautolimitmargin = (0, 0),
            yautolimitmargin = (0, 0),
            xzoomlock = true,
            yzoomlock = true,
            xpanlock = true,
            ypanlock = true,
            xrectzoom = false,
            yrectzoom = false,
        )
        GLMakie.xlims!(ax, flim)
        ax.titlesize = 18
        ax.xlabelsize = 18
        ax.ylabelsize = 18
        ax.xticklabelsize = 12
        ax.yticklabelsize = 12

        # plot powers
        Makie.lines!(sf, sp[:, 1], color = :black)

    else
        if db
            xl, yl, tt = _set_defaults(
                xlabel, ylabel, title, "Frequency [Hz]", "Power [dB $units^2/Hz]", "ERP power-spectrum"
            )
        else
            xl, yl, tt = _set_defaults(
                xlabel, ylabel, title, "Frequency [Hz]", "Power [$units^2/Hz]", "ERP power-spectrum"
            )
        end

        # prepare plot
        GLMakie.activate!(title = "plot_erop()")
        plot_size = (1200, 800)
        p = GLMakie.Figure(size = plot_size)
        ax1 = GLMakie.Axis(
            p[1, 1],
            xlabel = xl,
            ylabel = yl,
            title = tt,
            xminorticksvisible = true,
            xminorticks = IntervalsBetween(10),
            xscale = frq === :lin ? identity : log,
            xautolimitmargin = (0, 0),
            yautolimitmargin = (0, 0),
            xzoomlock = true,
            yzoomlock = true,
            xpanlock = true,
            ypanlock = true,
            xrectzoom = false,
            yrectzoom = false,
        )
        GLMakie.xlims!(ax1, flim)
        ax1.titlesize = 18
        ax1.xlabelsize = 18
        ax1.ylabelsize = 18
        ax1.xticklabelsize = 12
        ax1.yticklabelsize = 12

        # plot powers
        Makie.lines!(ax1, sf, sp[:, 1], color = :black)

        if db
            xl, yl, tt = _set_defaults(
                xlabel, ylabel, title, "Frequency [Hz]", "Power [dB $units^2/Hz]", "Averaged power-spectra of epochs"
            )
        else
            xl, yl, tt = _set_defaults(
                xlabel, ylabel, title, "Frequency [Hz]", "Power [$units^2/Hz]", "Averaged power-spectra of epochs"
            )
        end

        ax2 = GLMakie.Axis(
            p[2, 1],
            xlabel = xl,
            ylabel = yl,
            title = tt,
            xminorticksvisible = true,
            xminorticks = IntervalsBetween(10),
            xscale = frq === :lin ? identity : log,
            xautolimitmargin = (0, 0),
            yautolimitmargin = (0, 0),
            xzoomlock = true,
            yzoomlock = true,
            xpanlock = true,
            ypanlock = true,
            xrectzoom = false,
            yrectzoom = false,
        )
        GLMakie.xlims!(ax2, flim)
        ax2.titlesize = 18
        ax2.xlabelsize = 18
        ax2.ylabelsize = 18
        ax2.xticklabelsize = 12
        ax2.yticklabelsize = 12

        # plot powers
        Makie.lines!(ax2, sf, sp[:, 2], color = :black)
    end

    return p

end

"""
    plot_icatopo(obj; <keyword arguments>)

Topographical plot of external ICA components.

# Arguments

  - `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
  - `ic::Matrix{Float64}`: components IC(1)..IC(n)
  - `ic_mw::Matrix{Float64}`: weighting matrix IC(1)..IC(n)
  - `ch::Union{String, Vector{String}, Regex}`: channel name(s)
  - `ic_idx::Union{Int64, Vector{Int64}, AbstractRange}=axes(ic_idx, 1)`: component(s) to plot, default is all components
  - `tpos::Union{Nothing, Real, AbstractVector}=nothing`: time point in seconds to plot, ignored if `data` is provided
  - `imethod::Symbol=:sh`: interpolation method:
      + `:sh`: Shepard
      + `:mq`: Multiquadratic
      + `:imq`: InverseMultiquadratic
      + `:tp`: ThinPlate
      + `:nn`: NearestNeighbour
      + `:ga`: Gaussian
  - `nmethod::Symbol=:minmax`: method for normalization, see `normalize()`
  - `contours::Int64=0`: number of contour levels to plot
  - `electrodes::Bools=true`: plot electrodes over topo plot
  - `ps::Symbol=:l`: plot size (`:l`: large (800×800 px), `:m`: medium (300×300 px), `:s`: small (100×100 px))

# Returns

  - `p::GLMakie.figure`
"""
function plot_icatopo(
        obj::NeuroAnalyzer.NEURO;
        ch::Union{String, Vector{String}, Regex},
        ic::Matrix{Float64},
        ic_mw::Matrix{Float64},
        ic_idx::Union{Int64, Vector{Int64}, AbstractRange} = axes(ic_idx, 1),
        tpos::Union{Nothing, Real, AbstractVector},
        imethod::Symbol = :sh,
        nmethod::Symbol = :minmax,
        contours::Int64 = 0,
        electrodes::Bool = true,
        ps::Symbol = :l,
    )::GLMakie.Figure

    p_topo = GLMakie.Figure[]
    for idx in eachindex(ic_idx)
        obj_tmp = ica_reconstruct(obj, ch = ch, ic = ic, ic_mw = ic_mw, ic_idx = idx, keep = true)
        p_tmp = plot_topo(
            obj_tmp;
            ch = ch,
            tpos = tpos,
            title = "IC $idx",
            imethod = imethod,
            nmethod = nmethod,
            contours = contours,
            electrodes = electrodes,
            ps = ps,
            cb = true,
        )
        push!(p_topo, p_tmp)
    end

    p = plot_compose(p_topo, layout = (1, length(ic_idx)))

    return p

end

"""
    plot_ci(s, s_ci_l, s_ci_h, t; <keyword arguments>)

Confidence interval plot.

# Arguments

  - `s::AbstractVector`: signal
  - `s_l::AbstractVector`: CI lower bound
  - `s_u::AbstractVector`: CI upper bound
  - `t::AbstractVector`: time points
  - `xlabel::String=""`: x-axis label
  - `ylabel::String=""`: y-axis label
  - `title::String=""`: plot title
  - `mono::Bool=false`: use color or gray palette

# Returns

  - `p::GLMakie.Figure`
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

    @assert length(s) == length(s_l) == length(s_u) "All input signals must be of the same length."

    pal = mono ? :grays : :darktest

    yl = (floor(minimum(s_l), digits = 0), ceil(maximum(s_u), digits = 0))
    yl = _tuple_max(yl)
    yticks = [yl[1], 0, yl[2]]

    # prepare plot
    GLMakie.activate!(title = "plot_ci()")
    plot_size = (800, 500)
    p = GLMakie.Figure(size = plot_size)
    ax = GLMakie.Axis(
        p[1, 1],
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        xticks = LinearTicks(10),
        yticks = yticks,
        xminorticksvisible = true,
        xminorticks = IntervalsBetween(10),
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        xzoomlock = true,
        yzoomlock = true,
        xpanlock = true,
        ypanlock = true,
        xrectzoom = false,
        yrectzoom = false,
    )
    GLMakie.ylims!(ax, yl)
    ax.titlesize = 18
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    # draw 95% CI
    Makie.band!(t, s_u, s_l; alpha = 0.25, color = :grey, strokewidth = 0.5)

    # draw signal
    Makie.lines!(t, s; color = :black, linewidth = 2)

    return p

end

"""
    plot_heatmap(m; <keyword arguments>)

Plot heatmap.

# Arguments

  - `m::AbstractMatrix`
  - `x::AbstractVector`
  - `y::AbstractVector`
  - `xlabel::String=""`: x-axis label
  - `ylabel::String=""`: y-axis label
  - `title::String=""`: plot title
  - `mono::Bool=false`: use color or gray palette
  - `cb::Bool=true`: draw colorbar
  - `cb_title::String=""`: colorbar title
  - `threshold::Union{Nothing, Real}=nothing`: if set, use threshold to mark a region
  - `threshold_type::Symbol=:neq`: rule for thresholding:
      + `:eq`: draw region is values are equal to threshold
      + `:neq`: draw region is values are not equal to threshold
      + `:geq`: draw region is values are ≥ to threshold
      + `:leq`: draw region is values are ≤ to threshold
      + `:g`: draw region is values are > to threshold
      + `:l`: draw region is values are < to threshold

# Returns

  - `p::GLMakie.Figure`
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
        threshold::Union{Nothing, Real} = nothing,
        threshold_type::Symbol = :neq,
    )::GLMakie.Figure

    @assert size(m, 1) == length(y) "Number of m rows ($(size(m, 1))) and y length ($(length(y))) must be equal."
    @assert size(m, 2) == length(x) "Number of m columns ($(size(m, 2))) and x length ($(length(x))) must be equal."

    pal = mono ? :grays : :bluesreds

    # prepare plot
    GLMakie.activate!(title = "plot_heatmap()")
    plot_size = (800, 500)
    p = GLMakie.Figure(size = plot_size)
    ax = GLMakie.Axis(
        p[1, 1],
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        xticks = LinearTicks(10),
        yticks = LinearTicks(10),
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
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

    hm = GLMakie.heatmap!(x, y, m'; colormap = pal)
    if cb
        Colorbar(p[1, 2], hm; label = cb_title, labelsize = 16)
    end

    if !isnothing(threshold)
        _, bm = seg_extract(m; threshold = threshold, threshold_type = threshold_type)
        reg = ones(size(m)) .* minimum(m)
        reg[bm] .= maximum(m)
        GLMakie.contour!(ax, x, y, reg'; levels = 1, color = :black, linewidth = 2)
    end

    return p

end

"""
    plot_imf(imf; <keyword arguments>)

Plot intrinsic mode functions (IMF), the residual and reconstructed signal.

# Arguments

  - `imf::Matrix{Float64}`: IMFs
  - `n::Int64=size(imf, 1) - 1`: number of IMFs to plot
  - `t::AbstractVector`: time points

# Returns

  - `p::GLMakie.Figure`
"""
function plot_imf(
        imf::Matrix{Float64};
        n::Int64 = size(imf, 1) - 1,
        t::AbstractVector
    )::GLMakie.Figure

    @assert n > 0 "n must be ≥ 1."
    @assert n + 1 <= size(imf, 1) "n must be ≤ $(size(imf, 1) - 1)."
    @assert size(imf, 2) == length(t) "Length of t $(size(imf, 2)) and number of imf columns ($(size(m, 2))) must be equal."

    s_restored = sum(imf; dims = 1)[:]
    imf = vcat(imf, s_restored')

    ylim = (floor(minimum(imf), digits = 0), ceil(maximum(imf), digits = 0))
    ylim = _tuple_max(ylim)
    yticks = [ylim[1], 0, ylim[2]]

    # prepare plot
    GLMakie.activate!(title = "plot_imf()")
    plot_size = (1200, 800)
    p = GLMakie.Figure(size = plot_size)

    nr = ceil(Int64, (n + 1) / 2)

    idx = 1
    cidx = 1
    for idx1 in 1:nr
        cidx = 1
        for idx2 in 1:2
            if idx <= n + 1
                ax = GLMakie.Axis(
                    p[idx1, idx2],
                    xlabel = "Time [s]",
                    ylabel = "",
                    title = idx == n + 1 ? "Residual" : "IMF: $idx",
                    xticks = LinearTicks(10),
                    xminorticksvisible = true,
                    xminorticks = IntervalsBetween(10),
                    yticks = yticks,
                    xautolimitmargin = (0, 0),
                    yautolimitmargin = (0, 0),
                    xzoomlock = true,
                    yzoomlock = true,
                    xpanlock = true,
                    ypanlock = true,
                    xrectzoom = false,
                    yrectzoom = false,
                )
                GLMakie.ylims!(ax, ylim)
                ax.titlesize = 18
                ax.xlabelsize = 18
                ax.ylabelsize = 18
                ax.xticklabelsize = 12
                ax.yticklabelsize = 12

                GLMakie.lines!(ax, t, imf[idx, :], color = :black)
                idx += 1
                cidx += 1
            end
        end
    end

    if cidx == 1
        ax = GLMakie.Axis(
            p[nr, 1:2],
            xlabel = "Time [s]",
            ylabel = "",
            title = "Reconstructed signal",
            xticks = LinearTicks(10),
            xminorticksvisible = true,
            xminorticks = IntervalsBetween(10),
            yticks = yticks,
            xautolimitmargin = (0, 0),
            yautolimitmargin = (0, 0),
            xzoomlock = true,
            yzoomlock = true,
            xpanlock = true,
            ypanlock = true,
            xrectzoom = false,
            yrectzoom = false,
        )
    else
        ax = GLMakie.Axis(
            p[nr + 1, 1:2],
            xlabel = "Time [s]",
            ylabel = "",
            title = "Reconstructed signal",
            xticks = LinearTicks(10),
            xminorticksvisible = true,
            xminorticks = IntervalsBetween(10),
            xautolimitmargin = (0, 0),
            yautolimitmargin = (0, 0),
            xzoomlock = true,
            yzoomlock = true,
            xpanlock = true,
            ypanlock = true,
            xrectzoom = false,
            yrectzoom = false,
        )
    end
    GLMakie.ylims!(ax, ylim)
    ax.titlesize = 18
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    GLMakie.lines!(ax, t, s_restored; color = :black)

    return p

end

"""
    plot_hs(sp, st; <keyword arguments>)

Plot Hilbert spectrum.

# Arguments

  - `sp::Vector{Float64}`: Hilbert transform powers
  - `st::Vector{Float64}`: time
  - `xlabel::String="default"`: x-axis label, default is Time [s]
  - `ylabel::String="default"`: y-axis label, default is Power [μV^2/Hz]
  - `title::String="default"`: plot title
  - `mono::Bool=false`: use color or gray palette

# Returns

  - `p::GLMakie.Figure`
"""
function plot_hs(
        sp::Vector{Float64},
        st::Vector{Float64};
        xlabel::String = "default",
        ylabel::String = "default",
        title::String = "default",
        mono::Bool = false,
    )::GLMakie.Figure

    @assert length(sp) == length(st) "Length of powers ($(length(sp))) and time points ($(length(st))) must be equal."

    pal = mono ? :grays : :darktest

    xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "Power [μV^2/Hz]", "")

    # prepare plot
    GLMakie.activate!(title = "plot_hs()")
    plot_size = (900, 450)
    p = GLMakie.Figure(size = plot_size)
    ax = GLMakie.Axis(
        p[1, 1],
        xlabel = xl,
        ylabel = yl,
        title = tt,
        xticks = LinearTicks(10),
        xminorticksvisible = true,
        xminorticks = IntervalsBetween(10),
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        xzoomlock = true,
        yzoomlock = true,
        xpanlock = true,
        ypanlock = true,
        xrectzoom = false,
        yrectzoom = false,
    )
    GLMakie.xlims!(ax, _xlims(st))
    ax.titlesize = 18
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    # plot powers
    GLMakie.lines!(st, sp, linewidth = 1, color = :black)

    return p

end

"""
    plot_fi(fi, st; <keyword arguments>)

Plot instantaneous frequencies.

# Arguments

  - `fi::Vector{Float64}`: instantaneous frequencies
  - `st::Vector{Float64}`: time
  - `xlabel::String="default"`: x-axis label, default is Time [s]
  - `ylabel::String="default"`: y-axis label, default is Power [μV^2/Hz]
  - `title::String="default"`: plot title

# Returns

  - `p::GLMakie.Figure`
"""
function plot_fi(
        fi::Vector{Float64},
        st::Vector{Float64};
        xlabel::String = "default",
        ylabel::String = "default",
        title::String = "default",
    )::GLMakie.Figure

    @assert length(fi) == length(st) "Length of frequencies ($(length(fi))) and time points ($(length(st))) must be equal."

    xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [s]", "Frequency [Hz]", "")

    # prepare plot
    GLMakie.activate!(title = "plot_fi()")
    plot_size = (900, 450)
    p = GLMakie.Figure(size = plot_size)
    ax = GLMakie.Axis(
        p[1, 1],
        xlabel = xl,
        ylabel = yl,
        title = tt,
        xticks = LinearTicks(10),
        xminorticksvisible = true,
        xminorticks = IntervalsBetween(10),
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0, 0),
        xzoomlock = true,
        yzoomlock = true,
        xpanlock = true,
        ypanlock = true,
        xrectzoom = false,
        yrectzoom = false,
    )
    GLMakie.xlims!(ax, _xlims(st))
    ax.titlesize = 18
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    # plot frequencies
    GLMakie.lines!(st, fi; linewidth = 1, color = :black)

    return p

end

"""
    plot_phase(ph, sf; <keyword arguments>)

Plot phases.

# Arguments

  - `ph::Vector{Float64}`: phases
  - `sf::Vector{Float64}`: frequencies
  - `unit::Symbol=:rad`: phase unit (radians or degrees)
  - `type::Symbol=:line`: plot type (`:line`: line, `:stem`: stems)
  - `xlabel::String="default"`: x-axis label, default is Time [s]
  - `ylabel::String="default"`: y-axis label, default is Power [μV^2/Hz]
  - `title::String="default"`: plot title

# Returns

  - `p::GLMakie.Figure`
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

    _check_var(unit, [:rad, :deg], "unit")
    _check_var(type, [:line, :stem], "type")
    @assert length(ph) == length(sf) "Length of phases ($(length(fi))) and frequencies ($(length(st))) must be equal."

    xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Frequency [Hz]", unit === :rad ? "Phase [rad]" : "Phase [°]", "")

    # prepare plot
    GLMakie.activate!(title = "plot_phase()")
    plot_size = (900, 450)
    p = GLMakie.Figure(size = plot_size)
    ax = GLMakie.Axis(
        p[1, 1],
        xlabel = xl,
        ylabel = yl,
        title = tt,
        xminorticksvisible = true,
        xminorticks = IntervalsBetween(10),
        xautolimitmargin = (0, 0),
        yautolimitmargin = (0.05, 0.05),
        xzoomlock = true,
        yzoomlock = true,
        xpanlock = true,
        ypanlock = true,
        xrectzoom = false,
        yrectzoom = false,
    )
    GLMakie.xlims!(ax, _xlims(sf))
    ax.titlesize = 18
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    # plot phases
    if type === :line
        GLMakie.lines!(sf, ph, linewidth = 1, color = :black)
    else
        GLMakie.stem!(sf, ph, markersize = 10, color = :black)
    end

    return p

end

"""
    plot_polezero(pol, zer; <keyword arguments>)

Polar pole-zero map.

# Arguments

  - `p::Vector{Complex{Float64}}`: vector of poles
  - `z::Vector{Complex{Float64}}`: vector of zeros
  - `m::Tuple{Real, Real}=(0, 0)`: major value to plot
  - `title::String=""`: plot title
  - `ticks::Bool=false`: draw x- and y-axis ticks
  - `ms::Symbol=:circle`: marker shape for drawing complex numbers (`:circle` or `:xcross`)
  - `mono::Bool=false`: use color or gray palette

# Returns

  - `p::GLMakie.Figure`
"""
function plot_polezero(
        pol::Vector{Complex{Float64}},
        zer::Vector{Complex{Float64}};
        title::String = "default",
        mono::Bool = false
    )::GLMakie.Figure

    # prepare plot
    GLMakie.activate!(title = "plot_polezero()")
    plot_size = (600, 600)
    p = GLMakie.Figure(size = plot_size)
    ax = GLMakie.Axis(
        p[1, 1],
        xlabel = "Real",
        ylabel = "Imag",
        aspect = 1,
        title = title == "default" ? "Pole-zero map" : title,
        xzoomlock = true,
        yzoomlock = true,
        xpanlock = true,
        ypanlock = true,
        xrectzoom = false,
        yrectzoom = false,
    )
    GLMakie.scatter!(ax, real.(pol), imag.(pol), markersize = 15, color = mono ? :black : :blue, marker = :xcross)
    GLMakie.scatter!(
        ax,
        real.(zer),
        imag.(zer),
        markersize = 15,
        strokecolor = mono ? :black : :blue,
        strokewidth = 2,
        color = :transparent,
        marker = :circle,
    )
    GLMakie.arc!(Point2f(0), 1, -pi, pi; linestyle = :dot, linewidth = 0.5, color = :black)

    return p

end


"""
    plot_dwc(dc; <keyword arguments>)

Plot discrete wavelet decomposition coefficients.

# Arguments

  - `dc::Matrix{Float64}`: coefficients
  - `n::Int64=size(dc, 1) - 1`: number of coefficients to plot
  - `t::AbstractVector`: time points

# Returns

  - `p::GLMakie.Figure`
"""
function plot_dwc(
        dc::Matrix{Float64};
        n::Int64 = size(dc, 1) - 1,
        t::AbstractVector
    )::GLMakie.Figure

    @assert n > 1 "n must be > 1."
    @assert n <= size(dc, 1) - 1 "n must be ≤ $(size(dc, 1) - 1)."
    @assert size(dc, 2) == length(t) "Length of t $(size(dc, 2)) and number of dc columns ($(size(m, 2))) must be equal."

    ylim = (floor(minimum(dc), digits = 0), ceil(maximum(dc), digits = 1))
    ylim = _tuple_max(ylim)
    yticks = unique([ylim[1], 0, ylim[2]])

    # prepare plot
    GLMakie.activate!(title = "plot_dwc()")
    plot_size = (1200, 800)
    p = GLMakie.Figure(size = plot_size)

    nr = ceil(Int64, (n + 1) / 2)

    idx = 2
    cidx = 1
    for idx1 in 1:nr
        cidx = 1
        for idx2 in 1:2
            if idx < n + 2
                ax = GLMakie.Axis(
                    p[idx1, idx2],
                    xlabel = "Time [s]",
                    ylabel = "",
                    title = "Coefficient #$(idx - 1)",
                    xticks = LinearTicks(10),
                    xminorticksvisible = true,
                    xminorticks = IntervalsBetween(10),
                    yticks = yticks,
                    xautolimitmargin = (0, 0),
                    yautolimitmargin = (0, 0),
                    xzoomlock = true,
                    yzoomlock = true,
                    xpanlock = true,
                    ypanlock = true,
                    xrectzoom = false,
                    yrectzoom = false,
                )
                GLMakie.ylims!(ax, ylim)
                ax.titlesize = 18
                ax.xlabelsize = 18
                ax.ylabelsize = 18
                ax.xticklabelsize = 12
                ax.yticklabelsize = 12

                GLMakie.lines!(ax, t, dc[idx, :], color = :black)
                idx += 1
                cidx += 1
            end
        end
    end

    if cidx == 1
        ax = GLMakie.Axis(
            p[nr, 1:2],
            xlabel = "Time [s]",
            ylabel = "",
            title = "Original signal",
            xticks = LinearTicks(10),
            xminorticksvisible = true,
            xminorticks = IntervalsBetween(10),
            yticks = yticks,
            xautolimitmargin = (0, 0),
            yautolimitmargin = (0, 0),
            xzoomlock = true,
            yzoomlock = true,
            xpanlock = true,
            ypanlock = true,
            xrectzoom = false,
            yrectzoom = false,
        )
        GLMakie.ylims!(ax, ylim)
        ax.titlesize = 18
        ax.xlabelsize = 18
        ax.ylabelsize = 18
        ax.xticklabelsize = 12
        ax.yticklabelsize = 12

        GLMakie.lines!(ax, t, dc[1, :], color = :black)
    else
        ax = GLMakie.Axis(
            p[nr + 1, 1:2],
            xlabel = "Time [s]",
            ylabel = "",
            title = "Original signal",
            xticks = LinearTicks(10),
            xminorticksvisible = true,
            xminorticks = IntervalsBetween(10),
            xautolimitmargin = (0, 0),
            yautolimitmargin = (0, 0),
            xzoomlock = true,
            yzoomlock = true,
            xpanlock = true,
            ypanlock = true,
            xrectzoom = false,
            yrectzoom = false,
        )
        GLMakie.xlims!(ax, _xlims(t))
        GLMakie.ylims!(ax, ylim)
        ax.titlesize = 18
        ax.xlabelsize = 18
        ax.ylabelsize = 18
        ax.xticklabelsize = 12
        ax.yticklabelsize = 12

        GLMakie.lines!(ax, t, dc[1, :], color = :black)
    end

    return p

end
