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
- `cb::Bool=true`: draw color bar
- `cb_title::String=""`: color bar title
- `xrot::Int64=90`: rotate xlabels (in degrees)
- `mono::Bool=false`: use color or gray palette
- `kwargs`: optional arguments for plotting

# Returns

- `p::GLMakie.Figure`
"""
function plot_matrix(m::Matrix{<:Real}; xlabels::Vector{String}, ylabels::Vector{String}, xlabel::String="", ylabel::String="", title::String="", cb::Bool=true, cb_title::String="", xrot::Int64=90, mono::Bool=false, kwargs...)::GLMakie.Figure

    @assert size(m, 1) == size(m, 2) "Matrix must be square."
    @assert length(xlabels) == length(ylabels) "Lengths of xlabels ($(length(xlabels))) and ylabels ($(length(ylabels))) must be equal."
    @assert length(xlabels) == size(m, 1) "Length of xlabels ($(length(xlabels))) and matrix size $(size(m)) must be equal."
    @assert length(ylabels) == size(m, 2) "Length of ylabels ($(length(xlabels))) and matrix size $(size(m)) must be equal."

    n = size(m, 1)
    pal = mono ? :grays : :bluesreds

    # prepare plot
    plot_size = (800, 800)
    p = GLMakie.Figure(size=plot_size)
    ax = GLMakie.Axis(p[1, 1],
                      xlabel=xlabel,
                      ylabel=ylabel,
                      title=title,
                      xticks=(1:n, xlabels),
                      xticklabelrotation=deg2rad(xrot),
                      xticksvisible=false,
                      yticks=(1:n, ylabels),
                      yticksvisible=false,
                      xautolimitmargin=(0, 0),
                      yautolimitmargin=(0, 0);
                      kwargs...)
    ax.titlesize = 20
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    hm = GLMakie.heatmap!(m,
                          colormap=pal)
    if cb
        Colorbar(p[1, 2],
                 hm,
                 label=cb_title,
                 labelsize=18)
    end

    return p

end

"""
    plot_xac(m, lags; <keyword arguments>)

Plot cross/auto-covariance/correlation.

# Arguments

- `m::Abstractvector`: covariance matrix
- `lags::AbstractVector`: covariance lags
- `xlabel::String="lag"`
- `ylabel::String=""`
- `title::String=""`
- `cb_title::String=""`: color bar title
- `mono::Bool=false`: use color or gray palette
- `kwargs`: optional arguments for plotting

# Returns

- `p::GLMakie.Figure`
"""
function plot_xac(m::AbstractVector, lags::AbstractVector; xlabel::String="lag [s]", ylabel::String="", title::String="", mono::Bool=false, kwargs...)::GLMakie.Figure

    if minimum(m) >= -1.0 && maximum(m) <= 1.0
        ylim = (-1.0, 1.0)
    else
        ylim = extrema(m)
    end

    pal = mono ? :grays : :darktest
    r = length(string(lags[1])) > 4 ? 90 : 0

    # prepare plot
    plot_size = (800, 300)
    p = GLMakie.Figure(size=plot_size)
    ax = GLMakie.Axis(p[1, 1],
                      xlabel=xlabel,
                      ylabel=ylabel,
                      title=title,
                      xticks=lags,
                      xticklabelrotation=deg2rad(r),
                      xautolimitmargin=(0, 0),
                      yautolimitmargin=(0, 0);
                      kwargs...)
    ax.titlesize = 20
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12
    GLMakie.ylims!(ax, ylim)

    GLMakie.lines!(lags,
                   m,
                   linewidth=0.5,
                   color=:black)

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
- `xlabel::String=""`: X axis label
- `ylabel::String=""`: Y axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `draw_mean::Bool=true`
- `draw_median::Bool=true`
- `kwargs`: optional arguments for plotting

# Returns

- `p::GLMakie.Figure`
"""
function plot_histogram(s::AbstractVector, x::Union{Nothing, Real}=nothing; type::Symbol=:hist, bins::Int64=15, xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, draw_mean::Bool=true, draw_median::Bool=true, kwargs...)::GLMakie.Figure

    _check_var(type, [:hist, :kd], "type")

    type === :kd && (type = :density)

    pal = mono ? :grays : :darktest

    if !isnothing(x)
        xticks = [round(minimum(s), digits=2), round(mean(s), digits=2), round(median(s), digits=2), round(x, digits=2), round(maximum(s), digits=2)]
    else
        xticks = [round(minimum(s), digits=2), round(mean(s), digits=2), round(median(s), digits=2), round(maximum(s), digits=2)]
    end

    !draw_median && deleteat!(xticks, 3)
    !draw_mean && deleteat!(xticks, 2)
    sort!(unique(xticks))

    # prepare plot
    plot_size = (800, 800)
    p = GLMakie.Figure(size=plot_size)
    ax = GLMakie.Axis(p[1, 1],
                      xlabel=xlabel,
                      ylabel=ylabel,
                      title=title,
                      xticks=xticks,
                      xticklabelrotation=pi/2,
                      xautolimitmargin=(0, 0),
                      yautolimitmargin=(0, 0);
                      kwargs...)
    GLMakie.xlims!(ax, extrema(xticks))
    ax.titlesize = 20
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12
    GLMakie.hist!(s,
                  bins=bins,
                  colormap=pal,
                  strokecolor=:black,
                  color=:grey,
                  alpha=0.5)

    draw_mean && (GLMakie.vlines!(round(mean(s), digits=2), linestyle=:dot, color=:black, label="mean"))
    draw_median && (GLMakie.vlines!(round(median(s), digits=2), linestyle=:dash, color=:grey, label="median"))

    if isnothing(x) != true
        if mono
            GLMakie.vlines!(x, linewidth=2, color=:black, label="test value")
        else
            GLMakie.vlines!([x], linewidth=2, color=:red, label="test value")
        end
        prop = round(cmp_stat(s, x), digits=3)
        _info("Proportion of values > $x: $prop")
        _info("Proportion of values < $x: $(1 - prop)")
    end

    (draw_median || draw_mean || !isnoting(x)) && axislegend(position = :rt)

    return p

end

"""
    plot_bar(s; <keyword arguments>)

Bar plot.

# Arguments

- `s::AbstractVector`
- `xlabels::Vector{String}`: x-ticks labels
- `xlabel::String=""`: X axis label
- `ylabel::String=""`: Y axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `kwargs`: optional arguments for plotting

# Returns

- `p::GLMakie.Figure`
"""
function plot_bar(s::AbstractVector; xlabels::Vector{String}, xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)::GLMakie.Figure

    @assert length(s) == length(xlabels) "Lengths of signal ($(length(s))) and xlabels ($(length(xlabels))) must be equal."

    pal = mono ? :grays : :darktest
    color = mono ? :lightgrey : :lightblue

    yl = minimum(s) > 0 ? (0, ceil(Int64, round(maximum(s) * 1.5, digits=1))) : ylims = (floor(Int64, round(minimum(s) * 1.5, digits=1)), ceil(Int64, round(maximum(s) * 1.5, digits=1)))

    # prepare plot
    plot_size = (800, 500)
    p = GLMakie.Figure(size=plot_size)
    ax = GLMakie.Axis(p[1, 1],
                      xlabel=xlabel,
                      ylabel=ylabel,
                      title=title,
                      xticks=(eachindex(xlabels), xlabels),
                      xautolimitmargin=(0.01, 0.01),
                      yautolimitmargin=(0, 0);
                      kwargs...)
    GLMakie.ylims!(ax, yl)
    ax.titlesize = 20
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    GLMakie.barplot!(s,
                     color=color,
                     colormap=pal)

    return p

end

"""
    plot_line(s; <keyword arguments>)

Line plot.

# Arguments

- `s::AbstractVector`
- `xlabels::Vector{String}`: x-ticks labels
- `xlabel::String=""`: X axis label
- `ylabel::String=""`: Y axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `kwargs`: optional arguments for plotting

# Returns

- `p::GLMakie.Figure`
"""
function plot_line(s::AbstractVector; xlabels::Vector{String}, xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)::GLMakie.Figure

    @assert length(s) == length(xlabels) "Lengths of signal ($(length(s))) and xlabels ($(length(xlabels))) must be equal."

    pal = mono ? :grays : :darktest

    yl = minimum(s) > 0 ? (0, ceil(Int64, round(maximum(s) * 1.5, digits=1))) : ylims = (floor(Int64, round(minimum(s) * 1.5, digits=1)), ceil(Int64, round(maximum(s) * 1.5, digits=1)))

    # prepare plot
    plot_size = (800, 500)
    p = GLMakie.Figure(size=plot_size)
    ax = GLMakie.Axis(p[1, 1],
                      xlabel=xlabel,
                      ylabel=ylabel,
                      title=title,
                      xticks=(eachindex(xlabels), xlabels),
                      xautolimitmargin=(0.1, 0.1),
                      yautolimitmargin=(0.1, 0.1);
                      kwargs...)
    GLMakie.ylims!(ax, yl)
    ax.titlesize = 20
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    GLMakie.lines!(eachindex(xlabels),
                   s,
                   color=:black,
                   colormap=pal)

    return p

end

"""
    plot_line(s; <keyword arguments>)

Line plot.

# Arguments

- `s::AbstractArray`
- `rlabels::Vector{String}`: signal rows labels
- `xlabels::Vector{String}`: x-ticks labels
- `xlabel::String=""`: X axis label
- `ylabel::String=""`: Y axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `kwargs`: optional arguments for plotting

# Returns

- `p::GLMakie.Figure`
"""
function plot_line(s::AbstractArray; rlabels::Vector{String}, xlabels::Vector{String}, xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)::GLMakie.Figure

    _chk2d(s)
    @assert size(s, 1) == length(rlabels) "Number of s columns ($(size(s, 1))) and length or rlabels ($(length(rlabels))) must be equal."
    @assert size(s, 2) == length(xlabels) "Number of s columns ($(size(s, 2))) and length of xlabels ($(length(xlabels))) must be equal."

    pal = mono ? :grays : :darktest

    yl = minimum(s) > 0 ? (0, ceil(Int64, round(maximum(s) * 1.5, digits=1))) : ylims = (floor(Int64, round(minimum(s) * 1.5, digits=1)), ceil(Int64, round(maximum(s) * 1.5, digits=1)))

    # prepare plot
    plot_size = (800, 500)
    p = GLMakie.Figure(size=plot_size)
    ax = GLMakie.Axis(p[1, 1],
                      xlabel=xlabel,
                      ylabel=ylabel,
                      title=title,
                      xticks=(eachindex(xlabels), xlabels),
                      xautolimitmargin=(0.1, 0.1),
                      yautolimitmargin=(0.1, 0.1);
                      kwargs...)
    GLMakie.ylims!(ax, yl)
    ax.titlesize = 20
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    cmap = GLMakie.resample_cmap(pal, length(xlabels))
    for idx in axes(s, 1)
        GLMakie.lines!(eachindex(xlabels),
                       s[idx, :],
                       label=rlabels[idx],
                       color=cmap[idx],
                       colormap=pal,
                       colorrange=eachindex(xlabels))
    end

    axislegend(position = :rt)

    return p

end

"""
    plot_box(s; <keyword arguments>)

Box plot.

# Arguments

- `s::AbstractArray`
- `xlabels::Vector{String}`: group labels (X ticks)
- `xlabel::String=""`: X axis label
- `ylabel::String=""`: Y axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `kwargs`: optional arguments for plotting

# Returns

- `p::GLMakie.Figure`
"""
function plot_box(s::AbstractArray; xlabels::Vector{String}, xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)::GLMakie.Figure

    _chk2d(s)
    @assert size(s, 1) == length(xlabels) "Number of signal columns ($(size(s, 1))) and length of xlabels ($(length(xlabels))) must be equal."

    pal = mono ? :grays : :darktest
    color = mono ? :lightgrey : :lightblue

    yl = minimum(s) > 0 ? (0, ceil(Int64, round(maximum(s) * 1.5, digits=1))) : ylims = (floor(Int64, round(minimum(s) * 1.5, digits=1)), ceil(Int64, round(maximum(s) * 1.5, digits=1)))

    # prepare plot
    plot_size = (800, 500)
    p = GLMakie.Figure(size=plot_size)
    ax = GLMakie.Axis(p[1, 1],
                      xlabel=xlabel,
                      ylabel=ylabel,
                      title=title,
                      xticks=(eachindex(xlabels), xlabels),
                      xautolimitmargin=(0.01, 0.01),
                      yautolimitmargin=(0, 0);
                      kwargs...)
    GLMakie.ylims!(ax, yl)
    ax.titlesize = 20
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    GLMakie.boxplot!(repeat(eachindex(xlabels), size(s, 2)),
                     s[:],
                     color=color,
                     colormap=pal)

    return p

end

"""
    plot_violin(s; <keyword arguments>)

Violin plot.

# Arguments

- `s::AbstractArray`
- `glabels::Vector{String}`: group labels (X ticks)
- `xlabel::String=""`: X axis label
- `ylabel::String=""`: Y axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `kwargs`: optional arguments for plotting

# Returns

- `p::GLMakie.Figure`
"""
function plot_violin(s::AbstractArray; xlabels::Vector{String}, xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)::GLMakie.Figure

    _chk2d(s)
    @assert size(s, 1) == length(xlabels) "Number of s columns ($(size(s, 1))) and length of xlabels ($(length(xlabels))) must be equal."

    pal = mono ? :grays : :darktest
    color = mono ? :lightgrey : :lightblue

    yl = minimum(s) > 0 ? (0, ceil(Int64, round(maximum(s) * 1.5, digits=1))) : ylims = (floor(Int64, round(minimum(s) * 1.5, digits=1)), ceil(Int64, round(maximum(s) * 1.5, digits=1)))

    # prepare plot
    plot_size = (800, 500)
    p = GLMakie.Figure(size=plot_size)
    ax = GLMakie.Axis(p[1, 1],
                      xlabel=xlabel,
                      ylabel=ylabel,
                      title=title,
                      xticks=(eachindex(xlabels), xlabels),
                      xautolimitmargin=(0.01, 0.01),
                      yautolimitmargin=(0, 0);
                      kwargs...)
    GLMakie.ylims!(ax, yl)
    ax.titlesize = 20
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    GLMakie.violin!(repeat(eachindex(xlabels), size(s, 2)),
                    s[:],
                    strokecolor=:black,
                    strokewidth=0.25,
                    color=color)

    return p

end

"""
    plot_dots(s; <keyword arguments>)

Dots plot.

# Arguments

- `s::AbstractArray`
- `xlabels::Vector{String}`: group labels (X ticks)
- `xlabel::String=""`: X axis label
- `ylabel::String=""`: Y axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `kwargs`: optional arguments for plotting

# Returns

- `p::GLMakie.Figure`
"""
function plot_dots(s::AbstractArray; xlabels::Vector{String}, xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)::GLMakie.Figure

    @assert size(s, 1) == length(xlabels) "Number of signal columns ($(size(s, 1))) and length of xlabels ($(length(xlabels))) must be equal."

    pal = mono ? :grays : :darktest

    yl = minimum(s) > 0 ? (0, ceil(Int64, round(maximum(s) * 1.5, digits=1))) : ylims = (floor(Int64, round(minimum(s) * 1.5, digits=1)), ceil(Int64, round(maximum(s) * 1.5, digits=1)))

    # prepare plot
    plot_size = (800, 500)
    p = GLMakie.Figure(size=plot_size)
    ax = GLMakie.Axis(p[1, 1],
                      xlabel=xlabel,
                      ylabel=ylabel,
                      title=title,
                      xticks=(eachindex(xlabels), xlabels),
                      xautolimitmargin=(0.25, 0.25),
                      yautolimitmargin=(0, 0);
                      kwargs...)
    GLMakie.ylims!(ax, yl)
    ax.titlesize = 20
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    cmap = GLMakie.resample_cmap(pal, length(xlabels))
    for idx in eachindex(xlabels)
        if !mono
            GLMakie.scatter!(repeat([idx], size(s, 2)),
                             s[idx, :],
                             color=cmap[idx],
                             colormap=pal,
                             colorrange=eachindex(xlabels))
        else
            GLMakie.scatter!(repeat([idx], size(s, 2)),
                             s[idx, :],
                             color=:black)
        end
    end

    return p

end

"""
    plot_paired(signal; <keyword arguments>)

Plot paired data.

# Arguments

- `s::AbstractArray`
- `xlabels::Vector{String}`: group labels (X ticks)
- `xlabel::String=""`: X axis label
- `ylabel::String=""`: Y axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `kwargs`: optional arguments for plotting

# Returns

- `p::GLMakie.Figure`
"""
function plot_paired(s::AbstractArray; xlabels::Vector{String}, xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)::GLMakie.Figure

    @assert size(s, 1) == length(xlabels) "Number of signal columns ($(size(s, 1))) and length of xlabels ($(length(xlabels))) must be equal."

    pal = mono ? :grays : :darktest
    yl = minimum(s) > 0 ? (0, ceil(Int64, round(maximum(s) * 1.5, digits=1))) : ylims = (floor(Int64, round(minimum(s) * 1.5, digits=1)), ceil(Int64, round(maximum(s) * 1.5, digits=1)))

    # prepare plot
    plot_size = (800, 500)
    p = GLMakie.Figure(size=plot_size)
    ax = GLMakie.Axis(p[1, 1],
                      xlabel=xlabel,
                      ylabel=ylabel,
                      title=title,
                      xticks=(eachindex(xlabels), xlabels),
                      xautolimitmargin=(0.25, 0.25),
                      yautolimitmargin=(0, 0);
                      kwargs...)
    GLMakie.ylims!(ax, yl)
    ax.titlesize = 20
    ax.xlabelsize = 18
    ax.ylabelsize = 18
    ax.xticklabelsize = 12
    ax.yticklabelsize = 12

    cmap = GLMakie.resample_cmap(pal, length(xlabels))
    for idx in eachindex(xlabels)
        if !mono
            GLMakie.scatter!(repeat([idx], size(s, 2)),
                             s[idx, :],
                             color=cmap[idx],
                             colormap=pal,
                             colorrange=eachindex(xlabels))
        else
            GLMakie.scatter!(repeat([idx], size(s, 2)),
                             s[idx, :],
                             color=:black)
        end
    end

    cmap = GLMakie.resample_cmap(pal, length(xlabels))
    for idx in eachindex(xlabels)
        if !mono
            GLMakie.scatter!(repeat([idx], size(s, 2)),
                             s[idx, :],
                             color=cmap[idx],
                             colormap=pal,
                             colorrange=eachindex(xlabels))
        else
            GLMakie.scatter!(repeat([idx], size(s, 2)),
                             s[idx, :],
                             color=:black)
        end
    end
    for idx in axes(s, 2)
        GLMakie.lines!(eachindex(xlabels),
                       s[:, idx],
                       color=:black,
                       linewidth=0.5)
    end

    return p

end

"""
    plot_polar(s; <keyword arguments>)

Polar plot.

# Arguments

- `s::Union{AbstractVector, AbstractArray}`
- `m::Tuple{Real, Real}=(0, 0)`: major value to plot
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `ticks::Bool=false`: draw X and Y ticks
- `kwargs`: optional arguments for plotting

# Returns

- `p::GLMakie.Figure`
"""
function plot_polar(s::Union{AbstractVector, AbstractArray}; m::Tuple{Real, Real}=(0, 0), title::String="", mono::Bool=false, ticks::Bool=true, kwargs...)::GLMakie.Figure

    @assert length(m) == 2 "m must have exactly 2 values: phases and lengths."
    ndims(s) > 1 && @assert size(s, 2) == 2 "signal must have exactly 2 columns: phases and lengths."

    pal = mono ? :grays : :darktest

    # prepare plot
    plot_size = (800, 800)
    p = GLMakie.Figure(size=plot_size)
    ax = GLMakie.PolarAxis(p[1, 1],
                           title=title;
                           kwargs...)
    !ticks && hidespines!(ax)

    if ndims(s) == 1
        GLMakie.lines!([0, s[1]],
                       [0, 1],
                       linewidth=2,
                       color=:black)
        for idx in eachindex(s)[(begin + 1):end]
            GLMakie.lines!([0, s[idx]],
                           [0, 1],
                           linewidth=2,
                           color=:black)
        end
    else
        GLMakie.lines!([0, s[1, 1]],
                       [0, s[1, 2]],
                       linewidth=2,
                       color=:black)
        for idx in axes(s, 1)[(begin + 1):end]
            GLMakie.lines!([0, s[idx, 1]],
                           [0, s[idx, 2]],
                           linewidth=2,
                           color=:black)
        end

    end

    if m != (0, 0)
            GLMakie.lines!([0, m[1]],
                           [0, m[2]],
                           linewidth=2,
                           color=:red)
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
- `frq_lim::Tuple{Real, Real}=(f[1], f[end])`: frequency limit for the Y axis
- `tm::Union{Int64, Vector{Int64}}=0`: time markers (in milliseconds) to be plot as vertical lines, useful for adding topoplots at these time points
- `xlabel::String="default"`
- `ylabel::String="default"`
- `title::String="default"`
- `cb::Bool=true`: draw color bar
- `mono::Bool=false`: use color or gray palette
- `units::String="μV"`
- `smooth::Bool=false`: smooth the image using Gaussian blur
- `n::Int64=3`: kernel size of the Gaussian blur (larger kernel means more smoothing)
- `kwargs`: optional arguments for plotting

# Returns

- `p::GLMakie.Figure`
"""
function plot_eros(sp::AbstractArray, sf::AbstractVector, st::AbstractVector; db::Bool=true, frq::Symbol=:lin, frq_lim::Tuple{Real, Real}=(sf[1], sf[end]), tm::Union{Int64, Vector{Int64}}=0, xlabel::String="default", ylabel::String="default", title::String="default", cb::Bool=true, mono::Bool=false, units::String="μV", smooth::Bool=false, n::Int64=3, kwargs...)::GLMakie.Figure

    @assert size(sp, 1) == length(sf) "Length of sf ($(length(sf))) and number of spectrogram rows ($(size(sp, 1))) must be equal."
    @assert size(sp, 2) == length(st) "Length of st ($(length(st))) and number of spectrogram columns ($(size(sp, 2))) must be equal."
    @assert ndims(sp) == 3 "sp must have 3 dimensions."
    @assert size(sp, 3) <= 2 "sp must contain ≤ 2 epochs."
    @assert n > 0 "n must be ≥ 1."

    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(frq_lim, "frq_lim")

    pal = mono ? :grays : :darktest
    cb_title = db ? "[dB $units^2/Hz]" : "[$units^2/Hz]"

    if frq === :lin
        if frq_lim[2] > 100
            yt = frq_lim[1]:10:frq_lim[2]
        else
            yt = frq_lim[1]:5:frq_lim[2]
        end
    else
        if frq_lim[1] == 0
            _warn("Lower frequency bound truncated to $(sf[2]) Hz")
            frq_lim = (sf[2], frq_lim[2])
        end
        yt = round.(log10space(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1)
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
        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [ms]", "Frequency [Hz]", "Averaged spectrograms of epochs")

        # prepare plot
        plot_size = (1200, 800)
        p = GLMakie.Figure(size=plot_size)
        ax = GLMakie.Axis(p[1, 1],
                          xlabel=xl,
                          ylabel=yl,
                          title=tt,
                          xticks=NeuroAnalyzer._ticks(st),
                          yticks=yt,
                          xminorticksvisible=true,
                          xminorticks=IntervalsBetween(10),
                          yscale=frq===:lin ? identity : log10,
                          xautolimitmargin=(0, 0),
                          yautolimitmargin=(0, 0);
                          kwargs...)
        GLMakie.ylims!(ax, frq_lim)
        ax.titlesize = 20
        ax.xlabelsize = 18
        ax.ylabelsize = 18
        ax.xticklabelsize = 12
        ax.yticklabelsize = 12

        hm = GLMakie.heatmap!(ax,
                              st,
                              sf,
                              sp[:, :, 1]',
                              colormap=pal)
        if cb
            Colorbar(p[1, 2],
                     hm,
                     label=cb_title,
                     labelsize=18)
        end

        # draw time markers
        if tm != 0
            for tm_idx in eachindex(tm)
                GLMakie.vlines!(p[1, 1],
                                [st[tm[tm_idx]]],
                                color=:black,
                                linewidth=1)
            end
        end
    else
        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [ms]", "Frequency [Hz]", "ERP spectrogram")

        # prepare plot
        plot_size = (1200, 800)
        p = GLMakie.Figure(size=plot_size)
        ax1 = GLMakie.Axis(p[1, 1],
                           xlabel=xl,
                           ylabel=yl,
                           title=tt,
                           xticks=_ticks(st),
                           yticks=yt,
                           xminorticksvisible=true,
                           xminorticks=IntervalsBetween(10),
                           yscale=frq===:lin ? identity : log10,
                           xautolimitmargin=(0, 0),
                           yautolimitmargin=(0, 0);
                           kwargs...)
        GLMakie.ylims!(ax1, frq_lim)
        ax1.titlesize = 20
        ax1.xlabelsize = 18
        ax1.ylabelsize = 18
        ax1.xticklabelsize = 12
        ax1.yticklabelsize = 12

        hm1 = GLMakie.heatmap!(ax1,
                               st,
                               sf,
                               sp[:, :, 1]',
                               colormap=pal)

        if cb
            Colorbar(p[1, 2],
                     hm1,
                     label=cb_title,
                     labelsize=18)
        end

        # plot 0 v-line
        GLMakie.vlines!(ax1,
                         [0],
                         linestyle=:dash,
                         linewidth=0.5,
                         color=:black)

        # draw time markers
        if tm != 0
            for tm_idx in eachindex(tm)
                GLMakie.vlines!(p[1, 1],
                                [st[tm[tm_idx]]],
                                color=:black,
                                linewidth=1)
            end
        end

        xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Time [ms]", "Frequency [Hz]", "Averaged spectrograms of ERP epochs")
        ax2 = GLMakie.Axis(p[2, 1],
                           xlabel=xl,
                           ylabel=yl,
                           title=tt,
                           xticks=_ticks(st),
                           yticks=yt,
                           xminorticksvisible=true,
                           xminorticks=IntervalsBetween(10),
                           yscale=frq===:lin ? identity : log10,
                           xautolimitmargin=(0, 0),
                           yautolimitmargin=(0, 0);
                           kwargs...)
        GLMakie.ylims!(ax2, frq_lim)
        ax2.titlesize = 20
        ax2.xlabelsize = 18
        ax2.ylabelsize = 18
        ax2.xticklabelsize = 12
        ax2.yticklabelsize = 12

        hm2 = GLMakie.heatmap!(ax2,
                               st,
                               sf,
                               sp[:, :, 2]',
                               colormap=pal)

        if cb
            Colorbar(p[2, 2],
                     hm2,
                     label=cb_title,
                     labelsize=18)
        end

        # plot 0 v-line
        GLMakie.vlines!(ax2,
                        [0],
                        linestyle=:dash,
                        linewidth=0.5,
                        color=:black)

        # draw time markers
        if tm != 0
            for tm_idx in eachindex(tm)
                GLMakie.vlines!(p[2, 1],
                                [st[tm[tm_idx]]],
                                color=:black,
                                linewidth=1)
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
- `frq_lim::Tuple{Real, Real}=(f[1], f[end])`: frequency limit for the Y axis
- `frq::Symbol=:lin`: linear (`:lin`) or logarithmic (`:log`) frequencies scaling
- `units::String="μV"`
- `mono::Bool=false`: use color or gray palette
- `kwargs`: optional arguments for plotting

# Returns

- `p::GLMakie.Figure`
"""
function plot_erop(sp::AbstractArray, sf::AbstractVector; db::Bool=true, xlabel::String="default", ylabel::String="default", title::String="default", frq_lim::Tuple{Real, Real}=(sf[1], sf[end]), frq::Symbol=:lin, units::String="μV", mono::Bool=false, kwargs...)::GLMakie.Figure

    _in(frq_lim[1], (sf[1], sf[end]), "frq_lim")
    _in(frq_lim[2], (sf[1], sf[end]), "frq_lim")
    @assert size(sp, 1) == length(sf) "Length of sf ($(length(sf))) and number of powers rows ($(size(sp, 1)))) must be equal."
    @assert ndims(sp) == 2 "sp must have 2 dimensions."
    @assert size(sp, 2) <= 2 "sp must contain ≤ 2 epochs."

    _check_var(frq, [:lin, :log], "frq")
    _check_tuple(frq_lim, "frq_lim")

    pal = mono ? :grays : :darktest

    if frq === :lin
        if frq_lim[2] > 100
            xt = frq_lim[1]:10:frq_lim[2]
        else
            xt = frq_lim[1]:5:frq_lim[2]
        end
    else
        if frq_lim[1] == 0
            _warn("Lower frequency bound truncated to $(sf[2]) Hz")
            frq_lim = (sf[2], frq_lim[2])
        end
        xt = round.(log10space(log10(frq_lim[1]), log10(frq_lim[2]), 10), digits=1)
    end

    if size(sp, 2) == 1
        if db
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Frequency [Hz]", "Power [dB $units^2/Hz]", "Averaged power-spectra of epochs")
        else
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Frequency [Hz]", "Power [$units^2/Hz]", "Averaged power-spectra of epochs")
        end

        # prepare plot
        plot_size = (1200, 600)
        p = GLMakie.Figure(size=plot_size)
        ax = GLMakie.Axis(p[1, 1],
                          xlabel=xl,
                          ylabel=yl,
                          title=tt,
                          xticks=xt,
                          xminorticksvisible=true,
                          xminorticks=IntervalsBetween(10),
                          xscale=frq===:lin ? identity : log10,
                          xautolimitmargin=(0, 0),
                          yautolimitmargin=(0, 0);
                          kwargs...)
        GLMakie.xlims!(ax, frq_lim)
        ax.titlesize = 20
        ax.xlabelsize = 18
        ax.ylabelsize = 18
        ax.xticklabelsize = 12
        ax.yticklabelsize = 12

        # plot powers
        Makie.lines!(sf,
                     sp[:, 1],
                     color=:black)

    else
        if db
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Frequency [Hz]", "Power [dB $units^2/Hz]", "ERP power-spectrum")
        else
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Frequency [Hz]", "Power [$units^2/Hz]", "ERP power-spectrum")
        end

        # prepare plot
        plot_size = (1200, 600)
        p = GLMakie.Figure(size=plot_size)
        ax1 = GLMakie.Axis(p[1, 1],
                           xlabel=xl,
                           ylabel=yl,
                           title=tt,
                           xticks=xt,
                           xminorticksvisible=true,
                           xminorticks=IntervalsBetween(10),
                           xscale=frq===:lin ? identity : log10,
                           xautolimitmargin=(0, 0),
                           yautolimitmargin=(0, 0);
                           kwargs...)
        GLMakie.xlims!(ax1, frq_lim)
        ax1.titlesize = 20
        ax1.xlabelsize = 18
        ax1.ylabelsize = 18
        ax1.xticklabelsize = 12
        ax1.yticklabelsize = 12

        # plot powers
        Makie.lines!(ax1,
                     sf,
                     sp[:, 1],
                     color=:black)

        if db
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Frequency [Hz]", "Power [dB $units^2/Hz]", "Averaged power-spectra of epochs")
        else
            xl, yl, tt = _set_defaults(xlabel, ylabel, title, "Frequency [Hz]", "Power [$units^2/Hz]", "Averaged power-spectra of epochs")
        end

        ax2 = GLMakie.Axis(p[2, 1],
                           xlabel=xl,
                           ylabel=yl,
                           title=tt,
                           xticks=xt,
                           xminorticksvisible=true,
                           xminorticks=IntervalsBetween(10),
                           xscale=frq===:lin ? identity : log10,
                           xautolimitmargin=(0, 0),
                           yautolimitmargin=(0, 0);
                           kwargs...)
        GLMakie.xlims!(ax2, frq_lim)
        ax2.titlesize = 20
        ax2.xlabelsize = 18
        ax2.ylabelsize = 18
        ax2.xticklabelsize = 12
        ax2.yticklabelsize = 12

        # plot powers
        Makie.lines!(ax2,
                     sf,
                     sp[:, 2],
                     color=:black)
    end

    return p

end

"""
    plot_icatopo(obj; <keyword arguments>)

Topographical plot of embedded ICA components.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `ic_idx::Union{Int64, Vector{Int64}, AbstractRange}=0`: component(s) to plot, default is all components
- `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
- `cb::Bool=false`: plot color bar
- `cb_label::String="[A.U.]"`: color bar label
- `amethod::Symbol=:mean`: averaging method:
    - `:mean`
    - `:median`
- `imethod::Symbol=:sh`: interpolation method:
    - `:sh`: Shepard
    - `:mq`: Multiquadratic
    - `:imq`: InverseMultiquadratic
    - `:tp`: ThinPlate
    - `:nn`: NearestNeighbour
    - `:ga`: Gaussian
- `nmethod::Symbol=:minmax`: method for normalization, see `normalize()`
- `plot_contours::Bools=true`: plot contours over topo plot
- `plot_electrodes::Bools=true`: plot electrodes over topo plot
- `kwargs`: optional arguments for plotting

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_icatopo(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, ic_idx::Union{Int64, Vector{Int64}, AbstractRange}=0, seg::Tuple{Real, Real}=(0, 10), cb::Bool=false, cb_label::String="default", amethod::Symbol=:mean, imethod::Symbol=:sh, nmethod::Symbol=:minmax, plot_contours::Bool=true, plot_electrodes::Bool=true, kwargs...)::Plots.Plot{Plots.GRBackend}

    @assert :ic in keys(obj.components) "OBJ does not contain :ic component. Perform ica_decompose() first."
    @assert :ic_mw in keys(obj.components) "OBJ does not contain :ic_mw component. Perform ica_decompose() first."

    ic = obj.components[:ic]
    ic_mw = obj.components[:ic_mw]

    # select component channels, default is all channels
    ic_idx == 0 && (ic_idx = _select_cidx(ic, ic_idx))
    _check_cidx(ic, ic_idx)

    p_topo = Vector{Plots.Plot{Plots.GRBackend}}()
    for idx in eachindex(ic_idx)
        obj_tmp = ica_reconstruct(obj, ic, ic_mw, ch=ch, ic_idx=ic_idx[idx], keep=true)
        loc_ch = get_channel(obj_tmp, ch=ch)
        loc_idx = _loc_idx(obj_tmp, loc_ch)
        p_tmp = plot_topo(obj_tmp, ch=loc_idx, title="IC $(ic_idx[idx])", cb=cb, cb_label=cb_label, amethod=amethod, imethod=imethod, nmethod=nmethod, plot_contours=plot_contours, plot_electrodes=plot_electrodes, seg=seg, kwargs...)
        push!(p_topo, p_tmp)
    end

    if length(ic_idx) <= 4
        p = plot_compose(p_topo, layout=(1, 4))
    elseif length(ic_idx) <= 8
        p = plot_compose(p_topo, layout=(2, ceil(Int64, length(ic_idx) / 2)))
    elseif length(ic_idx) <= 12
        p = plot_compose(p_topo, layout=(3, ceil(Int64, length(ic_idx) / 3)))
    else
        p = plot_compose(p_topo, layout=(4, ceil(Int64, length(ic_idx) / 4)))
    end

    return p

end

"""
    plot_icatopo(obj, ic, ic_mw; <keyword arguments>)

Topographical plot of external ICA components.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: NeuroAnalyzer NEURO object
- `ic::Matrix{Float64}`: components IC(1)..IC(n)
- `ic_mw::Matrix{Float64}`: weighting matrix IC(1)..IC(n)
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `ic_idx::Union{Int64, Vector{Int64}, AbstractRange}=0`: component(s) to plot, default is all components
- `seg::Tuple{Real, Real}=(0, 10)`: segment (from, to) in seconds to display, default is 10 seconds or less if single epoch is shorter
- `cb::Bool=false`: plot color bar
- `cb_label::String="[A.U.]"`: color bar label
- `amethod::Symbol=:mean`: averaging method:
    - `:mean`
    - `:median`
- `imethod::Symbol=:sh`: interpolation method:
    - `:sh`: Shepard
    - `:mq`: Multiquadratic
    - `:imq`: InverseMultiquadratic
    - `:tp`: ThinPlate
    - `:nn`: NearestNeighbour
    - `:ga`: Gaussian
- `nmethod::Symbol=:minmax`: method for normalization, see `normalize()`
- `plot_contours::Bools=true`: plot contours over topo plot
- `plot_electrodes::Bools=true`: plot electrodes over topo plot
- `kwargs`: optional arguments for plotting

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_icatopo(obj::NeuroAnalyzer.NEURO, ic::Matrix{Float64}, ic_mw::Matrix{Float64}; ch::Union{String, Vector{String}, Regex}, ic_idx::Union{Int64, Vector{Int64}, AbstractRange}=0, seg::Tuple{Real, Real}=(0, 10), cb::Bool=false, cb_label::String="default", amethod::Symbol=:mean, imethod::Symbol=:sh, nmethod::Symbol=:minmax, plot_contours::Bool=true, plot_electrodes::Bool=true, kwargs...)::Plots.Plot{Plots.GRBackend}

    # select component channels, default is all channels
    ic_idx == 0 && (ic_idx = _select_cidx(ic, ic_idx))
    _check_cidx(ic, ic_idx)

    p_topo = Vector{Plots.Plot{Plots.GRBackend}}()
    for idx in eachindex(ic_idx)
        obj_tmp = ica_reconstruct(obj, ic, ic_mw, ch=ch, ic_idx=ic_idx[idx], keep=true)
        p_tmp = plot_topo(obj_tmp, ch=ch, title="IC $(ic_idx[idx])", cb=cb, cb_label=cb_label, amethod=amethod, imethod=imethod, nmethod=nmethod, plot_contours=plot_contours, plot_electrodes=plot_electrodes, seg=seg, kwargs...)
        push!(p_topo, p_tmp)
    end

    p = Plots.plot(p_topo...)

#=
    if length(ic_idx) <= 4
        p = plot_compose(p_topo, layout=(1, 4))
    elseif length(ic_idx) <= 8
        p = plot_compose(p_topo, layout=(2, ceil(Int64, length(ic_idx) / 2)))
    elseif length(ic_idx) <= 12
        p = plot_compose(p_topo, layout=(3, ceil(Int64, length(ic_idx) / 3)))
    else
        p = plot_compose(p_topo, layout=(4, ceil(Int64, length(ic_idx) / 4)))
    end
=#

    return p

end

"""
    plot_ci(s, s_ci_l, s_ci_h; <keyword arguments>)

Confidence interval plot.

# Arguments

- `s::AbstractVector`: signal
- `s_l::AbstractVector`: CI lower bound
- `s_u::AbstractVector`: CI upper bound
- `t::AbstractVector`: time points
- `xlabel::String=""`: X axis label
- `ylabel::String=""`: Y axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `kwargs`: optional arguments for plotting

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_ci(s::AbstractVector, s_l::AbstractVector, s_u::AbstractVector, t::AbstractVector; xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, kwargs...)::Plots.Plot{Plots.GRBackend}

    @assert length(s) == length(s_l) == length(s_u) "All input signals must be of the same length."

    pal = mono ? :grays : :darktest

    ylim = (floor(minimum(s_l), digits=0), ceil(maximum(s_u), digits=0))
    ylim = _tuple_max(ylim)
    yticks = [ylim[1], 0, ylim[2]]

    # prepare plot
    p = Plots.plot(xlabel=xlabel,
                   ylabel=ylabel,
                   xlims=_xlims(t),
                   xticks=_ticks(t),
                   ylims=ylim,
                   yticks=yticks,
                   ytick_direction=:out,
                   xtick_direction=:out,
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
    # plot upper bound
    Plots.plot!(t,
                s_u,
                fillrange=s_l,
                fillalpha=0.35,
                label=false,
                t=:line,
                c=:grey,
                lw=0.5)
    # plot lower bound
    Plots.plot!(t,
                s_l,
                label=false,
                t=:line,
                c=:grey,
                lw=0.5)
    # plot signal
    Plots.plot!(t,
                s,
                label=false,
                t=:line,
                c=:black,
                lw=0.5)

    return p

end

"""
    plot_heatmap(m; <keyword arguments>)

Plot heatmap.

# Arguments

- `m::AbstractMatrix`
- `x::AbstractVector`
- `y::AbstractVector`
- `xlabel::String=""`: X axis label
- `ylabel::String=""`: Y axis label
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `cb::Bool=true`: draw color bar
- `cb_title::String=""`: color bar title
- `threshold::Union{Nothing, Real}=nothing`: if set, use threshold to mark a region
- `threshold_type::Symbol=:neq`: rule for thresholding:
    - `:eq`: draw region is values are equal to threshold
    - `:neq`: draw region is values are not equal to threshold
    - `:geq`: draw region is values are ≥ to threshold
    - `:leq`: draw region is values are ≤ to threshold
    - `:g`: draw region is values are > to threshold
    - `:l`: draw region is values are < to threshold
- `kwargs`: optional arguments for plotting

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_heatmap(m::AbstractMatrix; x::AbstractVector, y::AbstractVector, xlabel::String="", ylabel::String="", title::String="", mono::Bool=false, cb::Bool=true, cb_title::String="", threshold::Union{Nothing, Real}=nothing, threshold_type::Symbol=:neq, kwargs...)::Plots.Plot{Plots.GRBackend}

    @assert size(m, 1) == length(y) "Number of m rows ($(size(m, 1))) and y length ($(length(y))) must be equal."
    @assert size(m, 2) == length(x) "Number of m columns ($(size(m, 2))) and x length ($(length(x))) must be equal."

    pal = mono ? :grays : :bluesreds

    p = Plots.heatmap(x,
                      y,
                      m,
                      size=(1200, 800),
                      margins=20Plots.px,
                      legend=false,
                      xlabel=xlabel,
                      ylabel=ylabel,
                      title=title,
                      seriescolor=pal,
                      xlims=_xlims(x),
                      ylims=_xlims(y),
                      xticks=_ticks(x),
                      ytick_direction=:out,
                      xtick_direction=:out,
                      cb=cb,
                      colorbar_title=cb_title,
                      titlefontsize=8,
                      xlabelfontsize=8,
                      ylabelfontsize=8,
                      colorbar_titlefontsize=8,
                      xtickfontsize=8,
                      ytickfontsize=8;
                      kwargs...)

    if !isnothing(threshold)
        _, bm = seg_extract(m, threshold=threshold, threshold_type=threshold_type)
        reg = ones(size(m)) .* minimum(m)
        reg[bm] .= maximum(m)
        Plots.plot!(x,
                    y,
                    reg,
                    seriestype=:contour,
                    levels=1,
                    linecolor=:black,
                    colorbar_entry=false,
                    linewidth=2;
                    kwargs...)
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
- `mono::Bool=false`: use color or gray palette
- `kwargs`: optional arguments for plotting

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_imf(imf::Matrix{Float64}; n::Int64=size(imf, 1) - 1, t::AbstractVector, mono::Bool=false, kwargs...)::Plots.Plot{Plots.GRBackend}

    @assert n > 0 "n must be ≥ 1."
    @assert n + 1 <= size(imf, 1) "n must be ≤ $(size(imf, 1) - 1)."
    @assert size(imf, 2) == length(t) "Length of t $(size(imf, 2)) and number of imf columns ($(size(m, 2))) must be equal."

    pal = mono ? :grays : :darktest

    s_restored = sum(imf, dims=1)[:]
    imf = vcat(imf, s_restored')

    ylim = (floor(minimum(imf), digits=0), ceil(maximum(imf), digits=0))
    ylim = _tuple_max(ylim)
    yticks = [ylim[1], 0, ylim[2]]

    # prepare plots
    p_imf = Vector{Plots.Plot{Plots.GRBackend}}()
    for idx in 1:n
        p = Plots.plot(t,
                       imf[idx, :],
                       xlabel="time [s]",
                       ylabel="",
                       xlims=_xlims(t),
                       xticks=_ticks(t),
                       ylims=ylim,
                       yticks=yticks,
                       ytick_direction=:out,
                       xtick_direction=:out,
                       title="IMF: $idx",
                       palette=pal,
                       size=(500, 250),
                       margins=10Plots.px,
                       label=false,
                       line_width=0.5,
                       titlefontsize=6,
                       xlabelfontsize=6,
                       ylabelfontsize=6,
                       xtickfontsize=4,
                       ytickfontsize=4;
                       kwargs...)
        push!(p_imf, p)
    end
    p = Plots.plot(t,
                   imf[end - 1, :],
                   xlabel="time [s]",
                   ylabel="",
                   xlims=_xlims(t),
                   xticks=_ticks(t),
                   ylims=ylim,
                   yticks=yticks,
                   ytick_direction=:out,
                   xtick_direction=:out,
                   title="Residual",
                   palette=pal,
                   size=(500, 250),
                   margins=10Plots.px,
                   label=false,
                   line_width=0.5,
                   titlefontsize=6,
                   xlabelfontsize=6,
                   ylabelfontsize=6,
                   xtickfontsize=4,
                   ytickfontsize=4;
                   kwargs...)
    push!(p_imf, p)
    p = Plots.plot(t,
                   s_restored,
                   xlabel="time [s]",
                   ylabel="",
                   xlims=_xlims(t),
                   xticks=_ticks(t),
                   ylims=ylim,
                   yticks=yticks,
                   ytick_direction=:out,
                   xtick_direction=:out,
                   title="Reconstruced signal",
                   palette=pal,
                   size=(500, 250),
                   margins=10Plots.px,
                   label=false,
                   line_width=0.5,
                   titlefontsize=6,
                   xlabelfontsize=6,
                   ylabelfontsize=6,
                   xtickfontsize=4,
                   ytickfontsize=4;
                   kwargs...)
    push!(p_imf, p)
    mod(n + 2, 2) != 0 && push!(p_imf, plot_empty())

    p = plot_compose(p_imf, layout=(ceil(Int64, (n + 2) / 2), 2))

    return p

end

"""
    plot_hs(sp, st; <keyword arguments>)

Plot Hilbert spectrum.

# Arguments

- `sp::Vector{Float64}`: Hilbert transform powers
- `st::Vector{Float64}`: time
- `xlabel::String="default"`: X axis label, default is Time [s]
- `ylabel::String="default"`: Y axis label, default is Power [μV^2/Hz]
- `title::String="default"`: plot title
- `mono::Bool=false`: use color or gray palette
- `kwargs`: optional arguments for plotting

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_hs(sp::Vector{Float64}, st::Vector{Float64}; xlabel::String="default", ylabel::String="default", title::String="default", mono::Bool=false, kwargs...)::Plots.Plot{Plots.GRBackend}

    @assert length(sp) == length(st) "Length of powers ($(length(sp))) and time points ($(length(st))) must be equal."

    pal = mono ? :grays : :darktest

    xl, yl, tt = _set_defaults(xlabel,
                               ylabel,
                               title,
                               "Time [s]",
                               "Power [μV^2/Hz]",
                               "")

    # prepare plot
    p = Plots.plot(xlabel=xl,
                   ylabel=yl,
                   legend=false,
                   xlims=_xlims(st),
                   xticks=_ticks(st),
                   ytick_direction=:out,
                   xtick_direction=:out,
                   title=tt,
                   palette=pal,
                   t=:line,
                   c=:black,
                   size=(1200, 500),
                   margins=20Plots.px,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=6,
                   ytickfontsize=6;
                   kwargs...)

    # plot powers
    Plots.plot!(st,
                sp,
                linewidth=1,
                label="",
                color=:black)

    return p

end

"""
    plot_fi(fi, st; <keyword arguments>)

Plot instantaneous frequencies.

# Arguments

- `fi::Vector{Float64}`: instantaneous frequencies
- `st::Vector{Float64}`: time
- `xlabel::String="default"`: X axis label, default is Time [s]
- `ylabel::String="default"`: Y axis label, default is Power [μV^2/Hz]
- `title::String="default"`: plot title
- `mono::Bool=false`: use color or gray palette
- `kwargs`: optional arguments for plotting

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_fi(fi::Vector{Float64}, st::Vector{Float64}; xlabel::String="default", ylabel::String="default", title::String="default", mono::Bool=false, kwargs...)::Plots.Plot{Plots.GRBackend}

    @assert length(fi) == length(st) "Length of frequencies ($(length(fi))) and time points ($(length(st))) must be equal."

    pal = mono ? :grays : :darktest

    xl, yl, tt = _set_defaults(xlabel,
                               ylabel,
                               title,
                               "Time [s]",
                               "Frequency [Hz]",
                               "")

    # prepare plot
    p = Plots.plot(xlabel=xl,
                   ylabel=yl,
                   legend=false,
                   xlims=_xlims(st),
                   xticks=_ticks(st),
                   ytick_direction=:out,
                   xtick_direction=:out,
                   title=tt,
                   palette=pal,
                   t=:line,
                   c=:black,
                   size=(1200, 500),
                   margins=20Plots.px,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=6,
                   ytickfontsize=6;
                   kwargs...)

    # plot frequencies
    Plots.plot!(st,
                fi,
                linewidth=1,
                label="",
                color=:black)

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
- `mono::Bool=false`: use color or gray palette
- `kwargs`: optional arguments for plotting

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_phase(ph::Vector{Float64}, sf::Vector{Float64}; unit::Symbol=:rad, type::Symbol=:line, xlabel::String="default", ylabel::String="default", title::String="default", mono::Bool=false, kwargs...)::Plots.Plot{Plots.GRBackend}

    _check_var(unit, [:rad, :deg], "unit")
    _check_var(type, [:line, :stem], "type")
    @assert length(ph) == length(sf) "Length of phases ($(length(fi))) and frequencies ($(length(st))) must be equal."

    pal = mono ? :grays : :darktest

    xl, yl, tt = _set_defaults(xlabel,
                               ylabel,
                               title,
                               "Frequency [Hz]",
                               unit === :rad ? "Phase [rad]" : "Phase [°]",
                               "")

    # prepare plot
    p = Plots.plot(xlabel=xl,
                   ylabel=yl,
                   legend=false,
                   xlims=_xlims(sf),
                   xticks=_ticks(sf),
                   ytick_direction=:out,
                   xtick_direction=:out,
                   title=tt,
                   palette=pal,
                   t=:line,
                   c=:black,
                   size=(1200, 500),
                   margins=20Plots.px,
                   titlefontsize=8,
                   xlabelfontsize=8,
                   ylabelfontsize=8,
                   xtickfontsize=6,
                   ytickfontsize=6;
                   kwargs...)

    # plot phases
    if type === :line
        Plots.plot!(sf,
                    ph,
                    linewidth=1,
                    label="",
                    color=:black)
    else
        Plots.plot!(sf,
                    ph,
                    seriestype=:stem,
                    linewidth=1,
                    label="",
                    color=:black)
        Plots.scatter!(sf,
                       ph,
                       linewidth=1,
                       label="",
                       mc=:black,
                       ms=2.0)
    end

    return p

end

"""
    plot_polezero(s; <keyword arguments>)

Polar pole-zero map.

# Arguments

- `p::Vector{Complex{Float64}}`: vector of poles
- `z::Vector{Complex{Float64}}`: vector of zeros
- `m::Tuple{Real, Real}=(0, 0)`: major value to plot
- `title::String=""`: plot title
- `mono::Bool=false`: use color or gray palette
- `ticks::Bool=false`: draw X and Y ticks
- `ms::Symbol=:circle`: marker shape for drawing complex numbers (`:circle` or `:xcross`)
- `kwargs`: optional arguments for plotting

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_polezero(pol::Vector{Complex{Float64}}, zer::Vector{Complex{Float64}}; title::String="default", mono::Bool=false, kwargs...)::Plots.Plot{Plots.GRBackend}

    pal = mono ? :grays : :darktest

    p = Plots.plot(size=(600, 600),
                   aspect_ratio=:equal,
                   framestyle=:origin,
                   xlims=(-1.1, 1.1),
                   ylims=(-1.1, 1.1),
                   left_margin=-20Plots.px,
                   right_margin=0Plots.px,
                   bottom_margin=0Plots.px,
                   title=title == "default" ? "Pole-zero map" : title,
                   xtitle="Real",
                   ytitle="Imag",
                   ytick_direction=:out,
                   xtick_direction=:out,
                   titlefontsize=8,
                   xtickfontsize=7,
                   ytickfontsize=7,
                   palette=pal,
                   legend=false;
                   kwargs...)
    Plots.scatter!([real.(pol)], [imag.(pol)],
                   markersize=4,
                   markercolor=mono ? :black : :blue,
                   markershape=:xcross)
    Plots.scatter!([real.(zer)], [imag.(zer)],
                   markersize=5,
                   markerstrokecolor=mono ? :black : :blue,
                   markercolor=:white,
                   markershape=:circle)

    Plots.plot!(Plots.partialcircle(0, 2π, 64, 1.0),
                ls=:dot,
                lw=0.5,
                lc=:black)

    return p

end


"""
    plot_dwc(dc; <keyword arguments>)

Plot discrete wavelet decomposition coefficients.

# Arguments

- `dc::Matrix{Float64}`: coefficients
- `n::Int64=size(dc, 1) - 1`: number of coefficients to plot
- `t::AbstractVector`: time points
- `mono::Bool=false`: use color or gray palette
- `kwargs`: optional arguments for plotting

# Returns

- `p::Plots.Plot{Plots.GRBackend}`
"""
function plot_dwc(dc::Matrix{Float64}; n::Int64=size(dc, 1) - 1, t::AbstractVector, mono::Bool=false, kwargs...)::Plots.Plot{Plots.GRBackend}

    @assert n > 1 "n must be > 1."
    @assert n <= size(dc, 1) - 1 "n must be ≤ $(size(dc, 1) - 1)."
    @assert size(dc, 2) == length(t) "Length of t $(size(dc, 2)) and number of dc columns ($(size(m, 2))) must be equal."

    pal = mono ? :grays : :darktest

    ylim = (floor(minimum(dc), digits=0), ceil(maximum(dc), digits=0))
    ylim = _tuple_max(ylim)
    yticks = [ylim[1], 0, ylim[2]]

    # prepare plots
    p_dc = Vector{Plots.Plot{Plots.GRBackend}}()
    for idx in 2:(n + 1)
        p = Plots.plot(t,
                       dc[idx, :],
                       xlabel="time [s]",
                       ylabel="",
                       xlims=_xlims(t),
                       xticks=_ticks(t),
                       ylims=ylim,
                       yticks=yticks,
                       ytick_direction=:out,
                       xtick_direction=:out,
                       title="Coefficient #$(idx - 1)",
                       palette=pal,
                       size=(500, 250),
                       margins=10Plots.px,
                       label=false,
                       line_width=0.5,
                       titlefontsize=6,
                       xlabelfontsize=6,
                       ylabelfontsize=6,
                       xtickfontsize=4,
                       ytickfontsize=4;
                       kwargs...)
        push!(p_dc, p)
    end

    p = Plots.plot(t,
                   dc[1, :],
                   xlabel="time [s]",
                   ylabel="",
                   xlims=_xlims(t),
                   xticks=_ticks(t),
                   ylims=ylim,
                   yticks=yticks,
                   ytick_direction=:out,
                   xtick_direction=:out,
                   title="Original signal",
                   palette=pal,
                   size=(500, 250),
                   margins=10Plots.px,
                   label=false,
                   line_width=0.5,
                   titlefontsize=6,
                   xlabelfontsize=6,
                   ylabelfontsize=6,
                   xtickfontsize=4,
                   ytickfontsize=4;
                   kwargs...)
    push!(p_dc, p)

    mod(n + 2, 2) != 0 && push!(p_dc, plot_empty())

    p = plot_compose(p_dc, layout=(ceil(Int64, (n + 2) / 2), 2))

    return p

end