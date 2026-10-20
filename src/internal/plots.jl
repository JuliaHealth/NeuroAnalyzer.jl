"""
    _draw_head_labels!(ax; <keyword arguments>)

Draw fiducial point labels (NAS, IN, LPA, RPA) onto a 2-D axis.

# Arguments

- `ax::GLMakie.Axis`: target axis
- `plane::Symbol`: which anatomical plane is displayed (`:xy`, `:xz`, or `:yz`); determines which coordinates of each fiducial point are used
- `font_size::Int64=8`: label font size in points
"""
function _draw_head_labels!(
    ax::GLMakie.Axis;
    plane::Symbol = :xy,
    font_size::Int64 = 8,
)::Nothing
    _check_var(plane, [:xy, :xz, :yz], "plane")
    fid_names = ["NAS", "IN", "LPA", "RPA"]
    # (eachindex of a NamedTuple gives its keys, not positions)
    for (idx, pt) in enumerate(values(NeuroAnalyzer.fiducial_points))
        fid_loc_x, fid_loc_y = if plane === :xy
            pt[1], pt[2]
        elseif plane === :xz
            pt[1], pt[3]
        elseif plane === :yz
            pt[2], pt[3]
        end
        GLMakie.text!(
            ax,
            fid_loc_x, fid_loc_y;
            text     = fid_names[idx],
            fontsize = font_size,
            align    = (:center, :center),
        )
    end
    return nothing
end

"""
    _smooth_path(x, y; <keyword arguments>)

Interpolate a polyline with a Catmull-Rom spline (passes through every input point).

# Arguments

- `x::Vector{<:Real}`: x coordinates
- `y::Vector{<:Real}`: y coordinates
- `n::Int64=8`: number of interpolated points per segment

# Returns

- `Vector{Point2f}`
"""
function _smooth_path(x::Vector{<:Real}, y::Vector{<:Real}; n::Int64 = 8)
    # Catmull-Rom on a single coordinate vector
    function cr(v::Vector{Float32})
        w = [v[1]; v; v[end]]  # duplicate endpoints
        out = Float32[]
        for k in 1:(length(v) - 1)
            v0, v1, v2, v3 = w[k], w[k + 1], w[k + 2], w[k + 3]
            for t in range(0.0f0, 1.0f0; length = n + 1)[1:n]
                push!(
                    out,
                    0.5f0 * (
                        2.0f0 * v1 +
                        (v2 - v0) * t +
                        (2.0f0 * v0 - 5.0f0 * v1 + 4.0f0 * v2 - v3) * t^2 +
                        (3.0f0 * v1 - v0 - 3.0f0 * v2 + v3) * t^3
                    ),
                )
            end
        end
        push!(out, v[end])
        return out
    end
    return Point2f.(cr(Float32.(x)), cr(Float32.(y)))
end

"""
    _draw_head_outline!(ax; <keyword arguments>)

Draw a schematic head outline (circle, nose, ears) onto a 2-D axis.

The head is assumed to have unit radius centered at the origin. Set `ax.aspect = DataAspect()` on the axis, otherwise the head is drawn as an ellipse.

# Arguments

- `ax::GLMakie.Axis`: target axis
- `lw::Real=1`: line width
- `color=:black`: line color
"""
function _draw_head_outline!(ax::GLMakie.Axis; lw::Real = 1, color = :black)::Nothing
    kw = (linewidth = lw, color = color)

    # head circle
    θ = range(0, 2π; length = 181)
    GLMakie.lines!(ax, Point2f.(cos.(θ), sin.(θ)); kw...)

    # nose: one polyline so the tip gets a proper join
    GLMakie.lines!(ax, Point2f[(-0.2, 0.98), (0.0, 1.1), (0.2, 0.98)]; kw...)

    # ears: right ear mirrors the left, so they cannot drift apart
    ear_x = [0.995, 1.03, 1.06, 1.1, 1.12, 1.1, 1.13, 1.09, 1.02, 0.98, 0.975]
    ear_y = [0.1, 0.15, 0.16, 0.14, 0.05, -0.1, -0.3, -0.37, -0.39, -0.33, -0.22]
    GLMakie.lines!(ax, _smooth_path(-ear_x, ear_y); kw...)
    GLMakie.lines!(ax, _smooth_path(ear_x, ear_y); kw...)

    return nothing
end

"""
    _xlims(t)

Return `(floor(t[1], digits=2), ceil(t[end], digits=2))` as axis x-limits.
"""
_xlims(t::Union{AbstractVector, AbstractRange})::Tuple{Real, Real} =
    floor(t[1]; digits = 2), ceil(t[end]; digits = 2)

"""
    _ylims(s)

Return symmetric y-limits `(-m, m)` suitable for displaying signal `s`.

The magnitude `m` is derived from the signal's peak absolute value, rounded to a precision that scales with the signal's dynamic range.
"""
function _ylims(s::Union{AbstractVector, AbstractMatrix})::Tuple{Real, Real}
    peak = maximum(abs, s)

    n = peak > 100 ? 2 : peak >= 10 ? 1 : 0

    hi = ceil(Int64, round(maximum(s); digits = n))
    lo = floor(Int64, round(minimum(s); digits = n))

    # all-zero signal: provide a unit range
    hi == 0 && lo == 0 && return (-1.0, 1.0)

    # one bound is zero: make the range symmetric around zero
    lo == 0 && return (-abs(hi), abs(hi))
    hi == 0 && return (-abs(lo), abs(lo))

    # general case: use the larger absolute bound
    m = max(abs(hi), abs(lo))
    return (-m, m)
end

"""
    _ticks(t)

Return a vector of axis tick positions covering the range of `t`.

The step size is chosen based on the total range:
- ≤ 10: step 1
- ≤ 20: step 2
- ≤ 50: step 5
- > 50: step 10
"""
function _ticks(t::Union{AbstractVector, AbstractRange, Tuple{Real, Real}})::AbstractVector
    t1 = t isa Tuple ? t[1] : t[1]
    t2 = t isa Tuple ? t[2] : t[end]
    rng = t2 - t1
    step = rng > 50 ? 10 : rng > 20 ? 5 : rng > 10 ? 2 : 1
    return collect(t1:step:t2)
end

"""
    _erpticks(t)

Return a tick vector suitable for ERP plots: 3 ticks from `floor(t[1])` to 0, then 8 more ticks from 0 to `ceil(t[end])`.

Accepts either a vector/range or a `Tuple{Real, Real}`.
"""
function _erpticks(
    t::Union{AbstractVector, AbstractRange, Tuple{Real, Real}},
)::AbstractVector
    t1 = t isa Tuple ? t[1] : t[1]
    t2 = t isa Tuple ? t[2] : t[end]
    return vcat(
        collect(range(floor(t1; digits = 2), 0; length = 3)),
        collect(range(0, ceil(t2; digits = 2); length = 9))[2:end],
    )
end

"""
    _set_defaults(xl, yl, tt, x, y, t)

Replace any `"default"` placeholder in `xl`, `yl`, `tt` with the corresponding fallback string (`x`, `y`, `t`).

# Returns

- `Tuple{String, String, String}`: `(xlabel, ylabel, title)`
"""
function _set_defaults(
    xl::String, yl::String, tt::String,
    x::String, y::String, t::String,
)::Tuple{String, String, String}
    xl == "default" && (xl = x)
    yl == "default" && (yl = y)
    tt == "default" && (tt = t)
    return xl, yl, tt
end

"""
    _bernstein(i, n; steps=50)

Return the values of the `i`-th Bernstein basis polynomial of degree `n` evaluated at `steps` evenly spaced points in [0, 1].
"""
_bernstein(i::Int, n::Int; steps::Int = 50) =
    [binomial(n, i) * t^i * (1 - t)^(n - i) for t in LinRange(0, 1, steps)]

"""
    _bernstein_poly(px, py; steps=50)

Evaluate a Bézier curve defined by control points `(px[i], py[i])` at `steps` parameter values, using Bernstein polynomial basis functions.

# Arguments

- `px::AbstractVector`: x-coordinates of control points
- `py::AbstractVector`: y-coordinates of control points
- `steps::Int=50`: number of evaluation points

# Returns

- `Tuple{Vector{Float64}, Vector{Float64}}`: `(x_vals, y_vals)` along the curve

# Reference

Based on https://opensourc.es/blog/bezier-curve/
"""
function _bernstein_poly(
    px::AbstractVector,
    py::AbstractVector;
    steps::Int = 50,
)::Tuple{Vector{Float64}, Vector{Float64}}
    n = length(px) - 1
    length(px) == length(py) ||
        throw(ArgumentError("px and py must have the same length."))

    b = [_bernstein(i, n; steps = steps) for i = 0:n]
    x_vals = [sum(px[k] * b[k][t] for k = 1:(n + 1)) for t = 1:steps]
    y_vals = [sum(py[k] * b[k][t] for k = 1:(n + 1)) for t = 1:steps]
    return x_vals, y_vals
end

# position and moment of a dipole at time point t
function _dipole_at(d::NeuroAnalyzer.DIPOLE, t::Int64)::Tuple{Vector{Float64}, Vector{Float64}}
    1 <= t <= size(d.q, 2) || throw(ArgumentError("t must be in [1, $(size(d.q, 2))]."))
    pos = size(d.pos, 1) == 1 ? d.pos[1, :] : d.pos[t, :]
    return pos, d.q[:, t]
end

"""
    _log_flim(f, flim, frq)

For a logarithmic frequency axis (`frq = :log`) the lower frequency bound must be > 0: a bound ≤ 0 is replaced (with a warning) by the lowest positive frequency in `f`.
"""
function _log_flim(f::AbstractVector, flim::Tuple{Real, Real}, frq::Symbol)::Tuple{Real, Real}
    (frq === :log && flim[1] <= 0) || return flim
    idx = findfirst(>(0), f)
    isnothing(idx) &&
        throw(ArgumentError("Logarithmic frequency axis requires positive frequencies."))
    _warn("Lower frequency bound truncated to $(f[idx]) Hz.")
    return (f[idx], flim[2])
end

"""
    _trace_colors(n; <keyword arguments>)

Colors of the traces of a multi-channel line plot (`plot_psd()`, `plot_phs()`, `plot_erp()`).

# Arguments

- `n::Int64`: number of traces
- `mono::Bool`: monochrome palette
- `black::Bool`: all traces black (single channel, or only the mean ± 95% CI is drawn)
"""
function _trace_colors(n::Int64; mono::Bool, black::Bool)::Vector
    return black ? fill(:black, n) : GLMakie.resample_cmap(mono ? :grays : :darktest, max(2, n))[1:n]
end

"""
    _ch_locs_inset!(fig, obj, ch, colors; <keyword arguments>)

Add a head outline with the locations of the plotted channels (`plot_locs(; ps = :s)`), drawn in the colors of their traces, next to or inside the axis of `fig`. Nothing is drawn (with a warning) if locations of the channels are not available.

# Arguments

- `fig::GLMakie.Figure`: plot (axis at `fig[1, 1]`)
- `obj::NeuroAnalyzer.NEURO`: plotted object
- `ch::Union{Integer, AbstractVector{<:Integer}}`: plotted channel(s) (in the order of the traces)
- `colors::AbstractVector`: trace colors
- `pos::Symbol=:outside`: inset position:
    - `:outside`: to the right of the axis (top-aligned)
    - `:rt`, `:lt`, `:rb`, `:lb`: inside the axis, in the right-top, left-top, right-bottom or left-bottom corner
"""
function _ch_locs_inset!(
    fig::GLMakie.Figure,
    obj::NeuroAnalyzer.NEURO,
    ch::Union{Integer, AbstractVector{<:Integer}},
    colors::AbstractVector;
    pos::Symbol = :outside,
)::Nothing
    _check_var(pos, [:outside, :rt, :lt, :rb, :lb], "pos")
    ch = ch isa Integer ? [ch] : ch
    loc_labels = lowercase.(String.(obj.locs[!, :label]))
    loc_idx = [findfirst(==(lowercase(l)), loc_labels) for l in labels(obj)[ch]]
    if isempty(loc_labels) || any(isnothing, loc_idx)
        _warn("Locations of the plotted channels are not available: head inset is not drawn.")
        return nothing
    end

    pl = plot_locs(
        obj.locs;
        ch = Vector{Int64}(loc_idx),
        ch_colors = colors,
        ps = :s,
        ch_labels = false,
        head_labels = false,
        grid = false,
        gui = false,
    )
    _figure_inset!(fig, pl; pos = pos)

    return nothing
end

"""
    _figure_inset!(fig, pl; <keyword arguments>)

Draw figure `pl` (rendered to an image, white background made transparent) as an inset next to or inside the axis at `fig[1, 1]`.

# Arguments

- `fig::GLMakie.Figure`: main figure (axis at `fig[1, 1]`)
- `pl::GLMakie.Figure`: figure to draw as an inset
- `pos::Symbol=:outside`: inset position:
    - `:outside`: to the right of the axis (top-aligned)
    - `:rt`, `:lt`, `:rb`, `:lb`: inside the axis, in the right-top, left-top, right-bottom or left-bottom corner
- `width::Union{Nothing, Int64}=nothing`: inset width in pixels (aspect ratio is kept); `nothing`: the size of `pl`
"""
function _figure_inset!(
    fig::GLMakie.Figure,
    pl::GLMakie.Figure;
    pos::Symbol = :outside,
    width::Union{Nothing, Int64} = nothing,
)::Nothing
    _check_var(pos, [:outside, :rt, :lt, :rb, :lb], "pos")
    isnothing(width) || width > 0 || throw(ArgumentError("width must be > 0."))

    img = GLMakie.colorbuffer(pl)
    # white background → transparent
    # (colorbuffer gives 8-bit N0f8 colors, whose sums wrap around: convert to Float64 first)
    img = map(img) do c
        r, g, b = Float64(red(c)), Float64(green(c)), Float64(blue(c))
        return RGBA(r, g, b, r + g + b > 2.97 ? 0.0 : 1.0)
    end
    w = isnothing(width) ? size(img, 2) : width
    h = round(Int64, w * size(img, 1) / size(img, 2))

    if pos === :outside
        # to the right of the plot, at the top of the side column (above the legend)
        side = _side_layout(fig)
        ax = GLMakie.Axis(
            side[1, 1];
            width = w,
            height = h,
            tellheight = false,
            valign = :top,
            aspect = DataAspect(),
            backgroundcolor = :transparent,
            _AXIS_LOCK_KWARGS...,
        )
        GLMakie.rowsize!(side, 1, GLMakie.Fixed(h))
        GLMakie.rowgap!(side, _SIDE_ROWGAP)
    else
        # in a corner of the plot area (same cell as the main axis, 10 px from its edges)
        ax = GLMakie.Axis(
            fig[1, 1];
            width = w,
            height = h,
            tellwidth = false,
            tellheight = false,
            halign = pos in (:rt, :rb) ? :right : :left,
            valign = pos in (:rt, :lt) ? :top : :bottom,
            alignmode = Outside(10),
            aspect = DataAspect(),
            backgroundcolor = :transparent,
            _AXIS_LOCK_KWARGS...,
        )
    end
    hidedecorations!(ax)
    hidespines!(ax)
    GLMakie.image!(ax, rotr90(img))

    # less room is left for the legend below an outside inset
    pos === :outside && _fit_side_legend!(fig)

    return nothing
end

# ---------------------------------------------------------------------------
# legend of multi-channel plots: outside the plot, on the right
# ---------------------------------------------------------------------------

# gap between the outside inset and the legend [px]
const _SIDE_ROWGAP = 10

# side column to the right of the main axis (fig[1, 2]): an outside inset (row 1) above
# the legend (row 2); an unused row takes no space
function _side_layout(fig::GLMakie.Figure)::GLMakie.GridLayout
    gl = [c for c in GLMakie.contents(fig[1, 2]) if c isa GLMakie.GridLayout]
    isempty(gl) || return gl[1]
    side = GLMakie.GridLayout()
    fig[1, 2] = side
    # until an inset is added, its row takes no space (the legend starts at the plot's top)
    GLMakie.rowsize!(side, 1, GLMakie.Fixed(0))
    GLMakie.rowgap!(side, 0)
    return side
end

"""
    _side_legend!(fig, ax)

Draw the legend of the labelled plots of `ax` outside the plot, on the right (below the head inset if it is drawn outside). The legend is drawn in 1 or 2 columns so that it is not higher than the plot; if it does not fit (too many channels), it is not drawn.
"""
function _side_legend!(fig::GLMakie.Figure, ax::GLMakie.Axis)::Nothing
    side = _side_layout(fig)
    GLMakie.Legend(
        side[2, 1],
        ax;
        tellheight = false,
        valign = :top,
        labelsize = 12,
        rowgap = 0,
    )
    _fit_side_legend!(fig)
    return nothing
end

# number of legend columns (1 or 2) so that the legend is not higher than the plot,
# minus the outside inset above it; the legend is removed if it does not fit
function _fit_side_legend!(fig::GLMakie.Figure)::Nothing
    legs = [c for c in fig.content if c isa GLMakie.Legend]
    isempty(legs) && return nothing
    leg = legs[1]
    h_av = Float64(_main_axis(fig).scene.viewport[].widths[2])
    side = _side_layout(fig)
    inset = [c for c in GLMakie.contents(side[1, 1]) if c isa GLMakie.Axis]
    isempty(inset) || (h_av -= inset[1].height[] + _SIDE_ROWGAP)
    for nb in (1, 2)
        leg.nbanks = nb
        leg.layoutobservables.autosize[][2] <= h_av && return nothing
    end
    delete!(leg)
    _info("Too many channels to fit the legend: legend is not drawn")
    return nothing
end

# main axis of a plot: the first axis created (an inset may share its layout cell)
_main_axis(fig::GLMakie.Figure)::GLMakie.Axis = first(c for c in fig.content if c isa GLMakie.Axis)
