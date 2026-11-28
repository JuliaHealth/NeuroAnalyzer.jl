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
    # only the points visible in the view (otherwise NAS/IN or LPA/RPA would be drawn at the same position):
    # :xz is seen from behind (no nasion), :yz from the right (no left preauricular point)
    visible =
        plane === :xy ? fid_names :
        plane === :xz ? ["IN", "LPA", "RPA"] : ["NAS", "IN", "RPA"]
    # (eachindex of a NamedTuple gives its keys, not positions)
    for (idx, pt) in enumerate(values(NeuroAnalyzer.fiducial_points))
        fid_names[idx] in visible || continue
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
        for k = 1:(length(v) - 1)
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

Draw a schematic head outline onto a 2-D axis.

The head is assumed to have unit radius centered at the origin. Set `ax.aspect = DataAspect()` on the axis, otherwise the head is drawn as an ellipse.

- `:xy` (top view, nose up): circle, nose and both ears
- `:xz` (coronal view from behind, subject's right on the right): upper head outline and both ears; the nose is not visible from behind
- `:yz` (sagittal view from the right, nose on the right): head circle, nose (bridge, tip and a short part of the underside) and the right ear (traced from a drawing); the nose extends to y ≈ 1.21, so the axis limits must include it

# Arguments

- `ax::GLMakie.Axis`: target axis
- `plane::Symbol=:xy`: `:xy`, `:xz` or `:yz`
- `lw::Real=1`: line width
- `color=:black`: line color
"""
function _draw_head_outline!(
    ax::GLMakie.Axis;
    plane::Symbol = :xy,
    lw::Real = 1,
    color = :black,
)::Nothing
    _check_var(plane, [:xy, :xz, :yz], "plane")
    kw = (linewidth = lw, color = color)

    # ear profile (top view and coronal view): right ear mirrors the left, so they cannot drift apart
    ear_x = [0.995, 1.03, 1.06, 1.1, 1.12, 1.1, 1.13, 1.09, 1.02, 0.98, 0.975]
    ear_y = [0.1, 0.15, 0.16, 0.14, 0.05, -0.1, -0.3, -0.37, -0.39, -0.33, -0.22]

    if plane === :xy
        # head circle
        θ = range(0, 2π; length = 181)
        GLMakie.lines!(ax, Point2f.(cos.(θ), sin.(θ)); kw...)
        # nose: one polyline so the tip gets a proper join
        GLMakie.lines!(ax, Point2f[(-0.2, 0.98), (0.0, 1.1), (0.2, 0.98)]; kw...)
        GLMakie.lines!(ax, _smooth_path(-ear_x, ear_y); kw...)
        GLMakie.lines!(ax, _smooth_path(ear_x, ear_y); kw...)
    elseif plane === :xz
        # upper part of the head, extended below the horizontal plane down to the ears
        θ = range(-0.6, π + 0.6; length = 181)
        GLMakie.lines!(ax, Point2f.(cos.(θ), sin.(θ)); kw...)
        # ears at the preauricular points, z ≈ -0.2
        GLMakie.lines!(ax, _smooth_path(-ear_x, ear_y); kw...)
        GLMakie.lines!(ax, _smooth_path(ear_x, ear_y); kw...)
    else
        # head circle from the nasion (z ≈ -0.19) over the top and down the back of the head, so that it
        # matches the electrode positions on the unit sphere
        θ = range(asin(-0.189), π + 0.6; length = 181)
        GLMakie.lines!(ax, Point2f.(cos.(θ), sin.(θ)); kw...)
        # nose and right ear traced from a drawing of a head profile (scaled to the head model)
        # bridge, rounded tip and a short part of the underside
        nose_y = [0.982, 1.012, 1.076, 1.15, 1.2, 1.208, 1.18, 1.12]
        nose_z = [-0.189, -0.262, -0.345, -0.444, -0.53, -0.565, -0.59, -0.615]
        GLMakie.lines!(ax, _smooth_path(nose_y, nose_z; n = 4); kw...)
        # right ear
        ear_y = [
            0.02,
            -0.043,
            -0.129,
            -0.192,
            -0.225,
            -0.235,
            -0.215,
            -0.17,
            -0.1,
            -0.02,
            0.039,
        ]
        ear_z = [
            -0.051,
            0.007,
            0.025,
            -0.029,
            -0.145,
            -0.291,
            -0.371,
            -0.455,
            -0.545,
            -0.6,
            -0.625,
        ]
        GLMakie.lines!(ax, _smooth_path(ear_y, ear_z; n = 4); kw...)
    end

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
    _cb_ticks(crange; n)

Return colorbar ticks: the minimum, the maximum and `n - 2` equally spaced values in between.

Labels use a fixed number of decimals chosen from the tick step, so that neighbouring labels differ (e.g. 0.23, 0.36, 0.49 for a range of 0.232–0.746); very small or very large values use 3 significant digits.

# Arguments

- `crange::Tuple{Real, Real}`: color range (minimum, maximum)
- `n::Int64=5`: number of ticks (≥ 2)

# Returns

- `Tuple{Vector{Float64}, Vector{String}}`: tick positions and labels
"""
function _cb_ticks(crange::Tuple{Real, Real}; n::Int64 = 5)::Tuple{Vector{Float64}, Vector{String}}
    n < 2 && throw(ArgumentError("n must be ≥ 2."))
    lo, hi = Float64.(crange)
    # constant data: a single tick
    hi == lo && return ([lo], [_cb_label(lo, 3)])
    v = collect(range(lo, hi; length = n))
    step = (hi - lo) / (n - 1)
    m = max(abs(lo), abs(hi))
    lab = if m >= 1e5 || m < 1e-3
        [(r = round(x; sigdigits = 3) + 0.0;
          isinteger(r) && abs(r) < 1e15 ? string(round(Int64, r)) : string(r))
         for x in v]
    else
        d = clamp(ceil(Int64, -log10(step)) + 1, 0, 6)
        [_cb_label(x, d) for x in v]
    end
    # the colorbar limits are compared in Float32: ticks placed exactly at the
    # limits may fall outside them and be dropped, so the end ticks are moved
    # inwards by a negligible amount (labels are not changed)
    v[1] += (hi - lo) * 1e-6
    v[end] -= (hi - lo) * 1e-6
    return (v, lab)
end

"""
    _cb_tickfun(; n)

Return a tick function for `GLMakie.Colorbar(...; ticks = _cb_tickfun())`: Makie calls it with the colorbar limits, so the ticks (see `_cb_ticks()`) always cover the actual color range, also when it is set automatically or changes.

# Arguments

- `n::Int64=5`: number of ticks (≥ 2)

# Returns

- `Function`: `(vmin, vmax) -> (positions, labels)`
"""
_cb_tickfun(; n::Int64 = 5)::Function = (vmin, vmax) -> _cb_ticks((vmin, vmax); n = n)

# format x with exactly d decimals (no Printf dependency); -0.0 is shown as 0
function _cb_label(x::Real, d::Int64)::String
    d == 0 && return string(round(Int64, x))
    s = string(round(x; digits = d) + 0.0)
    occursin('e', s) && return string(round(x; sigdigits = 3) + 0.0)
    i, f = occursin('.', s) ? split(s, '.') : (s, "")
    return i * "." * rpad(f, d, '0')
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
    _bernstein(i, n; <keyword arguments>)

Return the values of the `i`-th Bernstein basis polynomial of degree `n` evaluated at `steps` evenly spaced points in [0, 1].

# Arguments

- `i::Int`: index of the basis polynomial
- `n::Int`: degree of the basis polynomial
- `steps::Int=50`: number of evaluation points

# Returns

- `Vector{Float64}`
"""
_bernstein(i::Int, n::Int; steps::Int = 50) =
    [binomial(n, i) * t^i * (1 - t)^(n - i) for t in LinRange(0, 1, steps)]

"""
    _bernstein_poly(px, py; <keyword arguments>)

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
function _dipole_at(
    d::NeuroAnalyzer.DIPOLE,
    t::Int64,
)::Tuple{Vector{Float64}, Vector{Float64}}
    1 <= t <= size(d.q, 2) || throw(ArgumentError("t must be in [1, $(size(d.q, 2))]."))
    pos = size(d.pos, 1) == 1 ? d.pos[1, :] : d.pos[t, :]
    return pos, d.q[:, t]
end

"""
    _log_flim(f, flim, frq)

For a logarithmic frequency axis (`frq = :log`) the lower frequency bound must be > 0: a bound ≤ 0 is replaced (with a warning) by the lowest positive frequency in `f`.
"""
function _log_flim(
    f::AbstractVector,
    flim::Tuple{Real, Real},
    frq::Symbol,
)::Tuple{Real, Real}
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
function _trace_colors(n::Int64; mono::Bool, black::Bool = false)::Vector
    black && return fill(:black, n)
    # mono: grays from black to mid-gray (never white, which is invisible on the background)
    mono && return [GLMakie.RGBf(g, g, g) for g in range(0, _MONO_MAXGRAY; length = n)]
    return GLMakie.resample_cmap(:darktest, max(2, n))[1:n]
end

# monochrome palettes, used by all plots with `mono = true`:
# - continuous color scales (heatmaps, topographical maps, spectrograms): _MONO_CMAP
# - colors of overlaid traces (e.g. butterfly plots): _trace_colors(n; mono = true),
#   grays from black to _MONO_MAXGRAY
# (traces drawn in separate rows, e.g. signal plots, are black)
const _MONO_CMAP = :grays
const _MONO_MAXGRAY = 0.6

# continuous color scale `name`, or the monochrome color scale
_palette(name::Symbol; mono::Bool) = mono ? _MONO_CMAP : name

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
        _warn(
            "Locations of the plotted channels are not available: head inset is not drawn.",
        )
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

# figure `pl` rendered to an image, with the white background made transparent
function _figure_image(pl::GLMakie.Figure)::Matrix{RGBA{Float64}}
    img = GLMakie.colorbuffer(pl)
    # (colorbuffer gives 8-bit N0f8 colors, whose sums wrap around: convert to Float64 first)
    return map(img) do c
        r, g, b = Float64(red(c)), Float64(green(c)), Float64(blue(c))
        return RGBA(r, g, b, r + g + b > 2.97 ? 0.0 : 1.0)
    end
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

    img = _figure_image(pl)
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

# ------------------------------------------------------------------------------
# legend of multi-channel plots: outside the plot, on the right
# ------------------------------------------------------------------------------

# gap between the outside inset and the legend [px]
const _SIDE_ROWGAP = 10

# side column to the right of the main axis (fig[1, 2]): an outside inset (row 1) above
# the legend (row 2); an unused row takes no space
function _side_layout(fig::GLMakie.Figure)::GLMakie.GridLayout
    gl = [c for c in GLMakie.contents(fig[1, 2]) if c isa GLMakie.GridLayout]
    isempty(gl) || return gl[1]
    # tellheight = false: the side column (inset, legend) must not set the height of the
    # plot row, otherwise a small inset squashes the plot to its own height
    side = GLMakie.GridLayout(; tellheight = false, valign = :top)
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
_main_axis(fig::GLMakie.Figure)::GLMakie.Axis =
    first(c for c in fig.content if c isa GLMakie.Axis)

"""
    _plot_scale(s)

Common amplitude scale of a group of channels of the same type (rows of `s`) for the multi-channel signal plots (one row = the scale), the larger of:

- 6 × the median (over channels) of the robust standard deviation (1.4826 × MAD) of each channel: a typical channel spans about ±3 SD in its row,
- 1.2 × the median (over channels) of the range between the 1st and the 99th percentile of each channel: signals with frequent large one-sided peaks (e.g. movements in accelerometer recordings) fit in their rows (see `_plot_center()`).

For noise-like signals the first term is larger (for Gaussian noise the second is 5.6 SD). The scale is not affected by artifacts lasting < 1% of the signal, which may extend into the neighboring rows. Flat channels and non-finite values are ignored.

# Arguments

- `s::AbstractMatrix`: signal matrix, shape (channels, samples)

# Returns

- `Float64`
"""
function _plot_scale(s::AbstractMatrix)::Float64
    sd = Float64[]
    rng = Float64[]
    for r in eachrow(s)
        x = Base.filter(isfinite, r)
        isempty(x) && continue
        push!(sd, 1.4826 * median(abs.(x .- median(x))))
        push!(rng, quantile(x, 0.99) - quantile(x, 0.01))
    end
    filter!(>(0), sd)
    filter!(>(0), rng)
    sc = max(isempty(sd) ? 0.0 : 6 * median(sd), isempty(rng) ? 0.0 : 1.2 * median(rng))
    return sc > 0 ? sc : 1.0
end

"""
    _plot_center(s)

Center of each channel (rows of `s`) in the multi-channel signal plots: the midpoint between its 1st and 99th percentile. For signals symmetric around the baseline (e.g. EEG) it is close to the median; signals with one-sided peaks (e.g. accelerometer) use the whole height of their rows. Not affected by artifacts lasting < 1% of the signal; non-finite values are ignored.

# Arguments

- `s::AbstractMatrix`: signal matrix, shape (channels, samples)

# Returns

- `Matrix{Float64}`: shape (channels, 1)
"""
function _plot_center(s::AbstractMatrix)::Matrix{Float64}
    c = zeros(size(s, 1), 1)
    for (idx, r) in enumerate(eachrow(s))
        x = Base.filter(isfinite, r)
        isempty(x) || (c[idx] = (quantile(x, 0.01) + quantile(x, 0.99)) / 2)
    end
    return c
end

# channel types whose channels are unrelated signals of different amplitudes (e.g. NIRS
# auxiliary channels): each channel is scaled separately, without a common scale bar
const _PLOT_SCALE_EACH = ["nirs_aux"]

# readable names of channel types (row labels of butterfly plots)
const _CTYPE_LABELS = Dict(
    "nirs_int" => "Intensity", "nirs_od" => "OD", "nirs_hbo" => "HbO",
    "nirs_hbr" => "HbR", "nirs_hbt" => "HbT", "nirs_h2o" => "H₂O",
    "nirs_lipid" => "Lipid", "nirs_bfi" => "BFi", "nirs_aux" => "AUX",
    "nirs_dmean" => "dMean", "nirs_dvar" => "dVar", "nirs_dskew" => "dSkew",
    "nirs_mua" => "μa", "nirs_musp" => "μs′",
)
function _ctype_label(t::AbstractString)::String
    haskey(_CTYPE_LABELS, t) && return _CTYPE_LABELS[t]
    if startswith(t, "nirs_hrf_")
        b = t[10:end]
        return "HRF " * get(_CTYPE_LABELS, "nirs_" * b, uppercase(b))
    end
    return uppercase(t)
end

# scale and center the rows `g_idx` of `s` (channels of one type) in place, negating them for
# the reversed y-axis; return the common scale (`NaN` for types scaled channel by channel)
function _plot_scale_group!(
    s::AbstractMatrix,
    g_idx::AbstractVector{Bool},
    ctype::String,
    k::Real,
)::Float64
    if ctype in _PLOT_SCALE_EACH
        # signals of different devices (e.g. slow drifts, pulse, spikes): each channel is
        # scaled by its own 1st–99th percentile range (the MAD-based scale of `_plot_scale`
        # makes slow or spiky signals look flat)
        for i in findall(g_idx)
            row = s[i:i, :]
            x = Base.filter(isfinite, vec(row))
            rng = isempty(x) ? 0.0 : quantile(x, 0.99) - quantile(x, 0.01)
            sc = rng > 0 ? 1.2 * rng : 1.0
            s[i, :] = .-(row .- _plot_center(row)) ./ (sc * k)
        end
        return NaN
    end
    group = s[g_idx, :]
    sc = _plot_scale(group) * k
    s[g_idx, :] = .-(group .- _plot_center(group)) ./ sc

    return sc
end

"""
    _log_freq_ticks(lo, hi)

Return ticks for a log-transformed frequency axis (e.g. of `GLMakie.Axis3`, which has no log scale): positions are `log10` of the frequencies, labels are the frequencies in Hz. Frequencies of the 1-2-5 sequence within `[lo, hi]` are used; if fewer than 3 of them fit, 1-9 multiples of powers of 10 are used.

# Arguments

- `lo::Real`: lowest frequency (> 0)
- `hi::Real`: highest frequency

# Returns

- `Tuple{Vector{Float64}, Vector{String}}`: tick positions (`log10` of frequency) and labels
"""
function _log_freq_ticks(lo::Real, hi::Real)::Tuple{Vector{Float64}, Vector{String}}
    lo > 0 || throw(ArgumentError("lo must be > 0."))
    hi >= lo || throw(ArgumentError("hi must be ≥ lo."))
    tol = 1e-9 * hi
    dec = floor(Int64, log10(lo)):ceil(Int64, log10(hi))
    pick(m) = [x * 10.0^d for d in dec for x in m if lo - tol <= x * 10.0^d <= hi + tol]
    v = pick((1, 2, 5))
    length(v) < 3 && (v = pick(1:9))
    # very narrow range: the limits only
    length(v) < 2 && (v = unique([Float64(lo), Float64(hi)]))
    lab = [(r = round(x; sigdigits = 3);
            isinteger(r) ? string(round(Int64, r)) : string(r)) for x in v]
    return (log10.(v), lab)
end
