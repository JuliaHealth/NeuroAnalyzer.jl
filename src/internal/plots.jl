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
