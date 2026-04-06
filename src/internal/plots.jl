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
    font_size::Int64 = 8
)::Nothing
    _check_var(plane, [:xy, :xz, :yz], "plane")
    fid_names = ["NAS", "IN", "LPA", "RPA"]
    for idx in eachindex(NeuroAnalyzer.fiducial_points)
        pt = NeuroAnalyzer.fiducial_points[idx]
        fid_loc_x, fid_loc_y = if plane === :xy
            pt[1], pt[2]
        elseif plane === :xz
            pt[1], pt[3]
        elseif plabe === :yz
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
    _draw_head_outline!(ax; <keyword arguments>)

Draw a schematic head outline (circle, nose, ears) onto a 2-D axis.

The head is assumed to have unit radius centered at the origin.

# Arguments

- `ax::GLMakie.Axis`: target axis
- `lw::Int64=1`: line width in points
"""
function _draw_head_outline!(ax::GLMakie.Axis; lw::Int64 = 1)
    # nose
    GLMakie.lines!(ax, [-0.2, 0.0], [0.98, 1.08]; linewidth = lw, color = :black)
    GLMakie.lines!(ax, [0.2,  0.0], [0.98, 1.08]; linewidth = lw, color = :black)

    # left ear
    left_ear_x  = [-0.995, -1.03, -1.06, -1.1, -1.12, -1.1, -1.13, -1.09, -1.02, -0.98, -0.975]
    left_ear_y  = [0.1, 0.15, 0.16, 0.14, 0.05, -0.1, -0.3, -0.37, -0.39, -0.33, -0.22]
    GLMakie.lines!(ax, left_ear_x, left_ear_y; linewidth = lw, color = :black)

    # right ear
    right_ear_x = [0.995, 1.03, 1.06, 1.1, 1.12, 1.1, 1.13, 1.09, 1.02, 0.98, 0.975]
    right_ear_y = [0.1, 0.15, 0.16, 0.14, 0.05, -0.1, -0.3, -0.37, -0.39, -0.33, -0.22]
    GLMakie.lines!(ax, right_ear_x, right_ear_y; linewidth = lw, color = :black)

    # head circle
    GLMakie.arc!(ax, Point2f(0, 0), 1, 0, 2pi; linewidth = lw, color = :black)

    return nothing
end

"""
    _xlims(t)

Return `(floor(t[1], digits=2), ceil(t[end], digits=2))` as axis x-limits.
"""
_xlims(t::Union{AbstractVector, AbstractRange})::Tuple{Real, Real} =
    floor(t[1], digits = 2), ceil(t[end], digits = 2)

"""
    _ylims(s)

Return symmetric y-limits `(-m, m)` suitable for displaying signal `s`.

The magnitude `m` is derived from the signal's peak absolute value, rounded to a precision that scales with the signal's dynamic range.
"""
function _ylims(s::Union{AbstractVector, AbstractMatrix})::Tuple{Real, Real}
    peak = maximum(abs, s)

    n = peak > 100 ? 2 : peak >= 10 ? 1 : 0

    hi = ceil(Int64,  round(maximum(s); digits = n))
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
    t1  = t isa Tuple ? t[1] : t[1]
    t2  = t isa Tuple ? t[2] : t[end]
    rng = t2 - t1
    step = rng > 50 ? 10 : rng > 20 ? 5 : rng > 10 ? 2 : 1
    return collect(t1:step:t2)
end

"""
    _erpticks(t)

Return a tick vector suitable for ERP plots: 3 ticks from `floor(t[1])` to 0, then 8 more ticks from 0 to `ceil(t[end])`.

Accepts either a vector/range or a `Tuple{Real, Real}`.
"""
function _erpticks(t::Union{AbstractVector, AbstractRange, Tuple{Real, Real}})::AbstractVector
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
    x::String,  y::String,  t::String,
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

    b = [_bernstein(i, n; steps = steps) for i in 0:n]
    x_vals = [sum(px[k] * b[k][t] for k in 1:(n + 1)) for t in 1:steps]
    y_vals = [sum(py[k] * b[k][t] for k in 1:(n + 1)) for t in 1:steps]
    return x_vals, y_vals
end