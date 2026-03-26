function _draw_head_outline!(ax::GLMakie.Axis; lw::Int64)
    # nose
    GLMakie.lines!(ax, [-0.2, 0], [0.98, 1.08], linewidth=lw, color=:black)
    GLMakie.lines!(ax, [0.2, 0], [0.98, 1.08], linewidth=lw, color=:black)
    # ears
    ear_coords = [(-1.03,0.15), (-1.06,0.16), (-1.1,0.14), (-1.12,0.05),
                 (-1.1, -0.1), (-1.13,-0.3), (-1.09,-0.37), (-1.02,-0.39),
                 (-0.98,-0.33), (-0.975,-0.22), (1.03,0.15), (1.06,0.16),
                 (1.1,0.14), (1.12,0.05), (1.1, -0.1), (1.13,-0.3),
                 (1.09,-0.37), (1.02,-0.39), (0.98,-0.33), (0.975,-0.22)]
    GLMakie.lines!(ax, first.(ear_coords), last.(ear_coords), linewidth=lw, color=:black)
    # head outline
    GLMakie.arc!(ax, (0,0), 1, 0, 2π, linewidth=lw, color=:black)
end

_xlims(t::Union{AbstractVector, AbstractRange})::Tuple{Real, Real} = floor(t[1], digits = 2), ceil(t[end], digits = 2)

function _ylims(s::Union{AbstractVector, AbstractMatrix})::Tuple{Real, Real}
    if maximum(abs.(s)) > 100
        n = 2
    elseif maximum(abs.(s)) >= 10
        n = 1
    elseif maximum(abs.(s)) < 10
        n = 0
    end
    max = ceil(Int64, round(maximum(s), digits = n))
    min = floor(Int64, round(minimum(s), digits = n))
    if abs(min) == 0 && abs(max) == 0
        max = 1.0
        min = -1.0
    end
    if min == 0
        max = abs(max)
        min = -abs(max)
    elseif max == 0
        max = abs(min)
        min = -abs(min)
    end
    if abs(max) > abs(min)
        return (-abs(max), abs(max))
    else
        return (-abs(min), abs(min))
    end
end

function _ticks(t::Union{AbstractVector, AbstractRange, Tuple{Real, Real}})::AbstractVector
    tc = collect(t[1]:1:t[end])
    if t[end] - t[1] > 10
        tc = collect(t[1]:2:t[end])
    elseif t[end] - t[1] > 20
        tc = collect(t[1]:2:t[end])
    elseif t[end] - t[1] > 50
        tc = collect(t[1]:5:t[end])
    end
    return tc
end

_erpticks(t::Union{AbstractVector, AbstractRange})::AbstractVector = vcat(
    collect(range(floor(t[1], digits = 2), 0, 3)), collect(range(0, ceil(t[end], digits = 2), 9))[2:end]
)

_erpticks(t::Tuple{Real, Real})::AbstractVector = vcat(
    collect(range(floor(t[1], digits = 2), 0, 3)), collect(range(0, ceil(t[2], digits = 2), 9))[2:end]
)

function _set_defaults(
        xl::String, yl::String, tt::String, x::String, y::String, t::String
    )::Tuple{String, String, String}
    yl == "default" && (yl = y)
    xl == "default" && (xl = x)
    tt == "default" && (tt = t)
    return xl, yl, tt
end

_bernstein(i, n; steps = 50) = [binomial(n, i) * t^i * (1 - t)^(n - i) for t in LinRange(0, 1, steps)]

function _bernstein_poly(px, py; steps = 50)
    # the code is based on https://opensourc.es/blog/bezier-curve/
    n = length(px) - 1
    b = [_bernstein(i, n) for i in 0:n]
    x_vals = [sum(px[k] * b[k][t] for k in 1:(n + 1)) for t in 1:steps]
    y_vals = [sum(py[k] * b[k][t] for k in 1:(n + 1)) for t in 1:steps]
    return x_vals, y_vals
end
