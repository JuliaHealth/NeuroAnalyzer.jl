_xlims(t::Union{AbstractVector, AbstractRange})::Tuple{Real, Real} = floor(t[1], digits=2), ceil(t[end], digits=2)

function _ylims(s::AbstractVector)::Tuple{Real, Real}
    if maximum(abs.(s)) > 100
        n = 2
    elseif maximum(abs.(s)) >= 0
        n = 1
    elseif maximum(abs.(s)) < 0
        n = 0
    end
    max = ceil(Int64, round(maximum(s) * 1.5, digits=n))
    min = floor(Int64, round(minimum(s) * 1.5, digits=n))
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

function _ticks(t::Union{AbstractVector, AbstractRange})::AbstractVector
#    if length(t) >= 3
#        if t[2] - t[1] == t[3] - t[2]
#            tc = linspace(round(t[1]), round(t[end]), length(round.(Int64, t[1]:t[end])))
#        else
#            tc = collect(floor(t[1], digits=3):((ceil(t[end]) - floor(t[1])) / 10):ceil(t[end], digits=3))
#        end
#    end
    if t[end] == round(Int64, t[end])
        tc = linspace(round(t[1]), round(t[end]), length(round.(Int64, t[1]:((t[end] - t[1])/10):t[end])))
    else
        tc = linspace(t[1], t[end], length(round.(t[1]:((t[end] - t[1])/10):t[end], digits=2)))
    end
    tc = round.(tc, digits=2)
    return tc
end

function _ticks(t::Tuple{Real, Real})::AbstractVector
    if typeof(t[1]) <: Int && typeof(t[2]) <: Int
        if length(t[1]:t[2]) <= 30
            return collect(t[1]:t[2])
        elseif length(t[1]:t[2]) <= 100
            return collect(t[1]:5:t[2])
        elseif length(t[1]:t[2]) <= 1000
            return collect(t[1]:10:t[2])
        end
    else
        if length(collect(t[1]:t[2])) > (1 / (collect(t[1]:t[2])[2] - collect(t[1]:t[2])[1]))
            return floor(t[1], digits=3):((ceil(t[end]) - floor(t[1])) / 10):ceil(t[end], digits=2)
        else
            return floor(t[1], digits=3):((ceil(t[end]) - floor(t[1])) / 20):ceil(t[end], digits=2)
        end
    end
end

_erpticks(t::Union{AbstractVector, AbstractRange})::AbstractVector = vcat(collect(range(floor(t[1], digits=3), 0, 3)), collect(range(0, ceil(t[end], digits=2), 9))[2:end])

_erpticks(t::Tuple{Real, Real})::AbstractVector = vcat(collect(range(floor(t[1], digits=3), 0, 3)), collect(range(0, ceil(t[2], digits=2), 9))[2:end])

function _set_defaults(xl::String, yl::String, tt::String, x::String, y::String, t::String)::Tuple{String, String, String}
    yl == "default" && (yl = y)
    xl == "default" && (xl = x)
    tt == "default" && (tt = t)
    return xl, yl, tt
end

_bernstein(i, n; steps=50) = [binomial(n, i) * t^i * (1 - t)^(n - i) for t in LinRange(0, 1, steps)]

function _bernstein_poly(px, py; steps=50)
    # the code is based on https://opensourc.es/blog/bezier-curve/
    n = length(px) - 1
    b = [_bernstein(i, n) for i in 0:n]
    x_vals = [sum(px[k] * b[k][t] for k in 1:(n + 1)) for t in 1:steps]
    y_vals = [sum(py[k] * b[k][t] for k in 1:(n + 1)) for t in 1:steps]
    return x_vals, y_vals
end

