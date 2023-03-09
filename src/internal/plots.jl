_xlims(t::Union{Vector{<:Real}, AbstractRange}) = floor(t[1], digits=2), ceil(t[end], digits=2)

_ticks(t::Union{Vector{<:Real}, AbstractRange}) = floor(t[1], digits=2):((ceil(t[end]) - floor(t[1])) / 10):ceil(t[end], digits=2)

_ticks(t::Tuple{Real, Real}) = floor(t[1], digits=2):((ceil(t[2]) - floor(t[1])) / 10):ceil(t[2], digits=2)

_erpticks(t::Union{Vector{<:Real}, AbstractRange}) = vcat(collect(range(floor(t[1], digits=2), 0, 3)), collect(range(0, ceil(t[end], digits=2), 9))[2:end])

_erpticks(t::Tuple{Real, Real}) = vcat(collect(range(floor(t[1], digits=2), 0, 3)), collect(range(0, ceil(t[2], digits=2), 9))[2:end])

function _set_defaults(xl::String, yl::String, tt::String, x::String, y::String, t::String)
    xl == "default" && (xl = x)
    yl == "default" && (yl = y)
    tt == "default" && (tt = t)
    return xl, yl, tt
end
