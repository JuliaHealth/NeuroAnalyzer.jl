_info(s::String) = verbose == true && @info s

_warn(s::String) = verbose == true && @warn s

_pl(x::Union{AbstractRange, AbstractVector}) = length(collect(x)) > 1 ? "s" : ""

_pl(x::Real) = x > 1 ? "s" : ""

_get_range(signal::Union{AbstractVector, AbstractArray}) = round(abs(minimum(signal)) + abs(maximum(signal)), digits=0)

_c(n) = collect(1:n)

function _tuple_max(t::Tuple{Real, Real})
    abs(t[1]) > abs(t[2]) && (t = (-abs(t[1]), abs(t[1])))
    abs(t[1]) < abs(t[2]) && (t = (-abs(t[2]), abs(t[2])))
    return t
end

function _s2v(s::Union{<:Number, Vector{<:Number}})
    if typeof(s) <: Number
        return [s]
    else
        return s
    end
end

