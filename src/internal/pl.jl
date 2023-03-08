_pl(x::Union{AbstractRange, AbstractVector}) = length(collect(x)) > 1 ? "s" : ""

_pl(x::Real) = x > 1 ? "s" : ""
