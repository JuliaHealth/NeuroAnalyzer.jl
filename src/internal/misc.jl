_info(s::String) = verbose == true && @info s

_warn(s::String) = verbose == true && @warn s

_pl(x::Union{AbstractRange, AbstractVector}) = length(collect(x)) > 1 ? "s" : ""

_pl(x::Real) = x > 1 ? "s" : ""

# _get_range(signal::Union{AbstractVector, AbstractArray}) = round(abs(minimum(signal)) + abs(maximum(signal)), digits=0)
_get_range(signal::Union{AbstractVector, AbstractArray}) = round(rng(signal), digits=0)

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

function _copy_lt2ut(m::AbstractArray)
    if ndims(m) == 2
        return m + m' - diagm(diag(m))
    else
        Threads.@threads for ep_idx in 1:size(m, 3)
            @inbounds m[:, :, ep_idx] = m[:, :, ep_idx] + m[:, :, ep_idx]' - diagm(diag(m[:, :, ep_idx]))
        end
        return m
    end
end

_tlength(t::Tuple{Real, Real}) = length(t[1]:1:t[2])

function _s2i(x::String)
    if occursin(":", x)
        return collect(parse(Int64, split(x, ":")[1]):parse(Int64, split(x, ":")[2]))
    elseif occursin(", ", x)
        x = replace(x, "["=>"")
        x = replace(x, "]"=>"")
        return parse.(Int64, split(x, ", "))
    elseif occursin(",", x)
        x = replace(x, "["=>"")
        x = replace(x, "]"=>"")
        return parse.(Int64, split(x, ","))
    end
end

function _i2s(x::Union{AbstractRange, Vector{Int64}})
    x = collect(x)
    x = string(x)
    x = replace(x, "["=>"")
    x = replace(x, "]"=>"")
    return x
end