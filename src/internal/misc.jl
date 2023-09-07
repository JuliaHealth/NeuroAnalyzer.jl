_info(s::String) = verbose == true && @info s

_warn(s::String) = verbose == true && @warn s

_deprecated(s1::String, s2::String) = verbose == true && @error "Function $s1() is deprecated, please use $s2() instead."

_wip() = allow_wip == true ? (@warn "This function has the WIP (Work In Progress) status and is not ready for production use.") : (@error "This function has the WIP (Work In Progress) status and is not ready for production use.")

_pl(x::Union{AbstractRange, AbstractVector}) = length(collect(x)) > 1 ? "s" : ""

_pl(x::Real) = x > 1 ? "s" : ""

# _get_range(signal::Union{AbstractVector, AbstractArray}) = round(abs(minimum(signal)) + abs(maximum(signal)), digits=0)
_get_range(signal::Union{AbstractVector, AbstractArray}) = round(rng(signal), digits=0)

_c(n) = collect(1:n)

_tuple_max(t::Tuple{Real, Real}) = abs(t[1]) > abs(t[2]) ? (-abs(t[1]), abs(t[1])) : (-abs(t[2]), abs(t[2]))

_s2v(s::Union{<:Number, Vector{<:Number}}) = typeof(s) <: Number ? [s] : s

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

function _s2i(s::String)
    s = replace(s, " "=>"")
    if occursin(":", s)
        return collect(parse(Int64, split(s, ":")[1]):parse(Int64, split(s, ":")[2]))
    elseif occursin(",", s)
        s = replace(s, "["=>"")
        s = replace(s, "]"=>"")
        return parse.(Int64, split(s, ","))
    elseif _check_sint(s)
        return parse.(Int64, s)
    end
end

function _i2s(s::Union{Int64, Vector{Int64}, AbstractRange})
    !isa(s, Int64) && (s = collect(s))
    s = string(s)
    s = replace(s, "["=>"")
    s = replace(s, "]"=>"")
    return s
end

function _s2tf(s::String)
    s = replace(s, " "=>"")
    s = replace(s, "("=>"")
    s = replace(s, ")"=>"")
    return (parse(Float64, split(s, ",")[1]), parse(Float64, split(s, ",")[2]))
end

function _s2ti(s::String)
    s = replace(s, " "=>"")
    s = replace(s, "("=>"")
    s = replace(s, ")"=>"")
    return (parse(Int64, split(s, ",")[1]), parse(Int64, split(s, ",")[2]))
end
