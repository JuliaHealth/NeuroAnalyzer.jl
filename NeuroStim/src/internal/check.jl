function _check_tuple(t::Tuple{Real, Real}, name::String, range::Union{Nothing, Tuple{Real, Real}}=nothing)
    @assert t == tuple_order(t) "$name must contain two values in ascending order."
    @assert t[1] < t[2] "$name must contain two different values in ascending order."
    if range !== nothing 
        @assert !(t[1] < range[1] || t[2] < range[1] || t[1] > range[2] || t[2] > range[2]) "$name must be in [$(range[1]), $(range[2])]."
    end
    return nothing
end

function _check_var(s1::Symbol, s2::Vector{Symbol}, var::String)
    if length(s2) > 1
        m = var * " must be "
        for idx in 1:(length(s2) - 2)
            m *= ":" * string(s2[idx]) * ", "
        end
        m *= ":" * string(s2[end - 1]) * " or :" * string(s2[end]) * "."
        @assert s1 in s2 "$m"
    else
        m = var * " must be :" * string(s2[1])
        @assert s1 in s2 "$m"
    end
    return nothing
end

function _check_var(s1::String, s2::Vector{String}, var::String)
    if length(s2) > 1
        m = var * " must be "
        for idx in 1:(length(s2) - 2)
            m *= s2[idx] * ", "
        end
        m *= s2[end - 1] * " or " * s2[end] * "."
        @assert s1 in s2 "$m"
    else
        m = var * " must be " * string(s2[1])
        @assert s1 in s2 "$m"
    end
    return nothing
end
