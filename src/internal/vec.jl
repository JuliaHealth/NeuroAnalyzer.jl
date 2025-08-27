function _flipx(s::AbstractVector)::Vector{Float64}
    # center at Y=0
    m = mean(s)
    s_new = s .- m
    # flip the signal along the X axis
    s_new = .-s_new
    s_new .+= m
    return s_new
end

_zeros(s::AbstractVector)::Int64 = count(abs.(diff(sign.(s))) .!= 0)
