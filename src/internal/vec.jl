function _flipx(s::AbstractVector)::Vector{Float64}
    # center at Y=0
    m = mean(s)
    s_new = s .- m
    # flip the signal along the X axis
    s_new = .-s_new
    s_new .+= m
    return s_new
end

function _zeros(s::AbstractVector)::Vector{Int64}
    z = Int64[]
    for idx in 2:length(s)
        @views s[idx - 1] > 0 && s[idx] < 0 && push!(z, idx)
        @views s[idx - 1] < 0 && s[idx] > 0 && push!(z, idx)
    end
    return z
end
