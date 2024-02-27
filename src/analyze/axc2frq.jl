export axc2frq

"""
   axc2frq(c, l)

Detect peaks in auto-/cross- correlation/covariance and transform them into frequencies.

# Arguments

- `c::AbstractVector`: auto-/cross- correlation/covariance values
- `l::AbstractVector`: lags

# Returns

- `frq::Vector{Float64}`: list of frequencies dominating in the auto-/cross- correlation/covariance
"""
function axc2frq(c::AbstractVector, l::AbstractVector)

    p_idx = findpeaks(c, d=2)
    l_pts = l[p_idx]
    l_pts .+= abs(minimum(l_pts))
    frq = sort(unique(round.(1 ./ round.(diff(l_pts), digits=3), digits=2)))
    deleteat!(frq, isinf.(frq))

    return frq

end