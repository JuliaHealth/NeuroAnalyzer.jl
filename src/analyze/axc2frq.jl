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

    p_idx = findpeaks1d(c)[1]
    l_pts = l[p_idx]
    frq = sort(unique(round.(1 ./ diff(l_pts), digits=0)))
    deleteat!(frq, isinf.(frq))

    return frq

end