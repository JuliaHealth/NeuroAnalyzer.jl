export axc2frq

"""
    axc2frq(c, l)

Finds peaks in an auto-/cross-correlation or covariance vector, converts their lag positions to inter-peak intervals, and returns the corresponding frequencies (1 / interval).

# Arguments

- `c::AbstractVector`: auto-/cross- correlation/covariance values
- `l::AbstractVector`: lags corresponding to each value in `c`

# Returns

- `frq::Vector{Float64}`: frequencies (Hz) dominating in the auto-/cross-correlation/covariance, sorted in ascending order
"""
function axc2frq(c::AbstractVector, l::AbstractVector)::Vector{Float64}

    # find indices of local peaks in the correlation/covariance vector
    # d=2 sets the minimum distance between accepted peaks
    p_idx = findpeaks(c, d = 2)

    # extract the lag values at the peak positions
    l_pts = l[p_idx]

    # shift all lags so the minimum is zero, turning negative symmetric lags
    # into a non-negative interval sequence suitable for diff
    l_pts .+= abs(minimum(l_pts))

    # compute inter-peak intervals, convert to frequencies, round and de-duplicate
    frq = sort(unique(round.(1 ./ round.(diff(l_pts), digits = 3), digits = 2)))
    Base.filter!(!isinf, frq)

    return frq

end
