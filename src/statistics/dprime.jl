export dprime

"""
    dprime(p1::Real, p2::Real)

Calculate d' and response bias for two proportions.

# Arguments

- `p1::Real`
- `p2::Real`

# Returns

Named tuple containing:
- `dprime::Float64`
- `rb::Float64`: response bias
"""
function dprime(p1::Real, p2::Real)

    p1 in [0, 1] && throw(ArgumentError("p1 must be > 0 and < 1."))
    p2 in [0, 1] && throw(ArgumentError("p2 must be > 0 and < 1."))
    
    p1_zscore = quantile(Distributions.Normal(), p1)
    p2_zscore = quantile(Distributions.Normal(), p2)

    return (dprime=(p1_zscore - p2_zscore), rb=(-(p1_zscore + p2_zscore) / 2))

end
