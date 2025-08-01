export dprime

"""
    dprime(p1::Real, p2::Real)

Calculate d' and response bias for two proportions.

# Arguments

- `p1::Real`
- `p2::Real`

# Returns

Named tuple containing:
- `dp::Float64`
- `rb::Float64`: response bias
"""
function dprime(p1::Real, p2::Real)::@NamedTuple{dp::Float64, rb::Float64}

    _in(p1, (0.0, 1.0), "p1")
    _in(p2, (0.0, 1.0), "p2")

    p1_zscore = quantile(Distributions.Normal(), p1)
    p2_zscore = quantile(Distributions.Normal(), p2)

    return (dp=(p1_zscore - p2_zscore), rb=(-(p1_zscore + p2_zscore) / 2))

end
