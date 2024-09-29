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

    @assert !(p1 <= 0 || p1 >= 1) "p1 must be in [0, 1]."
    @assert !(p2 <= 0 || p2 >= 1) "p2 must be in [0, 1]."

    p1_zscore = quantile(Distributions.Normal(), p1)
    p2_zscore = quantile(Distributions.Normal(), p2)

    return (dp=(p1_zscore - p2_zscore), rb=(-(p1_zscore + p2_zscore) / 2))

end
