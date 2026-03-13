export dprime

"""
    dprime(p1::Real, p2::Real)

Calculate the sensitivity index d′ and response bias for two proportions.

`p1` is typically the hit rate and `p2` the false-alarm rate. Both are converted to Z-scores using the inverse normal CDF, then combined as:

- `d′ = Z(p1) − Z(p2)`
- `bias = −(Z(p1) + Z(p2)) / 2`

Boundary values `p ∈ {0, 1}` produce ±Inf; callers should apply a correction (e.g. the Hautus log-linear rule) before passing boundary proportions.

# Arguments

Named tuple:

- `dp::Float64`: sensitivity index d′
- `rb::Float64`: response bias (negative = liberal, positive = conservative)

# Throws

- `ArgumentError`: if `p1` or `p2` are outside `(0, 1)`

# References

Green DM, Swets JA. Signal Detection Theory and Psychophysics. Wiley; 1966.
"""
function dprime(p1::Real, p2::Real)::@NamedTuple{dp::Float64, rb::Float64}

    _in(p1, (0.0, 1.0), "p1")
    _in(p2, (0.0, 1.0), "p2")

    z1 = quantile(Distributions.Normal(), p1)
    z2 = quantile(Distributions.Normal(), p2)

    return (dp=z1 - z2, rb=-(z1 + z2) / 2)

end
