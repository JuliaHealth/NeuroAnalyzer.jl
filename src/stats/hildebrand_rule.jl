export hildebrand_rule

"""
    hildebrand_rule(x; verbose)

Calculate the Hildebrand rule statistic for a vector.

Computed as `H = (mean(x) − median(x)) / std(x)`. Values of `|H| < 0.2` indicate approximate symmetry of the distribution.

# Arguments

- `x::AbstractVector`: input vector; must contain at least 2 elements and have non-zero standard deviation
- `verbose::Bool=true`: if `true`, print a message when the symmetry criterion is met

# Returns

- `Float64`: Hildebrand statistic H

# Throws

- `ArgumentError`: if `length(x) < 2` or `std(x) == 0`

# Notes

- The criterion `|H| < 0.2` is checked (not just `H < 0.2`) because a negative mean–median difference also indicates asymmetry if its magnitude is ≥ 0.2.
"""
function hildebrand_rule(x::AbstractVector; verbose::Bool = true)::Float64

    @assert length(x) >= 2 "x must contain at least 2 elements."
    s = std(x)
    @assert s != 0 "std(x) must not be zero."

    h = (mean(x) - median(x)) / s
    abs(h) < 0.2 && _info("H < 0.2: x is approximately symmetric")

    return h

end
