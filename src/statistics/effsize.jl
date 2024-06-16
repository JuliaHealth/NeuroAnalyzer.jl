export effsize
export effsize_p2g

"""
    effsize(x1, x2)

Calculate Cohen's d and Hedges g effect sizes.

# Arguments

- `x1::AbstractVector`
- `x2::AbstractVector`

# Returns

Named tuple containing:
- `d::Float64`: Cohen's d
- `g::Float64`: Hedges g
"""
function effsize(x1::AbstractVector, x2::AbstractVector)

    d = (mean(x2) - mean(x1)) / sqrt((std(x1)^2 + std(x2)^2) / 2)
    g = (mean(x2) - mean(x1)) / sqrt((((length(x1) - 1) * (std(x1)^2)) + ((length(x2) - 1) * (std(x2)^2))) / (length(x1) + length(x2) - 2))

    return (cohen=d, hedges=g)

end

"""
    effsize_p2g(p1, p2)

Calculate effect size for two proportions `p1` and `p2`.

# Arguments

- `p1::Float64`: 1st proportion, e.g. 0.7
- `p2::Float64`: 2nd proportion, e.g. 0.3

# Returns

- `e::Float64`
"""
function effsize_p2g(p1::Float64, p2::Float64)

    @assert p1 + p2 == 1.0 "Proportions must add to 1.0."

    e = 2 * asin(sqrt(p1)) - 2 * asin(sqrt(p2))

    return e

end

