export friedman

"""
    friedman(m)

Estimate Friedman's non-parametric two-way analysis of variance and the associated Kendall's coefficient of concordance W.

Rankings are computed within each group (column) across observations (rows). The Q statistic approximately follows a χ²(k − 1) distribution under H₀.

# Arguments

- `m::AbstractArray`: data matrix of shape `(observations × groups)`; must have ≥ 2 groups (columns) and ≥ 2 observations (rows)

# Returns

Named tuple:

- `q::Float64`: Friedman's Q statistic
- `w::Float64`: Kendall's coefficient of concordance W ∈ [0, 1]
- `p::Float64`: p-value from the χ²(k − 1) distribution

# Throws

- `ArgumentError`: if the matrix has fewer than 2 groups or fewer than 2 observations

# Notes

- H₀ (Friedman): all treatment groups have the same distribution
- H₀ (Kendall): there is no agreement between rankings (W = 0)
- W = 0 indicates no agreement; W = 1 indicates perfect agreement

# References

Friedman M. The use of ranks to avoid the assumption of normality implicit in the analysis of variance. Journal of the American Statistical Association. 1937;32(200):675–701.
"""
function friedman(m::AbstractMatrix)::@NamedTuple{q::Float64, w::Float64, p::Float64}

    # number of observations (blocks)
    n = size(m, 1)
        # number of groups (treatments)
    k = size(m, 2)
    @assert k >= 2 "m must have at least 2 groups (columns)."
    @assert n >= 2 "m must have at least 2 observations (rows)."

    # accumulate rank sums per observation across groups
    rs = zeros(Float64, n)
    for col in axes(m, 2)
        rs .+= ordinalrank(m[:, col])
    end

    # compute the sum-of-squared-deviations from the expected rank sum
    s = sum(rs .^ 2) - (k^2 * n * (n + 1)^2) / 4

    # Friedman's Q statistic
    q = (12 * s) / (k * n * (n + 1))

    # Kendall's W (normalized Q)
    w = (12 * s) / (k^2 * n * (n^2 - 1))

    # p-value: Q ~ χ²(k − 1) under H₀
    p = ccdf(Distributions.Chisq(k - 1), q)

    return (q=q, w=w, p=p)

end
