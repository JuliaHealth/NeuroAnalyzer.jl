export friedman

"""
    friedman(m)

Estimate Friedman's nonparametric two-way analysis of variance (and Kendall's coefficient of concordance, normalized Friedman's Q statistics).

# Arguments

- `m::AbstractArray`: values × groups

# Returns

Named tuple containing:
- `q::Float64`: Friedman's Q statistics
- `w::Float64`: Kendall's coefficient of concordance
- `p::Float64`: p value

# Notes

- H0 (Friedman) is that the treatments are equal
- H0 (Kendall) is that there is agreement between rankings or test results
- Kendall's coefficient of concordance ranges from 0 to 1, with 0 meaning no agreement across raters (judges)
"""
function friedman(m::AbstractMatrix)::@NamedTuple{q::Float64, w::Float64, p::Float64}

    rs = zeros(eltype(m), size(m, 1))
    k = size(m, 2)
    n = size(m, 1)
    @assert k >= 2 "Data must have at least two groups."
    for idx in axes(m, 2)
        r = ordinalrank(m[:, idx])
        rs .+= r
    end
    rs2 = rs.^2

    s = sum(rs2) - (((k^2) * n * (n + 1)^2) / 4)
    q = (12 * s) / (k * n * (n + 1))
    w = (12 * s) / ((k^2) * n * (n^2 - 1))
    p = k * (n - 1) * w
    p = pdf(Distributions.Chisq(n - 1), p) * 2

    return (q=q, w=w, p=p)

end
