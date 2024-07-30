export friedman

"""
    friedman(m)

Estimate Friedman's nonparametric two-way analysis of variance (and Kendall's coefficient of concordance).

# Arguments

- `m::AbstractArray`: values × groups

# Returns

Named tuple containing:
- `f::Float64`: Friedman
- `k::Float64`: Kendall
- `p::Float64`: P-value

# Notes

- H0 (Friedman) is that the treatments are equal
- H0 (Kendall) is that there is agreement between rankings or test results
- Kendall's coefficient of concordance ranges from 0 to 1, with 0 meaning no agreement across raters (judges)
"""
function friedman(m::AbstractMatrix)
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
    f = (12 * s) / (k * n * (n + 1))
    w = (12 * s) / ((k^2) * n * (n^2 - 1))
    p = k * (n - 1) * w
    p = pdf(Distributions.Chisq(n - 1), p) * 2

    return(f=f, k=w, p=p)
end
