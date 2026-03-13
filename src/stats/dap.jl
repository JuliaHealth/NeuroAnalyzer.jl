export dap

"""
    dap(x)

Calculate the D'Agostino–Pearson omnibus test for normality.

Combines standardized skewness (`zs`) and standardized kurtosis (`zk`) into a single chi-squared test statistic `d = zs² + zk²`, which follows a `χ²(2)` distribution under the null hypothesis of normality.

# Arguments

- `x::AbstractVector`: input vector; must contain at least 8 elements (`ses` and `sek` require sufficient sample size); a warning is issued for `length(x) < 20` as the approximation is unreliable for small samples

# Returns

Named tuple:

- `zs::Float64`: standardized skewness statistic (`skewness(x) / ses(x)`)
- `zk::Float64`: standardized kurtosis statistic (`kurtosis(x) / sek(x)`)
- `d::Float64`: omnibus test statistic `zs² + zk²`; follows `χ²(2)` under H₀
- `p::Float64`: p-value from the `χ²(2)` survival function

# Throws

- `ArgumentError`: if `length(x) < 8`

# Notes

- The test requires `n ≥ 20` for reliable results; a diagnostic message is printed for smaller samples.
- H₀: the data are normally distributed.
- H₁: the data are not normally distributed.

# References

D'Agostino RB, Belanger A, D'Agostino RB Jr. A suggestion for using powerful and informative tests of normality. The American Statistician. 1990;44(4):316–21.
"""
function dap(x::AbstractVector)::@NamedTuple{zs::Float64, zk::Float64, d::Float64, p::Float64}

    @assert length(x) >= 8 "x must contain at least 8 elements."
    length(x) < 20 && _info("DAP test is unreliable for length(x) < 20; got $(length(x)).")

    # standardize skewness and kurtosis by their respective standard errors
    zs = skewness(x) / ses(x)
    zk = kurtosis(x) / sek(x)

    # omnibus statistic: sum of squared Z-scores ~ χ²(2) under normality
    d = zs^2 + zk^2

    # p-value from the upper tail of χ²(2)
    p = ccdf(Distributions.Chisq(2), d)

    return (zs=zs, zk=zk, d=d, p=p)

end
