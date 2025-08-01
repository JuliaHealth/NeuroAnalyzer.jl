export dap

"""
    dap(x)

Calculate D'Agostino-Pearson Omnibus Test for normality.

# Arguments

- `x::AbstractVector`

# Returns

Named tuple containing:
- `zs::Float64`: skewness test
- `zk::Float64`: kurtosis test
- `d::Float64`: test statistic
- `p::Float64`: p value
"""
function dap(x::AbstractVector)::@NamedTuple{zs::Float64, zk::Float64, d::Float64, p::Float64}

    length(x) < 20 && _info("For DAP test length of x should be ≥ 20.")

    zs = skewness(x) / ses(x)
    zk = kurtosis(x) / sek(x)

    # calculate test statistic
    d = zs^2 + zk^2

    # calculate p value for d from Chi2(df=2) distribution
    p = ccdf(Distributions.Chisq(2), d)

    return (zs=zs, zk=zk, d=d, p=p)

end
