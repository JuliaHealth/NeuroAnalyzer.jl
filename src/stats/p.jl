export p2z
export z2p
export t2p
export chi2p
export f2p
export norminv
export p2o
export o2p

"""
    p2z(p; <keyword arguments>)

Convert a p-value to the corresponding Z-score.

# Arguments

- `p::Float64=0.05`: p-value; must be in `(0, 1)`
- `twotailed::Bool=false`: if `true`, compute the two-tailed Z-score; if `false` (default), compute the one-tailed Z-score

# Returns

- `Float64`: Z-score such that `P(Z > z) = p` (one-tailed) or `P(|Z| > z) = p` (two-tailed)
"""
function p2z(p::Float64 = 0.05; twotailed::Bool = false)::Float64

    # validate
    _in(p, (0.0, 1.0), "p")
    d = Distributions.Normal(0.0, 1.0)

    return twotailed ? quantile(d, 1 - p / 2) : quantile(d, 1 - p)

end

"""
    z2p(z; <keyword arguments>)

Convert a Z-score to a p-value.

# Arguments

- `z::Real`: Z-score
- `twotailed::Bool=false`: if `true`, compute the two-tailed probability `P(|Z| > |z|)`; if `false`, compute the one-tailed probability `P(Z > z)` for `z ≥ 0` or `P(Z < z)` for `z < 0`.

# Returns

- `Float64`: p-value ∈ `(0, 1)`

# Notes

Derived probabilities:

- `P(Z < z)` (left-tailed): `1 - z2p(z; twotailed=false)`
- `P(|Z| < z)`: `1 - z2p(z; twotailed=true)`
- `P(0 < Z < z)`: `(1 - z2p(z; twotailed=true)) / 2`
"""
function z2p(z::Real; twotailed::Bool = false)::Float64

    d = Distributions.Normal(0.0, 1.0)
    if twotailed
        return 2 * ccdf(d, abs(z))
    else
        # one-tailed: P(Z > z) for z ≥ 0; P(Z < z) for z < 0
        return z >= 0 ? ccdf(d, z) : cdf(d, z)   # was: 1 - ccdf(d, abs(z)) == cdf(d, -|z|) == cdf(d, z)
    end

end

"""
    t2p(t; <keyword arguments>)

Convert a t-score to a p-value using the Student's t-distribution.

# Arguments

- `t::Real`: t-statistic
- `df::Real`: degrees of freedom; must be > 0
- `twotailed::Bool=false`: if `true`, compute `P(|T| > |t|)`; if `false`, compute `P(T > t)` for `t ≥ 0` or `P(T < t)` for `t < 0`

# Returns

- `Float64`: p-value ∈ `(0, 1)`

# Notes

Derived probabilities:

- `P(T < t)` (left-tailed): `1 - t2p(t; df=df, twotailed=false)`
- `P(|T| < t)`: `1 - t2p(t; df=df, twotailed=true)`
- `P(0 < T < t)`: `(1 - t2p(t; df=df, twotailed=true)) / 2`
"""
function t2p(t::Real; df::Real, twotailed::Bool = false)::Float64

    # validate
    df > 0 || throw(ArgumentError("df must be > 0."))

    d = Distributions.TDist(df)
    if twotailed
        return 2 * ccdf(d, abs(t))
    else
        return t >= 0 ? ccdf(d, t) : cdf(d, t)
    end

end

"""
    chi2p(chi; <keyword arguments>)

Convert a χ² statistic to a right-tailed p-value `P(χ² > chi)`.

# Arguments

- `chi::Real`: χ² statistic; must be ≥ 0
- `df::Real`: degrees of freedom; must be > 0

# Returns

- `Float64`: right-tailed p-value

# Notes

To obtain the left-tailed probability `P(χ² < chi)` use `1 - chi2p(chi; df=df)`.
"""
function chi2p(chi::Real; df::Real)::Float64

    # validate
    chi >= 0 || throw(ArgumentError("chi must be ≥ 0."))
    df > 0 || throw(ArgumentError("df must be > 0."))

    return ccdf(Distributions.Chisq(df), chi)

end


"""
    f2p(t; <keyword arguments>)

Convert an F-statistic to a right-tailed p-value `P(F > f)`.

# Arguments

- `f::Real`: F-statistic; must be ≥ 0
- `df1::Real`: numerator degrees of freedom; must be > 0
- `df2::Real`: denominator degrees of freedom; must be > 0

# Returns

- `Float64`: right-tailed p-value

# Notes

To obtain the left-tailed probability `P(F < f)` use `1 - f2p(f; df1=df1, df2=df2)`.
"""
function f2p(f::Real; df1::Real, df2::Real)::Float64

    # negative F silently returns > 1
    f >= 0 || throw(ArgumentError("f must be ≥ 0."))
    df1 > 0 || throw(ArgumentError("df1 must be > 0."))
    df2 > 0 || throw(ArgumentError("df2 must be > 0."))

    return ccdf(Distributions.FDist(df1, df2), f)

end

export norminv

"""
    norminv(x::Real)

Return the quantile of the standard normal distribution at probability `x`.

Equivalent to `Φ⁻¹(x)` where `Φ` is the standard normal CDF.

# Arguments

- `x::Real`: probability; must be in `(0, 1)`

# Returns

- `Float64`: normal quantile at `x`
"""
function norminv(x::Real)::Float64

    # norminv(0) = -Inf, norminv(1) = Inf
    0 < x < 1 || throw(ArgumentError("x must be in (0, 1)."))

    return quantile(Distributions.Normal(), x)

end

"""
    p2o(p::Float64)

Convert a probability to odds.

Computed as `o = p / (1 − p)`.

# Arguments

- `p::Real`: probability; must be in `[0, 1)` (odds undefined at `p = 1`)

# Returns

- `Float64`: odds
"""
function p2o(p::Real)::Float64

    # validate
    p >= 0 || throw(ArgumentError("p must be ≥ 0."))
    p <  1 || throw(ArgumentError("p must be < 1 (odds undefined at p = 1)."))

    return p / (1 - p)

end

"""
    o2p(o::Float64)

Convert odds to a probability.

Computed as `p = o / (1 + o)`.

# Arguments

- `o::Real`: odds; must be ≥ 0

# Returns

- `Float64`: probability ∈ `[0, 1)`
"""
function o2p(o::Real)::Float64

    # negative odds produce p > 1
    o >= 0 || throw(ArgumentError("o must be ≥ 0."))

    return o / (1 + o)

end
