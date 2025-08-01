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

Convert p value to z score.

# Arguments

- `p::Float64=0.05`: p value
- `twotailed::Bool=false`: one- or two-tailed probability

# Returns

- `z::Float64`
"""
function p2z(p::Float64=0.05; twotailed::Bool=false)::Float64

    _in(p, (0.0, 1.0), "p")

    d = Distributions.Normal(0.0, 1.0)
    if twotailed
        z = quantile(d, 1 - p / 2)
    else
        z = quantile(d, 1 - p)
    end

    return z

end

"""
    z2p(z; <keyword arguments>)

Convert z score to p value.

# Arguments

- `z::Real`: z value
- `twotailed::Bool=false`: one- (`P(Z > z)`) or two-tailed (`P(Z < -z or Z > z)`) probability

# Returns

- `p::Float64`

# Notes

To calculate:
- `P(Z < z)` (left-tailed) use `1 - z2p(z, twotailed=false)`
- `P(-z < Z < z)` use `1 - z2p(z, twotailed=true)`
- `P(0 < Z < z)` use `(1 - z2p(z, twotailed=true))/2`
"""
function z2p(z::Real; twotailed::Bool=false)::Float64

    d = Distributions.Normal(0.0, 1.0)
    if twotailed
        p = 2 * ccdf(d, abs(z))
    else
        if z >= 0
            p = ccdf(d, z)
        else
            p = 1 - ccdf(d, abs(z))
        end
    end

    return p

end

"""
    t2p(t; <keyword arguments>)

Convert t score to p value.

# Arguments

- `t::Real`: t score
- `df::Real`: degrees of freedom
- `twotailed::Bool=false`: one- (`P(Z > t)`) or two-tailed (`P(T < -t or T > t)`) probability

# Returns

- `p::Float64`

# Notes

To calculate:
- `P(T < t)` (left-tailed) use `1 - t2p(t, df=df, twotailed=false)`
- `P(-t < T < t)` use `1 - t2p(t, df=df, twotailed=true)`
- `P(0 < T < t)` use `(1 - t2p(t, df=df, twotailed=true))/2`
"""
function t2p(t::Real; df::Real, twotailed::Bool=false)::Float64

    d = Distributions.TDist(df)
    if twotailed
        p = 2 * ccdf(d, abs(t))
    else
        if t >= 0
            p = ccdf(d, t)
        else
            p = 1 - ccdf(d, abs(t))
        end
    end

    return p

end

"""
    chi2p(chi; <keyword arguments>)

Convert Χ² score to right-tailed (`P(Chi > chi)`) p value.

# Arguments

- `chi::Real`: Χ² score
- `df::Real`: degrees of freedom

# Returns

- `p::Float64`

# Notes

To calculate `P(Chi < chi)` (left-tailed) use `1 - chi2p(chi, df=df)`.
"""
function chi2p(chi::Real; df::Real)::Float64

    p = ccdf(Distributions.Chisq(df), chi)

    return p

end


"""
    f2p(t; <keyword arguments>)

Convert F score to right-tailed (`P(F > f)`) p value.

# Arguments

- `f::Real`: F score (`F = var(x) / var(y)`)
- `df1::Real`: numerator degrees of freedom (DF1)
- `df2::Real`: denominator degrees of freedom (DF2)

# Returns

- `p::Float64`

# Notes

To calculate `P(F < f)` (left-tailed) use `1 - f2p(f, df1=df1, df2=df2)`.
"""
function f2p(f::Real; df1::Real, df2::Real)::Float64

    p = ccdf(Distributions.FDist(df1, df2), f)

    return p

end

export norminv

"""
    norminv(x::Real)

Convert probability to a normal distribution with a peak at 0.5.

# Arguments

- `x::Real`

# Returns

- `n::Float64`
"""
function norminv(x::Real)::Float64

    n = quantile(Distributions.Normal(), x)

    return n

end

"""
    p2o(p::Float64)

Convert probability to odds.

# Arguments

- `p::Float64`

# Returns

- `o::Float64`
"""
function p2o(p::Real)::Float64

    o = p / (1 - p)

    return o

end

"""
    o2p(o::Float64)

Convert probability to odds.

# Arguments

- `o::Float64`

# Returns

- `p::Float64`
"""
function o2p(o::Real)::Float64

    p = o / (1 + o)

    return p

end
