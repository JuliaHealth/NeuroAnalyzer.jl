export crit_z
export crit_t
export crit_chi

"""
    crit_z(alpha; <keyword arguments>)

Calculate critical z score.

# Arguments

- `alpha::Float64=0.05`: alpha value
- `twotailed::Bool=true`: one- or two-tailed probability

# Returns

- `z::Float64`

# Notes

Critical region for one- tailed probability:
- left: `(-∞ , -z]`
- right: `[z , ∞)`

Critical region for two-tailed probability: `(-∞ , -z] ∪ [z, ∞)`
"""
function crit_z(alpha::Float64=0.05; twotailed::Bool=true)::Float64

    @assert alpha > 0.0 "alpha must be > 0.0."
    @assert alpha < 1.0 "alpha must be < 1.0."

    z = cl2z(1 - alpha, twotailed=twotailed)

    return z

end

"""
    crit_t(df, alpha; <keyword arguments>)

Calculate critical t value.

# Arguments

- `df::Real`: degrees of freedom (usually df = n - 1)
- `alpha::Float64=0.05`: alpha value
- `twotailed::Bool=true`: one- or two-tailed probability

# Returns

- `t::Float64`

# Notes

Critical region for one- tailed probability:
- left: `(-∞ , -t]`
- right: `[t , ∞)`

Critical region for two-tailed probability: `(-∞ , -t] ∪ [t, ∞)`
"""
function crit_t(df::Real, alpha::Float64=0.05; twotailed::Bool=true)::Float64

    @assert alpha > 0.0 "alpha must be > 0.0."
    @assert alpha < 1.0 "alpha must be < 1.0."

    if twotailed
        t = quantile(TDist(df), 1 - (alpha / 2))
    else
        t = quantile(TDist(df), 1 - alpha)
    end

    return t

end

"""
    crit_chi(df, alpha)

Calculate critical Χ² score.

# Arguments

- `df::Real`: degrees of freedom (usually df = n - 1)
- `alpha::Float64=0.05`: alpha value

# Returns

- `chi::Float64`
"""
function crit_chi(df::Real, alpha::Float64=0.05)::Float64

    @assert alpha > 0.0 "alpha must be > 0.0."
    @assert alpha < 1.0 "alpha must be < 1.0."

    chi = quantile(Distributions.Chisq(df), alpha)

    return chi

end
