export crit_z
export crit_t
export crit_chi

"""
    crit_z(alpha; <keyword arguments>)

Calculate the critical Z-score for a given significance level.

# Arguments

- `alpha::Float64=0.05`: significance level (upper or lower tail probability); must be in `(0, 1)`
- `twotailed::Bool=true`: if `true`, compute the two-tailed critical value; if `false`, compute the one-tailed critical value

# Returns

- `Float64`: critical Z-score (always positive; the rejection region is symmetric around zero)

# Throws

- `ArgumentError`: if `alpha ∉ (0, 1)`

# Notes

Critical regions:

- One-tailed left:  `(−∞, −z]`
- One-tailed right: `[z, +∞)`
- Two-tailed:       `(−∞, −z] ∪ [z, +∞)`

# See also

[`crit_t`](@ref), [`crit_chi`](@ref), [`cl2z`](@ref)
"""
function crit_z(alpha::Float64 = 0.05; twotailed::Bool = true)::Float64

    !(alpha > 0.0) && throw(ArgumentError("alpha must be > 0."))
    !(alpha < 1.0) && throw(ArgumentError("alpha must be < 1."))

    return cl2z(1 - alpha; twotailed=twotailed)

end

"""
    crit_t(df, alpha; <keyword arguments>)

Calculate the critical t-value for a given degrees of freedom and significance level.

# Arguments

- `df::Real`: degrees of freedom; must be > 0; typically `df = n − 1`
- `alpha::Float64=0.05`: significance level (upper or lower tail probability); must be in `(0, 1)`
- `twotailed::Bool=true`: if `true`, compute the two-tailed critical value; if `false`, compute the one-tailed critical value

# Returns

- `Float64`: critical t-value (always positive)

# Throws

- `ArgumentError`: if `alpha ∉ (0, 1)` or `df ≤ 0`

# Notes

Critical regions:

- One-tailed left:  `(−∞, −t]`
- One-tailed right: `[t, +∞)`
- Two-tailed:       `(−∞, −t] ∪ [t, +∞)`

# See also

[`crit_z`](@ref), [`crit_chi`](@ref)
"""
function crit_t(df::Real, alpha::Float64 = 0.05; twotailed::Bool = true)::Float64

    !(alpha > 0.0) && throw(ArgumentError("alpha must be > 0."))
    !(alpha < 1.0) && throw(ArgumentError("alpha must be < 1."))
    !(df > 0  ) && throw(ArgumentError("df must be > 0."))

    return twotailed ? quantile(TDist(df), 1 - alpha / 2) :
                       quantile(TDist(df), 1 - alpha)

end

"""
    crit_chi(df, alpha)

Calculate the critical χ² value for a given degrees of freedom and significance level.

# Arguments

- `df::Real`: degrees of freedom; must be > 0; typically `df = n − 1`
- `alpha::Float64=0.05`: significance level (upper or lower tail probability); must be in `(0, 1)`

# Returns

- `Float64`: critical χ² value such that `P(X ≤ chi) = alpha` under `χ²(df)`

# Throws

- `ArgumentError`: if `alpha ∉ (0, 1)` or `df ≤ 0`

# Notes

To obtain the upper-tail critical value (i.e. `P(X > chi) = alpha`) pass `1 − alpha` as the `alpha` argument.

# See also

[`crit_z`](@ref), [`crit_t`](@ref)
"""
function crit_chi(df::Real, alpha::Float64 = 0.05)::Float64

    !(alpha > 0.0) && throw(ArgumentError("alpha must be > 0."))
    !(alpha < 1.0) && throw(ArgumentError("alpha must be < 1."))
    !(df > 0  ) && throw(ArgumentError("df must be > 0."))

    return quantile(Distributions.Chisq(df), alpha)

end
