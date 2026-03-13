export cl2z
export cim
export cimd
export cip
export cir
export cis
export civ

"""
    cl2z(ci; <keyword arguments>)

Convert a confidence level to the corresponding Z-score.

# Arguments

- `cl::Float64`: confidence level; must be in `(0, 1)`
- `twotailed::Bool=true`: if `true`, returns the Z-score for a two-tailed interval, i.e. `z` such that `P(−z ≤ X ≤ z) = cl`; if `false`, returns the Z-score for a one-tailed interval, i.e. `z` such that `P(X ≤ z) = cl`

# Returns

- `Float64`: critical Z-score

# Throws

- `ArgumentError`: if `cl ∉ (0, 1)`

# Notes

- Two-tailed (`twotailed=true`): the CI is `(−z, +z)`
- One-tailed (`twotailed=false`): the CI is `(−∞, z)` (upper bound) or equivalently `(−z, +∞)` (lower bound) depending on direction

# See also

[`cim`](@ref), [`cip`](@ref), [`cir`](@ref)
"""
function cl2z(cl::Float64; twotailed::Bool = true)::Float64

    _bin(cl, (0.0, 1.0), "cl")
    d = Distributions.Normal(0, 1)

    return twotailed ? quantile(d, 1 - (1 - cl) / 2) : quantile(d, cl)

end

"""
    cim(x; <keyword arguments>)

Calculate the confidence interval for the mean.

# Arguments

- `x::AbstractVector`: data vector; must contain at least 2 elements
- `cl::Float64=0.95`: confidence level; must be in `(0, 1)`
- `d::Symbol=:t`: istribution used for the critical value; `:t` (Student's t, recommended for small samples) or `:z` (standard normal)
- `twotailed::Bool=true`: if `true`, compute a two-sided interval

# Returns

- `Tuple{Float64, Float64}`: `(lower_bound, upper_bound)`

# Throws

- `ArgumentError`: if `cl ∉ (0, 1)`, `d ∉ {:t, :z}`, or `length(x) < 2`

# See also

[`cimd`](@ref), [`cis`](@ref), [`civ`](@ref)
"""
function cim(x::AbstractVector; cl::Float64 = 0.95, d::Symbol = :t, twotailed::Bool = true)::Tuple{Float64, Float64}

    _bin(cl, (0.0, 1.0), "cl")
    _check_var(d, [:t, :z], "d")
    @assert length(x) >= 2 "x must contain at least 2 elements."

    n  = length(x)
    m  = mean(x)
    s  = sem(x)
    df = n - 1

    e = if d === :t
        crit_t(df, 1 - cl; twotailed=twotailed) * s
    else
        crit_z(1 - cl; twotailed=twotailed) * s
    end

    return (m - e, m + e)

end

"""
    cimd(x; <keyword arguments>)

Calculate the confidence interval for the median of a 1-D vector.

Uses the order-statistic method: the CI bounds are `x[j]` and `x[k]` where `j` and `k` are quantile-derived indices from the sorted data.

# Arguments

- `x::AbstractVector`: data vector; must contain enough elements so that the derived indices `j` and `k` fall within `[1, n]`
- `cl::Float64=0.95`: confidence level; must be in `(0, 1)`

# Returns

- `Tuple{Float64, Float64}`: `(lower_bound, upper_bound)`

# Throws
- `ArgumentError`: if `cl ∉ (0, 1)` or the sample is too small for the requested confidence level

# See also

[`cim`](@ref), [`cimd(::AbstractArray)`](@ref)
"""
function cimd(x::AbstractVector; cl::Float64 = 0.95)::Tuple{Float64, Float64}

    _bin(cl, (0.0, 1.0), "cl")

    x_sorted = sort(x)
    n  = length(x)
    # median quantile
    q  = 0.5
    z  = cl2z(cl)
    # half-width of the index interval
    hw = z * sqrt(n * q * (1 - q))

    # clamp to avoid index < 1
    j  = max(1, ceil(Int64, n * q - hw))
    # clamp to avoid index > n
    k  = min(n, ceil(Int64, n * q + hw))

    return (x_sorted[j], x_sorted[k])

end

"""
    cimd(x; <keyword arguments>)

Calculate the confidence interval for the median of a multi-dimensional array.

Column medians are computed, sorted, and the order-statistic CI method is applied across columns.

# Arguments

- `x::AbstractArray`: data array; medians are taken along `dims=1` (across rows per column); must have at least 2 columns
- `cl::Float64=0.95`: confidence level; must be in `(0, 1)`

# Returns

- `Tuple{Float64, Float64}`: `(lower_bound, upper_bound)`

# Throws

- `ArgumentError`: if `cl ∉ (0, 1)` or `size(x, 2) < 2`

# See also

[`cimd(::AbstractVector)`](@ref)
"""
function cimd(x::AbstractArray; cl::Float64 = 0.95)::Tuple{Float64, Float64}

    _bin(cl, (0.0, 1.0), "cl")
    @assert size(x, 2) >= 2 "x must have at least 2 columns."

    x_sorted = sort(vec(median(x; dims=1)))
    n = size(x, 2)
    # the quantile of interest; for a median, we will use q = 0.5
    q = 0.5
    z = cl2z(cl)
    hw = z * sqrt(n * q * (1 - q))

    j = max(1, ceil(Int64, n * q - hw))
    k = min(n, ceil(Int64, n * q + hw))

    return (x_sorted[j], x_sorted[k])
end

"""
    cip(p, n; <keyword arguments>)

Calculate the confidence interval for a proportion using the normal approximation (Wald interval).

# Arguments

- `p::Float64`: sample proportion; must be in `[0, 1]`
- `n::Int64`: sample size; must be ≥ 1
- `cl::Float64=0.95`: confidence level; must be in `(0, 1)`

# Returns

- `Tuple{Float64, Float64}`: `(lower_bound, upper_bound)`

# Throws
- `ArgumentError`: if `p ∉ [0, 1]`, `n < 1`, or `cl ∉ (0, 1)`

# See also

[`cim`](@ref), [`cir`](@ref)
"""
function cip(p::Float64, n::Int64; cl::Float64 = 0.95)::Tuple{Float64, Float64}

    _bin(cl, (0.0, 1.0), "cl")
    _in(p, (0.0, 1.0), "p")
    @assert n >= 1 "n must be ≥ 1."

    z= cl2z(cl)
    hw= z * sqrt((p * (1 - p)) / n)
    return (p - hw, p + hw)

end

"""
    cir(x, y; <keyword arguments>)

Calculate the confidence interval for a Pearson correlation coefficient computed from two vectors, using Fisher's Z transformation.

# Arguments

- `x::AbstractVector`: first data vector; must have the same length as `y` and length > 3
- `y::AbstractVector`: second data vector
- `cl::Float64=0.95`: confidence level; must be in `(0, 1)`

# Returns

- `Tuple{Float64, Float64}`: `(lower_bound, upper_bound)`

# Throws

- `ArgumentError`: if lengths differ, `length(x) ≤ 3`, or `cl ∉ (0, 1)`

# See also

[`cir(; r, n, cl)`](@ref), [`cim`](@ref)
"""
function cir(x::AbstractVector, y::AbstractVector; cl::Float64 = 0.95)::Tuple{Float64, Float64}

    _bin(cl, (0.0, 1.0), "cl")
    @assert length(x) == length(y) "x and y must have the same length."
    @assert length(x) > 3 "length(x) must be > 3 for Fisher's Z transform."

    return cir(; r=cor(x, y), n=length(x), cl=cl)

end

"""
    cir(; <keyword arguments>)

Calculate the confidence interval for a Pearson correlation coefficient using Fisher's Z transformation.

Transforms `r` to `z = arctanh(r)`, applies the normal CI, then back-transforms with `tanh`.

# Arguments

- `r::Float64`: Pearson correlation coefficient; must be in `(−1, 1)`
- `n::Int64`: number of observations; must be > 3
- `cl::Float64=0.95`: confidence level; must be in `(0, 1)`

# Returns

- `Tuple{Float64, Float64}`: `(lower_bound, upper_bound)` in correlation units

# Throws

- `ArgumentError`: if `r ∉ (−1, 1)`, `n ≤ 3`, or `cl ∉ (0, 1)`

# See also

[`cir(::AbstractVector, ::AbstractVector)`](@ref)
"""
function cir(; r::Float64, n::Int64, cl::Float64 = 0.95)::Tuple{Float64, Float64}

    _bin(cl, (0.0, 1.0), "cl")
    _in(r, (-1.0, 1.0), "r")
    @assert n > 3 "n must be > 3 for Fisher's Z transform."

    # standard error of Fisher's Z
    se = 1 / sqrt(n - 3)
    # Fisher Z transform: arctanh(r)
    z_score = rfz(r)
    z_crit  = cl2z(cl)

    ci_l = tanh(z_score - z_crit * se)
    ci_u = tanh(z_score + z_crit * se)

    return (ci_l, ci_u)

end

"""
    cis(x; <keyword arguments>)

Calculate the confidence interval for the standard deviation using the chi-squared distribution.

# Arguments

- `x::AbstractVector`: data vector; must contain at least 2 elements
- `cl::Float64=0.95`: confidence level; must be in `(0, 1)`

# Returns

- `Tuple{Float64, Float64}`: `(lower_bound, upper_bound)`

# Throws

- `ArgumentError`: if `cl ∉ (0, 1)` or `length(x) < 2`

# See also

[`civ`](@ref), [`cim`](@ref)
"""
function cis(x::AbstractVector; cl::Float64 = 0.95)::Tuple{Float64, Float64}

    _bin(cl, (0.0, 1.0), "cl")
    @assert length(x) >= 2 "x must contain at least 2 elements."

    α = 1 - cl
    s = std(x)
    n = length(x)
    df = n - 1

    # upper chi² critical value (lower SD bound)
    chi_l = crit_chi(df, 1 - α / 2)
    # lower chi² critical value (upper SD bound)
    chi_u = crit_chi(df, α / 2)

    return (sqrt((df * s^2) / chi_l), sqrt((df * s^2) / chi_u))

end

"""
    civ(x; <keyword arguments>)

Calculate the confidence interval for the variance using the chi-squared distribution.

# Arguments

- `x::AbstractVector`: data vector; must contain at least 2 elements
- `cl::Float64=0.95`: confidence level; must be in `(0, 1)`

# Returns

- `Tuple{Float64, Float64}`: `(lower_bound, upper_bound)`

# Throws

- `ArgumentError`: if `cl ∉ (0, 1)` or `length(x) < 2`

# See also

[`cis`](@ref), [`cim`](@ref)
"""
function civ(x::AbstractVector; cl::Float64 = 0.95)::Tuple{Float64, Float64}

    _bin(cl, (0.0, 1.0), "cl")
    @assert length(x) >= 2 "x must contain at least 2 elements."

    α = 1 - cl
    v = var(x)
    n = length(x)
    df = n - 1

    # upper chi² critical value (lower variance bound)
    chi_l = crit_chi(df, 1 - α / 2)
    # lower chi² critical value (upper variance bound)
    chi_u = crit_chi(df, α / 2)

    return ((df * v) / chi_l, (df * v) / chi_u)

end
