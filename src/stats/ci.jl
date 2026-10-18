export cl2z
export cim
export cimd
export cip
export cir
export civ

"""
    cl2z(ci; <keyword arguments>)

Calculate Z-score of confidence level.

# Arguments

- `cl::Float64`: confidence level; must be in `(0, 1)`
- `twotailed::Bool=true`: if `true`, returns the Z-score for a two-tailed interval, i.e. `z` such that `P(−z ≤ X ≤ z) = cl`; if `false`, returns the Z-score for a one-tailed interval, i.e. `z` such that `P(X ≤ z) = cl`

# Returns

- `Float64`: corresponding Z-score
"""
function cl2z(cl::Float64; twotailed::Bool = true)::Float64
    # validate
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
- `d::Symbol=:t`: distribution used for the critical value; `:t` (Student's t, recommended for small samples) or `:z` (standard normal)
- `twotailed::Bool=true`: if `true`, compute a two-sided interval; if `false`, a one-sided interval (see `tail`)
- `tail::Symbol=:lower`: one-sided interval (`twotailed=false`):
    - `:lower`: lower confidence bound, interval `(m − e, Inf)`
    - `:upper`: upper confidence bound, interval `(-Inf, m + e)`

# Returns

- `Tuple{Float64, Float64}`: `(lower_bound, upper_bound)`; the interval covers the mean with probability `cl`
"""
function cim(
    x::AbstractVector;
    cl::Float64 = 0.95,
    d::Symbol = :t,
    twotailed::Bool = true,
    tail::Symbol = :lower,
)::Tuple{Float64, Float64}
    # validate
    _bin(cl, (0.0, 1.0), "cl")
    _check_var(d, [:t, :z], "d")
    _check_var(tail, [:lower, :upper], "tail")
    length(x) >= 2 || throw(ArgumentError("x must contain at least 2 elements."))

    n = length(x)
    m = mean(x)
    s = sem(x)
    df = n - 1

    e = if d === :t
        crit_t(df, 1 - cl; twotailed = twotailed) * s
    else
        crit_z(1 - cl; twotailed = twotailed) * s
    end

    twotailed && return (m - e, m + e)

    return tail === :lower ? (m - e, Inf) : (-Inf, m + e)
end

# order-statistic confidence interval for the median of n sorted values:
# ranks n/2 − z√n/2 (rounded down) and 1 + n/2 + z√n/2 (rounded up)
function _cimd_sorted(x_sorted::AbstractVector, cl::Float64)::Tuple{Float64, Float64}
    n = length(x_sorted)
    q = 0.5
    hw = cl2z(cl) * sqrt(n * q * (1 - q))
    j = clamp(floor(Int64, n * q - hw), 1, n)
    k = clamp(ceil(Int64, 1 + n * q + hw), 1, n)
    return (Float64(x_sorted[j]), Float64(x_sorted[k]))
end

"""
    cimd(x; <keyword arguments>)

Calculate the confidence interval for the median of a 1-D vector.

Uses the order-statistic (binomial, normal approximation) method: the bounds are the sorted values of rank `⌊n/2 − z√n/2⌋` and `⌈1 + n/2 + z√n/2⌉`, clamped to `[1, n]`.

# Arguments

- `x::AbstractVector`: data vector; must not be empty
- `cl::Float64=0.95`: confidence level; must be in `(0, 1)`

# Returns

- `Tuple{Float64, Float64}`: `(lower_bound, upper_bound)`
"""
function cimd(x::AbstractVector; cl::Float64 = 0.95)::Tuple{Float64, Float64}
    # validate
    _bin(cl, (0.0, 1.0), "cl")
    isempty(x) && throw(ArgumentError("x must not be empty."))

    return _cimd_sorted(sort(x), cl)
end

"""
    cimd(x; <keyword arguments>)

Calculate the confidence interval for the median of a multi-dimensional array.

Medians are taken along the first dimension; the order-statistic CI method (see the vector method) is applied to all resulting medians.

# Arguments

- `x::AbstractArray`: data array; medians are taken along `dims=1`; must yield at least 2 medians
- `cl::Float64=0.95`: confidence level; must be in `(0, 1)`

# Returns

- `Tuple{Float64, Float64}`: `(lower_bound, upper_bound)`
"""
function cimd(x::AbstractArray; cl::Float64 = 0.95)::Tuple{Float64, Float64}
    # validate
    _bin(cl, (0.0, 1.0), "cl")
    md = vec(median(x; dims = 1))
    length(md) >= 2 || throw(ArgumentError("x must have at least 2 columns."))

    return _cimd_sorted(sort(md), cl)
end

"""
    cip(p, n; <keyword arguments>)

Calculate the confidence interval for a proportion using the normal approximation (Wald interval), truncated to `[0, 1]`.

# Arguments

- `p::Float64`: sample proportion; must be in `[0, 1]`
- `n::Int64`: sample size; must be ≥ 1
- `cl::Float64=0.95`: confidence level; must be in `(0, 1)`

# Returns

- `Tuple{Float64, Float64}`: `(lower_bound, upper_bound)`
"""
function cip(p::Float64, n::Int64; cl::Float64 = 0.95)::Tuple{Float64, Float64}
    # validate
    _bin(cl, (0.0, 1.0), "cl")
    _in(p, (0.0, 1.0), "p")
    n >= 1 || throw(ArgumentError("n must be ≥ 1."))

    z = cl2z(cl)
    hw = z * sqrt((p * (1 - p)) / n)

    return (max(0.0, p - hw), min(1.0, p + hw))
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
"""
function cir(
    x::AbstractVector,
    y::AbstractVector;
    cl::Float64 = 0.95,
)::Tuple{Float64, Float64}
    # validate
    _bin(cl, (0.0, 1.0), "cl")
    length(x) == length(y) || throw(ArgumentError("x and y must have the same length."))
    length(x) > 3 || throw(ArgumentError("length(x) must be > 3 for Fisher's Z transform."))

    return cir(; r = cor(x, y), n = length(x), cl = cl)
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
"""
function cir(; r::Float64, n::Int64, cl::Float64 = 0.95)::Tuple{Float64, Float64}
    # validate
    _bin(cl, (0.0, 1.0), "cl")
    _in(r, (-1.0, 1.0), "r")
    n > 3 || throw(ArgumentError("n must be > 3 for Fisher's Z transform."))

    # standard error of Fisher's Z
    se = 1 / sqrt(n - 3)
    # Fisher Z transform: arctanh(r)
    z_score = rfz(r)
    z_crit = cl2z(cl)

    ll = tanh(z_score - z_crit * se)
    ul = tanh(z_score + z_crit * se)

    return (ll, ul)
end

"""
    cis(x; <keyword arguments>)

Calculate the confidence interval for the standard deviation using the chi-squared distribution.

This is a method of `Base.cis()` (for a number, `cis(x)` is still `exp(im * x)`).

# Arguments

- `x::AbstractVector`: data vector; must contain at least 2 elements
- `cl::Float64=0.95`: confidence level; must be in `(0, 1)`

# Returns

- `Tuple{Float64, Float64}`: `(lower_bound, upper_bound)`
"""
function Base.cis(x::AbstractVector; cl::Float64 = 0.95)::Tuple{Float64, Float64}
    # validate
    _bin(cl, (0.0, 1.0), "cl")
    length(x) >= 2 || throw(ArgumentError("x must contain at least 2 elements."))

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
"""
function civ(x::AbstractVector; cl::Float64 = 0.95)::Tuple{Float64, Float64}
    # validate
    _bin(cl, (0.0, 1.0), "cl")
    length(x) >= 2 || throw(ArgumentError("x must contain at least 2 elements."))

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
