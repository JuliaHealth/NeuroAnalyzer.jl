export cl2z
export cim
export cimd
export cip
export cir
export cis
export civ

"""
    cl2z(ci; <keyword arguments>)

Convert confidence level to z score.

# Arguments

- `cl::Float64`: confidence level
- `twotailed::Bool=true`: one- or two-tailed probability

# Returns

- `z::Float64`

# Notes

The confidence interval is (-z, +z) if `twotailed=true`; otherwise it is (-∞, -z) on the left and (z, +∞) on the right.
"""
function cl2z(cl::Float64; twotailed::Bool=true)::Float64

    _bin(cl, (0.0, 1.0), "cl")
    
    d = Distributions.Normal(0, 1)
    if twotailed
        z = quantile(d, 1 - ((1 - cl) / 2))
    else
        z = quantile(d, cl)
    end

    return z

end

"""
    cim(x; <keyword arguments>)

Calculate confidence interval for the mean.

# Arguments

- `x::AbstractVector`
- `cl::Float64=0.95`: confidence level
- `d::Symbol=:t`: distribution used for critical value calculation (`:z` or `:t`)
- `twotailed::Bool=true`: interval type, use `twotailed=false` to calculate lower and upper bound separately (both will be returned)

# Returns

- `cim::Tuple{Float64, Float64}`: lower and upper bound
"""
function cim(x::AbstractVector; cl::Float64=0.95, d::Symbol=:t, twotailed::Bool=true)::Tuple{Float64, Float64}

    _bin(cl, (0.0, 1.0), "cl")
    _check_var(d, [:t, :z], "d")

    n = length(x)
    m = mean(x)
    df = n - 1
    s = sem(x)

    if d === :t
        tc = crit_t(df, 1-cl, twotailed=twotailed)
        e = tc * s
    else
        zc = crit_z(1 - cl, twotailed=twotailed)
        e = zc * s
    end

    return (m - e, m + e)

end

"""
    cimd(x; <keyword arguments>)

Calculate confidence interval for the median.

# Arguments

- `x::AbstractVector`
- `cl::Float64=0.95`: confidence level

# Returns

- `cimd::Tuple{Float64, Float64}`
"""
function cimd(x::AbstractVector; cl::Float64=0.95)::Tuple{Float64, Float64}

    _bin(cl, (0.0, 1.0), "cl")

    x_new = sort(x)
    n = length(x)
    q = 0.5 # the quantile of interest; for a median, we will use q = 0.5
    z = cl2z(cl)
    j = ceil(Int64, (n * q) - (z * sqrt((n * q) * (1 - q))))
    k = ceil(Int64, (n * q) + (z * sqrt((n * q) * (1 - q))))

    return (x_new[j], x_new[k])

end

"""
    cimd(x; <keyword arguments>)

Calculate confidence interval for the median.

# Arguments

- `x::AbstractArray`
- `cl::Float64=0.95`: confidence level

# Returns

- `cimd::Tuple{Float64, Float64}`
"""
function cimd(x::AbstractArray; cl::Float64=0.95)::Tuple{Float64, Float64}

    _bin(cl, (0.0, 1.0), "cl")

    x_new = sort(vec(median(x, dims=1)))
    n = size(x, 2)
    q = 0.5 # the quantile of interest; for a median, we will use q = 0.5
    z = cl2z(cl)
    j = ceil(Int64, (n * q) - (z * sqrt((n * q) * (1 - q))))
    k = ceil(Int64, (n * q) + (z * sqrt((n * q) * (1 - q))))

    return (x_new[j], x_new[k])

end

"""
    cip(p, n; <keyword arguments>)

Calculate confidence interval for the proportion.

# Arguments

- `p::Float64`: proportion
- `n::Int64`: sample size
- `cl::Float64=0.95`: confidence level

# Returns

- `cip::Tuple{Float64, Float64}`
"""
function cip(p::Float64, n::Int64; cl::Float64=0.95)::Tuple{Float64, Float64}

    _bin(cl, (0.0, 1.0), "cl")

    z = cl2z(cl)
    ci_l = (p - z * sqrt((p * (1 - p)) / n))
    ci_u = (p + z * sqrt((p * (1 - p)) / n))

    return (ci_l, ci_u)

end

"""
    cir(x, y; <keyword arguments>)

Calculate confidence interval for the correlation coefficient.

# Arguments

- `x::AbstractVector`
- `y::AbstractVector`
- `cl::Float64=0.95`: confidence level

# Returns

- `cir::Tuple{Float64, Float64}`
"""
function cir(x::AbstractVector, y::AbstractVector; cl::Float64=0.95)::Tuple{Float64, Float64}

    _bin(cl, (0.0, 1.0), "cl")
    @assert length(x) == length(y) "Lengths of x and y must be equal."
    @assert length(x) > 3 "Lengths of x and y must be > 3."

    n = length(x)
    r = cor(x, y)

    return cir(r=r, n=n, cl=cl)

end

"""
    cir(; <keyword arguments>)

Calculate confidence interval for the correlation coefficient.

# Arguments

- `r::Float64`: correlation coefficient
- `n::Int64`: number of observations
- `cl::Float64=0.95`: confidence level

# Returns

- `cir::Tuple{Float64, Float64}`
"""
function cir(; r::Float64, n::Int64, cl::Float64=0.95)::Tuple{Float64, Float64}

    _bin(cl, (0.0, 1.0), "cl")
    _in(r, (-1.0, 1.0), "r")
    @assert n > 0 "n must be > 0."

    z_r = 1 / sqrt(n - 3)
    z_score = rfz(r)
    ci_l = tanh(z_score - (z_r * cl2z(cl)))
    ci_u = tanh(z_score + (z_r * cl2z(cl)))

    return (ci_l, ci_u)

end

"""
    cis(x; <keyword arguments>)

Calculate confidence interval for the standard deviation.

# Arguments

- `x::AbstractVector`
- `cl::Float64=0.95`: confidence level

# Returns

- `cis::Tuple{Float64, Float64}`: lower and upper bound
"""
function cis(x::AbstractVector; cl::Float64=0.95)::Tuple{Float64, Float64}

    _bin(cl, (0.0, 1.0), "cl")

    α = 1 - cl
    s = std(x)
    n = length(x)
    df = n - 1

    chi1_crit = crit_chi(df, 1-α/2)
    chi2_crit = crit_chi(df, α/2)

    return (sqrt(((n - 1) * s^2)/chi1_crit), sqrt(((n - 1) * s^2)/chi2_crit))

end

"""
    civ(x; <keyword arguments>)

Calculate confidence interval for the variance.

# Arguments

- `x::AbstractVector`
- `cl::Float64=0.95`: confidence level

# Returns

- `civ::Tuple{Float64, Float64}`: lower and upper bound
"""
function civ(x::AbstractVector; cl::Float64=0.95)::Tuple{Float64, Float64}

    _bin(cl, (0.0, 1.0), "cl")

    α = 1 - cl
    v = var(x)
    n = length(x)
    df = n - 1

    chi1_crit = crit_chi(df, 1-α/2)
    chi2_crit = crit_chi(df, α/2)

    return (((n - 1) * v)/chi1_crit, ((n - 1) * v)/chi2_crit)

end
