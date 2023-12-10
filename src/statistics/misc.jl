export z_score
export k_categories
export se
export rng
export moe
export cvar
export effsize
export binom_prob
export binom_stat
export cvar_mean
export cvar_median
export ci_median
export ci_r
export ci2z
export r_test
export slope
export distance

"""
    z_score(x)

Calculate Z-scores for each value of the vector `x`.

# Arguments

- `x::AbstractVector`

# Returns

- `z_score::Vector{Float64}`
"""
function z_score(x::AbstractVector)

    return (x .- mean(x)) ./ std(x)

end

"""
    k_categories(n)

Calculate number of categories for a given sample size `n`.

# Arguments

- `n::Int64`

# Returns

Named tuple containing:
- `k1::Float64`: sqrt(n)
- `k2::Float64`: 1 + 3.222 * log10(n)
"""
function k_categories(n::Int64)

    return (k1=sqrt(n), k2=(1 + 3.222 * log10(n)))

end

"""
    se(x)

Calculate standard error.

# Arguments

- `x::AbstractVector`

# Returns

- `se::Float64`
"""
function se(x::AbstractVector)

    return std(x) / sqrt(length(x))

end

"""
    rng(x)

Calculate range.

# Arguments

- `x::AbstractVector`

# Returns

- `rng::Float64`
"""
function rng(x::AbstractVector)

    return maximum(x) - minimum(x)

end

"""
    rng(x)

Calculate range.

# Arguments

- `x::AbstractArray`

# Returns

- `rng::Float64`
"""
function rng(x::AbstractArray)

    return maximum(x) - minimum(x)

end

"""
    moe(n)

Calculate margin of error for given sample size `n`.

# Arguments

- `n::Int64`

# Returns

- `moe::Float64`
"""
function moe(n::Int64)

    return 1 / sqrt(n)

end

"""
    cvar(se, s)

Calculate coefficient of variation for statistic `s`.

# Arguments

- `se::Real`: standard error
- `s::Real`: statistics, e.g. mean value

# Returns

- `cvar::Float64`
"""
function cvar(se::Real, s::Real)

    return 100 * (se / s)

end

"""
    effsize(p1, p2)

Calculate effect size for two proportions `p1` and `p2`.

# Arguments

- `p1::Float64`: 1st proportion, e.g. 0.7
- `p2::Float64`: 2nd proportion, e.g. 0.3

# Returns

- `e::Float64`
"""
function effsize(p1::Float64, p2::Float64)

    @assert p1 + p2 == 1.0 "Proportions must add to 1.0."

    e = 2 * asin(sqrt(p1)) - 2 * asin(sqrt(p2))

    return e

end

"""
    binom_prob(p, r, n)

Calculate probability of exactly `r` successes in `n` trials.

# Arguments

- `p::Float64`: proportion of successes
- `r::Int64`: number of successes
- `n::Int64`: number of trials

# Returns

- `binom_prob::Float64`: probability
"""
function binom_prob(p::Float64, r::Int64, n::Int64)

    return binomial(n, r) * (p^r) * (1 - p)^(n - r)

end

"""
    binom_stat(p, n)

Calculate mean and standard deviation for probability `p`.

# Arguments

- `p::Float64`: proportion of successes
- `n::Int64`: number of trials

# Returns

Named tuple containing:
- `m::Float64`: mean
- `s::Float64`: standard deviation
"""
function binom_stat(p::Float64, n::Int64)

    m = n * p
    s = sqrt(n * p * (1 - p))

    return (m=m, s=s)

end

"""
    cvar_mean(x)

Calculate coefficient of variation for a mean.

# Arguments

- `x::AbstractVector`

# Returns

- `cvar_mean::Float64`
"""
function cvar_mean(x::AbstractVector)

    return std(x) / mean(x)

end

"""
    cvar_median(x)

Calculate coefficient of variation for a median.

# Arguments

- `x::AbstractVector`

# Returns

- `cvar_median::Float64`
"""
function cvar_median(x::AbstractVector)

    return ((quantile(x, 0.75) - quantile(x, 0.25)) / 2) / median(x)
    
end

"""
    ci_median(x; ci_level)

Calculate confidence interval for a median.

# Arguments

- `x::AbstractVector`
- `ci_level::Float64=0.95`: confidence level

# Returns

- `ci_median::Tuple(Float64, Float64)`
"""
function ci_median(x::AbstractVector; ci_level::Float64=0.95)

    x_new = sort(x)
    n = length(x)
    q = 0.5 # the quantile of interest; for a median, we will use q = 0.5
    z = ci2z(ci_level)

    j = ceil(Int64, (n * q) - (z * sqrt((n * q) * (1 - q))))
    k = ceil(Int64, (n * q) + (z * sqrt((n * q) * (1 - q))))

    return (x_new[j], x_new[k])

end

"""
    ci_r(x, y; ci_level)

Calculate confidence interval for a correlation coefficient.

# Arguments

- `x::AbstractVector`
- `y::AbstractVector`
- `ci_level::Float64=0.95`: confidence level

# Returns

- `ci_r::Tuple(Float64, Float64)`
"""
function ci_r(x::AbstractVector, y::AbstractVector; ci_level::Float64=0.95)

    @assert length(x) == length(y) "Both vectors must have the same length."

    n = length(x)
    r = cor(x, y)

    z_r = 1 / sqrt(n - 3)

    r_values = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99]
    z_values = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.121, 0.131, 0.141, 0.151, 0.161, 0.172, 0.182, 0.192, 0.203, 0.213, 0.224, 0.234, 0.245, 0.255, 0.266, 0.277, 0.288, 0.299, 0.31, 0.321, 0.332, 0.343, 0.354, 0.365, 0.377, 0.389, 0.4, 0.412, 0.424, 0.436, 0.448, 0.46, 0.472, 0.485, 0.497, 0.51, 0.523, 0.536, 0.549, 0.563, 0.577, 0.59, 0.604, 0.618, 0.633, 0.648, 0.663, 0.678, 0.693, 0.709, 0.725, 0.741, 0.758, 0.775, 0.793, 0.811, 0.829, 0.848, 0.867, 0.887, 0.908, 0.929, 0.95, 0.973, 0.996, 1.02, 1.045, 1.071, 1.099, 1.127, 1.157, 1.188, 1.221, 1.256, 1.293, 1.333, 1.376, 1.422, 1.472, 1.528, 1.589, 1.658, 1.738, 1.832, 1.946, 2.092, 2.298, 2.647]

    r_idx = vsearch(r, r_values)
    z_score = z_values[r_idx]

    ci_h = z_score + z_r * ci2z(ci_level)
    ci_l = z_score - z_r * ci2z(ci_level)

    ci_h_idx = vsearch(ci_h, z_values)
    ci_l_idx = vsearch(ci_l, z_values)

    return (r_values[ci_l_idx], r_values[ci_h_idx])

end

"""
    ci_r(; r, n, ci_level)

Calculate confidence interval for a correlation coefficient.

# Arguments

- `r::Float64`
- `n::Int64`: number of observations
- `ci_level::Float64=0.95`: confidence level

# Returns

- `ci_r::Tuple(Float64, Float64)`
"""
function ci_r(; r::Float64, n::Int64, ci_level::Float64=0.95)

    z_r = 1 / sqrt(n - 3)

    r_values = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99]
    z_values = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.121, 0.131, 0.141, 0.151, 0.161, 0.172, 0.182, 0.192, 0.203, 0.213, 0.224, 0.234, 0.245, 0.255, 0.266, 0.277, 0.288, 0.299, 0.31, 0.321, 0.332, 0.343, 0.354, 0.365, 0.377, 0.389, 0.4, 0.412, 0.424, 0.436, 0.448, 0.46, 0.472, 0.485, 0.497, 0.51, 0.523, 0.536, 0.549, 0.563, 0.577, 0.59, 0.604, 0.618, 0.633, 0.648, 0.663, 0.678, 0.693, 0.709, 0.725, 0.741, 0.758, 0.775, 0.793, 0.811, 0.829, 0.848, 0.867, 0.887, 0.908, 0.929, 0.95, 0.973, 0.996, 1.02, 1.045, 1.071, 1.099, 1.127, 1.157, 1.188, 1.221, 1.256, 1.293, 1.333, 1.376, 1.422, 1.472, 1.528, 1.589, 1.658, 1.738, 1.832, 1.946, 2.092, 2.298, 2.647]

    r_idx = vsearch(r, r_values)
    z_score = z_values[r_idx]

    ci_h = z_score + z_r * ci2z(ci_level)
    ci_l = z_score - z_r * ci2z(ci_level)

    ci_h_idx = vsearch(ci_h, z_values)
    ci_l_idx = vsearch(ci_l, z_values)

    return (r_values[ci_l_idx], r_values[ci_h_idx])

end

"""
    ci2z(ci_level)

Convert Confidence Interval level to z score.

# Arguments

- `ci_level::Float64=0.95`: confidence level

# Returns

- `z_score::Float64`
"""
function ci2z(ci_level::Float64=0.95)

    return quantile(Distributions.Normal(0.0, 1.0), 1 - (1 - ci_level) / 2)

end

"""
    r_test(; r1, r2, n1, n2)

Test if two correlation coefficients are significantly different.

# Arguments

- `r1::Float64`: correlation coefficient, group 1
- `r2::Float64`: correlation coefficient, group 2
- `n1::Int64`: number of observations, group 1
- `n2::Int64`: number of observations, group 2

# Returns

- `z_r1r2::Float64`
"""
function r_test(; r1::Float64, r2::Float64, n1::Int64, n2::Int64)

    r_values = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99]
    z_values = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.121, 0.131, 0.141, 0.151, 0.161, 0.172, 0.182, 0.192, 0.203, 0.213, 0.224, 0.234, 0.245, 0.255, 0.266, 0.277, 0.288, 0.299, 0.31, 0.321, 0.332, 0.343, 0.354, 0.365, 0.377, 0.389, 0.4, 0.412, 0.424, 0.436, 0.448, 0.46, 0.472, 0.485, 0.497, 0.51, 0.523, 0.536, 0.549, 0.563, 0.577, 0.59, 0.604, 0.618, 0.633, 0.648, 0.663, 0.678, 0.693, 0.709, 0.725, 0.741, 0.758, 0.775, 0.793, 0.811, 0.829, 0.848, 0.867, 0.887, 0.908, 0.929, 0.95, 0.973, 0.996, 1.02, 1.045, 1.071, 1.099, 1.127, 1.157, 1.188, 1.221, 1.256, 1.293, 1.333, 1.376, 1.422, 1.472, 1.528, 1.589, 1.658, 1.738, 1.832, 1.946, 2.092, 2.298, 2.647]

    r1_idx = vsearch(r1, r_values)
    z1_score = z_values[r1_idx]

    r2_idx = vsearch(r2, r_values)
    z2_score = z_values[r2_idx]

    z_r1r2 = (z1_score - z2_score) / sqrt((1 / (n1 - 3)) + (1 / (n2 - 3)))

    return z_r1r2

end

slope(p1::Tuple{Real, Real}, p2::Tuple{Real, Real}) = (p2[2] - p1[2]) / (p2[1] - p1[1])
distance(p1::Tuple{Real, Real}, p2::Tuple{Real, Real}) = sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)