export z_score
export k_categories
export sem
export rng
export moe
export binom_prob
export binom_stat
export cvar_mean
export cvar_median
export ci_median
export ci_prop
export ci_r
export r1r2_test
export slope
export distance
export count_thresh
export crit_t
export ci2z
export p2z
export z2p
export cmp_stat
export fwhm
export cosine_similarity
export std
export permute

"""
    z_score(x)

Calculate Z-scores for each value of the vector `x`.

# Arguments

- `x::AbstractVector`

# Returns

- `z::Vector{Float64}`
"""
function z_score(x::AbstractVector)

    m = mean(x)
    s = std(x)
    z = (x .- m) ./ s

    return z

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

    k1 = sqrt(n)
    k2 = 1 + 3.222 * log10(n)

    return (k1=k1, k2=k2)

end

"""
    sem(x)

Calculate standard error of the mean.

# Arguments

- `x::AbstractVector`

# Returns

- `s::Float64`
"""
function sem(x::AbstractVector)

    s = std(x) / sqrt(length(x))

    return s

end

"""
    rng(x)

Calculate range.

# Arguments

- `x::AbstractArray`

# Returns

- `r::Float64`
"""
function rng(x::AbstractArray)

    r = maximum(x) - minimum(x)

    return r

end

"""
    moe(n)

Calculate margin of error for given sample size `n`.

# Arguments

- `n::Int64`

# Returns

- `m::Float64`
"""
function moe(n::Int64)

    m = 1 / sqrt(n)

    return m

end

"""
    moe(x)

Calculate margin of error.

# Arguments

- `x::AbstractArray`

# Returns

- `m::Float64`
"""
function moe(x::AbstractArray)

    n = length(x)
    m = 1 / sqrt(n)

    return m

end

"""
    binom_prob(p, r, n)

Calculate probability of exactly `r` successes in `n` trials.

# Arguments

- `p::Float64`: proportion of successes
- `r::Int64`: number of successes
- `n::Int64`: number of trials

# Returns

- `bp::Float64`: probability
"""
function binom_prob(p::Float64, r::Int64, n::Int64)

    bp = binomial(n, r) * (p^r) * (1 - p)^(n - r)

    return bp

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
    ci_median(x; <keyword arguments>)

Calculate confidence interval for a median.

# Arguments

- `x::AbstractVector`
- `ci::Float64=0.95`: confidence level

# Returns

- `ci_median::Tuple(Float64, Float64)`
"""
function ci_median(x::AbstractVector; ci::Float64=0.95)

    x_new = sort(x)
    n = length(x)
    q = 0.5 # the quantile of interest; for a median, we will use q = 0.5
    z = ci2z(ci)

    j = ceil(Int64, (n * q) - (z * sqrt((n * q) * (1 - q))))
    k = ceil(Int64, (n * q) + (z * sqrt((n * q) * (1 - q))))

    return (x_new[j], x_new[k])

end

"""
    ci_median(x; <keyword arguments>)

Calculate confidence interval for a median.

# Arguments

- `x::AbstractArray`
- `ci::Float64=0.95`: confidence level

# Returns

- `ci_median::Tuple(Float64, Float64)`
"""
function ci_median(x::AbstractArray; ci::Float64=0.95)

    x_new = sort(vec(median(x, dims=1)))
    n = size(x, 2)
    q = 0.5 # the quantile of interest; for a median, we will use q = 0.5
    z = ci2z(ci)

    j = ceil(Int64, (n * q) - (z * sqrt((n * q) * (1 - q))))
    k = ceil(Int64, (n * q) + (z * sqrt((n * q) * (1 - q))))

    return (x_new[j], x_new[k])

end

"""
    ci_prop(p, n; <keyword arguments>)

Calculate confidence interval for a proportion.

# Arguments

- `p::Float64`: proportion
- `n::Int64`: sample size
- `ci::Float64=0.95`: confidence level

# Returns

- `ci_prop::Tuple(Float64, Float64)`
"""
function ci_prop(p::Float64, n::Int64; ci::Float64=0.95)

    z = ci2z(ci)
    q = 1 - p
    ci_l = (p - z * sqrt((p * q) / n))
    ci_u = (p + z * sqrt((p * q) / n))

    return (ci_l, ci_u)

end

"""
    ci_r(x, y; <keyword arguments>)

Calculate confidence interval for a correlation coefficient.

# Arguments

- `x::AbstractVector`
- `y::AbstractVector`
- `ci::Float64=0.95`: confidence level

# Returns

- `ci_r::Tuple(Float64, Float64)`
"""
function ci_r(x::AbstractVector, y::AbstractVector; ci::Float64=0.95)

    @assert length(x) == length(y) "Both vectors must have the same length."
    @assert length(x) > 3 "Length of both vectors must be > 3."

    n = length(x)
    r = cor(x, y)

    z_r = 1 / sqrt(n - 3)

    r_values = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99]
    z_values = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.121, 0.131, 0.141, 0.151, 0.161, 0.172, 0.182, 0.192, 0.203, 0.213, 0.224, 0.234, 0.245, 0.255, 0.266, 0.277, 0.288, 0.299, 0.31, 0.321, 0.332, 0.343, 0.354, 0.365, 0.377, 0.389, 0.4, 0.412, 0.424, 0.436, 0.448, 0.46, 0.472, 0.485, 0.497, 0.51, 0.523, 0.536, 0.549, 0.563, 0.577, 0.59, 0.604, 0.618, 0.633, 0.648, 0.663, 0.678, 0.693, 0.709, 0.725, 0.741, 0.758, 0.775, 0.793, 0.811, 0.829, 0.848, 0.867, 0.887, 0.908, 0.929, 0.95, 0.973, 0.996, 1.02, 1.045, 1.071, 1.099, 1.127, 1.157, 1.188, 1.221, 1.256, 1.293, 1.333, 1.376, 1.422, 1.472, 1.528, 1.589, 1.658, 1.738, 1.832, 1.946, 2.092, 2.298, 2.647]

    r_idx = vsearch(r, r_values)
    z = z_values[r_idx]

    ci_h = z + z_r * ci2z(ci)
    ci_l = z - z_r * ci2z(ci)

    ci_h_idx = vsearch(ci_h, z_values)
    ci_l_idx = vsearch(ci_l, z_values)

    return (r_values[ci_l_idx], r_values[ci_h_idx])

end

"""
    ci_r(; <keyword arguments>)

Calculate confidence interval for a correlation coefficient.

# Arguments

- `r::Float64`: correlation coefficient
- `n::Int64`: number of observations
- `ci::Float64=0.95`: confidence level

# Returns

- `ci_r::Tuple(Float64, Float64)`
"""
function ci_r(; r::Float64, n::Int64, ci::Float64=0.95)

    z_r = 1 / sqrt(n - 3)

    r_values = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99]
    z_values = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.121, 0.131, 0.141, 0.151, 0.161, 0.172, 0.182, 0.192, 0.203, 0.213, 0.224, 0.234, 0.245, 0.255, 0.266, 0.277, 0.288, 0.299, 0.31, 0.321, 0.332, 0.343, 0.354, 0.365, 0.377, 0.389, 0.4, 0.412, 0.424, 0.436, 0.448, 0.46, 0.472, 0.485, 0.497, 0.51, 0.523, 0.536, 0.549, 0.563, 0.577, 0.59, 0.604, 0.618, 0.633, 0.648, 0.663, 0.678, 0.693, 0.709, 0.725, 0.741, 0.758, 0.775, 0.793, 0.811, 0.829, 0.848, 0.867, 0.887, 0.908, 0.929, 0.95, 0.973, 0.996, 1.02, 1.045, 1.071, 1.099, 1.127, 1.157, 1.188, 1.221, 1.256, 1.293, 1.333, 1.376, 1.422, 1.472, 1.528, 1.589, 1.658, 1.738, 1.832, 1.946, 2.092, 2.298, 2.647]

    r_idx = vsearch(r, r_values)
    z_score = z_values[r_idx]

    ci_h = z_score + z_r * ci2z(ci)
    ci_l = z_score - z_r * ci2z(ci)

    ci_h_idx = vsearch(ci_h, z_values)
    ci_l_idx = vsearch(ci_l, z_values)

    return (r_values[ci_l_idx], r_values[ci_h_idx])

end

"""
    r1r2_test(; <keyword arguments>)

Test if two correlation coefficients are significantly different.

# Arguments

- `r1::Float64`: correlation coefficient, group 1
- `r2::Float64`: correlation coefficient, group 2
- `n1::Int64`: number of observations, group 1
- `n2::Int64`: number of observations, group 2

# Returns

- `z::Float64`: z score
"""
function r1r2_test(; r1::Float64, r2::Float64, n1::Int64, n2::Int64)

    r_values = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99]
    z_values = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.121, 0.131, 0.141, 0.151, 0.161, 0.172, 0.182, 0.192, 0.203, 0.213, 0.224, 0.234, 0.245, 0.255, 0.266, 0.277, 0.288, 0.299, 0.31, 0.321, 0.332, 0.343, 0.354, 0.365, 0.377, 0.389, 0.4, 0.412, 0.424, 0.436, 0.448, 0.46, 0.472, 0.485, 0.497, 0.51, 0.523, 0.536, 0.549, 0.563, 0.577, 0.59, 0.604, 0.618, 0.633, 0.648, 0.663, 0.678, 0.693, 0.709, 0.725, 0.741, 0.758, 0.775, 0.793, 0.811, 0.829, 0.848, 0.867, 0.887, 0.908, 0.929, 0.95, 0.973, 0.996, 1.02, 1.045, 1.071, 1.099, 1.127, 1.157, 1.188, 1.221, 1.256, 1.293, 1.333, 1.376, 1.422, 1.472, 1.528, 1.589, 1.658, 1.738, 1.832, 1.946, 2.092, 2.298, 2.647]

    r1_idx = vsearch(r1, r_values)
    z1_score = z_values[r1_idx]

    r2_idx = vsearch(r2, r_values)
    z2_score = z_values[r2_idx]

    z = (z1_score - z2_score) / sqrt((1 / (n1 - 3)) + (1 / (n2 - 3)))

    return z

end

"""
    slope(p1, p2)

Calculate slope of the line crossing two points.

# Arguments

- `p1::Tuple{Real, Real}`
- `p2::Tuple{Real, Real}`

# Returns

- `s::Float64`: slope
"""
function slope(p1::Tuple{Real, Real}, p2::Tuple{Real, Real})

    s = (p2[2] - p1[2]) / (p2[1] - p1[1])

    return s

end

"""
    slope(p1, p2)

Calculate distance between two points.

# Arguments

- `p1::Tuple{Real, Real}`
- `p2::Tuple{Real, Real}`

# Returns

- `d::Float64`: distance
"""
function distance(p1::Tuple{Real, Real}, p2::Tuple{Real, Real})

    d = sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)

    return d

end

"""
    count_thresh(x; <keyword arguments>)

Collect thresholded elements, e.g. in a topographical map.

# Arguments

- `x::AbstractMatrix`
- `t::Real`: threshold value
- `t_type::Symbol=:g`: rule for thresholding:
    - `:eq`: =
    - `:geq`: ≥
    - `:leq`: ≤
    - `:g`: >
    - `:l`: <

# Returns

Named tuple containing:
- `x_t::Int64`: thresholded matrix
- `n::Int64`: number of elements
"""
function count_thresh(x::AbstractMatrix; t::Real, t_type::Symbol=:g)

    _check_var(t_type, [:eq, :geq, :leq, :g, :l], "t_type")

    x_t = zeros(Bool, size(x))

    if t_type === :eq
        x_t[x .== t] .= true
    elseif t_type === :g
        x_t[x .> t] .= true
    elseif t_type === :geq
        x_t[x .>= t] .= true
    elseif t_type === :l
        x_t[x .< t] .= true
    elseif t_type === :leq
        x_t[x .<= t] .= true
    end

    n = count(==(true), x_t)

    return (x_t=x_t, n=n)

end

"""
    crit_t(df, alpha; <keyword arguments>)

Calculate critical t value.

# Arguments

- `df::Real`: degrees of freedom (usually df = n - 1)
- `alpha::Float64=0.05`: alpha value
- `twosided::Bool=false`: one or two tailed probability

# Returns

- `t::Float64`

# Notes

Critical region for one tailed probability:
- left: `(-∞ , -t]`
- right: `[t , ∞)`

Critical region for two tailed probability: `(-∞ , -t] ∪ [t, ∞)`
"""
function crit_t(df::Real, alpha::Float64=0.05; twosided::Bool=false)

    if twosided
        t = quantile(TDist(df), 1 - (alpha / 2))
    else
        t = quantile(TDist(df), 1 - alpha)
    end

    return t

end

"""
    ci2z(ci)

Calculate critical z score.

# Arguments

- `ci::Float64`: confidence level

# Returns

- `z::Float64`
"""
function ci2z(ci::Float64)

    z = quantile(Distributions.Normal(0, 1), ci)

    return z

end

"""
    p2z(p; <keyword arguments>)

Calculate z score for p value.

# Arguments

- `p::Float64=0.05`: confidence level
- `twosided::Bool=false`: one or two tailed probability

# Returns

- `z::Float64`
"""
function p2z(p::Float64=0.05; twosided::Bool=false)

    if twosided
        z = quantile(Distributions.Normal(0.0, 1.0), 1 - p / 2)
    else
        z = quantile(Distributions.Normal(0.0, 1.0), 1 - p)
    end

    return z

end

"""
    z2p(z; <keyword arguments>)

Calculate probability for a given z value.

# Arguments

- `z::Real`: z value
- `twosided::Bool=false`: one or two tailed probability

# Returns

- `p::Float64`
"""
function z2p(z::Real; twosided::Bool=false)

    if twosided
        p = 2 * ccdf(Distributions.Normal(0.0, 1.0), z)
    else
        p = ccdf(Distributions.Normal(0.0, 1.0), z)
    end

    return p

end

"""
    cmp_stat(stat_dist, v)

Calculate proportion of elements below or above a given statistic value.

# Arguments

- `stat_dist::AbstractVector`: statistic values distribution
- `v::Real`: statistic value
- `type::Symbol=:g`: calculation proportion of elements greater (`:g`) or lesser (`:l`) than `v`

# Returns

- `p::Float64`
"""
function cmp_stat(stat_dist::AbstractVector, v::Real; type::Symbol=:g)

    _check_var(type, [:g, :l], "type")

    type === :g && return count(stat_dist .> v) / length(stat_dist)
    type === :l && return count(stat_dist .< v) / length(stat_dist)

end

"""
    fwhm(s)

Calculate indices of full-width half-maximum points of a Gaussian-like distribution.

# Arguments

- `s::AbstractVector`

# Returns

- `p1_idx::Int64`: pre-peak half-maximum point
- `p_idx::Int64`: peak
- `p2_idx::Int64`: post-peak half-maximum point
"""
function fwhm(s::AbstractVector)

    s = normalize_n(s)

    # peak
    p_idx = vsearch(maximum(s), s)
    # pre-peak half-maximum width
    p1_idx = vsearch(0.5, s[1:p_idx])
    # post-peak half-maximum width
    p2_idx = p_idx + vsearch(0.5, s[p_idx:end]) - 1

    return p1_idx, p_idx, p2_idx

end

"""
    cosine_similarity(s1, s2)

Calculate cosine similarity.

# Arguments

- `s1::AbstractVector`
- `s2::AbstractVector`

# Returns

- `cs::Float64`
"""
function cosine_similarity(s1::AbstractVector, s2::AbstractVector)

    @assert length(s1) == length(s2) "Both vectors must have the same length."

    cs = sum(s1 .* s2) / (sqrt(sum(s1.^2)) * sqrt(sum(s2.^2)))

    return cs

end

"""
    std(obj)

Calculate standard deviation of the signal data (along epochs).

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

- `s::Matrix{Float64}`
"""
function Statistics.std(obj::NeuroAnalyzer.NEURO)

    @assert nepochs(obj) > 1 "OBJ must have > 1 epoch."

    if datatype(obj) == "erp"
        s = @views std(obj.data[:, :, 2:end], dims=3)
    else
        s = @views std(obj.data[:, :, :], dims=3)
    end
    s = reshape(s, size(s, 1), size(s, 2))

    return s

end

"""
    permute(s, n)

Permute signal data.

# Arguments

- `s::AbstractVector`
- `n::Int64`: number of permutations

# Returns

- `s_new::Matrix{Float64}`
"""
function permute(s::AbstractVector, n::Int64)

    @assert n > 0 "n must have > 0 epoch."

    s_new = zeros(n, length(s))
    for idx in 1:n
        x = rand(2:length(s))
        s1 = s[x:end]
        s2 = s[1:(x - 1)]
        s_new[idx, :] = @views vcat(s1, s2)
    end

    return s_new

end

"""
    permute(s, n)

Permute signal data.

# Arguments

- `s::AbstractArray`
- `n::Int64`: number of permutations

# Returns

- `s_new::Matrix{Float64}`
"""
function permute(s::AbstractArray, n::Int64)

    @assert n > 0 "n must have > 0 epoch."
    @assert ndims(s) <= 3 "permute() only works for arrays of ≤ 3 dimensions."

    if ndims(s) == 2
        s_new = zeros(n, size(s,1 ), size(s,2 ))
        @inbounds for idx1 in 1:n
            Threads.@threads for idx2 in axes(s, 1)
                x = rand(2:size(s, 2))
                s1 = s[idx2, x:end]
                s2 = s[idx2, 1:(x - 1)]
                s_new[idx1, idx2, :] = @views vcat(s1, s2)
            end
        end
    else
        s_new = zeros(n, size(s, 1), size(s, 2), size(s, 3))
        for idx1 in 1:n
            @inbounds for idx2 in axes(s, 1)
                Threads.@threads for idx3 in axes(s, 3)
                    x = rand(2:size(s, 2))
                    s1 = s[idx2, x:end, idx3]
                    s2 = s[idx2, 1:(x - 1), idx3]
                    s_new[idx1, idx2, :, idx3] = @views vcat(s1, s2)
                end
            end
        end
    end

    return s_new

end