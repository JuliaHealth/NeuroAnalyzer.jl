export efs
export efs_p1g
export efs_p2g
export stdp

"""
    efs(x1, x2)

Calculate effect size measures for two independent samples.

# Arguments

- `x1::AbstractVector`: first sample; must contain ≥ 2 elements with non-zero SD
- `x2::AbstractVector`: second sample; must contain ≥ 2 elements with non-zero SD

# Returns

Named tuple:

- `d::Float64`: Cohen's d (pooled SD denominator, unbiased for equal n)
- `g::Float64`: Hedges' g (maximum-likelihood pooled SD; less biased for small n)
- `delta::Float64`: Glass' Δ (uses SD of `x2` as the denominator; preferred when groups have different variances)
"""
function efs(
    x1::AbstractVector,
    x2::AbstractVector,
)::@NamedTuple{d::Float64, g::Float64, Δ::Float64}

    # validate
    length(x1) >= 2 || throw(ArgumentError("x1 must contain at least 2 elements."))
    length(x2) >= 2 || throw(ArgumentError("x2 must contain at least 2 elements."))
    std(x1) != 0 || throw(ArgumentError("std(x1) must not be zero."))
    std(x2) != 0 || throw(ArgumentError("std(x2) must not be zero."))

    Δm = mean(x1) - mean(x2)

    d = Δm / stdp(x1, x2; type = :cohen)
    g = Δm / stdp(x1, x2; type = :hedges)
    delta = Δm / std(x2)

    return (; d, g, delta)
end

"""
    efs_p1g(p)

Calculate Cohen's h arcsine transformation for a single proportion.

Computes the arcsine-transformed proportion `φ = 2 × arcsin(√p)`, which stabilizes the variance of a binomial proportion.

# Arguments

- `p::Float64`: proportion; must be in `[0, 1]`

# Returns

- `Float64`: transformed value `φ = 2 × arcsin(√p)`
"""
function efs_p1g(p::Float64)::Float64
    _in(p, (0.0, 1.0), "p")

    return 2 * asin(sqrt(p))
end

"""
    efs_p2g(p1, p2; nd)

Calculate Cohen's h effect size for two independent proportions.

Computed as `h = 2arcsin(√p1) − 2arcsin(√p2)`. The two proportions are **independent** and do not need to sum to 1.

# Arguments

- `p1::Float64`: first proportion; must be in `[0, 1]`
- `p2::Float64`: second proportion; must be in `[0, 1]`
- `nd::Bool=false`: if `true`, return `|h|` (non-directional effect size)

# Returns

- `Float64`: Cohen's h (signed, or absolute if `nd=true`)
"""
function efs_p2g(p1::Float64, p2::Float64; nd::Bool = false)::Float64
    _in(p1, (0.0, 1.0), "p1")
    _in(p2, (0.0, 1.0), "p2")

    h = 2 * asin(sqrt(p1)) - 2 * asin(sqrt(p2))

    return nd ? abs(h) : h
end

"""
    stdp(x1, x2; <keyword arguments>)

Calculate the pooled standard deviation from two sample vectors.

# Arguments
- `x1::AbstractVector`: first sample; must contain ≥ 2 elements
- `x2::AbstractVector`: second sample; must contain ≥ 2 elements
- `type::Symbol=:cohen`: pooling method:
    - `:cohen`: `√(((n1−1)s1² + (n2−1)s2²) / (n1+n2−2))`; requires `n1+n2 > 2`
    - `:hedges`: `√(((n1−1)s1² + (n2−1)s2²) / (n1+n2))` (MLE; slightly biased downward)

# Returns

- `Float64`: pooled standard deviation
"""
function stdp(x1::AbstractVector, x2::AbstractVector; type::Symbol = :cohen)::Float64

    # validate
    length(x1) > 0 || throw(ArgumentError("Length of x1 cannot be 0."))
    length(x1) >= 2 || throw(ArgumentError("x1 must contain at least 2 elements."))
    length(x2) >= 2 || throw(ArgumentError("x2 must contain at least 2 elements."))
    _check_var(type, [:cohen, :hedges], "type")

    n1, n2 = length(x1), length(x2)
    s1, s2 = std(x1), std(x2)

    if type === :cohen
        n1 + n2 > 2 || throw(ArgumentError("For :cohen, n1 + n2 must be > 2."))
        return sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
    elseif type === :hedges
        return sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2))
    end
end

"""
    stdp(s1, s2, n1, n2; <keyword arguments>)

Calculate the pooled standard deviation from summary statistics when group sizes differ.

# Arguments

- `s1::Real`: standard deviation of group 1; must be ≥ 0
- `s2::Real`: standard deviation of group 2; must be ≥ 0
- `n1::Int64`: sample size of group 1; must be ≥ 2
- `n2::Int64`: sample size of group 2; must be ≥ 2
- `type::Symbol=:cohen`: pooling method:
    - `:cohen`: `√(((n1−1)s1² + (n2−1)s2²) / (n1+n2−2))`; requires `n1+n2 > 2`
    - `:hedges`: `√(((n1−1)s1² + (n2−1)s2²) / (n1+n2))` (MLE; slightly biased downward)

# Returns

- `Float64`: pooled standard deviation
"""
function stdp(s1::Real, s2::Real, n1::Int64, n2::Int64; type::Symbol = :cohen)::Float64

    # validate
    s1 >= 0 || throw(ArgumentError("s1 must be ≥ 0."))
    s2 >= 0 || throw(ArgumentError("s2 must be ≥ 0."))
    n1 >= 2 || throw(ArgumentError("n1 must be ≥ 2."))
    n2 >= 2 || throw(ArgumentError("n2 must be ≥ 2."))
    _check_var(type, [:cohen, :hedges], "type")

    if type === :cohen
        n1 + n2 > 2 || throw(ArgumentError("For :cohen, n1 + n2 must be > 2."))
        return sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
    elseif type === :hedges
        return sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2))
    end
end

"""
    stdp(s1, s2)

Calculate the pooled standard deviation when both groups have equal size.

Computed as `√((s1² + s2²) / 2)`.

# Arguments

- `s1::Real`: standard deviation of group 1; must be ≥ 0
- `s2::Real`: standard deviation of group 2; must be ≥ 0

# Returns

- `Float64`: pooled standard deviation
"""
function stdp(s1::Real, s2::Real)::Float64

    # validate
    s1 >= 0 || throw(ArgumentError("s1 must be ≥ 0."))
    s2 >= 0 || throw(ArgumentError("s2 must be ≥ 0."))

    return sqrt((s1^2 + s2^2) / 2)
end
