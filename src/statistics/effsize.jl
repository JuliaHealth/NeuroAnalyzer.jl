export effsize
export effsize_p2g
export pooledstd

"""
    effsize(x1, x2)

Calculate Cohen's d and Hedges g effect sizes.

# Arguments

- `x1::AbstractVector`
- `x2::AbstractVector`

# Returns

Named tuple containing:
- `d::Float64`: Cohen's d
- `g::Float64`: Hedges g, uses maximum likelihood estimator by Hedges and Olkin
- `Δ::Float64`: Glass' Δ
"""
function effsize(x1::AbstractVector, x2::AbstractVector)::@NamedTuple{d::Float64, g::Float64, Δ::Float64}

    d = (mean(x2) - mean(x1)) / pooledstd(x1, x2, type=:cohen)
    g = (mean(x2) - mean(x1)) / pooledstd(x1, x2, type=:hedges)
    Δ = (mean(x2) - mean(x1)) / std(x2)

    return (d=d, g=g, Δ=Δ)

end

"""
    effsize_p2g(p1, p2)

Calculate effect size for two proportions `p1` and `p2`.

# Arguments

- `p1::Float64`: 1st proportion, e.g. 0.7
- `p2::Float64`: 2nd proportion, e.g. 0.3

# Returns

- `e::Float64`
"""
function effsize_p2g(p1::Float64, p2::Float64)::Float64

    @assert p1 + p2 == 1.0 "Proportions must add to 1.0."

    e = 2 * asin(sqrt(p1)) - 2 * asin(sqrt(p2))

    return e

end

"""
    pooledstd(x1, x2; <keyword arguments>)

Calculate pooled standard deviation

# Arguments

- `x1::AbstractVector`
- `x2::AbstractVector`
- `type::Symbol=:cohen`: use Cohen's equation (`:cohen`) or maximum likelihood estimator by Hedges and Olkin (`:hedges`)

# Returns

- `ps::Float64`
"""
function pooledstd(x1::AbstractVector, x2::AbstractVector; type::Symbol=:cohen)::Float64

    _check_var(type, [:cohen, :hedges], "type")

    v1 = var(x1)
    v2 = var(x2)
    n1 = length(x1)
    n2 = length(x2)

    if type === :cohen
        ps = sqrt((((n1 - 1) * v1) + ((n2 - 1) * v2))/(n1 + n2 - 2))
    else
        ps = sqrt((((n1 - 1) * v1) + ((n2 - 1) * v2))/(n1 + n2))
    end

    return ps

end
