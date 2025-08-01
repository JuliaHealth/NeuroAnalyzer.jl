export efs
export efs_p1g
export efs_p2g
export stdp

"""
    efs(x1, x2)

Calculate effect sizes.

# Arguments

- `x1::AbstractVector`
- `x2::AbstractVector`

# Returns

Named tuple containing:
- `d::Float64`: Cohen's d
- `g::Float64`: Hedges g, uses maximum likelihood estimator by Hedges and Olkin
- `Δ::Float64`: Glass' Δ
"""
function efs(x1::AbstractVector, x2::AbstractVector)::@NamedTuple{d::Float64, g::Float64, Δ::Float64}

    @assert std(x1) != 0 "std(x1) must not be equal 0."
    @assert std(x2) != 0 "std(x2) must not be equal 0."

    d = (mean(x1) - mean(x2)) / stdp(x1, x2, type=:cohen)
    g = (mean(x1) - mean(x2)) / stdp(x1, x2, type=:hedges)
    Δ = (mean(x1) - mean(x2)) / std(x2)

    return (d=d, g=g, Δ=Δ)

end

"""
    efs_p1g(p)

Calculate Cohen's h effect size for one proportion.

# Arguments

- `p::Float64`: proportion

# Returns

- `h::Float64`
"""
function efs_p1g(p::Float64)::Float64

    _in(p, (0.0, 1.0), "p")

    h = 2 * asin(sqrt(p))

    return h

end

"""
    efs_p2g(p1, p2; nd)

Calculate Cohen's h effect size for two proportions `p1` and `p2`.

# Arguments

- `p1::Float64`: 1st proportion, e.g. 0.7
- `p2::Float64`: 2nd proportion, e.g. 0.3
- `nd::Bool=false`: if true, calculate non-directional h value

# Returns

- `h::Float64`
"""
function efs_p2g(p1::Float64, p2::Float64; nd::Bool=false)::Float64

    _in(p1, (0.0, 1.0), "p1")
    _in(p2, (0.0, 1.0), "p2")
    @assert p1 + p2 == 1.0 "p1 and p1 must add to 1.0."

    h = 2 * asin(sqrt(p1)) - 2 * asin(sqrt(p2))

    return nd ? abs(h) : h

end

"""
    stdp(x1, x2; <keyword arguments>)

Calculate pooled standard deviation.

# Arguments

- `x1::AbstractVector`
- `x2::AbstractVector`
- `type::Symbol=:cohen`: use Cohen's equation (`:cohen`) or maximum likelihood estimator by Hedges and Olkin (`:hedges`)

# Returns

- `sp::Float64`
"""
function stdp(x1::AbstractVector, x2::AbstractVector; type::Symbol=:cohen)::Float64

    @assert length(x1) > 0 "Length of x1 cannot be 0."
    @assert length(x2) > 0 "Length of x2 cannot be 0."
    _check_var(type, [:cohen, :hedges], "type")
    type === :cohen && (@assert length(x1) + length(x2) > 2 "For :cohen equation length of x1 + length of x2 must be > 2.")

    s1 = std(x1)
    s2 = std(x2)
    n1 = length(x1)
    n2 = length(x2)

    if type === :cohen
        sp = sqrt((((n1 - 1) * s1^2) + ((n2 - 1) * s2^2))/(n1 + n2 - 2))
    else
        sp = sqrt((((n1 - 1) * s1^2) + ((n2 - 1) * s2^2))/(n1 + n2))
    end

    return sp

end

"""
    stdp(s1, s2, n1, n2; <keyword arguments>)

Calculate pooled standard deviation when number of subjects in groups are different.

# Arguments

- `s1::Real`
- `s2::Real`
- `n1::Int64`
- `n2::Int64`
- `type::Symbol=:cohen`: use Cohen's equation (`:cohen`) or maximum likelihood estimator by Hedges and Olkin (`:hedges`)

# Returns

- `ps::Float64`
"""
function stdp(s1::Real, s2::Real, n1::Int64, n2::Int64; type::Symbol=:cohen)::Float64

    @assert n1 > 0 "n1 must be > 0."
    @assert n2 > 0 "n2 must be > 0."
    _check_var(type, [:cohen, :hedges], "type")
    type === :cohen && (@assert n1 + n2 > 2 "For :cohen equation n1 + n2 must be > 2.")

    if type === :cohen
        ps = sqrt((((n1 - 1) * s1^2) + ((n2 - 1) * s2^2))/(n1 + n2 - 2))
    else
        ps = sqrt((((n1 - 1) * s1^2) + ((n2 - 1) * s2^2))/(n1 + n2))
    end

    return ps

end


"""
    stdp(s1, s2)

Calculate pooled standard deviation when number of subjects in groups are equal.

# Arguments

- `s1::Real`
- `s2::Real`

# Returns

- `ps::Float64`
"""
function stdp(s1::Real, s2::Real)::Float64

    ps = sqrt((s1^2 + s2^2)/2)

    return ps

end
