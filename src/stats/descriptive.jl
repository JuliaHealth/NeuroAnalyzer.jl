export varp
export varc
export stdp
export stdc
export rng
export mrng
export moe
export arf

"""
    varp(p, n)

Calculate variance of the proportion.

# Arguments

- `p::Float64`: proportion
- `n::Int64`: number of observations

# Returns

- `σ2::Float64`
"""
function varp(p::Float64, n::Int64)::Float64

    @assert n > 0 "n must be > 0."
    _in(p, (0.0, 1.0), "p")

    σ2 = (p * (1 - p)) / n

    return σ2

end

"""
    stdp(p, n)

Calculate standard deviation of the proportion.

# Arguments

- `p::Float64`: proportion
- `n::Int64`: number of observations

# Returns

- `σ::Float64`
"""
function stdp(p::Float64, n::Int64)::Float64

    σ = sqrt(varp(p, n))

    return σ

end

"""
    varc(g, x)

Calculate variance of categorical data.

# Arguments

- `g::Vector{Int64}`: group (e.g. [0, 1, 2, 3])
- `x::Vector{Int64}`: amount of subject per group (e.g. [2, 8, 27, 45])

# Returns

- `σ2::Float64`
"""
function varc(g::Vector{Int64}, x::Vector{Int64})::Float64

    @assert length(g) == length(x) "Length of g and length of x must be equal."
    @assert length(g) > 0 "Length of g must be > 0."
    @assert length(x) > 0 "Length of x must be > 0."

    σ2 = ((sum(g.^2 .* x) - sum(g .* x))^2 / sum(x)) / (sum(x) - 1)

    return σ2

end

"""
    stdc(g, x)

Calculate standard deviation of categorical data.

# Arguments

- `g::Vector{Int64}`: group (e.g. [0, 1, 2, 3])
- `x::Vector{Int64}`: amount of subject per group (e.g. [2, 8, 27, 45])

# Returns

- `σ::Float64`
"""
function stdc(g::Vector{Int64}, x::Vector{Int64})::Float64

    σ = sqrt(varc(g, x))

    return σ

end

"""
    rng(x)

Calculate range.

# Arguments

- `x::AbstractArray`

# Returns

- `r::Float64`
"""
function rng(x::AbstractArray)::Float64

    r = maximum(x) - minimum(x)

    return r

end

"""
    mrng(x)

Calculate midrange.

# Arguments

- `x::AbstractArray`

# Returns

- `mr::Float64`
"""
function mrng(x::AbstractArray)::Float64

    mr = (maximum(x) - minimum(x)) / 2

    return mr

end

"""
    moe(n)

Calculate margin of error for given sample size `n`.

# Arguments

- `n::Int64`

# Returns

- `m::Float64`
"""
function moe(n::Int64)::Float64

    @assert n > 0 "n must be > 0."

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
function moe(x::AbstractArray)::Float64

    @assert length(x) > 0 "Length of x must be > 0."

    n = length(x)
    m = 1 / sqrt(n)

    return m

end

"""
    arf(df, var)

Calculate absolute and relative frequencies.

# Arguments

- `df::DataFrame`
- `var::Union{Symbol, String}`: variable name

# Returns

- `m::Matrix{Float64}`: first row: absolute frequencies, second row: relative frequencies (proportion), third row: second row: relative frequencies (percentage); last column: totals
"""
function arf(df::DataFrame, var::Union{Symbol, String})::Matrix{Float64}

    _check_var(string(var), names(df), "var")
    x = df[!, var]
    @assert length(unique(x)) >= 2 "var must contain at least 2 different values."
    n = length(x)
    m = zeros(3, length(unique(x)) + 1)
    for idx in 1:length(unique(x))
        m[1, idx] = count(z -> z==unique(x)[idx], x)
        m[2, idx] = round(m[1, idx] / n, digits=3)
        m[3, idx] = round(m[2, idx] * 100, digits=2)
    end
    m[1, end] = n
    m[2, end] = 1.0
    m[3, end] = 100.0

    return m

end
