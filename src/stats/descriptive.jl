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

Calculate the variance of a proportion.

Computed as `p × (1 − p) / n` (the variance of a binomial proportion estimator).

# Arguments

- `p::Float64`: proportion; must be in `[0, 1]`
- `n::Int64`: number of observations; must be ≥ 1

# Returns

- `Float64`: variance of the proportion

# Throws

- `ArgumentError`: if `p ∉ [0, 1]` or `n < 1`

# See also

[`stdp`](@ref), [`varp`](@ref)
"""
function varp(p::Float64, n::Int64)::Float64

    _in(p, (0.0, 1.0), "p")
    @assert n >= 1 "n must be ≥ 1."

    return (p * (1 - p)) / n

end

"""
    stdp(p, n)

Calculate the standard deviation of a proportion.

Computed as `√(p × (1 − p) / n)`.

# Arguments

- `p::Float64`: proportion; must be in `[0, 1]`
- `n::Int64`: number of observations; must be ≥ 1

# Returns

- `Float64`: standard deviation of the proportion

# Throws

- `ArgumentError`: if `p ∉ [0, 1]` or `n < 1`

# See also

[`varp`](@ref), [`stdc`](@ref)
"""
function stdp(p::Float64, n::Int64)::Float64

    return sqrt(varp(p, n))

end

"""
    varc(g, x)

Calculate the variance of categorical data using group labels and counts.

# Arguments

- `g::Vector{Int64}`: group labels (e.g. `[0, 1, 2, 3]`)
- `x::Vector{Int64}`: number of subjects per group (e.g. `[2, 8, 27, 45]`); must be the same length as `g`, non-empty, and sum to > 1

# Returns

- `Float64`: variance of the categorical variable

# Throws

- `ArgumentError`: if `length(g) ≠ length(x)`, either is empty, or `sum(x) ≤ 1`

# Notes

Formula: `(Σ(g² × x) − (Σ(g × x))² / Σx) / (Σx − 1)`

# See also

[`stdc`](@ref), [`varp`](@ref)
"""
function varc(g::Vector{Int64}, x::Vector{Int64})::Float64

    @assert length(g) > 0 "g must not be empty."
    @assert length(g) == length(x) "g and x must have the same length."
    n = sum(x)
    @assert n > 1 "sum(x) must be > 1 (at least two observations needed)."

    σ2 = (sum(g .^ 2 .* x) - sum(g .* x)^2 / n) / (n - 1)

    return σ2

end

"""
    stdc(g, x)

Calculate the standard deviation of categorical data.

# Arguments

- `g::Vector{Int64}`: group labels (e.g. `[0, 1, 2, 3]`)
- `x::Vector{Int64}`: number of subjects per group; same length as `g`

# Returns

- `Float64`: standard deviation of the categorical variable

# Throws

- `ArgumentError`: if `length(g) ≠ length(x)`, either is empty, or `sum(x) ≤ 1`

# See also

[`varc`](@ref), [`stdp`](@ref)
"""
function stdc(g::Vector{Int64}, x::Vector{Int64})::Float64

    return sqrt(varc(g, x))

end

"""
    rng(x)

Calculate the range of an array (maximum − minimum).

# Arguments

- `x::AbstractArray`: input array; must not be empty

# Returns

- `Float64`: range of `x`

# Throws

- `ArgumentError`: if `x` is empty

# See also

[`mrng`](@ref)
"""
function rng(x::AbstractArray)::Float64

    @assert length(x) > 0 "x must not be empty."

    return Float64(maximum(x) - minimum(x))

end

"""
    mrng(x)

Calculate the midrange of an array: `(maximum(x) − minimum(x)) / 2`.

# Arguments

- `Float64`: midrange of `x`.

# Throws

- `ArgumentError`: if `x` is empty

# See also

[`rng`](@ref)
"""
function mrng(x::AbstractArray)::Float64

    @assert length(x) > 0 "x must not be empty."

    return Float64((maximum(x) - minimum(x)) / 2)

end

"""
    moe(n)

Calculate the margin of error for a given sample size.

Computed as `1 / √n` (the standard error of a proportion at `p = 0.5`).

# Arguments

- `n::Int64`: sample size; must be ≥ 1

# Returns

- `Float64`: margin of error

# Throws

- `ArgumentError`: if `n < 1`

# See also

[`moe(::AbstractArray)`](@ref), [`varp`](@ref)
"""
function moe(n::Int64)::Float64

    @assert n >= 1 "n must be ≥ 1."

    return 1 / sqrt(n)

end

"""
    moe(x)

Calculate the margin of error for a data array based on its sample size.

Computed as `1 / √length(x)`.

# Arguments

- `x::AbstractArray`: input array; must not be empty

# Returns

- `Float64`: margin of error

# Throws

- `ArgumentError`: if `x` is empty

# See also

[`moe(::Int64)`](@ref), [`varp`](@ref)
"""
function moe(x::AbstractArray)::Float64

    @assert length(x) > 0 "x must not be empty."

    return 1 / sqrt(length(x))

end

"""
    arf(df, var)

Calculate absolute and relative frequencies for a categorical variable.

# Arguments

- `df::DataFrame`: input DataFrame
- `var::Union{Symbol, String}`: column name; must be present in `df` and contain at least 2 distinct values

# Returns

- `Matrix{Float64}`: a `(3 × (k + 1))` matrix where `k = length(unique(x))`:
    - row 1: absolute frequencies per category, total in last column
    - row 2: relative frequencies as proportions (rounded to 3 d.p.), total = 1.0
    - row 3: relative frequencies as percentages (rounded to 2 d.p.), total = 100.0

# Throws

- `ArgumentError`: if `var` is not a column of `df` or contains fewer than 2 distinct values
"""
function arf(df::DataFrame, var::Union{Symbol, String})::Matrix{Float64}

    x = df[!, var]
    uvals = unique(x)
    k = length(uvals)
    @assert k >= 2 "var must contain at least 2 distinct values."

    n = length(x)
    m = zeros(3, k + 1)

    for (idx, val) in enumerate(uvals)
        abs_freq = count(==(val), x)
        m[1, idx] = abs_freq
        m[2, idx] = round(abs_freq / n, digits=3)
        m[3, idx] = round(m[2, idx] * 100, digits=2)
    end

    # totals column
    m[1, end] = n
    m[2, end] = 1.0
    m[3, end] = 100.0

    return m

end
