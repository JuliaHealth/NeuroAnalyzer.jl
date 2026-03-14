export k_categories
export slope
export distance
export count_thresh
export cmp_stat
export permute
export logit
export sumsq
export rmna
export df
export center

"""
    k_categories(n)

Calculate the recommended number of histogram categories for a sample of size `n`.

Returns two common rules:

- Square-root choice: `k1 = √n`
- Sturges' extended rule: `k2 = 1 + 3.222 × log₁₀(n)`

# Arguments

- `n::Int64`: sample size; must be ≥ 1

# Returns

Named tuple:

- `k1::Float64`: square-root estimate
- `k2::Float64`: sturges' extended estimate

# Throws

- `ArgumentError`: if `n < 1`
"""
function k_categories(n::Int64)::@NamedTuple{k1::Float64, k2::Float64}

    !(n >= 1) && throw(ArgumentError("n must be ≥ 1."))

    return (k1=sqrt(n), k2=1 + 3.222 * log10(n))

end

"""
    slope(p1, p2)

Calculate the slope of the line passing through two points.

# Arguments

- `p1::Tuple{Real, Real}`: first point `(x₁, y₁)`
- `p2::Tuple{Real, Real}`: second point `(x₂, y₂)`; `x₂ ≠ x₁`

# Returns

- `Float64`: slope `(y₂ − y₁) / (x₂ − x₁)`

# Throws

- `ArgumentError`: if `p2[1] == p1[1]` (vertical line; slope undefined)

# See also

[`distance`](@ref)
"""
function slope(p1::Tuple{Real, Real}, p2::Tuple{Real, Real})::Float64

    !(p2[1] != p1[1]) && throw(ArgumentError("p2[1] and p1[1] must not be equal (vertical line has no slope)."))

    return (p2[2] - p1[2]) / (p2[1] - p1[1])

end

"""
    distance(p1, p2)

Calculate the Euclidean distance between two points.

# Arguments
- `p1::Tuple{Real, Real}`: first point `(x₁, y₁)`
- `p2::Tuple{Real, Real}`: second point `(x₂, y₂)`

# Returns

- `Float64`: Euclidean distance `√((x₂−x₁)² + (y₂−y₁)²)`

# See also

[`slope`](@ref)
"""
function distance(p1::Tuple{Real, Real}, p2::Tuple{Real, Real})::Float64

    return sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)

end

"""
    count_thresh(x; <keyword arguments>)

Threshold a matrix and count the number of elements satisfying the condition.

# Arguments

- `x::AbstractMatrix`: input matrix
- `t::Real`: threshold value
- `t_type::Symbol=:g`: thresholding rule:
  - `:eq`:  `x == t`
  - `:geq`: `x ≥ t`
  - `:leq`: `x ≤ t`
  - `:g`:   `x > t`
  - `:l`:   `x < t`

# Returns

Named tuple:

- `x_t::Matrix{Bool}`: boolean mask; `true` where the condition holds
- `n::Int64`: number of elements satisfying the condition

# Throws

- `ArgumentError`: if `t_type` is not a recognized symbol

# See also

[`cmp_stat`](@ref)
"""
function count_thresh(x::AbstractMatrix; t::Real, t_type::Symbol = :g)::@NamedTuple{x_t::Matrix{Bool}, n::Int64}

    _check_var(t_type, [:eq, :geq, :leq, :g, :l], "t_type")

    x_t = if t_type === :eq
        x .== t
    elseif t_type === :g
        x .> t
    elseif t_type === :geq
        x .>= t
    elseif t_type === :l
        x .< t
    elseif t_type === :leq
        x .<= t
    end

    return (x_t=Matrix{Bool}(x_t), n=count(x_t))

end

"""
    cmp_stat(stat_dist, v)

Calculate the proportion of elements in a statistic distribution that are greater or lesser than a given value.

# Arguments

- `stat_dist::AbstractVector`: distribution of statistic values; must not be empty
- `v::Real`: reference statistic value
- `type::Symbol=:g`: Comparison direction:
    - `:g`: proportion of elements > `v`
    - `:l`: proportion of elements < `v`

# Returns

- `Float64`: proportion of elements satisfying the condition

# Throws

- `ArgumentError`: if `stat_dist` is empty or `type ∉ {:g, :l}`

# See also

[`count_thresh`](@ref)
"""
function cmp_stat(stat_dist::AbstractVector, v::Real; type::Symbol = :g)::Float64

    _check_var(type, [:g, :l], "type")
    !(length(stat_dist) > 0) && throw(ArgumentError("stat_dist must not be empty."))

    if type === :g
        return count(>(v), stat_dist) / length(stat_dist)
    elseif type === :l
        return count(<(v), stat_dist) / length(stat_dist)
    end

end

"""
    permute(s, n)

Generate `n` cyclic permutations of a signal vector.

Each permutation randomly selects a split point and rotates the vector (moves the tail before the head).

# Arguments

- `s::AbstractVector`: input signal vector; must contain at least 2 elements
- `n::Int64`: number of permutations; must be ≥ 1

# Returns

- `Matrix{Float64}`: matrix of shape `(n × length(s))`

# Throws

- `ArgumentError`: if `n < 1` or `length(s) < 2`

# See also

[`permute(::AbstractArray, ::Int64)`](@ref)
"""
function permute(s::AbstractVector, n::Int64)::Matrix{Float64}

    !(n >= 1) && throw(ArgumentError("n must be ≥ 1."))
    !(length(s) >= 2) && throw(ArgumentError("s must contain at least 2 elements."))

    # pre-allocate output
    s_new = zeros(n, length(s))
    for idx in 1:n
        x = rand(2:length(s))
        s_new[idx, :] = vcat(s[x:end], s[1:(x - 1)])
    end

    return s_new

end

"""
    permute(s, n)

Generate `n` cyclic permutations of each row of a 2- or 3-D array.

Each permutation randomly selects a split point along the second axis and rotates that row cyclically.

# Arguments

- `s::AbstractArray`: input array; must be 2- or 3-dimensional, with at least 2 elements along the second axis
- `n::Int64`: number of permutations; must be ≥ 1

# Returns

- `Array{Float64, 3}`: shape `(n × size(s,1) × size(s,2))` for 2-D input
- `Array{Float64, 4}`: shape `(n × size(s,1) × size(s,2) × size(s,3))` for 3-D input

# Throws

- `ArgumentError`: if `n < 1`, `ndims(s) ∉ {2, 3}`, or `size(s, 2) < 2`

# See also

[`permute(::AbstractVector, ::Int64)`](@ref)
"""
function permute(s::AbstractArray, n::Int64)::Union{Array{Float64, 3}, Array{Float64, 4}}

    !(n >= 1) && throw(ArgumentError("n must be ≥ 1."))
    !(2 <= ndims(s) <= 3) && throw(ArgumentError("permute() only supports 2- and 3-dimensional arrays."))
    !(size(s, 2) >= 2) && throw(ArgumentError("Second dimension of s must be ≥ 2."))

    ncols = size(s, 2)

    if ndims(s) == 2
        nrows = size(s, 1)
        s_new = zeros(n, nrows, ncols)

        # Thread over permutations (outer loop): each idx1 writes to a unique
        # s_new[idx1, :, :] slice → zero contention between threads
        Threads.@threads :static for idx1 in 1:n
            @inbounds for idx2 in 1:nrows
                x = rand(2:ncols)
                s_new[idx1, idx2, 1:(ncols - x + 1)]   .= @view s[idx2, x:end]
                s_new[idx1, idx2, (ncols - x + 2):end] .= @view s[idx2, 1:(x - 1)]
            end
        end

    elseif ndims(s) == 3
        nrows  = size(s, 1)
        nepochs = size(s, 3)
        s_new  = zeros(n, nrows, ncols, nepochs)

        # Thread over permutations: each idx1 writes to s_new[idx1, :, :, :]
        # the two inner loops are sequential (short; threading would add overhead)
        Threads.@threads :static for idx1 in 1:n
            @inbounds for idx3 in 1:nepochs
                for idx2 in 1:nrows
                    x = rand(2:ncols)
                    s_new[idx1, idx2, 1:(ncols - x + 1), idx3]   .= @view s[idx2, x:end, idx3]
                    s_new[idx1, idx2, (ncols - x + 2):end, idx3] .= @view s[idx2, 1:(x - 1), idx3]
                end
            end
        end
    end

    return s_new

end

"""
    logit(p)

Convert a proportion to its logit (log-odds).

Computed as `log(p / (1 − p))`. Returns `−Inf` for `p = 0` and `+Inf` for `p = 1`.

# Arguments

- `p::Float64`: proportion; must be in `[0, 1]`

# Returns

- `Float64`: log-odds `log(p / (1 − p))`

# Throws

- `ArgumentError`: if `p ∉ [0, 1]`
"""
function logit(p::Float64)::Float64

    _in(p, (0.0, 1.0), "p")
    return log(p / (1 - p))

end

"""
    sumsq(x)

Calculate the sum of squared deviations from the mean.

Computed as `Σ(xᵢ − x̄)²`.

# Arguments

- `x::AbstractVector`: input vector; must contain at least 2 elements

# Returns

- `Float64`: sum of squared deviations

# Throws

- `ArgumentError`: if `length(x) < 2`

# See also

[`df`](@ref)
"""
function sumsq(x::AbstractVector)::Float64

    !(length(x) >= 2) && throw(ArgumentError("x must contain at least 2 elements."))
    m = mean(x)

    return sum((x .- m) .^ 2)

end

"""
    rmna(x)

Return a copy of `x` with all `NaN` and `Missing` values removed.

# Arguments

- `x::AbstractVector`: input vector; may contain `Float64`, `Missing`, or `NaN`

# Returns
- `Vector{Float64}`: filtered vector converted to `Float64`
"""
function rmna(x::AbstractVector)::Vector{Float64}

    # filter missing first (type-level), then NaN (value-level)
    x_clean = collect(skipmissing(x))

    return Float64.(Base.filter(!isnan, x_clean))

end


"""
    df(x)

Calculate the degrees of freedom for a vector (`length(x) − 1`).

# Arguments

- `x::AbstractVector`: input vector; must contain at least 1 element

# Returns

- `Int64`: degrees of freedom `length(x) − 1`

# Throws

- `ArgumentError`: if `x` is empty

# See also

[`sumsq`](@ref)
"""
function df(x::AbstractVector)::Int64

    !(length(x) >= 1) && throw(ArgumentError("x must not be empty."))

    return length(x) - 1

end

"""
    center(x)

Center a vector by subtracting its mean.

# Arguments

- `x::AbstractVector`: input vector; must contain at least 1 element

# Returns

- `Vector{Float64}`: mean-centered vector `x .- mean(x)`

# Throws

- `ArgumentError`: if `x` is empty
"""
function center(x::AbstractVector)::Vector{Float64}

    !(length(x) >= 1) && throw(ArgumentError("x must not be empty."))

    return x .- mean(x)

end
