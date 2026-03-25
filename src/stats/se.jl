export sem
export semd
export sep
export sen
export sem_diff
export sep_diff
export sen_diff
export ses
export sek

"""
    sem(x)

Calculate the standard error of the mean.

Computed as `std(x) / √n`.

# Arguments

- `x::AbstractVector`: input vector; must contain at least 2 elements

# Returns

- `Float64`: standard error of the mean
"""
function sem(x::AbstractVector)::Float64

    # validate
    length(x) >= 2 || throw(ArgumentError("x must contain at least 2 elements."))

    return std(x) / sqrt(length(x))

end

"""
    semd(x)

Calculate the standard error of the median.

Approximated as `1.253 × std(x) / √n` (valid for large normal samples).

# Arguments

- `x::AbstractVector`: input vector; must contain at least 2 elements

# Returns

- `Float64`: standard error of the median
"""
function semd(x::AbstractVector)::Float64

    # validate
    length(x) >= 2 || throw(ArgumentError("x must contain at least 2 elements."))

    return 1.253 * std(x) / sqrt(length(x))

end

"""
    sep(p, n)

Calculate the standard error of a proportion.

Computed as `√(p(1 − p) / n)`.

# Arguments
- `p::Float64`: proportion; must be in `[0, 1]`
- `n::Int64`: number of observations; must be ≥ 1

# Returns

- `Float64`: standard error of the proportion
"""
function sep(p::Float64, n::Int64)::Float64

    # validate
    _in(p, (0.0, 1.0), "p")
    n >= 1 || throw(ArgumentError("n must be ≥ 1."))

    return sqrt((p * (1 - p)) / n)

end

"""
    sen(n)

Calculate the standard error of a count (`√n`).

# Arguments

- `n::Int64`: number of observations; must be ≥ 1

# Returns

- `Float64`: √n.
"""
function sen(n::Int64)::Float64

    # validate
    n >= 1 || throw(ArgumentError("n must be ≥ 1."))

    return sqrt(n)

end

"""
    sem_diff(x, y)

Calculate the standard error of the difference between two means.

For equal-length vectors: `√(SEM(x)² + SEM(y)²)`.

For unequal-length vectors: pooled SD × `√(1/n1 + 1/n2)`.

# Arguments

- `x::AbstractVector`: first sample; must contain at least 2 elements
- `y::AbstractVector`: second sample; must contain at least 2 elements

# Returns

- `Float64`: standard error of the mean difference
"""
function sem_diff(x::AbstractVector, y::AbstractVector)::Float64

    # validate
    length(x) >= 2 || throw(ArgumentError("x must contain at least 2 elements."))
    length(y) >= 2 || throw(ArgumentError("y must contain at least 2 elements."))

    if length(x) == length(y)
        return sqrt(sem(x)^2 + sem(y)^2)
    else
        return stdp(x, y) * sqrt(1 / length(x) + 1 / length(y))
    end

end

"""
    sep_diff(p1, p2, n1, n2)

Calculate the standard error of the difference between two proportions.

Computed as `√(p1(1−p1)/n1 + p2(1−p2)/n2)`.

# Arguments
- `p1::Float64`: proportion of group 1; must be in `[0, 1]`
- `p2::Float64`: proportion of group 2; must be in `[0, 1]`
- `n1::Int64`: group 1 sample size; must be ≥ 1
- `n2::Int64`: group 2 sample size; must be ≥ 1

# Returns

- `Float64`: standard error of the difference in proportions
"""
function sep_diff(
    p1::Float64,
    p2::Float64,
    n1::Int64,
    n2::Int64
)::Float64

    # validate
    _in(p1, (0.0, 1.0), "p1")
    _in(p2, (0.0, 1.0), "p2")
    n1 >= 1 || throw(ArgumentError("n1 must be ≥ 1."))
    n2 >= 1 || throw(ArgumentError("n2 must be ≥ 1."))

    return sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)

end

"""
    sen_diff(n1, n2)

Calculate the standard error of the difference between two counts.

Computed as `√(n1 + n2)`.

# Arguments

- `n1::Int64`: group 1 count; must be ≥ 1
- `n2::Int64`: group 2 count; must be ≥ 1

# Returns

- `Float64`: `√(n1 + n2)`
"""
function sen_diff(n1::Int64, n2::Int64)::Float64

    # validate
    n1 >= 1 || throw(ArgumentError("n1 must be ≥ 1."))
    n2 >= 1 || throw(ArgumentError("n2 must be ≥ 1."))

    return sqrt(n1 + n2)

end

"""
    ses(x)

Calculate the standard error of skewness.

Computed as `√(6n(n−1) / ((n−2)(n+1)(n+3)))`.

Requires `n ≥ 3` so that the denominator is non-zero.

# Arguments

- `x::AbstractVector`: input vector; must contain at least 3 elements

# Returns

- `Float64`: standard error of skewness
"""
function ses(x::AbstractVector)::Float64

    # validate
    n = length(x)
    n >= 3 || throw(ArgumentError("x must contain at least 3 elements."))

    return sqrt((6 * n * (n - 1)) / ((n - 2) * (n + 1) * (n + 3)))

end

"""
    ses(n)

Calculate the standard error of skewness for a sample of size `n`.

# Arguments

- `n::Int64`: sample size; must be ≥ 3

# Returns

- `Float64`: standard error of skewness
"""
function ses(n::Int64)::Float64

    # validate
    n >= 3 || throw(ArgumentError("n must be ≥ 3."))

    return sqrt((6 * n * (n - 1)) / ((n - 2) * (n + 1) * (n + 3)))

end

"""
    sek(x)

Calculate the standard error of kurtosis.

Computed as `2 × (n−1) × √(6n / ((n−2)(n−3)(n+3)(n+5)))`.

Requires `n ≥ 4` so that the `(n−3)` term in the denominator is non-zero.

# Arguments

- `x::AbstractVector`: input vector; must contain at least 4 elements.

# Returns

- `Float64`: standard error of kurtosis
"""
function sek(x::AbstractVector)::Float64

    # validate
    n = length(x)
    n >= 4 || throw(ArgumentError("x must contain at least 4 elements."))

    return 2 * (n - 1) * sqrt((6 * n) / ((n - 2) * (n - 3) * (n + 3) * (n + 5)))

end

"""
    sek(n)

Calculate the standard error of kurtosis for a sample of size `n`.

Computed as `2 × (n−1) × √(6n / ((n−2)(n−3)(n+3)(n+5)))`.

# Arguments

- `n::Int64`: sample size; must be ≥ 4

# Returns

- `Float64`: standard error of kurtosis
"""
function sek(n::Int64)::Float64

    # validate
    n >= 4 || throw(ArgumentError("n must be ≥ 4."))

    return 2 * (n - 1) * sqrt((6 * n) / ((n - 2) * (n - 3) * (n + 3) * (n + 5)))

end
