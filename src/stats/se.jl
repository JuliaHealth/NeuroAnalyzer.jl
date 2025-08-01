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

Calculate standard error of the mean.

# Arguments

- `x::AbstractVector`

# Returns

- `s::Float64`
"""
function sem(x::AbstractVector)::Float64

    s = std(x) / sqrt(length(x))

    return s

end

"""
    semd(x)

Calculate standard error of the median.

# Arguments

- `x::AbstractVector`

# Returns

- `s::Float64`
"""
function semd(x::AbstractVector)::Float64

    s = 1.253 * std(x) / sqrt(length(x))

    return s

end

"""
    sep(p, n)

Calculate standard error of the proportion.

# Arguments

- `p::Float64`: proportion
- `n::Int64`: number of observations

# Returns

- `s::Float64`
"""
function sep(p::Float64, n::Int64)::Float64

    @assert n > 0 "n must be > 0."
    _in(p, (0.0, 1.0), "p")

    s = sqrt((p * (1 - p)) / n)

    return s

end

"""
    sen(n)

Calculate standard error of the number.

# Arguments

- `n::Int64`: number of observations

# Returns

- `s::Float64`
"""
function sen(n::Int64)::Float64

    @assert n > 0 "n must be > 0."

    s = sqrt(n)

    return s

end

"""
    sem_diff(x, y)

Calculate SEM (standard error of the mean) of difference between two means.

# Arguments

- `x::AbstractVector`
- `y::AbstractVector`

# Returns

- `sd::Float64`
"""
function sem_diff(x::AbstractVector, y::AbstractVector)::Float64

    if length(x) == length(y)
        sd = sqrt(sem(x)^2 + sem(y)^2)
    else
        sd = stdp(x, y) * sqrt(1/legth(x) + 1/length(y))
    end

    return sd

end

"""
    sep_diff(p1, p2, n1, n2)

Calculate standard error of the difference of two proportions.

# Arguments

- `p1::Float64`: proportion 1
- `p2::Float64`: proportion 1
- `n1::Int64`: number of observations in group 1
- `n2::Int64`: number of observations in group 2

# Returns

- `s::Float64`
"""
function sep_diff(p1::Float64, p2::Float64, n1::Int64, n2::Int64)::Float64

    @assert n1 > 0 "n1 must be > 0."
    @assert n2 > 0 "n2 must be > 0."
    _in(p1, (0.0, 1.0), "p1")
    _in(p2, (0.0, 1.0), "p2")

    s = sqrt(((p1 * (1 - p1)) / n1) + ((p2 * (1 - p2)) / n2))

    return s

end

"""
    sen_diff(n1, n2)

Calculate standard error of the difference between two numbers.

# Arguments

- `n1::Int64`: number of observations in group 1
- `n2::Int64`: number of observations in group 2

# Returns

- `s::Float64`
"""
function sen_diff(n1::Int64, n2::Int64)::Float64

    @assert n1 > 0 "n1 must be > 0."
    @assert n2 > 0 "n2 must be > 0."

    s = sqrt(n1 + n2)

    return s

end

"""
    ses(x)

Calculate standard error of the skewness.

# Arguments

- `x::AbstractVector`

# Returns

- `s::Float64`
"""
function ses(x::AbstractVector)::Float64

    n = length(x)
    @assert n > 0 "x length must be > 0."

    s = sqrt((6 * n * (n - 1)) / ((n - 2) * (n + 1) * (n + 3)))

    return s

end

"""
    ses(n)

Calculate standard error of the skewness.

# Arguments

- `n::Int64`: number of observations

# Returns

- `s::Float64`
"""
function ses(n::Int64)::Float64

    @assert n > 0 "n must be > 0."

    s = sqrt((6 * n * (n - 1)) / ((n - 2) * (n + 1) * (n + 3)))

    return s

end

"""
    sek(x)

Calculate standard error of the kurtosis.

# Arguments

- `x::AbstractVector`

# Returns

- `s::Float64`
"""
function sek(x::AbstractVector)::Float64

    n = length(x)
    @assert n > 0 "x length must be > 0."

    s = 2 * (n - 1) * sqrt((6 * n) / ((n - 2) * (n - 3) * (n + 3) * (n + 5)))

    return s

end

"""
    ses(n)

Calculate standard error of the skewness.

# Arguments

- `n::Int64`: number of observations

# Returns

- `s::Float64`
"""
function sek(n::Int64)::Float64

    @assert n > 0 "n must be > 0."
    
    s = sqrt((6 * n * (n - 1)) / ((n - 2) * (n + 1) * (n + 3)))

    return s

end
