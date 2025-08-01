export cvm
export cvmd
export fano

"""
    cvm(x)

Calculate coefficient of variation for a mean.

# Arguments

- `x::AbstractVector`

# Returns

- `cvm::Float64`
"""
function cvm(x::AbstractVector)::Float64

    @assert mean(x) != 0 "mean(x) must not be equal 0."

    return std(x) / mean(x)

end

"""
    cvmd(x)

Calculate coefficient of variation for a median.

# Arguments

- `x::AbstractVector`

# Returns

- `cvmd::Float64`
"""
function cvmd(x::AbstractVector)::Float64

    @assert median(x) != 0 "median(x) must not be equal 0."

    return ((quantile(x, 0.75) - quantile(x, 0.25)) / 2) / median(x)

end

"""
    fano(x)

Calculate Fano factor.

# Arguments

- `x::AbstractVector`

# Returns

- `fano::Float64`
"""
function fano(x::AbstractVector)::Float64

    @assert mean(x) != 0 "mean(x) must not be equal 0."

    return var(x) / mean(x)

end
