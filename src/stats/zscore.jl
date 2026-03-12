export zscore

"""
    zscore(x)

Calculate z-scores for each value of the vector `x`.

# Arguments

- `x::AbstractVector`

# Returns

- `z::Vector{Float64}`
"""
function zscore(x::AbstractVector)::Vector{Float64}

    m = mean(x)
    s = std(x)
    z = (x .- m) ./ s

    return z

end

"""
    zscore(x, m, sd)

Calculate z-score for `x`.

# Arguments

- `x::Real`
- `m::Real`: mean
- `sd::Real`: standard deviation

# Returns

- `z::Float64`
"""
function zscore(x::Real, m::Real, sd::Real)::Float64

    z = (x - m) / sd

    return z

end
