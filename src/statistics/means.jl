export meang
export meanh
export meanw

"""
    meang(x)

Calculate geometric mean.

# Arguments

- `x::AbstractVector`

# Returns

- `m::Float64`
"""
function meang(x::AbstractVector)
    return exp(mean(log.(x[x .> 0])))
end

"""
    meanh(x)

Calculate harmonic mean.

# Arguments

- `x::AbstractVector`

# Returns

- `m::Float64`
"""
function meanh(x::AbstractVector)
    return length(x) / sum(1 ./ x)
end

"""
    meanw(x, w)

Calculate weighted mean.

# Arguments

- `x::AbstractVector`
- `w::AbstractVector`: weights

# Returns

- `m::Float64`
"""
function meanw(x::AbstractVector, w::AbstractVector)
    length(x) == length(w) || throw(ArgumentError("Weights and values vectors must have the same length."))
    return length(x) / sum(1 ./ x)
end
