export sem_diff

"""
    sem_diff(x::AbstractVector, y::AbstractVector)

Calculate SEM (standard error of the mean) for the difference of two means.

# Arguments

- `x::AbstractVector`
- `y::AbstractVector`

# Returns

- `sem_diff::Float64`
"""
function sem_diff(x::AbstractVector, y::AbstractVector)

    length(x) == length(y) || throw(ArgumentError("Both vectors must have the same length."))

    return sqrt((std(x)^2 / sqrt(length(x))) + (std(y)^2 / sqrt(length(y))))
    
end
