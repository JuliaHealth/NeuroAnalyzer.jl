export op

"""
    op(x, y)

Calculate the outer product of two vectors.

Computed as `x × yᵀ`, producing a matrix of shape `(length(x) × length(y))`.

# Arguments

- `x::AbstractVector`: left vector; must not be empty
- `y::AbstractVector`: right vector; must not be empty

# Returns

- `Matrix`: outer product matrix of shape `(length(x), length(y))`

# Throws

- `ArgumentError`: If either vector is empty.
"""
function op(x::AbstractVector, y::AbstractVector)::AbstractMatrix

    @assert length(x) > 0 "x must not be empty."
    @assert length(y) > 0 "y must not be empty."

    return x * y'

end
