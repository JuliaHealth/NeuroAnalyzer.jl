export op
export angle

"""
    op(x, y)

Calculate outer product of two vectors.

# Arguments

- `x::AbstractVector`
- `y::AbstractVector`

# Returns

- `op::AbstractMatrix`
"""
function op(x::AbstractVector, y::AbstractVector)::AbstractMatrix
    return x * y'
end

"""
    angle(x)

Return the phase angles (in radians) of a vector with complex elements.

# Arguments

- `x::Vector{ComplexF64}`

# Returns

- `angle::Vector{Float64}`
"""
function angle(x::Vector{ComplexF64})::Vector{Float64}
    return atan.(imag.(x), real.(x))
end
