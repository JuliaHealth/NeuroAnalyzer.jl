export jaccsim
export sdi

"""
    jaccsim(x, y)

Calculate the Jaccard similarity between two vectors.

Computed as `|x ∩ y| / |x ∪ y|` where `|x ∪ y| = |x| + |y| − |x ∩ y|`.

Values range from 0 (no common elements) to 1 (identical element sets).

# Arguments

- `x::AbstractVector`: first vector; must not be empty
- `y::AbstractVector`: second vector; must not be empty

# Returns

- `Float64`: Jaccard similarity ∈ `[0, 1]`

# Throws

- `ArgumentError`: if both vectors are empty (union is empty → division by zero)

# Notes

Jaccard distance = `1 − jaccsim(x, y)`.

# See also

[`sdi`](@ref)
"""
function jaccsim(x::AbstractVector, y::AbstractVector)::Float64

    i = length(intersect(x, y))
    u = length(x) + length(y) - i # |x ∪ y|
    @assert u > 0 "Union of x and y is empty; Jaccard similarity is undefined."

    return i / u

end

"""
    sdi(x, y)

Calculate the Sørensen–Dice similarity index between two vectors.

Computed as `2|x ∩ y| / (|x| + |y|)`. Values range from 0 to 1.

# Arguments

- `x::AbstractVector`: first vector; must not be empty
- `y::AbstractVector`: second vector; must not be empty

# Returns

- `Float64`: Sørensen–Dice index ∈ `[0, 1]`, rounded to 2 decimal places

# Throws

- `ArgumentError`: if both vectors are empty (denominator is zero)

# Notes

Interpretation guide:

- `0.80–1.00`: very high similarity
- `0.60–0.79`: high similarity
- `0.40–0.59`: moderate similarity
- `0.20–0.39`: low similarity
- `0.00–0.19`: very low similarity

# See also

[`jaccsim`](@ref)
"""
function sdi(x::AbstractVector, y::AbstractVector)::Float64

    a = length(x)
    b = length(y)
    @assert a + b > 0 "x and y must not both be empty (denominator is zero)."
    c = length(intersect(x, y))

    return round(2c / (a + b), digits=2)

end
