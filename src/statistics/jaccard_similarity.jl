export jaccard_similarity

"""
    jaccard_similarity(x, y)

Calculate Jaccard similarity between two vectors.

# Arguments

- `x::AbstractVector`
- `y::AbstractVector`

# Returns

- `j::Float64`
"""
function jaccard_similarity(x::AbstractVector, y::AbstractVector)

    i = float(length(intersect(x, y)))
    u = length(x) + length(y) - i

    return i / u

end

