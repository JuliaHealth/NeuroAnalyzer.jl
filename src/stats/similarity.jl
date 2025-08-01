export jaccsim
export sdi
export cosim

"""
    jaccsim(x, y)

Calculate Jaccard similarity between two vectors.

# Arguments

- `x::AbstractVector`
- `y::AbstractVector`

# Returns

- `j::Float64`

# Notes

To compute Jaccard distance, use `1 - jaccard(x, y)`
"""
function jaccsim(x::AbstractVector, y::AbstractVector)::Float64

    i = float(length(intersect(x, y)))
    u = length(x) + length(y) - i
    u == 0 && @error "Length of x + length of y - $i must not be 0."
    j = i / u

    return j

end

"""
    sdi(x, y)

Calculate Sorensen-Dice similarity index between two vectors.

# Arguments

- `x::AbstractVector`
- `y::AbstractVector`

# Returns

- `sdi::Float64`

# Notes

Sorensen-Dice Index values range from 0 to 1:
- 0.80-1.00: very high similarity
- 0.60-0.79: high similarity
- 0.40-0.59: moderate similarity
- 0.20-0.39: low similarity
- 0.00-0.19: very low similarity
"""
function sdi(x::AbstractVector, y::AbstractVector)::Float64

    a = length(x)
    b = length(y)
    c = float(length(intersect(x, y)))
    
    return round((2 * c) / (a + b), digits=2)

end

"""
    cosim(x, y)

Calculate cosine similarity between two vectors.

# Arguments

- `x::AbstractVector`
- `y::AbstractVector`

# Returns

- `cs::Float64`
"""
function cosim(x::AbstractVector, y::AbstractVector)::Float64

    @assert length(x) == length(y) "Lengths of x and y must be equal."

    cs = sum(x .* y) / (sqrt(sum(x.^2)) * sqrt(sum(y.^2)))

    return cs

end
