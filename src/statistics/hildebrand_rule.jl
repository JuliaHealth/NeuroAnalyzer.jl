export hildebrand_rule

"""
    hildebrand_rule(x)

Calculate Hildebrand rule for vector `x`.
If H < 0.2 then the vector `x` is symmetrical.

# Arguments

- `x::AbstractVector`

# Returns

- `h::Float64`
"""
function hildebrand_rule(x::AbstractVector)

    return (mean(x) - median(x)) ./ std(x)
    
end
