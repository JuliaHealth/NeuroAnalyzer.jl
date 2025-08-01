export hildebrand_rule

"""
    hildebrand_rule(x; verbose)

Calculate Hildebrand rule for vector `x`.

If H < 0.2 then the vector `x` is symmetric.

# Arguments

- `x::AbstractVector`
- `verbose::Bool=true`: print detailed output

# Returns

- `h::Float64`
"""
function hildebrand_rule(x::AbstractVector; verbose::Bool=true)::Float64

    @assert std(x) !=0 "Standard deviation of x must not be 0."

    h = (mean(x) - median(x)) / std(x)
    h < 0.2 && (verbose && println("H < 0.2, x is symmetric"))

    return h

end
