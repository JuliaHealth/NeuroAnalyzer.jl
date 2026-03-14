export paired_labels

"""
    paired_labels(l; <keyword arguments>)

Generate all ordered pairwise label combinations from a single label vector.

# Arguments

- `l::Vector{String}`: input label vector; must not be empty
- `unq::Bool=true`: if `true`, exclude self-pairs (`"label1-label1"`)

# Returns

- `Vector{String}`: ordered paired labels

# Throws

- `ArgumentError`: if `l` is empty

# See also

[`paired_labels(::Vector{String}, ::Vector{String})`](@ref)
"""
function paired_labels(l::Vector{String}; unq::Bool = true)::Vector{String}

    !(length(l) > 0) && throw(ArgumentError("l must not be empty."))

    if unq
        # exclude diagonal (self-pairs): n × (n − 1) ordered pairs
        return [l[i] * "-" * l[j] for i in eachindex(l) for j in eachindex(l) if i != j]
    else
        # include all n² ordered pairs
        return [l[i] * "-" * l[j] for i in eachindex(l) for j in eachindex(l)]
    end

end

"""
    paired_labels(l1, l2)

Generate element-wise paired labels from two label vectors of equal length.

# Arguments

- `l1::Vector{String}`: first label vector; must not be empty
- `l2::Vector{String}`: second label vector; must have the same length as `l1`

# Returns

- `Vector{String}`: element-wise paired labels of length `length(l1)`

# Throws
- `ArgumentError`: if `l1` and `l2` have different lengths, or if either is empty

# See also

[`paired_labels(::Vector{String})`](@ref)
"""
function paired_labels(l1::Vector{String}, l2::Vector{String})::Vector{String}

    !(length(l1) > 0) && throw(ArgumentError("l1 must not be empty."))
    !(length(l1) == length(l2)) && throw(ArgumentError("l1 and l2 must have the same length."))

    return l1 .* "-" .* l2

end
