"""
    _dict2labeled_matrix(d; rev)

Convert a dictionary to labeled matrix format (labels and values vectors).

# Arguments

- `d::Dict{String, Vector{Float64}}`: input dictionary with string keys and Float64 vectors as values
- `rev::Bool=true`: if `true`, reverse the order of labels and values

# Returns

- `Tuple{Vector{String}, Vector{Vector{Float64}}}`: a tuple containing:
    - vector of labels (strings)
    - vector of values (vectors of Float64)
"""
function _dict2labeled_matrix(
        d::Dict;
        rev::Bool = true,
    )::Tuple{Vector{String}, Vector{Vector{Float64}}}
    isempty(d) && throw(ArgumentError("Dictionary cannot be empty."))

    # extract labels and values
    l = collect(keys(d))
    v = collect(values(d))

    # validate all values are vectors of Float64
    all(x -> x isa Vector{Float64}, v) ||
        throw(ArgumentError("All dictionary values must be vectors of Float64"))

    # return in requested order
    if rev
        return reverse(l), reverse(v)
    else
        return l, v
    end
end

"""
    _labeled_matrix2dict(l, v)

Convert labeled matrix format (labels and values vectors) back to a dictionary.

# Arguments

- `l::Vector{String}`: vector of labels (strings)
- `v::Vector{Vector{Float64}}`: Vector of values (vectors of Float64)

# Returns

`Dict{String, Vector{Float64}}`: dictionary mapping labels to their corresponding values

"""
function _labeled_matrix2dict(l::Vector{String}, v::Vector{Vector{Float64}})::Dict
    length(l) == length(v) ||
        throw(ArgumentError("Labels and values vectors must have equal lengths."))
    isempty(l) && "Input vectors cannot be empty."
    isempty(v) && "Input vectors cannot be empty."
    # convert to dictionary
    return Dict(zip(l, v))
end
