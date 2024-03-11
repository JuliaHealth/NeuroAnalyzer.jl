export paired_labels

"""
    paired_labels(l; unq)

Return paired labels.

# Arguments

- `l::Vector{String}`
- `unq::Bool=true`: if true, do not add pairs of the same labels, e.g. "label1-label1"

# Returns

- `l_paired::Vector{String}`: paired labels
"""
function paired_labels(l::Vector{String}; unq::Bool=true)

    if unq == true
        l_paired = repeat([""], length(l)^2 - length(l))
    else
        l_paired = repeat([""], length(l)^2)
    end

    idx = 1
    for idx1 in 1:length(l), idx2 in 1:length(l)
        if unq == true
            if idx1 != idx2
                l_paired[idx] = l[idx1] * "-" * l[idx2]
                idx += 1
            end
        else
            l_paired[idx] = l[idx1] * "-" * l[idx2]
            idx += 1
        end
    end

    return l_paired

end

"""
    paired_labels(l1, l2)

Return paired labels.

# Arguments

- `l1::Vector{String}`
- `l2::Vector{String}`

# Returns

- `l_paired::Vector{String}`: paired labels
"""
function paired_labels(l1::Vector{String}, l2::Vector{String})

    @assert length(l1) == length(l2) "l1 and l2 length must be equal."
    l_paired = l1 .* "-" .* l2

    return l_paired

end
