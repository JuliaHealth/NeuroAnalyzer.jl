export fwhm

"""
    fwhm(s)

Calculate indices of full-width half-maximum points of a Gaussian-like distribution.

# Arguments

- `s::AbstractVector`

# Returns

- `p1_idx::Int64`: pre-peak half-maximum point
- `p_idx::Int64`: peak
- `p2_idx::Int64`: post-peak half-maximum point
"""
function fwhm(s::AbstractVector)::Tuple{Int64, Int64, Int64}

    s = normalize_n(s)

    # peak
    p_idx = vsearch(maximum(s), s)
    # pre-peak half-maximum width
    p1_idx = vsearch(0.5, s[1:p_idx])
    # post-peak half-maximum width
    p2_idx = p_idx + vsearch(0.5, s[p_idx:end]) - 1

    return p1_idx, p_idx, p2_idx

end
