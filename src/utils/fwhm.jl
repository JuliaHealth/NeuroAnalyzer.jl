export fwhm

"""
    fwhm(s)

Calculate the indices of the full-width at half-maximum (FWHM) points of a Gaussian-like signal.

# Arguments

- `s::AbstractVector`: signal vector; must contain at least 2 elements

# Returns

- `p1_idx::Int64`: index of the pre-peak half-maximum point
- `p_idx::Int64`: index of the signal peak
- `p2_idx::Int64`: index of the post-peak half-maximum point

# Throws

- `ArgumentError`: if `length(s) < 2`.

# Notes

- The input `s` is normalized internally; the original vector is not modified.
- For noisy or non-unimodal signals, `vsearch` may return the index of the closest sample to 0.5 rather than a true half-maximum crossing.

# See also

[`normalize_n`](@ref), [`vsearch`](@ref)
"""
function fwhm(s::AbstractVector)::Tuple{Int64, Int64, Int64}

    @assert length(s) >= 2 "s must contain at least 2 elements."

    # normalize to [0, 1] so the half-maximum level is always 0.5
    s = normalize_n(s)

    # index of the global peak
    p_idx  = vsearch(maximum(s), s)

    # nearest sample to 0.5 in the pre-peak segment [1 … p_idx]
    p1_idx = vsearch(0.5, s[1:p_idx])

    # nearest sample to 0.5 in the post-peak segment [p_idx … end];
    # offset by p_idx - 1 to convert the local index back to global
    p2_idx = p_idx + vsearch(0.5, s[p_idx:end]) - 1

    return p1_idx, p_idx, p2_idx

end
