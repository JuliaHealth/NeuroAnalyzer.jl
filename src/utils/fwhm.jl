export fwhm

"""
    fwhm(s)

Calculate the indices of the full-width at half-maximum (FWHM) points of a Gaussian-like signal.

# Arguments

- `s::AbstractVector`: signal vector; must contain at least 2 elements

# Returns

- `Int64`: index of the pre-peak half-maximum point
- `Int64`: index of the signal peak
- `Int64`: index of the post-peak half-maximum point

# Notes

- The input `s` is normalized internally; the original vector is not modified.
- For noisy or non-unimodal signals, `vsearch` may return the index of the closest sample to 0.5 rather than a true half-maximum crossing.
"""
function fwhm(s::AbstractVector)::Tuple{Int64, Int64, Int64}

    # validate
    length(s) >= 2 || throw(ArgumentError("s must contain at least 2 elements."))

    # normalize to [0, 1] so the half-maximum level is always 0.5
    s = normalize_n(s)

    # index of the global peak
    signal_peak_idx = vsearch(maximum(s), s)

    # nearest sample to 0.5 in the pre-peak segment [1 … signal_peak_idx]
    prepeak_hmp = vsearch(0.5, s[1:signal_peak_idx])

    # nearest sample to 0.5 in the post-peak segment [signal_peak_idx … end]
    # offset by signal_peak_idx - 1 to convert the local index back to global
    postpeak_hmp = signal_peak_idx + vsearch(0.5, s[signal_peak_idx:end]) - 1

    return prepeak_hmp, signal_peak_idx, postpeak_hmp
end
