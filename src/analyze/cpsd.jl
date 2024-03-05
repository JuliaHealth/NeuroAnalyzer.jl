export cpsd

"""
   cpsd(s1, s2)

Calculate cross power spectral density.

# Arguments

- `s1::AbstractVector`
- `s2::AbstractVector`

# Returns

- `pxy::AbstractVector`
"""
function cpsd(s1::AbstractVector, s2::AbstractVector; fs::Int64)

    s1, s2 = _veqlen(s1, s2)

    # multitaper
    # demean::Bool
    # freq_range::T5,
    s = hcat(s1, s2)'
    pxy = mt_cross_power_spectra(s, fs=fs)
    p = power(pxy)
    f = freq(pxy)

    abs.(p)

    # fft
    pxy = conj.(rfft(s1) .* rfft(s2))
    p = pxy 
    f, _ = freqs(s1, fs)

end
