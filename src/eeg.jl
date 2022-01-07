"""
    signal_derivative(x)

Returns the derivative of the signal with length same as the signal
"""
signal_derivative(x::Vector) = vcat(diff(x), diff(x)[end-1])

"""
    band_power(psd, f1, f2)

Calculates absolute band power between frequencies 'f1' and 'f2'
"""
function band_power(psd, f1, f2)
    frq_idx = [vsearch(psd.freq, f1), vsearch(psd.freq, f2)]
    # dx: frequency resolution
    dx = psd.freq[2] - psd.freq[1]
    result = simpson(psd.power[frq_idx[1]:frq_idx[2]], start=frq_idx[1], stop=frq_idx[1], dx=dx)
    return result
end

"""
    make_spectrum(y, fs)

Return 'y' signal FFT and DFT sample frequencies for a DFT of the 'y' length 
"""
function make_spectrum(y, fs)
    hs = fft(y)
    n = length(y)               # number of samples
    d = 1/fs                    # time between samples
    fs = fftfreq(n, d)
    return hs, fs
end