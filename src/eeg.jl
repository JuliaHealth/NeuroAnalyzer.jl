using FFTW
using Statistics

"""
Gets the position value n in array x
"""
function vsearch(x, n)
    _, n_idx = findmin(abs.(x .- n))
    return n_idx
end

"""
Converts cartographic coordinates to polar
"""
function cart2pol(x, y)
    rho = hypot(x, y)
    theta = atan(y, x)
    return rho, theta
end

"""
Converts polar coordinates to cartographic
"""
function pol2cart(theta, rho)
    x = rho * cos(theta)
    y = rho * sin(theta)
    return x, y
end

"""
Unity based (min..max) data scaling: 0..1
"""
function min_max(x)
    return (x - findmin(x)) / (x - findmax(x))
end

"""
Returns the phase angles, in radians, of a vector with complex elements
"""
function cvangle(v)
    return atan.(imag(v), real(v))
end

"""
Returns the n-point long symmetric Hann window column vector
"""
function hann(n)
    return 0.5 .* (1 .- cos.(2 .* pi .* range(0, 1, length = n)))
end

"""
Calculates Hildebrand rule for symmetry
H < 0.2 means symmetry
"""
function hildebrand_rule(x)
    return (mean(x) - median(x)) ./ std(x)
end

"""
Calculates Jaccard similarity between two vectors
"""
function jaccard_similarity(x, y)
    intersection = length(intersect(x, y))
    union = length(x) + length(y) - intersection
    return intersection / union
end

"""
Calculates Manhattan distance between two vectors
"""
function manhattan_distance(x, y)
    return sum(x .- y)
end

"""
Zero-padded FFT
"""
function fft0(x, padlength)
    return fft(vcat(x, zeros(padlength)))
end

"""
Returns the derivative of the signal with length same as the signal
"""
function signal_derivative(x)
    return np.append(np.diff(x), np.diff(x)[-1])
end

"""
Calculates absolute band power
"""
function band_power(psd, freqs, f1, f2)
    # find intersecting values in frequency vector
    idx_delta = np.logical_and(freqs >= f1, freqs <= f2)
    # dx: frequency resolution
    p = simpson(psd[idx_delta], dx=freqs[1]-freqs[0])
    return p
end

function simpson(x, dx)
    return int
end