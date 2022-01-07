"""
    linspace(start, stop, length)

Generates 'length' of evenly spaced numbers between 'start' and 'end'.
"""
linspace(start, stop, length) = collect(range(start, stop, length))

"""
    logspace(start, stop, length)

Generates 'length' of log10-spaced numbers between 'start' and 'end'.
"""
logspace(start, stop, length) = collect(exp10.(range(start, stop, length)))

"""
    zero_pad(m)

Pads the matrix 'm' with zeros to make 'm' square
"""
function zero_pad(m)
    nr, nc = size(m)
    if nr > nc
        mp = repeat([0], nr, nr - nc)
        mp = hcat(M, Mp)
    elseif nr < nc
        mp = repeat([0], nc - nr, nc)
        mp = vcat(m, mp)
    elseif nr == nc
        mp = m
    end
    return mp
end

"""
    vsearch(x, n) 

Gets the position of value 'n' in array 'x'.
"""
function vsearch(x, n)
    _, n_idx = findmin(abs.(x .- n))
    return n_idx
end

"""
    cart2pol(x, y)

Converts cartographic coordinates 'x' and 'y' to polar.
"""
function cart2pol(x, y)
    rho = hypot(x, y)
    theta = atan(y, x)
    return rho, theta
end

"""
    pol2cart(theta, rho)

Converts polar coordinates 'theta' and 'rho' to cartographic.
"""
function pol2cart(theta, rho)
    x = rho * cos(theta)
    y = rho * sin(theta)
    return x, y
end

"""
    minmax_scaler(x)

Apply unity based (min - max) data scaling to vector 'x'
'x' will have values between 0 and 1.
"""
minmax_scaler(x::Vector) = (x - findmin(x)) / (x - findmax(x))

"""
    cvangle(v)

Returns the phase angles, in radians, of vector 'v' with complex elements.
"""
cvangle(v::Vector) = atan.(imag(v), real(v))

"""
    hann(n)

Returns the 'n'-point long symmetric Hann window column vector.
"""
hann(n::Int) = 0.5 .* (1 .- cos.(2 .* pi .* range(0, 1, length = n)))

"""
    hildebrand_rule(x)

Calculates Hildebrand rule for vector 'x'.
If H < 0.2 then vector is symmetrical.
"""
hildebrand_rule(x::Vector) = (mean(x) - median(x)) ./ std(x)

"""
    jaccard_similarity(x, y)

Calculates Jaccard similarity between two vectors 'x' and 'y'.
"""
function jaccard_similarity(x::Vector, y::Vector)
    intersection = length(intersect(x, y))
    union = length(x) + length(y) - intersection
    return intersection / union
end

"""
    fft0(x, n)

Calculates FFT for vector 'x' padded with 'n' zeros at the end.
"""
fft0(x::Vector, padlength::Int) = fft(vcat(x, zeros(n)))