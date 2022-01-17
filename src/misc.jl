"""
    linspace(start, stop, length)

Generates `length`-long sequence of evenly spaced numbers between `start` and `stop`.
"""
linspace(start, stop, length) = collect(range(start, stop, length))

"""
    logspace(start, stop, length)

Generates `length`-long sequence of log10-spaced numbers between `start` and `stop`.
"""
logspace(start, stop, length) = collect(exp10.(range(start, stop, length)))

"""
    zero_pad(m)

Pads the matrix `m` with zeros to make it square.
"""
function zero_pad(m::Matrix)
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
    vsearch(x::Vector, y::Number; return_distance=false)

Returns the positions of the `y` value in the vector `x`.
"""
function vsearch(x::Vector, y::Number; return_distance=false)
    y_dist, y_idx = findmin(abs.(x .- y))
    return_distance == false && return y_idx
    return_distance == true && return y_idx, y_dist
end
"""
    vsearch(x::Vector, y::Vector)

Returns the positions of the `y` vector in the vector `x`.
"""
function vsearch(x::Vector, y::Vector)
    length(y) > length(x) && throw(ArgumentError("Length of 'y' cannot be larger than length 'x'"))
    y_idx = zeros(length(y))
    y_dist = zeros(length(y))
    for idx in 1:length(y)
        y_dist[idx], y_idx[idx] = findmin(abs.(x .- y[idx]))
    end
    return convert.(Int, y_idx), y_dist
end

"""
    cart2pol(x, y)

Converts cartographic coordinates `x` and `y` to polar.
"""
function cart2pol(x, y)
    rho = hypot(x, y)
    theta = atan(y, x)
end

"""
    pol2cart(theta, rho)

Converts polar coordinates `theta` and `rho` to cartographic.
"""
function pol2cart(theta, rho)
    x = rho * cos(theta)
    y = rho * sin(theta)
    return x, y
end

"""
    minmax_scaler(x)

Apply unity based (min - max) data scaling to the vector `x`. After scaling, `x` will have values between 0 and 1.
"""
minmax_scaler(x::Vector) = (x - findmin(x)[1]) / (x - findmax(x)[1])

"""
    cvangle(x)

Returns the phase angles, in radians, of the vector `x` with complex elements.
"""
cvangle(x::Vector) = atan.(imag(x), real(x))

"""
    hann(n)

Returns the `n`-point long symmetric Hanning window.
"""
hann(n::Int) = 0.5 .* (1 .- cos.(2 .* pi .* range(0, 1, length = n)))

"""
    hildebrand_rule(x)

Calculates Hildebrand rule for vector `x`.
If H < 0.2 then the vector `x` is symmetrical.
"""
hildebrand_rule(x::Vector) = (mean(x) - median(x)) ./ std(x)

"""
    jaccard_similarity(x, y)

Calculates Jaccard similarity between two vectors `x` and `y`.
"""
function jaccard_similarity(x::Vector, y::Vector)
    intersection = length(intersect(x, y))
    union = length(x) + length(y) - intersection
    return intersection / union
end

"""
    fft0(x, n)

Calculates FFT for the vector `x` padded with `n` zeros at the end.
"""
fft0(x::Vector, padlength::Int) = fft(vcat(x, zeros(padlength)))

"""
    nexpow2(x)

Returns the next power of 2 for given number `x`.
"""
function nexpow2(x)
    x == 0 && return 1
    x == 0 || return 2 ^ ndigits(x - 1, base=2)
end

"""
    vsplit(x, n)

Splits the vector `x` into `n`-long pieces.
"""
function vsplit(x::Vector, n::Int=1)
    length(x) % n == 0 || throw(ArgumentError("""Length of "x" must be a multiple of "n"."""))
    x_m = reshape(x, length(x) ÷ n, n)
    result = [x_m[1, :]]
    for idx in 2:size(x_m, 1)
        result = vcat(result, [x_m[idx, :]])
    end
    return result
end

"""
    rms(x)

Calculates Root Mean Square of the vector `x`.
"""
rms(x::Vector) = norm(x) / sqrt(length(x))

"""
    db(x)

Converts values of the vector `x` to dB.
"""
db(x::Vector) = 10 .* log10.(x ./ findmax(x)[1])

"""
    sine(f, t, a, p)

Generates sine wave of `f` frequency over `t` time; optional arguments are: `a` amplitude and  `p` phase.
"""
sine(f, t, a=1, p=0) = a .* sin.(2 * pi .* f * t .+ p)

"""
    frequencies(t)

Generates vector of frequencies for given time vector `t`.
"""
function frequencies(t::Vector)
    # sampling interval
    dt = t[2] - t[1]
    # sampling rate
    fs = 1 / dt
    # frequency step size
    df = 1 / (length(t) * dt)
    # Nyquist frequency
    nyquist_freq = fs / 2
    # frequency array
    hz = collect(0:df:nyquist_freq)
    return hz
end