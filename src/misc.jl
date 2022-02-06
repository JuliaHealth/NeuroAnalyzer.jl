"""
    linspace(start, stop, length)

Generates `length`-long sequence of evenly spaced numbers between `start` and `stop`.

# Arguments

- `start::Union{Int64, Float64}`
- `stop::Union{Int64, Float64}`
- `length::Int64`

# Returns

- `range::Number`
"""
linspace(start::Union{Int64, Float64}, stop::Union{Int64, Float64}, length::Int64) = collect(range(start, stop, length))

"""
    logspace(start, stop, length)

Generates `length`-long sequence of log10-spaced numbers between `start` and `stop`.

# Arguments

- `start::Union{Int64, Float64}`
- `stop::Union{Int64, Float64}`
- `length::Int64`

# Returns

- `range::Number`
"""
logspace(start::Union{Int64, Float64}, stop::Union{Int64, Float64}, length::Int64) = collect(exp10.(range(start, stop, length)))

"""
    zero_pad(m)

Pads the matrix `m` with zeros to make it square.

# Arguments

- `m::Union{Matrix{Int64}, Matrix{Float64}, Matrix{ComplexF64}}`

# Returns

- `m::Union{Matrix{Int64}, Matrix{Float64}, Matrix{ComplexF64}}`
"""
function zero_pad(m::Union{Matrix{Int64}, Matrix{Float64}, Matrix{ComplexF64}})
    nr, nc = size(m)
    if nr > nc
        mp = repeat([0], nr, nr - nc)
        mp = hcat(m, mp)
    elseif nr < nc
        mp = repeat([0], nc - nr, nc)
        mp = vcat(m, mp)
    elseif nr == nc
        mp = m
    end
    return mp
end

"""
    vsearch(x, y; return_distance=false)

Returns the positions of the `y` value in the vector `x` and the difference between `y` and `x[vsearch(x, y)].

# Arguments

- `x::Union{Vector{Int64}, Vector{Float64}}
- `y::Union{Int64, Float64}`
- `return_distance::Bool`

# Returns

y_idx::Int64
y_dist::Union{Int64, Float64}
"""
function vsearch(x::Union{Vector{Int64}, Vector{Float64}}, y::Union{Int64, Float64}; return_distance::Bool=false)
    y_dist, y_idx = findmin(abs.(x .- y))
    return_distance == false && return y_idx
    return_distance == true && return y_idx, y_dist
end

"""
    vsearch(x, y)

Returns the positions of the `y` vector in the vector `x`.

# Arguments

- `x::Union{Vector{Int64}, Vector{Float64}}`
- `y::Union{Vector{Int64}, Vector{Float64}}
- `return_distance::Bool`

# Returns

y_idx::Int64
y_dist::Union{Int64, Float64}
"""
function vsearch(x::Union{Vector{Int64}, Vector{Float64}}, y::Union{Vector{Int64}, Vector{Float64}}; return_distance=false)
    length(y) > length(x) && throw(ArgumentError("Length of 'y' cannot be larger than length 'x'"))

    y_idx = zeros(length(y))
    y_dist = zeros(length(y))
    for idx in 1:length(y)
        y_dist[idx], y_idx[idx] = findmin(abs.(x .- y[idx]))
    end
    return_distance == false && return convert.(Int64, y_idx)
    return_distance == true && return convert.(Int64, y_idx), y_dist
end

"""
    cart2pol(x, y)

Converts cartographic coordinates `x` and `y` to polar.

# Arguments

- `x::Union{Int64, Float64}`
- `y::Union{Int64, Float64}`

# Returns

- `rho::Float64`
- `theta::Float64`
"""
function cart2pol(x::Union{Int64, Float64}, y::Union{Int64, Float64})
    rho = hypot(x, y)
    theta = atan(y, x)
    return rho, theta
end

"""
    pol2cart(theta, rho)

Converts polar coordinates `theta` and `rho` to cartographic.

# Arguments

- `rho::Float64`
- `theta::Float64`

# Returns

- `x::Float64`
- `y::Float64`
"""
function pol2cart(theta::Float64, rho::Float64)
    x = rho * cos(theta)
    y = rho * sin(theta)
    return x, y
end

"""
    generate_hanning(n)

Returns the `n`-point long symmetric Hanning window.

# Arguments

- `n::Int64`

# Returns

- `hanning::Vector{Float64}`
"""
generate_hanning(n::Int64) = 0.5 .* (1 .+ cos.(2 .* pi .* range(0, 1, length = n)))

"""
    hildebrand_rule(x)

Calculates Hildebrand rule for vector `x`.
If H < 0.2 then the vector `x` is symmetrical.

# Arguments

- `n::Int64`

# Returns

- `h::Float64`
"""
hildebrand_rule(x::Union{Vector{Int64}, Vector{Float64}}) = (mean(x) - median(x)) ./ std(x)

"""
    jaccard_similarity(x, y)

Calculates Jaccard similarity between two vectors `x` and `y`.

# Arguments

- `n::Int64`

# Returns

- `j::Float64`
"""
function jaccard_similarity(x::Union{Vector{Int64}, Vector{Float64}}, y::Union{Vector{Int64}, Vector{Float64}})
    intersection = length(intersect(x, y))
    union = length(x) + length(y) - intersection
    j = intersection / union
    return j
end

"""
    fft0(x, n)

Calculates FFT for the vector `x` padded with `n` or `n - length(x)` zeros at the end.

# Arguments

- `x::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`
- `n::Int64`

# Returns

- `fft0::Vector{ComplexF64}`
"""
function fft0(x::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}, n::Int64)
    n < 0 && throw(ArgumentError("Pad must be positive."))

    n > length(x) && (n = n - length(x))
    return fft(vcat(x, zeros(eltype(x), n)))
end

"""
    ifft0(x, n)

Calculates IFFT for the vector `x` padded with `n` or `n - length(x)` zeros at the end.

# Arguments

- `x::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`
- `n::Int64`

# Returns

- `ifft0::Vector{ComplexF64}`
"""
function ifft0(x::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}, n::Int64)
    n < 0 && throw(ArgumentError("Pad must be positive."))

    n > length(x) && (n = n - length(x))
    return ifft(vcat(x, zeros(eltype(x), n)))
end

"""
    nextpow2(x)

Returns the next power of 2 for given number `x`.

# Argument

 - `x::Int64`

# Returns

- `nextpow::Int64`
"""
function nextpow2(x::Int64)
    x == 0 && return 1
    x == 0 || return 2 ^ ndigits(x - 1, base=2)
end

"""
    vsplit(x, n)

Splits the vector `x` into `n`-long pieces.

# Argument

- `x::Union{Vector{Int64}, Vector{Float64}}`
 - `n::Int64`

# Returns

- `x::Vector{Union{Vector{Int64}, Vector{Float64}}}`
"""
function vsplit(x::Union{Vector{Int64}, Vector{Float64}}, n::Int64=1)
    n < 0 && throw(ArgumentError("n must be positive."))
    length(x) % n == 0 || throw(ArgumentError("Length of x must be a multiple of n."))

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

# Arguments

- `x::Union{Vector{Int64}, Vector{Float64}}`

# Returns

- rms::Float64`
"""
rms(x::Union{Vector{Int64}, Vector{Float64}}) = norm(x) / sqrt(length(x))

"""
    generate_sine(f, t, a, p)

Generates sine wave of `f` frequency over `t` time; optional arguments are: `a` amplitude and  `p` phase.

# Arguments

- `f::Int64`
- `t::Union{Vector{Int64}, Vector{Float64}}`
- `a::Union{Int64, Float64}`
- `p::Union{Int64, Float64}`

# Returns

- sine::Vector{Float64}`
"""
generate_sine(f::Int64, t::Union{Vector{Int64}, Vector{Float64}}, a::Union{Int64, Float64}=1, p::Union{Int64, Float64}=0) = @. a * sin(2 * pi * f * t + p)

"""
    freqs(t)

Returns vector of frequencies and Nyquist frequency for given time vector `t`.

# Arguments

- `t::Union{Vector{Int64}, Vector{Float64}, AbstractRange}`

# Returns

- `hz::Vector{Float64}
- `nyquist_freq::Float64`
"""
function freqs(t::Union{Vector{Int64}, Vector{Float64}, AbstractRange})
    if typeof(t) <: AbstractRange
        t = collect(t)
    end

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

    return hz, nyquist_freq
end

"""
    freqs(signal, fs)

Returns vector of frequencies and Nyquist frequency for given `signal` and `fs`.

# Arguments

- `signal::Vector{Float64}`
- `fs::Union{Int64, Float64}`

# Returns

- `hz::Vector{Float64}
- `nyquist_freq::Float64`
"""
function freqs(signal::Vector{Float64}, fs::Union{Int64, Float64})
    fs < 0 && throw(ArgumentError("Sampling rate must be positive."))

    # Nyquist frequency
    nyquist_freq = fs / 2
    # frequency array
    hz = linspace(0, nyquist_freq, floor(Int64, length(signal) / 2) + 1)

    return hz, nyquist_freq
end

"""
    matrix_sortperm(m; rev=false, dims=1)

Generates sorting index for matrix `m` by columns (`dims` = 1) or by rows (`dims` = 2).

# Arguments

- `m::Matrix`
- `rev::Bool`
- `dims::Int64`

# Returns

- `m_idx::Matrix{Int64}`
"""
function matrix_sortperm(m::Matrix; rev::Bool=false, dims::Int64=1)
    (dims < 1 || dims > 2) && throw(ArgumentError("dims must be 1 or 2."))
    m_idx = zeros(Int, size(m))
    idx=1
    if dims == 1
        for idx = 1:size(m, 2)
            # sort by columns
            m_idx[:, idx] = sortperm(m[:, idx], rev=rev)
        end
    else
        for idx = 1:size(m, 1)
            # sort by rows
            m_idx[idx, :] = sortperm(m[idx, :], rev=rev)'
        end     
    end

    return m_idx
end

"""
    matrix_sort(m, m_idx; rev=false, dims=1)

Sorts matrix `m` using sorting index `m_idx` by columns (`dims` = 1) or by rows (`dims` = 2).

# Arguments

- `m::Matrix`
- `m_idx::Vector{Int64}`
- `rev::Bool`
- `dims::Int64`

# Returns

- `m_sorted::Matrix`
"""
function matrix_sort(m::Matrix, m_idx::Vector{Int64}; rev::Bool=false, dims::Int64=1)
    (dims < 1 || dims > 2) && throw(ArgumentError("dims must be 1 or 2."))

    m_sorted = zeros(eltype(m), size(m))
    if dims == 1
        for idx = 1:size(m, 2)
            # sort by columns
            tmp = m[:, idx]
            tmp = tmp[m_idx]
            m_sorted[:, idx] = tmp
        end
    else
        for idx = 1:size(m, 1)
            # sort by rows
            tmp = m[idx, :]
            tmp = tmp[m_idx]
            m_sorted[idx, :] = tmp
        end
    end
    return m_sorted
end

"""
    pad0(x, n)

Pads the vector `x` with `n` zeros at the beginning and at the end.

# Arguments

- `x::Union{Vector{Int64}, Vector{Float64}}`
- `n::Int64`

# Returns

- `v_pad::Union{Vector{Int64}, Vector{Float64}}`
"""
# to do: check if x is numeric vector
function pad0(x::Union{Vector{Int64}, Vector{Float64}}, n)
    n < 1 && throw(ArgumentError("n must be positive."))

    v_pad = vcat(zeros(eltype(x), n), x, zeros(eltype(x), n))
end
"""
    hz2rads(f)

Converts frequency `f` in Hz to rad/s.

# Arguments

- `f::Union{Int64, Float64}`

# Returns

- `f_rads::Float64`
"""
hz2rads(f::Union{Int64, Float64}) = 2 * pi * f

"""
    rads2hz(f)

Converts frequency `f` in rad/s to Hz.

# Arguments

- `f::Union{Int64, Float64}`

# Returns

- `f_rads::Float64`
"""
rads2hz(f::Union{Int64, Float64}) = f / 2 * pi

"""
    z_score(x)

Calculates Z-scores for each value of the vector `x`.

# Arguments

- `x::Union{Vector{Int64}, Vector{Float64}}`

# Returns

- `z_score::Vector{Float64}`
"""
z_score(x::Union{Vector{Int64}, Vector{Float64}}) = (x .- mean(x)) ./ std(x)

"""
    k(n)

Calculates number of categories for a given sample size `n`.

# Arguments

- `n::Int64`

# Returns

- `k::Float64`
"""
k(n::Int64) = (sqrt(n), (1 + 3.222 * log10(n)))

"""
    cmax(x)

Returns maximum value of the complex vector`x`.

# Arguments

- `x::Vector{ComplexF64}`

# Returns

- `cmax::ComplexF64`
"""
cmax(x::Vector{ComplexF64}) = argmax(abs, x)

"""
    cmin(x)

Returns minimum value of the complex vector`x`.

# Arguments

- `x::Vector{ComplexF64}`

# Returns

- `cmin::ComplexF64`
"""
cmin(x::Vector{ComplexF64}) = argmin(abs, x)

"""
    generate_sinc(t=0:0.01:10, f=10, peak=4)

Generates sinc function.

# Arguments

- `t::AbstractRange` - time
- `f::Union{Int64, Float64}` - frequency
- `peak::Union{Int64, Float64}` - sinc peak time

# Returns

- `sinc::Vector{Float64}
"""
function generate_sinc(t::AbstractRange=-2:0.01:2, f::Union{Int64, Float64}=10.0, peak::Union{Int64, Float64}=0)
    y_sinc = @. sin(2 * pi * f * (t - peak)) / (t - peak)
    nan_idx = y_sinc[y_sinc .== NaN]
    y_sinc[findall(isnan, y_sinc)[1]] = (y_sinc[findall(isnan, y_sinc)[1] - 1] + y_sinc[findall(isnan, y_sinc)[1] + 1]) / 2
    return y_sinc
end

"""
    generate_morlet(fs, wt, wf)

Generates Morlet wavelet.

# Arguments

- `fs::Int64` - sampling rate
- `wt::Union{Int64, Float64}` - length = -wt:1/fs:wt
- `wf::Union{Int64, Float64}` - frequency
- `ncyc::Int64` - number of cycles
- `complex::Bool` - generate complex Morlet

# Returns

- `morlet::Union{Vector{Float64}, Vector{ComplexF64}}
"""
function generate_morlet(fs::Int64, wt::Union{Int64, Float64}, wf::Union{Int64, Float64}; ncyc::Int64=5, complex::Bool=false)
    wt = -wt:1/fs:wt
    complex == false && (sin_wave = @. cos(2 * pi * wf * wt))           # for symmetry at x = 0
    complex == true && (sin_wave = @. exp(im * 2 * pi * wf * wt))       # for symmetry at x = 0
    w = 2 * (ncyc / (2 * pi * wf))^2                                    # ncyc: time-frequency precision
    gaussian = @. exp((-wt.^2) / w)
    morlet_wavelet = sin_wave .* gaussian
    return morlet_wavelet
end