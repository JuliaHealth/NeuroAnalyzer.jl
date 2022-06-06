################################
#                              #
# Low-level internal functions #
#                              #
################################

_reflect(signal::AbstractArray) = vcat(signal[end:-1:1], signal, signal[end:-1:1])
_chop(signal::AbstractArray) = signal[(length(signal) ÷ 3 + 1):(length(signal) ÷ 3) * 2]

################################

"""
    linspace(start, stop, length)

Generates `length`-long sequence of evenly spaced numbers between `start` and `stop`.

# Arguments

- `start::Real`
- `stop::Real`
- `length::Int64`

# Returns

- `range::Number`
"""
function linspace(start::Real, stop::Real, length::Int64)

    return collect(range(start, stop, length))
end

"""
    logspace(start, stop, length)

Generates `length`-long sequence of log10-spaced numbers between `start` and `stop`.

# Arguments

- `start::Real`
- `stop::Real`
- `length::Int64`

# Returns

- `range::Number`
"""
function logspace(start::Real, stop::Real, length::Int64)

    return collect(exp10.(range(start, stop, length)))
end

"""
    m_pad0(m)

Pad the matrix `m` with zeros to make it square.

# Arguments

- `m::Matrix{<:Number}`

# Returns

- `m::Matrix{Number}`
"""
function m_pad0(m::Matrix{<:Number})

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
    vsearch(y, x; return_distance=false)

Return the positions of the `y` value in the vector `x` and the difference between `y` and `x[vsearch(x, y)].

# Arguments

- `y::Real`
- `x::Vector{<:Real}`
- `return_distance::Bool`

# Returns

- `y_idx::Int64`
-` y_dist::Real`
"""
function vsearch(y::Real, x::Vector{<:Real}; return_distance::Bool=false)

    y_dist, y_idx = findmin(abs.(x .- y))

    return_distance == false && (return y_idx)
    return_distance == true && (return y_idx, y_dist)
end

"""
    vsearch(y, x; return_distance=false)

Return the positions of the `y` vector in the vector `x`.

# Arguments

- `x::Vector{<:Real}`
- `y::Vector{<:Real}`
- `return_distance::Bool`

# Returns

- `y_idx::Int64`
- `y_dist::Real`
"""
function vsearch(y::Vector{<:Real}, x::Vector{<:Real}; return_distance=false)

    length(y) > length(x) && throw(ArgumentError("Length of 'y' cannot be larger than length 'x'"))

    y_idx = zeros(length(y))
    y_dist = zeros(length(y))

    for idx in 1:length(y)
        y_dist[idx], y_idx[idx] = findmin(abs.(x .- y[idx]))
    end

    return_distance == false && (return convert.(Int64, y_idx))
    return_distance == true && (return convert.(Int64, y_idx), y_dist)
end

"""
    cart2pol(x, y)

Convert cartographic coordinates `x` and `y` to polar.

# Arguments

- `x::Real`
- `y::Real`

# Returns

- `phi::Float64`
- `theta::Float64`
"""
function cart2pol(x::Real, y::Real)

    return hypot(x, y), atan(y, x)
end

"""
    pol2cart(theta, phi)

Convert polar coordinates `theta` and `phi` to cartographic.

# Arguments

- `phi::Real`
- `theta::Real`

# Returns

- `x::Float64`
- `y::Float64`
"""
function pol2cart(theta::Real, phi::Real)

    return phi * cos(theta), phi * sin(theta)
end

"""
    pol2cart_sph(rho, theta, phi=0)

Convert spherical coordinates `theta` and `phi` and `rho` to cartographic.

# Arguments

- `phi::Real`: the angle with respect to the z-axis (elevation)
- `theta::Real`: the angle in the xy plane with respect to the x-axis (azimuth)
- `rho::Real`: the distance from the origin to the point

# Returns

- `x::Float64`
- `y::Float64`
- `z::Float64`
"""
function sph2cart(rho::Real, theta::Real, phi::Real=0)

    return rho * cos(phi) * cos(theta), rho * cos(phi) * sin(theta), rho * sin(phi)
end

"""
    generate_window(type, n; even)

Return the `n`-point long symmetric window `type`.

# Arguments

- `type::Symbol`: window type:
    - `:hann`: Hann
    - `:bh`: Blackman-Harris
    - `:bohman`: Bohman
    - `:flat`: Flat-top window
    - `:bn`: Blackman-Nuttall
    - `:nutall`: Nuttall
    - `:triangle`: symmetric triangle (left half ↑, right half ↓)
    - `:exp`: symmetric exponential (left half ↑, right half ↓)
- `n::Int64`: window length
- `even::Bool=false`: if true, make the window of even length (+1 for odd n)

# Returns

- `w::Vector{Float64}`:: generated window
"""
function generate_window(type::Symbol, n::Int64; even::Bool=false)
    n < 1 && throw(ArgumentError("n must be ≥ 1."))
    even == true && mod(n, 2) != 0 && (n += 1)
    t = range(0, 1, n)
    if type === :hann
        w = @. 0.5 * (1 - cos.(2 * pi * t))
    elseif type === :bh
        w = @. 0.35875 - 0.48829 * cos.(2 * pi * t) + 0.14128 * cos.(4 * pi * t) - 0.01168 * cos.(6 * pi * t)
    elseif type === :bohman
        w = @. (1 - abs.(t * 2 - 1)) * cos.(pi * abs.(t * 2 - 1)) + (1 / pi) * sin.(pi * abs.(t * 2 - 1))
    elseif type === :flat
        w = @. 0.21557 - 0.41663 * cos.(2 * pi * t) + 0.27726 * cos.(4 * pi * t) - 0.08357 * cos.(6 * pi * t) + 0.00694 * cos.(8 * pi * t)
    elseif type === :bn
        w = @. 0.3635819 - 0.4891775 * cos(2 * pi * t) + 0.1365995 * cos(4 * pi * t) - 0.0106411 * cos(6 * pi * t)
    elseif type === :nutall
        w = @. 0.355768 - 0.487396 * cos(2 * pi * t) + 0.144232 * cos(4 * pi * t) - 0.012604 * cos(6 * pi * t)
    elseif type === :triangle
        mod(n, 2) == 0 && (n += 1)
        w = zeros(n)
        for idx in 1:((n ÷ 2) + 1)
            w[idx] = @. (idx * (idx + 1)) / 2
        end
        w[((n ÷ 2) + 2):n] = reverse(w)[((n ÷ 2) + 2):n]
        w .= w ./ maximum(w)
    elseif type === :exp
        mod(n, 2) == 0 && (n += 1)
        w = ones(n)
        for idx in 1:((n ÷ 2) + 1)
            w[idx] = 1 / idx
        end
        w[1:((n ÷ 2) + 1)] = reverse(w[1:((n ÷ 2) + 1)])
        w[((n ÷ 2) + 2):n] = reverse(w[1:(n ÷ 2)])
    else
        throw(ArgumentError("Window type must be :hann, :bh, :bohman, :flat, :bn, :nutall, :triangle, :exp."))
    end

    return w
end

"""
    hildebrand_rule(x)

Calculate Hildebrand rule for vector `x`.
If H < 0.2 then the vector `x` is symmetrical.

# Arguments

- `x::Vector{<:Real}`

# Returns

- `h::Float64`
"""
function hildebrand_rule(x::Vector{<:Real})

    return (mean(x) - median(x)) ./ std(x)
end

"""
    jaccard_similarity(x, y)

Calculate Jaccard similarity between two vectors `x` and `y`.

# Arguments

- `n::Int64`

# Returns

- `j::Float64`
"""
function jaccard_similarity(x::Vector{<:Real}, y::Vector{<:Real})

    i = length(intersect(x, y))
    u = length(x) + length(y) - i
    j = i / u

    return j
end

"""
    fft0(x, n)

Calculate FFT for the vector `x` padded with `n` or `n - length(x)` zeros at the end.

# Arguments

- `x::AbstractArray`
- `n::Int64`

# Returns

- `fft0::Vector{ComplexF64}`
"""
function fft0(x::AbstractArray, n::Int64)

    n < 0 && throw(ArgumentError("Pad must be positive."))

    n > length(x) && (n = n - length(x))

    return fft(vcat(x, zeros(eltype(x), n)))
end

"""
    ifft0(x, n)

Calculate IFFT for the vector `x` padded with `n` or `n - length(x)` zeros at the end.

# Arguments

- `x::AbstractArray`
- `n::Int64`

# Returns

- `ifft0::Vector{ComplexF64}`
"""
function ifft0(x::AbstractArray, n::Int64)

    n < 0 && throw(ArgumentError("Pad must be positive."))

    n > length(x) && (n = n - length(x))

    return ifft(vcat(x, zeros(eltype(x), n)))
end

"""
    nextpow2(x)

Return the next power of 2 for given number `x`.

# Argument

 - `x::Int64`

# Returns

- `nextpow::Int64`
"""
function nextpow2(x::Int64)

    if x == 0
        return 1
    else
        return 2 ^ ndigits(x - 1, base=2)
    end
end

"""
    vsplit(x, n)

Splits the vector `x` into `n`-long pieces.

# Argument

- `x::Vector{<:Real}`
 - `n::Int64`

# Returns

- `x::Vector{Vector{<:Real}}`
"""
function vsplit(x::Vector{<:Real}, n::Int64=1)

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
    s_rms(signal)

Calculate Root Mean Square of `signal`.

# Arguments

- `signal::Vector{<:Real}`

# Returns

- rms::Float64`
"""
function s_rms(signal::Vector{<:Real})

    # rms = sqrt(mean(signal.^2))    
    rms = norm(signal) / sqrt(length(signal))

    return rms
end

"""
    generate_sine(f, t, a, p)

Generates sine wave of `f` frequency over `t` time; optional arguments are: `a` amplitude and  `p` phase.

# Arguments

- `f::Real`: frequency
- `t::Union{Vector{<:Real}, AbstractRange}`: time vector
- `a::Real`: amplitude
- `p::Real`: initial phase

# Returns

- sine::Vector{Float64}`
"""
function generate_sine(f::Real, t::Union{Vector{<:Real}, AbstractRange}, a::Real=1, p::Real=0)

    return @. a * sin(2 * pi * f * t + p)
end

"""
    s_freqs(t)

Return vector of frequencies and Nyquist frequency for given time vector `t`.

# Arguments

- `t::Vector{<:Real}, AbstractRange}`

# Returns

- `hz::Vector{Float64}`
- `nyquist_freq::Float64`
"""
function s_freqs(t::Union{Vector{<:Real}, AbstractRange})

    typeof(t) <: AbstractRange && (t = collect(t))

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
    s_freqs(signal, fs)

Return vector of frequencies and Nyquist frequency for given `signal` and `fs`.

# Arguments

- `signal::Vector{Float64}`
- `fs::Real`

# Returns

- `hz::Vector{Float64`
- `nyquist_freq::Float64`
"""
function s_freqs(signal::Vector{Float64}, fs::Real)

    fs < 0 && throw(ArgumentError("Sampling rate must be >0 Hz."))

    # Nyquist frequency
    nyquist_freq = fs / 2
    # frequency array
    hz = linspace(0, nyquist_freq, floor(Int64, length(signal) / 2) + 1)

    return hz, nyquist_freq
end

"""
    m_sortperm(m; rev=false, dims=1)

Generates sorting index for matrix `m` by columns (`dims` = 1) or by rows (`dims` = 2).

# Arguments

- `m::Matrix`
- `rev::Bool`
- `dims::Int64`

# Returns

- `m_idx::Matrix{Int64}`
"""
function m_sortperm(m::Matrix; rev::Bool=false, dims::Int64=1)

    (dims < 1 || dims > 2) && throw(ArgumentError("dims must be 1 or 2."))
    
    m_idx = zeros(Int, size(m))
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
    m_sort(m, m_idx; rev=false, dims=1)

Sorts matrix `m` using sorting index `m_idx` by columns (`dims` = 1) or by rows (`dims` = 2).

# Arguments

- `m::Matrix`
- `m_idx::Vector{Int64}`
- `rev::Bool`
- `dims::Int64`

# Returns

- `m_sorted::Matrix`
"""
function m_sort(m::Matrix, m_idx::Vector{Int64}; rev::Bool=false, dims::Int64=1)

    (dims < 1 || dims > 2) && throw(ArgumentError("dims must be 1 or 2."))

    m_sorted = zeros(eltype(m), size(m))
    if dims == 1
        for idx = 1:size(m, 2)
            # sort by columns
            tmp = @view m[:, idx]
            tmp = tmp[m_idx]
            m_sorted[:, idx] = tmp
        end
    else
        for idx = 1:size(m, 1)
            # sort by rows
            tmp = @view m[idx, :]
            tmp = tmp[m_idx]
            m_sorted[idx, :] = tmp
        end
    end
    return m_sorted
end

"""
    pad0(x, n, sym)

Pad the vector `x` with `n` zeros.

# Arguments

- `x::Vector{<:Real}`
- `n::Int64`
- `sym::Bool=false`: if true, than pad at the beginning and at the end, otherwise only at the end.

# Returns

- `v_pad::Vector{<:Real}`
"""
function pad0(x::AbstractArray, n::Int64, sym::Bool=false)

    n < 0 && throw(ArgumentError("n must be ≥ 0."))

    sym == false && (v_pad = vcat(x, zeros(eltype(x), n)))
    sym == true && (v_pad = vcat(zeros(eltype(x), n), x, zeros(eltype(x), n)))

    return v_pad
end
"""
    hz2rads(f)

Convert frequency `f` in Hz to rad/s.

# Arguments

- `f::Real`

# Returns

- `f_rads::Float64`
"""
function hz2rads(f::Real)

    return 2 * pi * f
end

"""
    rads2hz(f)

Convert frequency `f` in rad/s to Hz.

# Arguments

- `f::Real`

# Returns

- `f_rads::Float64`
"""
function rads2hz(f::Real)

    return f / 2 * pi
end

"""
    z_score(x)

Calculate Z-scores for each value of the vector `x`.

# Arguments

- `x::Vector{<:Real}`

# Returns

- `z_score::Vector{Float64}`
"""
function z_score(x::Vector{<:Real})

    return (x .- mean(x)) ./ std(x)
end

"""
    k_categories(n)

Calculate number of categories for a given sample size `n`.

# Arguments

- `n::Int64`

# Returns

- `k::Float64`
"""
function k_categories(n::Int64)

    return (sqrt(n), (1 + 3.222 * log10(n)))
end

"""
    cmax(x)

Return maximum value of the complex vector`x`.

# Arguments

- `x::Vector{ComplexF64}`

# Returns

- `cmax::ComplexF64`
"""
function cmax(x::Vector{ComplexF64})

    return argmax(abs, x)
end

"""
    cmin(x)

Return minimum value of the complex vector`x`.

# Arguments

- `x::Vector{ComplexF64}`

# Returns

- `cmin::ComplexF64`
"""
function cmin(x::Vector{ComplexF64})

    return argmin(abs, x)
end

"""
    generate_sinc(t, f, peak, norm)

Generate normalized or unnormalized sinc function.

# Arguments

- `t::AbstractRange=-2:0.01:2`: time
- `f::Real=10.0`: frequency
- `peak::Real=0`: sinc peak time
- `norm::Bool=true`: generate normalized function
# Returns

- `sinc::Vector{Float64}`
"""
function generate_sinc(t::AbstractRange=-2:0.01:2; f::Real=1, peak::Real=0, norm::Bool=true)

    norm == true && (y_sinc = @. sin(2 * pi * f * (t - peak)) / (pi * (t - peak)))
    norm == false && (y_sinc = @. sin(2 * f * (t - peak)) / (t - peak))
    nan_idx = isnan.(y_sinc)
    sum(nan_idx) != 0 && (y_sinc[findall(isnan, y_sinc)[1]] = (y_sinc[findall(isnan, y_sinc)[1] - 1] + y_sinc[findall(isnan, y_sinc)[1] + 1]) / 2)
    
    return y_sinc
end

"""
    generate_morlet(fs, f, t; ncyc, complex)

Generate Morlet wavelet.

# Arguments

- `fs::Int64`: sampling rate
- `f::Real`: frequency
- `t::Real=1`: length = -t:1/fs:t
- `ncyc::Int64=5`: number of cycles
- `complex::Bool=false`: generate complex Morlet

# Returns

- `morlet::Union{Vector{Float64}, Vector{ComplexF64}}`
"""
function generate_morlet(fs::Int64, f::Real, t::Real=1; ncyc::Int64=5, complex::Bool=false)

    t = -t:1/fs:t
    complex == true && (sin_wave = @. exp(im * 2 * pi * f * t))       # for symmetry at x = 0
    complex == false && (sin_wave = @. sin(2 * pi * f * t))           # for symmetry at x = 0
    g = generate_gaussian(fs, f, t[end], ncyc=ncyc)
    m = sin_wave .* g
    
    return m
end

"""
    generate_gaussian(fs, f, t; ncyc, a)

Generate Gaussian wave.

# Arguments

- `fs::Int64`: sampling rate
- `f::Real`: frequency
- `t::Real=1`: length = -t:1/fs:t
- `ncyc::Int64`: : number of cycles, width, SD of the Gaussian
- `a::Real=1`: peak amp

# Returns

- `gaussian::Vector{Float64}`
"""
function generate_gaussian(fs::Int64, f::Real, t::Real=1; ncyc::Int64=5, a::Real=1.0)

    t = -t:1/fs:t
    s = ncyc / (2 * pi * f)             # Gaussian width (standard deviation)
    g = @. a * exp(-(t/s)^2 / 2)        # Gaussian

    return g
end

"""
    tuple_order(t, rev)

Order tuple elements in ascending or descending (rev=true) order.

# Arguments

- `t::Tuple{Real, Real}`
- `rev::Bool=false`

# Returns

- `t::Tuple{Real, Real}`
"""
function tuple_order(t::Tuple{Real, Real}, rev::Bool=false)
    
    (rev == false && t[1] > t[2]) && (t = (t[2], t[1]))
    (rev == true && t[1] < t[2]) && (t = (t[2], t[1]))

    return t
end

"""
    s2_rmse(signal1, signal2)

Calculate RMSE between `signal1` and `signal2`.

# Arguments

- `signal1::Vector{Float64}`
- `signal2::Vector{Float64}`

# Returns

- `r::Float64`
"""
function s2_rmse(signal1::Vector{Float64}, signal2::Vector{Float64})

    # r = sum(signal1 .* signal2) ./ (sqrt(sum(signal1.^2)) .* sqrt(sum(signal2.^2)))
    r = sqrt(mean(signal2 - signal1)^2)

    return r
end

"""
    m_norm(m)

Normalize matrix `m`.

# Arguments

- `m::Matrix{Float64}`

# Returns

- `m_norm::Matrix{Float64}`
"""
function m_norm(m::Array{Float64, 3})
    m_norm = m ./ (size(m, 2) - 1)
    
    return m_norm
end

"""
   s_cov(signal; norm=true)

Calculates covariance between all channels of the `signal`.

# Arguments

- `signal::AbstractArray`
- `norm::Bool`: normalize covariance

# Returns

- `cov_mat::Matrix{Float64}`
"""
function s_cov(signal::AbstractArray; norm::Bool=false)

    signal = signal'
    
    # channels-vs-channels
    cov_mat = cov(signal)

    # normalize
    norm == true && (cov_mat = m_norm(cov_mat))

    return cov_mat
end

"""
   s2_cov(signal1, signal2; norm=true)

Calculates covariance between `signal1` and `signal2`.

# Arguments

- `signal1::AbstractArray`
- `signal2::AbstractArray`
- `norm::Bool`: normalize covariance

# Returns

- `cov_mat::Matrix{Float64}`
"""
function s2_cov(signal1::AbstractArray, signal2::AbstractArray; norm::Bool=false)

    length(signal1) == length(signal2) || throw(ArgumentError("Both signals must be of the same length."))

    # channels-vs-channels
    cov_mat = cov(signal1 * signal2')

    # normalize
    norm == true && (cov_mat = cov_mat ./ (size(cov_mat, 2) - 1))

    return cov_mat
end

"""
    s_dft(signal; fs)

Return FFT and DFT sample frequencies for a DFT for the `signal`.

# Arguments

- `signal::AbstractArray`
- `fs::Int64`: sampling rate

# Returns

Named tuple containing:
- `s_fft::Vector{ComplexF64}`
- `s_sf::Vector{Float64}`
"""
function s_dft(signal::AbstractArray; fs::Int64)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))

    s_fft = fft(signal)
    # number of samples
    n = length(signal)
    # time between samples
    d = 1 / fs
    s_sf = Vector(fftfreq(n, d))

    return (s_fft=s_fft, s_sf=s_sf)
end

"""
    s_msci95(signal)

Calculate mean, std and 95% confidence interval for `signal`.

# Arguments

- `signal::Vector{Float64}`

# Returns

- `s_m::Float64`: mean
- `s_s::Float64`: standard deviation
- `s_u::Float64`: upper 95% CI
- `s_l::Float64`: lower 95% CI
"""
function s_msci95(signal::Vector{Float64})

    s_m = mean(signal)
    s_s = std(signal) / sqrt(length(signal))
    s_u = s_m + 1.96 * s_s
    s_l = s_m - 1.96 * s_s

    return s_m, s_s, s_u, s_l
end

"""
    s_msci95(signal; n=3, method=:normal)

Calculates mean, std and 95% confidence interval for each the `signal` channels.

# Arguments

- `signal::AbstractArray`
- `n::Int64`: number of bootstraps
- `method::Symbol[:normal, :boot]`: use normal method or `n`-times boostrapping

# Returns

- `s_m::Vector{Float64}`: mean
- `s_s::Vector{Float64}`: standard deviation
- `s_u::Vector{Float64}`: upper 95% CI
- `s_l::Vector{Float64}`: lower 95% CI
"""
function s_msci95(signal::AbstractArray; n::Int64=3, method::Symbol=:normal)

    method in [:normal, :boot] || throw(ArgumentError("method must be :normal or :boot."))

    if method === :normal
        s_m = mean(signal, dims=1)'
        s_s = std(signal, dims=1)' / sqrt(size(signal, 1))
        s_u = s_m + 1.96 * s_s
        s_l = s_m - 1.96 * s_s
    else
        s_tmp1 = zeros(size(signal, 1) * n, size(signal, 2))
        Threads.@threads for idx1 in 1:size(signal, 1) * n
            s_tmp2 = zeros(size(signal))
            sample_idx = rand(1:size(signal, 1), size(signal, 1))
            for idx2 in 1:size(signal, 1)
                s_tmp2[idx2, :] = signal[sample_idx[idx2], :]'
            end
            s_tmp1[idx1, :] = mean(s_tmp2, dims=1)
        end

        s_m = mean(s_tmp1, dims=1)'
        s_s = std(s_tmp1, dims=1)' / sqrt(size(s_tmp1, 1))
        s_sorted = sort(s_tmp1, dims=1)
        s_l = s_sorted[round(Int, 0.025 * size(s_tmp1, 1)), :]
        s_u = s_sorted[round(Int, 0.975 * size(s_tmp1, 1)), :]
    end

    return vec(s_m[:, 1]), vec(s_s[:, 1]), vec(s_u[:, 1]), vec(s_l[:, 1])
end

"""
    s2_mean(signal1, signal2)

Calculates mean and 95% confidence interval for 2 signals.

# Arguments

- `signal1::Vector{Float64}`
- `signal2:Vector{Float64}`

# Returns

- `s_m::Float64`: mean
- `s_s::Float64`: standard deviation
- `s_u::Float64`: upper 95% CI
- `s_l::Float64`: lower 95% CI
"""
function s2_mean(signal1::Vector{Float64}, signal2::Vector{Float64})

    length(signal1) != length(signal2) && throw(ArgumentError("Both signals must be of the same as size."))

    s_m = zeros(length(signal1))
    s_s = zeros(length(signal1))
    s_u = zeros(length(signal1))
    s_l = zeros(length(signal1))

    s1_mean = mean(signal1)
    s2_mean = mean(signal2)
    s_m = s1_mean - s2_mean
    s1_sd = std(signal1) / sqrt(length(signal1))
    s2_sd = std(signal2) / sqrt(length(signal2))
    s_s = sqrt(s1_sd^2 + s2_sd^2)
    s_u = s_m + 1.96 * s_s
    s_l = s_m - 1.96 * s_s

    return s_m, s_s, s_u, s_l
end

"""
    s2_difference(signal1, signal2; n=3, method=:absdiff)

Calculates mean difference and 95% confidence interval for 2 signals.

# Arguments

- `signal1::AbstractArray`
- `signal2::AbstractArray`
- `n::Int64`: number of bootstraps
- `method::Symbol[:absdiff, :diff2int]`
    - `:absdiff`: maximum difference
    - `:diff2int`: integrated area of the squared difference

# Returns

- `s_statistic::Vector{Float64}`
- `s_statistic_single::Float64`
- `p::Float64`
"""
function s2_difference(signal1::AbstractArray, signal2::AbstractArray; n::Int64=3, method::Symbol=:absdiff)

    size(signal1) != size(signal2) && throw(ArgumentError("Both signals must be of the same size."))
    method in [:absdiff, :diff2int] || throw(ArgumentError("method must be :absdiff or :diff2int."))

    s1_mean = mean(signal1, dims=1)'
    s2_mean = mean(signal2, dims=1)'

    if method === :absdiff
        # statistic: maximum difference
        s_diff = s1_mean - s2_mean
        s_statistic_single = maximum(abs.(s_diff))
    else
        # statistic: integrated area of the squared difference
        s_diff_squared = (s1_mean - s2_mean).^2
        s_statistic_single = simpson(s_diff_squared)
    end

    signals = [signal1; signal2]
    s_statistic = zeros(size(signal1, 1) * n)

    Threads.@threads for idx1 in 1:(size(signal1, 1) * n)
        s_tmp1 = zeros(size(signal1, 1), size(signal1, 2))
        sample_idx = rand(1:size(signals, 1), size(signals, 1))
        # sample_idx = sample_idx[1:1000]
        for idx2 in 1:size(signal1, 1)
            s = @view signals[sample_idx[idx2], :]
            s_tmp1[idx2, :] = s'
        end
        s1_mean = mean(s_tmp1, dims=1)
        s_tmp1 = zeros(size(signal1, 1), size(signal1, 2))
        sample_idx = rand(1:size(signals, 1), size(signals, 1))
        for idx2 in 1:size(signal1, 1)
            s = @view signals[sample_idx[idx2], :]
            s_tmp1[idx2, :] = s'
        end
        s2_mean = mean(s_tmp1, dims=1)
        if method === :absdiff
            # statistic: maximum difference
            s_diff = s1_mean - s2_mean
            s_statistic[idx1] = maximum(abs.(s_diff))
        else
            # statistic: integrated area of the squared difference
            s_diff_squared = (s1_mean - s2_mean).^2
            s_statistic[idx1] = simpson(s_diff_squared)
        end
    end

    p = length(s_statistic[s_statistic .> s_statistic_single]) / size(signal1, 1) * n
    p > 1 && (p = 1.0)

    return s_statistic, s_statistic_single, p
end

"""
   s_acov(signal; lag=1, demean=false, norm=false)

Calculate autocovariance of the `signal`.

# Arguments

- `signal::AbstractArray`
- `lag::Int64`: lags range is `-lag:lag`
- `demean::Bool`: demean `signal` prior to calculations
- `norm::Bool`: normalize autocovariance

# Returns

- `acov::Vector{Float64}`
- `lags::Vector{Int64}`
"""
function s_acov(signal::AbstractArray; lag::Int64=1, demean::Bool=false, norm::Bool=false)

    lag < 1 && throw(ArgumentError("lag must be ≥ 1."))

    lags = collect(-lag:lag)

    if demean == true
        s_demeaned = s_demean(signal)
    else
        s_demeaned = signal
    end

    acov = zeros(length(lags))

    for idx in 1:length(lags)
        if lags[idx] == 0
            # no lag
            s_lagged = s_demeaned
            s_mul = s_demeaned .* s_lagged
        elseif lags[idx] > 0
            # positive lag
            s1 = @view s_demeaned[(1 + lags[idx]):end]
            s2 = @view s_demeaned[1:(end - lags[idx])]
            s_mul =  s1 .* s2
        elseif lags[idx] < 0
            # negative lag
            s1 = @view s_demeaned[1:(end - abs(lags[idx]))]
            s2 = @view s_demeaned[(1 + abs(lags[idx])):end]
            s_mul = s1 .* s2
        end
        s_sum = sum(s_mul)
        if norm == true
            acov[idx] = s_sum / length(signal)
        else
            acov[idx] = s_sum
        end
    end

    return acov, lags
end

"""
   s_xcov(signal1, signal2; lag=1, demean=false, norm=false)

Calculates cross-covariance between `signal1` and `signal2`.

# Arguments

- `signal1::AbstractArray`
- `signal2::AbstractArray`
- `lag::Int64`: lags range is `-lag:lag`
- `demean::Bool`: demean signal prior to analysis
- `norm::Bool`: normalize cross-covariance

# Returns

- `ccov::Vector{Float64}`
- `lags::Vector{Int64}`
"""
function s_xcov(signal1::AbstractArray, signal2::AbstractArray; lag::Int64=1, demean::Bool=false, norm::Bool=false)

    length(signal1) != length(signal2) && throw(ArgumentError("Both signals must be of the same as length."))
    lag < 1 && throw(ArgumentError("lag must be ≥ 1."))

    lags = collect(-lag:lag)

    if demean == true
        s_demeaned1 = s_demean(signal1)
        s_demeaned2 = s_demean(signal2)
    else
        s_demeaned1 = signal1
        s_demeaned2 = signal2
    end

    xcov = zeros(length(lags))

    for idx in 1:length(lags)
        if lags[idx] == 0
            # no lag
            s_lagged = s_demeaned2
            s_mul = s_demeaned1 .* s_lagged
        elseif lags[idx] > 0
            # positive lag
            s1 = @view s_demeaned1[(1 + lags[idx]):end]
            s2 = @view s_demeaned2[1:(end - lags[idx])]
            s_mul = s1 .* s2
        elseif lags[idx] < 0
            # negative lag
            s1 = @view s_demeaned1[1:(end - abs(lags[idx]))] 
            s2 = @view s_demeaned2[(1 + abs(lags[idx])):end]
            s_mul = s1 .* s2
        end
        s_sum = sum(s_mul)
        if norm == true
            xcov[idx] = s_sum / length(signal1)
        else
            xcov[idx] = s_sum
        end
    end

    return xcov, lags
end

"""
    s_spectrum(signal; pad)

Calculates FFT, amplitudes, powers and phases of the `signal`.

# Arguments

- `signal::AbstractArray`
- `pad::Int64=0`: pad the `signal` with `pad` zeros

# Returns

Named tuple containing:
- `s_fft::Vector{ComplexF64}`
- `s_amplitudes::Vector{Float64}`
- `s_powers::Vector{Float64}`
- `s_phases::Vector{Float64}`
"""
function s_spectrum(signal::AbstractArray; pad::Int64=0)

    pad < 0 && throw(ArgumentError("pad must be ≥ 0."))
    s_fft = fft0(signal, pad)

    # normalize
    s_fft ./= length(signal)
    # amplitudes
    s_amplitudes = @. 2 * abs(s_fft)
    # power
    s_powers = s_amplitudes.^2
    # phases
    s_phases = angle.(s_fft)

    return (s_fft=s_fft, s_amplitudes=s_amplitudes, s_powers=s_powers, s_phases=s_phases)
end

"""
    s_total_power(signal; fs)

Calculates `signal` total power.

# Arguments

- `signal::AbstractArray`
- `fs::Int64`: sampling rate
- `mt::Bool=false`: if true use multi-tapered periodogram

# Returns

- `stp::Float64`: signal total power
"""
function s_total_power(signal::AbstractArray; fs::Int64, mt::Bool=false)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))
    mt == false && (psd = welch_pgram(signal, 4*fs, fs=fs))
    mt == true && (psd = mt_pgram(signal, fs=fs))
    # dx: frequency resolution
    dx = psd.freq[2] - psd.freq[1]
    stp = simpson(psd.power, dx=dx)

    return stp
end

"""
    s_band_power(signal; fs, f)

Calculates `signal` power between `f[1]` and `f[2]`.

# Arguments

- `signal::AbstractArray`
- `fs::Int64`: sampling rate
- `f::Tuple{Real, Real}`: lower and upper frequency bounds

# Returns

- `sbp::Float64`: signal band power
"""
function s_band_power(signal::AbstractArray; fs::Int64, f::Tuple{Real, Real}, mt::Bool=false)

    mt == false && (psd = welch_pgram(signal, 4*fs, fs=fs))
    mt == true && (psd = mt_pgram(signal, fs=fs))

    psd_freq = Vector(psd.freq)
    
    f1_idx = vsearch(f[1], psd_freq)
    f2_idx = vsearch(f[2], psd_freq)
    frq_idx = [f1_idx, f2_idx]

    # dx: frequency resolution
    dx = psd_freq[2] - psd_freq[1]
    sbp = simpson(psd.power[frq_idx[1]:frq_idx[2]], psd_freq[frq_idx[1]:frq_idx[2]], dx=dx)

    return sbp
end

"""
    s_taper(signal; taper)

Taper the `signal` with `taper`.

# Arguments

- `signal::AbstractArray`
- `taper::Union{Vector{<:Real}, Vector{ComplexF64}}`

# Returns

- `s_tapered::Vector{Union{Float64, ComplexF64}}`
"""
function s_taper(signal::AbstractArray; taper::Union{Vector{<:Real}, Vector{ComplexF64}})

    length(taper) == length(signal) || throw(ArgumentError("Taper and signal lengths must be equal."))
    s_tapered = signal .* taper

    return s_tapered
end

"""
    s_detrend(signal; type, offset, order, span)
Perform piecewise detrending of `eeg`.

# Arguments

- `eeg::NeuroJ.EEG`
- `type::Symbol`, optional
    - `:ls`: the result of a linear least-squares fit to `signal` is subtracted from `signal`
    - `:linear`: linear trend is subtracted from `signal`
    - `:constant`: `offset` or the mean of `signal` (if `offset` = 0) is subtracted
    - `:poly`: polynomial of `order` is subtracted
    - `:loess`: fit and subtract loess approximation
- `offset::Real=0`: constant for :constant detrending
- `order::Int64=1`: polynomial fitting order
- `span::Float64=0.5`: smoothing of loess

# Returns
- `s_det::Vector{Float64}`
"""
function s_detrend(signal::AbstractArray; type::Symbol=:linear, offset::Real=0, order::Int64=1, span::Float64=0.5)

    type in [:ls, :linear, :constant, :poly, :loess] || throw(ArgumentError("type must be :ls, :linear, :constant, :poly, :loess."))

    if type === :loess
        t = collect(1.0:1:length(signal))
        model = loess(t, signal, span=span)
        trend = Loess.predict(model, t)
        s_det = signal .- trend

        return s_det
    end

    if type === :poly
        t = collect(1:1:length(signal))        
        p = Polynomials.fit(t, signal, order)
        trend = zeros(length(signal))
        for idx in 1:length(signal)
            trend[idx] = p(t[idx])
        end
        s_det = signal .- trend

        return s_det
    end

    if type === :constant
        offset == 0 && (offset = mean(signal))
        s_det = signal .- mean(signal)

        return s_det
    end

    if type === :ls
        T = eltype(signal)
        N = size(signal, 1)
        # create linear trend matrix
        A = similar(signal, T, N, 2)
        A[:,2] .= T(1)
        A[:,1] .= range(T(0),T(1),length=N)
        # create linear trend matrix
        R = transpose(A) * A
        # do the matrix inverse for 2×2 matrix
        Rinv = inv(Array(R)) |> typeof(R)
        factor = Rinv * transpose(A)
        s_det = signal .- A * (factor * signal)

        return s_det
    end

    if type == :linear
        trend = linspace(signal[1], signal[end], length(signal))
        s_det = signal .- trend

        return s_det
    end
end

"""
    s_demean(signal)

Remove mean value (DC offset) from the `signal`.

# Arguments

- `signal::AbstractArray`

# Returns

- `s_demeaned::Vector{Float64}`
"""
function s_demean(signal::AbstractArray)

    m = mean(signal)
    s_demeaned = signal .- m

    return s_demeaned
end

"""
    s_normalize_zscore(signal)

Normalize (by z-score) `signal`.

# Arguments

- `signal::AbstractArray`

# Returns

- `s_normalized::Vector{Float64}`
"""
function s_normalize_zscore(signal::AbstractArray)

    m = mean(signal)
    s = std(signal)
    s_normalized = (signal .- m) ./ s

    return s_normalized
end

"""
    s_normalize_minmax(signal)

Normalize `signal` in [-1, +1].

# Arguments

- `signal::AbstractArray`

# Returns

- `s_normalized::Vector{Float64}`
"""
function s_normalize_minmax(signal::AbstractArray)

    mi = minimum(signal)
    mx = maximum(signal)
    s_normalized = 2 .* (signal .- mi) ./ (mx - mi) .- 1

    return s_normalized
end

"""
    s_add_noise(signal)

Adds random noise to the `signal`.

# Arguments

- `signal::AbstractArray`

# Returns

- `s_noisy::Vector{Float64}`
"""
function s_add_noise(signal::AbstractArray)

    s_noise = signal .+ rand(length(signal))

    return s_noise
end

"""
    s_resample(signal; t, new_sr)

Resample `signal` to `new_sr` sampling frequency.

# Arguments

- `signal::AbstractArray`
- `t::AbstractRange`: time
- `new_sr::Int64`: new sampling rate

# Returns

- `s_resampled::Vector{Float64}`
- `t_resampled::AbstractRange`
"""
function s_resample(signal::AbstractArray; t::AbstractRange, new_sr::Int64)

    new_sr < 1 && throw(ArgumentError("New sampling rate must be positive."))

    # sampling interval
    dt = t[2] - t[1]
    # sampling rate
    sr = 1 / dt
    new_sr == sr && return(signal)

    # interpolate
    sr_ratio = new_sr / sr
    s_upsampled = resample(signal, sr_ratio)
    # s_interpolation = CubicSplineInterpolation(t, signal)
    t_upsampled = t[1]:1/new_sr:t[end]
    # s_upsampled = s_interpolation(t_upsampled)

    return s_upsampled, t_upsampled
end

"""
    s_resample(signal; t, new_sr)

Resamples all channels of the`signal` and time vector `t` to `new_sr` sampling frequency.

# Arguments

- `signal::Array{Float64, 3}`
- `t::AbstractRange`
- `new_sr::Int64`: new sampling rate

# Returns

- `s_downsampled::Array{Float64, 3}`
- `t_downsampled::AbstractRange`
"""
function s_resample(signal::Array{Float64, 3}; t::AbstractRange, new_sr::Int64)

    new_sr < 1 && throw(ArgumentError("New sampling rate must be positive."))

    channel_n, _, epoch_n = size(signal)

    s_resampled_len = length(s_resample(signal[1, :, 1], t=t, new_sr=new_sr)[1])
    s_resampled = zeros(channel_n, s_resampled_len, epoch_n) 

    t_resampled = nothing
    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            s = @view signal[channel_idx, :, epoch_idx]
            s_resampled[channel_idx, :, epoch_idx], t_resampled = s_resample(s, t=t, new_sr=new_sr)
        end
    end

    return s_resampled, t_resampled
end

"""
    s_derivative(signal)

Return derivative of `signal` of the same length.
"""
function s_derivative(signal::AbstractArray)

    s_der = diff(signal)
    s_der = vcat(s_der, s_der[end])
    
    return s_der
end

"""
    s_tconv(signal; kernel)

Performs convolution in the time domain between `signal` and `kernel`.

# Arguments

- `signal::AbstractArray`
- `kernel::Union{Vector{<:Real}, Vector{ComplexF64}}`

# Returns

- `s_conv::Union{Vector{Float64}, Vector{ComplexF64}}`
"""
function s_tconv(signal::AbstractArray; kernel::Union{Vector{<:Real}, Vector{ComplexF64}})

    signal = Vector(signal)
    s_conv = conv(signal, kernel)
    half_kernel = floor(Int, length(kernel) / 2)

    # remove in- and out- edges
    if mod(length(kernel), 2) == 0 
        s_conv = s_conv[half_kernel:(end - half_kernel)]
    else
        s_conv = s_conv[(half_kernel + 1):(end - half_kernel)]
    end

    return s_conv
end

"""
    s_filter(signal; <keyword arguments>)

Filter `signal`.

# Arguments

- `signal::AbstractArray`
- `fprototype::Symbol`: filter prototype:
    - `:butterworth`
    - `:chebyshev1`
    - `:chebyshev2`
    - `:elliptic`
    - `:fir`
    - `:iirnotch`
    - `:remez`
    - `:mavg`: moving average (with threshold)
    - `:mmed`: moving median (with threshold)
    - `:poly`: polynomial of `order` order
- `ftype::Symbol`: filter type:
    - `:lp`: low pass
    - `:hp`: high pass
    - `:bp`: band pass
    - `:bs`: band stop
- `cutoff::Union{Real, Tuple}`: filter cutoff in Hz (vector for `:bp` and `:bs`)
- `order::Int64=8`: filter order, number of taps for :remez filter, k-value for :mavg and :mmed (window length = 2 × k + 1)
- `rp::Real=-1`: ripple amplitude in dB in the pass band; default: 0.0025 dB for :elliptic, 2 dB for others
- `rs::Real=-1`: ripple amplitude in dB in the stop band; default: 40 dB for :elliptic, 20 dB for others
- `bw::Real=-1`: bandwidth for :iirnotch and :remez filters
- `dir:Symbol=:twopass`: filter direction (:onepass, :onepass_reverse, :twopass), for causal filter use :onepass
- `t::Real`: threshold for :mavg and :mmed filters; threshold = threshold * std(signal) + mean(signal) for :mavg or threshold = threshold * std(signal) + median(signal) for :mmed filter
- `window::Union{Vector{<:Real}, Nothing} - window, required for FIR filter, weighting window for :mavg and :mmed 

# Returns

- `s_filtered::Vector{Float64}`
"""
function s_filter(signal::AbstractArray; fprototype::Symbol, ftype::Union{Symbol, Nothing}=nothing, cutoff::Union{Real, Tuple}=0, fs::Int64=0, order::Int64=8, rp::Real=-1, rs::Real=-1, bw::Real=-1, dir::Symbol=:twopass, t::Real=0, window::Union{Vector{<:Real}, Nothing}=nothing)

    fprototype in [:mavg, :mmed, :poly, :butterworth, :chebyshev1, :chebyshev2, :elliptic, :fir, :iirnotch, :remez] || throw(ArgumentError("fprototype must be :mavg, :mmed, :poly, :butterworth, :chebyshev1, :chebyshev2, :elliptic, :fir, :iirnotch or :remez."))

    if fprototype === :fir
        if window === nothing
            @warn "Using default window for :fir filter: hanning($(3 * floor(Int64, fs / cutoff[1])))."
            window = hanning(3 * floor(Int64, fs / cutoff[1]))
        end
        if ftype === :hp || ftype === :bp || ftype === :bs
            mod(length(window), 2) == 0 && (window = vcat(window[1:((length(window) ÷ 2) - 1)], window[((length(window) ÷ 2) + 1):end]))
        end
        fprototype = FIRWindow(window)
    end

    order < 1 && throw(ArgumentError("order must be > 1."))
    window !== nothing && length(window) > length(signal) && throw(ArgumentError("For :fir filter window must be shorter than signal."))
    (fprototype !== :mavg && fprototype !== :mmed && fprototype !== :poly && fprototype !== :conv && fprototype !== :iirnotch && fprototype !== :remez) && (ftype in [:lp, :hp, :bp, :bs] || throw(ArgumentError("ftype must be :bp, :hp, :bp or :bs.")))
    (fprototype !== :mavg && fprototype !== :conv && fprototype !== :mmed) && (fs < 1 && throw(ArgumentError("fs must be > 0.")))
    dir in [:onepass, :onepass_reverse, :twopass] || throw(ArgumentError("direction must be :onepass, :onepass_reverse or :twopass."))
    ((order < 2 && fprototype !== :poly && fprototype !== :remez) && mod(order, 2) != 0) && throw(ArgumentError("order must be even and ≥ 2."))
    (window !== nothing && length(window) != (2 * order + 1) && (fprototype === :mavg || fprototype === :mmed)) && throw(ArgumentError("For :mavg and :mmed window length must be 2 × order + 1 ($(2 * order + 1))."))
    order > length(signal) && throw(ArgumentError("order must be ≤ signal length ($length(signal))."))
    
    if rp == -1
        if fprototype === :elliptic
            rp = 0.0025
        else
            rp = 2
        end
    end
    
    if rs == -1
        if fprototype === :elliptic
            rp = 40
        else
            rp = 20
        end
    end

    if fprototype === :mavg
        signal = _reflect(signal)
        s_filtered = zeros(length(signal))
        window === nothing && (window = ones(2 * order + 1))
        for idx in (1 + order):(length(signal) - order)
            if t > 0
                if signal[idx] > t * std(signal) + mean(signal)
                    s_filtered[idx] = mean(signal[(idx - order):(idx + order)] .* window)
                end
            else
                s_filtered[idx] = mean(signal[(idx - order):(idx + order)] .* window)
            end
        end
        s_filtered = _chop(s_filtered)

        return s_filtered
    end

    if fprototype === :mmed
        signal = _reflect(signal)
        s_filtered = zeros(length(signal))
        window === nothing && (window = ones(2 * order + 1))
        for idx in (1 + order):(length(signal) - order)
            if t > 0
                if signal[idx] > t * std(signal) + median(signal)
                    s_filtered[idx] = median(signal[(idx - order):(idx + order)] .* window)
                end
            else
                s_filtered[idx] = median(signal[(idx - order):(idx + order)] .* window)
            end
        end
        s_filtered = _chop(s_filtered)

        return s_filtered
    end

    if fprototype === :poly
        t = collect(0:1/fs:(length(signal) - 1) / fs)        
        p = Polynomials.fit(t, signal, order)
        s_filtered = zeros(length(signal))
        for idx in 1:length(signal)
            s_filtered[idx] = p(t[idx])
        end

        return s_filtered
    end

    if fprototype === :conv
        s_filtered = s_tconv(signal, window)

        return s_filtered
    end

    if ftype === :lp
        length(cutoff) != 1 && throw(ArgumentError("For :lp filter one frequency must be given."))
        responsetype = Lowpass(cutoff; fs=fs)
    elseif ftype === :hp
        length(cutoff) != 1 && throw(ArgumentError("For :hp filter one frequency must be given."))
        responsetype = Highpass(cutoff; fs=fs)
    elseif ftype === :bp
        length(cutoff) != 2 && throw(ArgumentError("For :bp filter two frequencies must be given."))
        responsetype = Bandpass(cutoff[1], cutoff[2]; fs=fs)
    elseif ftype === :bs
        length(cutoff) != 2 && throw(ArgumentError("For :bs filter two frequencies must be given."))
        responsetype = Bandstop(cutoff[1], cutoff[2]; fs=fs)
    end

    fprototype === :butterworth && (prototype = Butterworth(order))
    if fprototype === :fir
        prototype = FIRWindow(window)
    end
    if fprototype === :chebyshev1
        (rs < 0 || rs > fs / 2) && throw(ArgumentError("For :chebyshev1 filter rs must be > 0 and ≤ $(fs / 2)."))
        prototype = Chebyshev1(order, rs)
    end
    if fprototype === :chebyshev2
        (rp < 0 || rp > fs / 2) && throw(ArgumentError("For :chebyshev2 filter rp must be > 0 and ≤ $(fs / 2)."))
        prototype = Chebyshev2(order, rp)
    end
    if fprototype === :elliptic
        (rs < 0 || rs > fs / 2) && throw(ArgumentError("For :elliptic filter rs must be > 0 and ≤ $(fs / 2)."))
        (rp < 0 || rp > fs / 2) && throw(ArgumentError("For :elliptic filter rp must be > 0 and ≤ $(fs / 2)."))
        prototype = Elliptic(order, rp, rs)
    end

    if fprototype === :iirnotch
        flt = iirnotch(cutoff, bw, fs=fs)
    elseif fprototype === :remez
        ftype === :lp && (window = [(0, cutoff - bw) => 1, (cutoff + bw, fs / 2) => 0])
        ftype === :hp && (window = [(0, cutoff - bw) => 0, (cutoff + bw, fs / 2) => 1])
        ftype === :bp && (window = [(0, cutoff[1] - bw / 2) => 0, (cutoff[1] + bw / 2, cutoff[2] - bw / 2) => 1, (cutoff[2] + bw / 2, fs / 2) => 0])
        ftype === :bs && (window = [(0, cutoff[1] - bw / 2) => 1, (cutoff[1] + bw / 2, cutoff[2] - bw / 2) => 0, (cutoff[2] + bw / 2, fs / 2) => 1])
        flt = remez(order, window, Hz=fs)
    else
        flt = digitalfilter(responsetype, prototype)
    end

    if fprototype !== :mavg && fprototype !== :mmed && fprototype !== :conv
        dir === :twopass && (s_filtered = filtfilt(flt, signal))
        dir === :onepass && (s_filtered = filt(flt, signal))
        dir === :onepass_reverse && (s_filtered = filt(flt, reverse(signal)))
    end

    return s_filtered
end

"""
    s_psd(signal; fs, norm=false)

Calculate power spectrum density of the `signal`.

# Arguments
- `signal::AbstractArray`
- `fs::Int64`: sampling rate
- `norm::Bool`: normalize do dB
- `mt::Bool=false`: if true use multi-tapered periodogram

# Returns

Named tuple containing:
- `psd_pow::Vector{Float64}`
- `psd_frq::Vector{Float64}`
"""
function s_psd(signal::AbstractArray; fs::Int64, norm::Bool=false, mt::Bool=false)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))
    
    mt == false && (psd = welch_pgram(signal, 4*fs, fs=fs))
    mt == true && (psd = mt_pgram(signal, fs=fs))
    psd_pow = power(psd)
    psd_frq = freq(psd)
    norm == true && (psd_pow = pow2db.(psd_pow))

    return (psd_pow=psd_pow, psd_frq=Vector(psd_frq))
end

"""
    s_psd(signal; fs, norm=false)

Calculate power spectrum density of the `signal`.

# Arguments
- `signal::Matrix{Float64}`
- `fs::Int64`: sampling rate
- `norm::Bool`: normalize do dB
- `mt::Bool=false`: if true use multi-tapered periodogram

# Returns

named tuple containing:
- `psd_pow::Matrix{Float64}`
- `psd_frq::Matrix{Float64}`
"""
function s_psd(signal::Matrix{Float64}; fs::Int64, norm::Bool=false, mt::Bool=false)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))

    channel_n = size(signal, 1)
    psd_tmp, frq_tmp = s_psd(signal[1, :], fs=fs, norm=norm, mt=mt)
    psd_pow = zeros(channel_n, length(psd_tmp))
    psd_frq = zeros(channel_n, length(frq_tmp))

    @inbounds @simd for channel_idx in 1:channel_n
        s = @view signal[channel_idx, :]
        psd_pow[channel_idx, :], psd_frq[channel_idx, :] = s_psd(s, fs=fs, norm=norm, mt=mt)
    end
    
    return (psd_pow=psd_pow, psd_frq=psd_frq)
end

"""
    s_psd(signal; fs, norm=false)

Calculate power spectrum density of the `signal`.

# Arguments
- `signal::Array{Float64, 3}`
- `fs::Int64`: sampling rate
- `norm::Bool`: normalize do dB
- `mt::Bool=false`: if true use multi-tapered periodogram

# Returns

Named tuple containing:
- `psd_pow::Array{Float64, 3}`
- `psd_frq::Array{Float64, 3}`
"""
function s_psd(signal::Array{Float64, 3}; fs::Int64, norm::Bool=false, mt::Bool=false)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))

    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)
    psd_tmp, frq_tmp = s_psd(signal[1, :, 1], fs=fs, norm=norm, mt=mt)
    psd_pow = zeros(channel_n, length(psd_tmp), epoch_n)
    psd_frq = zeros(channel_n, length(frq_tmp), epoch_n)

    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            s = @view signal[channel_idx, :, epoch_idx]
            psd_pow[channel_idx, :, epoch_idx], psd_frq[channel_idx, :, epoch_idx] = s_psd(s, fs=fs, norm=norm, mt=mt)
        end
    end
    
    return (psd_pow=psd_pow, psd_frq=psd_frq)
end

"""
    s_stationarity_hilbert(signal::Vector{Float64})

Calculate phase stationarity using Hilbert transformation.

# Arguments

- `signal::AbstractArray`

# Returns

- `phase_stationarity::Vector{Float64}`
"""
function s_stationarity_hilbert(signal::AbstractArray)
    
    phase_stationarity = diff(DSP.unwrap(angle.(hilbert(signal))))
    
    return phase_stationarity
end

"""
    s_stationarity_mean(signal)

Calculate mean stationarity.

# Arguments

- `signal::AbstractArray`
- `window::Int64`: time window in samples

# Returns

- `mean_stationarity::Vector{Float64}`
"""
function s_stationarity_mean(signal::AbstractArray; window::Int64)

    signal = signal[1:(window * floor(Int64, length(signal) / window))]
    signal = reshape(signal, Int(length(signal) / window), window)
    mean_stationarity = mean(signal, dims=1)

    return mean_stationarity
end

"""
    s_stationarity_var(signal)

Calculate variance stationarity.

# Arguments

- `signal::AbstractArray`
- `window::Int64`: time window in samples

# Returns

- `var_stationarity::Vector{Float64}`
"""
function s_stationarity_var(signal::AbstractArray; window::Int64)

    signal = signal[1:(window * floor(Int64, length(signal) / window))]
    signal = reshape(signal, Int(length(signal) / window), window)
    var_stationarity = var(signal, dims=1)

    return var_stationarity
end

"""
    s_trim(signal; len::Int64)

Remove `len` samples from the beginning (`from` = :start, default) or end (`from` = :end) of the `signal`.

# Arguments

- `signal::AbstractArray`
- `len::Int64`: trimming length in samples
- `offset::Int64`: offset from which trimming starts, only works for `from` = :start
- `from::Symbol[:start, :end]

# Returns

- `s_trimmed::Vector{Float64}`
"""
function s_trim(signal::AbstractArray; len::Int64, offset::Int64=1, from::Symbol=:start)

    from in [:start, :end] || throw(ArgumentError("from must be :start or :end."))
    len < 1 && throw(ArgumentError("len must be ≥ 1."))
    len >= length(signal) && throw(ArgumentError("len must be < $(length(signal))."))
    offset < 1 && throw(ArgumentError("offset must be ≥ 1."))
    offset >= length(signal) - 1 && throw(ArgumentError("offset must be < $(length(signal))."))
    (from ===:start && 1 + offset + len > length(signal)) && throw(ArgumentError("offset + len must be < $(length(signal))."))
    
    offset == 1 && (from === :start && (s_trimmed = signal[(offset + len):end]))
    offset > 1 && (from === :start && (s_trimmed = vcat(signal[1:offset], signal[(1 + offset + len):end])))
    from === :end && (s_trimmed = signal[1:(end - len)])
    
    return s_trimmed
end

"""
    s2_mi(signal1::AbstractArray, signal2::AbstractArray)

Calculate mutual information between `signal1` and `signal2`.

# Arguments

- `signal1::AbstractArray`
- `signal2::AbstractArray`

# Returns

- `mi::Float64`
"""
function s2_mi(signal1::AbstractArray, signal2::AbstractArray)

    mi = get_mutual_information(signal1, signal2)
    
    return mi
end

"""
    s_entropy(signal)

Calculate entropy of `signal`.

# Arguments

- `signal::AbstractArray`

# Returns

- `ent::Float64`
"""
function s_entropy(signal::AbstractArray)

    n = length(signal)
    maxmin_range = maximum(signal) - minimum(signal)
    fd_bins = ceil(Int64, maxmin_range/(2.0 * iqr(signal) * n^(-1/3))) # Freedman-Diaconis

    # recompute entropy with optimal bins for comparison
    h = StatsKit.fit(Histogram, signal, nbins=fd_bins)
    hdat1 = h.weights ./ sum(h.weights)

    # convert histograms to probability values
    ent = -sum(hdat1 .* log2.(hdat1 .+ eps()))    

    return ent
end

"""
    s_negentropy(signal)

Calculate negentropy of `signal`.

# Arguments

- `signal::AbstractArray`

# Returns

- `ent::Float64`
"""
function s_negentropy(signal::AbstractArray)

    s = s_demean(signal)
    ne = 0.5 * log(2 * pi * exp(1) * var(s)) - s_entropy(s)

    return ne
end

"""
    s_average(signal)

Average all channels of `signal`.

# Arguments

- `signal::Array{Float64, 3}`

# Returns

- `s_averaged::Array{Float64, 3}`
"""
function s_average(signal::Array{Float64, 3})

    return mean(signal, dims=1)
end

"""
    s2_average(signal1, signal2)

Averages `signal1` and `signal2`.

# Arguments

- `signal1::AbstractArray`
- `signal2::AbstractArray`

# Returns

- `s_averaged::Vector{Float64}`
"""
function s2_average(signal1::AbstractArray, signal2::AbstractArray)

    return mean(hcat(signal1, signal2), dims=2)
end

"""
    s2_tcoherence(signal1, signal2)

Calculate coherence (mean over time), IC (imaginary coherence) and MSC (magnitude-squared coherence) between `signal1` and `signal2`.

# Arguments

- `signal1::AbstractArray`
- `signal2::AbstractArray`

# Returns

Named tuple containing:
- `c::Vector{Float64}`: coherence
- `msc::Vector{Float64}`: magnitude-squares coherence
- `ic::Vector{Float64}`: imaginary part of coherence
"""
function s2_tcoherence(signal1::AbstractArray, signal2::AbstractArray)

    length(signal1) == length(signal2) || throw(ArgumentError("Both signals must have the same length."))

    signal1 = _reflect(signal1)
    signal2 = _reflect(signal2)

    s1_fft = fft(signal1) ./ length(signal1)
    s2_fft = fft(signal2) ./ length(signal2)

    coh = (abs.((s1_fft) .* conj.(s2_fft)).^2) ./ (s1_fft .* s2_fft)
    coh = _chop(coh)
    msc = @. abs(coh)^2

    return (c=real.(coh), msc=msc, ic=imag.(coh))
end

"""
    s_pca(signal, n)

Calculates `n` first PCs for `signal`.

# Arguments

- `signal::Array{Float64, 3}`
- `n::Int64`: number of PCs

# Returns

- `pc::Array{Float64, 3}:`: PC(1)..PC(n) × epoch
- `pc_var::Matrix{Float64}`: variance of PC(1)..PC(n) × epoch
- `pc_m::PCA{Float64}`: PC mean
"""
function s_pca(signal::Array{Float64, 3}; n::Int64)

    n < 0 && throw(ArgumentError("n must be ≥ 1."))
    n > size(signal, 1) && throw(ArgumentError("Number of PCs must be ≤ $(size(signal, 1))."))

    channel_n, _, epoch_n = size(signal)
    pc_m = []

    # check maximum n
    n_tmp = n
    Threads.@threads for epoch_idx in 1:epoch_n
        s = @view signal[:, :, epoch_idx]
        pc_m = MultivariateStats.fit(PCA, s, maxoutdim=n)
        size(pc_m, 2) < n_tmp && (n_tmp = size(pc_m, 2))
    end
    n_tmp < n && @warn "Only $n_tmp PC components were generated."
    n = n_tmp
    
    pc = zeros(n, size(signal, 2), epoch_n)
    pc_var = zeros(n, epoch_n)
    pc_reconstructed = zeros(size(signal))

    Threads.@threads for epoch_idx in 1:epoch_n
        s = @view signal[:, :, epoch_idx]
        # m_cov = s_cov(s)
        # eig_val, eig_vec = eigen(m_cov)
        # eig_val_idx = sortperm(eig_val, rev=true)
        # eig_val = eig_val[eig_val_idx]
        # eig_vec = m_sort(eig_vec, eig_val_idx)
        # eig_val = 100 .* eig_val / sum(eig_val) # convert to %

        pc_m = MultivariateStats.fit(PCA, s, maxoutdim=n)
        v = principalvars(pc_m) ./ var(pc_m) * 100

        for idx in 1:n
            pc_var[idx, epoch_idx] = v[idx]
            # pc[idx, :, epoch_idx] = (eig_vec[:, idx] .* s)[idx, :]
            pc[idx, :, epoch_idx] = MultivariateStats.predict(pc_m, s)[idx, :]
        end
    end

    return pc, pc_var, pc_m
end


"""
    s_pca_reconstruct(signal, pc, pcm)

Reconstructs `signal` using PCA components.

# Arguments

- `signal::Array{Float64, 3}`
- `pc::Array{Float64, 3}:`: IC(1)..IC(n) × epoch
- `pc_m::PCA{Float64}:`: IC(1)..IC(n) × epoch

# Returns

- `s_reconstructed::Array{Float64, 3}`
"""
function s_pca_reconstruct(signal::Array{Float64, 3}; pc::Array{Float64, 3}, pc_m::PCA{Float64})

    s_reconstructed = zeros(size(signal))

    _, _, epoch_n = size(signal)

    @inbounds @simd for epoch_idx in 1:epoch_n
        s_reconstructed[:, :, epoch_idx] = reconstruct(pc_m, pc[:, :, epoch_idx])
    end

    return s_reconstructed
end

"""
    s_fconv(signal; kernel, norm)

Perform convolution in the frequency domain between `signal` and `kernel`.

# Arguments

- `signal::AbstractArray`
- `kernel::Union{Vector{<:Real}, Vector{ComplexF64}}`
- `norm::Bool=false`: normalize kernel

# Returns

- `s_conv::Vector{ComplexF64}`
"""
function s_fconv(signal::AbstractArray; kernel::Union{Vector{<:Real}, Vector{ComplexF64}}, norm::Bool=false)

    n_signal = length(signal)
    n_kernel = length(kernel)
    n_conv = n_signal + n_kernel - 1
    half_kernel = floor(Int64, n_kernel / 2)
    s_fft = fft0(signal, n_conv)
    kernel_fft = fft0(kernel, n_conv)
    norm == true && (kernel_fft ./= cmax(kernel_fft))
    s_conv = ifft(s_fft .* kernel_fft)
    
    # remove in- and out- edges
    if mod(n_kernel, 2) == 0 
        s_conv = s_conv[half_kernel:(end - half_kernel)]
    else
        s_conv = s_conv[(half_kernel + 1):(end - half_kernel)]
    end

    return s_conv
end

"""
    s_ica(signal, n, tol=1.0e-6, iter=100, f=:tanh)

Calculate `n` first ICs for `signal`.

# Arguments

- `signal::Array{Float64, 3}`
- `n::Int64`: number of PCs
- `tol::Float64`: tolerance for ICA
- `iter::Int64`: maximum number of iterations
- `f::Symbol[:tanh, :gaus]`: neg-entropy functor

# Returns

- `ic::Array{Float64, 3}:`: IC(1)..IC(n) × epoch
- `ic_mw::Array{Float64, 3}:`: IC(1)..IC(n) × epoch
"""
function s_ica(signal::Array{Float64, 3}; n::Int64, tol::Float64=1.0e-6, iter::Int64=100, f::Symbol=:tanh)

    f in [:tanh, :gaus] || throw(ArgumentError("f must be :tanh or :gaus."))
    n < 0 && throw(ArgumentError("n must be ≥ 1."))
    n > size(signal, 1) && throw(ArgumentError("Number of ICs must be ≤ $(size(signal, 1))."))
    channel_n, _, epoch_n = size(signal)
    ic = zeros(n, size(signal, 2), epoch_n)
    ic_mw = zeros(channel_n, n, epoch_n)

    @inbounds @simd for epoch_idx in 1:epoch_n
        s = @view signal[:, :, epoch_idx]

        f === :tanh && (M = MultivariateStats.fit(ICA, s, n, tol=tol, maxiter=iter, fun=MultivariateStats.Tanh(1.0)))
        f === :gaus && (M = MultivariateStats.fit(ICA, s, n, tol=tol, maxiter=iter, fun=MultivariateStats.Gaus()))

        n == size(signal, 1) && (mw = inv(M.W)')
        n < size(signal, 1) && (mw = pinv(M.W)')

        for idx in 1:n
            ic[idx, :, epoch_idx] = MultivariateStats.predict(M, s)[idx, :]
        end

        ic_mw[:, :, epoch_idx] = mw
    end

    return ic, ic_mw
end

"""
    s_ica_reconstruct(signal, ic, ic_mw, ic_v)

Reconstructs `signal` using removal of `ic_v` ICA components.

# Arguments

- `signal::Array{Float64, 3}`
- `ic::Array{Float64, 3}:`: IC(1)..IC(n) × epoch
- `ic_mw::Array{Float64, 3}:`: IC(1)..IC(n) × epoch
- `ic_v::Union{Int64, Vector{Int64}, AbstractRange} - list of ICs to remove

# Returns

- `s_reconstructed::Array{Float64, 3}`
"""
function s_ica_reconstruct(signal::Array{Float64, 3}; ic::Array{Float64, 3}, ic_mw::Array{Float64, 3}, ic_v::Union{Int64, Vector{Int64}, AbstractRange})

    typeof(ic_v) <: AbstractRange && (ic_v = collect(ic_v))
    if typeof(ic_v) == Vector{Int64}
        sort!(ic_v)
        for idx in 1:length(ic_v)
            (ic_v[idx] < 1 || ic_v[idx] > size(ic_mw, 2)) && throw(ArgumentError("ic_v must be ≥ 1 and ≤ $(size(ic_mw, 2))"))
        end
    else
        (ic_v < 1 || ic_v > size(ic_mw, 2)) && throw(ArgumentError("ic_v must be ≥ 1 and ≤ $(size(ic_mw, 2))"))
    end
    ic_removal = setdiff(1:size(ic_mw, 2), ic_v)

    s_reconstructed = zeros(size(signal))

    _, _, epoch_n = size(signal)

    @inbounds @simd for epoch_idx in 1:epoch_n
        s_reconstructed[:, :, epoch_idx] = ic_mw[:, ic_removal, epoch_idx] * ic[ic_removal, :, epoch_idx]
    end

    return s_reconstructed
end

"""
    s_spectrogram(signal; fs, norm=true, demean=true)

Calculate spectrogram of `signal`.

# Arguments

- `signal::AbstractArray`
- `fs::Int64`: sampling frequency
- `norm::Bool=true`: normalize powers to dB
- `mt::Bool=false`: if true use multi-tapered spectrogram
- `demean::Bool=true`: demean signal prior to analysis

# Returns

Named tuple containing:
- `s_pow::Matrix{Float64}`: powers
- `s_frq::Vector{Float64}`: frequencies
- `s_t::Vector{Float64}`: time
"""
function s_spectrogram(signal::AbstractArray; fs::Int64, norm::Bool=true, mt::Bool=false, demean::Bool=true)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1 Hz."))

    demean == true && (signal = s_demean(signal))

    nfft = length(signal)
    interval = fs
    overlap = round(Int64, fs * 0.85)

    mt == false && (spec = spectrogram(signal, interval, overlap, nfft=nfft, fs=fs, window=hanning))
    mt == true && (spec = mt_spectrogram(signal, fs=fs))
    s_t = collect(spec.time)
    s_frq = Vector(spec.freq)
    if norm == true
        s_pow = pow2db.(spec.power)
    else
        s_pow = Vector(spec.power)
    end

    return (s_pow=s_pow, s_frq=s_frq, s_t=s_t)
end

"""
    s_detect_epoch_flat(signal)

Detect bad `signal` epochs based on: flat channel(s)

# Arguments

- `signal::Array{Float64, 3}`

# Returns

- `bad_epochs_score::Vector{Int64}`: percentage of bad channels per epoch
"""
function s_detect_epoch_flat(signal::Array{Float64, 3})
    
    channel_n, _, epoch_n = size(signal)

    bad_epochs_score = zeros(epoch_n)

    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            s = @view signal[channel_idx, :, epoch_idx]
            c_tmp_diff = abs.(diff(s))
            # add tolerance around zero μV
            sum(c_tmp_diff) < eps() && (bad_epochs_score[epoch_idx] += 1)
        end
    end

    bad_epochs_score = round.(bad_epochs_score ./ channel_n, digits=1)

    return bad_epochs_score
end

"""
    s_detect_epoch_rmse(signal)

Detect bad `signal` epochs based on: RMSE vs average channel > 95%CI.

# Arguments

- `signal::Array{Float64, 3}`

# Returns

- `bad_epochs_score::Vector{Int64}`: percentage of bad channels per epoch
"""
function s_detect_epoch_rmse(signal::Array{Float64, 3})
    
    channel_n, _, epoch_n = size(signal)

    bad_epochs_score = zeros(epoch_n)

    @inbounds @simd for epoch_idx in 1:epoch_n
        ch_m = vec(median(signal[:, :, epoch_idx], dims=1))
        rmse_ch = zeros(channel_n)
        Threads.@threads for channel_idx in 1:channel_n
            rmse_ch[channel_idx] = s2_rmse(signal[channel_idx, :, epoch_idx], ch_m)
        end
        Threads.@threads for channel_idx in 1:channel_n
            rmse_ch[channel_idx] > HypothesisTests.confint(OneSampleTTest(rmse_ch))[2] && (bad_epochs_score[epoch_idx] += 1)
        end
    end

    bad_epochs_score = round.(bad_epochs_score ./ channel_n, digits=1)

    return bad_epochs_score
end

"""
    detect_epoch_rmsd(signal)

Detect bad `signal` epochs based on: RMSD vs average channel > 95%CI.
# Arguments

- `signal::Array{Float64, 3}`

# Returns

- `bad_epochs_score::Vector{Int64}`: percentage of bad channels per epoch
"""
function s_detect_epoch_rmsd(signal::Array{Float64, 3})
    
    channel_n, _, epoch_n = size(signal)

    bad_epochs_score = zeros(epoch_n)

    @inbounds @simd for epoch_idx in 1:epoch_n
        ch_m = median(signal[:, :, epoch_idx], dims=1)
        rmsd_ch = zeros(channel_n)
        Threads.@threads for channel_idx in 1:channel_n
            rmsd_ch[channel_idx] = Distances.rmsd(signal[channel_idx, :, epoch_idx], ch_m)
        end
        Threads.@threads for channel_idx in 1:channel_n
            rmsd_ch[channel_idx] > HypothesisTests.confint(OneSampleTTest(rmsd_ch))[2] && (bad_epochs_score[epoch_idx] += 1)
        end
    end

    bad_epochs_score = round.(bad_epochs_score ./ channel_n, digits=1)

    return bad_epochs_score
end

"""
    s_detect_epoch_euclid(signal)

Detect bad `signal` epochs based on: Euclidean distance vs median channel > 95% CI.

# Arguments

- `signal::Array{Float64, 3}`

# Returns

- `bad_epochs_score::Vector{Int64}`: percentage of bad channels per epoch
"""
function s_detect_epoch_euclid(signal::Array{Float64, 3})
    
    channel_n, _, epoch_n = size(signal)

    bad_epochs_score = zeros(epoch_n)

    @inbounds @simd for epoch_idx in 1:epoch_n
        ch_m = median(signal[:, :, epoch_idx], dims=1)
        ed_ch = zeros(channel_n)
        Threads.@threads for channel_idx in 1:channel_n
            ed_ch[channel_idx] = euclidean(signal[channel_idx, :, epoch_idx], ch_m)
        end
        Threads.@threads for channel_idx in 1:channel_n
            ed_ch[channel_idx] > HypothesisTests.confint(OneSampleTTest(ed_ch))[2] && (bad_epochs_score[epoch_idx] += 1)
        end
    end

    bad_epochs_score = round.(bad_epochs_score ./ channel_n, digits=1)

    return bad_epochs_score
end

"""
    s_detect_epoch_p2p(signal)

Detect bad `signal` epochs based on: p2p amplitude > upper 95% CI p2p amplitude.

# Arguments

- `signal::Array{Float64, 3}`

# Returns

- `bad_epochs_score::Vector{Int64}`: percentage of bad channels per epoch
"""
function s_detect_epoch_p2p(signal::Array{Float64, 3})
    
    channel_n, _, epoch_n = size(signal)

    bad_epochs_score = zeros(epoch_n)

    @inbounds @simd for epoch_idx in 1:epoch_n
        p2p = zeros(channel_n)
        Threads.@threads for channel_idx in 1:channel_n
            p2p[channel_idx] = maximum(signal[channel_idx, :, epoch_idx]) + abs(minimum(signal[channel_idx, :, epoch_idx]))
        end
        Threads.@threads for channel_idx in 1:channel_n
            p2p[channel_idx] > HypothesisTests.confint(OneSampleTTest(p2p))[2] && (bad_epochs_score[epoch_idx] += 1)
        end
    end

    bad_epochs_score = round.(bad_epochs_score ./ channel_n, digits=1)

    return bad_epochs_score
end

"""
    s_snr(signal)

Calculate SNR of `signal`.

# Arguments

- `signal::AbstractArray`

# Returns

- `snr::Float64`: SNR

# Source

D. J. Schroeder (1999). Astronomical optics (2nd ed.). Academic Press. ISBN 978-0-12-629810-9, p.278
"""
function s_snr(signal::AbstractArray)

    # make signal positive
    signal .+= abs(minimum(signal))
    snr = mean(signal) / std(signal)

    return snr
end

"""
    s_findpeaks(signal; d)

Find peaks in `signal`.

# Arguments

- `signal::AbstractArray`
- `d::Int64=32`: distance between peeks in samples

# Returns

- `p_idx::Vector{Int64}`

"""
function s_findpeaks(signal::AbstractArray; d::Int64=32)
    p_idx, = findpeaks1d(signal, distance=d)
    
    return p_idx
end

"""
    s_wdenoise(signal; wt)

Perform wavelet denoising.

# Arguments

- `signal::AbstractArray`
- `wt::Symbol=:db4`: wavelet type: :db2, :db4, :db8, :db10, :haar, :coif2, :coif4, :coif8

# Returns

- `signal_denoised::Vector{Float64}`
"""
function s_wdenoise(signal::AbstractArray; wt::Symbol=:db4)
    
    wt in [:db2, :db4, :db8, :db10, :haar, :coif2, :coif4, :coif8] || throw(ArgumentError("wt must be :db2, :db4, :db8, :db10, :haar, :coif2, :coif4, :coif8"))

    wt === :db2 && (wt = wavelet(WT.db2))
    wt === :db4 && (wt = wavelet(WT.db4))
    wt === :db8 && (wt = wavelet(WT.db8))
    wt === :db10 && (wt = wavelet(WT.db10))
    wt === :haar && (wt = wavelet(WT.haar))
    wt === :coif2 && (wt = wavelet(WT.coif2))
    wt === :coif4 && (wt = wavelet(WT.coif4))
    wt === :coif8 && (wt = wavelet(WT.coif8))

    signal_denoised = denoise(signal, wt)

    return signal_denoised
end

"""
    effsize(x1, x2)

Calculate Cohen's d and Hedges g effect sizes.

# Arguments

- `x1::Vector{Float64}`
- `x2::Vector{Float64}`

# Returns

Named tuple containing:
- `d::Float64`: Cohen's d
- `g::Float64`: Hedges g
"""
function effsize(x1::Vector{<:Real}, x2::Vector{<:Real})
    d = (mean(x2) - mean(x1)) / sqrt((std(x1)^2 + std(x2)^2) / 2)
    g = (mean(x2) - mean(x1)) / sqrt((((length(x1) - 1) * (std(x1)^2)) + ((length(x2) - 1) * (std(x2)^2))) / (length(x1) + length(x2) - 2))
    return (cohen=d, hedges=g)
end

"""
    s_ispc(signal1, signal2)

Calculate ISPC (Inter-Site-Phase Clustering) between `signal1` and `signal2`.

# Arguments

- `signal1::AbstractArray`
- `signal2::AbstractArray`

# Returns

Named tuple containing:
- `ispc::Float64`: ISPC value
- `ispc_angle::Float64`: ISPC angle
- `signal_diff::Vector{Float64}`: signal difference (signal2 - signal1)
- `phase_diff::Vector{Float64}`: phase difference (signal2 - signal1)
- `s1_phase::Vector{Float64}`: signal 1 phase
- `s2_phase::Vector{Float64}`: signal 2 phase
"""
function s_ispc(signal1::AbstractArray, signal2::AbstractArray)

    length(signal1) == length(signal2) || throw(ArgumentError("Both signals must have the same length."))

    _, _, _, s1_phase = s_hspectrum(signal1)
    _, _, _, s2_phase = s_hspectrum(signal2)

    signal_diff = signal2 - signal1
    phase_diff = s2_phase - s1_phase

    ispc = abs(mean(exp.(1im .* phase_diff)))
    ispc_angle = angle(mean(exp.(1im .* phase_diff)))

    return (ispc=ispc, ispc_angle=ispc_angle, signal_diff=signal_diff, phase_diff=phase_diff, s1_phase=s1_phase, s2_phase=s2_phase)
end

"""
    s_itpc(signal; t)

Calculate ITPC (Inter-Trial-Phase Clustering) over epochs/trials at time `t` of `signal`.

# Arguments

- `signal::AbstractArray`
- `t::Int64`: time point (sample number) at which ITPC is calculated
- `w::Union{Vector{<:Real}, Nothing}`: optional vector of epochs/trials weights for wITPC calculation

# Returns

Named tuple containing:
- `itpc::Float64`: ITPC value
- `itpcz::Float64`: Rayleigh's ITPC Z value
- `itpc_angle::Float64`: ITPC angle
- `itpc_phases::Vector{Float64}`: phases at time `t` averaged across trials/epochs
"""
function s_itpc(signal::AbstractArray; t::Int64, w::Union{Vector{<:Real}, Nothing}=nothing)

    t < 1 && throw(ArgumentError("t must be ≥ 1."))
    size(signal, 1) == 1 || throw(ArgumentError("signals must have 1 channel."))
    t > size(signal, 2) && throw(ArgumentError("t must be ≤ $(size(signal, 2))."))
    epoch_n = size(signal, 3)

    w === nothing && (w = ones(epoch_n))
    # scale w if w contains negative values
    any(i -> i < 0, w) && (w .+= abs(minimum(w)))
    length(w) == epoch_n || throw(ArgumentError("Length of w should be equal to number of epochs ($epoch_n)."))
    
    s_phase = zeros(size(signal, 2), epoch_n)
    @inbounds @simd for epoch_idx in 1:epoch_n
        _, _, _, s_phase[:, epoch_idx] = s_hspectrum(signal[1, :, epoch_idx])
    end
 
    itpc_phases = s_phase[t, :]
    itpc = abs.(mean(exp.(1im .* itpc_phases .* w)))
    itpc_angle = angle.(mean(exp.(1im .* itpc_phases .* w)))
    itpcz = epoch_n * itpc^2

    return (itpc=itpc, itpcz=itpcz, itpc_angle=itpc_angle, itpc_phases=itpc_phases)
end

"""
    s_pli(signal1, signal2)

Calculate PLI (Phase-Lag Index) between `signal1` and `signal2`.

# Arguments

- `signal1::AbstractArray`
- `signal2::AbstractArray`

# Returns

Named tuple containing:
- `pli::Float64`: PLI value
- `signal_diff::Vector{Float64}`: signal difference (signal2 - signal1)
- `phase_diff::Vector{Float64}`: phase difference (signal2 - signal1)
- `s1_phase::Vector{Float64}`: signal 1 phase
- `s2_phase::Vector{Float64}`: signal 2 phase
"""
function s_pli(signal1::AbstractArray, signal2::AbstractArray)

    length(signal1) == length(signal2) || throw(ArgumentError("Both signals must have the same length."))

    _, _, _, s1_phase = s_hspectrum(signal1)
    _, _, _, s2_phase = s_hspectrum(signal2)

    signal_diff = signal2 - signal1
    phase_diff = s2_phase - s1_phase

    pli = abs(mean(sign.(imag.(exp.(1im .* phase_diff)))))

    return (pli=pli, signal_diff=signal_diff, phase_diff=phase_diff, s1_phase=s1_phase, s2_phase=s2_phase)
end

"""
    s_ged(signal1, signal2)

Perform generalized eigendecomposition between `signal1` and `signal2`.

# Arguments

- `signal1::AbstractArray`: signal to be analyzed
- `signal2::AbstractArray`: original signal

# Returns

Named tuple containing:
- `sged::AbstractArray`
- `ress::AbstractArray`
- `ress_normalized::AbstractArray`: RESS normalized to -1..1
"""
function s_ged(signal1::AbstractArray, signal2::AbstractArray)

    size(signal1) == size(signal2) || throw(ArgumentError("signal1 and signal2 must have the same size."))

    channel_n = size(signal1, 1)

    s1cov = cov(signal1')
    s2cov = cov(signal2')
    eig_val, eig_vec = eigen(s1cov, s2cov)
    eig_val_idx = sortperm(eig_val, rev=true)
    eig_val = eig_val[eig_val_idx]
    eig_vec = m_sort(eig_vec, eig_val_idx, dims=2)
    sged = signal2 .* eig_vec[:, 1]
    ress = pinv(eig_vec[:, 1]')
    ress_normalized = ress ./ maximum(abs.(ress))

    return (sged=sged, ress=ress, ress_normalized=ress_normalized)
end

"""
    s_frqinst(signal; fs)

Calculate instantaneous frequency `signal`.

# Arguments

- `signal::AbstractArray`
- `fs::Int64`

# Returns

- `frqinst::Vector{Float64}`
"""
function s_frqinst(signal::AbstractArray; fs::Int64)

    fs < 0 && throw(ArgumentError("fs must be > 0."))

    _, _, _, h_phases = s_hspectrum(signal)
    frqinst = 256 * s_derivative(h_phases) / (2*pi)

    return frqinst
end

"""
    s_hspectrum(signal; pad=0)

Calculate amplitudes, powers and phases of the `signal` using Hilbert transform.

# Arguments

- `signal::AbstractArray`
- `pad::Int64`: pad the `signal` with `pad` zeros

# Returns

Named tuple containing:
- `h::Vector(ComplexF64}`: Hilbert components
- `h_amplitudes::Vector{Float64}`
- `h_powers::Vector{Float64}`
- `h_phases::Vector{Float64}`
"""
function s_hspectrum(signal::AbstractArray; pad::Int64=0)

    pad < 0 && throw(ArgumentError("pad must be ≥ 0."))
    h = hilbert(pad0(signal, pad))

    # amplitudes
    h_amplitudes = @. abs(h)
    # power
    h_powers = h_amplitudes.^2
    # phases
    h_phases = angle.(h)

    return (h=h, h_amplitudes=h_amplitudes, h_powers=h_powers, h_phases=h_phases)
end

"""
    t2f(t)

Convert cycle length in ms `t` to frequency.

# Arguments

- `t::Real`: cycle length in ms

# Returns

- `f::Float64`: frequency in Hz
"""
function t2f(t::Real)

    t < 0 && throw(ArgumentError("t must be > 0."))
    f = round(1000 / t, digits=2)

    return f
end

"""
    f2t(f)

Convert frequency `f` to cycle length in ms.

# Arguments

- `f::Real`: frequency in Hz

# Returns

- `f::Float64`: cycle length in ms
"""
function f2t(f::Real)

    f < 0 && throw(ArgumentError("t must be > 0."))
    f = round(1000 / f, digits=2)
    
    return f
end

"""
    s_wspectrogram(signal; pad, norm, frq_lim, frq_n, frq, fs, ncyc, demean)

Calculate spectrogram of the `signal` using wavelet convolution.

# Arguments

- `signal::AbstractArray`
- `pad::Int64`: pad the `signal` with `pad` zeros
- `norm::Bool=true`: normalize powers to dB
- `frq_lim::Tuple{Real, Real}`: frequency bounds for the spectrogram
- `frq_n::Int64`: number of frequencies
- `frq::Symbol=:log`: linear (:lin) or logarithmic (:log) frequencies
- `fs::Int64`: sampling rate
- `ncyc::Int64=6`: number of cycles for Morlet wavelet
- `demean::Bool=true`: demean signal prior to analysis

# Returns

Named tuple containing:
- `w_conv::Matrix(ComplexF64}`: convoluted signal
- `w_powers::Matrix{Float64}`
- `w_phases::Matrix{Float64}`
- `frq_list::Vector{Float64}`
"""
function s_wspectrogram(signal::AbstractArray; pad::Int64=0, norm::Bool=true, frq_lim::Tuple{Real, Real}, frq_n::Int64, frq::Symbol=:lin, fs::Int64, ncyc::Int64=6, demean::Bool=true)

    fs <= 0 && throw(ArgumentError("fs must be > 0."))
    pad < 0 && throw(ArgumentError("pad must be ≥ 0."))
    frq in [:log, :lin] || throw(ArgumentError("frq must be :log or :lin."))
    frq_lim = tuple_order(frq_lim)
    frq_lim[1] < 0 && throw(ArgumentError("Lower frequency bound must be ≥ 0."))
    frq_lim[2] > fs ÷ 2 && throw(ArgumentError("Upper frequency bound must be ≤ $(fs ÷ 2)."))
    frq_n < 2 && throw(ArgumentError("frq_n frequency bound must be ≥ 2."))
    frq_lim[1] == 0 && (frq_lim = (0.1, frq_lim[2]))
    if frq === :log
        frq_lim = (frq_lim[1], frq_lim[2])
        frq_list = round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), frq_n), digits=1)
    else
        frq_list = linspace(frq_lim[1], frq_lim[2], frq_n)
    end

    demean == true && (signal = s_demean(signal))
    pad > 0 && (signal = pad0(signal, pad))
    w_conv = zeros(ComplexF64, length(frq_list), length(signal))
    w_powers = zeros(length(frq_list), length(signal))
    w_phases = zeros(length(frq_list), length(signal))
    @inbounds @simd for frq_idx in 1:frq_n
        kernel = generate_morlet(fs, frq_list[frq_idx], 1, ncyc=ncyc, complex=true)
        w_conv[frq_idx, :] = s_fconv(signal, kernel=kernel, norm=true)
        w_powers[frq_idx, :] = @. abs(w_conv[frq_idx, :])^2
        w_phases[frq_idx, :] = @. angle(w_conv[frq_idx, :])
    end
    norm == true && (w_powers = pow2db.(w_powers))

    return (w_conv=w_conv, w_powers=w_powers, w_phases=w_phases, frq_list=frq_list)
end

"""
   s_fftdenoise(signal::AbstractArray; pad::Int64=0, threshold::Int64=100) 

Perform FFT denoising.

# Arguments

- `signal::AbstractArray`
- `pad::Int64=0`: pad the `signal` with `pad` zeros
- `threshold::Int64=100`: PSD threshold for keeping frequency components

# Returns

- `signal_denoised::Vector{Float64}`
"""
function s_fftdenoise(signal::AbstractArray; pad::Int64=0, threshold::Int64=100)

    signal_fft = fft0(signal, pad)
    signal_psd = real.(signal_fft .* conj.(signal_fft)) / length(signal)
    # signal_psd = abs2.(signal_fft) / length(signal)

    # zero frequencies
    signal_idx = signal_psd .> threshold
    signal_psd .*= signal_idx
    signal_fft .*= signal_idx
    signal_denoised = real.(ifft(signal_fft))

    return signal_denoised
end

"""
    s_gfilter(signal, fs, f, gw)

Filter `signal` using Gaussian in the frequency domain.

# Arguments

- `signal::AbstractArray`
- `fs::Int64`: sampling rate
- `f::Real`: filter frequency
- `gw::Real=5`: Gaussian width in Hz

# Returns

Named tuple containing:
- `s_f::Vector{Float64}`
"""
function s_gfilter(signal::Vector{Float64}; fs::Int64, f::Real, gw::Real=5)

    # add reflected signals to reduce edge artifacts
    s_r = _reflect(signal)
    
    # create Gaussian in frequency domain
    gf = linspace(0, fs, length(s_r))
    gs = (gw * (2 * pi - 1)) / (4 * pi)     # normalized width
    gf .-= f                                # shifted frequencies
    # g = @. exp(-0.5 * (gf / gs)^2)          # Gaussian
    g = @. exp((-gf^2 ) / 2 * gs^2)         # Gaussian
    g ./= abs(maximum(g))                   # gain-normalized

    # filter
    s_f = 2 .* real.(ifft(fft(s_r).*g))
    
    # remove reflected signals
    s_f = _chop(s_f)

    return s_f
end

"""
    s_ghspectrogram(signal; pad, norm, frq_lim, frq_n, frq, fs, gw, demean)

Calculate spectrogram of the `signal` using Gaussian and Hilbert transform.

# Arguments

- `signal::AbstractArray`
- `fs::Int64`: sampling rate
- `norm::Bool=true`: normalize powers to dB
- `frq_lim::Tuple{Real, Real}`: frequency bounds for the spectrogram
- `frq_n::Int64`: number of frequencies
- `frq::Symbol=:log`: linear (:lin) or logarithmic (:log) frequencies
- `gw::Real=5`: Gaussian width in Hz
- `demean::Bool`=true: demean signal prior to analysis

# Returns

Named tuple containing:
- `h_powers::Matrix{Float64}`
- `frq_list::Vector{Float64}`
"""
function s_ghspectrogram(signal::AbstractArray; fs::Int64, norm::Bool=true, frq_lim::Tuple{Real, Real}, frq_n::Int64, frq::Symbol=:lin, gw::Real=5, demean::Bool=true)

    fs <= 0 && throw(ArgumentError("fs must be > 0."))
    frq in [:log, :lin] || throw(ArgumentError("frq must be :log or :lin."))
    frq_lim = tuple_order(frq_lim)
    frq_lim[1] < 0 && throw(ArgumentError("Lower frequency bound must be ≥ 0."))
    frq_lim[2] > fs ÷ 2 && throw(ArgumentError("Upper frequency bound must be ≤ $(fs ÷ 2)."))
    frq_n < 2 && throw(ArgumentError("frq_n frequency bound must be ≥ 2."))
    frq_lim[1] == 0 && (frq_lim = (0.1, frq_lim[2]))
    if frq === :log
        frq_lim = (frq_lim[1], frq_lim[2])
        frq_list = round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), frq_n), digits=1)
    else
        frq_list = linspace(frq_lim[1], frq_lim[2], frq_n)
    end

    demean == true && (signal = s_demean(signal))
    h_powers = zeros(length(frq_list), length(signal))
    h_phases = zeros(length(frq_list), length(signal))
    @inbounds @simd for frq_idx in 1:length(frq_list)
        s = s_gfilter(signal, fs=fs, f=frq_list[frq_idx], gw=gw)
        h_powers[frq_idx, :] = abs.(hilbert(s)).^2
        h_phases[frq_idx, :] = angle.(hilbert(s))
    end
    norm == true && (h_powers = pow2db.(h_powers))

    return (h_powers=h_powers, frq_list=frq_list)
end

"""
    s_tkeo(signal)

Calculate Teager-Kaiser energy-tracking operator: y(t) = x(t)^2 - x(t-1)x(t+1)

# Arguments

- `signal::AbstractArray`

# Returns

- `s_new::Vector{Float64}`
"""
function s_tkeo(signal::AbstractArray)
    tkeo = zeros(length(signal))
    tkeo[1] = signal[1]
    tkeo[end] = signal[end]
    for idx in 2:(length(signal) - 1)
        tkeo[idx] = signal[idx]^2 - (signal[idx - 1] * signal[idx + 1])
    end

    return tkeo
end

"""
    s_wspectrum(signal; pad, norm, frq_lim, frq_n, frq, fs, ncyc)

Calculate power spectrum of the `signal` using wavelet convolution.

# Arguments

- `signal::AbstractArray`
- `pad::Int64`: pad the `signal` with `pad` zeros
- `norm::Bool=true`: normalize powers to dB
- `frq_lim::Tuple{Real, Real}`: frequency bounds for the spectrogram
- `frq_n::Int64`: number of frequencies
- `frq::Symbol=:log`: linear (:lin) or logarithmic (:log) frequencies
- `fs::Int64`: sampling rate
- `ncyc::Int64=6`: number of cycles for Morlet wavelet

# Returns

Named tuple containing:
- `w_powers::Matrix{Float64}`
- `frq_list::Vector{Float64}`
"""
function s_wspectrum(signal::AbstractArray; pad::Int64=0, norm::Bool=true, frq_lim::Tuple{Real, Real}, frq_n::Int64, frq::Symbol=:lin, fs::Int64, ncyc::Int64=6)

    fs <= 0 && throw(ArgumentError("fs must be > 0."))
    pad < 0 && throw(ArgumentError("pad must be ≥ 0."))
    frq in [:log, :lin] || throw(ArgumentError("frq must be :log or :lin."))
    frq_lim = tuple_order(frq_lim)
    frq_lim[1] < 0 && throw(ArgumentError("Lower frequency bound must be ≥ 0."))
    frq_lim[2] > fs ÷ 2 && throw(ArgumentError("Upper frequency bound must be ≤ $(fs ÷ 2)."))
    frq_n < 2 && throw(ArgumentError("frq_n frequency bound must be ≥ 2."))
    frq_lim[1] == 0 && (frq_lim = (0.1, frq_lim[2]))
    if frq === :log
        frq_lim = (frq_lim[1], frq_lim[2])
        frq_list = round.(logspace(log10(frq_lim[1]), log10(frq_lim[2]), frq_n), digits=1)
    else
        frq_list = linspace(frq_lim[1], frq_lim[2], frq_n)
    end

    pad > 0 && (signal = pad0(signal, pad))
    w_powers = zeros(length(frq_list))
    @inbounds @simd for frq_idx in 1:frq_n
        kernel = generate_morlet(fs, frq_list[frq_idx], 1, ncyc=ncyc, complex=true)
        w_conv = s_fconv(signal, kernel=kernel, norm=true)
        w_powers[frq_idx] = mean(@. abs(w_conv)^2)
    end
    norm == true && (w_powers = pow2db.(w_powers))

    return (w_powers=w_powers, frq_list=frq_list)
end

"""
    a2_cmp(a1, a2; p, perm_n)

Compare two 3-dimensional arrays `a1` and `a2` (e.g. two spectrograms), using permutation based statistic.

# Arguments

- `a1::Array{<:Real, 3}`: first array
- `a2::Array{<:Real, 3}`: second array
- `p::Float64=0.05`: p-value
- `perm_n::Int64=1000`: number of permutations

# Returns

Named tuple containing:
- `zmap::Array{Float64, 3}`: array of Z-values
- `zmap_b::Array{Float64, 3}`: binarized mask of statistically significant positions
"""
function a2_cmp(a1::Array{<:Real, 3}, a2::Array{<:Real, 3}; p::Float64=0.05, perm_n::Int64=1000)
    size(a1) == size(a2) || throw(ArgumentError("Both arrays must have the same size"))

    spec_diff = dropdims(mean(a2, dims=3) .- mean(a1, dims=3), dims=3)
    zval = abs(norminvcdf(p))
    perm_n = 1000
    spec_all = cat(a1, a2, dims=3)
    perm_maps = zeros(size(a1, 1), size(a1, 2), perm_n)
    epoch_n = size(spec_all, 3)
    @inbounds @simd for perm_idx in 1:perm_n
        rand_idx = sample(1:epoch_n, epoch_n, replace=false)
        rand_spec = spec_all[:, :, rand_idx]
        a2 = @view rand_spec[:, :, (epoch_n ÷ 2 + 1):end]
        a1 = @view rand_spec[:, :, 1:(epoch_n ÷ 2)]
        perm_maps[:, :, perm_idx] = dropdims(mean(a2, dims=3) .- mean(a1, dims=3), dims=3)
    end
    mean_h0 = dropdims(mean(perm_maps, dims=3), dims=3)
    std_h0 = dropdims(std(perm_maps, dims=3), dims=3)

    # threshold real data
    zmap = @. (spec_diff - mean_h0) / std_h0

    # threshold at p-value
    zmap_b = deepcopy(zmap)
    @. zmap_b[abs(zmap) < zval] = 0
    @. zmap_b[zmap_b != 0] = 1
    zmap_b = Bool.(zmap_b)
    zmap_b = map(x -> !x, zmap_b)

    return (zmap=zmap, zmap_b=zmap_b)
end

"""
    s_fcoherence(signal; fs, frq)

Calculate coherence (mean over all frequencies) and MSC (magnitude-squared coherence) between channels of `signal`.

# Arguments

- `signal::AbstractArray`
- `fs::Int64`
- `frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing`: return coherence only for the given frequency range

# Returns

- `c::Array{Float64, 3}`: coherence
- `msc::Array{Float64, 3}`: MSC
- `f::Vector{Float64}`: frequencies
"""
function s_fcoherence(signal::AbstractArray; fs::Int64, frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing)

    fs < 0 && throw(ArgumentError("fs must be > 0."))
    c = mt_coherence(signal, fs=fs)
    f = Vector(c.freq)
    c = c.coherence
    if frq_lim !== nothing
        frq_lim = tuple_order(frq_lim)
        frq_lim[1] < 0 && throw(ArgumentError("Lower frequency bound must be ≥ 0."))
        frq_lim[2] > fs / 2 && throw(ArgumentError("Upper frequency bound must be ≤ $fs."))
        idx1 = vsearch(frq_lim[1], f)
        idx2 = vsearch(frq_lim[2], f)
        c = c[:, :, idx1:idx2]
        f = f[idx1:idx2]
    end
    msc = c.^2

    return (c=c, msc=msc, f=f)
end

"""
    s2_fcoherence(signal1, signal2; fs, frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing)

Calculate coherence (mean over all frequencies) and MSC (magnitude-squared coherence) between channels of `signal1` and `signal2`.

# Arguments

- `signal1::AbstractArray`
- `signal2::AbstractArray`
- `fs::Int64`
- `frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing`: return coherence only for the given frequency range

# Returns

- `c::Array{Float64, 3}`: coherence
- `msc::Array{Float64, 3}`: MSC
- `f::Vector{Float64}`: frequencies
"""
function s2_fcoherence(signal1::AbstractArray, signal2::AbstractArray; fs::Int64, frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing)

    length(signal1) == length(signal2) || throw(ArgumentError("Both signals must have the same length."))
    fs < 0 && throw(ArgumentError("fs must be > 0."))

    signal = hcat(signal1, signal2)'
    c = mt_coherence(signal, fs=fs)
    f = Vector(c.freq)
    c = c.coherence
    if frq_lim !== nothing
        frq_lim = tuple_order(frq_lim)
        frq_lim[1] < 0 && throw(ArgumentError("Lower frequency bound must be ≥ 0."))
        frq_lim[2] > fs / 2 && throw(ArgumentError("Upper frequency bound must be ≤ $fs."))
        idx1 = vsearch(frq_lim[1], f)
        idx2 = vsearch(frq_lim[2], f)
        c = c[:, :, idx1:idx2]
        f = f[idx1:idx2]
    end
    c = c[1, 2, :]
    msc = @. abs(c)^2
    
    return (c=c, msc=msc, f=f)
end

"""
    a2_l1(a1, a2)

Compare two 3-dimensional arrays `a1` and `a2` (e.g. two spectrograms), using L1 (Manhattan) distance.

# Arguments

- `a1::AbstractArray`: first array
- `a2::AbstractArray`: second array

# Returns

- `l1::Float64`
"""
function a2_l1(a1::AbstractArray, a2::AbstractArray)

    size(a1) == size(a2) || throw(ArgumentError("a1 and a2 mast have the same size."))

    l1 = sum(abs.(a1 .- a2))

    return l1
end

"""
    a2_l2(a1, a2)

Compare two 3-dimensional arrays `a1` and `a2` (e.g. two spectrograms), using L2 (Euclidean) distance.

# Arguments

- `a1::AbstractArray`: first array
- `a2::AbstractArray`: second array

# Returns

- `l2::Float64`
"""
function a2_l2(a1::AbstractArray, a2::AbstractArray)

    size(a1) == size(a2) || throw(ArgumentError("a1 and a2 mast have the same size."))

    l2 = sqrt(sum((a1 .- a2).^2))

    return l2
end

"""
    s_cums(signal)

Calculate cumulative sum of the `signal`.

# Arguments

- `signal::Vector{<:Real}`

# Returns

- `signal_cs::Vector{Float64}`
"""
function s_cums(signal::Vector{<:Real})
    
    signal_cs = cumsum(signal)

    return signal_cs
end

"""
    s_cums(signal)

Calculate cumulative sum of the `signal`.

# Arguments

- `signal::Array{<:Real, 3}`

# Returns

- `signal_cs::Array{Float64, 3}`
"""
function s_cums(signal::Array{<:Real, 3})
    
    channel_n, _, epoch_n = size(signal)
    signal_cs = similar(signal)

    @inbounds @simd for epoch_idx in 1:epoch_n
        Threads.@threads for channel_idx in 1:channel_n
            s = @view signal[channel_idx, :, epoch_idx]
            signal_cs[channel_idx, :, epoch_idx] = cumsum(s)
        end
    end

    return signal_cs
end

"""
    s_dfa(signal)

Perform detrended fluctuation analysis of the `signal`.

# Arguments

- `signal::Array{<:Real, 3}`

# Returns

- `signal_cs::Vector{Float64}`
"""
function s_dfa(signal::Array{<:Real, 3}; fs::Int64)
    
    # INCOMPLETE, DOES NOT WORK, DO NOT USE
    # also, eeg_dfa() is missing

    signal_dm = s_demean(signal)
    signal_cs = s_cums(signal_dm)

    scale_order = collect(1:20)
    scale_duration = log.(scale_order)
    scatter(scale_order, scale_duration)

    for idx in length(scale_duration):-1:1
        if scale_duration[idx] > length(signal) ÷ 2
            popat!(scale_duration, idx)
            popat!(scale_order, idx)
        end
    end
    for idx in length(scale_duration):-1:2
        if scale_duration[idx] == scale_duration[idx - 1]
            popat!(scale_duration, idx)
            popat!(scale_order, idx)
        end
    end

    epoch_rms = zeros(length(scale_duration))
    for order_idx in 1:length(scale_duration)
        epochs_n = Int(length(signal_cs) ÷ scale_duration[order_idx])
        epochs = zeros(epochs_n, scale_duration[order_idx])
        for epoch_idx in 1:epochs_n
            epochs[epoch_idx, :] = signal_cs[(((epoch_idx - 1) * scale_duration[order_idx]) + 1):(epoch_idx * scale_duration[order_idx])]
        end
        epochs_dt = s_detrend(epochs, type=:ls)
        epochs_rms = zeros(epochs_n)
        for idx in 1:size(epochs, 1)
            epochs_rms[idx] = s_rms(epochs_dt[idx, :])
        end
        epoch_rms[order_idx] = mean(epochs_rms)
    end
    scale_duration = log.(scale_duration)
    epoch_rms = log.(epoch_rms)
    scatter(scale_order, epoch_rms)
    atilde = pinv(log.(scale_order)) * log.(epoch_rms)
    plot!(atilde .* log.(scale_order))

    return signal_cs
end

"""
    s_gfp(signal)

Calculate GFP (Global Field Power) of the `signal`.

# Arguments

- `signal::Vector{<:Real}`

# Returns

- `gfp::Float64`
"""
function s_gfp(signal::Vector{<:Real})
    
    gfp = sum(signal.^2) / length(signal)

    return gfp
end

"""
    s_gfp_norm(signal)

Calculate `signal` values normalized for GFP (Global Field Power) of that signal.

# Arguments

- `signal::Vector{<:Real}`

# Returns

- `gfp_norm::Float64`
"""
function s_gfp_norm(signal::Vector{<:Real})
    
    gfp = s_gfp(signal)
    gfp_norm = signal ./ gfp

    return gfp_norm
end

"""
    s2_diss(signal1, signal2)

Calculate DISS (global dissimilarity) and spatial correlation between `signal1` and `signal2`.

# Arguments

- `signal1::Vector{<:Real}`
- `signal2::Vector{<:Real}`

# Returns

Named tuple containing:
- `diss::Float64`: global dissimilarity
- `c::Float64`: spatial correlation
"""
function s2_diss(signal1::Vector{<:Real}, signal2::Vector{<:Real})
    
    length(signal1) == length(signal2) || throw(ArgumentError("Both signals must have the same length."))
    gfp_norm1 = s_gfp_norm(signal1)
    gfp_norm2 = s_gfp_norm(signal2)
    diss = sqrt(sum((gfp_norm1 .- gfp_norm2).^2) / length(signal1))
    c = 0.5 * (2 - diss^2)

    return (diss=diss, c=c)
end