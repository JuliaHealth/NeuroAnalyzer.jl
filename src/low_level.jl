################################
#                              #
# Low-level external functions #
#                              #
################################

################################

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
function linspace(start::Union{Int64, Float64}, stop::Union{Int64, Float64}, length::Int64)

    return collect(range(start, stop, length))
end

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
function logspace(start::Union{Int64, Float64}, stop::Union{Int64, Float64}, length::Int64)

    return collect(exp10.(range(start, stop, length)))
end

"""
    m_pad0(m)

Pad the matrix `m` with zeros to make it square.

# Arguments

- `m::Union{Matrix{Int64}, Matrix{Float64}, Matrix{ComplexF64}}`

# Returns

- `m::Union{Matrix{Int64}, Matrix{Float64}, Matrix{ComplexF64}}`
"""
function m_pad0(m::Union{Matrix{Int64}, Matrix{Float64}, Matrix{ComplexF64}})

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

- `y::Union{Int64, Float64}`
- `x::Union{Vector{Int64}, Vector{Float64}}
- `return_distance::Bool`

# Returns

- `y_idx::Int64`
-` y_dist::Union{Int64, Float64}`
"""
function vsearch(y::Union{Int64, Float64}, x::Union{Vector{Int64}, Vector{Float64}}; return_distance::Bool=false)

    y_dist, y_idx = findmin(abs.(x .- y))

    return_distance == false && (return y_idx)
    return_distance == true && (return y_idx, y_dist)
end

"""
    vsearch(y, x; return_distance=false)

Return the positions of the `y` vector in the vector `x`.

# Arguments

- `x::Union{Vector{Int64}, Vector{Float64}}`
- `y::Union{Vector{Int64}, Vector{Float64}}
- `return_distance::Bool`

# Returns

- `y_idx::Int64`
- `y_dist::Union{Int64, Float64}`
"""
function vsearch(y::Union{Vector{Int64}, Vector{Float64}}, x::Union{Vector{Int64}, Vector{Float64}}; return_distance=false)

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

- `x::Union{Int64, Float64}`
- `y::Union{Int64, Float64}`

# Returns

- `phi::Float64`
- `theta::Float64`
"""
function cart2pol(x::Union{Int64, Float64}, y::Union{Int64, Float64})

    return hypot(x, y), atan(y, x)
end

"""
    pol2cart(theta, phi)

Convert polar coordinates `theta` and `phi` to cartographic.

# Arguments

- `phi::Union{Float64, Int64}`
- `theta::Union{Float64, Int64}`

# Returns

- `x::Float64`
- `y::Float64`
"""
function pol2cart(theta::Union{Float64, Int64}, phi::Union{Float64, Int64})

    return phi * cos(theta), phi * sin(theta)
end

"""
    pol2cart_sph(rho, theta, phi=0)

Convert spherical coordinates `theta` and `phi` and `rho` to cartographic.

# Arguments

- `phi::Union{Float64, Int64}`: the angle with respect to the z-axis (elevation)
- `theta::Union{Float64, Int64}`: the angle in the xy plane with respect to the x-axis (azimuth)
- `rho::Union{Float64, Int64}`: the distance from the origin to the point

# Returns

- `x::Float64`
- `y::Float64`
- `z::Float64`
"""
function sph2cart(rho::Union{Float64, Int64}, theta::Union{Float64, Int64}, phi::Union{Float64, Int64}=0)

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
    else
        throw(ArgumentError("Window type must be :hann, :bh, :bohman, :flat, :bn, :nutall."))
    end

    return w
end

"""
    hildebrand_rule(x)

Calculate Hildebrand rule for vector `x`.
If H < 0.2 then the vector `x` is symmetrical.

# Arguments

- `x::Union{Vector{Int64}, Vector{Float64}}`

# Returns

- `h::Float64`
"""
function hildebrand_rule(x::Union{Vector{Int64}, Vector{Float64}})

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
function jaccard_similarity(x::Union{Vector{Int64}, Vector{Float64}}, y::Union{Vector{Int64}, Vector{Float64}})

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
    s_rms(signal)

Calculate Root Mean Square of `signal`.

# Arguments

- `signal::Union{Vector{Int64}, Vector{Float64}}`

# Returns

- rms::Float64`
"""
function s_rms(signal::Union{Vector{Int64}, Vector{Float64}})

    return norm(signal) / sqrt(length(signal))
end

"""
    generate_sine(f, t, a, p)

Generates sine wave of `f` frequency over `t` time; optional arguments are: `a` amplitude and  `p` phase.

# Arguments

- `f::Union{Int64, Float64}`
- `t::Union{Vector{Int64}, Vector{Float64}}`
- `a::Union{Int64, Float64}`
- `p::Union{Int64, Float64}`

# Returns

- sine::Vector{Float64}`
"""
function generate_sine(f::Union{Int64, Float64}, t::Union{Vector{Int64}, Vector{Float64}}, a::Union{Int64, Float64}=1, p::Union{Int64, Float64}=0)

    return @. a * sin(2 * pi * f * t + p)
end

"""
    s_freqs(t)

Return vector of frequencies and Nyquist frequency for given time vector `t`.

# Arguments

- `t::Union{Vector{Int64}, Vector{Float64}, AbstractRange}`

# Returns

- `hz::Vector{Float64}
- `nyquist_freq::Float64`
"""
function s_freqs(t::Union{Vector{Int64}, Vector{Float64}, AbstractRange})

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
- `fs::Union{Int64, Float64}`

# Returns

- `hz::Vector{Float64}
- `nyquist_freq::Float64`
"""
function s_freqs(signal::Vector{Float64}, fs::Union{Int64, Float64})

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

- `x::Union{Vector{Int64}, Vector{Float64}}`
- `n::Int64`
- `sym::Bool=false`: if true, than pad at the beginning and at the end, otherwise only at the end.

# Returns

- `v_pad::Union{Vector{Int64}, Vector{Float64}}`
"""
function pad0(x::Union{Vector{Int64}, Vector{Float64}}, n::Int64, sym::Bool=false)

    n < 1 && throw(ArgumentError("n must be positive."))

    sym == false && (v_pad = vcat(x, zeros(eltype(x), n)))
    sym == true && (v_pad = vcat(zeros(eltype(x), n), x, zeros(eltype(x), n)))

    return v_pad
end
"""
    hz2rads(f)

Convert frequency `f` in Hz to rad/s.

# Arguments

- `f::Union{Int64, Float64}`

# Returns

- `f_rads::Float64`
"""
function hz2rads(f::Union{Int64, Float64})

    return 2 * pi * f
end

"""
    rads2hz(f)

Convert frequency `f` in rad/s to Hz.

# Arguments

- `f::Union{Int64, Float64}`

# Returns

- `f_rads::Float64`
"""
function rads2hz(f::Union{Int64, Float64})

    return f / 2 * pi
end

"""
    z_score(x)

Calculate Z-scores for each value of the vector `x`.

# Arguments

- `x::Union{Vector{Int64}, Vector{Float64}}`

# Returns

- `z_score::Vector{Float64}`
"""
function z_score(x::Union{Vector{Int64}, Vector{Float64}})

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
- `f::Union{Int64, Float64}=10.0`: frequency
- `peak::Union{Int64, Float64}=0`: sinc peak time
- `norm::Bool=true`: generate normalzied function
# Returns

- `sinc::Vector{Float64}
"""
function generate_sinc(t::AbstractRange=-2:0.01:2; f::Union{Int64, Float64}=1, peak::Union{Int64, Float64}=0, norm::Bool=true)

    norm == true && (y_sinc = @. sin(2 * pi * f * (t - peak)) / (pi * (t - peak)))
    norm == false && (y_sinc = @. sin(2 * f * (t - peak)) / (t - peak))
    nan_idx = y_sinc[y_sinc .== NaN]
    length(nan_idx) !=0 && (y_sinc[findall(isnan, y_sinc)[1]] = (y_sinc[findall(isnan, y_sinc)[1] - 1] + y_sinc[findall(isnan, y_sinc)[1] + 1]) / 2)
    
    return y_sinc
end

"""
    generate_morlet(fs, wt, wf)

Generate Morlet wavelet.

# Arguments

- `fs::Int64`: sampling rate
- `wt::Union{Int64, Float64}`: length = -wt:1/fs:wt
- `wf::Union{Int64, Float64}`: frequency
- `ncyc::Int64=5`: number of cycles
- `complex::Bool=false`: generate complex Morlet

# Returns

- `morlet::Union{Vector{Float64}, Vector{ComplexF64}}`
"""
function generate_morlet(fs::Int64, wt::Union{Int64, Float64}, wf::Union{Int64, Float64}; ncyc::Int64=5, complex::Bool=false)

    wt = -wt:1/fs:wt
    complex == false && (sin_wave = @. cos(2 * pi * wf * wt))           # for symmetry at x = 0
    complex == true && (sin_wave = @. exp(im * 2 * pi * wf * wt))       # for symmetry at x = 0
    w = 2 * (ncyc / (2 * pi * wf))^2                                    # ncyc: time-frequency precision
    g = generate_gaussian(fs, wt[end], w)
    m = sin_wave .* g

    return m
end

"""
    generate_gaussian(fs, gt, gw, pt, pa)

Generate Gaussian wave.

# Arguments

- `fs::Int64`: sampling rate
- `gt::Union{Int64, Float64}`: length = 0:1/fs:gt
- `gw::Union{Int64, Float64}=1`: width
- `pt::Union{Int64, Float64}=0`: peak time
- `pa::Union{Int64, Float64}=1`: peak amp
- 
# Returns

- `gaussian::Vector{Float64}`
"""
function generate_gaussian(fs::Int64, gt::Union{Int64, Float64}, gw::Union{Int64, Float64}=1, pt::Union{Int64, Float64}=0, pa::Union{Int64, Float64}=1.0)

    t = -gt:1/fs:gt
    g = @. pa * exp(-(t - pt)^2 / gw)

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

    return s_fft, s_sf
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
    s_spectrum(signal; pad=0)

Calculates FFT, amplitudes, powers and phases of the `signal`.

# Arguments

- `signal::AbstractArray`
- `pad::Int64`: pad the `signal` with `pad` zeros

# Returns

- `fft::Vector(ComplexF64}`
- `amplitudes::Vector{Float64}`
- `powers::Vector{Float64}`
- `phases::Vector{Float64}
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

    return s_fft, s_amplitudes, s_powers, s_phases
end

"""
    s_total_power(signal; fs)

Calculates `signal` total power.

# Arguments

- `signal::AbstractArray`
- `fs::Int64`: sampling rate

# Returns

- `stp::Float64`: signal total power
"""
function s_total_power(signal::AbstractArray; fs::Int64)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))

    psd = welch_pgram(signal, 4*fs, fs=fs)
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
- `f::Tuple`: lower and upper frequency bounds

# Returns

- `stp::Float64`: signal total power
"""
function s_band_power(signal::AbstractArray; fs::Int64, f::Tuple{Real, Real})

    psd = welch_pgram(signal, 4*fs, fs=fs)
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
- `taper::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`

# Returns

- `s_tapered::Vector{Union{Float64, ComplexF64}}`
"""
function s_taper(signal::AbstractArray; taper::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})

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
- `offset::Union{Int64, Float64}=0`: constant for :constant detrending
- `order::Int64=1`: polynomial fitting order
- `span::Float64=0.5`: smoothing of loess

# Returns
- `s_det::Vector{Float64}`
"""
function s_detrend(signal::AbstractArray; type::Symbol=:linear, offset::Union{Int64, Float64}=0, order::Int64=1, span::Float64=0.5)

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
        # do the matrix inverse for 2x2 matrix
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
    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s = @view signal[idx, :, epoch]
            s_resampled[idx, :, epoch], t_resampled = s_resample(s, t=t, new_sr=new_sr)
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
- `kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`

# Returns

- `s_conv::Union{Vector{Float64}, Vector{ComplexF64}}`
"""
function s_tconv(signal::AbstractArray; kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})

    signal = Vector(signal)
    s_conv = conv(signal, kernel)
    half_kernel = floor(Int, length(kernel) / 2)

    # remove in- and out- edges
    if mod(length(kernel), 2) == 0 
        s_conv = s_conv[half_kernel:(end - half_kernel)]
    else
        s_conv = s_conv[half_kernel:(end - half_kernel - 1)]
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
    - `:mavg`: moving average (with threshold and/or weight window)
    - `:mmed`: moving median (with threshold and/or weight window)
    - `:poly`: polynomial of `order` order
- `ftype::Symbol`: filter type:
    - `:lp`: low pass
    - `:hp`: high pass
    - `:bp`: band pass
    - `:bs`: band stop
- `cutoff::Union{Int64, Float64, Tuple}`: filter cutoff in Hz (vector for `:bp` and `:bs`)
- `order::Int64=8`: filter order
- `rp::Union{Int64, Float64}=-1`: ripple amplitude in dB in the pass band; default: 0.0025 dB for :elliptic, 2 dB for others
- `rs::Union{Int64, Float64}=-1`: ripple amplitude in dB in the stop band; default: 40 dB for :elliptic, 20 dB for others
- `dir:Symbol=:twopass`: filter direction (:onepass, :onepass_reverse, :twopass), for causal filter use :onepass
- `d::Int64=1`: window length for mean average and median average filter
- `t::Union{Int64, Float64}`: threshold for :mavg and :mmed filters; threshold = threshold * std(signal) + mean(signal) for :mavg or threshold = threshold * std(signal) + median(signal) for :mmed filter
- `window::Union{Vector{Float64}, Nothing} - window, required for FIR filter

# Returns

- `s_filtered::Vector{Float64}`
"""
function s_filter(signal::AbstractArray; fprototype::Symbol, ftype::Union{Symbol, Nothing}=nothing, cutoff::Union{Int64, Float64, Tuple}=0, fs::Int64=0, order::Int64=8, rp::Union{Int64, Float64}=-1, rs::Union{Int64, Float64}=-1, dir::Symbol=:twopass, d::Int64=1, t::Union{Int64, Float64}=0, window::Union{Vector{Float64}, Nothing}=nothing)

    fprototype in [:mavg, :mmed, :poly, :butterworth, :chebyshev1, :chebyshev2, :elliptic, :fir] || throw(ArgumentError("fprototype must be :mavg, :mmed,:butterworth, :chebyshev1, :chebyshev2, :elliptic or :fir."))
    (fprototype === :fir && (window === nothing || length(window) > length(signal))) && throw(ArgumentError("For :fir filter window must be shorter than signal."))
    (fprototype !== :mavg && fprototype !== :mmed && fprototype !== :poly) && (ftype in [:lp, :hp, :bp, :bs] || throw(ArgumentError("ftype must be :bp, :hp, :bp or :bs.")))
    (fprototype !== :mavg && fprototype !== :mmed) && (fs < 1 && throw(ArgumentError("fs must be > 0.")))
    dir in [:onepass, :onepass_reverse, :twopass] || throw(ArgumentError("direction must be :onepass, :onepass_reverse or :twopass."))
    (order < 2 && fprototype !== :poly) && (mod(order, 2) != 0 && throw(ArgumentError("order must be even and ≥ 2.")))
    (order < 1 && (fprototype !== :mavg && fprototype !== :mmed)) && throw(ArgumentError("order must be > 0."))
    d > length(signal) && throw(ArgumentError("d must be ≤ signal length."))
    
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
        if window === nothing
            s_filtered = signal
            for idx in d:-1:1
                if t > 0
                    if signal[idx] > t * std(signal) + mean(signal)
                        s_filtered[idx] = mean(signal[idx:idx+1])
                    end
                else
                    s_filtered[idx] = mean(signal[idx:idx+1])
                end
            end
            for idx in (1 + d):(length(signal) - d)
                if t > 0
                    if signal[idx] > t * std(signal) + mean(signal)
                        s_filtered[idx] = mean(signal[(idx - d):(idx + d)])
                    end
                else
                    s_filtered[idx] = mean(signal[(idx - d):(idx + d)])
                end
            end
            for idx in (length(signal) - d + 1):length(signal)
                if t > 0
                    if signal[idx] > t * std(signal) + mean(signal)
                        s_filtered[idx] = mean(signal[idx-1:idx])
                    end
                else
                    s_filtered[idx] = mean(signal[idx-1:idx])
                end
            end
        else
            s_filtered = s_tconv(signal, kernel=window)
        end

        return s_filtered
    end

    if fprototype === :mmed
        s_filtered = signal
        for idx in d:-1:1
            if t > 0
                if signal[idx] > t * std(signal) + median(signal)
                    s_filtered[idx] = median(signal[idx:idx+1])
                end
            else
                s_filtered[idx] = median(signal[idx:idx+1])
            end
        end
        for idx in (1 + d):(length(signal) - d)
            if t > 0
                if signal[idx] > t * std(signal) + median(signal)
                    s_filtered[idx] = median(signal[(idx - d):(idx + d)])
                end
            else
                s_filtered[idx] = median(signal[(idx - d):(idx + d)])
            end
        end
        for idx in (length(signal) - d + 1):length(signal)
            if t > 0
                if signal[idx] > t * std(signal) + median(signal)
                    s_filtered[idx] = median(signal[idx-1:idx])
                end
            else
                s_filtered[idx] = median(signal[idx-1:idx])
            end
        end

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
        if window === nothing
            @warn "Using default window for :fir filter: hanning($(3 * floor(Int64, fs / cutoff[1])))."
            window = hanning(3 * floor(Int64, fs / cutoff[1]))
        end
        if ftype === :hp || ftype === :bp || ftype === :bs
            mod(length(window), 2) == 0 && (window = vcat(window[1:((length(window) ÷ 2) - 1)], window[((length(window) ÷ 2) + 1):end]))
        end
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

    eeg_filter = digitalfilter(responsetype, prototype)

    dir === :twopass && (s_filtered = filtfilt(eeg_filter, signal))
    dir === :onepass && (s_filtered = filt(eeg_filter, signal))
    dir === :onepass_reverse && (s_filtered = filt(eeg_filter, reverse(signal)))

    return s_filtered
end

"""
    s_psd(signal; fs, norm=false)

Calculate power spectrum density of the `signal`.

# Arguments
- `signal::AbstractArray`
- `fs::Int64`: sampling rate
- `norm::Bool`: normalize do dB

# Returns

- `psd_pow::Vector{Float64}`
- `psd_frq::Vector{Float64}`
"""
function s_psd(signal::AbstractArray; fs::Int64, norm::Bool=false)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))
    
    psd = welch_pgram(signal, 4*fs, fs=fs)
    psd_pow = power(psd)
    psd_frq = freq(psd)
    norm == true && (psd_pow = pow2db.(psd_pow))

    return psd_pow, Vector(psd_frq)
end

"""
    s_psd(signal; fs, norm=false)

Calculate power spectrum density of the `signal`.

# Arguments
- `signal::Matrix{Float64}`
- `fs::Int64`: sampling rate
- `norm::Bool`: normalize do dB

# Returns

- `psd_pow::Matrix{Float64}`
- `psd_frq::Matrix{Float64}`
"""
function s_psd(signal::Matrix{Float64}; fs::Int64, norm::Bool=false)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))

    channel_n = size(signal, 1)
    psd_tmp, frq_tmp = s_psd(signal[1, :], fs=fs, norm=norm)
    psd_pow = zeros(channel_n, length(psd_tmp))
    psd_frq = zeros(channel_n, length(frq_tmp))

    @inbounds @simd for idx in 1:channel_n
        s = @view signal[idx, :]
        psd_pow[idx, :], psd_frq[idx, :] = s_psd(s, fs=fs, norm=norm)
    end
    
    return psd_pow, psd_frq
end

"""
    s_psd(signal; fs, norm=false)

Calculate power spectrum density of the `signal`.

# Arguments
- `signal::Array{Float64, 3}`
- `fs::Int64`: sampling rate
- `norm::Bool`: normalize do dB

# Returns

- `psd_pow::Array{Float64, 3}`
- `psd_frq::Array{Float64, 3}`
"""
function s_psd(signal::Array{Float64, 3}; fs::Int64, norm::Bool=false)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))

    channel_n = size(signal, 1)
    epoch_n = size(signal, 3)
    psd_tmp, frq_tmp = s_psd(signal[1, :, 1], fs=fs, norm=norm)
    psd_pow = zeros(channel_n, length(psd_tmp), epoch_n)
    psd_frq = zeros(channel_n, length(frq_tmp), epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            s = @view signal[idx, :, epoch]
            psd_pow[idx, :, epoch], psd_frq[idx, :, epoch] = s_psd(s, fs=fs, norm=norm)
        end
    end
    
    return psd_pow, psd_frq
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
    s2_coherence(signal1, signal2)

Calculate coherence between `signal1` and `signal2`.

# Arguments

- `signal1::AbstractArray`
- `signal2::AbstractArray`

# Returns

- `coherence::Vector{ComplexF64}`
"""
function s2_coherence(signal1::AbstractArray, signal2::AbstractArray)

    length(signal1) == length(signal2) || throw(ArgumentError("Both signals must have the same length."))

    s1_fft = fft(signal1) ./ length(signal1)
    s2_fft = fft(signal2) ./ length(signal2)

    coherence = (abs.((s1_fft) .* conj.(s2_fft)).^2) ./ (s1_fft .* s2_fft)

    return coherence
end

"""
    s_pca(signal, n)

Calculates `n` first PCs for `signal`.

# Arguments

- `signal::Array{Float64, 3}`
- `n::Int64`: number of PCs

# Returns

- `pc::Array{Float64, 3}:`: PC(1)..PC(n) × epoch
- `pc_var::Matrix{Float64}`: PC_VAR(1)..PC_VAR(n) × epoch
"""
function s_pca(signal::Array{Float64, 3}; n::Int64)

    n < 0 && throw(ArgumentError("n must be ≥ 1."))
    n > size(signal, 1) && throw(ArgumentError("Number of PCs must be ≤ $(size(signal, 1))."))

    channel_n, _, epoch_n = size(signal)
    pc = zeros(n, size(signal, 2), epoch_n)
    pc_var = zeros(n, epoch_n)
    pc_reconstructed = zeros(size(signal))
    M = []
    
    Threads.@threads for epoch in 1:epoch_n
        s = @view signal[:, :, epoch]
        # m_cov = s_cov(s)
        # eig_val, eig_vec = eigen(m_cov)
        # eig_val_idx = sortperm(eig_val, rev=true)
        # eig_val = eig_val[eig_val_idx]
        # eig_vec = m_sort(eig_vec, eig_val_idx)
        # eig_val = 100 .* eig_val / sum(eig_val) # convert to %

        M = MultivariateStats.fit(PCA, s, maxoutdim=n)
        v = principalvars(M) ./ var(M) * 100

        for idx in 1:n
            pc_var[idx, epoch] = v[idx]
            # pc[idx, :, epoch] = (eig_vec[:, idx] .* s)[idx, :]
            pc[idx, :, epoch] = MultivariateStats.predict(M, s)[idx, :]
        end
    end

    return pc, pc_var, M
end

"""
    s_fconv(signal; kernel)

Perform convolution in the frequency domain between `signal` and `kernel`.

# Arguments

- `signal::AbstractArray`
- `kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}}`

# Returns

- `s_conv::Vector{ComplexF64}`
"""
function s_fconv(signal::AbstractArray; kernel::Union{Vector{Int64}, Vector{Float64}, Vector{ComplexF64}})

    n_signal = length(signal)
    n_kernel = length(kernel)
    n_conv = n_signal + n_kernel - 1
    half_kernel = floor(Int64, n_kernel / 2)
    s_fft = fft0(signal, n_conv)
    kernel_fft = fft0(kernel, n_conv)
    s_conv = ifft(s_fft .* kernel_fft)
    
    # remove in- and out- edges
    if mod(n_kernel, 2) == 0 
        s_conv = s_conv[half_kernel:(end - half_kernel)]
    else
        s_conv = s_conv[half_kernel:(end - half_kernel - 1)]
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

    @inbounds @simd for epoch in 1:epoch_n
        s = @view signal[:, :, epoch]

        f === :tanh && (M = MultivariateStats.fit(ICA, s, n, tol=tol, maxiter=iter, fun=MultivariateStats.Tanh(1.0)))
        f === :gaus && (M = MultivariateStats.fit(ICA, s, n, tol=tol, maxiter=iter, fun=MultivariateStats.Gaus()))

        n == size(signal, 1) && (mw = inv(M.W)')
        n < size(signal, 1) && (mw = pinv(M.W)')

        for idx in 1:n
            ic[idx, :, epoch] = MultivariateStats.predict(M, s)[idx, :]
        end

        ic_mw[:, :, epoch] = mw
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

    @inbounds @simd for epoch in 1:epoch_n
        s_reconstructed[:, :, epoch] = ic_mw[:, ic_removal, epoch] * ic[ic_removal, :, epoch]
    end

    return s_reconstructed
end

"""
    s_spectrogram(signal; fs, norm=true, demean=true)

Calculate spectrogram of `signal`.

# Arguments

- `signal::AbstractArray`
- `fs::Int64`: sampling frequency
- `norm::Bool`: normalize powers to dB
- `demean::Bool`: demean signal prior to analysis

# Returns

- `s_pow::Matrix{Float64}`: powers
- `s_frq::Vector{Float64}`: frequencies
- `s_t::Vector{Float64}`: time
"""
function s_spectrogram(signal::AbstractArray; fs::Int64, norm::Bool=true, demean::Bool=true)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1 Hz."))

    demean == true && (signal = s_demean(signal))

    nfft = length(signal)
    interval = fs
    overlap = round(Int64, fs * 0.85)

    spec = spectrogram(signal, interval, overlap, nfft=nfft, fs=fs, window=hanning)
    s_t = collect(spec.time)
    s_frq = Vector(spec.freq)
    if norm == true
        s_pow = pow2db.(spec.power)
    else
        s_pow = Vector(spec.power)
    end

    return s_pow, s_frq, s_t
end

"""
    s_detect_epoch_flat(signal::Array{Float64, 3}, threshold=0.1)

Detect bad `signal` epochs based on: flat channel(s)

# Arguments

- `signal::Array{Float64, 3}`

# Returns

- `bad_epochs_score::Vector{Int64}`: percentage of bad channels per epoch
"""
function s_detect_epoch_flat(signal::Array{Float64, 3})
    
    channel_n, _, epoch_n = size(signal)

    bad_epochs_score = zeros(epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        Threads.@threads for idx in 1:channel_n
            idx=1
            epoch=1
            c_tmp_diff = abs.(diff(signal[idx, :, epoch]))
            # add tolerance around zero μV
            sum(c_tmp_diff) < eps() && (bad_epochs_score[epoch] += 1)
        end
    end

    bad_epochs_score = round.(bad_epochs_score ./ channel_n, digits=1)

    return bad_epochs_score
end

"""
    s_detect_epoch_rmse(signal::Array{Float64, 3})

Detect bad `signal` epochs based on: RMSE vs average channel > 95%CI.

# Arguments

- `signal::Array{Float64, 3}`

# Returns

- `bad_epochs_score::Vector{Int64}`: percentage of bad channels per epoch
"""
function s_detect_epoch_rmse(signal::Array{Float64, 3})
    
    channel_n, _, epoch_n = size(signal)

    bad_epochs_score = zeros(epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        ch_m = vec(median(signal[:, :, epoch], dims=1))
        rmse_ch = zeros(channel_n)
        Threads.@threads for idx in 1:channel_n
            rmse_ch[idx] = s2_rmse(signal[idx, :, epoch], ch_m)
        end
        Threads.@threads for idx in 1:channel_n
            rmse_ch[idx] > HypothesisTests.confint(OneSampleTTest(rmse_ch))[2] && (bad_epochs_score[epoch] += 1)
        end
    end

    bad_epochs_score = round.(bad_epochs_score ./ channel_n, digits=1)

    return bad_epochs_score
end

"""
    detect_epoch_rmsd(signal::Array{Float64, 3})

Detect bad `signal` epochs based on: RMSD vs average channel > 95%CI.
# Arguments

- `signal::Array{Float64, 3}`

# Returns

- `bad_epochs_score::Vector{Int64}`: percentage of bad channels per epoch
"""
function s_detect_epoch_rmsd(signal::Array{Float64, 3})
    
    channel_n, _, epoch_n = size(signal)

    bad_epochs_score = zeros(epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        ch_m = median(signal[:, :, epoch], dims=1)
        rmsd_ch = zeros(channel_n)
        Threads.@threads for idx in 1:channel_n
            rmsd_ch[idx] = Distances.rmsd(signal[idx, :, epoch], ch_m)
        end
        Threads.@threads for idx in 1:channel_n
            rmsd_ch[idx] > HypothesisTests.confint(OneSampleTTest(rmsd_ch))[2] && (bad_epochs_score[epoch] += 1)
        end
    end

    bad_epochs_score = round.(bad_epochs_score ./ channel_n, digits=1)

    return bad_epochs_score
end

"""
    s_detect_epoch_euclid(signal::Array{Float64, 3})

Detect bad `signal` epochs based on: Euclidean distance vs median channel > 95% CI.

# Arguments

- `signal::Array{Float64, 3}`

# Returns

- `bad_epochs_score::Vector{Int64}`: percentage of bad channels per epoch
"""
function s_detect_epoch_euclid(signal::Array{Float64, 3})
    
    channel_n, _, epoch_n = size(signal)

    bad_epochs_score = zeros(epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        ch_m = median(signal[:, :, epoch], dims=1)
        ed_ch = zeros(channel_n)
        Threads.@threads for idx in 1:channel_n
            ed_ch[idx] = euclidean(signal[idx, :, epoch], ch_m)
        end
        Threads.@threads for idx in 1:channel_n
            ed_ch[idx] > HypothesisTests.confint(OneSampleTTest(ed_ch))[2] && (bad_epochs_score[epoch] += 1)
        end
    end

    bad_epochs_score = round.(bad_epochs_score ./ channel_n, digits=1)

    return bad_epochs_score
end

"""
    s_detect_epoch_p2p(signal::Array{Float64, 3})

Detect bad `signal` epochs based on: p2p amplitude > upper 95% CI p2p amplitude.

# Arguments

- `signal::Array{Float64, 3}`

# Returns

- `bad_epochs_score::Vector{Int64}`: percentage of bad channels per epoch
"""
function s_detect_epoch_p2p(signal::Array{Float64, 3})
    
    channel_n, _, epoch_n = size(signal)

    bad_epochs_score = zeros(epoch_n)

    @inbounds @simd for epoch in 1:epoch_n
        p2p = zeros(channel_n)
        Threads.@threads for idx in 1:channel_n
            p2p[idx] = maximum(signal[idx, :, epoch]) + abs(minimum(signal[idx, :, epoch]))
        end
        Threads.@threads for idx in 1:channel_n
            p2p[idx] > HypothesisTests.confint(OneSampleTTest(p2p))[2] && (bad_epochs_score[epoch] += 1)
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