export generate_window
export generate_sine
export generate_cosine
export generate_csine
export generate_sinc
export generate_morlet
export generate_gaussian
export generate_noise
export generate_signal
export generate_morlet_fwhm
export generate_square
export generate_triangle

"""
    generate_window(type, n; <keyword arguments>)

Return an `n`-point symmetric window of the given type.

# Arguments

- `type::Symbol`: window type:
    - `:hann`: Hann
    - `:bh`: Blackman-Harris
    - `:bohman`: Bohman
    - `:flat`: Flat-top
    - `:bn`: Blackman-Nuttall
    - `:nutall`: Nuttall
    - `:triangle`: symmetric triangle (left half ↑, right half ↓)
    - `:exp`: symmetric exponential (left half ↑, right half ↓)
- `n::Int64`: window length in samples; must be ≥ 1
- `even::Bool=false`: if `true` and `n` is odd, increment `n` by 1 to enforce even length

# Returns

- `w::Vector{Float64}`: generated window of length `n` (or `n + 1` if `even=true` and `n` was odd)

# Throws
- `ArgumentError`: if `n < 1` or `type` is not a recognised symbol

# See also

[`generate_sine`](@ref), [`generate_gaussian`](@ref)
"""
function generate_window(type::Symbol, n::Int64; even::Bool = false)::Vector{Float64}

    _check_var(type, [:hann, :bh, :bohman, :flat, :bn, :nutall, :triangle, :exp], "type")
    @assert n >= 1 "n must be ≥ 1."

    even && mod(n, 2) != 0 && (n += 1)
    t = range(0, 1, n)

    if type === :hann

        # Standard Hann window: 0.5 × (1 − cos(2πt))
        w = @. 0.5 * (1 - cos(2 * pi * t))

    elseif type === :bh

        # Blackman-Harris 4-term window
        w = @. 0.35875 - 0.48829 * cos(2 * pi * t) + 0.14128 * cos(4 * pi * t) -
            0.01168 * cos(6 * pi * t)

    elseif type === :bohman

        # Bohman window (product of triangle and sinc)
        w = @. (1 - abs(t * 2 - 1)) * cos(pi * abs(t * 2 - 1)) +
            (1 / pi) * sin(pi * abs(t * 2 - 1))

    elseif type === :flat

        # Flat-top 5-term window (minimises amplitude error)
        w = @. 0.21557 - 0.41663 * cos(2 * pi * t) + 0.27726 * cos(4 * pi * t) -
            0.08357 * cos(6 * pi * t) + 0.00694 * cos(8 * pi * t)

    elseif type === :bn

        # Blackman-Nuttall 4-term window
        w = @. 0.3635819 - 0.4891775 * cos(2 * pi * t) + 0.1365995 * cos(4 * pi * t) -
            0.0106411 * cos(6 * pi * t)

    elseif type === :nutall

        # Nuttall 4-term window
        w = @. 0.355768 - 0.487396 * cos(2 * pi * t) + 0.144232 * cos(4 * pi * t) -
            0.012604 * cos(6 * pi * t)

    elseif type === :triangle

        # symmetric triangle: ramp up to the midpoint, then ramp down
        w = zeros(n)
        @inbounds for idx in 1:((n ÷ 2) + 1)
            w[idx] = @. (idx * (idx + 1)) / 2
        end
        w[((n ÷ 2) + 2):n] = reverse(w)[((n ÷ 2) + 2):n]
        w .= w ./ maximum(w)

    elseif type === :exp

        # symmetric exponential: decaying from center outwards
        w = ones(n)
        if mod(n, 2) == 0
            @inbounds for idx in 1:(n ÷ 2)
                w[idx] = 1 / idx
            end
            w[1:((n ÷ 2))] = reverse(w[1:((n ÷ 2))])
            w[((n ÷ 2) + 1):n] = reverse(w[1:(n ÷ 2)])
        else
            @inbounds for idx in 1:((n ÷ 2) + 1)
                w[idx] = 1 / idx
            end
            w[1:((n ÷ 2) + 1)] = reverse(w[1:((n ÷ 2) + 1)])
            w[((n ÷ 2) + 2):n] = reverse(w[1:(n ÷ 2)])
        end

    end

    return w

end

"""
    generate_sine(f, t, a, p)

Generate a sine wave.

Computes `a × sin(2πft + φ)` where `φ = deg2rad(p)`.

# Arguments

- `f::Real`: frequency in Hz
- `t::AbstractVector`: time vector in seconds
- `a::Real`: amplitude; the signal ranges over `[-a, +a]`
- `p::Real`: phase shift in degrees

# Returns

- `Vector{Float64}`: sine wave sampled at the points in `t`

# See also

[`generate_cosine`](@ref), [`generate_csine`](@ref)
"""
function generate_sine(
    f::Real,
    t::AbstractVector,
    a::Real = 1,
    p::Real = 0
)::Vector{Float64}

    return @. a * sin(2 * pi * f * t + deg2rad(p))

end

"""
    generate_cosine(f, t, a, p)

Generate a cosine wave.

Computes `a × cos(2πft + φ)` where `φ = deg2rad(p)`.

# Arguments

- `f::Real`: frequency in Hz
- `t::AbstractVector`: time vector in seconds
- `a::Real`: amplitude; the signal ranges over `[-a, +a]`
- `p::Real`: phase shift in degrees

# Returns

- `Vector{Float64}`: cosine wave sampled at the points in `t`

# See also

[`generate_sine`](@ref), [`generate_csine`](@ref)
"""
function generate_cosine(
    f::Real,
    t::AbstractVector,
    a::Real = 1,
    p::Real = 0
)::Vector{Float64}

    return @. a * cos(2 * pi * f * t + deg2rad(p))

end

"""
    generate_csine(f, t, a)

Generate a complex sine wave (complex exponential).

Computes `a × exp(i × 2πft)`.

# Arguments

- `f::Real`: frequency in Hz
- `t::AbstractVector`: time vector in seconds
- `a::Real`: amplitude; the signal ranges over `[-a, +a]`

# Returns

- cs::Vector{ComplexF64}`: complex exponential sampled at the points in `t`

# See also

[`generate_sine`](@ref), [`generate_cosine`](@ref)
"""
function generate_csine(f::Real, t::AbstractVector, a::Real = 1)::Vector{ComplexF64}

    return @. a * exp(1im * 2 * pi * f * t)

end

"""
    generate_sinc(t, f, peak, norm)

Generate a sinc function.

Produces either the normalized (`sin(2πf(t−peak)) / (π(t−peak))`) or unnormalized (`sin(2f(t−peak)) / (t−peak)`) sinc. The singularity at `t = peak` is resolved by linear interpolation from the two neighboring samples.

# Arguments

- `t::AbstractVector=-2:0.01:2`: time vector
- `f::Real=1.0: frequency in Hz
- `peak::Real=0`: time of the sinc peak (location of the singularity)
- `norm::Bool=true`: if `true`, generate the normalized sinc; otherwise generate the unnormalized sinc

# Returns

- `Vector{Float64}`: sinc function sampled at the points in `t`

# See also

[`generate_sine`](@ref)
"""
function generate_sinc(
    t::AbstractVector = -2:0.01:2;
    f::Real = 1.0,
    peak::Real = 0,
    norm::Bool = true
)::Vector{Float64}

    s = if norm
        (@. sin(2 * pi * f * (t - peak)) / (pi * (t - peak)))
    else
        (@. sin(2 * f * (t - peak)) / (t - peak))
    end

    # resolve the singularity at t == peak by linear interpolation
    nan_idx = findfirst(isnan, s)
    if nan_idx !== nothing
        s[nan_idx] = (s[nan_idx - 1] + s[nan_idx + 1]) / 2
    end

    return s

end

"""
    generate_morlet(fs, f, t; <keyword arguments>)

Generate a Morlet wavelet.

The wavelet is the product of a complex (or real) sine wave at frequency `f` and a Gaussian envelope whose width is determined by `ncyc`. The time axis spans `−t : 1/fs : t`.

# Arguments

- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `f::Real`: centre frequency in Hz
- `t::Real=1`: half-length of the wavelet in seconds; must be > 0
- `ncyc::Int64=5`: number of cycles; controls the Gaussian width; must be ≥ 1
- `complex::Bool=false`: if `true`, return a complex Morlet; otherwise return the real part only

# Returns

- `Union{Vector{Float64}, Vector{ComplexF64}}`: Morlet wavelet of length `length(-t:1/fs:t)`

# Throws

- `ArgumentError`: if `fs < 1`, `ncyc < 1`, or `t ≤ 0`.

# See also

[`generate_morlet_fwhm`](@ref), [`generate_gaussian`](@ref)
"""
function generate_morlet(
    fs::Int64,
    f::Real,
    t::Real = 1;
    ncyc::Int64 = 5,
    complex::Bool = false
)::Union{Vector{Float64}, Vector{ComplexF64}}

    @assert fs >= 1 "fs must be ≥ 1."
    @assert ncyc >= 1 "ncyc must be ≥ 1."
    @assert t > 0 "t must be > 0."

    tvec = (-t):(1 / fs):t
    sin_wave = @. exp(im * 2 * pi * f * tvec)
    g = generate_gaussian(fs, f, t; ncyc=ncyc)
    
    return @. complex ? sin_wave * g : real(sin_wave) * g

end

"""
    generate_gaussian(fs, f, t; <keyword arguments>)

Generate a Gaussian envelope.

The standard deviation of the Gaussian is `σ = ncyc / (2πf)`. The time axis spans `−t : 1/fs : t`.

# Arguments

- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `f::Real`: frequency in Hz; determines the Gaussian width via `σ = ncyc / (2πf)`
- `t::Real=1`: half-length of the time axis in seconds; must be > 0
- `ncyc::Int64`: : number of cycles; controls the width (SD) of the Gaussian; must be ≥ 1
- `a::Real=1`: peak amplitude

# Returns

- `Vector{Float64}`: Gaussian envelope of length `length(-t:1/fs:t)`

# Throws
- `ArgumentError`: if `fs < 1`, `ncyc < 1`, or `t ≤ 0`

# See also

[`generate_morlet`](@ref), [`generate_morlet_fwhm`](@ref)
"""
function generate_gaussian(
    fs::Int64,
    f::Real,
    t::Real = 1;
    ncyc::Int64 = 5,
    a::Real = 1.0
)::Vector{Float64}

    @assert fs >= 1 "fs must be ≥ 1."
    @assert ncyc >= 1 "ncyc must be ≥ 1."
    @assert t > 0 "t must be > 0."

    tvec = (-t):(1 / fs):t
    # Gaussian SD: wider for more cycles
    s = ncyc / (2 * pi * f)
    return @. a * exp(-0.5 * (tvec / s)^2)

end

"""
    generate_noise(n, a; <keyword arguments>)

Generate a noise signal.

The raw noise is normalized to `[−1, 1]` and then scaled by `a`.

# Arguments

- `n::Int64`: number of samples; must be ≥ 1
- `a::Real=1.0`: amplitude; the output ranges over `[-a, +a]`
- `type::Symbol=:whiten`: noise type:
    - `:whiten`: white noise (normally distributed)
    - `:whiteu`: white noise (uniformly distributed)
    - `:pink`: pink (1/f²) noise.

# Returns

- `Vector{Float64}`: noise signal of length `n`

# Throws
- `ArgumentError`: if `n < 1` or `type` is not a recognised symbol

# See also

[`generate_signal`](@ref)
"""
function generate_noise(n::Int64, a::Real = 1.0; type::Symbol = :whiten)::Vector{Float64}

    _check_var(type, [:whiten, :whiteu, :pink], "type")
    @assert n >= 1 "n must be ≥ 1."

    if type === :whiten
        s = randn(n)
    elseif type === :whiteu
        s = rand(n)
    elseif type === :pink
        # 1/f² shaping via spectral coloring
        noise = randn(n)
        spec = fft(noise)
        # build a 1/f² shaping filter over the full spectrum length
        shaping = range(1.0, 2.0, length=length(spec)) .^ 2
        s = real(ifft(spec .* shaping))
    end

    return normalize_minmax(s) .* a

end

"""
    generate_signal(n, a)

Generate a random walk signal based on cumulative normally distributed noise.

The cumulative sum introduces temporal autocorrelation, producing a Brownian-motion–like signal. The result is normalised to `[−1, 1]` and scaled by `a`.

# Arguments

- `n::Int64`: number of samples; must be ≥ 1
- `a::Real=1.0`: amplitude; the output ranges over `[-a, +a]`

# Returns

- `Vector{Float64}`: random walk signal of length `n`

# Throws

- `ArgumentError`: Iif `n < 1`.

# See also

[`generate_noise`](@ref)
"""
function generate_signal(n::Int64, a::Real = 1.0)::Vector{Float64}

    @assert n >= 1 "n must be ≥ 1."

    s = cumsum(randn(n))
    return normalize_minmax(s) .* a

end

"""
    generate_morlet_fwhm(fs, f, t; <keyword arguments>)

Generate a complex Morlet wavelet parameterized by full-width at half-maximum (FWHM).

Uses the FWHM-based Gaussian envelope `exp(−4 ln 2 × t² / h²)` instead of a cycle-count parameter, providing more intuitive control over time–frequency resolution.

# Arguments

- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `f::Real`: centre frequency in Hz
- `t::Real=1`: half-length of the wavelet in seconds; must be > 0
- `h::Float64=0.25`: full-width at half-maximum of the Gaussian envelope in seconds; must be > 0

# Returns

- `Vector{ComplexF64}`: complex Morlet wavelet of length `length(-t:1/fs:t)`

# Throws

- `ArgumentError`: if `fs < 1`, `t ≤ 0`, or `h ≤ 0`

# Reference

Cohen MX. A better way to define and describe Morlet wavelets for time-frequency analysis. NeuroImage. 2019 Oct;199:81–6.

# See also

[`generate_morlet`](@ref), [`generate_gaussian`](@ref)
"""
function generate_morlet_fwhm(
    fs::Int64,
    f::Real,
    t::Real = 1;
    h::Float64 = 0.25
)::Vector{ComplexF64}

    @assert fs >= 1 "fs must be ≥ 1."
    @assert t > 0 "t must be > 0."
    @assert h > 0 "h must be > 0."

    tvec = (-t):(1 / fs):t
    return @. exp(2 * 1im * π * f * tvec) * exp((-4 * log(2) * tvec^2) / h^2)

end

"""
    generate_square(t, a, p, w, offset)

Generate a square wave.

# Arguments

- `t::AbstractVector`: time vector
- `a::Real=1`: amplitude scaling factor
- `p::Real`: phase (horizontal shift of the wave) = duty cycle
- `w::Real`: width = duty-cycle threshold; samples where `mod(p + t, 2) > w` are high
- `offset::Real=0`: vertical amplitude offset

# Returns

- `Vector{Float64}`: square wave sampled at the points in `t`

# See also

[`generate_triangle`](@ref)
"""
function generate_square(
    t::AbstractVector,
    a::Real = 1,
    p::Real,
    w::Real = 1,
    offset::Real = 0
)::Vector{Float64}

    return @. offset + a * (mod(p + t, 2) > w)

end

"""
    generate_triangle(t, a)

Generate a triangle wave.

Computes `a × |mod(t, 2) − 1|`, which produces a symmetric triangle wave with period 2 and amplitude `a`.

# Arguments

- `t::AbstractVector`: time vector
- `a::Real`: amplitude scaling factor

# Returns

- `Vector{Float64}`: triangle wave sampled at the points in `t`

# See also

[`generate_square`](@ref)
"""
function generate_triangle(t::AbstractVector, a::Real = 1)::Vector{Float64}

    return @. a * abs(mod(t, 2) - 1)

end
