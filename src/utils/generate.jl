export generate_window
export generate_sine
export generate_csine
export generate_sinc
export generate_morlet
export generate_gaussian
export generate_noise
export generate_morlet_fwhm
export generate_square
export generate_triangle

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

    _check_var(type, [:hann, :bh, :bohman, :flat, :bn, :nutall, :triangle, :exp], "type")
    @assert n >= 1 "n must be ≥ 1."

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
        w = zeros(n)
        @inbounds @simd for idx in 1:((n ÷ 2) + 1)
            w[idx] = @. (idx * (idx + 1)) / 2
        end
        w[((n ÷ 2) + 2):n] = reverse(w)[((n ÷ 2) + 2):n]
        w .= w ./ maximum(w)
    elseif type === :exp
        w = ones(n)
        if mod(n, 2) == 0
            @inbounds @simd for idx in 1:(n ÷ 2)
                w[idx] = 1 / idx
            end
            w[1:((n ÷ 2))] = reverse(w[1:((n ÷ 2))])
            w[((n ÷ 2) + 1):n] = reverse(w[1:(n ÷ 2)])
        else
            @inbounds @simd for idx in 1:((n ÷ 2) + 1)
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

Generates sine wave.

# Arguments

- `f::Real`: frequency [Hz]
- `t::Union{AbstractVector, AbstractRange}`: time vector
- `a::Real`: amplitude
- `p::Real`: initial phase

# Returns

- s::Vector{Float64}`
"""
function generate_sine(f::Real, t::Union{AbstractVector, AbstractRange}, a::Real=1, p::Real=0)
    
    s = @. a * sin(2 * pi * f * t + p)
    
    return s

end

"""
    generate_csine(f, t, a)

Generates complex sine wave.

# Arguments

- `f::Real`: frequency [Hz]
- `t::Union{AbstractVector, AbstractRange}`: time vector
- `a::Real`: amplitude

# Returns

- cs::Vector{ComplexF64}`
"""
function generate_csine(f::Real, t::Union{AbstractVector, AbstractRange}, a::Real=1)

    cs = @. a * exp(1im * 2 * pi * f * t)
    
    return cs

end

"""
    generate_sinc(t, f, peak, norm)

Generate sinc function.

# Arguments

- `t::AbstractRange=-2:0.01:2`: time
- `f::Real=10.0`: frequency
- `peak::Real=0`: sinc peak time
- `norm::Bool=true`: generate normalized function

# Returns

- `s::Vector{Float64}`
"""
function generate_sinc(t::AbstractRange=-2:0.01:2; f::Real=1, peak::Real=0, norm::Bool=true)

    s = norm == true ? (@. sin(2 * pi * f * (t - peak)) / (pi * (t - peak))) : (@. sin(2 * f * (t - peak)) / (t - peak))
    nan_idx = isnan.(s)
    sum(nan_idx) != 0 && (s[findall(isnan, s)[1]] = (s[findall(isnan, s)[1] - 1] + s[findall(isnan, s)[1] + 1]) / 2)
    
    return s

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

    @assert fs >= 1 "fs must be ≥ 1."
    @assert ncyc >= 1 "ncyc must be ≥ 1."
    @assert t > 0 "t must be > 0."

    t = -t:1/fs:t
    sin_wave = complex == true ? (@. exp(im * 2 * pi * f * t)) : (@. sin(2 * pi * f * t))
    g = generate_gaussian(fs, f, t[end], ncyc=ncyc)
    mw = sin_wave .* g

    return mw

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

- `g::Vector{Float64}`
"""
function generate_gaussian(fs::Int64, f::Real, t::Real=1; ncyc::Int64=5, a::Real=1.0)

    @assert fs >= 1 "fs must be ≥ 1."
    @assert ncyc >= 1 "ncyc must be ≥ 1."
    @assert t > 0 "t must be > 0."

    t = -t:1/fs:t
    s = ncyc / (2 * pi * f)             # Gaussian width (standard deviation)
    g = @. a * exp(-(t/s)^2 / 2)        # Gaussian

    return g

end

"""
    generate_noise(n, amp; type)

Generate noise.

# Arguments

- `n::Int64`: length (in samples)
- `amp::Real=1.0`: amplitude, signal amplitude will be in `[-amp, +amp]`
- `type::Symbol=:whiten`: noise type:
    - `:whiten`: normal distributed
    - `:whiteu`: uniformly distributed
    - `:pink`

# Returns

- `s::Float64`
"""
function generate_noise(n::Int64, amp::Real=1.0; type::Symbol=:whiten)

    _check_var(type, [:whiten, :whiteu, :pink], "type")

    if type === :whiten
        s = randn(n)
    elseif type === :whiteu
        s = rand(n)
    elseif type === :pink
        s = real(ifft(fft(randn(n)) .* linspace(-1, 1, length(fft(randn(n)))).^2)) .* 2
    end

    s = normalize_minmax(s)
    s .*= amp

    return s

end

"""
    generate_morlet_fwhm(fs, f, t; h)

Generate Morlet wavelet using Mike X Cohen formula.

# Arguments

- `fs::Int64`: sampling rate
- `f::Real`: frequency
- `t::Real=1`: length = -t:1/fs:t
- `h::Float64=0.25`: full width at half-maximum in seconds (FWHM)

# Returns

- `mw::Vector{ComplexF64}`

# Source

Cohen MX. A better way to define and describe Morlet wavelets for time-frequency analysis. NeuroImage. 2019 Oct;199:81–6. 
"""
function generate_morlet_fwhm(fs::Int64, f::Real, t::Real=1; h::Float64=0.25)
    
    t = -t:1/fs:t
    mw = @. exp(2 * 1im * π * f * t) * exp((-4 * log(2) * t^2) / (h^2))
    
    return mw
    
end

"""
    generate_square(t, a, p, w, offset)

Generates square wave.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: time vector
- `a::Real`: amplitude
- `p::Real`: duty cycle
- `w::Real`: width
- `offset::Real`: amplitude offset

# Returns

- s::Vector{Float64}`
"""
function generate_square(t::Union{AbstractVector, AbstractRange}, p::Real, a::Real=1, offset::Real=0, w::Real=1)

    s = @. offset + a * (mod(p + t, 2) > w)
    return s

end

"""
    generate_triangle(t, a)

Generates triangle wave.

# Arguments

- `t::Union{AbstractVector, AbstractRange}`: time vector
- `a::Real`: amplitude

# Returns

- s::Vector{Float64}`
"""
function generate_triangle(t::Union{AbstractVector, AbstractRange}, a::Real=1)

    s = @. a * abs(mod(t, 2) - 1)
    return s

end
