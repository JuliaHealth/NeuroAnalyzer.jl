export transform
export ftransform
export htransform
export hanalytic

"""
    ftransform(s; <keyword arguments>)

Calculate Fourier transformation.

# Arguments

- `s::AbstractVector`
- `pad::Int64=0`: number of zeros padding the signal
- `db::Bool=false`: normalize powers to dB

# Returns

Named tuple containing:
- `ft::Vector{ComplexF64}`: Fourier transform
- `a::Vector{Float64}`: amplitudes
- `p::Vector{Float64}`: powers
- `ph::Vector{Float64}`: phases (in radians)

# Notes

To get frequencies for the signal, use `f, _ = freqs(s, fs)`.
"""
function ftransform(s::AbstractVector; pad::Int64=0, db::Bool=false)::@NamedTuple{ft::Vector{ComplexF64}, a::Vector{Float64}, p::Vector{Float64}, ph::Vector{Float64}}

    # this will only return positive frequencies
    ft = rfft0(s, pad)

    # normalize
    ft ./= length(s)

    # multiple by 2 to compensate removed negative frequencies
    ft[2:end] .*= 2

    # amplitudes per frequencies
    a = abs.(ft)

    # powers
    p = abs2.(ft)               # p = a .^ 2 = abs.(ft .* conj(ft))
    db && (p = pow2db.(p))

    # remove very small values of |ft|, since they will affect calculations
    ft[abs.(ft) .< eps()] .= 0

    # phase
    ph = DSP.angle.(ft)         # ph = atan.(imag(ft), real(ft))
    # ph = DSP.unwrap(ph)

    return (ft=ft, a=a, p=p, ph=ph)

end

"""
    htransform(s; <keyword arguments>)

Calculate Hilbert transformation.

# Arguments

- `s::AbstractVector`
- `pad::Int64`: number of zeros padding the signal
- `db::Bool=false`: normalize powers to dB

# Returns

Named tuple containing:
- `ht::Vector{ComplexF64}`: Hilbert transform
- `a::Vector{Float64}`: amplitudes
- `p::Vector{Float64}`: powers
- `ph::Vector{Float64}`: phases (in radians)
"""
function htransform(s::AbstractVector; pad::Int64=0, db::Bool=false)::@NamedTuple{ht::Vector{ComplexF64}, a::Vector{Float64}, p::Vector{Float64}, ph::Vector{Float64}}

    # Hilbert transform
    ht = DSP.hilbert(pad0(s, pad))

    # instantaneous amplitude
    a = abs.(ht)

    # instantaneous phases
    ph = DSP.angle.(ht)

    # instantaneous power / energy
    p = abs2.(ht)
    db && (p = pow2db.(ht))

    return (ht=ht, a=a, p=p, ph=ph)

end

"""
    ftransform(s; <keyword arguments>)

Calculate Fourier transformation.

# Arguments

- `s::AbstractArray`
- `pad::Int64`: number of zeros padding the signal
- `db::Bool=false`: normalize powers to dB

# Returns

Named tuple containing:
- `ft::Array{ComplexF64, 3}`: Fourier transform
- `a::Array{Float64, 3}`: amplitudes
- `p::Array{Float64, 3}`: powers
- `ph::Array{Float64, 3}`: phases (in radians)
"""
function ftransform(s::AbstractArray; pad::Int64=0, db::Bool=false)::@NamedTuple{ht::Array{ComplexF64, 3}, a::Array{Float64, 3}, p::Array{Float64, 3}, ph::Array{Float64, 3}}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_len = size(s, 2)
    ep_n = size(s, 3)

    ft = zeros(ComplexF64, ch_n, ep_len, ep_n)
    a = similar(s)
    p = similar(s)
    ph = similar(s)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads :greedy for ch_idx in 1:ch_n
            ft[ch_idx, :, ep_idx], a[ch_idx, :, ep_idx], p[ch_idx, :, ep_idx], ph[ch_idx, :, ep_idx] = @views NeuroAnalyzer.ftransform(s[ch_idx, :, ep_idx], pad=pad, db=db)
        end
    end

    return (ft=ft, a=a, p=p, ph=ph)

end

"""
    htransform(s; <keyword arguments>)

Calculate Hilbert transformation.

# Arguments

- `s::AbstractArray`
- `pad::Int64`: number of zeros padding the signal
- `db::Bool=false`: normalize powers to dB

# Returns

Named tuple containing:
- `ht::Array{ComplexF64, 3}`: Hilbert transform
- `a::Array{Float64, 3}`: amplitudes
- `p::Array{Float64, 3}`: powers
- `ph::Array{Float64, 3}`: phases (in radians)
"""
function htransform(s::AbstractArray; pad::Int64=0, db::Bool=false)::@NamedTuple{ht::Array{ComplexF64, 3}, a::Array{Float64, 3}, p::Array{Float64, 3}, ph::Array{Float64, 3}}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_len = size(s, 2)
    ep_n = size(s, 3)

    ht = zeros(ComplexF64, ch_n, ep_len, ep_n)
    a = similar(s)
    p = similar(s)
    ph = similar(s)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads :greedy for ch_idx in 1:ch_n
            ht[ch_idx, :, ep_idx], a[ch_idx, :, ep_idx], p[ch_idx, :, ep_idx], ph[ch_idx, :, ep_idx] = @views NeuroAnalyzer.htransform(s[ch_idx, :, ep_idx], pad=pad, db=db)
        end
    end

    return (ht=ht, a=a, p=p, ph=ph)

end

"""
    transform(s; <keyword arguments>)

Calculate Fourier/Hilbert transformation.

# Arguments

- `s::AbstractArray`
- `pad::Int64=0`: number of zeros padding the signal
- `h::Bool=false`: perform Hilbert transformation
- `db::Bool=false`: normalize powers to dB

# Returns

Named tuple containing:
- `t::Array{ComplexF64, 3}`: Fourier or Hilbert transform
- `a::Array{Float64, 3}`: amplitudes
- `p::Array{Float64, 3}`: powers
- `ph::Array{Float64, 3}: phases (in radians)
"""
function transform(s::AbstractArray; pad::Int64=0, h::Bool=false, db::Bool=false)::@NamedTuple{t::Array{ComplexF64, 3}, a::Array{Float64, 3}, p::Array{Float64, 3}, ph::Array{Float64, 3}}

    _chk3d(s)
    h && _warn("htransform() uses Hilbert transform, the signal should be narrowband for best results.")

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    if !h
        fft_size = div(size(s, 2) + pad, 2) + 1
    else
        fft_size = size(s, 2) + pad
    end

    t = zeros(ComplexF64, ch_n, fft_size, ep_n)
    ph = zeros(ch_n, fft_size, ep_n)
    a = zeros(ch_n, fft_size, ep_n)
    p = zeros(ch_n, fft_size, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads :greedy for ch_idx in 1:ch_n
            if h
                t[ch_idx, :, ep_idx], a[ch_idx, :, ep_idx], p[ch_idx, :, ep_idx], ph[ch_idx, :, ep_idx] = @views NeuroAnalyzer.htransform(s[ch_idx, :, ep_idx], pad=pad, db=db)
            else
                t[ch_idx, :, ep_idx], a[ch_idx, :, ep_idx], p[ch_idx, :, ep_idx], ph[ch_idx, :, ep_idx] = @views NeuroAnalyzer.transform(s[ch_idx, :, ep_idx], pad=pad, db=db)
            end
        end
    end

    return (t=t, a=a, p=p, ph=ph)
end

"""
    transform(obj; <keyword arguments>)

Calculate Fourier/Hilbert transformation.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `pad::Int64=0`: number of zeros to add signal for FFT
- `h::Bool=false`: use Hilbert transform for calculations instead of FFT
- `db::Bool=false`: normalize powers to dB

# Returns

Named tuple containing:
- `t::Array{ComplexF64, 3}`: Fourier or Hilbert transform
- `a::Array{Float64, 3}`: amplitudes
- `p::Array{Float64, 3}`: powers
- `ph::Array{Float64, 3}: phases (in radians)
"""
function transform(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, pad::Int64=0, h::Bool=false, db::Bool=false)::@NamedTuple{t::Array{ComplexF64, 3}, a::Array{Float64, 3}, p::Array{Float64, 3}, ph::Array{Float64, 3}}

    ch = exclude_bads ? get_channel(obj, ch=ch, exclude="bad") : get_channel(obj, ch=ch, exclude="")
    t, a, p, ph = NeuroAnalyzer.transform(obj.data[ch, :, :], pad=pad, h=h, db=db)

    return (t=t, a=a, p=p, ph=ph)

end

"""
    hanalytic(s)

Calculate analytic signal using Hilbert transformation.

# Arguments

- `s::AbstractVector`
- `pad::Int64`: number of zeros padding the signal

# Returns

- `ha::Vector{ComplexF64}`:
"""
function hanalytic(s::AbstractVector; pad::Int64)::Vector{ComplexF64}

    ha = s + im * DSP.hilbert(pad0(s, pad))

    return ha

end

"""
    hanalytic(s; <keyword arguments>)

Calculate analytic signal using Hilbert transformation.

# Arguments

- `s::AbstractArray`
- `pad::Int64`: number of zeros padding the signal

# Returns

- `ha::Vector{ComplexF64}`:
"""
function hanalytic(s::AbstractArray; pad::Int64=0)::Array{ComplexF64, 3}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_len = size(s, 2)
    ep_n = size(s, 3)

    ha = zeros(ComplexF64, ch_n, ep_len, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads :greedy for ch_idx in 1:ch_n
            ha = @views NeuroAnalyzer.hanalytic(s[ch_idx, :, ep_idx], pad=pad)
        end
    end

    return ha

end

"""
    hanalytic(obj; <keyword arguments>)

Calculate analytic signal using Hilbert transformation.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `pad::Int64=0`: number of zeros to add signal for FFT

# Returns

- `ha::Vector{ComplexF64}`:
"""
function hanalytic(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, pad::Int64=0, h::Bool=false, db::Bool=false)::Array{ComplexF64, 3}

    ch = exclude_bads ? get_channel(obj, ch=ch, exclude="bad") : get_channel(obj, ch=ch, exclude="")
    ha = NeuroAnalyzer.hanalytic(obj.data[ch, :, :], pad=pad)

    return ha

end