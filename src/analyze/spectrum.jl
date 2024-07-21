export spectrum
export hspectrum

"""
    spectrum(s; <keyword arguments>)

Calculate FFT, amplitudes, powers and phases.

# Arguments

- `s::AbstractVector`
- `pad::Int64=0`: number of zeros to add at the end of the signal
- `db::Bool=false`: normalize powers to dB

# Returns

Named tuple containing:
- `ft::Vector{ComplexF64}`: Fourier transforms
- `a::Vector{Float64}`: amplitudes
- `p::Vector{Float64}`: powers
- `ph::Vector{Float64}`: phases
"""
function spectrum(s::AbstractVector; pad::Int64=0, db::Bool=false)

    ft = rfft0(s, pad)

    # normalize
    ft ./= length(s)
    ft[2:end] .*= 2

    # amplitudes
    a = abs.(ft)

    # power
    p = abs.(ft .* conj(ft))       # p = a .^ 2;

    db && (p = pow2db.(p))

    # phases
    ph = angle.(ft)

    return (ft=ft, a=a, p=p, ph=ph)

end

"""
    hspectrum(s; <keyword arguments>)

Calculate amplitudes, powers and phases using Hilbert transform.

# Arguments

- `s::AbstractVector`
- `pad::Int64`: number of zeros to add at the end of the signal
- `db::Bool=false`: normalize powers to dB

# Returns

Named tuple containing:
- `hc::Vector(ComplexF64}`: Hilbert components
- `a::Vector{Float64}`: amplitudes
- `p::Vector{Float64}`: powers
- `ph::Vector{Float64}`: phases
"""
function hspectrum(s::AbstractVector; pad::Int64=0, db::Bool=false)

    hc = hilbert(pad0(s, pad))

    # amplitudes
    a = abs.(hc)

    # powers
    p = a.^2
    db && (p = pow2db.(p))

    # phases
    ph = angle.(hc)

    return (hc=hc, a=a, p=p, ph=ph)

end

"""
    hspectrum(s; <keyword arguments>)

Calculate amplitudes, powers and phases using Hilbert transform.

# Arguments

- `s::AbstractArray`
- `pad::Int64`: number of zeros to add at the end of the signal
- `db::Bool=false`: normalize powers to dB

# Returns

Named tuple containing:
- `hc::Array(ComplexF64, 3}`: Hilbert components
- `a::Array{Float64, 3}`: amplitudes
- `p::Array{Float64, 3}`: powers
- `ph::Array{Float64, 3}`: phases
"""
function hspectrum(s::AbstractArray; pad::Int64=0, db::Bool=false)

    ch_n = size(s, 1)
    ep_len = size(s, 2)
    ep_n = size(s, 3)

    hc = zeros(ComplexF64, ch_n, ep_len, ep_n)
    a = similar(s)
    p = similar(s)
    ph = similar(s)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            hc[ch_idx, :, ep_idx], a[ch_idx, :, ep_idx], p[ch_idx, :, ep_idx], ph[ch_idx, :, ep_idx] = @views hspectrum(s[ch_idx, :, ep_idx], pad=pad, db=db)
        end
    end

    return (hc=hc, a=a, p=p, ph=ph)

end

"""
    spectrum(s; <keyword arguments>)

Calculate FFT/Hilbert transformation components, amplitudes, powers and phases.

# Arguments

- `s::AbstractArray`
- `pad::Int64=0`: number of zeros to add signal for FFT
- `h::Bool=false`: use Hilbert transform for calculations instead of FFT
- `db::Bool=false`: normalize powers to dB

# Returns

Named tuple containing:
- `c::Array{ComplexF64, 3}`: Fourier or Hilbert components
- `a::Array{Float64, 3}`: amplitudes
- `p::Array{Float64, 3}`: powers
- `ph::Array{Float64, 3}: phase angles
"""
function spectrum(s::AbstractArray; pad::Int64=0, h::Bool=false, db::Bool=false)

    h && _warn("hspectrum() uses Hilbert transform, the signal should be narrowband for best results.")

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    if !h
        fft_size = div(size(s, 2) + pad, 2) + 1
    else
        fft_size = size(s, 2) + pad
    end

    c = zeros(ComplexF64, ch_n, fft_size, ep_n)
    ph = zeros(ch_n, fft_size, ep_n)
    a = zeros(ch_n, fft_size, ep_n)
    p = zeros(ch_n, fft_size, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            if h
                c[ch_idx, :, ep_idx], a[ch_idx, :, ep_idx], p[ch_idx, :, ep_idx], ph[ch_idx, :, ep_idx] = @views hspectrum(s[ch_idx, :, ep_idx], pad=pad, db=db)
            else
                c[ch_idx, :, ep_idx], a[ch_idx, :, ep_idx], p[ch_idx, :, ep_idx], ph[ch_idx, :, ep_idx] = @views spectrum(s[ch_idx, :, ep_idx], pad=pad, db=db)
            end
        end
    end

    return (c=c, a=a, p=p, ph=ph)
end

"""
    spectrum(obj; <keyword arguments>)

Calculate FFT/Hilbert transformation components, amplitudes, powers and phases.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}}`: list of channels
- `pad::Int64=0`: number of zeros to add signal for FFT
- `h::Bool=false`: use Hilbert transform for calculations instead of FFT
- `db::Bool=false`: normalize powers to dB

# Returns

Named tuple containing:
- `c::Array{ComplexF64, 3}`: Fourier or Hilbert components
- `a::Array{Float64, 3}`: amplitudes
- `p::Array{Float64, 3}`: powers
- `ph::Array{Float64, 3}: phase angles
"""
function spectrum(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}}, pad::Int64=0, h::Bool=false, db::Bool=false)

    ch = _ch_idx(obj, ch)
    c, a, p, ph = spectrum(obj.data[ch, :, :], pad=pad, h=h, db=db)

    return (c=c, a=a, p=p, ph=ph)

end
