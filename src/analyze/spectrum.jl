export spectrum
export hspectrum

"""
    spectrum(s; pad, norm)

Calculate FFT, amplitudes, powers and phases.

# Arguments

- `s::AbstractVector`
- `pad::Int64=0`: number of zeros to add at the end of the signal
- `norm::Bool=false`: normalize do dB

# Returns

Named tuple containing:
- `ft::Vector{ComplexF64}`: Fourier transforms
- `sa::Vector{Float64}`: amplitudes
- `sp::Vector{Float64}`: powers
- `sph::Vector{Float64}`: phases
"""
function spectrum(s::AbstractVector; pad::Int64=0, norm::Bool=false)

    ft = fft0(s, pad) / length(s)

    # amplitudes
    sa = @. abs(ft)                          # normalize
    sa = @views sa[1:(length(sa) ÷ 2)]       # remove negative frequencies
    sa[2:end] .*= 2                          # double positive frequencies

    # replace amplitudes at extreme frequencies
    # sa[1] = sa[2]
    # sa[end] = sa[end - 1]

    # power
    sp = sa.^2
    norm == true && (sp = pow2db.(sp))

    # phases
    sph = angle.(ft)

    return (ft=ft, sa=sa, sp=sp, sph=sph)

end

"""
    hspectrum(s; pad=0)

Calculate amplitudes, powers and phases using Hilbert transform.

# Arguments

- `s::AbstractVector`
- `pad::Int64`: number of zeros to add at the end of the signal
- `norm::Bool=false`: normalize do dB

# Returns

Named tuple containing:
- `hc::Vector(ComplexF64}`: Hilbert components
- `sa::Vector{Float64}`: amplitudes
- `sp::Vector{Float64}`: powers
- `sph::Vector{Float64}`: phases
"""
function hspectrum(s::AbstractVector; pad::Int64=0, norm::Bool=false)

    hc = hilbert(pad0(s, pad))

    # amplitudes
    sa = abs.(hc)

    # powers
    sp = sa.^2
    norm == true && (sp = pow2db.(sp))

    # phases
    sph = angle.(hc)

    return (hc=hc, sa=sa, sp=sp, sph=sph)
    
end

"""
    hspectrum(s; pad, norm)

Calculate amplitudes, powers and phases using Hilbert transform.

# Arguments

- `s::AbstractArray`
- `pad::Int64`: number of zeros to add at the end of the signal
- `norm::Bool=false`: normalize do dB

# Returns

Named tuple containing:
- `hc::Array(ComplexF64, 3}`: Hilbert components
- `sa::Array{Float64, 3}`: amplitudes
- `sp::Array{Float64, 3}`: powers
- `sph::Array{Float64, 3}`: phases
"""
function hspectrum(s::AbstractArray; pad::Int64=0, norm::Bool=false)

    ch_n = size(s, 1)
    ep_len = size(s, 2)
    ep_n = size(s, 3)

    hc = zeros(ComplexF64, ch_n, ep_len, ep_n)
    sa = similar(s)
    sp = similar(s)
    sph = similar(s)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            hc[ch_idx, :, ep_idx], sa[ch_idx, :, ep_idx], sp[ch_idx, :, ep_idx], sph[ch_idx, :, ep_idx] = @views hspectrum(s[ch_idx, :, ep_idx], pad=pad, norm=norm)
        end
    end  

    return (hc=hc, sa=sa, sp=sp, sph=sph)
    
end

"""
    spectrum(s; pad, h, norm)

Calculate FFT/Hilbert transformation components, amplitudes, powers and phases.

# Arguments

- `s::AbstractArray`
- `pad::Int64=0`: number of zeros to add signal for FFT
- `h::Bool=false`: use Hilbert transform for calculations instead of FFT
- `norm::Bool=false`: normalize do dB

# Returns

Named tuple containing:
- `c::Array{ComplexF64, 3}`: Fourier or Hilbert components
- `sa::Array{Float64, 3}`: amplitudes
- `sp::Array{Float64, 3}`: powers
- `sph::Array{Float64, 3}: phase angles
"""
function spectrum(s::AbstractArray; pad::Int64=0, h::Bool=false, norm::Bool=false)

    ch_n = size(s, 1)
    ep_n = size(s, 3)
    fft_size = size(s, 2) + pad

    c = zeros(ComplexF64, ch_n, fft_size, ep_n)
    sph = zeros(ch_n, fft_size, ep_n)

    if h == true
        _warn("hspectrum() uses Hilbert transform, the signal should be narrowband for best results.")
        sa = zeros(ch_n, fft_size, ep_n)
        sp = zeros(ch_n, fft_size, ep_n)
    else
        sa = zeros(ch_n, fft_size ÷ 2, ep_n)
        sp = zeros(ch_n, fft_size ÷ 2, ep_n)
    end        

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            if h == true
                c[ch_idx, :, ep_idx], sa[ch_idx, :, ep_idx], sp[ch_idx, :, ep_idx], sph[ch_idx, :, ep_idx] = @views hspectrum(s[ch_idx, :, ep_idx], pad=pad, norm=norm)
            else
                c[ch_idx, :, ep_idx], sa[ch_idx, :, ep_idx], sp[ch_idx, :, ep_idx], sph[ch_idx, :, ep_idx] = @views spectrum(s[ch_idx, :, ep_idx], pad=pad, norm=norm)
            end
        end
    end  

    return (c=c, sa=sa, sp=sp, sph=sph)
end


"""
    spectrum(obj; ch, pad, h, norm)

Calculate FFT/Hilbert transformation components, amplitudes, powers and phases.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `pad::Int64=0`: number of zeros to add signal for FFT
- `h::Bool=false`: use Hilbert transform for calculations instead of FFT
- `norm::Bool=false`: normalize do dB

# Returns

Named tuple containing:
- `c::Array{ComplexF64, 3}`: Fourier or Hilbert components
- `sa::Array{Float64, 3}`: amplitudes
- `sp::Array{Float64, 3}`: powers
- `sph::Array{Float64, 3}: phase angles
"""
function spectrum(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), pad::Int64=0, h::Bool=false, norm::Bool=false)

    _check_channels(obj, ch)

    c, sa, sp, sph = spectrum(obj.data[ch, :, :], pad=pad, h=h, norm=norm)

    return (c=c, sa=sa, sp=sp, sph=sph)

end
