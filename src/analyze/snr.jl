export snr
export snr2

"""
    snr(s)

Calculate mean-based SNR.

# Arguments

- `s::AbstractVector`

# Returns

- `snr::Float64`: SNR

# Source

D. J. Schroeder (1999). Astronomical optics (2nd ed.). Academic Press. ISBN 978-0-12-629810-9, p.278
"""
function snr(s::AbstractVector)

    return mean(s) / std(s)

end

"""
    snr2(s)

Calculate RMS-based SNR.

# Arguments

- `s::AbstractVector`

# Returns

- `snr2::Float64`: SNR
"""
function snr2(s::AbstractVector)
    
    return (maximum(s) - minimum(s)) / rms(s)

end

"""
    snr(s; t, type)

Calculate SNR.

# Arguments

- `s::AbstractArray`
- `t::Vector{Float64}`: epoch time
- `type::Symbol=:rms`: SNR type:
    - `:mean`: mean-based
    - `:rms`: RMS-based

# Returns

Named tuple containing:
- `s::Matrix(Float64)`: SNR for each channel over frequencies 1:Nyquist
- `f::Vector(Float64)`: frequencies
"""
function snr(s::AbstractArray; t::Vector{Float64}, type::Symbol=:rms)

    _check_var(type, [:mean, :rms], "type")

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    ep_n == 1 && throw(ArgumentError("OBJ must contain ≥ 2 epochs."))

    f, _ = freqs(t)
    _, amp, _, _ = @views spectrum(s[1, :, 1])
    amp = zeros(ch_n, length(amp), ep_n)
    sn = zeros(ch_n, length(f))

    # create spectrum for each channel
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            _, amp[ch_idx, :, ep_idx], _, _ = @views spectrum(s[ch_idx, :, ep_idx])
        end
    end

    # calculate SNR for each channel spectrum
    @inbounds @simd for f_idx in eachindex(f)
        Threads.@threads for ch_idx in 1:ch_n
            if type === :mean
                sn[ch_idx, f_idx] = @views snr(amp[ch_idx, f_idx, :])
            else
                sn[ch_idx, f_idx] = @views snr2(amp[ch_idx, f_idx, :])
            end
        end
    end

    return (sn=sn, f=f)
    
end


"""
    snr(obj; ch, type)

Calculate SNR.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `type::Symbol=:rms`: SNR type:
    - `:mean`: mean-based
    - `:rms`: RMS-based

# Returns

Named tuple containing:
- `sn::Matrix(Float64)`: SNR for each channel over frequencies 1:Nyquist
- `f::Vector(Float64)`: frequencies
"""
function snr(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), type::Symbol=:rms)

    _check_channels(obj, ch)

    sn, f = @views snr(obj.data[ch, :, :], t=obj.epoch_time, type=type)

    return (sn=sn, f=f)

end
