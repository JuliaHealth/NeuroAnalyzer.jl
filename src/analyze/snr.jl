export snr
export snr2

"""
    snr(s)

Calculate mean-based SNR.

# Arguments

- `s::Vector{<:Real}`

# Returns

- `snr::Float64`: SNR

# Source

D. J. Schroeder (1999). Astronomical optics (2nd ed.). Academic Press. ISBN 978-0-12-629810-9, p.278
"""
function snr(s::Vector{<:Real})::Float64

    return mean(s) / std(s)

end

"""
    snr2(s)

Calculate RMS-based SNR.

# Arguments

- `s::Vector{<:Real}`

# Returns

- `snr2::Float64`: SNR
"""
function snr2(s::Vector{<:Real})::Float64

    a = amp(s)
    return (maximum(s) - minimum(s)) / a.rms

end

"""
    snr(s; <keyword arguments>)

Calculate SNR.

# Arguments

- `s::Array{<:Real, 3}`
- `t::Vector{Float64}`: epoch time
- `type::Symbol=:rms`: SNR type:
    - `:mean`: mean-based
    - `:rms`: RMS-based

# Returns

Named tuple containing:
- `sn::Matrix{Float64}`: SNR for each channel over frequencies 1:Nyquist
- `f::Vector{Float64}`: frequencies
"""
function snr(s::Array{<:Real, 3}; t::Vector{Float64}, type::Symbol=:rms)::@NamedTuple{sn::Matrix{Float64}, f::Vector{Float64}}

    _check_var(type, [:mean, :rms], "type")

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    @assert ep_n >= 2 "OBJ must contain ≥ 2 epochs."

    f, _ = freqs(t)
    _, amp, _, _ = spectrum(s[1, :, 1])
    amp = zeros(ch_n, length(amp), ep_n)
    sn = zeros(ch_n, length(f))

    # create spectrum for each channel
    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            _, amp[ch_idx, :, ep_idx], _, _ = spectrum(s[ch_idx, :, ep_idx])
        end
    end

    # calculate SNR for each channel spectrum
    @inbounds for f_idx in eachindex(f)
        Threads.@threads for ch_idx in 1:ch_n
            sn[ch_idx, f_idx] = type === :mean ? snr(amp[ch_idx, f_idx, :]) : snr2(amp[ch_idx, f_idx, :])
        end
    end

    return (sn=sn, f=f)

end

"""
    snr(obj; <keyword arguments>)

Calculate SNR.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}}`: channel name or list of channel names
- `type::Symbol=:rms`: SNR type:
    - `:mean`: mean-based
    - `:rms`: RMS-based

# Returns

Named tuple containing:
- `sn::Matrix{Float64}`: SNR for each channel over frequencies 1:Nyquist
- `f::Vector{Float64}`: frequencies
"""
function snr(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}}, type::Symbol=:rms)::@NamedTuple{sn::Matrix{Float64}, f::Vector{Float64}}

    ch = get_channel(obj, ch=ch)
    sn, f = snr(obj.data[ch, :, :], t=obj.epoch_time, type=type)

    return (sn=sn, f=f)

end
