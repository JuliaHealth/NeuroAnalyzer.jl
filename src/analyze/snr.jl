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
function snr(s::AbstractVector)::Float64

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
function snr2(s::AbstractVector)::Float64

    a = amp(s)
    return (maximum(s) - minimum(s)) / a.rmsq

end

"""
    snr(s; <keyword arguments>)

Calculate SNR.

# Arguments

- `s::AbstractArray`
- `t::Vector{Float64}`: epoch time
- `type::Symbol=:rms`: SNR type:
    - `:mean`: mean-based
    - `:rms`: RMS-based

# Returns

Named tuple containing:
- `sn::Matrix{Float64}`: SNR for each channel over frequencies 1:Nyquist
- `f::Vector{Float64}`: frequencies
"""
function snr(s::AbstractArray; t::Vector{Float64}, type::Symbol=:rms)::@NamedTuple{sn::Matrix{Float64}, f::Vector{Float64}}

    _check_var(type, [:mean, :rms], "type")
    _chk3d(s)

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    @assert ep_n >= 2 "OBJ must contain ≥ 2 epochs."

    f, _ = freqs(t)
    sp = @views NeuroAnalyzer.transform(s[1, :, 1])
    amp = zeros(ch_n, length(sp.a), ep_n)
    sn = zeros(ch_n, length(f))

    # create spectrum for each channel
    @inbounds for ep_idx in 1:ep_n
        Threads.@threads :greedy for ch_idx in 1:ch_n
            _, amp[ch_idx, :, ep_idx], _, _ = @views NeuroAnalyzer.transform(s[ch_idx, :, ep_idx])
        end
    end

    # calculate SNR for each channel spectrum
    @inbounds for f_idx in eachindex(f)
        Threads.@threads :greedy for ch_idx in 1:ch_n
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
    snr(obj; <keyword arguments>)

Calculate SNR.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names
- `type::Symbol=:rms`: SNR type:
    - `:mean`: mean-based
    - `:rms`: RMS-based

# Returns

Named tuple containing:
- `sn::Matrix{Float64}`: SNR for each channel over frequencies 1:Nyquist
- `f::Vector{Float64}`: frequencies
"""
function snr(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, type::Symbol=:rms)::@NamedTuple{sn::Matrix{Float64}, f::Vector{Float64}}

    ch = exclude_bads ? get_channel(obj, ch=ch, exclude="bad") : get_channel(obj, ch=ch, exclude="")
    sn, f = @views snr(obj.data[ch, :, :], t=obj.epoch_time, type=type)

    return (sn=sn, f=f)

end
