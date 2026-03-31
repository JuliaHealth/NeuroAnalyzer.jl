export snr
export snr2

"""
    snr(s1, s2)

Calculate SNR between two 1-D signal vectors.

# Arguments

- `s1::AbstractVector`: signal vector
- `s2::AbstractVector`: signal vector

# Returns

- `Float64`: SNR
"""
function snr(s1::AbstractVector, s2::AbstractVector)::Float64
    return -20 * log10(norm(abs.(s2 - s1)) / norm(s2))
end

"""
    snr(s)

Calculate mean-based SNR for a 1-D signal vector.

# Arguments

- `s::AbstractVector`: signal vector

# Returns

- `Float64`: SNR

# References

D. J. Schroeder (1999). Astronomical optics (2nd ed.). Academic Press. ISBN 978-0-12-629810-9, p.278
"""
function snr(s::AbstractVector)::Float64
    return mean(s) / std(s)
end

"""
    snr2(s)

Calculate RMS-based SNR for a 1-D signal vector.

# Arguments

- `s::AbstractVector`: signal vector

# Returns

- `Float64`: SNR
"""
function snr2(s::AbstractVector)::Float64
    a = amp(s)
    return (maximum(s) - minimum(s)) / a.rmsq
end

"""
    snr(s; <keyword arguments>)

Calculate SNR for a 3-D signal array.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `t::Vector{Float64}`: epoch time
- `type::Symbol=:rms`: SNR type:
    - `:mean`: mean-based
    - `:rms`: RMS-based

# Returns

Named tuple:

- `sn::Matrix{Float64}`: SNR for each channel over frequencies 1:Nyquist
- `f::Vector{Float64}`: frequencies
"""
function snr(
    s::AbstractArray;
    t::Vector{Float64},
    type::Symbol = :rms,
)::@NamedTuple{
    sn::Matrix{Float64},
    f::Vector{Float64},
}

    # validate
    _check_var(type, [:mean, :rms], "type")

    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # validate
    ep_n >= 2 || throw(ArgumentError("OBJ must contain ≥ 2 epochs."))

    f, _ = freqs(t)
    sp = @views NeuroAnalyzer.ftransform(s[1, :, 1])
    amp = zeros(ch_n, length(sp.a), ep_n)
    sn = zeros(ch_n, length(f))

    # create spectrum for each channel
    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        amp[ch_idx, :, ep_idx] = NeuroAnalyzer.ftransform(@view(s[ch_idx, :, ep_idx])).a
    end

    # calculate SNR for each channel spectrum
    @inbounds Threads.@threads :static for idx in CartesianIndices((length(f), ch_n))
        f_idx, ch_idx = idx[1], idx[2]
        if type === :mean
            sn[ch_idx, f_idx] = snr(@view(amp[ch_idx, f_idx, :]))
        else
            sn[ch_idx, f_idx] = snr2(@view(amp[ch_idx, f_idx, :]))
        end
    end

    return (; sn, f)
end

"""
    snr(obj; <keyword arguments>)

Calculate SNR for a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `type::Symbol=:rms`: SNR type:
    - `:mean`: mean-based
    - `:rms`: RMS-based

# Returns

Named tuple:

- `sn::Matrix{Float64}`: SNR for each channel over frequencies 1:Nyquist
- `f::Vector{Float64}`: frequencies
"""
function snr(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    type::Symbol = :rms,
)::@NamedTuple{
    sn::Matrix{Float64},
    f::Vector{Float64},
}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch =
        exclude_bads ?
        get_channel(obj; ch = ch, exclude = "bad") :
        get_channel(obj; ch = ch, exclude = "")
    isempty(ch) && throw(ArgumentError("No channels selected."))

    return snr(@view(obj.data[ch, :, :]); t = obj.epoch_time, type = type)
end
