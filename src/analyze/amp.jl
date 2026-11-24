export amp

"""
    _amp(s)

Computes amplitude descriptors for a 1-D signal vector.

# Arguments

- `s::AbstractVector`: signal vector

# Returns

Named tuple:

- `peak_amp::Float64`: peak amplitude (`max(|s|)`)
- `rms_amp::Float64`: RMS amplitude (`√mean(s²)`)
- `p2p_amp::Float64`: peak-to-peak amplitude (`max(s) - min(s)`)
- `semi_p2p_amp::Float64`: half of the peak-to-peak amplitude (`p2p/2`)
- `ms_amp::Float64`: mean square amplitude (`mean(s²)`)
- `te_signal::Float64`: total signal energy (`Σ s²`)
"""
function _amp(
    s::AbstractVector,
)::@NamedTuple{
    peak_amp::Float64,
    rms_amp::Float64,
    p2p_amp::Float64,
    semi_p2p_amp::Float64,
    ms_amp::Float64,
    te_signal::Float64,
}
    peak_amp = maximum(abs, s)
    rms_amp = _rms(s)
    s_min, s_max = extrema(s)
    p2p_amp = s_max - s_min
    semi_p2p_amp = p2p_amp / 2
    ms_amp = sum(abs2, s) / length(s)
    te_signal = sum(abs2, s)

    return (; peak_amp, rms_amp, p2p_amp, semi_p2p_amp, ms_amp, te_signal)
end

"""
    _amp(s)

Computes amplitude descriptors for a 3-D signal array.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)

# Returns

Named tuple:

- `peak_amp::Matrix{Float64}`: peak amplitude (`max(|s|)`), shape (channels, epochs)
- `rms_amp::Matrix{Float64}`: RMS amplitude (`√mean(s²)`), shape (channels, epochs)
- `p2p_amp::Matrix{Float64}`: peak-to-peak amplitude (`max(s) - min(s)`), shape (channels, epochs)
- `semi_p2p_amp::Matrix{Float64}`: half of the peak-to-peak amplitude, shape (channels, epochs)
- `ms_amp::Matrix{Float64}`: mean square amplitude (`mean(s²)`), shape (channels, epochs)
- `te_signal::Matrix{Float64}`: total signal energy (`Σ s²`), shape (channels, epochs)
"""
function _amp(
    s::AbstractArray,
)::@NamedTuple{
    peak_amp::Matrix{Float64},
    rms_amp::Matrix{Float64},
    p2p_amp::Matrix{Float64},
    semi_p2p_amp::Matrix{Float64},
    ms_amp::Matrix{Float64},
    te_signal::Matrix{Float64},
}
    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate output
    peak_amp = zeros(ch_n, ep_n)
    rms_amp = zeros(ch_n, ep_n)
    p2p_amp = zeros(ch_n, ep_n)
    semi_p2p_amp = zeros(ch_n, ep_n)
    ms_amp = zeros(ch_n, ep_n)
    te_signal = zeros(ch_n, ep_n)

    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        amp_data = _amp(@view(s[ch_idx, :, ep_idx]))
        peak_amp[ch_idx, ep_idx] = amp_data.peak_amp
        rms_amp[ch_idx, ep_idx] = amp_data.rms_amp
        p2p_amp[ch_idx, ep_idx] = amp_data.p2p_amp
        semi_p2p_amp[ch_idx, ep_idx] = amp_data.semi_p2p_amp
        ms_amp[ch_idx, ep_idx] = amp_data.ms_amp
        te_signal[ch_idx, ep_idx] = amp_data.te_signal
    end

    return (; peak_amp, rms_amp, p2p_amp, semi_p2p_amp, ms_amp, te_signal)
end

"""
    amp(obj; <keyword arguments>)

Computes amplitude descriptors for a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `exclude_bad_segments::Bool=NeuroAnalyzer.exclude_bad_segments`: if `true`, bad segments (markers of type `"BAD"` of all channels or of the selected channels) of a continuous signal are not used: their samples are removed and the remaining parts are joined (default: the `exclude_bad_segments` preference, see `set_exclude_bad_segments()`)
- `exclude_bad_channels::Bool=NeuroAnalyzer.exclude_bad_channels`: if `true`, channels marked as bad are not used (default: the `exclude_bad_channels` preference, see `set_exclude_bad_channels()`)

# Returns

Named tuple:

- `peak_amp::Matrix{Float64}`: peak amplitude (`max(|s|)`), shape (channels, epochs)
- `rms_amp::Matrix{Float64}`: RMS amplitude (`√mean(s²)`), shape (channels, epochs)
- `p2p_amp::Matrix{Float64}`: peak-to-peak amplitude (`max(s) - min(s)`), shape (channels, epochs)
- `semi_p2p_amp::Matrix{Float64}`: half of the peak-to-peak amplitude, shape (channels, epochs)
- `ms_amp::Matrix{Float64}`: mean square amplitude (`mean(s²)`), shape (channels, epochs)
- `te_signal::Matrix{Float64}`: total signal energy (`Σ s²`), shape (channels, epochs)
"""
function amp(
    obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex},
    exclude_bad_channels::Bool = NeuroAnalyzer.exclude_bad_channels,
    exclude_bad_segments::Bool = NeuroAnalyzer.exclude_bad_segments,
)::@NamedTuple{
    peak_amp::Matrix{Float64},
    rms_amp::Matrix{Float64},
    p2p_amp::Matrix{Float64},
    semi_p2p_amp::Matrix{Float64},
    ms_amp::Matrix{Float64},
    te_signal::Matrix{Float64},
}
    # resolve channel names to integer indices, optionally skipping bad channels
    ch = _get_ch(obj, ch, exclude_bad_channels)
    isempty(ch) && throw(ArgumentError("No channels selected."))

    return _amp(_seg_data(obj, ch, exclude_bad_segments))
end
