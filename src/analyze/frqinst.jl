export frqinst

"""
    frqinst(signal; fs)

Calculate instantaneous frequency.

# Arguments

- `signal::AbstractVector`
- `fs::Int64`

# Returns

- `frqinst::Vector{Float64}`
"""
function frqinst(signal::AbstractVector; fs::Int64)
    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))
    _, _, _, h_phases = hspectrum(signal)
    return 256 * derivative(h_phases) / (2*pi)
end

"""
    frqinst(obj; channel)

Calculate instantaneous frequency.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

# Returns

- `f::Array{Float64, 3}`
"""
function frqinst(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj))

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    f = zeros(ch_n, epoch_len(obj), ep_n)
    fs = sr(obj)

    _info("frqinst() uses Hilbert transform, the signal should be narrowband for best results.")

    # initialize progress bar
    progress_bar == true && (p = Progress(ep_n * ch_n, 1))
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            f[ch_idx, :, ep_idx] = @views frqinst(obj.data[channel[ch_idx], :, ep_idx], fs=fs)
        end
        # update progress bar
        progress_bar == true && next!(p)
    end
    return f
end
