export phdiff

"""
    phdiff(s1, s2; pad, h)

Calculate phase difference between signals.

# Arguments

- `s1::AbstractVector`
- `s2::AbstractVector`
- `pad::Int64=0`: number of zeros to add
- `h::Bool=false`: use FFT or Hilbert transformation (if h=true)

# Returns

Named tuple containing:
- `ph_diff::Vector{Float64}`: phase differences in radians
"""
function phdiff(s1::AbstractVector, s2::AbstractVector; pad::Int64=0, h::Bool=false)

    if h
        _, _, _, ph1 = hspectrum(s1, pad=pad)
        _, _, _, ph2 = hspectrum(s2, pad=pad)
    else
        _, _, _, ph1 = spectrum(s1, pad=pad)
        _, _, _, ph2 = spectrum(s2, pad=pad)
    end

    return round.(ph1 - ph2, digits=2)
end

"""
    phdiff(obj; channel, pad, h)

Calculate phase difference between channels and mean phase of reference `channel`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of reference channels, default is all signal channels except the analyzed one
- `avg::Symbol=:phase`: method of averaging:
    - `:phase`: phase is calculated for each reference channel separately and then averaged
    - `:signal`: signals are averaged prior to phase calculation
- `pad::Int64=0`: pad signals with 0s
- `h::Bool=false`: use FFT or Hilbert transformation

# Returns
 
- `ph_diff::Array{Float64, 3}`
"""
function phdiff(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), avg::Symbol=:phase, pad::Int64=0, h::Bool=false)

    avg in [:phase, :signal] || throw(ArgumentError("avg must be :phase or :signal."))

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    ph_diff = zeros(ch_n, epoch_len(obj), ep_n)
    if avg === :phase
        @inbounds @simd for ep_idx in 1:ep_n
            Threads.@threads for ch_idx in 1:ch_n
                ref_channels = setdiff(channel, ch_idx)
                ph_ref = zeros(length(ref_channels), epoch_len(obj))
                for ref_idx in eachindex(ref_channels)
                    if h == true
                        _, _, _, ph = @views hspectrum(obj.data[ref_channels[ref_idx], :, ep_idx], pad=pad)
                    else
                        _, _, _, ph = @views spectrum(obj.data[ref_channels[ref_idx], :, ep_idx], pad=pad)
                    end
                    ph_ref[ref_idx, :] = ph
                end
                ph_ref = vec(mean(ph_ref, dims=1))
                if h == true
                    _, _, _, ph = @views hspectrum(obj.data[channel[ch_idx], :, ep_idx], pad=pad)
                else
                    _, _, _, ph = @views spectrum(obj.data[channel[ch_idx], :, ep_idx], pad=pad)
                end
                ph_diff[ch_idx, :, ep_idx] = ph - ph_ref
            end
        end
    else
        @inbounds @simd for ep_idx in 1:ep_n
            Threads.@threads for ch_idx in 1:ch_n
                ref_channels = setdiff(channel, ch_idx)
                signal_m = @views vec(mean(obj.data[ref_channels, :, ep_idx], dims=1))
                ph_diff[ch_idx, :, ep_idx] = @views phdiff(obj.data[channel[ch_idx], :, ep_idx], signal_m)
            end
        end
    end

    return ph_diff
end
