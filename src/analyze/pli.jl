export pli

"""
    pli(s1, s2)

Calculate PLI (Phase-Lag Index).

# Arguments

- `s1::AbstractVector`
- `s2::AbstractVector`

# Returns

Named tuple containing:
- `pv::Float64`: PLI value
- `sd::Vector{Float64}`: signal difference (s2 - s1)
- `phd::Vector{Float64}`: phase difference (s2 - s1)
- `s1ph::Vector{Float64}`: signal 1 phase
- `s2ph::Vector{Float64}`: signal 2 phase
"""
function pli(s1::AbstractVector, s2::AbstractVector)

    length(s1) == length(s2) || throw(ArgumentError("Both signals must have the same length."))

    _, _, _, s1ph = hspectrum(s1)
    _, _, _, s2ph = hspectrum(s2)

    sd = s2 - s1
    phd = s2ph - s1ph

    pv = abs(mean(sign.(imag.(exp.(1im .* phd)))))

    return (pv=pv, sd=sd, phd=phd, s1ph=s1ph, s2ph=s2ph)
end

"""
    pli(obj1, obj2; ch1, ch2, ep1, ep2)

Calculate PLI (Phase Lag Index).

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `ch1::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj1)`: index of channels, default is all signal channels
- `ch2::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj2)`: index of channels, default is all signal channels
- `ep1::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(epoch_n(obj1))`: default use all epochs
- `ep2::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(epoch_n(obj2))`: default use all epochs

# Returns

Named tuple containing:
- `pv::Array{Float64, 2}`: PLI value
- `sd::Array{Float64, 3}`: signal difference (s2 - s1)
- `phd::Array{Float64, 3}`: phase difference (s2 - s1)
- `s1ph::Array{Float64, 3}`: signal 1 phase
- `s2ph::Array{Float64, 3}`: signal 2 phase
"""
function pli(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch1::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj1), ch2::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj2), ep1::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(epoch_n(obj1)), ep2::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(epoch_n(obj2)))

    _check_channels(obj1, ch1)
    _check_channels(obj2, ch2)
    length(ch1) == length(ch2) || throw(ArgumentError("ch1 and ch2 must have the same length."))
    
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    length(ep1) == length(ep2) || throw(ArgumentError("ep1 and ep2 must have the same length."))
    epoch_len(obj1) == epoch_len(obj2) || throw(ArgumentError("OBJ1 and OBJ2 must have the same epoch lengths."))

    ep_n = length(ep1)
    ch_n = length(ch1)

    pv = zeros(ch_n, ep_n)
    sd = zeros(ch_n, epoch_len(obj1), ep_n)
    phd = zeros(ch_n, epoch_len(obj1), ep_n)
    s1ph = zeros(ch_n, epoch_len(obj1), ep_n)
    s2ph = zeros(ch_n, epoch_len(obj1), ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            pv[ch_idx, ep_idx], sd[ch_idx, :, ep_idx], phd[ch_idx, :, ep_idx], s1ph[ch_idx, :, ep_idx], s2ph[ch_idx, :, ep_idx] = @views pli(obj1.data[ch1[ch_idx], :, ep1[ep_idx]], obj2.data[ch2[ch_idx], :, ep2[ep_idx]])
        end
    end

    return (pv=pv, sd=sd, phase_dif=phd, s1ph=s1ph, s2ph=s2ph)
end

"""
    pli(obj; channel)

Calculate PLIs (Phase Lag Index).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

# Returns

- `pv::Array{Float64, 3}`: PLI value matrices over epochs
"""
function pli(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj))

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    pv = zeros(ch_n, ch_n, ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx1 in 1:ch_n
            for ch_idx2 in 1:ch_idx1
                pv[ch_idx1, ch_idx2, ep_idx], _, _, _, _ = @views pli(obj.data[channel[ch_idx1], :, ep_idx], obj.data[channel[ch_idx2], :, ep_idx])
            end
        end
    end

    pv = _copy_lt2ut(pv)

    return pv
end
