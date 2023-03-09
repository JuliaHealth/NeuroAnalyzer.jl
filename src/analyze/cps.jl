export cps

"""
    cps(signal1, signal2; fs, norm)

Calculate cross power spectrum between `signal1` and `signal2`.

# Arguments

- `signal1::AbstractVector`
- `signal2::AbstractVector`
- `fs::Int64`: sampling rate
- `norm::Bool=true`: normalize do dB

# Returns

Named tuple containing:
- `cps_pw::Vector{Float64}`: cross power spectrum power
- `cps_ph::Vector{Float64}`: cross power spectrum phase (in radians)
- `cps_fq::Vector{Float64}`: cross power spectrum frequencies
"""
function cps(signal1::AbstractVector, signal2::AbstractVector; fs::Int64, norm::Bool=true)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))
    s = hcat(signal1, signal2)'
    p = mt_cross_power_spectra(s, fs=fs)
    cps_pw = real.(p.power)[1, 2, :]
    cps_ph = angle.(imag.(p.power))[1, 2, :]
    cps_fq = Vector(p.freq)

    return (cps_pw=cps_pw, cps_ph=cps_ph, cps_fq=cps_fq)
end

"""
    cps(obj; channel, norm)

Calculate cross power spectrum.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `norm::Bool=true`: normalize do dB

# Returns

Named tuple containing:
- `cps_pw::Array{Float64, 4}`: cross power spectrum power
- `cps_ph::Array{Float64, 4}`: cross power spectrum phase (in radians)
- `cps_fq::Vector{Float64}`: cross power spectrum frequencies
"""
function cps(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), norm::Bool=true)

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    fs = sr(obj)
    
    cps_pw_tmp, cps_ph_tmp, cps_fq = @views s2_cps(obj.data[1, :, 1], obj.data[1, :, 1], fs=fs)
    cps_pw = zeros(ch_n, ch_n, length(cps_pw_tmp), ep_n)
    cps_ph = zeros(ch_n, ch_n, length(cps_ph_tmp), ep_n)

    # initialize progress bar
    progress_bar == true && (p = Progress(ep_n * ch_n, 1))

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx1 in 1:ch_n
           for ch_idx2 in 1:ch_idx1
                cps_pw[ch_idx1, ch_idx2, :, ep_idx], cps_ph[ch_idx1, ch_idx2, :, ep_idx], _ = @views cps(obj.data[channel[ch_idx1], :, ep_idx], obj.data[channel[ch_idx2], :, ep_idx], fs=fs, norm=norm)
            end

        # update progress bar
        progress_bar == true && next!(p)
        end
    end

    @inbounds @simd for time_idx in 1:size(cps_pw, 3)
        Threads.@threads for ep_idx in 1:ep_n
            for ch_idx1 in 1:(ch_n - 1)
                for ch_idx2 in (ch_idx1 + 1):ch_n
                    cps_pw[ch_idx1, ch_idx2, time_idx, ep_idx] = @views cps_pw[ch_idx2, ch_idx1, time_idx, ep_idx]
                    cps_ph[ch_idx1, ch_idx2, time_idx, ep_idx] = @views cps_ph[ch_idx2, ch_idx1, time_idx, ep_idx]
                end
            end
        end
    end

    return (cps_pw=cps_pw, cps_ph=cps_ph, cps_fq=cps_fq)
end

"""
    cps(obj1, obj2; channel1, channel2, epoch1, epoch2, norm)

Calculate cross power spectrum.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj1.header.recording[:data_type]))`: index of channels, default is all channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj2, type=Symbol(obj2.header.recording[:data_type]))`: index of channels, default is all channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2))`: default use all epochs
- `norm::Bool=true`: normalize do dB

# Returns

Named tuple containing:
- `cps_pw::Vector{Float64}`: cross power spectrum power
- `cps_ph::Vector{Float64}`: cross power spectrum phase (in radians)
- `cps_fq::Vector{Float64}`: cross power spectrum frequencies
"""
function cps(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; channel1::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj1.header.recording[:data_type])), channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj2, type=Symbol(obj2.header.recording[:data_type])), epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1)), epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2)), norm::Bool=true)

    sr(obj1) == sr(obj2) || throw(ArgumentError("OBJ1 and OBJ2 must have the same sampling rate."))

    _check_channels(obj1, channel1)
    _check_channels(obj2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))
    
    _check_epochs(obj1, epoch1)
    _check_epochs(obj2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    epoch_len(obj1) == epoch_len(obj2) || throw(ArgumentError("obj1 and obj2 epoch lengths must be equal."))

    ep_n = length(epoch1)
    ch_n = length(channel1)
    fs = sr(obj1)

    cps_pw, cps_ph, cps_fq = @views s2_cps(obj1.data[1, :, 1], obj2.data[1, :, 1], fs=fs, norm=norm)

    cps_pw = zeros(ch_n, length(cps_pw), ep_n)
    cps_ph = zeros(ch_n, length(cps_ph), ep_n)
    cps_fq = zeros(ch_n, length(cps_fq), ep_n)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            cps_pw[ch_idx, :, ep_idx], cps_ph[ch_idx, :, ep_idx], cps_fq[ch_idx, :, ep_idx] = @views cps(obj1.data[channel1[ch_idx], :, epoch1[ep_idx]], obj2.data[channel2[ch_idx], :, epoch2[ep_idx]], fs=fs, norm=norm)
        end
    end

    return (cps_pw=cps_pw, cps_ph=cps_ph, cps_fq=cps_fq)
end
