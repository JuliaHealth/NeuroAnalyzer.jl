export difference
export chdiff
export phdiff
export ampdiff

"""
    difference(signal1, signal2; n, method)

Calculate mean difference and 95% confidence interval for 2 signals.

# Arguments

- `signal1::AbstractArray`
- `signal2::AbstractArray`
- `n::Int64=3`: number of bootstraps
- `method::Symbol=:absdiff`: maximum difference (`:absdiff`), integrated area of the squared difference (`:diff2int`)

# Returns

Named tuple containing:
- `s_stat::Vector{Float64}`
- `s_stat_single::Float64`
- `p::Float64`
"""
function difference(signal1::AbstractArray, signal2::AbstractArray; n::Int64=3, method::Symbol=:absdiff)

    size(signal1) == size(signal2) || throw(ArgumentError("Both signals must be of the same size."))
    _check_var(method, [:absdiff, :diff2int], "method")
    n < 1 && throw(ArgumentError("n must be ≥ 1."))

    s1_mean = vec(mean(signal1, dims=1))
    s2_mean = vec(mean(signal2, dims=1))

    if method === :absdiff
        # statistic: maximum difference
        s_diff = s1_mean - s2_mean
        s_stat_single = maximum(abs.(s_diff))
    else
        # statistic: integrated area of the squared difference
        s_diff_squared = (s1_mean - s2_mean).^2
        s_stat_single = simpson(s_diff_squared)
    end

    signals = [signal1; signal2]
    s_stat = zeros(size(signal1, 1) * n)

    Threads.@threads for idx1 in 1:(size(signal1, 1) * n)
        s_tmp1 = zeros(size(signal1, 1), size(signal1, 2))
        sample_idx = rand(1:size(signals, 1), size(signals, 1))
        @inbounds @simd for idx2 in 1:size(signal1, 1)
            s_tmp1[idx2, :] = @views signals[sample_idx[idx2], :]'
        end
        s1_mean = vec(mean(s_tmp1, dims=1))
        s_tmp1 = zeros(size(signal1, 1), size(signal1, 2))
        sample_idx = rand(1:size(signals, 1), size(signals, 1))
        @inbounds @simd for idx2 in 1:size(signal1, 1)
            s_tmp1[idx2, :] = @views signals[sample_idx[idx2], :]'
        end
        s2_mean = vec(mean(s_tmp1, dims=1))
        if method === :absdiff
            # statistic: maximum difference
            s_diff = s1_mean - s2_mean
            @inbounds s_stat[idx1] = maximum(abs.(s_diff))
        else
            # statistic: integrated area of the squared difference
            s_diff_squared = (s1_mean - s2_mean).^2
            @inbounds s_stat[idx1] = simpson(s_diff_squared)
        end
    end

    p = length(s_stat[s_stat .> s_stat_single]) / size(signal1, 1) * n
    p > 1 && (p = 1.0)

    return (s_stat=s_stat, s_stat_single=s_stat_single, p=p)
end

"""
    difference(obj; channel, n, method)

Calculate mean difference and its 95% CI between channels.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `n::Int64=3`: number of bootstraps
- `method::Symbol=:absdiff`
    - `:absdiff`: maximum difference
    - `:diff2int`: integrated area of the squared difference

# Returns

Named tuple containing:
- `signals_statistic::Matrix{Float64}`
- `signals_statistic_single::Vector{Float64}`
- `p::Vector{Float64}`
"""
function difference(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj), n::Int64=3, method::Symbol=:absdiff)

    ep_n = epoch_n(obj)
    _check_channels(obj, channel)

    s_stat = zeros(ep_n, length(channel) * n)
    s_stat_single = zeros(ep_n)
    p = zeros(ep_n)

    @inbounds @simd for ep_idx in 1:ep_n
        s_stat[ep_idx, :], s_stat_single[ep_idx], p[ep_idx] = difference(obj.data[channel, :, ep_idx], obj.data[channel, :, ep_idx], n=n, method=method)
    end

    return (s_stat=s_stat, s_stat_single=s_stat_single, p=p)
end

"""
    difference(obj1, obj2; channel1, channel2, epoch1, epoch2, n, method)

Calculates mean difference and 95% confidence interval for two channels.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2:NeuroAnalyzer.NEURO`
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj1.header.recording[:data_type]))`: index of channels, default is all channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj2, type=Symbol(obj2.header.recording[:data_type]))`: index of channels, default is all channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2))`: default use all epochs
- `n::Int64`: number of bootstraps
- `method::Symbol[:absdiff, :diff2int]`
    - `:absdiff`: maximum difference
    - `:diff2int`: integrated area of the squared difference

# Returns

Named tuple containing:
- `s_stat::Matrix{Float64}`
- `s_stat_single::Vector{Float64}`
- `p::Vector{Float64}`
"""
function difference(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; channel1::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj1.header.recording[:data_type])), channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj2, type=Symbol(obj2.header.recording[:data_type])), epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1)), epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2)), n::Int64=3, method::Symbol=:absdiff)

    _check_channels(obj1, channel1)
    _check_channels(obj2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))
    
    _check_epochs(obj1, epoch1)
    _check_epochs(obj2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    epoch_len(obj1) == epoch_len(obj2) || throw(ArgumentError("OBJ1 and OBJ2 epoch lengths must be equal."))

    ep_n = length(epoch1)

    s_stat = zeros(ep_n, length(channel1) * n)
    s_stat_single = zeros(ep_n)
    p = zeros(ep_n)

    Threads.@threads for ep_idx in 1:ep_n
        s_stat[ep_idx, :], s_stat_single[ep_idx], p[ep_idx] = @views difference(obj1.signals[channel1, :, epoch1[ep_idx]], obj2.signals[channel2, :, epoch2[ep_idx]], n=n, method=method)
    end

    return (s_stat=s_stat, statsitic_single=s_stat_single, p=p)
end

"""
    chdiff(obj1, obj2; channel1, channel2, epoch1, epoch2)

Subtract channels.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj1.header.recording[:data_type]))`: index of channels, default is all channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj2, type=Symbol(obj2.header.recording[:data_type]))`: index of channels, default is all channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2))`: default use all epochs

# Returns

- `ch_diff::Matrix{Float64}`
"""
function chdiff(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; channel1::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj1.header.recording[:data_type])), channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj2, type=Symbol(obj2.header.recording[:data_type])), epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1)), epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2)))

    _check_channels(obj1, channel1)
    _check_channels(obj2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))
    
    _check_epochs(obj1, epoch1)
    _check_epochs(obj2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    epoch_len(obj1) == epoch_len(obj2) || throw(ArgumentError("obj1 and obj2 epoch lengths must be equal."))

    ep_n = length(epoch1)
    ch_n = length(channel1)

    ch_diff = zeros(ch_n, epoch_len(obj1), ep_n)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            ch_diff[ch_idx, :, ep_idx] = @views obj1.signals[channel1[ch_idx], :, epoch1[ep_idx]] .- obj2.signals[channel2[ch_idx], :, epoch2[ep_idx]]
        end
    end

    return ch_diff
end

"""
    phdiff(signal1, signal2; pad, h)

Calculate phase difference between signals.

# Arguments

- `signal1::AbstractVector`
- `signal2::AbstractVector`
- `pad::Int64=0`: number of zeros to add
- `h::Bool=false`: use FFT or Hilbert transformation (if h=true)

# Returns

Named tuple containing:
- `ph_diff::Vector{Float64}`: phase differences in radians
"""
function phdiff(signal1::AbstractVector, signal2::AbstractVector; pad::Int64=0, h::Bool=false)

    if h
        _, _, _, ph1 = hspectrum(signal1, pad=pad)
        _, _, _, ph2 = hspectrum(signal2, pad=pad)
    else
        _, _, _, ph1 = spectrum(signal1, pad=pad)
        _, _, _, ph2 = spectrum(signal2, pad=pad)
    end

    return round.(ph1 - ph2, digits=2)
end

"""
    phdiff(obj; channel, pad, h)

Calculate phase difference between channels and mean phase of reference `channel`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of reference channels, default is all OBJ/MEG channels except the analyzed one
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

"""
    ampdiff(obj; channel)

Calculate amplitude difference between each channel and mean amplitude.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of reference channels, default is all OBJ/MEG channels except the analyzed one

# Returns
 
- `amp_diff::Array{Float64, 3}`
"""
function ampdiff(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj))

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    amp_diff = zeros(ch_n, epoch_len(obj), ep_n)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            ref_channels = setdiff(channel, ch_idx)
            amp_ref = @views vec(mean(obj.data[ref_channels, :, ep_idx], dims=1))
            amp_diff[ch_idx, :, ep_idx] = @views obj.data[channel[ch_idx], :, ep_idx] - amp_ref
        end
    end

    return amp_diff
end
