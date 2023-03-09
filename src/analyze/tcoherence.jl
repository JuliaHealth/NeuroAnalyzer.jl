export tcoherence

"""
    tcoherence(signal1, signal2; pad)

Calculate coherence (mean over time), IC (imaginary coherence) and MSC (magnitude-squared coherence) between `signal1` and `signal2`.

# Arguments

- `signal1::AbstractVector`
- `signal2::AbstractVector`
- `pad::Int64=0`: number of zeros to add

# Returns

Named tuple containing:
- `c::Vector{Float64}`: coherence
- `msc::Vector{Float64}`: magnitude-squares coherence
- `ic::Vector{Float64}`: imaginary part of coherence
"""
function s2_tcoherence(signal1::AbstractVector, signal2::AbstractVector; pad::Int64=0)

    length(signal1) == length(signal2) || throw(ArgumentError("Both signals must have the same length."))

    s1_fft = fft0(signal1, pad) ./ length(signal1)
    s2_fft = fft0(signal2, pad) ./ length(signal2)

    coh = @. (abs((s1_fft) * conj.(s2_fft))^2) / (s1_fft * s2_fft)

    msc = @. abs(coh)^2

    return (c=real.(coh), msc=msc, ic=imag.(coh))
end

"""
    tcoherence(obj1, obj2; channel1, channel2, epoch1, epoch2)

Calculate coherence (mean over time), IC (imaginary coherence) and MSC (magnitude-squared coherence).

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj1.header.recording[:data_type]))`: index of channels, default is all signal channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj2, type=Symbol(obj2.header.recording[:data_type]))`: index of channels, default is all signal channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2))`: default use all epochs
- `pad::Int64=0`: number of zeros to add signal for FFT

# Returns

Named tuple containing:
- `c::Array{Float64, 3}`: coherence
- `msc::Array{Float64, 3}`: MSC
- `ic::Array{Float64, 3}`: imaginary part of coherence
"""
function tcoherence(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; channel1::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj1.header.recording[:data_type])), channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj2, type=Symbol(obj2.header.recording[:data_type])), epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1)), epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2)), pad::Int64=0)

    _check_channels(obj1, channel1)
    _check_channels(obj2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))

    _check_epochs(obj1, epoch1)
    _check_epochs(obj2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    epoch_len(obj1) == epoch_len(obj2) || throw(ArgumentError("OBJ1 and OBJ2 epoch lengths must be equal."))

    ep_n = length(epoch1)
    ch_n = length(channel1)

    c = zeros(length(channel1), epoch_len(obj1), length(epoch1))
    msc = zeros(length(channel1), epoch_len(obj1), length(epoch1))
    ic = zeros(length(channel1), epoch_len(obj1), length(epoch1))

    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            c[ch_idx, :, ep_idx], msc[ch_idx, :, ep_idx], ic[ch_idx, :, ep_idx] = @views s2_tcoherence(obj1.data[channel1[ch_idx], :, epoch1[ep_idx]], obj2.data[channel2[ch_idx], :, epoch2[ep_idx]], pad=pad)
        end
    end

    return (c=c, msc=msc, ic=ic)
end
