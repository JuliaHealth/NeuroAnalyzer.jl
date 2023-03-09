export fcoherence

"""
    fcoherence(signal; fs, frq_lim)

Calculate coherence (mean over all frequencies) and MSC (magnitude-squared coherence) between channels.

# Arguments

- `signal::AbstractArray`
- `fs::Int64`: sampling rate
- `frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing`: return coherence only for the given frequency range

# Returns

Named tuple containing:
- `c::Array{Float64, 3}`: coherence
- `msc::Array{Float64, 3}`: MSC
- `f::Vector{Float64}`: frequencies
"""
function fcoherence(signal::AbstractArray; fs::Int64, frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing)

    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))
    c = mt_coherence(signal, fs=fs)
    f = Vector(c.freq)
    c = c.coherence
    if frq_lim !== nothing
        frq_lim = tuple_order(frq_lim)
        frq_lim[1] < 0 && throw(ArgumentError("Lower frequency bound must be ≥ 0."))
        frq_lim[2] > fs / 2 && throw(ArgumentError("Upper frequency bound must be ≤ $fs."))
        idx1 = vsearch(frq_lim[1], f)
        idx2 = vsearch(frq_lim[2], f)
        c = c[:, :, idx1:idx2]
        f = f[idx1:idx2]
    end

    return (c=c, msc=c.^2, f=f)
end

"""
    fcoherence(signal1, signal2; fs, frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing)

Calculate coherence (mean over all frequencies) and MSC (magnitude-squared coherence) between channels of `signal1` and `signal2`.

# Arguments

- `signal1::AbstractArray`
- `signal2::AbstractArray`
- `fs::Int64`
- `frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing`: return coherence only for the given frequency range

# Returns

- `c::Array{Float64, 3}`: coherence
- `msc::Array{Float64, 3}`: MSC
- `f::Vector{Float64}`: frequencies
"""
function fcoherence(signal1::AbstractArray, signal2::AbstractArray; fs::Int64, frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing)

    length(signal1) == length(signal2) || throw(ArgumentError("Both signals must have the same length."))
    fs < 1 && throw(ArgumentError("fs must be ≥ 1."))

    signal = hcat(signal1, signal2)'
    c = mt_coherence(signal, fs=fs)
    f = Vector(c.freq)
    c = c.coherence
    if frq_lim !== nothing
        frq_lim = tuple_order(frq_lim)
        frq_lim[1] < 0 && throw(ArgumentError("Lower frequency bound must be ≥ 0."))
        frq_lim[2] > fs / 2 && throw(ArgumentError("Upper frequency bound must be ≤ $fs."))
        idx1 = vsearch(frq_lim[1], f)
        idx2 = vsearch(frq_lim[2], f)
        c = c[:, :, idx1:idx2]
        f = f[idx1:idx2]
    end
    c = c[1, 2, :]
    msc = @. abs(c)^2
    
    return (c=c, msc=msc, f=f)
end

"""
    fcoherence(obj1, obj2; channel1, channel2, epoch1, epoch2, frq_lim)

Calculate coherence (mean over frequencies) and MSC (magnitude-squared coherence).

# Arguments

- `obj1::NeuroAnalyzer.NEURO`
- `obj2::NeuroAnalyzer.NEURO`
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj1.header.recording[:data_type]))`: index of channels, default is all channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj2, type=Symbol(obj2.header.recording[:data_type]))`: index of channels, default is all channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2))`: default use all epochs
- `frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing`: return coherence only for the given frequency range

# Returns

Named tuple containing:
- `c::Array{Float64, 3}`: coherence
- `msc::Array{Float64, 3}`: MSC
- `f::Vector{Float64}`: frequencies
"""
function fcoherence(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; channel1::Union{Int64, Vector{Int64}, AbstractRange}=0, channel2::Union{Int64, Vector{Int64}, AbstractRange}=0, epoch1::Union{Int64, Vector{Int64}, AbstractRange}=0, epoch2::Union{Int64, Vector{Int64}, AbstractRange}=0, frq_lim::Union{Tuple{Real, Real}, Nothing}=nothing)

    _check_channels(obj1, channel1)
    _check_channels(obj2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))
    
    _check_epochs(obj1, epoch1)
    _check_epochs(obj2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    epoch_len(obj1) == epoch_len(obj2) || throw(ArgumentError("obj1 and obj2 epoch lengths must be equal."))

    sr(obj1) == sr(obj2) || throw(ArgumentError("OBJ1 and OBJ2 must have the same sampling rate."))

    c_tmp, _, f = @views s2_fcoherence(obj1.signals[1, :, 1], obj1.signals[1, :, 1], fs=sr(obj1), frq_lim=frq_lim)
    c = zeros(length(channel1), length(c_tmp), length(epoch1))
    msc = zeros(length(channel1), length(c_tmp), length(epoch1))
    f = zeros(length(channel1), length(c_tmp), length(epoch1))
    @inbounds @simd for ep_idx in eachindex(epoch1)
        Threads.@threads for ch_idx in eachindex(channel1)
            c[ch_idx, :, ep_idx], msc[ch_idx, :, ep_idx], _ = @views s2_fcoherence(obj1.signals[channel1[ch_idx], :, epoch1[ep_idx]], obj2.signals[channel2[ch_idx], :, epoch2[ep_idx]], fs=sr(obj1), frq_lim=frq_lim)
        end
    end

    return (c=c, msc=msc, f=f)
end
