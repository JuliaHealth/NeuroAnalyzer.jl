export entropy
export negentropy

"""
    entropy(signal)

Calculate entropy.

# Arguments

- `signal::AbstractVector`

# Returns

- `ent::Float64`
- `sent::Float64`: Shanon entropy
- `leent::Float64`: log energy entropy
"""
function entropy(signal::AbstractVector)

    n = length(signal)
    maxmin_range = maximum(signal) - minimum(signal)
    fd_bins = ceil(Int64, maxmin_range/(2.0 * iqr(signal) * n^(-1/3))) # Freedman-Diaconis

    # recompute entropy with optimal bins for comparison
    h = StatsKit.fit(Histogram, signal, nbins=fd_bins)
    hdat1 = h.weights ./ sum(h.weights)

    # convert histograms to probability values
    return (ent=-sum(hdat1 .* log2.(hdat1 .+ eps())),
            sent=coefentropy(signal, ShannonEntropy()),
            leent=coefentropy(signal, LogEnergyEntropy()))
end

"""
    negentropy(signal)

Calculate negentropy.

# Arguments

- `signal::AbstractVector`

# Returns

- `negent::Float64`
"""
function negentropy(signal::AbstractVector)
    s = demean(signal)
    return 0.5 * log(2 * pi * exp(1) * var(s)) - entropy(s)[1]
end

"""
    entropy(obj; channel)

Calculate entropy.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels

# Returns

Named tuple containing:
- `ent::Array{Float64, 2}`
- `s_ent::Array{Float64, 2}`: Shanon entropy
- `le_ent::Array{Float64, 2}`: log energy entropy
"""
function entropy(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj))

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    ent = zeros(ch_n, ep_n)
    sent = zeros(ch_n, ep_n)
    leent = zeros(ch_n, ep_n)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            ent[ch_idx, ep_idx], sent[ch_idx, ep_idx], leent[ch_idx, ep_idx] = @views s_entropy(obj.data[channel[ch_idx], :, ep_idx])
        end
    end

    return (ent=ent, s_ent=sent, le_ent=leent)
end

"""
    negentropy(obj; channel)

Calculate negentropy.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj)`: index of channels, default is all EEG channels

# Returns

- `ne::Matrix{Float64}`
"""
function negentropy(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj))

    _check_channels(obj, channel)
    ch_n = length(channel)
    ep_n = epoch_n(obj)

    ne = zeros(ch_n, ep_n)
    @inbounds @simd for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            ne[ch_idx, ep_idx] = @views s_negentropy(obj.data[channel[ch_idx], :, ep_idx])
        end
    end

    return ne
end
