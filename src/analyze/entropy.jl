export entropy
export negentropy

"""
    entropy(s)

Calculate entropy.

# Arguments

- `s::AbstractVector`

# Returns

- `ent::Float64`
- `sent::Float64`: Shanon entropy
- `leent::Float64`: log energy entropy
"""
function entropy(s::AbstractVector)

    n = length(s)
    maxmin_range = maximum(s) - minimum(s)
    fd_bins = ceil(Int64, maxmin_range/(2.0 * iqr(s) * n^(-1/3))) # Freedman-Diaconis

    # recompute entropy with optimal bins for comparison
    h = StatsKit.fit(Histogram, s, nbins=fd_bins)
    hdat1 = h.weights ./ sum(h.weights)

    # convert histograms to probability values
    return (ent=-sum(hdat1 .* log2.(hdat1 .+ eps())),
            sent=coefentropy(s, ShannonEntropy()),
            leent=coefentropy(s, LogEnergyEntropy()))
end

"""
    entropy(s)

Calculate entropy.

# Arguments

- `s::AbstractArray`

# Returns

Named tuple containing:
- `ent::Array{Float64, 2}`
- `sent::Array{Float64, 2}`: Shanon entropy
- `leent::Array{Float64, 2}`: log energy entropy
"""
function entropy(s::AbstractArray)

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    ent = zeros(ch_n, ep_n)
    sent = zeros(ch_n, ep_n)
    leent = zeros(ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            ent[ch_idx, ep_idx], sent[ch_idx, ep_idx], leent[ch_idx, ep_idx] = @views entropy(s[ch_idx, :, ep_idx])
        end
    end

    return (ent=ent, sent=sent, leent=leent)

end

"""
    entropy(obj; ch)

Calculate entropy.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

# Returns

Named tuple containing:
- `ent::Array{Float64, 2}`
- `sent::Array{Float64, 2}`: Shanon entropy
- `leent::Array{Float64, 2}`: log energy entropy
"""
function entropy(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj))

    _check_channels(obj, ch)
    isa(ch, Int64) && (ch = [ch])

    ent, sent, leent = @views entropy(obj.data[ch, :, :])

    return (ent=ent, sent=sent, leent=leent)

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

    s = remove_dc(signal)

    negent = 0.5 * log(2 * pi * exp(1) * var(s)) - entropy(s)[1]

    return negent

end

"""
    negentropy(s)

Calculate negentropy.

# Arguments

- `s::AbstractArray`

# Returns

- `ne::Array{Float64, 2}`
"""
function negentropy(s::AbstractArray)

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    ne = zeros(ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            ne[ch_idx, ep_idx] = @views negentropy(s[ch_idx, :, ep_idx])
        end
    end

    return ne

end

"""
    negentropy(obj; ch)

Calculate negentropy.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

# Returns

- `ne::Array{Float64, 2}`
"""
function negentropy(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj))

    _check_channels(obj, ch)
    isa(ch, Int64) && (ch = [ch])

    ne = @views negentropy(obj.data[ch, :, :])

    return ne

end
