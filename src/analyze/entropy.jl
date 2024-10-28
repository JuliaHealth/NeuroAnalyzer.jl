export entropy
export negentropy

"""
    entropy(s)

Calculate entropy.

# Arguments

- `s::Vector{<:Real}`

# Returns

Named tuple containing:
- `ent::Float64`
- `sent::Float64`: Shanon entropy
- `leent::Float64`: log energy entropy
"""
function entropy(s::Vector{<:Real})::@NamedTuple{ent::Float64, sent::Float64, leent::Float64}

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

- `s::Array{<:Real, 3}`

# Returns

Named tuple containing:
- `ent::Matrix{Float64}`
- `sent::Matrix{Float64}`: Shanon entropy
- `leent::Matrix{Float64}`: log energy entropy
"""
function entropy(s::Array{<:Real, 3})::@NamedTuple{ent::Matrix{Float64}, sent::Matrix{Float64}, leent::Matrix{Float64}}

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    ent = zeros(ch_n, ep_n)
    sent = zeros(ch_n, ep_n)
    leent = zeros(ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            ent[ch_idx, ep_idx], sent[ch_idx, ep_idx], leent[ch_idx, ep_idx] = entropy(s[ch_idx, :, ep_idx])
        end
    end

    return (ent=ent, sent=sent, leent=leent)

end

"""
    entropy(obj; <keyword arguments>)

Calculate entropy.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}}`: channel name or list of channel names

# Returns

Named tuple containing:
- `ent::Matrix{Float64}`
- `sent::Matrix{Float64}`: Shanon entropy
- `leent::Matrix{Float64}`: log energy entropy
"""
function entropy(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}})::@NamedTuple{ent::Matrix{Float64}, sent::Matrix{Float64}, leent::Matrix{Float64}}

    ch = get_channel(obj, ch=ch)
    ent, sent, leent = entropy(obj.data[ch, :, :])

    return (ent=ent, sent=sent, leent=leent)

end

"""
    negentropy(signal)

Calculate negentropy.

# Arguments

- `signal::Vector{<:Real}`

# Returns

- `negent::Float64`
"""
function negentropy(signal::Vector{<:Real})::Float64

    s = remove_dc(signal)

    negent = 0.5 * log(2 * pi * exp(1) * var(s)) - entropy(s)[1]

    return negent

end

"""
    negentropy(s)

Calculate negentropy.

# Arguments

- `s::Array{<:Real, 3}`

# Returns

- `ne::Matrix{Float64}`
"""
function negentropy(s::Array{<:Real, 3})::Matrix{Float64}

    ch_n = size(s, 1)
    ep_n = size(s, 3)

    ne = zeros(ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads for ch_idx in 1:ch_n
            ne[ch_idx, ep_idx] = negentropy(s[ch_idx, :, ep_idx])
        end
    end

    return ne

end

"""
    negentropy(obj; <keyword arguments>)

Calculate negentropy.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}}`: channel name or list of channel names

# Returns

- `ne::Matrix{Float64}`
"""
function negentropy(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}})::Matrix{Float64}

    ch = get_channel(obj, ch=ch)
    ne = negentropy(obj.data[ch, :, :])

    return ne

end
