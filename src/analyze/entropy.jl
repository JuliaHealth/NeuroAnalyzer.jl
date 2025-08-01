export entropy
export negentropy

"""
    entropy(s)

Calculate entropy.

# Arguments

- `s::AbstractVector`

# Returns

Named tuple containing:
- `ent::Float64`
- `shent::Float64`: Shanon entropy
- `leent::Float64`: log energy entropy
- `sent::Float64`: sample entropy
- `nsent::Float64`: normalized sample entropy

# Note

Entropy of the signal is calculated using its histogram bins (number of bins is calculated using the Freedman-Diaconis formula) using the formulas `p = n / sum(n)` and `entropy = -sum(p .* log2(p))`, where `p` is the probability of each bin and `n` are bins' weights.

Shannon entropy and log energy entropy are calculated using `Wavelets.coefentropy()`.

Completely regular signals should have sample entropy approaching zero, while less regular signals should have higher sample entropy.
"""
function entropy(s::AbstractVector)::@NamedTuple{ent::Float64, shent::Float64, leent::Float64, sent::Float64, nsent::Float64}

    n = length(s)
    maxmin_range = maximum(s) - minimum(s)
    fd_bins = ceil(Int64, maxmin_range/(2.0 * iqr(s) * n^(-1/3))) # Freedman-Diaconis

    # recompute entropy with optimal bins for comparison
    h = StatsKit.fit(Histogram, s, nbins=fd_bins)
    hdat1 = h.weights ./ sum(h.weights)

    # convert histograms to probability values
    return (ent=-sum(hdat1 .* log2.(hdat1 .+ eps())),
            shent=Wavelets.coefentropy(s, ShannonEntropy()),
            leent=Wavelets.coefentropy(s, LogEnergyEntropy()),
            sent=ComplexityMeasures.complexity(SampleEntropy(s), s),
            nsent=ComplexityMeasures.complexity_normalized(SampleEntropy(s), s))

end

"""
    entropy(s)

Calculate entropy.

# Arguments

- `s::AbstractArray`

# Returns

Named tuple containing:
- `ent::Matrix{Float64}`
- `shent::Matrix{Float64}`: Shanon entropy
- `leent::Matrix{Float64}`: log energy entropy
- `sent::Matrix{Float64}`: sample entropy
- `nsent::Matrix{Float64}`: normalized sample entropy
"""
function entropy(s::AbstractArray)::@NamedTuple{ent::Matrix{Float64}, shent::Matrix{Float64}, leent::Matrix{Float64}, sent::Matrix{Float64}, nsent::Matrix{Float64}}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)

    ent = zeros(ch_n, ep_n)
    shent = zeros(ch_n, ep_n)
    leent = zeros(ch_n, ep_n)
    sent = zeros(ch_n, ep_n)
    nsent = zeros(ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads :greedy for ch_idx in 1:ch_n
            ent[ch_idx, ep_idx], shent[ch_idx, ep_idx], leent[ch_idx, ep_idx], sent[ch_idx, ep_idx], nsent[ch_idx, ep_idx] = @views entropy(s[ch_idx, :, ep_idx])
        end
    end

    return (ent=ent, shent=shent, leent=leent, sent=sent, nsent=nsent)

end

"""
    entropy(obj; <keyword arguments>)

Calculate entropy.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names

# Returns

Named tuple containing:
- `ent::Matrix{Float64}`
- `shent::Matrix{Float64}`: Shanon entropy
- `leent::Matrix{Float64}`: log energy entropy
- `sent::Matrix{Float64}`: sample entropy
- `nsent::Matrix{Float64}`: normalized sample entropy
"""
function entropy(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::@NamedTuple{ent::Matrix{Float64}, shent::Matrix{Float64}, leent::Matrix{Float64}, sent::Matrix{Float64}, nsent::Matrix{Float64}}

    ch = exclude_bads ? get_channel(obj, ch=ch, exclude="bad") : get_channel(obj, ch=ch, exclude="")
    ent, shent, leent, sent, nsent = @views entropy(obj.data[ch, :, :])

    return (ent=ent, shent=shent, leent=leent, sent=sent, nsent=nsent)

end

"""
    negentropy(signal)

Calculate negentropy.

# Arguments

- `signal::AbstractVector`

# Returns

- `ne::Float64`
"""
function negentropy(signal::AbstractVector)::Float64

    s = remove_dc(signal)

    ne = 0.5 * log(2 * pi * exp(1) * var(s)) - entropy(s)[1]

    return ne

end

"""
    negentropy(s)

Calculate negentropy.

# Arguments

- `s::AbstractArray`

# Returns

- `ne::Matrix{Float64}`
"""
function negentropy(s::AbstractArray)::Matrix{Float64}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)

    ne = zeros(ch_n, ep_n)

    @inbounds for ep_idx in 1:ep_n
        Threads.@threads :greedy for ch_idx in 1:ch_n
            ne[ch_idx, ep_idx] = @views negentropy(s[ch_idx, :, ep_idx])
        end
    end

    return ne

end

"""
    negentropy(obj; <keyword arguments>)

Calculate negentropy.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{String, Vector{String}, Regex}`: channel name or list of channel names

# Returns

- `ne::Matrix{Float64}`
"""
function negentropy(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::Matrix{Float64}

    ch = exclude_bads ? get_channel(obj, ch=ch, exclude="bad") : get_channel(obj, ch=ch, exclude="")
    ne = @views negentropy(obj.data[ch, :, :])

    return ne

end
