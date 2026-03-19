export entropy
export negentropy

"""
    entropy(s)

Calculate signal entropy descriptors:

- histogram-based entropy in bits (Freedman-Diaconis binning)
- Shannon entropy (Wavelets.coefentropy)
- log energy entropy (Wavelets.coefentropy)
- sample entropy (ComplexityMeasures)
- normalised sample entropy (ComplexityMeasures)

# Arguments

- `s::AbstractVector`: signal vector

# Returns

Named tuple:

- `ent::Float64`: entropy in bits
- `shent::Float64`: Shanon entropy
- `leent::Float64`: log energy entropy
- `sent::Float64`: sample entropy
- `nsent::Float64`: normalized sample entropy

# Note

Histogram entropy uses Freedman-Diaconis binning: `p = n / sum(n)`, `ent = −Σ p·log₂(p)`. Shannon and log energy entropy use `Wavelets.coefentropy()`. Sample entropy approaches zero for completely regular signals and increases
with irregularity.
"""
function entropy(
    s::AbstractVector
)::@NamedTuple{
    ent::Float64,
    shent::Float64,
    leent::Float64,
    sent::Float64,
    nsent::Float64
}

    n = length(s)

    # Freedman-Diaconis rule: optimal bin width = 2·IQR·N^(−1/3).
    maxmin_range = maximum(s) - minimum(s)
    fd_bins = ceil(Int64, maxmin_range / (2.0 * iqr(s) * n^(-1 / 3)))

    # fit histogram and convert bin counts to probabilities
    h = StatsKit.fit(Histogram, s, nbins = fd_bins)
    p = h.weights ./ sum(h.weights)

    # histogram entropy in bits; eps() guards against log(0)
    ent = -sum(p .* log2.(p .+ eps()))

    # construct SampleEntropy estimator once and reuse for both sent and nsent
    se = SampleEntropy(s)

    return (
        ent = ent,
        shent = Wavelets.coefentropy(s, ShannonEntropy()),
        leent = Wavelets.coefentropy(s, LogEnergyEntropy()),
        sent = ComplexityMeasures.complexity(se, s),
        nsent = ComplexityMeasures.complexity_normalized(se, s),
    )

end

"""
    entropy(s)

Calculate signal entropy descriptors:

- histogram-based entropy in bits (Freedman-Diaconis binning)
- Shannon entropy (Wavelets.coefentropy)
- log energy entropy (Wavelets.coefentropy)
- sample entropy (ComplexityMeasures)
- normalised sample entropy (ComplexityMeasures)

# Arguments

- `s::AbstractArray`: signal array (channels, samples, epochs)

# Returns

Named tuple:

- `ent::Matrix{Float64}`: entropy in bits, shape `(channels, epochs)`
- `shent::Matrix{Float64}`: Shanon entropy, shape `(channels, epochs)`
- `leent::Matrix{Float64}`: log energy entropy, shape `(channels, epochs)`
- `sent::Matrix{Float64}`: sample entropy, shape `(channels, epochs)`
- `nsent::Matrix{Float64}`: normalized sample entropy, shape `(channels, epochs)`
"""
function entropy(
    s::AbstractArray
)::@NamedTuple{
    ent::Matrix{Float64},
    shent::Matrix{Float64},
    leent::Matrix{Float64},
    sent::Matrix{Float64},
    nsent::Matrix{Float64},
}

    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate outputs
    ent = zeros(ch_n, ep_n)
    shent = zeros(ch_n, ep_n)
    leent = zeros(ch_n, ep_n)
    sent = zeros(ch_n, ep_n)
    nsent = zeros(ch_n, ep_n)

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        entropy_data = entropy(@view(s[ch_idx, :, ep_idx]))
        ent[ch_idx, ep_idx]   = entropy_data.ent
        shent[ch_idx, ep_idx] = entropy_data.shent
        leent[ch_idx, ep_idx] = entropy_data.leent
        sent[ch_idx, ep_idx]  = entropy_data.sent
        nsent[ch_idx, ep_idx] = entropy_data.nsent
    end

    return (ent = ent, shent = shent, leent = leent, sent = sent, nsent = nsent)

end

"""
    entropy(obj; <keyword arguments>)

Calculate signal entropy descriptors:

- histogram-based entropy in bits (Freedman-Diaconis binning)
- Shannon entropy (Wavelets.coefentropy)
- log energy entropy (Wavelets.coefentropy)
- sample entropy (ComplexityMeasures)
- normalised sample entropy (ComplexityMeasures)

# Returns

Named tuple:

- `ent::Matrix{Float64}`: entropy in bits, shape `(channels, epochs)`
- `shent::Matrix{Float64}`: Shanon entropy, shape `(channels, epochs)`
- `leent::Matrix{Float64}`: log energy entropy, shape `(channels, epochs)`
- `sent::Matrix{Float64}`: sample entropy, shape `(channels, epochs)`
- `nsent::Matrix{Float64}`: normalized sample entropy, shape `(channels, epochs)`
"""
function entropy(
        obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}
    )::@NamedTuple{
        ent::Matrix{Float64}, shent::Matrix{Float64}, leent::Matrix{Float64}, sent::Matrix{Float64}, nsent::Matrix{Float64},
    }

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")
    entropy_data = entropy(@view(obj.data[ch, :, :]))

    return entropy_data

end

"""
    negentropy(signal)

Calculate negentropy. Negentropy measures how far a signal's distribution departs from Gaussian: `ne = 0.5·ln(2πe·var(s)) − H(s)`, where `H(s)` is the histogram entropy. ne ≈ 0 for Gaussian; ne > 0 for distributions that are more structured (peaky, multi-modal, etc.).

# Arguments

- `signal::AbstractVector`: signal vector

# Returns

- `ne::Float64`: negentropy (≥ 0; equals 0 for a Gaussian signal)
"""
function negentropy(signal::AbstractVector)::Float64

    # remove DC offset so variance reflects only signal variability
    s = remove_dc(signal)

    # Gaussian differential entropy: 0.5·ln(2πe·σ²).
    # ℯ is the built-in mathematical constant (more readable than exp(1)).
    gaussian_h = 0.5 * log(2 * π * ℯ * var(s))

    ne = gaussian_h - entropy(s).ent

    return ne

end

"""
    negentropy(s)

Calculate negentropy. Negentropy measures how far a signal's distribution departs from Gaussian: `ne = 0.5·ln(2πe·var(s)) − H(s)`, where `H(s)` is the histogram entropy. ne ≈ 0 for Gaussian; ne > 0 for distributions that are more structured (peaky, multi-modal, etc.).

# Arguments

- `s::AbstractArray`: signal array (channels, samples, epochs)

# Returns

- `ne::Matrix{Float64}`: negentropy (≥ 0; equals 0 for a Gaussian signal), shape `(channel, epochs)`
"""
function negentropy(s::AbstractArray)::Matrix{Float64}

    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate output
    ne = zeros(ch_n, ep_n)

    # initialize progress bar
    progbar = Progress(ep_n * ch_n, dt = 1, barlen = 20, color = :white, enabled = progress_bar)

    # calculate over channel and epochs
    @inbounds Threads.@threads :dynamic for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        ne[ch_idx, ep_idx] = negentropy(@view(s[ch_idx, :, ep_idx]))
        progress_bar && next!(progbar)
    end

    return ne

end

"""
    negentropy(obj; <keyword arguments>)

Calculate negentropy. Negentropy measures how far a signal's distribution departs from Gaussian: `ne = 0.5·ln(2πe·var(s)) − H(s)`, where `H(s)` is the histogram entropy. ne ≈ 0 for Gaussian; ne > 0 for distributions that are more structured (peaky, multi-modal, etc.).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)

# Returns

- `ne::Matrix{Float64}`: negentropy (≥ 0; equals 0 for a Gaussian signal), shape `(channel, epochs)`
"""
function negentropy(obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex})::Matrix{Float64}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch = exclude_bads ? get_channel(obj, ch = ch, exclude = "bad") : get_channel(obj, ch = ch, exclude = "")

    ne = negentropy(@view(obj.data[ch, :, :]))

    return ne

end
