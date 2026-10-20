export entropy
export negentropy

"""
    _entropy(s)

Calculate signal entropy descriptors for a 1-D signal vector:

- histogram-based entropy in bits (Freedman-Diaconis binning)
- Shannon entropy (Wavelets.coefentropy)
- log energy entropy (Wavelets.coefentropy)
- sample entropy (ComplexityMeasures)
- normalized sample entropy (ComplexityMeasures)
- differential entropy

# Arguments

- `s::AbstractVector`: signal vector

# Returns

Named tuple:

- `ent::Float64`: entropy in bits
- `shent::Float64`: Shanon entropy
- `leent::Float64`: log energy entropy
- `sent::Float64`: sample entropy
- `nsent::Float64`: normalized sample entropy
- `dent::Float64`: differential entropy

# Note

Histogram entropy uses Freedman-Diaconis binning: `p = n / sum(n)`, `ent = −Σ p·log₂(p)`. Shannon and log energy entropy use `Wavelets.coefentropy()`. Sample entropy approaches zero for completely regular signals and increases
with irregularity.
"""
function _entropy(
    s::AbstractVector,
)::@NamedTuple{
    ent::Float64,
    shent::Float64,
    leent::Float64,
    sent::Float64,
    nsent::Float64,
    dent::Float64,
}
    n = length(s)

    # Freedman-Diaconis rule: optimal bin width = 2·IQR·N^(−1/3).
    maxmin_range = maximum(s) - minimum(s)
    q = StatsKit.iqr(s)
    # IQR = 0 (e.g. flat or heavily quantized signal): fall back to √n bins
    fd_bins = q > 0 ? max(1, ceil(Int64, maxmin_range / (2.0 * q * n^(-1 / 3)))) :
        max(1, ceil(Int64, sqrt(n)))

    # fit histogram and convert bin counts to probabilities
    h = StatsKit.fit(Histogram, s; nbins = fd_bins)
    p = h.weights ./ sum(h.weights)

    # histogram entropy in bits; eps() guards against log(0)
    ent = -sum(p .* log2.(p .+ eps()))

    # construct SampleEntropy estimator once and reuse for both sent and nsent
    se = SampleEntropy(s)

    # differential entropy
    # estimate the PDF using Kernel Density Estimation (KDE)
    kde_model = kde(s)
    # integrate over the whole KDE grid (it includes the kernel tails beyond [min, max]);
    # truncating to the data range makes the estimate depend on the signal scale
    points = collect(kde_model.x)
    pdf_values = max.(kde_model.density, 0)
    # calculate differential entropy in bits
    dent = -trapz(points, pdf_values .* log2.(pdf_values .+ eps()))

    return (
        ent = ent,
        shent = Wavelets.coefentropy(s, ShannonEntropy()),
        leent = Wavelets.coefentropy(s, LogEnergyEntropy()),
        sent = ComplexityMeasures.complexity(se, s),
        nsent = ComplexityMeasures.complexity_normalized(se, s),
        dent = dent,
    )
end

"""
    _entropy(s)

Calculate signal entropy descriptors for a 3-D signal array:

- histogram-based entropy in bits (Freedman-Diaconis binning)
- Shannon entropy (Wavelets.coefentropy)
- log energy entropy (Wavelets.coefentropy)
- sample entropy (ComplexityMeasures)
- normalised sample entropy (ComplexityMeasures)
- differential entropy

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)

# Returns

Named tuple:

- `ent::Matrix{Float64}`: entropy in bits, shape (channels, epochs)
- `shent::Matrix{Float64}`: Shanon entropy, shape (channels, epochs)
- `leent::Matrix{Float64}`: log energy entropy, shape (channels, epochs)
- `sent::Matrix{Float64}`: sample entropy, shape (channels, epochs)
- `nsent::Matrix{Float64}`: normalized sample entropy, shape (channels, epochs)
- `dent::Matrix{Float64}`: differential entropy, shape (channels, epochs)
"""
function _entropy(
    s::AbstractArray,
)::@NamedTuple{
    ent::Matrix{Float64},
    shent::Matrix{Float64},
    leent::Matrix{Float64},
    sent::Matrix{Float64},
    nsent::Matrix{Float64},
    dent::Matrix{Float64},
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
    dent = zeros(ch_n, ep_n)

    # calculate over channel and epochs
    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        entropy_data = _entropy(@view(s[ch_idx, :, ep_idx]))
        ent[ch_idx, ep_idx] = entropy_data.ent
        shent[ch_idx, ep_idx] = entropy_data.shent
        leent[ch_idx, ep_idx] = entropy_data.leent
        sent[ch_idx, ep_idx] = entropy_data.sent
        nsent[ch_idx, ep_idx] = entropy_data.nsent
        dent[ch_idx, ep_idx] = entropy_data.dent
    end

    return (; ent, shent, leent, sent, nsent, dent)
end

"""
    entropy(obj; <keyword arguments>)

Calculate signal entropy descriptors:

- histogram-based entropy in bits (Freedman-Diaconis binning)
- Shannon entropy (Wavelets.coefentropy)
- log energy entropy (Wavelets.coefentropy)
- sample entropy (ComplexityMeasures)
- normalized sample entropy (ComplexityMeasures)
- differential entropy

# Returns

Named tuple:

- `ent::Matrix{Float64}`: entropy in bits, shape (channels, epochs)
- `shent::Matrix{Float64}`: Shanon entropy, shape (channels, epochs)
- `leent::Matrix{Float64}`: log energy entropy, shape (channels, epochs)
- `sent::Matrix{Float64}`: sample entropy, shape (channels, epochs)
- `nsent::Matrix{Float64}`: normalized sample entropy, shape (channels, epochs)
- `dent::Matrix{Float64}`: differential entropy, shape (channels, epochs)
"""
function entropy(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
)::@NamedTuple{
    ent::Matrix{Float64},
    shent::Matrix{Float64},
    leent::Matrix{Float64},
    sent::Matrix{Float64},
    nsent::Matrix{Float64},
    dent::Matrix{Float64},
}
    # resolve channel names to integer indices, optionally skipping bad channels
    ch =
        exclude_bads ?
        get_channel(obj; ch = ch, exclude = "bad") :
        get_channel(obj; ch = ch, exclude = "")
    isempty(ch) && throw(ArgumentError("No channels selected."))

    return _entropy(@view(obj.data[ch, :, :]))
end

"""
    _negentropy(s; <keyword arguments>)

Calculate negentropy for a 1-D signal vector.

Negentropy measures how far a signal's distribution departs from Gaussian: `ne = 0.5·log₂(2πe·var(s)) − H(s)` (bits), where `H(s)` is the entropy selected by `type` (differential entropy by default). ne ≈ 0 for Gaussian; ne > 0 for distributions that are more structured (peaky, multi-modal, etc.).

# Arguments

- `s::AbstractVector`: signal vector
- `demean::Bool=true`: if `true`, subtract DC before calculating negentropy
- `norm::Bool=true`: if `true`, normalize the signal by its total energy
- `type::Symbol=:diff`: entropy type used for calculations (`:diff` differential, `:shannon` Shannon, `:sample` sample)

# Returns

- `Float64`: negentropy (≥ 0; equals 0 for a Gaussian signal)
"""
function _negentropy(
    s::AbstractVector;
    demean::Bool = true,
    norm::Bool = true,
    type::Symbol = :diff,
)::Float64
    # validate
    _check_var(type, [:diff, :shannon, :sample], "type")

    # work on a copy (the input, e.g. a view of the object data, is not modified)
    s = Vector{Float64}(s)
    # remove DC offset so variance reflects only signal variability
    demean && (s = _remove_dc(s))
    # normalize the signal by its total energy
    norm && (s ./= sum(abs2, s))
    # Gaussian differential entropy in bits (as the entropies below): 0.5·log₂(2πe·σ²)
    gaussian_h = 0.5 * log2(2 * π * ℯ * var(s))

    if type === :diff
        # differential entropy in bits
        signal_h = NeuroAnalyzer._entropy(s).dent
    elseif type === :shannon
        signal_h = NeuroAnalyzer._entropy(s).ent
    elseif type === :sample
        signal_h = NeuroAnalyzer._entropy(s).sent
    end

    return gaussian_h - signal_h
end

"""
    _negentropy(s; <keyword arguments>)

Calculate negentropy for a 3-D signal array.

Negentropy measures how far a signal's distribution departs from Gaussian: `ne = 0.5·log₂(2πe·var(s)) − H(s)` (bits), where `H(s)` is the entropy selected by `type` (differential entropy by default). ne ≈ 0 for Gaussian; ne > 0 for distributions that are more structured (peaky, multi-modal, etc.).

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `demean::Bool=true`: if `true`, subtract DC before calculating negentropy
- `norm::Bool=true`: if `true`, normalize the signal by its total energy
- `type::Symbol=:diff`: entropy type used for calculations (`:diff` differential, `:shannon` Shannon, `:sample` sample)

# Returns

- `Matrix{Float64}`: negentropy (≥ 0; equals 0 for a Gaussian signal), shape (channel, epochs)
"""
function _negentropy(
    s::AbstractArray;
    demean::Bool = true,
    norm::Bool = true,
    type::Symbol = :diff,
)::Matrix{Float64}
    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate output
    ne = zeros(ch_n, ep_n)

    # initialize progress bar
    progbar = Progress(
        ch_n * ep_n;
        dt = 1,
        barlen = 20,
        color = :white,
        enabled = progress_bar,
    )

    # calculate over channel and epochs
    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        ne[ch_idx, ep_idx] = _negentropy(
            @view(s[ch_idx, :, ep_idx]),
            demean = demean,
            norm = norm,
            type = type,
        )
        # update progress bar
        progress_bar && next!(progbar)
    end
    progress_bar && finish!(progbar)

    return ne
end

"""
    negentropy(obj; <keyword arguments>)

Calculate negentropy for a NEURO object.

Negentropy measures how far a signal's distribution departs from Gaussian: `ne = 0.5·log₂(2πe·var(s)) − H(s)` (bits), where `H(s)` is the entropy selected by `type` (differential entropy by default). ne ≈ 0 for Gaussian; ne > 0 for distributions that are more structured (peaky, multi-modal, etc.).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `demean::Bool=true`: if `true`, subtract DC before calculating negentropy
- `norm::Bool=true`: if `true`, normalize the signal by its total energy
- `type::Symbol=:diff`: entropy type used for calculations (`:diff` differential, `:shannon` Shannon, `:sample` sample)

# Returns

- `Matrix{Float64}`: negentropy (≥ 0; equals 0 for a Gaussian signal), shape (channel, epochs)
"""
function negentropy(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    demean::Bool = true,
    norm::Bool = true,
    type::Symbol = :diff,
)::Matrix{Float64}
    # resolve channel names to integer indices, optionally skipping bad channels
    ch =
        exclude_bads ?
        get_channel(obj; ch = ch, exclude = "bad") :
        get_channel(obj; ch = ch, exclude = "")
    isempty(ch) && throw(ArgumentError("No channels selected."))

    return _negentropy(
        @view(obj.data[ch, :, :]);
        demean = demean,
        norm = norm,
        type = type,
    )
end
