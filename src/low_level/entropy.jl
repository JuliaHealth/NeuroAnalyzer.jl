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

