export ged

"""
    ged(s1, s2)

Perform Generalized Eigendecomposition.

Solves the generalised eigenvalue problem  S1·W = S2·W·Λ, where S1 and S2 are the covariance matrices of the "target" and "reference" signals.

The first (largest) eigenvector defines the spatial filter that maximally distinguishes S1 from S2. The RESS (Rhythmic Entrainment Source Separation) spatial filter is the pseudoinverse of that eigenvector.

# Arguments

- `s1::AbstractMatrix`: target signal (channels × samples)
- `s2::AbstractMatrix`: reference signal (channels × samples)

# Returns

Named tuple containing:

- `sged::Matrix{Float64}`: reference signal weighted by the first eigenvector, shape `(channels, samples)`
- `ress::Vector{Float64}`: RESS spatial filter (pseudoinverse of leading eigenvector)
- `ress_norm::Vector{Float64}`: RESS normalized to −1..1
"""
function ged(
    s1::AbstractMatrix,
    s2::AbstractMatrix
)::@NamedTuple{
    sged::Matrix{Float64},
    ress::Vector{Float64},
    ress_norm::Vector{Float64}
}

    @assert size(s1) == size(s2) "s1 and s2 must have the same size."

    # compute channels × channels covariance matrices
    # cov() expects observations in rows, so we transpose (channels × samples → samples × channels)
    s1cov = cov(s1')
    s2cov = cov(s2')

    # solve the generalised eigenvalue problem S1·W = S2·W·Λ
    eig_val, eig_vec = eigen(s1cov, s2cov)

    # sort by descending eigenvalue so the first component has maximum power ratio
    eig_val_idx = sortperm(eig_val, rev = true)
    eig_val = eig_val[eig_val_idx]
    eig_vec = m_sort(eig_vec, eig_val_idx, dims = 2)

    # weight each channel of s2 by its component in the leading eigenvector
    sged = s2 .* eig_vec[:, 1]

    # RESS spatial filter: pseudoinverse of the leading eigenvector (as a row vector)
    # pinv of a 1×ch row vector returns a ch×1 matrix; vec() converts to Vector
    ress = vec(pinv(eig_vec[:, 1]'))

    # normalise RESS to [−1, 1] by dividing by the largest absolute value
    ress_norm = ress ./ maximum(abs, ress)

    return (sged = sged, ress = ress, ress_norm = ress_norm)

end

"""
    ged(obj1, obj2; <keyword arguments>)

Perform Generalized Eigendecomposition.

Solves the generalised eigenvalue problem  S1·W = S2·W·Λ, where S1 and S2 are the covariance matrices of the "target" and "reference" signals.

The first (largest) eigenvector defines the spatial filter that maximally distinguishes S1 from S2. The RESS (Rhythmic Entrainment Source Separation) spatial filter is the pseudoinverse of that eigenvector.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: object to be analyzed (target)
- `obj2::NeuroAnalyzer.NEURO`: original object (reference)
- `ch1::Union{String, Vector{String}, Regex}`: channel name(s)
- `ch2::Union{String, Vector{String}, Regex}`: channel name(s)
- `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj1))`: epoch number(s)
- `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(nepochs(obj2))`: epoch number(s)

# Returns

Named tuple containing:

- `sged::Array{Float64, 3}`: GED output, shape `(channels, samples, epochs)`
- `ress::Matrix{Float64}`: RESS spatial filter, shape `(channels, epochs)`
- `ress_norm::Matrix{Float64}`: RESS normalized to −1..1, shape `(channels, epochs)`
"""
function ged(
    obj1::NeuroAnalyzer.NEURO,
    obj2::NeuroAnalyzer.NEURO;
    ch1::Union{String, Vector{String}, Regex},
    ch2::Union{String, Vector{String}, Regex},
    ep1::Union{Int64, Vector{Int64}, AbstractRange} = _c(nepochs(obj1)),
    ep2::Union{Int64, Vector{Int64}, AbstractRange} = _c(nepochs(obj2)),
)::@NamedTuple{
    sged::Array{Float64, 3},
    ress::Matrix{Float64},
    ress_norm::Matrix{Float64}
}

    # resolve channel names to integer indices, optionally skipping bad channels
    ch1 = exclude_bads ? get_channel(obj1, ch = ch1, exclude = "bad") : get_channel(obj1, ch = ch1, exclude = "")
    ch2 = exclude_bads ? get_channel(obj2, ch = ch2, exclude = "bad") : get_channel(obj2, ch = ch2, exclude = "")
    @assert length(ch1) == length(ch2) "Lengths of ch1 ($(length(ch1))) and ch2 ($(length(ch2))) must be equal."

    # validate epoch indices and ensure both objects have matching epoch structure
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    # normalize scalar epoch arguments to vectors so indexing is uniform
    isa(ep1, Int64) && (ep1 = [ep1])
    isa(ep2, Int64) && (ep2 = [ep2])
    @assert length(ep1) == length(ep2) "Lengths of ep1 ($(length(ep1))) and ep2 ($(length(ep2))) must be equal."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."

    # number of channels
    ch_n = length(ch1)
    # number of epochs
    ep_n = length(ep1)
    # epoch length
    ep_len = epoch_len(obj1)

    # pre-allocate output
    sged = zeros(ch_n, ep_len, ep_n)
    ress = zeros(ch_n, ep_n)
    ress_norm = zeros(ch_n, ep_n)

    @inbounds Threads.@threads :dynamic for ep_idx in 1:ep_n
        ged_data = ged(
            @view(obj1.data[ch1, :, ep1[ep_idx]]),
            @view(obj2.data[ch2, :, ep2[ep_idx]]),
        )
        sged[:, :, ep_idx] = ged_data.sged
        ress[:, ep_idx] = ged_data.ress
        ress_norm[:, ep_idx] = ged_data.ress_norm
    end

    return (sged = sged, ress = ress, ress_norm = ress_norm)

end
