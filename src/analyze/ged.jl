export ged

"""
    ged(s1, s2)

Perform generalized eigendecomposition.

# Arguments

- `s1::AbstractArray`: signal to be analyzed
- `s2::AbstractArray`: original signal

# Returns

Named tuple containing:
- `sged::Matrix{Float64}`
- `ress::Vector{Float64}`
- `ress_norm::Vector{Float64}`: RESS normalized to -1..1
"""
function ged(s1::AbstractArray, s2::AbstractArray)

    size(s1) == size(s2) || throw(ArgumentError("s1 and s2 must have the same size."))

    s1cov = cov(s1')
    s2cov = cov(s2')

    eig_val, eig_vec = eigen(s1cov, s2cov)
    eig_val_idx = sortperm(eig_val, rev=true)
    eig_val = eig_val[eig_val_idx]
    eig_vec = m_sort(eig_vec, eig_val_idx, dims=2)

    sged = s2 .* eig_vec[:, 1]
    ress = pinv(eig_vec[:, 1]')
    ress_norm = ress ./ maximum(abs.(ress))

    return (sged=sged, ress=ress, ress_norm=ress_norm)

end

"""
    ged(obj1, obj2; ch1, ch2, ep1, ep2)

Perform generalized eigendecomposition.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: signal data to be analyzed
- `obj2::NeuroAnalyzer.NEURO`: original signal data
- `ch1::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj1)`: index of channels, default is all signal channels
- `ch2::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj2)`: index of channels, default is all signal channels
- `ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1))`: default use all epochs
- `ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2))`: default use all epochs

# Returns

- `sged::Array{Float64, 3}`
- `ress::Matrix{Float64}`
- `ress_norm::Matrix{Float64}`
"""
function ged(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch1::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj1), ch2::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj2), ep1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1)), ep2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2)))

    _check_channels(obj1, ch1)
    _check_channels(obj2, ch2)
    length(ch1) == length(ch2) || throw(ArgumentError("ch1 and ch2 must have the same length."))
    
    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    length(ep1) == length(ep2) || throw(ArgumentError("ep1 and ep2 must have the same length."))
    epoch_len(obj1) == epoch_len(obj2) || throw(ArgumentError("OBJ1 and OBJ2 epoch lengths must be equal."))

    ep_n = length(ep1)
    ch_n = length(ch1)

    sged = zeros(ch_n, epoch_len(obj1), ep_n)
    ress = zeros(ch_n, ep_n)
    ress_norm = zeros(ch_n, ep_n)

    Threads.@threads for ep_idx in 1:ep_n
        sged[:, :, ep_idx], ress[:, ep_idx], ress_norm[:, ep_idx] = @views ged(obj1.data[ch1, :, ep1[ep_idx]], obj2.data[ch2, :, ep2[ep_idx]])
    end

    return (sged=sged, ress=ress, ress_norm=ress_norm)

end
