export ged

"""
    ged(signal1, signal2)

Perform generalized eigendecomposition.

# Arguments

- `signal1::AbstractArray`: signal to be analyzed
- `signal2::AbstractArray`: original signal

# Returns

Named tuple containing:
- `sged::Matrix{Float64}`
- `ress::Vector{Float64}`
- `resnormalized::Vector{Float64}`: RESS normalized to -1..1
"""
function ged(signal1::AbstractArray, signal2::AbstractArray)

    size(signal1) == size(signal2) || throw(ArgumentError("signal1 and signal2 must have the same size."))

    s1cov = cov(signal1')
    s2cov = cov(signal2')
    eig_val, eig_vec = eigen(s1cov, s2cov)
    eig_val_idx = sortperm(eig_val, rev=true)
    eig_val = eig_val[eig_val_idx]
    eig_vec = m_sort(eig_vec, eig_val_idx, dims=2)
    sged = signal2 .* eig_vec[:, 1]
    ress = pinv(eig_vec[:, 1]')
    resnormalized = ress ./ maximum(abs.(ress))

    return (sged=sged, ress=ress, resnormalized=resnormalized)
end

"""
    ged(obj1, obj2; channel1, channel2, epoch1, epoch2)

Perform generalized eigendecomposition.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: signal data to be analyzed
- `obj2::NeuroAnalyzer.NEURO`: original signal data
- `channel1::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj1)`: index of channels, default is all signal channels
- `channel2::Union{Int64, Vector{Int64}, AbstractRange}=signal_channels(obj2)`: index of channels, default is all signal channels
- `epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1))`: default use all epochs
- `epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2))`: default use all epochs

# Returns

- `sged::Array{Float64, 3}`
- `ress::Matrix{Float64}`
- `resnormalized::Matrix{Float64}`
"""
function ged(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; channel1::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj1, type=Symbol(obj1.header.recording[:data_type])), channel2::Union{Int64, Vector{Int64}, AbstractRange}=get_channel_bytype(obj2, type=Symbol(obj2.header.recording[:data_type])), epoch1::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj1)), epoch2::Union{Int64, Vector{Int64}, AbstractRange}=_c(epoch_n(obj2)))

    _check_channels(obj1, channel1)
    _check_channels(obj2, channel2)
    length(channel1) == length(channel2) || throw(ArgumentError("channel1 and channel2 lengths must be equal."))
    
    _check_epochs(obj1, epoch1)
    _check_epochs(obj2, epoch2)
    length(epoch1) == length(epoch2) || throw(ArgumentError("epoch1 and epoch2 lengths must be equal."))
    epoch_len(obj1) == epoch_len(obj2) || throw(ArgumentError("OBJ1 and OBJ2 epoch lengths must be equal."))

    ep_n = length(epoch1)
    ch_n = length(channel1)

    sged = zeros(ch_n, epoch_len(obj1), ep_n)
    ress = zeros(ch_n, ep_n)
    resnormalized = zeros(ch_n, ep_n)

    Threads.@threads for ep_idx in 1:ep_n
        sged[:, :, ep_idx], ress[:, ep_idx], resnormalized[:, ep_idx] = @views ged(obj1.data[channel1, :, epoch1[ep_idx]], obj2.data[channel2, :, epoch2[ep_idx]])
    end

    return (sged=sged, ress=ress, resnormalized=resnormalized)
end
