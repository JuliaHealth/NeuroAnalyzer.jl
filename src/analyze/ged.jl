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
function ged(s1::AbstractArray, s2::AbstractArray)::NamedTuple{(:sged, :ress, :ress_norm), Tuple{Matrix{Float64}, Vector{Float64}, Vector{Float64}}}

    @assert size(s1) == size(s2) "s1 and s2 must have the same size."
    _chk3d(s1)
    _chk3d(s2)

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
    ged(obj1, obj2; <keyword arguments>)

Perform generalized eigendecomposition.

# Arguments

- `obj1::NeuroAnalyzer.NEURO`: signal data to be analyzed
- `obj2::NeuroAnalyzer.NEURO`: original signal data
- `ch1::Union{String, Vector{String}}: list of channels
- `ch2::Union{String, Vector{String}}: list of channels
- `ep1::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj1))`: default use all epochs
- `ep2::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj2))`: default use all epochs

# Returns

Named tuple containing:
- `sged::Array{Float64, 3}`
- `ress::Matrix{Float64}`
- `ress_norm::Matrix{Float64}`
"""
function ged(obj1::NeuroAnalyzer.NEURO, obj2::NeuroAnalyzer.NEURO; ch1::Union{String, Vector{String}}, ch2::Union{String, Vector{String}}, ep1::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj1)), ep2::Union{Int64, Vector{Int64}, <:AbstractRange}=_c(nepochs(obj2)))::NamedTuple{(:sged, :ress, :ress_norm), Tuple{Array{Float64, 3}, Matrix{Float64}, Matrix{Float64}}}

    ch1 = get_channel(obj1, ch=ch1)
    ch2 = get_channel(obj2, ch=ch2)
    @assert length(ch1) == length(ch2) "ch1 and ch2 must have the same length."

    _check_epochs(obj1, ep1)
    _check_epochs(obj2, ep2)
    @assert length(ep1) == length(ep2) "ep1 and ep2 must have the same length."
    @assert epoch_len(obj1) == epoch_len(obj2) "OBJ1 and OBJ2 must have the same epoch lengths."
    length(ep1) == 1 && (ep1 = [ep1])
    length(ep2) == 1 && (ep2 = [ep2])

    ep_n = length(ep1)
    ch_n = length(ch1)

    sged = zeros(ch_n, epoch_len(obj1), ep_n)
    ress = zeros(ch_n, ep_n)
    ress_norm = zeros(ch_n, ep_n)

    Threads.@threads for ep_idx in 1:ep_n
        @inbounds sged[:, :, ep_idx], ress[:, ep_idx], ress_norm[:, ep_idx] = @views ged(obj1.data[ch1, :, ep1[ep_idx]], obj2.data[ch2, :, ep2[ep_idx]])
    end

    return (sged=sged, ress=ress, ress_norm=ress_norm)

end
