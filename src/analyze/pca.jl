export pca
export pca_reconstruct
export pca_reconstruct!

"""
    pca(signal, n)

Calculate `n` first PCs.

# Arguments

- `signal::AbstractArray`
- `n::Int64`: number of PCs

# Returns

Named tuple containing:
- `pc::Array{Float64, 3}:`: PC(1)..PC(n) × epoch
- `pc_var::Matrix{Float64}`: variance of PC(1)..PC(n) × epoch
- `pc_m::PCA{Float64}`: PC mean
- `pca::MultivariateStats.PCA{Float64}`: PC model
"""
function pca(signal::AbstractArray; n::Int64)

    n < 0 && throw(ArgumentError("n must be ≥ 1."))
    n > size(signal, 1) && throw(ArgumentError("n must be ≤ $(size(signal, 1))."))

    ep_n = size(signal, 3)
    pc_m = []

    # check maximum n
    n_tmp = n
    @inbounds @simd for ep_idx in 1:ep_n
        pc_m = @views MultivariateStats.fit(PCA, signal[:, :, ep_idx], maxoutdim=n)
        size(pc_m)[2] < n_tmp && (n_tmp = size(pc_m)[2])
    end
    (n_tmp < n && verbose == true) && @info "Only $n_tmp PC components were generated."
    n = n_tmp
    
    pc = zeros(n, size(signal, 2), ep_n)
    pc_var = zeros(n, ep_n)
    pca = nothing
    
    @inbounds @simd for ep_idx in 1:ep_n
        # m_cov = s_cov(s)
        # eig_val, eig_vec = eigen(m_cov)
        # eig_val_idx = sortperm(eig_val, rev=true)
        # eig_val = eig_val[eig_val_idx]
        # eig_vec = m_sort(eig_vec, eig_val_idx)
        # eig_val = 100 .* eig_val / sum(eig_val) # convert to %

        pca = @views MultivariateStats.fit(PCA, signal[:, :, ep_idx], maxoutdim=n)
        v = MultivariateStats.principalvars(pca) ./ MultivariateStats.var(pca) * 100

        for idx in 1:n
            pc_var[idx, ep_idx] = v[idx]
            # pc[idx, :, ep_idx] = (eig_vec[:, idx] .* s)[idx, :]
            pc[idx, :, ep_idx] = @views MultivariateStats.predict(pca, signal[:, :, ep_idx])[idx, :]
        end
    end

    return (pc=pc, pc_var=pc_var, pc_m=pca.mean, pca=pca)
end

"""
    pca_reconstruct(signal, pc, pca)

Reconstructs `signal` using PCA components.

# Arguments

- `signal::AbstractArray`
- `pc::AbstractArray:`: IC(1)..IC(n) × epoch
- `pca::MultivariateStats.PCA{Float64}:`: IC(1)..IC(n) × epoch

# Returns

- `s_reconstructed::Array{Float64, 3}`
"""
function pca_reconstruct(signal::AbstractArray; pc::AbstractArray, pca::MultivariateStats.PCA{Float64})

    s_reconstructed = similar(signal)
    ep_n = size(signal, 3)
    @inbounds @simd for ep_idx in 1:ep_n
        s_reconstructed[:, :, ep_idx] = @views MultivariateStats.reconstruct(pca, pc[:, :, ep_idx])
    end

    return s_reconstructed
end

"""
    pca(obj; channel, n)

Perform principal component analysis (PCA).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
- `n::Int64`: number of PCs to calculate

# Returns

Named tuple containing:
- `pc::Array{Float64, 3}:`: PC(1)..PC(n) × epoch
- `pc_var::Matrix{Float64}`: variance of PC(1)..PC(n) × epoch
- `pc_m::Vector{Float64}`: PC mean
- `pca::MultivariateStats.PCA{Float64}`: PC model
"""
function pca(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj)), n::Int64)

    _check_channels(obj, channel)

    pc, pc_var, pc_m, pca = @views pca(obj.data[channel, :, :], n=n)

    return (pc=pc, pc_var=pc_var, pc_m=pc_m, pca=pca)
end

"""
    pca_reconstruct(obj; channel)

Reconstruct signals using embedded PCA components (`:pc`) and model (`:pca`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function pca_reconstruct(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj)))

    :pc in obj.header.component_names || throw(ArgumentError("OBJ does not contain :pc component. Perform pca(OBJ) first."))
    :pca in obj.header.component_names || throw(ArgumentError("OBJ does not contain :pca component. Perform pca(OBJ) first."))

    _check_channels(obj, channel)

    obj_new = deepcopy(obj)
    pc_idx = component_idx(obj, c=:pc)
    pc_m_idx = component_idx(obj, c=:pca)
    obj_new.data[channel, :, :] = @views pca_reconstruct(obj_new.data[channel, :, :], pc=obj_new.components[pc_idx], pca=obj_new.components[pc_m_idx])
    reset_components!(obj_new)
    push!(obj_new.header.history, "pca_reconstruct(OBJ, channel=$channel)")

    return obj_new
end

"""
    pca_reconstruct!(obj; channel)

Reconstruct signals using embedded PCA components (`:pc`) and model (`:pca`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
"""
function pca_reconstruct!(obj::NeuroAnalyzer.NEURO; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj)))

    obj_tmp = pca_reconstruct(obj, channel=channel)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    obj.components = obj_tmp.components

    return nothing
end

"""
    pca_reconstruct(obj, pc, pca; channel)

Reconstruct signals using external PCA components (`pc` and `pca`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `pc::Array{Float64, 3}:`: PC(1)..PC(n) × epoch
- `pca::MultivariateStats.PCA{Float64}`: PC model
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels

# Returns

- `obj::NeuroAnalyzer.NEURO`
"""
function pca_reconstruct(obj::NeuroAnalyzer.NEURO, pc::Array{Float64, 3}, pca::MultivariateStats.PCA{Float64}; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj)))

    _check_channels(obj, channel)

    obj_new = deepcopy(obj)
    obj_new.data[channel, :, :] = @views pca_reconstruct(obj_new.data[channel, :, :], pc=pc, pca=pca)
    reset_components!(obj_new)
    push!(obj_new.header.history, "pca_reconstruct(OBJ, channel=$channel, pc=$pc, pca=$pca)")

    return obj_new
end

"""
    pca_reconstruct!(obj, pc, pca; channel)

Reconstruct signals using external PCA components (`pc` and `pca`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `pc::Array{Float64, 3}:`: PC(1)..PC(n) × epoch
- `pca::MultivariateStats.PCA{Float64}`: PC model
- `channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj))`: index of channels, default is all channels
"""
function pca_reconstruct!(obj::NeuroAnalyzer.NEURO, pc::Array{Float64, 3}, pca::MultivariateStats.PCA{Float64}; channel::Union{Int64, Vector{Int64}, AbstractRange}=_c(channel_n(obj)))

    tmp = pca_reconstruct(obj, pc, pca, channel=channel)
    obj.data = obj_tmp.data
    obj.header = obj_tmp.header
    reset_components!(obj)

    return nothing
end
