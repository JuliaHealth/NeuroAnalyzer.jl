export pca_decompose
export pca_reconstruct
export pca_reconstruct!

"""
    pca_decompose(s, n)

Calculate `n` first Primary Components (PCs).

# Arguments

- `s::AbstractArray`
- `n::Int64`: number of PCs

# Returns

Named tuple containing:
- `pc::Array{Float64, 3}:`: PC(1)..PC(n) × epoch
- `pcv::Matrix{Float64}`: variance of PC(1)..PC(n) × epoch (% of total variance)
- `pcm::PCA{Float64}`: PC mean
- `pc_model::MultivariateStats.PCA{Float64}`: PC model
"""
function pca_decompose(s::AbstractArray; n::Int64)

    @assert n >= 1 "n must be ≥ 1."
    @assert n <= size(s, 1) "n must be ≤ $(size(s, 1))."

    ep_n = size(s, 3)

    # check maximum n
    pc_tmp = []
    n_tmp = n
    @inbounds for ep_idx in 1:ep_n
        pc_tmp = @views MultivariateStats.fit(PCA, s[:, :, ep_idx], maxoutdim=n)
        size(pc_tmp)[2] < n_tmp && (n_tmp = size(pc_tmp)[2])
    end
    (n_tmp < n && verbose) && _warn("Only $n_tmp PCs were generated.")
    n = n_tmp

    pc = zeros(n, size(s, 2), ep_n)
    pcv = zeros(n, ep_n)
    pc_model = nothing

    @inbounds for ep_idx in 1:ep_n
        # m_cov = s_cov(s)
        # eig_val, eig_vec = eigen(m_cov)
        # eig_val_idx = sortperm(eig_val, rev=true)
        # eig_val = eig_val[eig_val_idx]
        # eig_vec = m_sort(eig_vec, eig_val_idx)
        # eig_val = 100 .* eig_val / sum(eig_val) # convert to %

        pc_model = @views MultivariateStats.fit(PCA, s[:, :, ep_idx], maxoutdim=n)
        v = MultivariateStats.principalvars(pc_model) ./ MultivariateStats.var(pc_model) * 100

        for idx in 1:n
            pcv[idx, ep_idx] = v[idx]
            # pc[idx, :, ep_idx] = (eig_vec[:, idx] .* s)[idx, :]
            pc[idx, :, ep_idx] = @views MultivariateStats.predict(pc_model, s[:, :, ep_idx])[idx, :]
        end
    end

    return (pc=pc, pcv=pcv, pcm=pc_model.mean, pc_model=pc_model)
end

"""
    pca_decompose(obj; ch, n)

Calculate `n` first Primary Components (PCs).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
- `n::Int64`: number of PCs to calculate

# Returns

Named tuple containing:
- `pc::Array{Float64, 3}:`: PC(1)..PC(n) × epoch
- `pcv::Matrix{Float64}`: variance of PC(1)..PC(n) × epoch (% of total variance)
- `pcm::Vector{Float64}`: PC mean
- `pc_model::MultivariateStats.PCA{Float64}`: PC model
"""
function pca_decompose(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj), n::Int64)

    _check_channels(obj, ch)
    isa(ch, Int64) && (ch = [ch])

    pc, pcv, pcm, pc_model = @views pca_decompose(obj.data[ch, :, :], n=n)

    return (pc=pc, pcv=pcv, pcm=pcm, pc_model=pc_model)
end

"""
    pca_reconstruct(s, pc, pca)

Reconstructs signal using PCA components.

# Arguments

- `s::AbstractArray`
- `pc::AbstractArray:`: IC(1)..IC(n) × epoch
- `pc_model::MultivariateStats.PCA{Float64}:`: PC model

# Returns

- `s_new::Array{Float64, 3}`
"""
function pca_reconstruct(s::AbstractArray; pc::AbstractArray, pc_model::MultivariateStats.PCA{Float64})

    s_new = similar(s)
    ep_n = size(s, 3)

    @inbounds for ep_idx in 1:ep_n
        s_new[:, :, ep_idx] = @views MultivariateStats.reconstruct(pc_model, pc[:, :, ep_idx])
    end

    return s_new
end

"""
    pca_reconstruct(obj; ch)

Reconstruct signal using embedded PCA components (`:pc`) and model (`:pc_model`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function pca_reconstruct(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj))

    @assert :pc in keys(obj.components) "OBJ does not contain :pc component. Perform pca_decompose() first."
    @assert :pc_model in keys(obj.components) "OBJ does not contain :pc_model component. Perform pca_decompose() first."

    _check_channels(obj, ch)
    isa(ch, Int64) && (ch = [ch])

    obj_new = deepcopy(obj)

    obj_new.data[ch, :, :] = @views pca_reconstruct(obj_new.data[ch, :, :], pc=obj_new.components[:pc], pc_model=obj_new.components[:pc_model])

    reset_components!(obj_new)
    push!(obj_new.history, "pca_reconstruct(OBJ, ch=$ch)")

    return obj_new

end

"""
    pca_reconstruct!(obj; ch)

Reconstruct signal using embedded PCA components (`:pc`) and model (`:pc_model`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
"""
function pca_reconstruct!(obj::NeuroAnalyzer.NEURO; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj))

    obj_new = pca_reconstruct(obj, ch=ch)
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components

    return nothing

end

"""
    pca_reconstruct(obj, pc, pc_model; ch)

Reconstruct signal using external PCA components (`pc` and `pca`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `pc::Array{Float64, 3}:`: PC(1)..PC(n) × epoch
- `pc_model::MultivariateStats.PCA{Float64}`: PC model
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels

# Returns

- `obj_new::NeuroAnalyzer.NEURO`
"""
function pca_reconstruct(obj::NeuroAnalyzer.NEURO, pc::Array{Float64, 3}, pc_model::MultivariateStats.PCA{Float64}; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj))

    _check_channels(obj, ch)
    isa(ch, Int64) && (ch = [ch])

    obj_new = deepcopy(obj)

    obj_new.data[ch, :, :] = @views pca_reconstruct(obj_new.data[ch, :, :], pc=pc, pc_model=pc_model)

    reset_components!(obj_new)
    push!(obj_new.history, "pca_reconstruct(OBJ, ch=$ch, pc, pc_model)")

    return obj_new

end

"""
    pca_reconstruct!(obj, pc, pc_model; ch)

Reconstruct signals using external PCA components (`pc` and `pc_model`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `pc::Array{Float64, 3}:`: PC(1)..PC(n) × epoch
- `pc_model::MultivariateStats.PCA{Float64}`: PC model
- `ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj)`: index of channels, default is all signal channels
"""
function pca_reconstruct!(obj::NeuroAnalyzer.NEURO, pc::Array{Float64, 3}, pc_model::MultivariateStats.PCA{Float64}; ch::Union{Int64, Vector{Int64}, <:AbstractRange}=signal_channels(obj))

    obj_new = pca_reconstruct(obj, pc, pc_model, ch=ch)
    obj.data = obj_new.data
    obj.history = obj_new.history
    obj.components = obj_new.components

    return nothing

end
