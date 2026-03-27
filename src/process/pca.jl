export pca_decompose
export pca_reconstruct
export pca_reconstruct!

"""
    pca_decompose(s, n)

Calculate `n` first Primary Components (PCs) for a 3-D signal array.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `n::Int64`: number of PCs

# Returns

Named tuple:

- `pc::Array{Float64, 3}`: PC(1)..PC(n) × epoch
- `pcv::Matrix{Float64}`: PC(1)..PC(n) variances (fraction of total variance explained)
- `pcm::Vector{Float64}`: PC means
- `pc_model::MultivariateStats.PCA{Float64}`: PC model
"""
function pca_decompose(
    s::AbstractArray;
    n::Int64,
)::@NamedTuple{
    pc::Array{Float64, 3},
    pcv::Matrix{Float64},
    pcm::Vector{Float64},
    pc_model::MultivariateStats.PCA{Float64},
}
    _chk3d(s)
    n >= 1 || throw(ArgumentError("n must be ≥ 1."))
    n <= size(s, 1) || throw(ArgumentError("n must be ≤ $(size(s, 1))."))

    ep_n = size(s, 3)

    # check maximum n
    pc_tmp = []
    n_tmp = n
    @inbounds for ep_idx in 1:ep_n
        pc_tmp =
            MultivariateStats.fit(PCA, @view(s[:, :, ep_idx]), maxoutdim = n, pratio = 1)
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

        pc_model =
            MultivariateStats.fit(PCA, @view(s[:, :, ep_idx]), maxoutdim = n, pratio = 1)
        v =
            MultivariateStats.principalvars(pc_model) ./ MultivariateStats.var(pc_model) *
            100

        for idx in 1:n
            pcv[idx, ep_idx] = v[idx]
            # pc[idx, :, ep_idx] = (eig_vec[:, idx] .* s)[idx, :]
            pc[idx, :, ep_idx] =
                MultivariateStats.predict(pc_model, @view(s[:, :, ep_idx]))[idx, :]
        end
    end

    pcm = pc_model.mean

    return (; pc, pcv, pcm, pc_model)
end

"""
    pca_decompose(obj; <keyword arguments>)

Calculate `n` first Primary Components (PCs).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `n::Int64`: number of PCs to calculate

# Returns

Named tuple:

- `pc::Array{Float64, 3}`: PC(1)..PC(n) × epoch
- `pcv::Matrix{Float64}`: PC variances (fraction of total variance explained)
- `pcm::Vector{Float64}`: PC means
- `pc_model::MultivariateStats.PCA{Float64}`: PC model
"""
function pca_decompose(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    n::Int64,
)::@NamedTuple{
    pc::Array{Float64, 3},
    pcv::Matrix{Float64},
    pcm::Vector{Float64},
    pc_model::MultivariateStats.PCA{Float64},
}

    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)

    return pca_decompose(@view(obj.data[ch, :, :]); n = n)
end

"""
    pca_reconstruct(s, pc, pca)

Reconstructs signal using PCA components for a 3-D signal array.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `pc::AbstractArray`: IC(1)..IC(n) × epoch
- `pc_model::MultivariateStats.PCA{Float64}`: PC model

# Returns

- `Array{Float64, 3}`
"""
function pca_reconstruct(
    s::AbstractArray;
    pc::AbstractArray,
    pc_model::MultivariateStats.PCA{Float64},
)::Array{Float64, 3}
    _chk3d(s)
    s_new = similar(s, Float64)
    ep_n = size(s, 3)

    @inbounds for ep_idx in 1:ep_n
        s_new[:, :, ep_idx] =
            MultivariateStats.reconstruct(pc_model, @view(pc[:, :, ep_idx]))
    end

    return s_new
end

"""
    pca_reconstruct(obj, pc, pc_model; <keyword arguments>)

Reconstruct signal using PCA components (`pc` and `pca`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `pc::Array{Float64, 3}`: PC(1)..PC(n) × epoch
- `pc_model::MultivariateStats.PCA{Float64}`: PC model

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function pca_reconstruct(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    pc::Array{Float64, 3},
    pc_model::MultivariateStats.PCA{Float64},
)::NeuroAnalyzer.NEURO

    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)

    # create new dataset
    obj_new = deepcopy(obj)

    obj_new.data[ch, :, :] =
        pca_reconstruct(@view(obj_new.data[ch, :, :]); pc = pc, pc_model = pc_model)
    push!(obj_new.history, "pca_reconstruct(obj; ch=$ch)")

    return obj_new
end

"""
    pca_reconstruct!(obj, pc, pc_model; <keyword arguments>)

Reconstruct signals using PCA components (`pc` and `pc_model`).

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `pc::Array{Float64, 3}`: PC(1)..PC(n) × epoch
- `pc_model::MultivariateStats.PCA{Float64}`: PC model

# Returns

- `Nothing`
"""
function pca_reconstruct!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    pc::Array{Float64, 3},
    pc_model::MultivariateStats.PCA{Float64},
)::Nothing
    obj_new = pca_reconstruct(obj; ch = ch, pc = pc, pc_model = pc_model)
    obj.data = obj_new.data
    obj.history = obj_new.history

    return nothing
end
