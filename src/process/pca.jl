export pca_decompose
export pca_reconstruct
export pca_reconstruct!

"""
    _pca_decompose(s, n)

Calculate `n` first Principal Components (PCs) for a 3-D signal array.

One PCA model is fitted to all epochs (channels × samples of all epochs), so the components and the model are shared by all epochs and `_pca_reconstruct()` can reconstruct every epoch.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `n::Int64`: number of PCs

# Returns

Named tuple:

- `pc::Array{Float64, 3}`: PC(1)..PC(n) × samples × epochs
- `pcv::Matrix{Float64}`: variance of PC(1)..PC(n) in each epoch (% of the total variance of the epoch)
- `pcm::Vector{Float64}`: channel means of the model
- `pc_model::MultivariateStats.PCA{Float64}`: PC model
"""
function _pca_decompose(
    s::AbstractArray;
    n::Int64,
)::@NamedTuple{
    pc::Array{Float64, 3},
    pcv::Matrix{Float64},
    pcm::Vector{Float64},
    pc_model::MultivariateStats.PCA{Float64},
}
    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    n >= 1 || throw(ArgumentError("n must be ≥ 1."))
    n <= size(s, 1) || throw(ArgumentError("n must be ≤ $(size(s, 1))."))

    ch_n, ep_len, ep_n = size(s)

    # one model for all epochs
    pc_model = MultivariateStats.fit(
        PCA,
        reshape(Float64.(s), ch_n, ep_len * ep_n);
        maxoutdim = n,
        pratio = 1,
    )
    n_pc = size(pc_model, 2)
    n_pc < n && _warn("Only $n_pc PCs were generated.")

    pc = zeros(n_pc, ep_len, ep_n)
    pcv = zeros(n_pc, ep_n)
    @inbounds for ep_idx = 1:ep_n
        pc[:, :, ep_idx] = MultivariateStats.predict(pc_model, @view(s[:, :, ep_idx]))
        # total variance of the epoch (sum of channel variances)
        v_tot = sum(var(@view(s[ch_idx, :, ep_idx])) for ch_idx = 1:ch_n)
        for idx = 1:n_pc
            pcv[idx, ep_idx] = 100 * var(@view(pc[idx, :, ep_idx])) / v_tot
        end
    end

    return (; pc, pcv, pcm = pc_model.mean, pc_model)
end

"""
    pca_decompose(obj; <keyword arguments>)

Calculate `n` first Principal Components (PCs); one model is fitted to all epochs.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `n::Int64`: number of PCs to calculate

# Returns

Named tuple:

- `pc::Array{Float64, 3}`: PC(1)..PC(n) × epoch
- `pcv::Matrix{Float64}`: variance of each PC in each epoch (% of the total variance of the epoch)
- `pcm::Vector{Float64}`: channel means of the model
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
    isempty(ch) && throw(ArgumentError("No channels selected."))

    return _pca_decompose(@view(obj.data[ch, :, :]); n = n)
end

"""
    _pca_reconstruct(s, pc, pca)

Reconstructs signal using PCA components for a 3-D signal array.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `pc::AbstractArray`: PC(1)..PC(n) × samples × epochs
- `pc_model::MultivariateStats.PCA{Float64}`: PC model

# Returns

- `Array{Float64, 3}`
"""
function _pca_reconstruct(
    s::AbstractArray;
    pc::AbstractArray,
    pc_model::MultivariateStats.PCA{Float64},
)::Array{Float64, 3}
    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate output
    s_new = similar(s, Float64)

    @inbounds for ep_idx = 1:ep_n
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
    isempty(ch) && throw(ArgumentError("No channels selected."))

    # create new dataset
    obj_tmp = deepcopy(obj)

    obj_tmp.data[ch, :, :] =
        _pca_reconstruct(@view(obj_tmp.data[ch, :, :]); pc = pc, pc_model = pc_model)
    push!(obj_tmp.history, "pca_reconstruct(obj; ch=$ch)")

    return obj_tmp
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
    obj_tmp = pca_reconstruct(obj; ch = ch, pc = pc, pc_model = pc_model)
    obj.data = obj_tmp.data
    obj.history = obj_tmp.history
    obj_tmp = nothing

    return nothing
end
