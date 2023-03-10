export l1
export l2
export perm_cmp

"""
    l1(a1, a2)

Compare two arrays (e.g. two spectrograms), using L1 (Manhattan) distance.

# Arguments

- `a1::AbstractArray`: first array
- `a2::AbstractArray`: second array

# Returns

- `l1::Float64`
"""
function l1(a1::AbstractArray, a2::AbstractArray)
    size(a1) == size(a2) || throw(ArgumentError("a1 and a2 mast have the same size."))
    return sum(abs.(a1 .- a2))
end

"""
    l2(a1, a2)

Compare two arrays (e.g. two spectrograms), using L2 (Euclidean) distance.

# Arguments

- `a1::AbstractArray`: first array
- `a2::AbstractArray`: second array

# Returns

- `l2::Float64`
"""
function l2(a1::AbstractArray, a2::AbstractArray)
    size(a1) == size(a2) || throw(ArgumentError("a1 and a2 mast have the same size."))
    # return sqrt(sum((a1 .- a2).^2))
    return euclidean(a1, a2)
end

"""
    perm_cmp(a1, a2; p, perm_n)

Compare two 3-dimensional arrays `a1` and `a2` (e.g. two spectrograms), using permutation based statistic.

# Arguments

- `a1::Array{<:Real, 3}`: first array
- `a2::Array{<:Real, 3}`: second array
- `p::Float64=0.05`: p-value
- `perm_n::Int64=1000`: number of permutations

# Returns

Named tuple containing:
- `zmap::Array{Float64, 3}`: array of Z-values
- `zmap_b::Array{Float64, 3}`: binarized mask of statistically significant positions
"""
function perm_cmp(a1::Array{<:Real, 3}, a2::Array{<:Real, 3}; p::Float64=0.05, perm_n::Int64=1000)

    size(a1) == size(a2) || throw(ArgumentError("Both arrays must have the same size"))
    perm_n <= 0 && throw(ArgumentError("perm_n must be > 0."))
    p < 0 && throw(ArgumentError("p must be ≥ 0."))
    p > 1 && throw(ArgumentError("p must be ≤ 1."))

    spec_diff = dropdims(mean(a2, dims=3) .- mean(a1, dims=3), dims=3)
    zval = abs(norminvcdf(p))
    perm_n = 1000
    spec_all = cat(a1, a2, dims=3)
    perm_maps = zeros(size(a1, 1), size(a1, 2), perm_n)
    ep_n = size(spec_all, 3)
    @inbounds @simd for perm_idx in 1:perm_n
        rand_idx = sample(1:ep_n, ep_n, replace=false)
        rand_spec = @view spec_all[:, :, rand_idx]
        perm_maps[:, :, perm_idx] = @views dropdims(mean(rand_spec[:, :, (ep_n ÷ 2 + 1):end], dims=3) .- mean(rand_spec[:, :, 1:(ep_n ÷ 2)], dims=3), dims=3)
    end
    mean_h0 = dropdims(mean(perm_maps, dims=3), dims=3)
    std_h0 = dropdims(std(perm_maps, dims=3), dims=3)

    # threshold real data
    zmap = @. (spec_diff - mean_h0) / std_h0

    # threshold at p-value
    zmap_b = deepcopy(zmap)
    @. zmap_b[abs(zmap) < zval] = 0
    @. zmap_b[zmap_b != 0] = 1
    zmap_b = Bool.(zmap_b)
    zmap_b = map(x -> !x, zmap_b)

    return (zmap=zmap, zmap_b=zmap_b)
end
