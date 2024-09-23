export l1
export l2
export perm_cmp
export tavg
export delmean
export areduce

"""
    l1(a1, a2)

Compare two arrays (e.g. two spectrograms), using L1 (Manhattan) distance.

# Arguments

- `a1::AbstractArray`: first array
- `a2::AbstractArray`: second array

# Returns

- `l1::Float64`
"""
function l1(a1::AbstractArray, a2::AbstractArray)::Float64

    @assert size(a1) == size(a2) "a1 and a2 mast have the same size."

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
function l2(a1::AbstractArray, a2::AbstractArray)::Float64

    @assert size(a1) == size(a2) "a1 and a2 mast have the same size."

    # return sqrt(sum((a1 .- a2).^2))
    return euclidean(a1, a2)

end

"""
    perm_cmp(a1, a2; <keyword arguments>)

Compare two 3-dimensional arrays (e.g. two spectrograms), using permutation based statistic.

# Arguments

- `a1::Array{<:Real, 3}`: first array
- `a2::Array{<:Real, 3}`: second array
- `p::Float64=0.05`: p-value
- `perm_n::Int64=1000`: number of permutations

# Returns

Named tuple containing:
- `zmap::Array{Float64, 3}`: array of Z-values
- `bm::Array{Float64, 3}`: binarized mask of statistically significant positions
"""
function perm_cmp(a1::Array{<:Real, 3}, a2::Array{<:Real, 3}; p::Float64=0.05, perm_n::Int64=1000)::NamedTuple{(:zmap, :bm), Tuple{Array{Float64, 3}, Array{Float64, 3}}}

    @assert size(a1) == size(a2) "Both arrays must have the same size"
    @assert perm_n > 0 "perm_n must be > 0."
    @assert p >= 0 "p must be ≥ 0."
    @assert p <= 1 "p must be ≤ 1."

    spec_diff = dropdims(mean(a2, dims=3) .- mean(a1, dims=3), dims=3)
    zval = abs(norminvcdf(p))
    perm_n = 1000
    spec_all = cat(a1, a2, dims=3)
    perm_maps = zeros(size(a1, 1), size(a1, 2), perm_n)
    ep_n = size(spec_all, 3)
    @inbounds for perm_idx in 1:perm_n
        rand_idx = sample(1:ep_n, ep_n, replace=false)
        rand_spec = @view spec_all[:, :, rand_idx]
        perm_maps[:, :, perm_idx] = @views dropdims(mean(rand_spec[:, :, (ep_n ÷ 2 + 1):end], dims=3) .- mean(rand_spec[:, :, 1:(ep_n ÷ 2)], dims=3), dims=3)
    end
    mean_h0 = dropdims(mean(perm_maps, dims=3), dims=3)
    std_h0 = dropdims(std(perm_maps, dims=3), dims=3)

    # threshold real data
    zmap = @. (spec_diff - mean_h0) / std_h0

    # threshold at p-value
    bm = deepcopy(zmap)
    @. bm[abs(zmap) < zval] = 0
    @. bm[bm != 0] = 1
    bm = Bool.(bm)
    bm = map(x -> !x, bm)

    return (zmap=zmap, bm=bm)

end

"""
    tavg(s)

Average signal across trials.

# Arguments

- `s::AbstractArray`

# Returns

- `s_new::AbstractArray`
"""
function tavg(s::AbstractArray)::AbstractArray

    @assert ndims(s) == 3 "Signal must have 3 dimensions."

    return mean(s, dims=3)

end

"""
    delmean(s)

Demean signal.

# Arguments

- `s::AbstractArray`

# Returns

- `s_new::AbstractArray`
"""
function delmean(s::AbstractArray; dims::Union{Int64, Nothing}=nothing)::AbstractArray

    ms = 0
    if isnothing(dims)
        ms = mean(s)
    else
        @assert dims <= ndims(s) "dims must be ≤ $(ndims(s))"
        ms = mean(s, dims=dims)
    end

    return s .- ms

end

"""
    areduce(a, f; <keyword arguments>)

Reduce an array at indices of a vector being multiplications of a constant. Useful e.g. for simplifying values across frequencies, when the number of frequencies (and thus values) is high.

# Arguments

- `a::AbstractArray`: e.g. signal data
- `f::AbstractVector`: e.g. frequencies
- `n::Float64=0.5`: reduce at multiplications of this value

# Returns

- `a_new::Array{eltype(a), ndims(a)}`
- `f_new::Vector{eltype(f)}`
"""
function areduce(a::AbstractArray, f::AbstractVector; n::Float64=0.5)::Tuple{AbstractArray, AbstractVector}

    @assert ndims(a) <= 3 "areduce() only works for 2- and 3-dimensional arrays."
    @assert size(a, 2) == length(f) "Length of both vectors must be equal."

    f1_idx = vsearch(round(f[1]), f)
    f2_idx = vsearch(round(f[end]), f)
    f1 = round(f[f1_idx])
    f2 = round(f[f2_idx])
    f_new = collect(f1:n:f2)

    if ndims(a) == 2
        a_new = zeros(size(a, 1), length(f_new))
        @inbounds for ch_idx in axes(a, 1)
            for idx in eachindex(f_new)
                f_idx = vsearch(f_new[idx], f)
                a_new[ch_idx, idx] = a[ch_idx, f_idx]
            end
        end
    else
        a_new = zeros(size(a, 1), length(f_new), size(a, 3))
        @inbounds for ep_idx in axes(a, 3)
            for ch_idx in axes(a, 1)
                for idx in eachindex(f_new)
                    f_idx = vsearch(f_new[idx], f)
                    a_new[ch_idx, idx, ep_idx] = a[ch_idx, f_idx, ep_idx]
                end
            end
        end
    end

    return a_new, f_new

end