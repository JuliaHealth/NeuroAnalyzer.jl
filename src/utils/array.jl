export l1
export l2
export perm_cmp
export tavg
export areduce

"""
    l1(a1, a2)

Compare two arrays using the L1 (Manhattan) distance.

Computes `∑|a1ᵢ - a2ᵢ|` element-wise.

Suitable for comparing spectrograms, feature maps, or any same-shaped numeric arrays.

# Arguments

- `a1::AbstractArray`: first array
- `a2::AbstractArray`: second array; must be the same size as `a1`.

# Returns

- `l1::Float64`: L1 distance between `a1` and `a2`
"""
function l1(a1::AbstractArray, a2::AbstractArray)::Float64

    @assert size(a1) == size(a2) "a1 and a2 must have the same size."

    return sum(abs.(a1 .- a2))

end

"""
    l2(a1, a2)

Compare two arrays using the L2 (Euclidean) distance.

Computes `√∑(a1ᵢ - a2ᵢ)²` element-wise via `Distances.euclidean`.

Suitable for comparing spectrograms, feature maps, or any same-shaped numeric arrays.

# Arguments

- `a1::AbstractArray`: first array
- `a2::AbstractArray`: second array; must be the same size as `a1`

# Returns

- `l2::Float64`: L2 distance between `a1` and `a2`.
"""
function l2(a1::AbstractArray, a2::AbstractArray)::Float64

    @assert size(a1) == size(a2) "a1 and a2 must have the same size."

    return euclidean(a1, a2)

end

"""
    perm_cmp(a1, a2; <keyword arguments>)

Compare two 3-dimensional arrays using a permutation-based statistic.

Randomly shuffles the combined pool of epochs `perm_n` times, splits each shuffle into two equal halves, and builds a null distribution of difference maps. The real difference map is then Z-scored against this null distribution
and thresholded at the Z-value corresponding to `p`.

# Arguments

- `a1::Array{<:Real, 3}`: first array (e.g. spectrogram), shape (freq, time, epochs)
- `a2::Array{<:Real, 3}`: second array; must match `a1` in size.
- `p::Float64=0.05`: two-tailed p-value threshold (must be in `(0, 1)`).
- `perm_n::Int64=1000`: number of permutations (must be > 0).

# Returns

Named tuple:

- `zmap::Matrix{Float64}`: Z-scored difference map `(a2 mean − a1 mean)` normalized by the permutation null distribution
- `bm::BitMatrix`: Boolean mask where `true` indicates a **statistically significant**
  position (`|z| ≥ zval`)
"""
function perm_cmp(
    a1::Array{<:Real, 3},
    a2::Array{<:Real, 3};
    p::Float64 = 0.05,
    perm_n::Int64 = 1000
)::@NamedTuple{zmap::Matrix{Float64}, bm::BitMatrix}

    @assert size(a1) == size(a2) "Both arrays must have the same size"
    @assert perm_n > 0 "perm_n must be > 0."
    @assert p >= 0 "p must be ≥ 0."
    @assert p <= 1 "p must be ≤ 1."

    # real observed difference (a2 − a1), averaged across epochs
    spec_diff = dropdims(mean(a2, dims = 3) .- mean(a1, dims = 3); dims = 3)

    # z-value threshold corresponding to the two-tailed p-value
    zval = abs(norminvcdf(p))

    # pool all epochs from both conditions
    spec_all = cat(a1, a2, dims=3)
    ep_n = size(spec_all, 3)
    half = ep_n ÷ 2

    # build null distribution via random epoch label permutations
    perm_maps = zeros(size(a1, 1), size(a1, 2), perm_n)
    @inbounds for perm_idx in 1:perm_n
        rand_idx  = sample(1:ep_n, ep_n; replace=false)
        rand_spec = @view spec_all[:, :, rand_idx]
        # difference between the two random halves → one null sample
        perm_maps[:, :, perm_idx] = @views dropdims(
            mean(rand_spec[:, :, (half + 1):end]; dims=3) .-
            mean(rand_spec[:, :, 1:half];         dims=3),
            dims=3,
        )
    end

    # H0 distribution statistics (mean and SD across permutations)
    mean_h0 = dropdims(mean(perm_maps; dims=3); dims=3)
    std_h0  = dropdims(std(perm_maps;  dims=3); dims=3)

    # z-score the real difference map against the null distribution
    zmap = @. (spec_diff - mean_h0) / std_h0

    # build binary significance mask: true = |z| ≥ threshold (significant)
    bm = BitMatrix(abs.(zmap) .>= zval)

    return (zmap=zmap, bm=bm)

end

"""
    tavg(s)

Average a 3-dimensional signal array across the trial (third) dimension.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)

# Returns

- `s_new::AbstractArray`: mean across epochs, shape (channels, samples, 1)
"""
function tavg(s::AbstractArray)::AbstractArray

    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    return mean(s, dims = 3)

end

"""
    areduce(a, f; <keyword arguments>)

Reduce an array by resampling along its second axis at regular intervals of `n`.

Useful for downsampling a frequency axis (and its associated data) when the number of frequency bins is very large. Only the nearest existing frequency bin to each target frequency is used (via `vsearch`); no interpolation is performed.

# Arguments

- `a::AbstractArray`: 2- or 3-dimensional data array whose second axis corresponds to `f`, shape (channels, frequencies) or (channels, frequencies, epochs)
- `f::AbstractVector`: frequency (or index) vector; `length(f)` must equal `size(a, 2)`
- `n::Float64=0.5`: step size between retained values (in the same units as `f`); smaller values retain more points; larger values reduce more aggressively

# Returns

- `a_new::Array{eltype(a), ndims(a)}`: deduced data array
- `f_new::Vector{eltype(f)}`: reduced frequency vector

# Throws

- `ArgumentError`: if `ndims(a) ∉ {2, 3}` or `size(a, 2) ≠ length(f)`
"""
function areduce(
    a::AbstractArray, f::AbstractVector; n::Float64 = 0.5
)::Tuple{AbstractArray, AbstractVector}

    @assert ndims(a) <= 3 "areduce() only works for 2- and 3-dimensional arrays."
    @assert size(a, 2) == length(f) "size(a, 2) ($(size(a, 2))) must equal length(f) ($(length(f)))."

    # build the reduced frequency grid from the rounded min/max frequencies
    f1 = round(f[vsearch(round(f[1]), f)])
    f2 = round(f[vsearch(round(f[end]), f)])
    f_new = collect(f1:n:f2)

    if ndims(a) == 2
        # allocate with matching element type to avoid silent precision loss
        a_new = zeros(eltype(a), size(a, 1), length(f_new))
        @inbounds for ch_idx in axes(a, 1)
            for (idx, freq) in enumerate(f_new)
                a_new[ch_idx, idx] = a[ch_idx, vsearch(freq, f)]
            end
        end
    else
        # allocate with matching element type to avoid silent precision loss
        a_new = zeros(eltype(a), size(a, 1), length(f_new), size(a, 3))
        @inbounds for ep_idx in axes(a, 3)
            for ch_idx in axes(a, 1)
                for (idx, freq) in enumerate(f_new)
                    a_new[ch_idx, idx, ep_idx] = a[ch_idx, vsearch(freq, f), ep_idx]
                end
            end
        end
    end

    return a_new, f_new

end
