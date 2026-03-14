export prank
export dranks

"""
    prank(x)

Calculate percentile ranks for each element of `x`.

The percentile rank of element `xᵢ` is the proportion of values in `x` that are strictly less than `xᵢ`, expressed as a value in `[0, 1)`: `PR(xᵢ) = count(x .< xᵢ) / length(x)`

Ties receive the same rank. Results are returned in the original element order.

# Arguments

- `x::AbstractVector`: input vector; must not be empty

# Returns

- `Vector{Float64}`: percentile ranks ∈ `[0, 1)`, in the same order as `x`

# Throws

- `ArgumentError`: if `x` is empty

# See also
[`dranks`](@ref)
"""
function prank(x::AbstractVector)::Vector{Float64}

    !(length(x) > 0) && throw(ArgumentError("x must not be empty."))
    n = length(x)

    return [count(<(xi), x) / n for xi in x]

end

"""
    dranks(x, nbins)

Scale tied ranks into discrete bins `1 … nbins`.

Tied ranks are computed with `StatsBase.tiedrank`, normalised to `(0, 1]`, then mapped to integers in `[1, nbins]` via `ceil`.

# Arguments

- `x::AbstractArray`: input array of continuous values; must not be empty
- `nbins::Int64`: number of bins; defaults to Sturges' formula `ceil(Int64, 1 + log2(length(x)))`; must be ≥ 1

# Returns

- `Array{Int64}`: rank-bin indices ∈ `[1, nbins]`, same shape as `x`

# Throws

- `ArgumentError`: if `x` is empty or `nbins < 1`

# See also

[`prank`](@ref)
"""
function dranks(
    x::AbstractArray,
    nbins::Int64=ceil(Int64, 1 + log2(length(x))),
)::Array{Int64}

    !(length(x) > 0) && throw(ArgumentError("x must not be empty."))
    !(nbins >= 1) && throw(ArgumentError("nbins must be ≥ 1."))

    # normalise tied ranks to (0, 1], then bin into 1..nbins
    r  = tiedrank(x) ./ length(x)

    return ceil.(Int64, r .* nbins)

end
