export prank
export dranks

"""
    prank(x)

Calculate percentile rank.

# Arguments

- `x::AbstractVector`: the vector to analyze

# Returns

- `prnk::Vector{Float64}`
"""
function prank(x::AbstractVector)
    xorder = sortperm(x)
    x = sort(x)
    prnk = zeros(length(x))
    for idx in eachindex(x)
        percentile = length(x[x .< x[idx]]) / length(x) * 100
        prnk[idx] = percentile / (100 * (length(x) + 1))
    end
    return prnk[xorder]
end

"""
    dranks(x, nbins)

Calculate ranks scaled in 0..nbins. Number of bins is calculated automatically using Sturges' formula.

# Arguments

- `x::AbstractArray`: some continuous variable such as reaction time (the time it takes to indicate the response)
- `nbins::Int64`: number of bins, default is Sturges' formula

# Returns

- `caf::Array{Float64}`
"""
function dranks(x::AbstractArray, nbins::Int64=round(Int64, 1 + log2(length(x))))
    # scale ranks in 0..1
    ranks = tiedrank(x) ./ length(x)
    # scale ranks in 0..nbins
    return ceil.(Int64, ranks .* nbins)
end
