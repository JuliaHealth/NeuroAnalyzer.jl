export prank
export dranks

"""
    prank(x)

Calculate percentile ranks.

# Arguments

- `x::AbstractVector`: the vector to analyze

# Returns

- `p::Vector{Float64}`: percentile ranks
"""
function prank(x::AbstractVector)::Vector{Float64}

    xorder = sortperm(x)
    x = sort(x)
    p = zeros(length(x))

    for idx in eachindex(x)
        percentile = length(x[x .< x[idx]]) / length(x) * 100
        p[idx] = percentile / (100 * (length(x) + 1))
    end

    return p[xorder]

end

"""
    dranks(x, nbins)

Calculate ranks scaled in 0..nbins.

# Arguments

- `x::AbstractArray`: some continuous variable, such as reaction time (the time it takes to indicate the response)
- `nbins::Int64`: number of bins, default is Sturges' formula (`nbins = 1 + log2(length(x)))`)

# Returns

- `sr::Array{Int64}`
"""
function dranks(x::AbstractArray, nbins::Int64=round(Int64, 1 + log2(length(x))))::Array{Int64}

    # scale ranks in 0..1
    r = tiedrank(x) ./ length(x)
    
    # scale ranks in 0..nbins
    sr = ceil.(Int64, r .* nbins)

    return sr

end
