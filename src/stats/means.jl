export meanp
export meanc
export meang
export meanh
export meanw
export meancirc
export meant

"""
    meanp(p, n)

Calculate mean of the proportion.

# Arguments

- `p::Float64`: proportion
- `n::Int64`: number of observations

- `m::Float64`
"""
function meanp(p::Float64, n::Int64)::Float64

    @assert n > 0 "n must be > 0."
    _in(p, (0.0, 1.0), "p")

    m = n * p

    return m

end

"""
    meanc(g, x)

Calculate mean of categorical data.

# Arguments

- `g::Vector{Int64}`: group (e.g. [0, 1, 2, 3])
- `x::Vector{Int64}`: amount of subject per group (e.g. [2, 8, 27, 45])

# Returns

- `m::Float64`
"""
function meanc(g::Vector{Int64}, x::Vector{Int64})::Float64

    @assert length(g) == length(x) "Length of g and length of x must be equal."
    @assert length(g) > 0 "Length of g must be > 0."
    @assert length(x) > 0 "Length of x must be > 0."

    m = sum(g .* x) / sum(x)

    return m

end

"""
    meang(x)

Calculate geometric mean.

# Arguments

- `x::AbstractVector`

# Returns

- `m::Float64`
"""
function meang(x::AbstractVector)::Float64

    # @assert length(x[x .> 0]) != 0 "All values of x must be > 0."
    @assert length(x) > 0 "Length of x must be > 0."

    # m = exp(mean(log.(x)))
    m = prod(x)^(1/length(x))

    return m

end

"""
    meanh(x)

Calculate harmonic mean.

# Arguments

- `x::AbstractVector`

# Returns

- `m::Float64`
"""
function meanh(x::AbstractVector)::Float64

    @assert !(0 in x) "x must not contain value(s) of 0."

    m = length(x) / sum(1 ./ x)

    return m

end

"""
    meanw(x, w)

Calculate weighted mean.

# Arguments

- `x::AbstractVector`
- `w::AbstractVector`: weights

# Returns

- `m::Float64`
"""
function meanw(x::AbstractVector, w::AbstractVector)::Float64

    @assert !(0 in x) "x must not contain value(s) of 0."
    @assert length(x) == length(w) "Weights and values vectors must have the same length."
    
    m = sum(x .* w) / sum(x)

    return m

end

"""
    meancirc(x; <keyword arguments>)

Calculate circular mean.

# Arguments

- `x::AbstractVector`: angles
- `rad::Bool=false`: angles are provided in radians (`rad=true`) or degrees (`rad=false`)

# Returns

- `m::Float64`
"""
function meancirc(x::AbstractVector; rad::Bool=false)::Float64

    if rad
        m = atan(sum(sin.(x)), sum(cos.(x)))
    else
        m = rad2deg(atan(sum(sind.(x)), sum(cosd.(x))))
    end

    return m

end

"""
    meant(x; <keyword arguments>)

Calculate trimmed mean.

# Arguments

- `x::AbstractVector`
- `n::Float64=0.1`: percentage of extreme values to trim from both ends

# Returns

- `m::Float64`
"""
function meant(x::AbstractVector; n::Float64=0.1)::Float64

    @assert n > 0 "n must be > 0."
    _in(n, (0.0, 1.0), "n")

    xs = sort(x)
    xn = round(Int64, length(xs) * n)
    @assert xn + 1 < length(x) "n does not match x length."
    @assert length(x) - xn > 0 "n does not match x length."
    m = mean(x[xn + 1:end - xn])

    return m

end
