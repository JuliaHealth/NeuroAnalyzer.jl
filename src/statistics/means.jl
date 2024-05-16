export meang
export meanh
export meanw
export meanc

"""
    meang(x)

Calculate geometric mean.

# Arguments

- `x::AbstractVector`

# Returns

- `m::Float64`
"""
function meang(x::AbstractVector)

    return exp(mean(log.(x[x .> 0])))

end

"""
    meanh(x)

Calculate harmonic mean.

# Arguments

- `x::AbstractVector`

# Returns

- `m::Float64`
"""
function meanh(x::AbstractVector)

    return length(x) / sum(1 ./ x)

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
function meanw(x::AbstractVector, w::AbstractVector)

    @assert length(x) == length(w) "Weights and values vectors must have the same length."

    return length(x) / sum(1 ./ x)

end

"""
    meanc(x; rad)

Calculate circular mean.

# Arguments

- `x::AbstractVector`: angles
- `rad::Bool=false`: angles in radians (`rad=true`) or degrees (`rad=false`)

# Returns

- `m::Float64`
"""
function meanc(x::AbstractVector; rad::Bool=false)

    if rad
        return atan(sum(sin.(x)), sum(cos.(x)))
    else
        return rad2deg(atan(sum(sind.(x)), sum(cosd.(x))))
    end

end
