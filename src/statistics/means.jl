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
function meang(x::AbstractVector)::Float64

    m = exp(mean(log.(x[x .> 0])))

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

    @assert length(x) == length(w) "Weights and values vectors must have the same length."

    m = length(x) / sum(1 ./ x)

    return m

end

"""
    meanc(x; <keyword arguments>)

Calculate circular mean.

# Arguments

- `x::AbstractVector`: angles
- `rad::Bool=false`: angles in radians (`rad=true`) or degrees (`rad=false`)

# Returns

- `m::Float64`
"""
function meanc(x::AbstractVector; rad::Bool=false)::Float64

    if rad
        m = atan(sum(sin.(x)), sum(cos.(x)))
    else
        m = rad2deg(atan(sum(sind.(x)), sum(cosd.(x))))
    end

    return m

end
