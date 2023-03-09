export gfp
export gfp_norm
export diss

"""
    gfp(signal)

Calculate GFP (Global Field Power).

# Arguments

- `signal::AbstractVector`

# Returns

- `gfp::Float64`
"""
function gfp(signal::AbstractVector)
    return sum(signal.^2) / length(signal)
end

"""
    gfp_norm(signal)

Calculate signal normalized for GFP (Global Field Power).

# Arguments

- `signal::AbstractVector`

# Returns

- `gfp_norm::Float64`
"""
function gfp_norm(signal::AbstractVector)
    return signal ./ s_gfp(signal)
end

"""
    diss(signal1, signal2)

Calculate DISS (global dissimilarity) and spatial correlation between `signal1` and `signal2`.

# Arguments

- `signal1::AbstractVector`
- `signal2::AbstractVector`

# Returns

Named tuple containing:
- `glob_diss::Float64`: global dissimilarity
- `c::Float64`: spatial correlation
"""
function diss(signal1::AbstractVector, signal2::AbstractVector)
    
    length(signal1) == length(signal2) || throw(ArgumentError("Both signals must have the same length."))
    gfp_norm1 = gfp_norm(signal1)
    gfp_norm2 = gfp_norm(signal2)
    glob_diss = sqrt(sum((gfp_norm1 .- gfp_norm2).^2) / length(signal1))
    c = 0.5 * (2 - glob_diss^2)

    return (glob_diss=glob_diss, c=c)
end

