export aff_mni2tal
export aff_tal2mni
export mni2tal
export tal2mni

"""
    aff_mni2tal(pts)

Convert MNI coordinates to Talairach coordinates using an affine transform.

Applies the Brett affine transformation:
- `x′ = 0.88x − 0.8`
- `y′ = 0.97y − 3.32`
- `z′ = 0.05y + 0.88z − 0.44`

# Arguments

- `pts::Vector{<:Number}`: MNI `[X, Y, Z]` coordinates

# Returns

- `Vector{Float64}`: Talairach `[X, Y, Z]` coordinates

# Throws

- `ArgumentError`: if `pts` does not contain exactly 3 elements

# References
Brett M. https://www.brainmap.org/training/BrettTransform.html

# See also

[`aff_tal2mni`](@ref), [`mni2tal`](@ref)
"""
function aff_mni2tal(pts::Vector{<:Number})::Vector{Float64}

    !(length(pts) == 3) && throw(ArgumentError("pts must contain exactly 3 coordinates (x, y, z)."))

    x′ = 0.88  * pts[1] - 0.8
    y′ = 0.97  * pts[2] - 3.32
    z′ = 0.05  * pts[2] + 0.88 * pts[3] - 0.44

    return [x′, y′, z′]

end

"""
    aff_tal2mni(pts)

Convert Talairach coordinates to MNI coordinates using the inverse affine transform.

Applies the inverse of the Brett affine transformation:
- `x = (x′ + 0.8)  / 0.88`
- `y = (y′ + 3.32) / 0.97`
- `z = (z′ − 0.05y + 0.44) / 0.88`

# Arguments

- `pts::Vector{<:Number}`: Talairach `[X, Y, Z]` coordinates

# Returns

- `Vector{Float64}`: MNI `[X, Y, Z]` coordinates

# Throws

- `ArgumentError`: if `pts` does not contain exactly 3 elements

# References

Brett M. https://www.brainmap.org/training/BrettTransform.html

# See also

[`aff_mni2tal`](@ref), [`tal2mni`](@ref)
"""
function aff_tal2mni(pts::Vector{<:Number})::Vector{Float64}

    !(length(pts) == 3) && throw(ArgumentError("pts must contain exactly 3 coordinates (x, y, z)."))
    x′ = (pts[1] + 0.8)  / 0.88
    y′ = (pts[2] + 3.32) / 0.97
    # yeuse y′ to avoid repeating the y inversion inline
    z′ = (pts[3] - 0.05 * y′ + 0.44) / 0.88
    return [x′, y′, z′]

end

"""
    mni2tal(pts)

Convert MNI coordinates to Talairach coordinates using the non-linear Brett transform.

The transform is piecewise-linear in z:
- z ≥ 0 (above the AC–PC plane): `x′ = 0.99x`, `y′ = 0.9688y + 0.046z`, `z′ = −0.0485y + 0.9189z`
- z < 0 (below the AC–PC plane): `x′ = 0.99x`, `y′ = 0.9688y + 0.042z`, `z′ = −0.0485y + 0.839z`

# Arguments

- `pts::Vector{<:Number}`: MNI `[X, Y, Z]` coordinates

# Returns

- `t::Vector{Float64}`: Talairach `[X, Y, Z]` coordinates

# Throws

- `ArgumentError`: if `pts` does not contain exactly 3 elements

# References

Brett M. https://www.brainmap.org/training/BrettTransform.html

# See also

[`tal2mni`](@ref), [`aff_mni2tal`](@ref)
"""
function mni2tal(pts::Vector{<:Number})::Vector{Float64}

    !(length(pts) == 3) && throw(ArgumentError("pts must contain exactly 3 coordinates (x, y, z)."))

    # x scaling is identical in both branches
    x′ = 0.99 * pts[1]
    if pts[3] >= 0
        # above (or on) the AC–PC plane
        y′ =  0.9688 * pts[2] + 0.046  * pts[3]
        z′ = -0.0485 * pts[2] + 0.9189 * pts[3]
    else
        # below the AC–PC plane
        y′ =  0.9688 * pts[2] + 0.042 * pts[3]
        z′ = -0.0485 * pts[2] + 0.839 * pts[3]
    end
    return [x′, y′, z′]

end

"""
    tal2mni(pts)

Convert Talairach coordinates to MNI coordinates using the inverse non-linear Brett transform.

This is the analytical inverse of [`mni2tal`](@ref). The transform is piecewise-linear in z, with the branching threshold applied to the Talairach z coordinate:

z ≥ 0 — inverse of the upper matrix `[[0.9688, 0.046], [−0.0485, 0.9189]]` (det ≈ 0.8925):
- `x = x′ / 0.99`
- `y = (0.9189 y′ − 0.046  z′) / det`
- `z = (0.0485 y′ + 0.9688 z′) / det`

z < 0 — inverse of the lower matrix `[[0.9688, 0.042], [−0.0485, 0.839]]` (det ≈ 0.8151):
- `x = x′ / 0.99`
- `y = (0.839  y′ − 0.042  z′) / det`
- `z = (0.0485 y′ + 0.9688 z′) / det`

# Arguments

- `pts::Vector{<:Number}`: Talairach `[X, Y, Z]` coordinates

# Returns

- `m::Vector{Float64}`: MNI `[X, Y, Z]` coordinates

# Throws

- `ArgumentError`: if `pts` does not contain exactly 3 elements

# References

Brett M. https://www.brainmap.org/training/BrettTransform.html

# See also

[`mni2tal`](@ref), [`aff_tal2mni`](@ref)
"""
function tal2mni(pts::Vector{<:Number})::Vector{Float64}

    !(length(pts) == 3) && throw(ArgumentError("pts must contain exactly 3 coordinates (x, y, z)."))

    # x scaling is identical in both branches (inverse of 0.99)
    x′ = pts[1] / 0.99
    if pts[3] >= 0
        # inverse of upper 2×2 block; det = 0.9688×0.9189 − 0.046×(−0.0485)
        det = 0.9688 * 0.9189 + 0.046 * 0.0485
        y′  = ( 0.9189 * pts[2] - 0.046  * pts[3]) / det
        z′  = ( 0.0485 * pts[2] + 0.9688 * pts[3]) / det
    else
        # inverse of lower 2×2 block; det = 0.9688×0.839 − 0.042×(−0.0485)
        det = 0.9688 * 0.839 + 0.042 * 0.0485
        y′  = ( 0.839  * pts[2] - 0.042  * pts[3]) / det
        z′  = ( 0.0485 * pts[2] + 0.9688 * pts[3]) / det
    end
    return [x′, y′, z′]

end
