export aff_mni2tal
export aff_tal2mni
export mni2tal
export tal2mni

"""
    aff_mni2tal(pts)

Convert MNI coordinates to Talairach coordinates: redo the affine transform of Talairach coordinates.

# Arguments

- `pts::Vector{<:Number}`: MNI X, Y, Z coordinates

# Returns

- `t::Vector{Float64}`: Talairach X, Y, Z coordinates

# Source

https://www.brainmap.org/training/BrettTransform.html
"""
function aff_mni2tal(pts::Vector{<:Number})::Vector{Float64}

    @assert length(pts) == 3 "pts must contain 3 coordinates (x, y, z)."

    x_prime = 0.88 * pts[1] - 0.8
    y_prime = 0.97 * pts[2] - 3.32
    z_prime = 0.05 * pts[2] + 0.88 * pts[3] - 0.44

    t = ([x_prime, y_prime, z_prime])

    return t

end

"""
    aff_tal2mni(pts)

Convert Talairach coordinates to MNI coordinates: do the affine transform of MNI coordinates.

# Arguments

- `pts::Vector{<:Number}`: Talairach X, Y, Z coordinates

# Returns

- `t::Vector{Float64}`: MNI X, Y, Z coordinates

# Source

https://www.brainmap.org/training/BrettTransform.html
"""
function aff_tal2mni(pts::Vector{<:Number})::Vector{Float64}

    @assert length(pts) == 3 "pts must contain 3 coordinates (x, y, z)."

    x_prime = (pts[1] + 0.8) / 0.88
    y_prime = (pts[2] + 3.32) / 0.97
    z_prime = (pts[3] - 0.05 * ((pts[2] + 3.32) / 0.97) + 0.44) / 0.88

    m = ([x_prime, y_prime, z_prime])

    return m

end

"""
    mni2tal(pts)

Convert MNI coordinates to Talairach coordinates: a non-linear transform of MNI to Talairach coordinates.

# Arguments

- `pts::Vector{<:Number}`: MNI X, Y, Z coordinates

# Returns

- `t::Vector{Float64}`: Talairach X, Y, Z coordinates

# Source

https://www.brainmap.org/training/BrettTransform.html
"""
function mni2tal(pts::Vector{<:Number})::Vector{Float64}

    @assert length(pts) == 3 "pts must contain 3 coordinates (x, y, z)."

    if pts[3] >= 0
        x_prime = 0.9900 * pts[1]
        y_prime = 0.9688 * pts[2] + 0.0460 * pts[3]
        z_prime = -0.0485 * pts[2] + 0.9189 * pts[3]
    else
        x_prime = 0.9900 * pts[1]
        y_prime = 0.9688 * pts[2] + 0.0420 * pts[3]
        z_prime = -0.0485 * pts[2] + 0.8390 * pts[3]
    end

    t = ([x_prime, y_prime, z_prime])

    return t

end

"""
Convert Talairach coordinates to MNI coordinates: a non-linear transform of MNI to Talairach coordinates.

# Arguments

- `pts::Vector{<:Number}`: Talairach X, Y, Z coordinates

# Returns

- `m::Vector{Float64}`: MNI X, Y, Z coordinates

# Source

https://www.brainmap.org/training/BrettTransform.html
"""
function tal2mni(pts::Vector{<:Number})::Vector{Float64}

    @assert length(pts) == 3 "pts must contain 3 coordinates (x, y, z)."

    if pts[3] >= 0
        x_prime = pts[1] / 0.9900
        y_prime = 4.4819869616310097*10^-8 * (22972500 * pts[2] - 1150000 * pts[3])
        z_prime = 4.4819869616310097*10^-8 * (24220000 * pts[3] + 1212500 * pts[2])
    else
        x_prime = 0.9900 * pts[1]
        y_prime = 2.454408743978415*10^-7 * (4195000 * pts[2] - 210000 * pts[2])
        z_prime = 2.454408743978415*10^-7 * (4844000 * pts[2] + 242500 * pts[2])
    end

    m = ([x_prime, y_prime, z_prime])

    return m

end
