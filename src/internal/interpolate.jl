"""
    _xy2sph(x, y)

Map planar (topographic) coordinates onto the unit sphere, inverting the azimuthal equidistant projection: the distance from the origin `r` is the polar angle from the vertex (`r = 1` at the equator, `r = 2` at the bottom pole) and the planar angle is the azimuth.

# Arguments

- `x::AbstractVector{<:Real}`: x-coordinates
- `y::AbstractVector{<:Real}`: y-coordinates

# Returns

- `Matrix{Float64}`: unit vectors (rows), shape (points, 3)
"""
function _xy2sph(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})::Matrix{Float64}
    polar = @. hypot(x, y) * π / 2
    azimuth = @. atan(y, x)
    return [(@. sin(polar) * cos(azimuth)) (@. sin(polar) * sin(azimuth)) cos.(polar)]
end

"""
    _interpolator(imethod, loc_x, loc_y, s; <keyword arguments>)

Scattered-data interpolator of values `s` at points `(loc_x, loc_y)`.

# Arguments

- `imethod::Symbol`: interpolation method:
    - ScatteredInterpolation.jl: `:sh` (Shepard), `:mq` (Multiquadratic), `:imq` (Inverse Multiquadratic), `:tp` (Thin Plate Spline), `:nn` (Nearest Neighbor), `:ga` (Gaussian)
    - NaturalNeighbours.jl: `:sib` (Sibson), `:sib1` (Sibson-1), `:lap` (Laplace), `:far` (Farin), `:hiy` (Hiyoshi), `:tri` (piecewise linear on Delaunay triangles)
    - `:sph`: spherical spline (Perrin et al. 1989; `m = 4`, 50 Legendre terms, `λ = 1e-8`, so the map passes through the data values), fitted on the unit sphere, see `_xy2sph()`
- `loc_x::AbstractVector{<:Real}`: x-coordinates of data points
- `loc_y::AbstractVector{<:Real}`: y-coordinates of data points
- `s::AbstractVector{<:Real}`: values at the data points
- `loc_xyz::Union{Nothing, AbstractMatrix{<:Real}}=nothing`: 3D positions (rows) of the data points, shape (points, 3); if given, `:sph` fits the spline at these positions (projected onto the unit sphere) instead of at `(loc_x, loc_y)` mapped onto the sphere by `_xy2sph()`; ignored by other methods

# Returns

- `Function`: `itp(x, y)` evaluates the interpolator at points `(x[i], y[i])` and returns `Vector{Float64}`

# Notes

Natural neighbour methods need at least 3 non-collinear data points; outside their convex hull values are extrapolated.
"""
function _interpolator(
    imethod::Symbol,
    loc_x::AbstractVector{<:Real},
    loc_y::AbstractVector{<:Real},
    s::AbstractVector{<:Real};
    loc_xyz::Union{Nothing, AbstractMatrix{<:Real}} = nothing,
)::Function
    _check_var(imethod, _IMETHODS, "imethod")

    if imethod === :sph
        # spherical spline with a constant term (as in sphinterpolate_channel(), but with
        # weaker regularization):
        # [G + λI 1; 1ᵀ 0] [c; c0] = [s; 0]
        if isnothing(loc_xyz)
            xyz = _xy2sph(loc_x, loc_y)
        else
            xyz = Matrix{Float64}(loc_xyz)
            r = sqrt.(sum(abs2, xyz; dims = 2))
            any(iszero, r) && throw(ArgumentError("loc_xyz must not contain points at the origin."))
            xyz ./= r
        end
        np = length(s)
        A = [(_sph_g(xyz, xyz; m = 4, n = 50) + 1e-8 * I(np)) ones(np);
             ones(1, np) 0.0]
        c = A \ [Vector{Float64}(s); 0.0]
        return (x, y) -> [_sph_g(_xy2sph(x, y), xyz; m = 4, n = 50) ones(length(x))] * c
    end

    if imethod in (:sib, :sib1, :lap, :far, :hiy, :tri)
        m = if imethod === :sib
            NaturalNeighbours.Sibson()
        elseif imethod === :sib1
            NaturalNeighbours.Sibson(1)
        elseif imethod === :lap
            NaturalNeighbours.Laplace()
        elseif imethod === :far
            NaturalNeighbours.Farin()
        elseif imethod === :hiy
            NaturalNeighbours.Hiyoshi(2)
        else
            NaturalNeighbours.Triangle()
        end
        # gradients (and Hessians) are needed by the smooth methods only
        itp_nn = NaturalNeighbours.interpolate(
            Vector{Float64}(loc_x),
            Vector{Float64}(loc_y),
            Vector{Float64}(s);
            derivatives = imethod in (:sib1, :far, :hiy),
        )
        # multithreading pays off only for many query points (e.g. a grid)
        return (x, y) -> itp_nn(
            Vector{Float64}(x), Vector{Float64}(y);
            method = m, parallel = length(x) > 1000,
        )
    end

    m = if imethod === :sh
        Shepard()
    elseif imethod === :mq
        Multiquadratic()
    elseif imethod === :imq
        InverseMultiquadratic()
    elseif imethod === :tp
        ThinPlate()
    elseif imethod === :nn
        NearestNeighbor()
    else
        Gaussian()
    end
    # data points as a (2 × n_points) matrix
    itp_si = ScatteredInterpolation.interpolate(m, Matrix{Float64}([loc_x loc_y]'), Vector{Float64}(s))
    # query points as a (2 × n_points) matrix
    return (x, y) -> vec(ScatteredInterpolation.evaluate(itp_si, Matrix{Float64}([x y]')))
end

"""
    _interpolate2d(s, loc_x, loc_y, ifactor, imethod, nmethod; <keyword arguments>)

Perform 2D scattered interpolation of signal data on a regular grid.

# Arguments

- `s::AbstractVector`: signal values to interpolate (length must match `loc_x` and `loc_y`)
- `loc_x::Vector{Float64}`: x-coordinates of data points
- `loc_y::Vector{Float64}`: y-coordinates of data points
- `ifactor::Int64=100`: interpolation factor determining grid density
- `imethod::Symbol=:tp`: interpolation method:
    - `:sh`: Shepard interpolation
    - `:mq`: Multiquadratic interpolation
    - `:imq`: Inverse Multiquadratic interpolation
    - `:tp`: Thin Plate Spline interpolation
    - `:nn`: Nearest Neighbor interpolation
    - `:ga`: Gaussian interpolation
    - `:sib`: Sibson (natural neighbour) interpolation
    - `:sib1`: Sibson-1 (smooth natural neighbour) interpolation
    - `:lap`: Laplace (non-Sibsonian natural neighbour) interpolation
    - `:far`: Farin (C1 natural neighbour) interpolation
    - `:hiy`: Hiyoshi (C2 natural neighbour) interpolation
    - `:tri`: piecewise linear interpolation on Delaunay triangles
    - `:sph`: spherical spline (Perrin et al. 1989) interpolation
- `nmethod::Symbol=:minmax`: normalization method for output
- `loc_xyz::Union{Nothing, AbstractMatrix{<:Real}}=nothing`: 3D positions of the data points for `:sph`, see `_interpolator()`

# Returns

- `Tuple{Matrix{Float64}, Vector{Float64}, Vector{Float64}}`: a tuple containing:
    - Interpolated signal matrix, shape (`ifactor`, `ifactor`)
    - x-coordinates of the interpolation grid
    - y-coordinates of the interpolation grid

# Notes

- The interpolation creates a regular grid covering the area containing all electrodes.
- The grid size is determined by `ifactor`.
- The output is normalized using the specified method.
"""
function _interpolate2d(
    s::AbstractVector,
    loc_x::Vector{Float64},
    loc_y::Vector{Float64},
    ifactor::Int64 = 100,
    imethod::Symbol = :tp,
    nmethod::Symbol = :minmax;
    loc_xyz::Union{Nothing, AbstractMatrix{<:Real}} = nothing,
)::Tuple{Matrix{Float64}, Vector{Float64}, Vector{Float64}}
    # validate
    ifactor > 0 || throw(ArgumentError("Interpolation factor (ifactor) must be positive"))
    length(s) == length(loc_x) == length(loc_y) ||
        throw(ArgumentError("Signal and location vectors must have equal lengths"))
    _check_var(imethod, _IMETHODS, "imethod")

    # calculate grid limits based on electrode positions
    max_abs_loc =
        max(ceil(maximum(abs, loc_x); digits = 1), ceil(maximum(abs, loc_y); digits = 1))
    # expand grid slightly beyond electrode positions
    extr = max_abs_loc > 1.2 ? 1.6 : 1.2

    # create interpolation grid
    interpolated_x = collect(range(-extr, extr; length = ifactor))
    interpolated_y = collect(range(-extr, extr; length = ifactor))

    # evaluate all grid points at once: s_interpolated[idx1, idx2] is the value at
    # (interpolated_x[idx1], interpolated_y[idx2]) (column-major order)
    itp = _interpolator(imethod, loc_x, loc_y, s; loc_xyz = loc_xyz)
    s_interpolated = reshape(
        itp(repeat(interpolated_x, ifactor), repeat(interpolated_y; inner = ifactor)),
        ifactor,
        ifactor,
    )

    # rotate matrix to match standard orientation (origin at top-left)
    s_interpolated = rotl90(s_interpolated)

    # normalize the interpolated data
    normalized_signal = NeuroAnalyzer._normalize(s_interpolated; method = nmethod)

    return normalized_signal, interpolated_x, interpolated_y
end
