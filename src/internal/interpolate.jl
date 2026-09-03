"""
    _interpolate2d(s, loc_x, loc_y, ifactor, imethod, nmethod)

Perform 2D scattered interpolation of signal data on a regular grid.

# Arguments

- `s::AbstractVector`: signal values to interpolate (length must match `loc_x` and `loc_y`)
- `loc_x::Vector{Float64}`: x-coordinates of data points
- `loc_y::Vector{Float64}`: y-coordinates of data points
- `ifactor::Int64=100`: interpolation factor determining grid density
- `imethod::Symbol=:sh`: interpolation method:
    - `:sh`: Shepard interpolation
    - `:mq`: Multiquadratic interpolation
    - `:imq`: Inverse Multiquadratic interpolation
    - `:tp`: Thin Plate Spline interpolation
    - `:nn`: Nearest Neighbor interpolation
    - `:ga`: Gaussian interpolation
- `nmethod::Symbol=:minmax`: normalization method for output

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
    imethod::Symbol = :sh,
    nmethod::Symbol = :minmax,
)::Tuple{Matrix{Float64}, Vector{Float64}, Vector{Float64}}
    # validate
    ifactor > 0 || throw(ArgumentError("Interpolation factor (ifactor) must be positive"))
    length(s) == length(loc_x) == length(loc_y) ||
        throw(ArgumentError("Signal and location vectors must have equal lengths"))
    _check_var(imethod, [:sh, :mq, :imq, :tp, :nn, :ga], "imethod")

    # calculate grid limits based on electrode positions
    max_abs_loc =
        max(ceil(maximum(abs, loc_x); digits = 1), ceil(maximum(abs, loc_y); digits = 1))
    # expand grid slightly beyond electrode positions
    extr = max_abs_loc > 1.2 ? 1.6 : 1.2

    # create interpolation grid
    x_lim_int = (-extr, extr)
    y_lim_int = (-extr, extr)
    interpolated_x = range(x_lim_int[1], x_lim_int[2]; length = ifactor) |> collect
    interpolated_y = range(y_lim_int[1], y_lim_int[2]; length = ifactor) |> collect

    # round for cleaner visualization (optional)
    interpolated_x = round.(interpolated_x; digits = 2)
    interpolated_y = round.(interpolated_y; digits = 2)

    # pre-allocate interpolation matrix
    interpolation_m = Matrix{NTuple{2, Float64}}(undef, ifactor, ifactor)
    s_interpolated = zeros(ifactor, ifactor)

    # prepare electrode locations as matrix (2 × n_points)
    electrode_locations = [loc_x loc_y]'

    # create appropriate interpolator based on method
    itp = if imethod === :sh
        ScatteredInterpolation.interpolate(Shepard(), electrode_locations, s)
    elseif imethod === :mq
        ScatteredInterpolation.interpolate(Multiquadratic(), electrode_locations, s)
    elseif imethod === :imq
        ScatteredInterpolation.interpolate(InverseMultiquadratic(), electrode_locations, s)
    elseif imethod === :tp
        ScatteredInterpolation.interpolate(ThinPlate(), electrode_locations, s)
    elseif imethod === :nn
        ScatteredInterpolation.interpolate(NearestNeighbor(), electrode_locations, s)
    elseif imethod === :ga
        ScatteredInterpolation.interpolate(Gaussian(), electrode_locations, s)
    end

    # perform interpolation on each grid point
    @inbounds for idx1 = 1:ifactor
        for idx2 = 1:ifactor
            x_val, y_val =
                interpolation_m[idx1, idx2] = (interpolated_x[idx1], interpolated_y[idx2])
            s_interpolated[idx1, idx2] =
                ScatteredInterpolation.evaluate(itp, [x_val, y_val])[1]
        end
    end

    # rotate matrix to match standard orientation (origin at top-left)
    s_interpolated = rotl90(s_interpolated)

    # normalize the interpolated data
    normalized_signal = NeuroAnalyzer.normalize(s_interpolated; method = nmethod)

    return normalized_signal, interpolated_x, interpolated_y
end
