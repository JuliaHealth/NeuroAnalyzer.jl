function _interpolate2d(s::AbstractVector, loc_x::Vector{Float64}, loc_y::Vector{Float64}, ifactor::Int64=100, imethod::Symbol=:sh, nmethod::Symbol=:minmax)::Tuple{Matrix{Float64}, Vector{Float64}, Vector{Float64}}
    # `imethod::Symbol=:sh`: interpolation method Shepard (`:sh`), Multiquadratic (`:mq`), InverseMultiquadratic (`:imq`), ThinPlate (`:tp`), NearestNeighbour (`:nn`), Gaussian (`:ga`)

    _check_var(imethod, [:sh, :mq, :imq, :tp, :nn, :ga], "imethod")

    x_lim_int = (-1.2, 1.2)
    y_lim_int = (-1.2, 1.2)

    interpolated_x = linspace(x_lim_int[1], x_lim_int[2], ifactor)
    interpolated_y = linspace(y_lim_int[1], y_lim_int[2], ifactor)
    interpolated_x = round.(interpolated_x, digits=2)
    interpolated_y = round.(interpolated_y, digits=2)
    interpolation_m = Matrix{Tuple{Float64, Float64}}(undef, ifactor, ifactor)

    @inbounds for idx1 in 1:ifactor
        for idx2 in 1:ifactor
            interpolation_m[idx1, idx2] = (interpolated_x[idx1], interpolated_y[idx2])
        end
    end

    s_interpolated = zeros(ifactor, ifactor)

    electrode_locations = [loc_x loc_y]'

    imethod === :sh && (itp = ScatteredInterpolation.interpolate(Shepard(), electrode_locations, s))
    imethod === :mq && (itp = ScatteredInterpolation.interpolate(Multiquadratic(), electrode_locations, s))
    imethod === :imq && (itp = ScatteredInterpolation.interpolate(InverseMultiquadratic(), electrode_locations, s))
    imethod === :tp && (itp = ScatteredInterpolation.interpolate(ThinPlate(), electrode_locations, s))
    imethod === :nn && (itp = ScatteredInterpolation.interpolate(NearestNeighbor(), electrode_locations, s))
    imethod === :ga && (itp = ScatteredInterpolation.interpolate(Gaussian(), electrode_locations, s))

    @inbounds for idx1 in 1:ifactor
        for idx2 in 1:ifactor
            s_interpolated[idx1, idx2] = ScatteredInterpolation.evaluate(itp, [interpolation_m[idx1, idx2][1]; interpolation_m[idx1, idx2][2]])[1]
        end
    end

    s_interpolated = rotl90(s_interpolated)

    return normalize(s_interpolated, method=nmethod), interpolated_x, interpolated_y

end
