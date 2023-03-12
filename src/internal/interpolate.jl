function _interpolate(s::AbstractVector, loc_x::Vector{Float64}, loc_y::Vector{Float64}, interpolation_factor::Int64=100, imethod::Symbol=:sh, nmethod::Symbol=:minmax)
    # `imethod::Symbol=:sh`: interpolation method Shepard (`:sh`), Multiquadratic (`:mq`), InverseMultiquadratic (`:imq`), ThinPlate (`:tp`), NearestNeighbour (`:nn`), Gaussian (`:ga`)
    _check_var(imethod, [:sh, :mq, :imq, :tp, :nn, :ga], "imethod")
    x_lim_int = (-1.4, 1.4)
    y_lim_int = (-1.4, 1.4)
    interpolated_x = linspace(x_lim_int[1], x_lim_int[2], interpolation_factor)
    interpolated_y = linspace(y_lim_int[1], y_lim_int[2], interpolation_factor)
    interpolated_x = round.(interpolated_x, digits=2)
    interpolated_y = round.(interpolated_y, digits=2)
    interpolation_m = Matrix{Tuple{Float64, Float64}}(undef, interpolation_factor, interpolation_factor)
    @inbounds @simd for idx1 in 1:interpolation_factor
        for idx2 in 1:interpolation_factor
            interpolation_m[idx1, idx2] = (interpolated_x[idx1], interpolated_y[idx2])
        end
    end
    s_interpolated = zeros(interpolation_factor, interpolation_factor)
    electrode_locations = [loc_x loc_y]'
    imethod === :sh && (itp = ScatteredInterpolation.interpolate(Shepard(), electrode_locations, s))
    imethod === :mq && (itp = ScatteredInterpolation.interpolate(Multiquadratic(), electrode_locations, s))
    imethod === :imq && (itp = ScatteredInterpolation.interpolate(InverseMultiquadratic(), electrode_locations, s))
    imethod === :tp && (itp = ScatteredInterpolation.interpolate(ThinPlate(), electrode_locations, s))
    imethod === :nn && (itp = ScatteredInterpolation.interpolate(NearestNeighbor(), electrode_locations, s))
    imethod === :ga && (itp = ScatteredInterpolation.interpolate(Gaussian(), electrode_locations, s))
    @inbounds @simd for idx1 in 1:interpolation_factor
        for idx2 in 1:interpolation_factor
            s_interpolated[idx1, idx2] = ScatteredInterpolation.evaluate(itp, [interpolation_m[idx1, idx2][1]; interpolation_m[idx1, idx2][2]])[1]
        end
    end
    s_interpolated = s_interpolated'
    
    return normalize(s_interpolated, method=nmethod), interpolated_x, interpolated_y
end
