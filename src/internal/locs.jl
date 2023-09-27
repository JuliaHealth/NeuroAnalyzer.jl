function _initialize_locs()
    return DataFrame(:channel=>Int64,
                     :labels=>String[],
                     :loc_theta=>Float64[],
                     :loc_radius=>Float64[],
                     :loc_x=>Float64[],
                     :loc_y=>Float64[],
                     :loc_z=>Float64[],
                     :loc_radius_sph=>Float64[],
                     :loc_theta_sph=>Float64[],
                     :loc_phi_sph=>Float64[])
end

function _initialize_locs!(obj::NeuroAnalyzer.NEURO)
    locs_ch = signal_channels(obj)
    obj.locs = DataFrame(:channel=>signal_channels(obj), :labels=>labels(obj)[locs_ch], :loc_theta=>zeros(length(locs_ch)), :loc_radius=>zeros(length(locs_ch)), :loc_x=>zeros(length(locs_ch)), :loc_y=>zeros(length(locs_ch)), :loc_z=>zeros(length(locs_ch)), :loc_radius_sph=>zeros(length(locs_ch)), :loc_theta_sph=>zeros(length(locs_ch)), :loc_phi_sph=>zeros(length(locs_ch)))
    return nothing
end

function _initialize_locs(obj::NeuroAnalyzer.NEURO)
    locs_ch = signal_channels(obj)
    return DataFrame(:channel=>signal_channels(obj), :labels=>labels(obj)[locs_ch], :loc_theta=>zeros(length(locs_ch)), :loc_radius=>zeros(length(locs_ch)), :loc_x=>zeros(length(locs_ch)), :loc_y=>zeros(length(locs_ch)), :loc_z=>zeros(length(locs_ch)), :loc_radius_sph=>zeros(length(locs_ch)), :loc_theta_sph=>zeros(length(locs_ch)), :loc_phi_sph=>zeros(length(locs_ch)))
end

function _has_locs(obj::NeuroAnalyzer.NEURO)
    return nrow(obj.locs) > 0 ? true : false
end

function _locnorm(x::Union{AbstractVector, Real}, y::Union{AbstractVector, Real})
    xy = normalize_minmax(hcat(x, y))
    x = xy[:, 1]
    y = xy[:, 2]
    return x, y
end

function _locnorm(x::Union{AbstractVector, Real}, y::Union{AbstractVector, Real}, z::Union{AbstractVector, Real})
    xyz = normalize_minmax(hcat(x, y, z))
    x = xyz[:, 1]
    y = xyz[:, 2]
    z = xyz[:, 3]
    return x, y, z
end

function _round_locs(locs::DataFrame)
    locs[!, :loc_x] = round.(locs[!, :loc_x], digits=3)
    locs[!, :loc_y] = round.(locs[!, :loc_y], digits=3)
    locs[!, :loc_z] = round.(locs[!, :loc_z], digits=3)
    locs[!, :loc_radius] = round.(locs[!, :loc_radius], digits=3)
    locs[!, :loc_theta] = round.(locs[!, :loc_theta], digits=3)
    locs[!, :loc_radius_sph] = round.(locs[!, :loc_radius_sph], digits=3)
    locs[!, :loc_theta_sph] = round.(locs[!, :loc_theta_sph], digits=3)
    locs[!, :loc_phi_sph] = round.(locs[!, :loc_phi_sph], digits=3)
    return locs
end

function _angle_quadrant(a::Real)
    if a >= 0
        a = mod(a, 360)
        a <= 90 && (q = 1)
        (a > 90 && a <= 180) && (q = 2)
        (a > 180 && a <= 270) && (q = 3)
        (a > 270 && a < 360) && (q = 4)
    else
        a = mod(a, -360)
        a >= -90 && (q = 4)
        (a < -90 && a >= -180) && (q = 3)
        (a < -180 && a >= -270) && (q = 2)
        (a < -270 && a > -360) && (q = 1)
    end
    return q    
end 
