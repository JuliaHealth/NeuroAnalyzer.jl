function _find_bylabel(locs::DataFrame, l::Union{String, Vector{String}, Vector{SubString{String}}})
    if typeof(l) == String
        if length(findall(occursin.(lowercase(l), lowercase.(locs[!, :labels])))) > 0
            return findall(occursin.(lowercase(l), lowercase.(locs[!, :labels])))[1]
        else
            return Int64[]
        end
    else
        l_idx = Vector{Int64}()
        for idx in l
            lowercase(idx) in lowercase.(locs[!, :labels]) && push!(l_idx, findall(occursin.(lowercase(idx), lowercase.(locs[!, :labels])))[1])
        end
        return l_idx
    end
end

function _initialize_locs()
    return DataFrame(:labels=>String[],
                     :loc_radius=>Float64[],
                     :loc_theta=>Float64[],
                     :loc_x=>Float64[],
                     :loc_y=>Float64[],
                     :loc_z=>Float64[],
                     :loc_radius_sph=>Float64[],
                     :loc_theta_sph=>Float64[],
                     :loc_phi_sph=>Float64[])
end

function _initialize_locs!(obj::NeuroAnalyzer.NEURO)
    locs_ch = signal_channels(obj)
    obj.locs = DataFrame(:labels=>labels(obj)[locs_ch], :loc_radius=>zeros(length(locs_ch)), :loc_theta=>zeros(length(locs_ch)), :loc_x=>zeros(length(locs_ch)), :loc_y=>zeros(length(locs_ch)), :loc_z=>zeros(length(locs_ch)), :loc_radius_sph=>zeros(length(locs_ch)), :loc_theta_sph=>zeros(length(locs_ch)), :loc_phi_sph=>zeros(length(locs_ch)))
    return nothing
end

function _initialize_locs(obj::NeuroAnalyzer.NEURO)
    locs_ch = signal_channels(obj)
    return DataFrame(:labels=>labels(obj)[locs_ch], :loc_radius=>zeros(length(locs_ch)), :loc_theta=>zeros(length(locs_ch)), :loc_x=>zeros(length(locs_ch)), :loc_y=>zeros(length(locs_ch)), :loc_z=>zeros(length(locs_ch)), :loc_radius_sph=>zeros(length(locs_ch)), :loc_theta_sph=>zeros(length(locs_ch)), :loc_phi_sph=>zeros(length(locs_ch)))
end

function _locs_round(locs::DataFrame)
    locs_new = deepcopy(locs)
    locs_new[!, :loc_radius] = round.(locs[!, :loc_radius], digits=2)
    locs_new[!, :loc_theta] = round.(locs[!, :loc_theta], digits=2)
    locs_new[!, :loc_x] = round.(locs[!, :loc_x], digits=2)
    locs_new[!, :loc_y] = round.(locs[!, :loc_y], digits=2)
    locs_new[!, :loc_z] = round.(locs[!, :loc_z], digits=2)
    locs_new[!, :loc_radius_sph] = round.(locs[!, :loc_radius_sph], digits=2)
    locs_new[!, :loc_theta_sph] = round.(locs[!, :loc_theta_sph], digits=2)
    locs_new[!, :loc_phi_sph] = round.(locs[!, :loc_phi_sph], digits=2)
    return locs_new
end

function _locs_round!(locs::DataFrame)
    locs[!, :loc_radius] = round.(locs[!, :loc_radius], digits=2)
    locs[!, :loc_theta] = round.(locs[!, :loc_theta], digits=2)
    locs[!, :loc_x] = round.(locs[!, :loc_x], digits=2)
    locs[!, :loc_y] = round.(locs[!, :loc_y], digits=2)
    locs[!, :loc_z] = round.(locs[!, :loc_z], digits=2)
    locs[!, :loc_radius_sph] = round.(locs[!, :loc_radius_sph], digits=2)
    locs[!, :loc_theta_sph] = round.(locs[!, :loc_theta_sph], digits=2)
    locs[!, :loc_phi_sph] = round.(locs[!, :loc_phi_sph], digits=2)
end

function _locs_remove_nans(locs::DataFrame)
    locs_new = deepcopy(locs)
    for idx in eachcol(locs_new)
        replace!(idx, NaN => 0)
    end
    return locs_new
end

function _locs_remove_nans!(locs::DataFrame)
    for idx in eachcol(locs)
        replace!(idx, NaN => 0)
    end
end

_locs_round(obj::NeuroAnalyzer.NEURO) = _locs_round(obj.locs)

function _locs_round!(obj::NeuroAnalyzer.NEURO)
    obj.locs = _locs_round(obj.locs)
end

_has_locs(obj::NeuroAnalyzer.NEURO) = nrow(obj.locs) > 0 ? true : false

function _locs_norm(x::Union{AbstractVector, Real}, y::Union{AbstractVector, Real})
    xy = normalize_minmax(hcat(x, y))
    x = xy[:, 1]
    y = xy[:, 2]
    return x, y
end

function _locs_norm(x::Union{AbstractVector, Real}, y::Union{AbstractVector, Real}, z::Union{AbstractVector, Real})
    xyz = normalize_minmax(hcat(x, y, z))
    x = xyz[:, 1]
    y = xyz[:, 2]
    z = xyz[:, 3]
    return x, y, z
end

function _locs_norm(locs::DataFrame)
    locs_new = deepcopy(locs)
    x, y, z = locs[!, :loc_x], locs[!, :loc_y], locs[!, :loc_z]
    xyz = normalize_minmax(hcat(x, y, z))
    x = xyz[:, 1]
    y = xyz[:, 2]
    z = xyz[:, 3]
    locs_new[!, :loc_x], locs_new[!, :loc_y], locs_new[!, :loc_z] = x, y, 
    return locs_new
end

function _locs_norm!(locs::DataFrame)
    x, y, z = locs[!, :loc_x], locs[!, :loc_y], locs[!, :loc_z]
    xyz = normalize_minmax(hcat(x, y, z))
    x = xyz[:, 1]
    y = xyz[:, 2]
    z = xyz[:, 3]
    locs[!, :loc_x], locs[!, :loc_y], locs[!, :loc_z] = x, y, z
    return nothing
end

function _locs_norm(obj::NeuroAnalyzer.NEURO)
    obj_new = deepcopy(obj)
    _locs_norm!(obj_new.locs)
    return obj_new
end

function _locs_norm!(obj::NeuroAnalyzer.NEURO)
    _locs_norm!(obj.locs)
    return nothing
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

_xyz2r(x, y, z) = sqrt(x^2 + y^2 + z^2)
