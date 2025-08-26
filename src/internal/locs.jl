function _sph_distance_sph(r1::Real, theta1::Real, phi1::Real, r2::Real, theta2::Real, phi2::Real)
    d = sqrt(r1^2 + r2^2 - (2 * r1 * r2) * cosd(theta1 - theta2) + (2 * r1 * r2) * sind(theta1) * sind(theta2) * (cosd(phi1 - phi2 - 1)))
    return d
end
function _sph_distance_cart(x1::Real, y1::Real, z1::Real, x2::Real, y2::Real, z2::Real)
    d = sqrt((x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2)
    return d
end

function _check_ch_locs(ch::Union{Int64, Vector{Int64}}, objl::Vector{String}, locsl::Vector{String})::Nothing
    for idx in ch
        @assert objl[idx] in locsl "Channel $(objl[idx]) does not have a location."
    end
    return nothing
end

_loc_idx(obj::NeuroAnalyzer.NEURO, ch::Union{Int64, Vector{Int64}})::Union{Int64, Vector{Int64}} = _find_bylabel(obj.locs, labels(obj)[ch])
_loc_idx(obj::NeuroAnalyzer.NEURO, ch::Union{String, Vector{String}, Regex})::Vector{Int64} = _find_bylabel(obj.locs, labels(obj)[get_channel(obj, ch=ch)])
_idx2lab(obj::NeuroAnalyzer.NEURO, ch::Union{Int64, Vector{Int64}})::Union{String, Vector{String}} = obj.locs[_loc_idx(obj, ch), :label]
_idx2lab(obj::NeuroAnalyzer.NEURO, ch::Union{String, Vector{String}, Regex})::Vector{String} = obj.locs[NeuroAnalyzer._loc_idx(obj, ch), :label]

function _ch_locs(obj::NeuroAnalyzer.NEURO, ch::Union{Int64, Vector{Int64}})::DataFrame
    chl = labels(obj)[ch]
    chs = intersect(obj.locs[!, :label], chl)
    locs = Base.filter(:label => in(chs), obj.locs)
    _check_ch_locs(ch, labels(obj), obj.locs[!, :label])
    return locs
end

function _ch_locs(obj::NeuroAnalyzer.NEURO, ch::Union{String, Vector{String}, Regex})::DataFrame
    return _ch_locs(obj, get_channel(obj, ch=ch))
end

function _find_bylabel(locs::DataFrame, l::Union{String, Vector{String}, Vector{SubString{String}}})::Union{Int64, Vector{Int64}}
    if isa(l, String)
        if !isnothing(findfirst(isequal.(lowercase(l), lowercase.(locs[!, :label]))))
            return findfirst(isequal.(lowercase(l), lowercase.(locs[!, :label])))
        else
            return Int64[]
        end
    else
        l_idx = Vector{Int64}()
        for idx in l
            lowercase(idx) in lowercase.(locs[!, :label]) && push!(l_idx, findfirst(isequal.(lowercase(idx), lowercase.(locs[!, :label]))))
        end
        return l_idx
    end
end

function _initialize_locs()::DataFrame
    return DataFrame(:label=>String[],
                     :loc_radius=>Float64[],
                     :loc_theta=>Float64[],
                     :loc_x=>Float64[],
                     :loc_y=>Float64[],
                     :loc_z=>Float64[],
                     :loc_radius_sph=>Float64[],
                     :loc_theta_sph=>Float64[],
                     :loc_phi_sph=>Float64[])
end

function _initialize_locs!(obj::NeuroAnalyzer.NEURO)::Nothing
    locs_ch = get_channel(obj, ch=get_channel(obj, type=["meg", "grad", "mag", "eeg", "ecog", "seeg", "ieeg", "nirs_int", "nirs_od", "eog", "ref"]))
    obj.locs = DataFrame(:label=>labels(obj)[locs_ch], :loc_radius=>zeros(length(locs_ch)), :loc_theta=>zeros(length(locs_ch)), :loc_x=>zeros(length(locs_ch)), :loc_y=>zeros(length(locs_ch)), :loc_z=>zeros(length(locs_ch)), :loc_radius_sph=>zeros(length(locs_ch)), :loc_theta_sph=>zeros(length(locs_ch)), :loc_phi_sph=>zeros(length(locs_ch)))
    return nothing
end

function _initialize_locs(obj::NeuroAnalyzer.NEURO)::DataFrame
    locs_ch = get_channel(obj, ch=get_channel(obj, type=datatype(obj)))
    return DataFrame(:label=>labels(obj)[locs_ch], :loc_radius=>zeros(length(locs_ch)), :loc_theta=>zeros(length(locs_ch)), :loc_x=>zeros(length(locs_ch)), :loc_y=>zeros(length(locs_ch)), :loc_z=>zeros(length(locs_ch)), :loc_radius_sph=>zeros(length(locs_ch)), :loc_theta_sph=>zeros(length(locs_ch)), :loc_phi_sph=>zeros(length(locs_ch)))
end

function _locs_round(locs::DataFrame)::DataFrame
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

function _locs_round!(locs::DataFrame)::Nothing
    locs[!, :loc_radius] = round.(locs[!, :loc_radius], digits=2)
    locs[!, :loc_theta] = round.(locs[!, :loc_theta], digits=2)
    locs[!, :loc_x] = round.(locs[!, :loc_x], digits=2)
    locs[!, :loc_y] = round.(locs[!, :loc_y], digits=2)
    locs[!, :loc_z] = round.(locs[!, :loc_z], digits=2)
    locs[!, :loc_radius_sph] = round.(locs[!, :loc_radius_sph], digits=2)
    locs[!, :loc_theta_sph] = round.(locs[!, :loc_theta_sph], digits=2)
    locs[!, :loc_phi_sph] = round.(locs[!, :loc_phi_sph], digits=2)
    return nothing
end

function _locs_remove_nans(locs::DataFrame)::DataFrame
    locs_new = deepcopy(locs)
    for idx in eachcol(locs_new)
        replace!(idx, NaN => 0)
    end
    return locs_new
end

function _locs_remove_nans!(locs::DataFrame)::Nothing
    for idx in eachcol(locs)
        replace!(idx, NaN => 0)
    end
    return nothing
end

_locs_round(obj::NeuroAnalyzer.NEURO)::DataFrame = _locs_round(obj.locs)

function _locs_round!(obj::NeuroAnalyzer.NEURO)::Nothing
    obj.locs = _locs_round(obj.locs)
    return nothing
end

function _has_locs(obj::NeuroAnalyzer.NEURO)::Nothing
    @assert nrow(obj.locs) > 0 "Electrode locations not available, use load_locs() or add_locs() first."
    return nothing
end

function _locs_norm(x::Union{AbstractVector, Real}, y::Union{AbstractVector, Real})::Tuple{Vector{Float64}, Vector{Float64}}
    xy = normalize_minmax(hcat(x, y))
    x = xy[:, 1]
    y = xy[:, 2]
    return x, y
end

function _locs_norm(x::Union{AbstractVector, Real}, y::Union{AbstractVector, Real}, z::Union{AbstractVector, Real})::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}
    xyz = normalize_minmax(hcat(x, y, z))
    x = xyz[:, 1]
    y = xyz[:, 2]
    z = xyz[:, 3]
    return x, y, z
end

function _locs_norm(locs::DataFrame)::DataFrame
    locs_new = deepcopy(locs)
    x, y, z = locs[!, :loc_x], locs[!, :loc_y], locs[!, :loc_z]
    xyz = normalize_minmax(hcat(x, y, z))
    x = xyz[:, 1]
    y = xyz[:, 2]
    z = xyz[:, 3]
    locs_new[!, :loc_x], locs_new[!, :loc_y], locs_new[!, :loc_z] = x, y, z
    return locs_new
end

function _locs_norm!(locs::DataFrame)::Nothing
    x, y, z = locs[!, :loc_x], locs[!, :loc_y], locs[!, :loc_z]
    xyz = normalize_minmax(hcat(x, y, z))
    x = xyz[:, 1]
    y = xyz[:, 2]
    z = xyz[:, 3]
    locs[!, :loc_x], locs[!, :loc_y], locs[!, :loc_z] = x, y, z
    return nothing
end

function _locs_norm(obj::NeuroAnalyzer.NEURO)::NeuroAnalyzer.NEURO
    obj_new = deepcopy(obj)
    _locs_norm!(obj_new.locs)
    return obj_new
end

function _locs_norm!(obj::NeuroAnalyzer.NEURO)::Nothing
    _locs_norm!(obj.locs)
    return nothing
end

function _angle_quadrant(a::Real)::Int64
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

_xyz2r(x::Real, y::Real, z::Real)::Float64 = sqrt(x^2 + y^2 + z^2)

_midxy(x1::Real, y1::Real, x2::Real, y2::Real)::Tuple{Float64, Float64} = (x1 + ((x2 - x1) / 2), y1 + ((y2 - y1) / 2))