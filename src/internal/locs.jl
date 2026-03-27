"""
    _sph_distance_sph(r1, theta1, phi1, r2, theta2, phi2)

Calculate spherical distance between two points in spherical coordinates.

# Arguments

- `r1`, `r2`: radial distances from origin
- `theta1`, `theta2`: polar angles (degrees) from z-axis (0-180°)
- `phi1`, `phi2`: azimuthal angles (degrees) in x-y plane (0-360°)

# Returns

- `Float64`: distance between points in same units as inputs

# Formula

Uses the spherical law of cosines:

    d = √(r₁² + r₂² - 2·r₁·r₂·[cos(θ₁)cos(θ₂) + sin(θ₁)sin(θ₂)cos(φ₁-φ₂)])
"""
function _sph_distance_sph(
    r1::Real,
    theta1::Real,
    phi1::Real,
    r2::Real,
    theta2::Real,
    phi2::Real,
)
    # convert angles to radians for calculation
    θ1, θ2 = deg2rad(theta1), deg2rad(theta2)
    φ1, φ2 = deg2rad(phi1), deg2rad(phi2)
    # calculate spherical distance using law of cosines
    d = sqrt(
        r1^2 + r2^2 -
        2 * r1 * r2 * (cos(θ1) * cos(θ2) + sin(θ1) * sin(θ2) * cos(φ1 - φ2)),
    )
    return d
end

"""Calculate Euclidean distance between two points in Cartesian coordinates."""
function _sph_distance_cart(x1::Real, y1::Real, z1::Real, x2::Real, y2::Real, z2::Real)
    return sqrt((x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2)
end

"""
    _check_ch_locs(ch, objl, locsl)

Validate that specified channels have corresponding location data available.

# Arguments
- `ch::Union{Int64, Vector{Int64}}`: channel index/indices to validate
- `objl::Vector{String}`: vector of channel labels from the NEURO object
- `locsl::Vector{String}`: vector of available location labels

# Returns

- `Nothing` if all channels have valid locations
"""
function _check_ch_locs(
    ch::Union{Int64, Vector{Int64}},
    objl::Vector{String},
    locsl::Vector{String},
)::Nothing
    isempty(ch) || throw(ArgumentError("Channel specification cannot be empty."))
    isempty(objl) || throw(ArgumentError("Channel labels vector cannot be empty."))
    isempty(locsl) || throw(ArgumentError("Location labels vector cannot be empty."))
    if ch isa Int64
        1 ≤ ch ≤ length(objl) ||
            throw(ArgumentError("Channel index $ch is out of bounds (1 to $(length(objl)))."))
        objl[ch] in locsl ||
            throw(ArgumentError("Channel $(objl[ch]) does not have a location"))
    else
        for idx in ch
            1 ≤ ch ≤ length(objl) ||
                throw(ArgumentError("Channel index $idx is out of bounds (1 to $(length(objl)))."))
            objl[idx] in locsl ||
                throw(ArgumentError("Channel $(objl[idx]) does not have a location"))
        end
    end
    return nothing
end

"""Find channel location using its index."""
_loc_idx(
    obj::NeuroAnalyzer.NEURO,
    ch::Union{Int64, Vector{Int64}},
)::Union{Int64, Vector{Int64}} = _find_bylabel(
    obj.locs, labels(obj)[ch],
)

"""Find channel location using its name."""
_loc_idx(
    obj::NeuroAnalyzer.NEURO,
    ch::Union{String, Vector{String}, Regex},
)::Vector{Int64} = _find_bylabel(
    obj.locs, labels(obj)[get_channel(obj; ch = ch)],
)

"""Find channel location in OBJ.locs using its index."""
_idx2lab(
    obj::NeuroAnalyzer.NEURO,
    ch::Union{Int64, Vector{Int64}},
)::Union{String, Vector{String}} = obj.locs[
    _loc_idx(obj, ch), :label,
]

_idx2lab(
    obj::NeuroAnalyzer.NEURO,
    ch::Union{String, Vector{String}, Regex},
)::Vector{String} = obj.locs[
    _loc_idx(obj, ch), :label,
]

function _ch_locs(obj::NeuroAnalyzer.NEURO, ch::Union{Int64, Vector{Int64}})::DataFrame
    chl = labels(obj)[ch]
    chs = intersect(obj.locs[!, :label], chl)
    locs = Base.filter(:label => in(chs), obj.locs)
    _check_ch_locs(ch, labels(obj), obj.locs[!, :label])
    return locs
end

function _ch_locs(
    obj::NeuroAnalyzer.NEURO,
    ch::Union{String, Vector{String}, Regex},
)::DataFrame
    return _ch_locs(obj, get_channel(obj; ch = ch))
end

"""
    _find_bylabel(locs, l)

Find location indices by label with case-insensitive matching.

# Arguments

- `locs::DataFrame`: channel locations DataFrame containing a `:label` column
- `l`:
    - `String`: single label to search for
    - `Vector{String}`: multiple labels to search for
    - `Vector{SubString{String}}`: multiple labels (substring type) to search for

# Returns

- if `l` is a `String`: returns the first matching index as `Int64`, or empty `Int64[]` if not found
- if `l` is a `Vector`: returns a `Vector{Int64}` of all matching indices
"""
function _find_bylabel(
    locs::DataFrame,
    l::Union{String, Vector{String}, Vector{SubString{String}}},
)::Union{Int64, Vector{Int64}}
    :label in names(locs) ||
        throw(ArgumentError("Channel locations DataFrame must contain a \":label\" column."))
    # convert labels to lowercase for case-insensitive comparison
    loc_labels = lowercase.(locs[!, ])
    if l isa String
        # single label case - return first match index or empty vector
        match_idx = findfirst(==(lowercase(l)), loc_labels)
        return isnothing(match_idx) ? Int64[] : Int64(match_idx)
    else
        # multiple labels case - return vector of all matches
        l_idx = Int64[]
        for label in l
            match_idx = findfirst(==(lowercase(label)), loc_labels)
            if !isnothing(match_idx)
                push!(l_idx, Int64(match_idx))
            end
        end
        return l_idx
    end
end

"""Initialize an empty channel locations DataFrame."""
function _initialize_locs()::DataFrame
    return DataFrame(
        :label => String[],
        :loc_radius => Float64[],
        :loc_theta => Float64[],
        :loc_x => Float64[],
        :loc_y => Float64[],
        :loc_z => Float64[],
        :loc_radius_sph => Float64[],
        :loc_theta_sph => Float64[],
        :loc_phi_sph => Float64[],
    )
end

"""
Initialize in-place an empty channel locations DataFrame for a NEURO object, containing entries for all channels of standard types.

Initializes locations for these channel types:

- EEG, MEG (magnetometers), MEG (gradometers)
- ECoG, sEEG, iEEG
- NIRS (intensity and optical density)
- EOG, reference channels

All location values (x, y, z, spherical coordinates) are initialized to 0.0.

Channel labels are preserved from the original object.
"""
function _initialize_locs(obj::NeuroAnalyzer.NEURO)::DataFrame
    locs_ch = get_channel(
        obj;
        ch = get_channel(
            obj;
            type = [
                "meg",
                "grad",
                "mag",
                "eeg",
                "ecog",
                "seeg",
                "ieeg",
                "nirs_int",
                "nirs_od",
                "eog",
                "ref",
            ],
        ),
    )
    return DataFrame(
        :label => labels(obj)[locs_ch],
        :loc_radius => zeros(length(locs_ch)),
        :loc_theta => zeros(length(locs_ch)),
        :loc_x => zeros(length(locs_ch)),
        :loc_y => zeros(length(locs_ch)),
        :loc_z => zeros(length(locs_ch)),
        :loc_radius_sph => zeros(length(locs_ch)),
        :loc_theta_sph => zeros(length(locs_ch)),
        :loc_phi_sph => zeros(length(locs_ch)),
    )
end

"""
Initialize in-place an empty channel locations DataFrame for a NEURO object, containing entries for all channels of standard types.

Initializes locations for these channel types:

- EEG, MEG (magnetometers), MEG (gradometers)
- ECoG, sEEG, iEEG
- NIRS (intensity and optical density)
- EOG, reference channels

All location values (x, y, z, spherical coordinates) are initialized to 0.0.

Channel labels are preserved from the original object.
"""
function _initialize_locs!(obj::NeuroAnalyzer.NEURO)::Nothing
    locs_ch = get_channel(
        obj;
        ch = get_channel(
            obj;
            type = [
                "meg",
                "grad",
                "mag",
                "eeg",
                "ecog",
                "seeg",
                "ieeg",
                "nirs_int",
                "nirs_od",
                "eog",
                "ref",
            ],
        ),
    )
    obj.locs = DataFrame(
        :label => labels(obj)[locs_ch],
        :loc_radius => zeros(length(locs_ch)),
        :loc_theta => zeros(length(locs_ch)),
        :loc_x => zeros(length(locs_ch)),
        :loc_y => zeros(length(locs_ch)),
        :loc_z => zeros(length(locs_ch)),
        :loc_radius_sph => zeros(length(locs_ch)),
        :loc_theta_sph => zeros(length(locs_ch)),
        :loc_phi_sph => zeros(length(locs_ch)),
    )
    return nothing
end

"""
Round all location coordinate values in a channel locations DataFrame to 2 decimal places.
"""
function _locs_round(locs::DataFrame)::DataFrame
    cols = ["loc_x", "loc_y", "loc_z", "loc_radius", "loc_theta", "loc_radius_sph", "loc_theta_sph", "loc_phi_sph"]
    for col in cols
        col in names(locs) ||
            throw(ArgumentError("Location DataFrame must contain \":$col\" column"))
    end
    locs_new = copy(locs)
    for col in cols
        locs_new[!, col] = round.(getproperty(locs, col); digits=2)
    end
    return locs_new
end

"""
Round all location coordinate values in-place in a channel locations DataFrame to 2 decimal places.
"""
function _locs_round!(locs::DataFrame)::Nothing
    cols = ["loc_x", "loc_y", "loc_z", "loc_radius", "loc_theta", "loc_radius_sph", "loc_theta_sph", "loc_phi_sph"]
    for col in cols
        col in names(locs) ||
            throw(ArgumentError("Location DataFrame must contain \":$col\" column"))
    end
    for col in cols
        locs[!, col] = round.(getproperty(locs, col); digits=2)
    end
    return nothing
end

"""
Round all location coordinate values to 2 decimal places.
"""
_locs_round(obj::NeuroAnalyzer.NEURO)::DataFrame = _locs_round(obj.locs)

"""
Round all location coordinate values in-place to 2 decimal places.
"""
function _locs_round!(obj::NeuroAnalyzer.NEURO)::Nothing
    obj.locs = _locs_round(obj.locs)
    return nothing
end

"""
Replace all NaN values with zeros in a channel locations DataFrame.
"""
function _locs_remove_nans(locs::DataFrame)::DataFrame
    isempty(locs) &&
        throw(ArgumentError("Channel locations DataFrame cannot be empty."))
    any(eltype.(eachcol(locs)) .<: Real) ||
        throw(ArgumentError("Channel locations DataFrame must contain at least one numeric column."))
    locs_new = copy(locs)
    # replace NaN with 0.0 in all columns
    for col in names(locs)[eltype.(eachcol(locs)) .<: Real]
        replace!(collect(skipmissing(locs_new[!, col])), NaN => 0.0)
    end
    return locs_new
end

"""
Replace in-place all NaN values with zeros in a channel locations DataFrame.
"""
function _locs_remove_nans!(locs::DataFrame)::Nothing
    isempty(locs) &&
        throw(ArgumentError("Channel locations DataFrame cannot be empty."))
    any(eltype.(eachcol(locs)) .<: Real) ||
        throw(ArgumentError("Channel locations DataFrame must contain at least one numeric column."))
    # replace NaN with 0.0 in all columns
    for col in names(locs)[eltype.(eachcol(locs)) .<: Real]
        replace!(collect(skipmissing(locs[!, col])), NaN => 0.0)
    end
    return nothing
end

"""
Check whether OBJ.locs is empty.
"""
function _has_locs(obj::NeuroAnalyzer.NEURO)::Nothing
    isempty(obj.locs) && throw(
        ArgumentError(
            "Channel locations is empty, use load_locs() or add_locs() first.",
        ),
    )
    return nothing
end

function _locs_norm(
    x::Union{AbstractVector, Real}, y::Union{AbstractVector, Real},
)::Tuple{Vector{Float64}, Vector{Float64}}
    xy = normalize_minmax(hcat(x, y))
    x = xy[:, 1]
    y = xy[:, 2]
    return x, y
end

function _locs_norm(
    x::Union{AbstractVector, Real}, y::Union{AbstractVector, Real},
    z::Union{AbstractVector, Real},
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}
    xyz = normalize_minmax(hcat(x, y, z))
    x = xyz[:, 1]
    y = xyz[:, 2]
    z = xyz[:, 3]
    return x, y, z
end

function _locs_norm(locs::DataFrame)::DataFrame
    locs_new = copy(locs)
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
    # create new dataset
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

_midxy(x1::Real, y1::Real, x2::Real, y2::Real)::Tuple{Float64, Float64} =
    (x1 + ((x2 - x1) / 2), y1 + ((y2 - y1) / 2))
