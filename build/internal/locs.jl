"""
    _sph_distance_sph(r1, theta1, phi1, r2, theta2, phi2)

Calculate the 3-D distance between two points given in spherical coordinates.

# Arguments

- `r1`, `r2::Real`: radial distances from the origin
- `theta1`, `theta2::Real`: polar angles in degrees from the z-axis (0–180°)
- `phi1`, `phi2::Real`: azimuthal angles in degrees in the x-y plane (0–360°)

# Returns

- `Float64`: Euclidean distance between the two points, in the same units as `r1`/`r2`

# Notes

Uses the spherical law of cosines:

    d = √(r₁² + r₂² − 2·r₁·r₂·[cos θ₁ cos θ₂ + sin θ₁ sin θ₂ cos(φ₁−φ₂)])
"""
function _sph_distance_sph(
    r1::Real,
    theta1::Real,
    phi1::Real,
    r2::Real,
    theta2::Real,
    phi2::Real,
)::Float64
    # convert angles to radians for calculation
    θ1, θ2 = deg2rad(theta1), deg2rad(theta2)
    φ1, φ2 = deg2rad(phi1), deg2rad(phi2)
    # calculate spherical distance using law of cosines
    return sqrt(
        r1^2 + r2^2 -
        2 * r1 * r2 * (cos(θ1) * cos(θ2) + sin(θ1) * sin(θ2) * cos(φ1 - φ2)),
    )
end

"""
    _sph_distance_cart(x1, y1, z1, x2, y2, z2)

Calculate the Euclidean distance between two points in Cartesian coordinates.

# Arguments

- `x1`, `y1`, `z1::Real`: coordinates of the first point
- `x2`, `y2`, `z2::Real`: coordinates of the second point

# Returns

- `Float64`: Euclidean distance between the two points
"""
function _sph_distance_cart(
    x1::Real,
    y1::Real,
    z1::Real,
    x2::Real,
    y2::Real,
    z2::Real,
)::Float64
    return sqrt((x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2)
end

"""
    _check_ch_locs(ch, objl, locsl)

Validate that specified channels all have corresponding location entries.

Throws `ArgumentError` if any channel index is out of bounds or its label is absent from the locations table.

# Arguments

- `ch::Union{Int64, Vector{Int64}}`: channel index or indices to validate
- `objl::Vector{String}`: channel labels from the NEURO object
- `locsl::Vector{String}`: labels present in the locations table
"""
function _check_ch_locs(
    ch::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}},
    objl::Vector{String},
    locsl::Vector{String},
)::Nothing
    isempty(ch) && throw(ArgumentError("Channel specification cannot be empty."))
    isempty(objl) && throw(ArgumentError("Channel labels vector cannot be empty."))
    isempty(locsl) && throw(ArgumentError("Location labels vector cannot be empty."))

    indices = ch isa Int64 ? [ch] : ch
    for idx in indices
        1 ≤ idx ≤ length(objl) || throw(
            ArgumentError("Channel index $idx is out of bounds (1 to $(length(objl)))."),
        )
        objl[idx] in locsl ||
            throw(ArgumentError("Channel $(objl[idx]) does not have a location."))
    end
    return nothing
end

"""
    _loc_idx(obj, ch)

Return the row index (or indices) in `obj.locs` for the given channel(s).

Accepts integer channel indices.
"""
_loc_idx(
    obj::NeuroAnalyzer.NEURO,
    ch::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}},
)::Union{Int64, Vector{Int64}} = _find_bylabel(obj.locs, labels(obj)[ch])

"""
    _loc_idx(obj, ch)

Return the row index (or indices) in `obj.locs` for the given channel(s).

Accepts channel name(s)/regex.
"""
_loc_idx(
    obj::NeuroAnalyzer.NEURO,
    ch::Union{String, Vector{String}, Regex},
)::Vector{Int64} = _find_bylabel(obj.locs, labels(obj)[get_channel(obj; ch = ch)])

"""
    _idx2lab(obj, ch)

Return the location label(s) in `obj.locs` corresponding to channel index/indices `ch`.
"""
_idx2lab(
    obj::NeuroAnalyzer.NEURO,
    ch::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}},
)::Union{String, Vector{String}} = obj.locs[
    _loc_idx(obj, ch), :label,
]

"""
    _idx2lab(obj, ch)

Return the location label(s) in `obj.locs` corresponding to channel name(s) `ch`.
"""
_idx2lab(
    obj::NeuroAnalyzer.NEURO,
    ch::Union{String, Vector{String}, Regex},
)::Vector{String} = obj.locs[_loc_idx(obj, ch), :label]

"""
    _ch_locs(obj, ch)

Return a filtered `obj.locs` DataFrame containing only the rows for channel(s) `ch`.

Validates that all requested channels have location data before returning.

Accepts integer indices.
"""
function _ch_locs(
    obj::NeuroAnalyzer.NEURO,
    ch::Union{Int64, Vector{Int64}, AbstractUnitRange{Int64}},
)::DataFrame
    chl  = labels(obj)[ch]
    chs  = intersect(obj.locs[!, :label], chl)
    locs = Base.filter(:label => in(chs), obj.locs)
    _check_ch_locs(ch, labels(obj), obj.locs[!, :label])
    return locs
end

"""
    _ch_locs(obj, ch)

Return a filtered `obj.locs` DataFrame containing only the rows for channel(s) `ch`.

Validates that all requested channels have location data before returning.

Accepts channel names/regex.
"""
_ch_locs(
    obj::NeuroAnalyzer.NEURO,
    ch::Union{String, Vector{String}, Regex},
)::DataFrame = _ch_locs(obj, get_channel(obj; ch = ch))

"""
    _find_bylabel(locs, l)

Find row index (or indices) in a locations DataFrame by label, case-insensitively.

# Arguments

- `locs::DataFrame`: channel locations DataFrame containing a `:label` column
- `l`:
    - `String`: returns the first matching index as `Int64`, or empty `Int64[]` if not found
    - `Vector{String}`: multiple labels to search for
    - `Vector{String}` / `Vector{SubString{String}}`: returns a `Vector{Int64}` of matching indices

# Returns

- if `l` is a `String`: returns the first matching index as `Int64`, or empty `Int64[]` if not found
- if `l` is a `Vector`: returns a `Vector{Int64}` of all matching indices
"""
function _find_bylabel(
    locs::DataFrame,
    l::Union{String, Vector{String}, Vector{SubString{String}}},
)::Union{Int64, Vector{Int64}}
    "label" in names(locs) ||
        throw(ArgumentError("Locations DataFrame must contain a \"label\" column."))

    # convert labels to lowercase for case-insensitive comparison
    loc_labels = lowercase.(locs.label)

    if l isa String
        # single label case - return first match index or empty vector
        match_idx = findfirst(==(lowercase(l)), loc_labels)
        return isnothing(match_idx) ? Int64[] : Int64(match_idx)
    else
        # multiple labels case - return vector of all matches
        l_idx = Int64[]
        for label in l
            match_idx = findfirst(==(lowercase(label)), loc_labels)
            isnothing(match_idx) || push!(l_idx, Int64(match_idx))
        end
        return l_idx
    end
end

"""
    _initialize_locs()

Return an empty channel locations DataFrame with the standard column schema: `label`, `loc_radius`, `loc_theta`, `loc_x`, `loc_y`, `loc_z`, `loc_radius_sph`, `loc_theta_sph`, `loc_phi_sph`.
"""
function _initialize_locs()::DataFrame
    return DataFrame(
        :label          => String[],
        :loc_radius     => Float64[],
        :loc_theta      => Float64[],
        :loc_x          => Float64[],
        :loc_y          => Float64[],
        :loc_z          => Float64[],
        :loc_radius_sph => Float64[],
        :loc_theta_sph  => Float64[],
        :loc_phi_sph    => Float64[],
    )
end

"""
    _initialize_locs(obj)

Return a channel locations DataFrame pre-populated for all locatable channels in `obj`.

Covers channel types: `eeg`, `meg`, `grad`, `mag`, `ecog`, `seeg`, `ieeg`, `nirs_int`, `nirs_od`, `eog`, `ref`.

All coordinate columns are initialised to 0.0.
"""
function _initialize_locs(obj::NeuroAnalyzer.NEURO)::DataFrame
    locs_ch = get_channel(
        obj;
        ch = get_channel(
            obj;
            type = ["meg", "grad", "mag", "eeg", "ecog", "seeg", "ieeg",
                "nirs_int", "nirs_od", "eog", "ref"],
        ),
    )
    n = length(locs_ch)
    return DataFrame(
        :label          => labels(obj)[locs_ch],
        :loc_radius     => zeros(n),
        :loc_theta      => zeros(n),
        :loc_x          => zeros(n),
        :loc_y          => zeros(n),
        :loc_z          => zeros(n),
        :loc_radius_sph => zeros(n),
        :loc_theta_sph  => zeros(n),
        :loc_phi_sph    => zeros(n),
    )
end

"""
    _initialize_locs!(obj)

Initialise `obj.locs` in-place with zero-valued location entries for all locatable
channels. See `_initialize_locs(obj)` for the list of covered channel types.
"""
function _initialize_locs!(obj::NeuroAnalyzer.NEURO)::Nothing
    obj.locs = _initialize_locs(obj)
    return nothing
end

"""
    _check_locs_cols(locs)
 
Throw `ArgumentError` if any expected coordinate column is absent from `locs`.
"""
function _check_locs_cols(locs::DataFrame)::Nothing
    for col in _LOCS_COORD_COLS
        col in names(locs) ||
            throw(ArgumentError("Locations DataFrame must contain a \":$col\" column."))
    end
    return nothing
end

"""
    _locs_round(locs)

Return a copy of `locs` with all coordinate columns rounded to 2 decimal places.
"""
function _locs_round(locs::DataFrame)::DataFrame
    _check_locs_cols(locs)
    locs_new = copy(locs)
    for col in _LOCS_COORD_COLS
        locs_new[!, col] = round.(locs[!, col]; digits = 2)
    end
    return locs_new
end

"""
    _locs_round!(locs)

Round all coordinate columns in `locs` to 2 decimal places in-place.
"""
function _locs_round!(locs::DataFrame)::Nothing
    _check_locs_cols(locs)
    for col in _LOCS_COORD_COLS
        locs[!, col] = round.(locs[!, col]; digits = 2)
    end
    return nothing
end

"""
    _locs_round(obj)

Return `obj.locs` with all coordinate columns rounded to 2 decimal places.
"""
_locs_round(obj::NeuroAnalyzer.NEURO)::DataFrame = _locs_round(obj.locs)

"""
    _locs_round!(obj)

Round all coordinate columns in `obj.locs` to 2 decimal places in-place.
"""
function _locs_round!(obj::NeuroAnalyzer.NEURO)::Nothing
    _locs_round!(obj.locs)
    return nothing
end

"""
    _locs_remove_nans(locs)

Return a copy of `locs` with all `NaN` values in numeric columns replaced by `0.0`.
"""
function _locs_remove_nans(locs::DataFrame)::DataFrame
    isempty(locs) &&
        throw(ArgumentError("Locations DataFrame cannot be empty."))
    locs_new = copy(locs)
    for col in names(locs_new)
        eltype(locs_new[!, col]) <: AbstractFloat || continue
        locs_new[!, col] = replace(locs_new[!, col], NaN => 0.0)
    end
    return locs_new
end

"""
    _locs_remove_nans!(locs)

Replace all `NaN` values in numeric columns of `locs` with `0.0` in-place.
"""
function _locs_remove_nans!(locs::DataFrame)::Nothing
    isempty(locs) &&
        throw(ArgumentError("Locations DataFrame cannot be empty."))
    for col in names(locs)
        eltype(locs[!, col]) <: AbstractFloat || continue
        locs[!, col] = replace(locs[!, col], NaN => 0.0)
    end
    return nothing
end

"""
    _has_locs(obj)

Throw `ArgumentError` if `obj.locs` is empty.

Use this as a guard before any function that requires location data to be loaded.
"""
function _has_locs(obj::NeuroAnalyzer.NEURO)::Nothing
    isempty(obj.locs) && throw(
        ArgumentError(
            "Channel locations is empty, use load_locs() or add_locs() first.",
        ),
    )
    return nothing
end

"""
    _locs_norm(x, y)

Normalize location coordinates to the [0, 1] range using min-max scaling.
"""
function _locs_norm(
    x::Union{AbstractVector, Real},
    y::Union{AbstractVector, Real},
)::Tuple{Vector{Float64}, Vector{Float64}}
    xy = normalize_minmax(hcat(x, y))
    return xy[:, 1], xy[:, 2]
end

"""
    _locs_norm(x, y, z)

Normalize location coordinates to the [0, 1] range using min-max scaling.
"""
function _locs_norm(
    x::Union{AbstractVector, Real},
    y::Union{AbstractVector, Real},
    z::Union{AbstractVector, Real},
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}
    xyz = normalize_minmax(hcat(x, y, z))
    return xyz[:, 1], xyz[:, 2], xyz[:, 3]
end

"""
    _locs_norm(locs)

Normalize location coordinates to the [0, 1] range using min-max scaling.
"""
function _locs_norm(locs::DataFrame)::DataFrame
    locs_new = copy(locs)
    xyz = normalize_minmax(hcat(locs.loc_x, locs.loc_y, locs.loc_z))
    x = xyz[:, 1]
    y = xyz[:, 2]
    z = xyz[:, 3]
    locs_new.loc_x, locs_new.loc_y, locs_new.loc_z = x, y, z
    return locs_new
end

"""
    _locs_norm!(locs)

Normalize location coordinates to the [0, 1] range using min-max scaling in-place.
"""
function _locs_norm!(locs::DataFrame)::Nothing
    xyz = normalize_minmax(hcat(locs.loc_x, locs.loc_y, locs.loc_z))
    x = xyz[:, 1]
    y = xyz[:, 2]
    z = xyz[:, 3]
    locs.loc_x, locs.loc_y, locs.loc_z = x, y, z
    return nothing
end

"""
    _locs_norm(obj)

Normalize location coordinates to the [0, 1] range using min-max scaling.
"""
function _locs_norm(obj::NeuroAnalyzer.NEURO)::NeuroAnalyzer.NEURO
    # create new dataset
    obj_new = deepcopy(obj)
    _locs_norm!(obj_new.locs)
    return obj_new
end

"""
    _locs_norm!(obj)

Normalize location coordinates to the [0, 1] range using min-max scaling in-place.
"""
function _locs_norm!(obj::NeuroAnalyzer.NEURO)::Nothing
    _locs_norm!(obj.locs)
    return nothing
end

"""
    _angle_quadrant(a)

Return the quadrant (1–4) of angle `a` given in degrees.

Quadrants follow the standard mathematical convention:
- Q1:   0° to  90°  (positive x, positive y)
- Q2:  90° to 180°  (negative x, positive y)
- Q3: 180° to 270°  (negative x, negative y)
- Q4: 270° to 360°  (positive x, negative y)

Angles outside [0°, 360°) are first normalized via `mod`.
"""
function _angle_quadrant(a::Real)::Int64
    # normalize to [0°, 360°) regardless of sign
    a = mod(Float64(a), 360.0)
    a <= 90 && return 1
    a <= 180 && return 2
    a <= 270 && return 3
    return 4
end

"""
    _xyz2r(x, y, z)

Return the spherical radial distance `r = √(x² + y² + z²)`.
"""
_xyz2r(x::Real, y::Real, z::Real)::Float64 = sqrt(x^2 + y^2 + z^2)

"""
    _midxy(x1, y1, x2, y2)

Return the midpoint between two 2-D points as `(xm, ym)`.
"""
_midxy(x1::Real, y1::Real, x2::Real, y2::Real)::Tuple{Float64, Float64} =
    (x1 + ((x2 - x1) / 2), y1 + ((y2 - y1) / 2))
