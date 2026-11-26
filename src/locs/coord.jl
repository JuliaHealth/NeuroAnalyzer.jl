export locs_coord
export set_locs_coord
export set_locs_coord!

# units of channel locations (BIDS: m, cm, mm; "": unknown, or normalized, e.g. unit sphere)
const _LOCS_UNITS = ("m", "cm", "mm", "")

"""
    locs_coord(obj)

Coordinate system and unit of the channel locations (`obj.locs`).

They are stored in the header (`:coord_system`, `:length_unit`) when they are known, e.g. by `import_bids()` (`coordsystem.json`), the FIFF, FieldTrip, EEGLAB, BrainVision and NIRS importers, `load_locs!()` (ELC files) or `set_locs_coord!()`, together with the locations they describe. When the locations are changed afterwards, the stored values may no longer apply:

- the unit, when the scale of the locations changes (the median distance of the channels from the origin changes more than 2×, e.g. `locs_scale!()`, `locs_normalize!()`, or locations of another unit loaded by `load_locs!()`)
- the coordinate system, when most channels (more than half) have moved (e.g. `locs_rotz!()`, `locs_center!()`, `load_locs!()`)

Editing the positions of a few channels, deleting, adding or renaming channels keeps the stored values. A coordinate system that no longer applies is unknown (`""`); a unit that is not stored, unknown or no longer applies is guessed from the size of the head (median distance of the channels from the origin, assuming a head radius of about 8–10 cm):

| Distance        | Unit                                |
| --------------- | ----------------------------------- |
| 0.02–0.3        | `"m"`                               |
| 2.5–25          | `"cm"`                              |
| 25–300          | `"mm"`                              |
| other (e.g. 1)  | `""` (unknown, e.g. unit sphere)    |

The coordinates of `obj.locs` are always in the NeuroAnalyzer orientation (x: right, y: nose, z: up); the coordinate system (e.g. `"CapTrak"`, `"EEGLAB"`, `"MNI305"`, BIDS names) describes their origin and scale, and the axes of the original data (e.g. `"EEGLAB"`: x towards the nose).

# Arguments

- `obj::NeuroAnalyzer.NEURO`

# Returns

Named tuple:

- `system::String`: coordinate system (`""`: unknown)
- `unit::String`: `"m"`, `"cm"`, `"mm"` or `""` (unknown)
- `guessed::Bool`: `true` if the unit was guessed (no valid stored unit, or the stored unit is unknown)
"""
function locs_coord(obj::NeuroAnalyzer.NEURO)::NamedTuple
    rec = obj.header.recording
    sys_ok, unit_ok = _locs_coord_valid(obj)
    sys = sys_ok ? string(get(rec, :coord_system, "")) : ""
    u = unit_ok ? string(get(rec, :length_unit, "")) : ""
    # unit not stored, unknown or no longer valid: guessed
    isempty(u) && return (system = sys, unit = _locs_guess_unit(obj.locs), guessed = true)
    return (system = sys, unit = u, guessed = false)
end

"""
    set_locs_coord(obj; <keyword arguments>)

Return a copy of `obj` with the coordinate system and unit of the channel locations set, see `set_locs_coord!()`.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `system::Union{Nothing, String}=nothing`: coordinate system, e.g. `"CapTrak"`, `"EEGLAB"`, `"MNI305"`, `"Other"` (BIDS names); `nothing`: unchanged (`""` if not set)
- `unit::Union{Nothing, String}=nothing`: `"m"`, `"cm"`, `"mm"` or `""` (unknown); `nothing`: unchanged, or guessed if not set (see `locs_coord()`)

# Returns

- `NeuroAnalyzer.NEURO`
"""
function set_locs_coord(
    obj::NeuroAnalyzer.NEURO;
    system::Union{Nothing, String} = nothing,
    unit::Union{Nothing, String} = nothing,
)::NeuroAnalyzer.NEURO
    obj_new = deepcopy(obj)
    set_locs_coord!(obj_new; system = system, unit = unit)

    return obj_new
end

"""
    set_locs_coord!(obj; <keyword arguments>)

Set the coordinate system and unit of the channel locations (`obj.locs`) in the header, see `locs_coord()`. The values apply to the current locations: when the locations are changed afterwards, they are no longer used.

# Arguments

- `obj::NeuroAnalyzer.NEURO`
- `system::Union{Nothing, String}=nothing`: coordinate system, e.g. `"CapTrak"`, `"EEGLAB"`, `"MNI305"`, `"Other"` (BIDS names); `nothing`: unchanged (`""` if not set)
- `unit::Union{Nothing, String}=nothing`: `"m"`, `"cm"`, `"mm"` or `""` (unknown); `nothing`: unchanged, or guessed if not set (see `locs_coord()`)

# Returns

- `Nothing`
"""
function set_locs_coord!(
    obj::NeuroAnalyzer.NEURO;
    system::Union{Nothing, String} = nothing,
    unit::Union{Nothing, String} = nothing,
)::Nothing
    DataFrames.nrow(obj.locs) > 0 || throw(ArgumentError("OBJ has no channel locations."))
    isnothing(unit) || _check_var(unit, collect(_LOCS_UNITS), "unit")
    c = locs_coord(obj)
    sys = isnothing(system) ? c.system : system
    u = isnothing(unit) ? c.unit : unit
    _set_locs_coord!(obj, sys, u)
    _history!(obj, "set_locs_coord!", (; system = sys, unit = u))
    return nothing
end

# store the coordinate system and unit of the current locations (with the locations they
# describe, to detect later changes)
function _set_locs_coord!(obj::NeuroAnalyzer.NEURO, system::AbstractString, unit::AbstractString)::Nothing
    rec = obj.header.recording
    rec[:coord_system] = string(system)
    rec[:length_unit] = string(unit)
    rec[:coord_ref] = Dict{String, NTuple{3, Float64}}(
        string(r.label) => (Float64(r.loc_x), Float64(r.loc_y), Float64(r.loc_z))
        for r in eachrow(obj.locs)
    )
    return nothing
end

# do the stored system and unit still apply to the current locations? (system, unit)
# - unit: the scale is the same (median distance from the origin within 2×)
# - system: at most half of the channels present in both have moved
function _locs_coord_valid(obj::NeuroAnalyzer.NEURO)::Tuple{Bool, Bool}
    rec = obj.header.recording
    DataFrames.nrow(obj.locs) == 0 && return (false, false)
    ref = get(rec, :coord_ref, nothing)
    # no stored locations (e.g. objects saved before, or headers created with a unit):
    # a stored unit is kept, the system is unknown
    (isnothing(ref) || isempty(ref)) && return (false, !isempty(string(get(rec, :length_unit, ""))))
    # scale
    dist(p) = sqrt(p[1]^2 + p[2]^2 + p[3]^2)
    cur = [(Float64(r.loc_x), Float64(r.loc_y), Float64(r.loc_z)) for r in eachrow(obj.locs)]
    d1 = Base.filter(x -> isfinite(x) && x > 0, dist.(cur))
    d0 = Base.filter(x -> isfinite(x) && x > 0, dist.(collect(values(ref))))
    unit_ok = if isempty(d1) || isempty(d0)
        isempty(d1) && isempty(d0)
    else
        0.5 <= median(d1) / median(d0) <= 2
    end
    # moved channels
    n = 0
    moved = 0
    for (k, r) in enumerate(eachrow(obj.locs))
        p = get(ref, string(r.label), nothing)
        isnothing(p) && continue
        n += 1
        all(isapprox.(cur[k], p; rtol = 1e-6, atol = 1e-9)) || (moved += 1)
    end
    sys_ok = n > 0 && moved <= n / 2
    return (sys_ok, unit_ok)
end

# unit of the locations guessed from the head size: median distance of the channels from
# the origin (head radius ~ 8–10 cm)
function _locs_guess_unit(locs::DataFrame)::String
    DataFrames.nrow(locs) == 0 && return ""
    d = [sqrt(r.loc_x^2 + r.loc_y^2 + r.loc_z^2) for r in eachrow(locs)]
    d = Base.filter(x -> isfinite(x) && x > 0, d)
    isempty(d) && return ""
    r = median(d)
    0.02 <= r < 0.3 && return "m"
    2.5 <= r < 25 && return "cm"
    25 <= r < 300 && return "mm"
    return ""
end

# unit name of a file (e.g. SNIRF LengthUnit "mm", "millimeter") → "m", "cm", "mm" or ""
function _locs_unit_name(u::AbstractString)::String
    v = lowercase(strip(u))
    v in ("m", "meter", "meters", "metre", "metres") && return "m"
    v in ("cm", "centimeter", "centimeters", "centimetre", "centimetres") && return "cm"
    v in ("mm", "millimeter", "millimeters", "millimetre", "millimetres") && return "mm"
    return ""
end
