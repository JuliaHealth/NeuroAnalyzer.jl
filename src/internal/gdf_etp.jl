# ---------------------------------------------------------------------------
# GDF event type code lookup table.
#
# Source: BioSig eventcodes.txt (biosig 3.9.5, table © 2004-2015 Alois Schlögl),
# which is part of the GDF v1.x / v2.x specification. The file is shipped
# unmodified next to this source file, as required by its (non-GPL) license -
# see the file header. It is parsed once, at precompilation; include_dependency
# makes the package recompile if the file changes.
#
# Format: "0x####<TAB>description"; lines starting with "#" are omitted from the
# official table. Commented code lines (e.g. "#0x8502 ecg:P-wave-end") are kept
# in a separate table of reserved names, used only for codes not in the official
# table.
# ---------------------------------------------------------------------------
const _GDF_EVENTCODES_FILE = joinpath(@__DIR__, "eventcodes.txt")
include_dependency(_GDF_EVENTCODES_FILE)

"""
    _read_gdf_eventcodes(file)

Parse a BioSig `eventcodes.txt` file.

# Returns

- `Dict{UInt16, String}`: official event codes (uncommented lines)
- `Dict{UInt16, String}`: reserved/proposed event codes (commented-out code lines)
"""
function _read_gdf_eventcodes(file::String)::Tuple{Dict{UInt16, String}, Dict{UInt16, String}}
    official = Dict{UInt16, String}()
    reserved = Dict{UInt16, String}()
    for line in eachline(file)
        m = match(r"^(#?)0x([0-9A-Fa-f]{4})\s+(.*?)\s*$", line)
        isnothing(m) && continue
        # display form only: some descriptions are column-aligned with tabs
        desc = replace(m[3], r"\s+" => " ")
        # range annotations such as "#0x031b - 0x037f reserved ..." are comments, not codes
        (isempty(desc) || startswith(desc, "- 0x")) && continue
        code = parse(UInt16, m[2]; base = 16)
        tbl = isempty(m[1]) ? official : reserved
        haskey(tbl, code) || (tbl[code] = desc)
    end
    return official, reserved
end

const _GDF_ETP_TABLE, _GDF_ETP_RESERVED = _read_gdf_eventcodes(_GDF_EVENTCODES_FILE)

# description of a single code: official table first, then reserved names
_gdf_etp_desc(code::UInt16)::Union{String, Nothing} =
    get(_GDF_ETP_TABLE, code, get(_GDF_ETP_RESERVED, code, nothing))

"""
    _gdf_etp(etp)

Decode a 2-byte GDF event type code into a human-readable description.

The code is a little-endian `UInt16`, as stored in the GDF event table, and is looked up in BioSig `eventcodes.txt`. If bit `0x8000` is set the event marks the end (offset) of the event `code - 0x8000`; unless the end code has its own name in the table (e.g. `0x8502` → "ecg:P-wave-end, MDC_ECG_TIME_END_P"), `" (end)"` is appended to the description of that event. Returns an "unknown code" string for unrecognized codes.

# Arguments

- `etp::Vector{UInt8}`: 2 bytes of the event type as read from the file; `etp[1]` is the low byte (specific event), `etp[2]` is the high byte (category)
"""
function _gdf_etp(etp::Vector{UInt8})::String
    length(etp) >= 2 || throw(ArgumentError("etp must contain at least 2 bytes."))

    # GDF is little-endian
    code = UInt16(etp[1]) | (UInt16(etp[2]) << 8)
    desc = _gdf_etp_desc(code)
    !isnothing(desc) && return desc

    # bit 15 flags the end of the event given by the remaining bits
    if (code & 0x8000) != 0
        desc = _gdf_etp_desc(code & 0x7fff)
        !isnothing(desc) && return desc * " (end)"
    end

    return "unknown code (0x$(string(code; base = 16, pad = 4)))"
end
