# ---------------------------------------------------------------------------
# FIFF (Elekta/MEGIN Neuromag Functional Image File Format) low-level writer
#
# All FIFF data are big-endian. Tag header = 4 × Int32: kind, type, size, next.
# Tags are written sequentially (next = 0); the file ends with a NOP tag whose
# next = -1. No tag directory is written: readers (MNE-Python, NeuroAnalyzer)
# rebuild it by scanning the file.
#
# Reference implementation used for cross-checking: MNE-Python (mne/_fiff/write.py).
# ---------------------------------------------------------------------------

# FIFF data types (FIFFT_*)
const _FIFFT_VOID = 0
const _FIFFT_INT = 3
const _FIFFT_FLOAT = 4
const _FIFFT_DOUBLE = 5
const _FIFFT_STRING = 10
const _FIFFT_CH_INFO_STRUCT = 30
const _FIFFT_ID_STRUCT = 31
const _FIFFT_DIG_POINT_STRUCT = 33
const _FIFFT_COORD_TRANS_STRUCT = 35

# FIFF file format version 1.3
const _FIFF_VERSION = Int32((1 << 16) | 3)

_fiff_be_bytes(x::AbstractVector{T}) where {T} = reinterpret(UInt8, hton.(x))

# tag header + payload
function _fiff_write_tag(io::IO, kind::Integer, type::Integer, payload::AbstractVector{UInt8}; next::Integer = 0)::Nothing
    write(io, hton(Int32(kind)), hton(Int32(type)), hton(Int32(length(payload))), hton(Int32(next)))
    write(io, payload)
    return nothing
end

_fiff_write_int(io::IO, kind::Integer, v) =
    _fiff_write_tag(io, kind, _FIFFT_INT, _fiff_be_bytes(Int32.(vcat(v))))
_fiff_write_float(io::IO, kind::Integer, v) =
    _fiff_write_tag(io, kind, _FIFFT_FLOAT, _fiff_be_bytes(Float32.(vcat(v))))
_fiff_write_double(io::IO, kind::Integer, v) =
    _fiff_write_tag(io, kind, _FIFFT_DOUBLE, _fiff_be_bytes(Float64.(vcat(v))))
_fiff_write_string(io::IO, kind::Integer, s::AbstractString) =
    _fiff_write_tag(io, kind, _FIFFT_STRING, Vector{UInt8}(codeunits(String(s))))

# colon-separated name list; literal colons are escaped as {COLON} (MNE-Python convention)
_fiff_write_name_list(io::IO, kind::Integer, names::AbstractVector{<:AbstractString}) =
    _fiff_write_string(io, kind, join([replace(String(n), ":" => "{COLON}") for n in names], ":"))

_fiff_start_block(io::IO, kind::Integer) =
    _fiff_write_int(io, _find_fiff_tag("block_start"), kind)
_fiff_end_block(io::IO, kind::Integer) =
    _fiff_write_int(io, _find_fiff_tag("block_end"), kind)

# id_t: version, machine id (2 × Int32), time (seconds, microseconds)
function _fiff_write_id(io::IO, kind::Integer, machid::NTuple{2, Int32}, t::NTuple{2, Int32})::Nothing
    v = Int32[_FIFF_VERSION, machid[1], machid[2], t[1], t[2]]
    return _fiff_write_tag(io, kind, _FIFFT_ID_STRUCT, _fiff_be_bytes(v))
end

# dense float matrix: data in row-major (C) order, then dims (reversed), then ndim
function _fiff_write_float_matrix(io::IO, kind::Integer, m::AbstractMatrix{<:Real})::Nothing
    payload = vcat(
        _fiff_be_bytes(Float32.(vec(permutedims(m)))),
        _fiff_be_bytes(Int32[size(m, 2), size(m, 1), 2]),
    )
    return _fiff_write_tag(io, kind, _FIFFT_FLOAT | _FIFFT_MATRIX, payload)
end

# ch_info_t (96 B): scan_no, log_no, kind, range, cal, coil_type, loc[12], unit, unit_mul, name[16]
function _fiff_write_ch_info(
    io::IO;
    scan_no::Integer,
    log_no::Integer,
    kind::Integer,
    range::Real,
    cal::Real,
    coil_type::Integer,
    loc::AbstractVector{<:Real},
    unit::Integer,
    unit_mul::Integer,
    ch_name::AbstractString,
)::Nothing
    length(loc) == 12 || throw(ArgumentError("loc must have 12 elements."))
    name = Vector{UInt8}(codeunits(String(ch_name)))
    length(name) <= 15 || throw(ArgumentError("ch_info name must have at most 15 bytes."))
    payload = vcat(
        _fiff_be_bytes(Int32[scan_no, log_no, kind]),
        _fiff_be_bytes(Float32[range, cal]),
        _fiff_be_bytes(Int32[coil_type]),
        _fiff_be_bytes(Float32.(loc)),
        _fiff_be_bytes(Int32[unit, unit_mul]),
        name,
        zeros(UInt8, 16 - length(name)),
    )
    return _fiff_write_tag(io, _find_fiff_tag("ch_info"), _FIFFT_CH_INFO_STRUCT, payload)
end

# dig_point_t (20 B): kind, ident, r[3]
function _fiff_write_dig_point(io::IO, kind::Integer, ident::Integer, r::AbstractVector{<:Real})::Nothing
    payload = vcat(_fiff_be_bytes(Int32[kind, ident]), _fiff_be_bytes(Float32.(r[1:3])))
    return _fiff_write_tag(io, _find_fiff_tag("dig_point"), _FIFFT_DIG_POINT_STRUCT, payload)
end

# coord_trans_t (104 B): from, to, rot (row-major), move, invrot (row-major), invmove
function _fiff_write_coord_trans(io::IO, from::Integer, to::Integer, t::AbstractMatrix{<:Real})::Nothing
    size(t) == (4, 4) || throw(ArgumentError("Coordinate transform must be a 4 × 4 matrix."))
    ti = inv(Matrix{Float64}(t))
    payload = vcat(
        _fiff_be_bytes(Int32[from, to]),
        _fiff_be_bytes(Float32.(vcat(vec(permutedims(t[1:3, 1:3])), t[1:3, 4]))),
        _fiff_be_bytes(Float32.(vcat(vec(permutedims(ti[1:3, 1:3])), ti[1:3, 4]))),
    )
    return _fiff_write_tag(io, _find_fiff_tag("coord_trans"), _FIFFT_COORD_TRANS_STRUCT, payload)
end

# file header: file_id, dir_pointer = -1 (no directory), free_list = -1
function _fiff_start_file(io::IO, machid::NTuple{2, Int32}, t::NTuple{2, Int32})::Nothing
    _fiff_write_id(io, _find_fiff_tag("file_id"), machid, t)
    _fiff_write_int(io, _find_fiff_tag("dir_pointer"), -1)
    _fiff_write_int(io, _find_fiff_tag("free_list"), -1)
    return nothing
end

# last tag of the file: NOP with next = -1
_fiff_end_file(io::IO) =
    _fiff_write_tag(io, _find_fiff_tag("nop"), _FIFFT_VOID, UInt8[]; next = -1)

# unique channel names of at most 15 bytes (ch_info_t limit), as in MNE-Python:
# names are cut to 15 bytes (whole UTF-8 characters), later duplicates get a "-N" suffix
function _fiff_short_names(names::AbstractVector{<:AbstractString})::Vector{String}
    function cut(s::String, nb::Int)
        r = s
        while ncodeunits(r) > nb
            r = chop(r)
        end
        return String(r)
    end
    short = [cut(String(n), 15) for n in names]
    seen = Set{String}()
    for i in eachindex(short)
        if short[i] in seen
            k = 1
            while true
                sfx = "-$k"
                cand = cut(short[i], 15 - ncodeunits(sfx)) * sfx
                if !(cand in seen) && !(cand in short)
                    short[i] = cand
                    break
                end
                k += 1
            end
        end
        push!(seen, short[i])
    end
    return short
end
