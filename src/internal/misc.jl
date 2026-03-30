"""
    _pl(x)

Return `"s"` if `x` represents a plural quantity, otherwise `""`.

Accepts a range/vector (plural when length > 1).
"""
_pl(x::Union{AbstractRange, AbstractVector})::String = length(x) > 1 ? "s" : ""

"""
    _pl(x)

Return `"s"` if `x` represents a plural quantity, otherwise `""`.

Accepts a scalar Real (plural when > 1).
"""
_pl(x::Real)::String = x > 1 ? "s" : ""

"""
    _get_range(s)

Return the peak-to-peak range of `s` (i.e. `rng(s)`), rounded to 0 decimal places.
"""
_get_range(s::Union{AbstractVector, AbstractArray})::Float64 = round(rng(s); digits = 0)

"""
    _c(n)

Return `collect(1:n)` as a `Vector{Int64}`.
"""
_c(n::Integer)::Vector{Int64} = collect(1:n)

"""
    _tuple_max(t)

Return a symmetric tuple `(-m, m)` where `m = max(|t[1]|, |t[2]|)`.

Useful for building balanced y-axis limits.
"""
function _tuple_max(t::Tuple{Real, Real})::Tuple{Real, Real}
    m = max(abs(t[1]), abs(t[2]))
    return (-m, m)
end

"""
    _n2v(s)

Wrap a scalar number in a single-element vector, or return the vector unchanged.
"""
_n2v(s::Union{<:Number, Vector{<:Number}})::Vector{<:Number} = s isa Number ? [s] : s

"""
    _v2s(v)

Concatenate all strings in `v` into a single `String`.
"""
_v2s(v::Vector{String})::String = join(v)

"""
    _v2s(x)

Convert a numeric vector to a comma-separated `String` (e.g. `[1, 2, 3]` → `"1, 2, 3"`).
"""
_v2s(x::Vector{<:Number})::String = join(string.(x), ", ")

"""
    _copy_lt2ut(m)

Mirror the lower triangle of a symmetric matrix (or batch of matrices) to the upper triangle, producing a fully symmetric result.

For a 2-D matrix returns a new matrix `m + m' - diag(m)`.

For a 3-D array operates in-place on each slice along the third dimension.

# Notes

The 2-D method returns a new array while the 3-D method mutates the input. This asymmetry is a known limitation; use `copy(m)` before calling if in-place behaviour is undesirable for the 3-D case.
"""
function _copy_lt2ut(m::AbstractArray)::AbstractArray
    if ndims(m) == 2
        return m + m' - diagm(diag(m))
    else
        Threads.@threads :static for ep_idx in axes(m, 3)
            @inbounds m[:, :, ep_idx] =
                m[:, :, ep_idx] + m[:, :, ep_idx]' - diagm(diag(m[:, :, ep_idx]))
        end
        return m
    end
end

"""
    _tlength(t)

Return the number of integer steps from `t[1]` to `t[2]` inclusive.
"""
_tlength(t::Tuple{Real, Real})::Int64 = length(t[1]:1:t[2])

"""
    _s2i(s)

Parse a string representation of an integer or integer range/list into `Int64` or `Vector{Int64}`.

Supported formats:
- `"3"` → `3`
- `"1:5"` → `[1, 2, 3, 4, 5]`
- `"1,3,7"` or `"[1,3,7]"` → `[1, 3, 7]`

Throws if the string does not match any recognized format.
"""
function _s2i(s::String)::Union{Int64, Vector{Int64}}
    s = replace(s, " " => "", "[" => "", "]" => "")
    if occursin(":", s)
        parts = split(s, ":")
        return collect(parse(Int64, parts[1]):parse(Int64, parts[2]))
    elseif occursin(",", s)
        return parse.(Int64, split(s, ","))
    elseif _check_sint(s)
        return parse(Int64, s)
    else
        throw(ArgumentError("Cannot parse \"$s\" as an integer or integer range/list."))
    end
end

"""
    _i2s(s)

Convert an integer, vector of integers, or range to a comma-separated `String`.
"""
function _i2s(s::Union{Int64, Vector{Int64}, AbstractRange})::String
    s_str = string(collect(s))
    return replace(s_str, "[" => "", "]" => "")
end

"""
    _s2tf(s)

Parse a string of the form `"(f1, f2)"` into a `Tuple{Float64, Float64}`.
"""
function _s2tf(s::String)::Tuple{Float64, Float64}
    s = replace(s, " " => "", "(" => "", ")" => "")
    parts = split(s, ",")
    return (parse(Float64, parts[1]), parse(Float64, parts[2]))
end

"""
    _s2ti(s)

Parse a string of the form `"(i1, i2)"` into a `Tuple{Int64, Int64}`.
"""
function _s2ti(s::String)::Tuple{Int64, Int64}
    s = replace(s, " " => "", "(" => "", ")" => "")
    parts = split(s, ",")
    return (parse(Int64, parts[1]), parse(Int64, parts[2]))
end

"""
    _detect_montage(clabels, ch_type, data_type)

Infer the EEG montage type from channel labels.

# Arguments

- `clabels::Vector{String}`: all channel labels
- `ch_type::Vector{String}`: corresponding channel type strings
- `data_type::String`: the data type to inspect (e.g. `"eeg"`)

# Returns

- `String`: one of `"common (A)"`, `"common (M)"`, `"common"`, `"bipolar"`, or `"physical"`
"""
function _detect_montage(
    clabels::Vector{String},
    ch_type::Vector{String},
    data_type::String,
)::String
    target = clabels[ch_type .== data_type]

    # check for bipolar / common reference pattern: "label-ref"
    m = match.(r"(.+)\-(.+)", lowercase.(target))
    if length(findall(!isnothing, m)) == length(target)
        refs = [m[idx].captures[2] for idx in eachindex(m)]
        if length(unique(refs)) == 1
            occursin("a", lowercase(refs[1])) && return "common (A)"
            occursin("m", lowercase(refs[1])) && return "common (M)"
            return "common"
        else
            return "bipolar"
        end
    end

    # check for paired label pattern: "label1num1label2num2"
    m = match.(r"([a-z]+)([0-9]+[0-9]?)([a-z]+)([0-9]+)", lowercase.(target))
    if length(findall(!isnothing, m)) == length(target)
        refs = [m[idx].captures[3] for idx in eachindex(m)]
        if length(unique(refs)) == 1
            occursin("a", lowercase(refs[1])) && return "common (A)"
            occursin("m", lowercase(refs[1])) && return "common (M)"
            return "common"
        else
            return "bipolar"
        end
    end

    return "physical"
end

"""
    _fread(fid, n, t)
 
Read `n` values of type `t` from an open IO stream `fid`, returning the result as `Int64`, `Float64`, or `Vector{Int64}`.
 
# Type symbols
 
| Symbol   | Meaning                          |
|----------|----------------------------------|
| `:s`     | string of `n` bytes (→ `Vector`) |
| `:c`     | single byte character            |
| `:l`     | 32-bit signed integer            |
| `:ul`    | 32-bit unsigned integer          |
| `:ui8`   | 8-bit unsigned integer           |
| `:ui16`  | 16-bit unsigned integer          |
| `:ui32`  | 32-bit unsigned integer          |
| `:ui64`  | 64-bit unsigned integer          |
| `:i`     | 32-bit signed integer            |
| `:i8`    | 8-bit signed integer             |
| `:i16`   | 16-bit signed integer            |
| `:i32`   | 32-bit signed integer            |
| `:i64`   | 64-bit signed integer            |
| `:f16`   | 16-bit float                     |
| `:f32`   | 32-bit float                     |
| `:f64`   | 64-bit float                     |
"""
function _fread(fid, n::Int64, t::Symbol)::Union{Int64, Float64, Vector{Int64}}
    (n > 1 && t === :c) && (t = :s)

    # compute byte count
    nbytes = if t in (:s, :c, :ui8, :i8)
        n
    elseif t in (:ui16, :i16, :f16)
        n * 2
    elseif t in (:l, :ul, :ui32, :i, :i32, :f32)
        n * 4
    elseif t in (:ui64, :i64, :f64)
        n * 8
    else
        throw(ArgumentError("Unknown type symbol :$t"))
    end

    header = zeros(UInt8, nbytes)
    readbytes!(fid, header, nbytes)

    # decode and return
    t === :s && return Int64.(map(ltoh, reinterpret(UInt8, header)))
    t === :c && return Int64(map(ltoh, reinterpret(UInt8, header))[1])
    t === :l && return Int64(map(ltoh, reinterpret(Int32, header))[1])
    t === :ul && return Int64(map(ltoh, reinterpret(UInt32, header))[1])
    t === :ui8 && return Int64(map(ltoh, reinterpret(UInt8, header))[1])
    t === :ui16 && return Int64(map(ltoh, reinterpret(UInt16, header))[1])
    t === :ui32 && return Int64(map(ltoh, reinterpret(UInt32, header))[1])
    t === :ui64 && return Int64(map(ltoh, reinterpret(UInt64, header))[1])
    t === :i && return Int64(map(ltoh, reinterpret(Int32, header))[1])
    t === :i8 && return Int64(map(ltoh, reinterpret(Int8, header))[1])
    t === :i16 && return Int64(map(ltoh, reinterpret(Int16, header))[1])
    t === :i32 && return Int64(map(ltoh, reinterpret(Int32, header))[1])
    t === :i64 && return Int64(map(ltoh, reinterpret(Int64, header))[1])
    t === :f16 && return Float64(map(ltoh, reinterpret(Float16, header))[1])
    t === :f32 && return Float64(map(ltoh, reinterpret(Float32, header))[1])
    return Float64(map(ltoh, reinterpret(Float64, header))[1])
end

"""
    _vint2str(x)

Convert a vector of integer code points to a `String`, stripping null characters.
"""
function _vint2str(x::Vector{Int64})::String
    return replace(strip(String(Char.(x))), "\0" => "")
end

"""
    _swap(x, y)

Return `(y, x)` — swap two values.
"""
_swap(x, y) = (y, x)

"""
    _veqlen(s1, s2)

Pad the shorter of two vectors with trailing zeros so both have equal length.

Returns `(s1, s2)` with the shorter one zero-padded.
"""
function _veqlen(
    s1::AbstractVector,
    s2::AbstractVector,
)::Tuple{AbstractVector, AbstractVector}
    if length(s1) > length(s2)
        return s1, pad0(s2, length(s1) - length(s2))
    elseif length(s2) > length(s1)
        return pad0(s1, length(s2) - length(s1)), s2
    else
        return s1, s2
    end
end

"""
    _fmem()

Return the amount of free system memory in megabytes.
"""
_fmem()::Float64 = Sys.free_memory() / 2^20
