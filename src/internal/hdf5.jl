# ---------------------------------------------------------------------------
# HDF5 encoding of Julia values (save_neuro(), load_neuro())
#
# values of HDF5 native types are written as plain datasets: numbers, Bool, String and
# arrays of them (HDF5 stores arrays in row-major order: a Julia array of size (a, b, c) is
# seen by other software as (c, b, a)); other values are written as datasets or groups
# with the `na_type` attribute:
#
# - "Symbol", "Char", "Date", "DateTime", "Time": string dataset
# - "Regex": string dataset (pattern) with the `flags` attribute (e.g. "i")
# - "BitArray": Bool array
# - "nothing", "missing": empty group
# - "Missings": array of Union{Missing, T}: group with the `values` (missing → zero or
#   "") and `mask` (true: missing) datasets
# - "Ragged": vector of vectors of a native type: group with the `values` (concatenated)
#   and `lengths` datasets
# - "Array": other arrays: group with one member per element ("1", "2", …, in
#   column-major order) and the `size` attribute
# - "Tuple": group with one member per element ("1", "2", …)
# - "NamedTuple": group with one member per field and the `keys` attribute (field order)
# - "Dict": group with one member per key, the `key_type` attribute ("Symbol", "String",
#   "Int64") and the `value_type` attribute (Julia type of the values, e.g. "Any",
#   "Tuple{Float64, Float64, Float64}"; types not defined in Base are read as the common
#   type of the values)
# - "DataFrame": group with one member per column and the `columns` attribute (column
#   order)
# - "Text": any other value, as its text representation (`repr()`), with the `julia_type`
#   attribute; read back as a string
#
# member names are the keys or column names with "%" and "/" escaped as "%25" and "%2F"
# ("" → "%", "." → "%2E")
# ---------------------------------------------------------------------------

# element types written natively
const _H5_NATIVE = Union{
    Bool, Int8, Int16, Int32, Int64, UInt8, UInt16, UInt32, UInt64, Float32, Float64,
    ComplexF32, ComplexF64,
}

_h5_native_eltype(T::Type)::Bool =
    isconcretetype(T) && (T <: _H5_NATIVE || T === String)

function _h5_name(s::AbstractString)::String
    s == "" && return "%"
    s == "." && return "%2E"
    return replace(string(s), "%" => "%25", "/" => "%2F")
end

function _h5_unname(s::AbstractString)::String
    s == "%" && return ""
    return replace(
        string(s),
        r"%([0-9A-Fa-f]{2})" => m -> string(Char(parse(UInt8, m[2:3]; base = 16))),
    )
end

_h5_tag!(obj, tag::String)::Nothing = (HDF5.attributes(obj)["na_type"] = tag; nothing)

function _h5_tag(obj)::String
    a = HDF5.attributes(obj)
    return haskey(a, "na_type") ? string(read(a["na_type"])) : ""
end

# zero value of a native element type (placeholder of missing values)
_h5_zero(::Type{String})::String = ""
_h5_zero(T::Type)::T = zero(T)

"""
    _h5_write(parent, name, x; <keyword arguments>)

Write the value `x` as member `name` of the HDF5 file or group `parent`.

# Arguments

- `parent::Union{HDF5.File, HDF5.Group}`
- `name::String`: member name (not escaped)
- `x`: value
- `path::String=name`: path of the value (warning messages)

# Returns

- `Nothing`
"""
function _h5_write(parent, name::AbstractString, x; path::String = string(name))::Nothing
    n = _h5_name(name)

    # native scalars and arrays
    if x isa _H5_NATIVE || x isa String
        parent[n] = x
    elseif x isa Array && _h5_native_eltype(eltype(x))
        parent[n] = x
    elseif x isa AbstractString
        parent[n] = string(x)
    elseif x isa Symbol || x isa Char || x isa Dates.Date || x isa Dates.DateTime ||
           x isa Dates.Time
        parent[n] = string(x)
        _h5_tag!(parent[n], string(nameof(typeof(x))))
    elseif x isa Regex
        parent[n] = x.pattern
        _h5_tag!(parent[n], "Regex")
        r = repr(x)
        HDF5.attributes(parent[n])["flags"] = r[(findlast('"', r) + 1):end]
    elseif x isa BitArray
        parent[n] = Array{Bool}(x)
        _h5_tag!(parent[n], "BitArray")
    elseif x === nothing || x === missing
        g = HDF5.create_group(parent, n)
        _h5_tag!(g, x === nothing ? "nothing" : "missing")
    elseif x isa AbstractArray && !(x isa Array) && _h5_native_eltype(eltype(x))
        # ranges, views, etc.
        parent[n] = collect(x)
    elseif x isa AbstractArray && (T = nonmissingtype(eltype(x))) !== eltype(x) &&
           _h5_native_eltype(T)
        g = HDF5.create_group(parent, n)
        _h5_tag!(g, "Missings")
        mask = Array{Bool}(ismissing.(x))
        g["values"] = T[ismissing(v) ? _h5_zero(T) : v for v in x]
        g["mask"] = mask
    elseif x isa AbstractVector && eltype(x) <: AbstractVector &&
           _h5_native_eltype(eltype(eltype(x)))
        T = eltype(eltype(x))
        g = HDF5.create_group(parent, n)
        _h5_tag!(g, "Ragged")
        g["values"] = isempty(x) ? T[] : T[v for e in x for v in e]
        g["lengths"] = Int64[length(e) for e in x]
    elseif x isa AbstractDict
        g = HDF5.create_group(parent, n)
        _h5_tag!(g, "Dict")
        K = keytype(x)
        kt = K === Symbol ? "Symbol" : K <: AbstractString ? "String" :
             K <: Integer ? "Int64" : ""
        if kt == ""
            # keys of other types: all keys must be of one of the supported types
            if all(k -> k isa Symbol, keys(x))
                kt = "Symbol"
            elseif all(k -> k isa AbstractString, keys(x))
                kt = "String"
            elseif all(k -> k isa Integer, keys(x))
                kt = "Int64"
            else
                throw(
                    ArgumentError(
                        "$path: dictionary keys must be symbols, strings or integers.",
                    ),
                )
            end
        end
        HDF5.attributes(g)["key_type"] = kt
        HDF5.attributes(g)["value_type"] = string(valtype(x))
        for (k, v) in x
            _h5_write(g, string(k), v; path = "$path/$k")
        end
    elseif x isa DataFrame
        g = HDF5.create_group(parent, n)
        _h5_tag!(g, "DataFrame")
        HDF5.attributes(g)["columns"] = names(x)
        for c in names(x)
            _h5_write(g, c, x[!, c]; path = "$path/$c")
        end
    elseif x isa NamedTuple
        g = HDF5.create_group(parent, n)
        _h5_tag!(g, "NamedTuple")
        HDF5.attributes(g)["keys"] = String[string(k) for k in keys(x)]
        for (k, v) in pairs(x)
            _h5_write(g, string(k), v; path = "$path/$k")
        end
    elseif x isa Tuple
        g = HDF5.create_group(parent, n)
        _h5_tag!(g, "Tuple")
        for (k, v) in enumerate(x)
            _h5_write(g, string(k), v; path = "$path/$k")
        end
    elseif x isa AbstractArray
        g = HDF5.create_group(parent, n)
        _h5_tag!(g, "Array")
        HDF5.attributes(g)["size"] = Int64[size(x)...]
        for (k, v) in enumerate(x)
            _h5_write(g, string(k), v; path = "$path/$k")
        end
    else
        _warn(
            "$path: values of type $(typeof(x)) are saved as text and loaded as strings.",
        )
        parent[n] = repr(x)
        _h5_tag!(parent[n], "Text")
        HDF5.attributes(parent[n])["julia_type"] = string(typeof(x))
    end

    return nothing
end

# narrowest common type of the values (the element type of a decoded container)
function _h5_common_type(v)::Type
    isempty(v) && return Any
    return mapreduce(typeof, typejoin, v)
end

# Julia type from its name (e.g. "Tuple{Float64, Float64, Float64}"); only types of Base and
# Core (and their parameters: types, integers, symbols) are resolved, nothing is evaluated;
# nothing: unknown type
function _h5_type(s::AbstractString)::Union{Nothing, Type}
    e = try
        Meta.parse(s)
    catch
        return nothing
    end
    function resolve(x)
        x isa Integer && return x
        x isa QuoteNode && x.value isa Symbol && return x.value
        if x isa Symbol
            isdefined(Base, x) || return nothing
            t = getfield(Base, x)
            return t isa Type ? t : nothing
        end
        if x isa Expr && x.head === :curly
            ps = Any[resolve(a) for a in x.args]
            any(isnothing, ps) && return nothing
            ps[1] isa Type || return nothing
            return try
                ps[1]{ps[2:end]...}
            catch
                nothing
            end
        end
        return nothing
    end
    t = resolve(e)
    return t isa Type ? t : nothing
end

# members of a group ordered by their number ("1", "2", …)
_h5_numbered(g)::Vector{String} = sort(collect(keys(g)); by = k -> parse(Int64, k))

"""
    _h5_read(obj)

Read a value written by `_h5_write()`.

# Arguments

- `obj::Union{HDF5.Dataset, HDF5.Group}`

# Returns

- the value
"""
function _h5_read(obj)
    tag = _h5_tag(obj)

    if obj isa HDF5.Dataset
        v = read(obj)
        tag == "" && return v
        tag == "Symbol" && return Symbol(v)
        tag == "Char" && return only(v)
        tag == "Date" && return Dates.Date(v)
        tag == "DateTime" && return Dates.DateTime(v)
        tag == "Time" && return Dates.Time(v)
        tag == "BitArray" && return BitArray(v)
        tag == "Regex" && return Regex(v, string(HDF5.read_attribute(obj, "flags")))
        # "Text" and unknown tags: the stored value
        return v
    end

    tag == "nothing" && return nothing
    tag == "missing" && return missing
    if tag == "Missings"
        v = read(obj["values"])
        mask = read(obj["mask"])
        out = Array{Union{Missing, eltype(v)}}(v)
        out[mask] .= missing
        return out
    elseif tag == "Ragged"
        v = read(obj["values"])
        lens = read(obj["lengths"])
        ends = cumsum(lens)
        return typeof(v)[v[(e - l + 1):e] for (l, e) in zip(lens, ends)]
    elseif tag == "Dict"
        kt = string(HDF5.read_attribute(obj, "key_type"))
        K = kt == "Symbol" ? Symbol : kt == "Int64" ? Int64 : String
        conv = kt == "Symbol" ? Symbol : kt == "Int64" ? (s -> parse(Int64, s)) : identity
        ks = K[conv(_h5_unname(k)) for k in keys(obj)]
        vs = Any[_h5_read(obj[k]) for k in keys(obj)]
        V = haskey(HDF5.attributes(obj), "value_type") ?
            _h5_type(string(HDF5.read_attribute(obj, "value_type"))) : nothing
        (isnothing(V) || !all(v -> v isa V, vs)) && (V = _h5_common_type(vs))
        return Dict{K, V}(zip(ks, vs))
    elseif tag == "DataFrame"
        cols = HDF5.read_attribute(obj, "columns")
        cols isa AbstractString && (cols = [cols])
        return DataFrame(
            [string(c) => _h5_read(obj[_h5_name(c)]) for c in cols]; copycols = false,
        )
    elseif tag == "NamedTuple"
        ks = HDF5.read_attribute(obj, "keys")
        ks isa AbstractString && (ks = [ks])
        return NamedTuple{Tuple(Symbol.(ks))}(
            Tuple(_h5_read(obj[_h5_name(k)]) for k in ks),
        )
    elseif tag == "Tuple"
        return Tuple(_h5_read(obj[k]) for k in _h5_numbered(obj))
    elseif tag == "Array"
        sz = Tuple(HDF5.read_attribute(obj, "size"))
        vs = Any[_h5_read(obj[k]) for k in _h5_numbered(obj)]
        out = Array{_h5_common_type(vs)}(undef, sz...)
        for k in eachindex(vs)
            out[k] = vs[k]
        end
        return out
    end

    # a group without a tag: dictionary of its members
    return Dict{String, Any}(_h5_unname(k) => _h5_read(obj[k]) for k in keys(obj))
end
