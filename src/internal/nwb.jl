# ------------------------------------------------------------------------------
# NWB (Neurodata Without Borders 2.x) on HDF5: reading and writing of the NWB building blocks,
# used by import_nwb() and export_nwb(); shares the HDF5 helpers of src/internal/hdf5.jl
#
# - objects: groups and datasets with the attributes `neurodata_type`, `namespace` (e.g.
#   "core", "hdmf-common", "ndx-events") and `object_id` (UUID)
# - DynamicTable: group with the `colnames` and `description` attributes, an `id` dataset
#   and one dataset per column (VectorData); ragged columns: the values and a `<col>_index`
#   dataset (VectorIndex: end position of each row, 0-based exclusive, `target` attribute:
#   reference to the values); DynamicTableRegion columns: row numbers (0-based) of another
#   table (`table` attribute: reference to the table)
# - TimeSeries: group with `data` (attributes `unit`, `conversion`, `offset`, `resolution`),
#   `starting_time` (attribute `rate`) or `timestamps`, optional `electrodes`
#   (DynamicTableRegion) and `channel_conversion`
# - Bool values: 8-bit enumeration FALSE = 0, TRUE = 1 (as written by h5py/pynwb)
# - array dimensions: NWB stores (time, channels); HDF5 row-major order is read by HDF5.jl in
#   reverse, i.e. as (channels, time)
# ------------------------------------------------------------------------------

# compression filters of HDF5 (decoded by HDF5.jl without plugins): deflate (gzip), shuffle,
# Fletcher32 checksum, scale-offset, N-bit, SZIP
const _NWB_FILTERS_OK = (1, 2, 3, 4, 5, 6)

# value of the attribute `name` of an HDF5 object (`default` if absent)
function _nwb_attr(obj, name::AbstractString, default = nothing)
    a = HDF5.attributes(obj)
    haskey(a, name) || return default
    return read(a[name])
end

# neurodata type and namespace of an HDF5 object ("" if none)
_nwb_type(obj)::Tuple{String, String} = (
    string(something(_nwb_attr(obj, "neurodata_type"), "")),
    string(something(_nwb_attr(obj, "namespace"), "")),
)

# identity of an HDF5 object (the same object reached through links has the same identity)
function _nwb_objid(obj)
    i = HDF5.API.h5o_get_info(obj)
    return (i.fileno, i.token)
end

"""
    _nwb_objects(f; <keyword arguments>)

List the NWB objects (groups and datasets with a `neurodata_type` attribute) of an NWB file. Each object is listed once, under the first path at which it is found (objects reached through links, e.g. the device of an electrode group, are not listed again); the cached specifications (`/specifications`) are not searched.

# Arguments

- `f::HDF5.File`
- `types::Union{Nothing, Vector{String}}=nothing`: neurodata types to list (e.g. `["ElectricalSeries"]`); `nothing`: all

# Returns

- `DataFrame`: columns `path`, `type` (neurodata type), `namespace`
"""
function _nwb_objects(f::HDF5.File; types::Union{Nothing, Vector{String}} = nothing)::DataFrame
    out = DataFrame(:path => String[], :type => String[], :namespace => String[])
    seen = Set{Any}()
    function walk(g, prefix::String)
        for k in keys(g)
            p = prefix * "/" * k
            p == "/specifications" && continue
            o = try
                g[k]
            catch
                # broken or external link
                continue
            end
            try
                id = _nwb_objid(o)
                id in seen && continue
                push!(seen, id)
                t, ns = _nwb_type(o)
                !isempty(t) && (isnothing(types) || t in types) && push!(out, (p, t, ns))
                o isa HDF5.Group && walk(o, p)
            finally
                close(o)
            end
        end
        return nothing
    end
    walk(f, "")
    return out
end

# path of the object an HDF5 object reference points to
_nwb_ref_path(f::HDF5.File, r::HDF5.Reference)::String = (o = f[r]; p = HDF5.name(o); close(o); p)

# is the dataset (or attribute) type an h5py Bool enumeration (FALSE = 0, TRUE = 1)?
function _nwb_isbool(obj)::Bool
    dt = HDF5.datatype(obj)
    try
        HDF5.API.h5t_get_class(dt) == HDF5.API.H5T_ENUM || return false
        HDF5.API.h5t_get_nmembers(dt) == 2 || return false
        names = [HDF5.API.h5t_get_member_name(dt, i) for i = 0:1]
        return uppercase.(names) == ["FALSE", "TRUE"]
    finally
        close(dt)
    end
end

"""
    _nwb_read(f, ds)

Read an NWB dataset: Bool enumerations as `Bool`, object references as the paths of the referenced objects (also in compound values, e.g. the `timeseries` column of trials), other values as stored. Datasets compressed with filters that HDF5.jl cannot decode (e.g. Blosc, Zstandard, LZF) raise an error naming the filter.

# Arguments

- `f::HDF5.File`
- `ds::HDF5.Dataset`

# Returns

- the values
"""
function _nwb_read(f::HDF5.File, ds::HDF5.Dataset)
    _nwb_check_filters(ds)
    _nwb_isbool(ds) && return Bool.(read(ds) .!= 0)
    v = read(ds)
    conv(x) = x isa HDF5.Reference ? _nwb_ref_path(f, x) :
              x isa NamedTuple ? map(conv, x) : x
    if v isa AbstractArray && !isempty(v) && (eltype(v) <: HDF5.Reference || eltype(v) <: NamedTuple)
        return map(conv, v)
    end
    return conv(v)
end

# error for compression filters that HDF5.jl cannot decode
function _nwb_check_filters(ds::HDF5.Dataset)::Nothing
    p = HDF5.get_create_properties(ds)
    try
        n = HDF5.API.h5p_get_nfilters(p)
        for i = 0:(n - 1)
            id = HDF5.API.h5p_get_filter(p, i, C_NULL, C_NULL, C_NULL, 0, C_NULL, C_NULL)
            Int(id) in _NWB_FILTERS_OK && continue
            name = Int(id) == 32001 ? "Blosc" : Int(id) == 32015 ? "Zstandard" :
                   Int(id) == 32000 ? "LZF" : Int(id) == 32004 ? "LZ4" :
                   Int(id) == 32008 ? "Bitshuffle" : "filter $id"
            throw(
                ArgumentError(
                    "$(HDF5.name(ds)) is compressed with $name, which cannot be read (supported: gzip/deflate, shuffle, Fletcher32, scale-offset, N-bit, SZIP).",
                ),
            )
        end
    finally
        close(p)
    end
    return nothing
end

"""
    _nwb_table(f, g)

Read an NWB DynamicTable (e.g. the electrodes table, trials, intervals, units) into a DataFrame, in the order of the `colnames` attribute, after the `id` column.

- ragged columns (with a `<column>_index` VectorIndex, also doubly indexed): a vector per row
- DynamicTableRegion columns: row numbers (1-based) of the target table, whose path is stored as the column metadata `table`
- references (e.g. the electrode `group` column): paths of the referenced objects
- Bool enumerations: `Bool`

# Arguments

- `f::HDF5.File`
- `g::HDF5.Group`: the table

# Returns

- `DataFrame`
"""
function _nwb_table(f::HDF5.File, g::HDF5.Group)::DataFrame
    cols = _nwb_attr(g, "colnames", String[])
    cols isa AbstractString && (cols = [cols])
    df = DataFrame()
    haskey(g, "id") && (df[!, :id] = Int64.(read(g["id"])))
    for c in string.(cols)
        haskey(g, c) || continue
        ds = g[c]
        v = _nwb_read(f, ds)
        # DynamicTableRegion: 0-based row numbers of another table
        tbl = _nwb_attr(ds, "table")
        if tbl isa HDF5.Reference
            v = v isa AbstractArray ? Int64.(v) .+ 1 : Int64(v) + 1
        end
        # ragged columns: one or two VectorIndex levels (`c_index`, `c_index_index`)
        idx = c * "_index"
        while haskey(g, idx)
            ends = Int64.(read(g[idx]))
            starts = vcat(0, ends[1:(end - 1)])
            v = [v[(s + 1):e] for (s, e) in zip(starts, ends)]
            idx *= "_index"
        end
        df[!, c] = v isa AbstractArray ? v : fill(v, max(DataFrames.nrow(df), 1))
        tbl isa HDF5.Reference && DataFrames.colmetadata!(df, c, "table", _nwb_ref_path(f, tbl); style = :note)
        close(ds)
    end
    return df
end

"""
    _nwb_series(f, g; <keyword arguments>)

Read an NWB TimeSeries (e.g. ElectricalSeries, TimeSeries of other signals).

# Arguments

- `f::HDF5.File`
- `g::HDF5.Group`: the series
- `scaled::Bool=true`: return the data in the unit of the series (`data × channel_conversion × conversion + offset`, `Float64`); `false`: as stored

# Returns

Named tuple:

- `type::String`: neurodata type
- `data`: values, `(channels, time points)` for 2-D series (NWB stores time × channels), a vector for 1-D series
- `unit::String`, `conversion::Float64`, `offset::Float64`, `resolution::Float64`
- `channel_conversion::Union{Nothing, Vector{Float64}}`
- `rate::Union{Nothing, Float64}`, `starting_time::Float64`
- `timestamps::Union{Nothing, Vector{Float64}}`: time of each sample (`nothing`: regular sampling, see `rate`)
- `electrodes::Union{Nothing, Vector{Int64}}`: rows (1-based) of the electrodes table
- `electrodes_table::String`: path of the electrodes table (`""`: none)
- `description::String`, `comments::String`, `filtering::String`
"""
function _nwb_series(f::HDF5.File, g::HDF5.Group; scaled::Bool = true)::NamedTuple
    haskey(g, "data") || throw(ArgumentError("$(HDF5.name(g)) has no data."))
    ds = g["data"]
    data = _nwb_read(f, ds)
    unit = string(something(_nwb_attr(ds, "unit"), ""))
    conversion = Float64(something(_nwb_attr(ds, "conversion"), 1.0))
    offset = Float64(something(_nwb_attr(ds, "offset"), 0.0))
    resolution = Float64(something(_nwb_attr(ds, "resolution"), -1.0))
    close(ds)
    cc = haskey(g, "channel_conversion") ? Float64.(read(g["channel_conversion"])) : nothing
    if scaled && eltype(data) <: Real
        x = Float64.(data)
        if !isnothing(cc) && ndims(x) == 2
            size(x, 1) == length(cc) ||
                throw(ArgumentError("$(HDF5.name(g)): channel_conversion does not match the channels."))
            x .*= cc
        end
        data = x .* conversion .+ offset
    end
    rate = nothing
    t0 = 0.0
    ts = nothing
    if haskey(g, "starting_time")
        st = g["starting_time"]
        t0 = Float64(read(st))
        r = _nwb_attr(st, "rate")
        isnothing(r) || (rate = Float64(r))
        close(st)
    end
    haskey(g, "timestamps") && (ts = Float64.(read(g["timestamps"])); t0 = isempty(ts) ? 0.0 : ts[1])
    el = nothing
    el_tbl = ""
    if haskey(g, "electrodes")
        e = g["electrodes"]
        el = Int64.(read(e)) .+ 1
        tbl = _nwb_attr(e, "table")
        tbl isa HDF5.Reference && (el_tbl = _nwb_ref_path(f, tbl))
        close(e)
    end
    txt(k) = haskey(g, k) ? string(read(g[k])) : string(something(_nwb_attr(g, k), ""))
    return (
        type = _nwb_type(g)[1],
        data = data,
        unit = unit,
        conversion = conversion,
        offset = offset,
        resolution = resolution,
        channel_conversion = cc,
        rate = rate,
        starting_time = t0,
        timestamps = ts,
        electrodes = el,
        electrodes_table = el_tbl,
        description = txt("description"),
        comments = txt("comments"),
        filtering = txt("filtering"),
    )
end

# ------------------------------------------------------------------------------
# writing
# ------------------------------------------------------------------------------

# random UUID (version 4), the `object_id` of NWB objects
function _nwb_uuid()::String
    b = rand(UInt8, 16)
    b[7] = (b[7] & 0x0f) | 0x40
    b[9] = (b[9] & 0x3f) | 0x80
    h = bytes2hex(b)
    return join([h[1:8], h[9:12], h[13:16], h[17:20], h[21:32]], "-")
end

"""
    _nwb_attrs!(obj, neurodata_type; <keyword arguments>)

Mark an HDF5 object as an NWB object: the attributes `neurodata_type`, `namespace` and `object_id` (a new UUID), and further attributes given as keyword arguments (e.g. `description`); strings are written as variable-length UTF-8 strings (as by pynwb).

# Arguments

- `obj`: HDF5 group or dataset
- `neurodata_type::String`
- `namespace::String="core"`
- `attrs...`: further attributes

# Returns

- `Nothing`
"""
function _nwb_attrs!(obj, neurodata_type::String; namespace::String = "core", attrs...)::Nothing
    # (strings: variable-length UTF-8, as written by pynwb)
    _h5_vstr_attr!(obj, "neurodata_type", neurodata_type)
    _h5_vstr_attr!(obj, "namespace", namespace)
    _h5_vstr_attr!(obj, "object_id", _nwb_uuid())
    for (k, v) in attrs
        if v isa AbstractString
            _h5_vstr_attr!(obj, string(k), v)
        else
            HDF5.attributes(obj)[string(k)] = v
        end
    end
    return nothing
end

"""
    _nwb_dataset!(parent, name, x; <keyword arguments>)

Write an array as a dataset; with `compress > 0`, chunked and compressed with gzip (level `compress`), optionally with the shuffle filter. A Julia array of size `(channels, time points)` is stored as NWB `(time, channels)`.

# Arguments

- `parent`: HDF5 file or group
- `name::String`
- `x::AbstractArray`
- `compress::Int64=0`: gzip level (0–9; 0: not compressed)
- `shuffle::Bool=false`: use the shuffle filter (with compression)
- `chunks::Union{Nothing, Tuple}=nothing`: chunk size (in Julia dimensions); `nothing`: chosen automatically (about 1 MB per chunk)

# Returns

- `HDF5.Dataset`
"""
function _nwb_dataset!(
    parent,
    name::String,
    x::AbstractArray;
    compress::Int64 = 0,
    shuffle::Bool = false,
    chunks::Union{Nothing, Tuple} = nothing,
)
    0 <= compress <= 9 || throw(ArgumentError("compress must be in [0, 9]."))
    x = collect(x)
    if compress == 0 || isempty(x)
        parent[name] = x
        return parent[name]
    end
    if isnothing(chunks)
        # whole first dimension (channels), ~1 MB along the last one
        per = max(1, sizeof(eltype(x)) * prod(size(x)[1:(end - 1)]))
        chunks = (size(x)[1:(end - 1)]..., clamp(2^20 ÷ per, 1, size(x, ndims(x))))
    end
    ds = HDF5.create_dataset(parent, name, eltype(x), size(x); chunk = chunks,
                             shuffle = shuffle, deflate = compress)
    ds[ntuple(_ -> Colon(), ndims(x))...] = x
    return ds
end

# write Bool values as the 8-bit enumeration FALSE = 0, TRUE = 1 (h5py/pynwb Bool type)
function _nwb_bool!(parent, name::String, x::AbstractArray{Bool})
    # (H5Tenum_create is not wrapped by HDF5.jl)
    id = ccall((:H5Tenum_create, HDF5.API.libhdf5), HDF5.API.hid_t, (HDF5.API.hid_t,),
               HDF5.API.H5T_NATIVE_INT8)
    id < 0 && error("H5Tenum_create failed.")
    dt = HDF5.Datatype(id)
    try
        for (nm, v) in (("FALSE", Int8(0)), ("TRUE", Int8(1)))
            HDF5.API.h5t_enum_insert(dt, nm, Ref(v))
        end
        ds = HDF5.create_dataset(parent, name, dt, HDF5.dataspace(size(x)))
        vals = Int8.(collect(x))
        HDF5.API.h5d_write(ds, dt, HDF5.API.H5S_ALL, HDF5.API.H5S_ALL, HDF5.API.H5P_DEFAULT, vals)
        return ds
    finally
        close(dt)
    end
end

# write the object references to the given paths (e.g. the electrode `group` column)
function _nwb_refs!(f::HDF5.File, parent, name::String, paths::AbstractVector{<:AbstractString})
    parent[name] = [HDF5.Reference(f, string(p)) for p in paths]
    return parent[name]
end

# attribute holding a reference to the object at `path` (e.g. `table` of a DynamicTableRegion)
function _nwb_ref_attr!(f::HDF5.File, obj, name::String, path::AbstractString)::Nothing
    HDF5.attributes(obj)[name] = HDF5.Reference(f, string(path))
    return nothing
end

"""
    _nwb_table!(f, parent, name, df; <keyword arguments>)

Write a DataFrame as an NWB DynamicTable (or a subtype, e.g. `ElectrodesTable`, `TimeIntervals`).

- columns of vectors (e.g. tags): ragged columns with a `<column>_index` VectorIndex
- columns listed in `regions` (column => path of the target table): DynamicTableRegion columns (1-based row numbers, stored 0-based)
- columns listed in `references`: references to the objects at the given paths (e.g. electrode groups)
- Bool columns: 8-bit enumeration; String columns: variable-length UTF-8 strings
- the `id` column (if present) is written as the ElementIdentifiers `id`, otherwise 0, 1, …

# Arguments

- `f::HDF5.File`
- `parent`: HDF5 file or group
- `name::String`
- `df::DataFrame`
- `neurodata_type::String="DynamicTable"`
- `namespace::String="hdmf-common"`
- `description::String=""`
- `col_descriptions::AbstractDict=Dict()`: descriptions of the columns
- `regions::AbstractDict=Dict()`: DynamicTableRegion columns and their target tables
- `references::Vector{String}=String[]`: columns of object paths written as references

# Returns

- `HDF5.Group`
"""
function _nwb_table!(
    f::HDF5.File,
    parent,
    name::String,
    df::DataFrame;
    neurodata_type::String = "DynamicTable",
    namespace::String = "hdmf-common",
    description::String = "",
    col_descriptions::AbstractDict = Dict{String, String}(),
    regions::AbstractDict = Dict{String, String}(),
    references::Vector{String} = String[],
)
    g = HDF5.create_group(parent, name)
    cols = [c for c in names(df) if c != "id"]
    _nwb_attrs!(g, neurodata_type; namespace = namespace, description = description,
                colnames = cols)
    ids = "id" in names(df) ? Int64.(df.id) : collect(0:(DataFrames.nrow(df) - 1))
    g["id"] = ids
    _nwb_attrs!(g["id"], "ElementIdentifiers"; namespace = "hdmf-common")
    for c in cols
        v = df[!, c]
        desc = string(get(col_descriptions, c, c))
        if haskey(regions, c)
            ds = _nwb_dataset!(g, c, Int64.(v) .- 1)
            _nwb_attrs!(ds, "DynamicTableRegion"; namespace = "hdmf-common", description = desc)
            _nwb_ref_attr!(f, ds, "table", regions[c])
        elseif c in references
            ds = _nwb_refs!(f, g, c, string.(v))
            _nwb_attrs!(ds, "VectorData"; namespace = "hdmf-common", description = desc)
        elseif eltype(v) <: AbstractVector
            # ragged column: concatenated values and end positions
            vals = reduce(vcat, v; init = eltype(eltype(v))[])
            ds = eltype(vals) <: Bool ? _nwb_bool!(g, c, vals) :
                 eltype(vals) <: AbstractString ? (g[c] = string.(vals); g[c]) : _nwb_dataset!(g, c, vals)
            _nwb_attrs!(ds, "VectorData"; namespace = "hdmf-common", description = desc)
            ix = _nwb_dataset!(g, c * "_index", cumsum(Int64.(length.(v))))
            _nwb_attrs!(ix, "VectorIndex"; namespace = "hdmf-common",
                        description = "Index for VectorData '$c'")
            _nwb_ref_attr!(f, ix, "target", HDF5.name(ds))
        else
            ds = eltype(v) <: Bool ? _nwb_bool!(g, c, Vector{Bool}(v)) :
                 eltype(v) <: AbstractString ? (g[c] = string.(v); g[c]) :
                 _nwb_dataset!(g, c, collect(v))
            _nwb_attrs!(ds, "VectorData"; namespace = "hdmf-common", description = desc)
        end
    end
    return g
end
