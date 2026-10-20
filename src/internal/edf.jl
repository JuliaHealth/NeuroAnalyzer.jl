# ---------------------------------------------------------------------------
# shared helpers for record-based formats (EDF/EDF+, BDF/BDF+, Alice 4 EDF, GDF),
# for channel units read from file headers (also used by BrainVision) and for
# converting non-integer sampling rates (NeuroAnalyzer stores integer rates)
# ---------------------------------------------------------------------------

# relative tolerance below which a sampling rate is treated as an integer
# (rates estimated from time stamps, e.g. 256.0016 Hz, are not resampled)
const _RATE_RTOL = 1.0e-4

_is_int_rate(fs::Real)::Bool = abs(fs - round(fs)) <= _RATE_RTOL * fs

"""
    _resample_mirror(x, m)

Fourier-resample `x` to `m` samples. The signal is mirrored before the FFT, so that its periodic extension has no jump between the last and the first sample (which would ring at both edges). Upsampling by an integer factor keeps the original samples.
"""
function _resample_mirror(x::AbstractVector{<:Real}, m::Int64)::Vector{Float64}
    n = length(x)
    n == m && return Float64.(x)
    n == 1 && return fill(Float64(x[1]), m)
    return FourierTools.resample(vcat(Float64.(x), reverse(Float64.(x))), 2 * m)[1:m]
end

"""
    _resample_nearest(x, m)

Resample `x` to `m` samples by taking the nearest original sample (for trigger / marker channels, whose codes must not be interpolated).
"""
function _resample_nearest(x::AbstractVector{<:Real}, m::Int64)::Vector{Float64}
    n = length(x)
    n == m && return Float64.(x)
    # sample-and-hold (floor, in integer arithmetic): new sample j takes the last original sample at or
    # before its time, so a trigger onset stays at its original time (rounding would move it up to
    # half an original sample earlier)
    return Float64.(x[clamp.(div.((0:(m - 1)) .* n, m) .+ 1, 1, n)])
end

"""
    _integer_rate(data, fs; <keyword arguments>)

Convert data sampled at `fs` Hz to an integer sampling rate.

- `fs` within 0.01% of an integer: data unchanged, rate rounded
- otherwise, if `k × fs` is an integer for some `k ≤ kmax` (e.g. 12.5 Hz × 2 = 25 Hz): exact upsampling by the smallest such `k`; the original samples are kept
- otherwise: resampling to the nearest integer rate (approximate)

Signals are Fourier-resampled (see [`_resample_mirror`](@ref)); channels in `nearest` (trigger / marker channels) are resampled by nearest sample.

# Arguments

- `data::AbstractArray{<:Real, 3}`: channels × samples × epochs
- `fs::Real`: sampling rate in Hz; must be > 0
- `kmax::Int64=16`: largest upsampling factor tried for exact upsampling
- `nearest::Vector{Int64}=Int64[]`: channels resampled by nearest sample

# Returns

- `Array{Float64, 3}`: data at the integer rate
- `Int64`: integer sampling rate
"""
function _integer_rate(
    data::AbstractArray{<:Real, 3},
    fs::Real;
    kmax::Int64 = 16,
    nearest::Vector{Int64} = Int64[],
)::Tuple{Array{Float64, 3}, Int64}
    fs > 0 || throw(ArgumentError("Sampling rate must be > 0."))
    _is_int_rate(fs) && return Float64.(data), round(Int64, fs)

    n = size(data, 2)
    ks = 2:kmax
    k_idx = findfirst(k -> _is_int_rate(k * fs), ks)
    if !isnothing(k_idx)
        k = ks[k_idx]
        new_fs = round(Int64, k * fs)
        m = n * k
        _info("Sampling rate $fs Hz is not an integer: data upsampled ×$k to $new_fs Hz (exact)")
    else
        new_fs = max(1, round(Int64, fs))
        m = round(Int64, n * new_fs / fs)
        _warn(
            "Sampling rate $fs Hz is not an integer: data resampled to $new_fs Hz (approximate).",
        )
    end

    out = zeros(size(data, 1), m, size(data, 3))
    @inbounds for ep_idx in axes(data, 3), ch_idx in axes(data, 1)
        x = @view data[ch_idx, :, ep_idx]
        out[ch_idx, :, ep_idx] =
            ch_idx in nearest ? _resample_nearest(x, m) : _resample_mirror(x, m)
    end

    return out, new_fs
end

"""
    _normalize_unit(u)

Normalise a physical unit string read from a file header (e.g. `"uV"`, `"µV"` → `"μV"`).
"""
function _normalize_unit(u::AbstractString)::String
    s = strip(replace(u, '\0' => ""))
    l = lowercase(s)
    l in ("uv", "µv", "μv") && return "μV"
    l == "mv" && return "mV"
    l == "nv" && return "nV"
    l == "v" && return "V"
    return String(s)
end

"""
    _resolve_units(units, ch_type)

Resolve channel units from the units stored in a file header.

- channels whose NeuroAnalyzer unit is μV (EEG, EOG, EMG, …) stored in V, mV or nV are converted to μV
- empty units fall back to the default unit of the channel type; marker channels are unitless
- any other unit is kept as stored in the file (e.g. an ECG channel recorded in μV stays in μV)

# Returns

- `Vector{String}`: channel units
- `Vector{Float64}`: factor by which the data of each channel must be multiplied
"""
function _resolve_units(
    units::AbstractVector{<:AbstractString},
    ch_type::AbstractVector{<:AbstractString},
)::Tuple{Vector{String}, Vector{Float64}}
    length(units) == length(ch_type) ||
        throw(ArgumentError("units and ch_type must have the same length."))
    to_uv = Dict("V" => 1.0e6, "mV" => 1.0e3, "μV" => 1.0, "nV" => 1.0e-3)
    u_out = Vector{String}(undef, length(units))
    scale = ones(length(units))
    for idx in eachindex(units)
        def = _ch_units(ch_type[idx])
        u = _normalize_unit(units[idx])
        if ch_type[idx] == "mrk" || isempty(u)
            u_out[idx] = def
        elseif def == "μV" && haskey(to_uv, u)
            scale[idx] = to_uv[u]
            u_out[idx] = "μV"
        else
            u_out[idx] = u
        end
    end
    return u_out, scale
end

"""
    _phys_scale!(x, dmin, dmax, pmin, pmax)

Convert digital values to physical values in place: `(x - dmin) × gain + pmin`, with `gain = (pmax - pmin) / (dmax - dmin)`.

# Returns

- `Float64`: gain
"""
function _phys_scale!(
    x::AbstractVector{Float64},
    dmin::Real,
    dmax::Real,
    pmin::Real,
    pmax::Real,
)::Float64
    # a degenerate digital range cannot be scaled; keep the raw values
    dmax == dmin && return 1.0
    gain = (pmax - pmin) / (dmax - dmin)
    @. x = (x - dmin) * gain + pmin
    return gain
end

"""
    _record_count(data_records, available_bytes, record_bytes, file_name)

Number of complete data records to read. A negative header value (unknown, e.g. Alice 4 or a recording that was not closed properly) is replaced by the number of records present in the file; a header value larger than that (truncated file) is reduced with a warning.
"""
function _record_count(
    data_records::Integer,
    available_bytes::Integer,
    record_bytes::Integer,
    file_name::String,
)::Int64
    record_bytes > 0 || throw(ArgumentError("$file_name: data record size is 0."))
    n_avail = available_bytes ÷ record_bytes
    data_records < 0 && return n_avail
    if data_records > n_avail
        _warn("$file_name is truncated: header declares $data_records data records, $n_avail found.")
        return n_avail
    end
    return data_records
end

# 24-bit little-endian two's-complement → signed integer
@inline _int24(b0::UInt8, b1::UInt8, b2::UInt8)::Int32 =
    reinterpret(Int32, (UInt32(b0) << 8) | (UInt32(b1) << 16) | (UInt32(b2) << 24)) >> 8

"""
    _read_edf_records(raw, spr, n_rec, bps; <keyword arguments>)

Split the data block of a record-based file into per-channel sample vectors.

Each data record holds, for every channel in turn, `spr[ch]` samples of `bps` bytes (2: little-endian Int16, as in EDF; 3: little-endian 24-bit, as in BDF).

# Arguments

- `raw::AbstractVector{UInt8}`: data block (`n_rec` complete records)
- `spr::Vector{Int64}`: samples per data record for each channel
- `n_rec::Int64`: number of data records
- `bps::Int64`: bytes per sample (2 or 3)
- `annotation_channels::Vector{Int64}=Int64[]`: channels holding TAL text; their bytes are returned as strings
- `trigger_channels::Vector{Int64}=Int64[]`: 24-bit channels whose lower 16 bits are trigger inputs (BDF Status)

# Returns

- `Vector{Vector{Float64}}`: digital samples of each channel (empty for annotation channels)
- `Vector{String}`: annotation text of each record and annotation channel
"""
function _read_edf_records(
    raw::AbstractVector{UInt8},
    spr::Vector{Int64},
    n_rec::Int64,
    bps::Int64;
    annotation_channels::Vector{Int64} = Int64[],
    trigger_channels::Vector{Int64} = Int64[],
)::Tuple{Vector{Vector{Float64}}, Vector{String}}
    bps in (2, 3) || throw(ArgumentError("bps must be 2 or 3."))
    ch_n = length(spr)
    sig = [ch in annotation_channels ? Float64[] : zeros(spr[ch] * n_rec) for ch = 1:ch_n]
    annotations = String[]

    pos = 0
    @inbounds for rec = 0:(n_rec - 1), ch = 1:ch_n
        n = spr[ch]
        nb = n * bps
        chunk = view(raw, (pos + 1):(pos + nb))
        pos += nb
        if ch in annotation_channels
            push!(annotations, String(Char.(chunk)))
            continue
        end
        dst = view(sig[ch], (rec * n + 1):((rec + 1) * n))
        if bps == 2
            dst .= ltoh.(reinterpret(Int16, chunk))
        elseif ch in trigger_channels
            # BDF Status: bits 0-15 are the trigger inputs, bits 16-23 status flags
            for i = 1:n
                j = 3 * (i - 1)
                dst[i] = Float64(UInt16(chunk[j + 1]) | (UInt16(chunk[j + 2]) << 8))
            end
        else
            for i = 1:n
                j = 3 * (i - 1)
                dst[i] = Float64(_int24(chunk[j + 1], chunk[j + 2], chunk[j + 3]))
            end
        end
    end

    return sig, annotations
end

"""
    _merge_channels(sig, skip)

Assemble per-channel sample vectors into a `(channels × samples × 1)` array. Channels in `skip` (e.g. annotation channels) are left as zeros; signal channels with fewer samples than the longest one are upsampled to its length (Fourier resampling of the whole channel, see [`_resample_mirror`](@ref); channels in `nearest`, e.g. trigger channels, by nearest sample).

# Returns

- `Array{Float64, 3}`: data
- `Bool`: `true` if any channel was resampled
"""
function _merge_channels(
    sig::Vector{Vector{Float64}},
    skip::Vector{Int64} = Int64[];
    nearest::Vector{Int64} = Int64[],
)::Tuple{Array{Float64, 3}, Bool}
    chs = setdiff(eachindex(sig), skip)
    isempty(chs) && throw(ArgumentError("No signal channels found."))
    n = maximum(length(sig[ch]) for ch in chs)
    d = zeros(length(sig), n, 1)
    resampled = false
    for ch in chs
        x = sig[ch]
        if length(x) != n
            x = ch in nearest ? _resample_nearest(x, n) : _resample_mirror(x, n)
            resampled = true
        end
        d[ch, :, 1] = x
    end
    return d, resampled
end
