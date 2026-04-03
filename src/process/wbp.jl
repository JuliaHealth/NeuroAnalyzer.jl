export wbp
export wbp!

"""
    wbp(s; <keyword arguments>)

Perform wavelet band-pass filtering.

# Arguments

- `s::AbstractVector`: signal vector
- `pad::Int64=0`: number of zeros to append
- `frq::Real`: filter frequency
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `ncyc::Int64=6`: Morlet wavelet cycles

# Returns

- `Vector{Float64}`
"""
function wbp(
    s::AbstractVector;
    pad::Int64 = 0,
    frq::Real,
    fs::Int64,
    ncyc::Int64 = 6,
)::Vector{Float64}
    # validate
    fs >= 1 || throw(ArgumentError("fs must be ≥ 1."))
    frq > 0 || throw(ArgumentError("frq must be > 0."))
    ncyc > 0 || throw(ArgumentError("ncyc must be > 0."))
    pad >= 0 || throw(ArgumentError("pad must be ≥ 0."))
    frq <= fs / 2 || throw(ArgumentError("frq must be ≤ $(fs / 2)."))

    pad > 0 && (s = pad0(s, pad))

    kernel = generate_morlet(fs, frq, 1; ncyc = ncyc, complex = true)

    return real.(fconv(s; kernel = kernel, norm = true))
end

"""
    wbp(s; <keyword arguments>)

Perform wavelet band-pass filtering for a 3-D signal array.

# Arguments

- `s::AbstractArray`: signal array, shape (channels, samples, epochs)
- `pad::Int64=0`: pad the `signal` with `pad` zeros
- `frq::Real`: filter frequency
- `fs::Int64`: sampling rate in Hz; must be ≥ 1
- `ncyc::Int64=6`: Morlet wavelet cycles

# Returns

- `Array{Float64, 3}`
"""
function wbp(
    s::AbstractArray;
    pad::Int64 = 0,
    frq::Real,
    fs::Int64,
    ncyc::Int64 = 6,
)::Array{Float64, 3}
    # validate that the input is a proper 3-D array (channels, samples, epochs)
    _chk3d(s)

    # number of channels
    ch_n = size(s, 1)
    # number of epochs
    ep_n = size(s, 3)

    # pre-allocate output
    s_new = similar(s, Float64)

    # calculate over channels and epochs
    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        s_new[ch_idx, :, ep_idx] = wbp(
            @view(s[ch_idx, :, ep_idx]),
            pad = pad,
            frq = frq,
            fs = fs,
            ncyc = ncyc,
        )
    end

    return s_new
end

"""
    wbp(obj; <keyword arguments>)

Perform wavelet band-pass filtering.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `pad::Int64=0`: pad the `signal` with `pad` zeros
- `frq::Real`: filter frequency
- `ncyc::Int64=6`: Morlet wavelet cycles

# Returns

- `NeuroAnalyzer.NEURO`: output NEURO object
"""
function wbp(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    pad::Int64 = 0,
    frq::Real,
    ncyc::Int64 = 6,
)::NeuroAnalyzer.NEURO
    # resolve channel names to integer indices
    ch = get_channel(obj; ch = ch)
    isempty(ch) && throw(ArgumentError("No channels selected."))

    # create new dataset
    obj_new = deepcopy(obj)

    obj_new.data[ch, :, :] =
        wbp(@view(obj.data[ch, :, :]); pad = pad, frq = frq, fs = sr(obj), ncyc = ncyc)
    push!(obj_new.history, "wbp(obj; ch=$ch, pad=$pad, frq=$frq, ncyc=$ncyc)")

    return obj_new
end

"""
    wbp!(obj; <keyword arguments>)

Perform wavelet band-pass filtering.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `pad::Int64=0`: pad the `signal` with `pad` zeros
- `frq::Real`: filter frequency
- `ncyc::Int64=6`: Morlet wavelet cycles

# Returns

- `Nothing`
"""
function wbp!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    pad::Int64 = 0,
    frq::Real,
    ncyc::Int64 = 6,
)::Nothing
    obj_new = wbp(obj; ch = ch, pad = pad, frq = frq, ncyc = ncyc)
    obj.data = obj_new.data
    obj.history = obj_new.history

    return nothing
end
