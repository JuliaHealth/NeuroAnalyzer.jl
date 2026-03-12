export wbp
export wbp!

"""
    wbp(s; <keyword arguments>)

Perform wavelet band-pass filtering.

# Arguments

- `s::AbstractVector`: signal vector
- `pad::Int64`: number of zeros to append
- `frq::Real`: filter frequency
- `fs::Int64`: sampling rate
- `ncyc::Int64=6`: Morlet wavelet cycles

# Returns

- `s_new::Vector{Float64}`
"""
function wbp(s::AbstractVector; pad::Int64 = 0, frq::Real, fs::Int64, ncyc::Int64 = 6)::Vector{Float64}

    @assert fs >= 1 "fs must be ≥ 1."
    @assert frq > 0 "frq must be > 0."
    @assert ncyc > 0 "ncyc must be > 0."
    @assert pad >= 0 "pad must be ≥ 0."
    @assert frq <= fs / 2 "frq must be ≤ $(fs / 2)."

    pad > 0 && (s = pad0(s, pad))

    kernel = generate_morlet(fs, frq, 1; ncyc = ncyc, complex = true)

    s_new = real.(fconv(s; kernel = kernel, norm = true))

    return s_new

end

"""
    wbp(s; <keyword arguments>)

Perform wavelet band-pass filtering.

# Arguments

- `s::AbstractArray`
- `pad::Int64`: pad the `signal` with `pad` zeros
- `frq::Real`: filter frequency
- `fs::Int64`: sampling rate
- `ncyc::Int64=6`: Morlet wavelet cycles

# Returns

- `s_new::Array{Float64, 3}`
"""
function wbp(s::AbstractArray; pad::Int64 = 0, frq::Real, fs::Int64, ncyc::Int64 = 6)::Array{Float64, 3}

    _chk3d(s)
    ch_n = size(s, 1)
    ep_n = size(s, 3)

    s_new = similar(s)
    @inbounds for ep_idx in 1:ep_n
        Threads.@threads :dynamic for ch_idx in 1:ch_n
            s_new[ch_idx, :, ep_idx] = @views wbp(s[ch_idx, :, ep_idx], pad = pad, frq = frq, fs = fs, ncyc = ncyc)
        end
    end

    return s_new

end

"""
    wbp(obj; <keyword arguments>)

Perform wavelet band-pass filtering.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `pad::Int64`: pad the `signal` with `pad` zeros
- `frq::Real`: filter frequency
- `ncyc::Int64=6`: Morlet wavelet cycles

# Returns

- `obj_new::NeuroAnalyzer.NEURO`: output NEURO object
"""
function wbp(
        obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, pad::Int64 = 0, frq::Real, ncyc::Int64 = 6
    )::NeuroAnalyzer.NEURO

    ch = get_channel(obj, ch = ch)
    obj_new = deepcopy(obj)
    obj_new.data[ch, :, :] = @views wbp(obj.data[ch, :, :], pad = pad, frq = frq, fs = sr(obj), ncyc = ncyc)
    push!(obj_new.history, "wbp(OBJ, ch=$ch, pad=$pad, frq=$frq, ncyc=$ncyc)")

    return obj_new

end

"""
    wbp!(obj; <keyword arguments>)

Perform wavelet band-pass filtering.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `pad::Int64`: pad the `signal` with `pad` zeros
- `frq::Real`: filter frequency
- `ncyc::Int64=6`: Morlet wavelet cycles

# Returns

- `Nothing`
"""
function wbp!(
        obj::NeuroAnalyzer.NEURO; ch::Union{String, Vector{String}, Regex}, pad::Int64 = 0, frq::Real, ncyc::Int64 = 6
    )::Nothing

    obj_new = wbp(obj, ch = ch, pad = pad, frq = frq, ncyc = ncyc)
    obj.data = obj_new.data
    obj.history = obj_new.history

    return nothing

end
