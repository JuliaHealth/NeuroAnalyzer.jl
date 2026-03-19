export cbp

"""
    cbp(s; <keyword arguments>)

Perform convolution band-pass filtering at a single frequency.

Generates a sine-wave kernel at `frq` Hz over a ±1 s window and convolves it with the signal. This acts as a narrow band-pass filter centred at `frq`.

# Arguments

- `s::AbstractVector`: signal vector
- `pad::Int64=0`: number of zeros to append before filtering; must be ≥ 0
- `frq::Real`: filter centre frequency in Hz; must be in `(0, fs/2]`
- `fs::Int64`: sampling rate in Hz; must be ≥ 1

# Returns

- `Vector{Float64}`: band-pass filtered signal of the same length as `s` (padding is removed after convolution)

# Throws

- `ArgumentError`: if `fs < 1`, `pad < 0`, `frq ≤ 0`, or `frq > fs/2`

# See also

[`cbp(::NeuroAnalyzer.NEURO)`](@ref), [`tconv`](@ref)
"""
function cbp(
    s::AbstractVector;
    pad::Int64 = 0,
    frq::Real,
    fs::Int64
)::Vector{Float64}

    fs >= 1 || throw(ArgumentError("fs must be ≥ 1."))
    pad >= 0 || throw(ArgumentError("pad must be ≥ 0."))
    frq > 0 || throw(ArgumentError("frq must be > 0."))
    frq <= fs / 2 || throw(ArgumentError("frq must be ≤ $(fs / 2) Hz."))

    pad > 0 && (s = pad0(s, pad))
    kernel = generate_sine(frq, -1:(1 / fs):1)

    return tconv(s, kernel=kernel)

end

"""
    cbp(obj; <keyword arguments>)

Perform convolution band-pass filtering on selected channels of a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `pad::Int64=0`: number of zeros to append before filtering; must be ≥ 0
- `frq::Real`: filter centre frequency in Hz; must be in `(0, fs/2]`

# Returns

- `NeuroAnalyzer.NEURO`: new object with filtered channels

# Throws

- `ArgumentError`: if `frq` or `pad` are out of range

# See also

[`cbp!`](@ref), [`cbp(::AbstractVector)`](@ref)
"""
function cbp(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    pad::Int64 = 0,
    frq::Real
)::NeuroAnalyzer.NEURO

    # resolve channel names to integer indices
    ch = get_channel(obj, ch=ch)

    # number of channels
    ch_n = length(ch)
    # number of epochs
    ep_n = nepochs(obj)
    # sampling rate
    fs = sr(obj)

    obj_new = deepcopy(obj)

    # calculate over channel and epochs
    @inbounds Threads.@threads :static for idx in CartesianIndices((ch_n, ep_n))
        ch_idx, ep_idx = idx[1], idx[2]
        obj_new.data[ch[ch_idx], :, ep_idx] = @views cbp(
            obj_new.data[ch[ch_idx], :, ep_idx], pad = pad, frq = frq, fs = fs
        )
    end

    push!(obj_new.history, "cbp(OBJ, ch=$ch, pad=$pad, frq=$frq)")

    return obj_new

end

"""
    cbp!(obj; <keyword arguments>)

Perform convolution band-pass filtering in-place on selected channels of a NEURO object.

# Arguments

- `obj::NeuroAnalyzer.NEURO`: input NEURO object; modified in-place
- `ch::Union{String, Vector{String}, Regex}`: channel name(s)
- `pad::Int64=0`: number of zeros to append before filtering; must be ≥ 0
- `frq::Real`: filter centre frequency in Hz; must be in `(0, fs/2]`

# Returns

- `Nothing`

# Throws

- `ArgumentError`: if `frq` or `pad` are out of range

# See also

[`cbp`](@ref)
"""
function cbp!(
    obj::NeuroAnalyzer.NEURO;
    ch::Union{String, Vector{String}, Regex},
    pad::Int64 = 0,
    frq::Real
)::Nothing

    obj_new = cbp(obj, ch = ch, pad = pad, frq = frq)
    obj.data = obj_new.data
    obj.history = obj_new.history

    return nothing

end
